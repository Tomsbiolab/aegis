"""
Peak-to-gene assignment for ChIP-seq / DAP-seq peaks.

Peaks already called by a peak caller (GEM, GPS, MACS2...) are read from a
GEM events table (``*.GEM_events.txt``), a narrowPeak file or a plain BED file
and linked to the genes of an :class:`~aegis.annotation.Annotation`:

* ``mode="nearest"``: every peak is linked to its closest gene(s). Distance is
  measured either from the peak summit to the TSS (``reference="tss"``, as the
  ChIPpeakAnno ``nearestLocation`` default and HOMER) or to the gene body
  (``reference="gene"``, 0 when the peak overlaps the gene).
* ``mode="window"``: every peak is linked to all genes whose regulatory window
  contains the summit, i.e. from ``upstream`` bp before the TSS to
  ``downstream`` bp after the TES (strand aware).

Besides the distances AEGIS uses, every link gets the location of the whole peak
relative to its gene and the shortest distance between their ends, as
ChIPpeakAnno's ``insideFeature`` and ``shortestDistance`` (``inside_feature``,
``shortest_distance``). Peak values given as -1 (the narrowPeak convention for
"not available") are read as missing. narrowPeak and BED files can carry a fold
enrichment and a binding motif (e.g. GEM's ``Fold`` and ``KmerGroup``) in extra
columns named ``fold_enrichment`` and ``binding_motif`` in a ``#`` header line.

All coordinates in the output are 1-based and inclusive, as in GFF.
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Iterable, Literal

if TYPE_CHECKING:
    from .annotation import Annotation

import gzip
import math
import re
from bisect import bisect_left, bisect_right
from dataclasses import dataclass, field

import pandas as pd


PEAK_FORMATS = ("auto", "gem", "narrowpeak", "bed")

# GEM event positions look like "00:5064740" or "chr1:1234"
_GEM_POSITION_RE = re.compile(r"^(\S+):(\d+)$")


@dataclass
class Peak:
    """A peak in 1-based inclusive coordinates."""
    chrom: str
    start: int
    end: int
    summit: int
    name: str = ""
    score: float = math.nan
    signal: float = math.nan
    fold: float = math.nan
    p_log10: float = math.nan
    q_log10: float = math.nan
    motif: str = ""
    peak_id: str = ""


@dataclass
class _GeneRecord:
    id: str
    ch: str
    start: int
    end: int
    strand: str
    tss: int = field(init=False)
    tes: int = field(init=False)

    def __post_init__(self):
        if self.strand == "-":
            self.tss, self.tes = self.end, self.start
        else:
            self.tss, self.tes = self.start, self.end


# ---------------------------------------------------------------------------
# Reading peaks
# ---------------------------------------------------------------------------

def _open(path:str):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return open(path, "rt", encoding="utf-8")


def _to_float(value:str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _peak_value(value:str) -> float:
    """A narrowPeak / BED value: -1 means not available."""
    number = _to_float(value)
    return math.nan if number == -1 else number


# Extra columns of narrowPeak / BED files, found by their name in a '#' header line
_EXTRA_COLUMNS = {"fold": ("fold_enrichment", "fold"), "motif": ("binding_motif", "motif")}


def _header_columns(path:str) -> list[str]:
    """Column names of a '#' header line before the first data line (empty without one)."""
    names = []
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith("#"):
                names = [c.strip().lower() for c in line.lstrip("#").split("\t")]
            elif line.strip() and not line.startswith(("track", "browser")):
                break
    return names


def _is_int(value:str) -> bool:
    return value.strip().lstrip("-").isdigit()


def _data_lines(path:str):
    with _open(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            yield line


def detect_peak_format(path:str) -> str:
    """Return 'gem', 'narrowpeak' or 'bed' for a peak file."""
    for line in _data_lines(path):
        cols = line.split("\t")
        if cols[0].strip() == "Position" and "IP" in [c.strip() for c in cols]:
            return "gem"
        if len(cols) >= 10 and all(_is_int(cols[i]) for i in (1, 2, 9)):
            return "narrowpeak"
        if len(cols) >= 3 and _is_int(cols[1]) and _is_int(cols[2]):
            return "bed"
        if _GEM_POSITION_RE.match(cols[0].strip()):
            return "gem"
        raise ValueError(f"Could not detect the format of peak file '{path}' from line: {line!r}")
    raise ValueError(f"Peak file '{path}' has no data lines.")


def _read_gem(path:str, event_width:int) -> list[Peak]:
    peaks = []
    columns = None
    # 200 bp around an event at P -> [P-99, P+100], the interval of GEM's own BED/narrowPeak files
    width = max(1, event_width)
    half_left = (width - 1) // 2
    half_right = width - 1 - half_left
    for line in _data_lines(path):
        cols = [c.strip() for c in line.split("\t")]
        if cols[0] == "Position":
            columns = {name: i for i, name in enumerate(cols)}
            continue
        m = _GEM_POSITION_RE.match(cols[0])
        if not m:
            raise ValueError(f"Unexpected GEM event position '{cols[0]}' in '{path}'.")
        chrom, pos = m.group(1), int(m.group(2))
        def text(name, default_index):
            i = columns.get(name, default_index) if columns else default_index
            return cols[i] if i is not None and i < len(cols) else ""
        def col(name, default_index):
            return _to_float(text(name, default_index))
        peaks.append(Peak(
            chrom=chrom,
            start=max(1, pos - half_left),
            end=pos + half_right,
            summit=pos,
            name=cols[0],
            signal=col("IP", 1),
            fold=col("Fold", 3),
            q_log10=col("Q_-lg10", 5),
            p_log10=col("P_-lg10", 6),
            motif=text("KmerGroup", 10),
        ))
    return peaks


def _read_bed_like(path:str, narrowpeak:bool) -> list[Peak]:
    peaks = []
    seen = set()
    # extra columns (fold enrichment, motif) named in a header line, after the narrowPeak / BED6 columns
    header = _header_columns(path)
    first_extra = 10 if narrowpeak else 6
    extra = {}
    for field_name, names in _EXTRA_COLUMNS.items():
        found = [i for i, c in enumerate(header) if c in names and i >= first_extra]
        if found:
            extra[field_name] = found[0]
    for line in _data_lines(path):
        cols = [c.strip() for c in line.split("\t")]
        if len(cols) < 3 or not (_is_int(cols[1]) and _is_int(cols[2])):
            # header lines of BED-like tables
            continue
        chrom = cols[0]
        start0, end = int(cols[1]), int(cols[2])
        name = cols[3] if len(cols) > 3 and cols[3] not in ("", ".") else ""
        key = (chrom, start0, end, name)
        # BED files that were already annotated repeat each peak once per gene
        if key in seen:
            continue
        seen.add(key)
        start = start0 + 1
        summit = (start + end) // 2
        score = _peak_value(cols[4]) if len(cols) > 4 else math.nan
        signal = p_log10 = q_log10 = math.nan
        fold = _peak_value(cols[extra["fold"]]) if "fold" in extra and extra["fold"] < len(cols) else math.nan
        motif = cols[extra["motif"]] if "motif" in extra and extra["motif"] < len(cols) else ""
        if motif == ".":
            motif = ""
        if narrowpeak:
            signal = _peak_value(cols[6])
            p_log10 = _peak_value(cols[7])
            q_log10 = _peak_value(cols[8])
            offset = int(cols[9])
            m = _GEM_POSITION_RE.match(name)
            if m and start - 1 <= int(m.group(2)) <= end + 1:
                # GEM narrowPeaks carry the exact event position in the name
                summit = int(m.group(2))
            elif offset >= 0:
                summit = start + offset
        peaks.append(Peak(chrom=chrom, start=start, end=end, summit=summit, name=name,
                          score=score, signal=signal, fold=fold, p_log10=p_log10, q_log10=q_log10, motif=motif))
    return peaks


def read_peaks(path:str, peak_format:str="auto", event_width:int=200) -> list[Peak]:
    """
    Read peaks from a GEM events table, a narrowPeak file or a BED file.

    Parameters
    ----------
    path: peak file (optionally gzipped).
    peak_format: 'auto', 'gem', 'narrowpeak' or 'bed'.
    event_width: width given to GEM events (GEM reports single positions); GEM's
        own BED/narrowPeak files use 200 bp. Set to 0 to treat events as points.
    """
    peak_format = peak_format.lower()
    if peak_format not in PEAK_FORMATS:
        raise ValueError(f"Invalid peak format '{peak_format}'. Choose from: {', '.join(PEAK_FORMATS)}.")
    if peak_format == "auto":
        peak_format = detect_peak_format(path)

    if peak_format == "gem":
        peaks = _read_gem(path, event_width)
    else:
        peaks = _read_bed_like(path, narrowpeak=(peak_format == "narrowpeak"))

    names = [p.name for p in peaks]
    unique_names = all(names) and len(set(names)) == len(names)
    width = len(str(len(peaks)))
    for n, p in enumerate(peaks, start=1):
        p.peak_id = p.name if unique_names else f"peak_{n:0{width}d}"
    return peaks


class ChromosomeMismatchError(ValueError):
    """Peak sequence names that are not present in the annotation."""


_SEQUENCE_REGION_RE = re.compile(r"^##sequence-region\s*(\S+)")


def annotation_sequences(annotation:Annotation) -> set[str]:
    """Sequences of an annotation: those with genes plus those declared in ##sequence-region header lines."""
    sequences = set(annotation.chrs)
    for line in getattr(annotation, "gff_header", []):
        m = _SEQUENCE_REGION_RE.match(line)
        if m:
            sequences.add(m.group(1))
    return sequences


def check_chromosomes(peak_chroms:Iterable[str], annotation_chroms:Iterable[str]) -> None:
    """Raise ChromosomeMismatchError unless every peak sequence name is an annotation sequence name."""
    annotation_chroms = set(annotation_chroms)
    missing = sorted(set(peak_chroms) - annotation_chroms)
    if missing:
        shown = ", ".join(missing[:10]) + (f" (+{len(missing) - 10} more)" if len(missing) > 10 else "")
        examples = ", ".join(sorted(annotation_chroms)[:10]) + (" ..." if len(annotation_chroms) > 10 else "")
        raise ChromosomeMismatchError(
            f"The sequence names of the peak file and the annotation are not the same. "
            f"{len(missing)} sequence name(s) in the peak file are not in the annotation: {shown}. "
            f"Annotation sequences: {examples}."
        )


# ---------------------------------------------------------------------------
# Peak to gene assignment
# ---------------------------------------------------------------------------

class _ChromIndex:
    """Genes of one chromosome sorted by start and by TSS, for fast nearest searches."""

    def __init__(self, genes:list[_GeneRecord]):
        self.by_start = sorted(genes, key=lambda g: (g.start, g.end, g.id))
        self.starts = [g.start for g in self.by_start]
        self.max_len = max((g.end - g.start + 1 for g in genes), default=0)
        self.by_tss = sorted(genes, key=lambda g: (g.tss, g.id))
        self.tss = [g.tss for g in self.by_tss]

    def nearest_by_gene(self, start:int, end:int) -> tuple[list[_GeneRecord], int]:
        best, best_d = [], math.inf
        i = bisect_right(self.starts, end)
        # genes starting at or before the peak end: overlapping or to the left
        j = i - 1
        while j >= 0:
            g = self.by_start[j]
            if g.start + self.max_len - 1 < start - best_d:
                break
            d = max(0, start - g.end)
            if d < best_d:
                best, best_d = [g], d
            elif d == best_d:
                best.append(g)
            j -= 1
        # genes starting after the peak end
        j = i
        while j < len(self.by_start):
            g = self.by_start[j]
            d = g.start - end
            if d > best_d:
                break
            if d < best_d:
                best, best_d = [g], d
            else:
                best.append(g)
            j += 1
        return best, best_d

    def nearest_by_tss(self, summit:int) -> tuple[list[_GeneRecord], int]:
        best, best_d = [], math.inf
        i = bisect_left(self.tss, summit)
        j = i - 1
        while j >= 0:
            d = summit - self.tss[j]
            if d > best_d:
                break
            if d < best_d:
                best, best_d = [self.by_tss[j]], d
            else:
                best.append(self.by_tss[j])
            j -= 1
        j = i
        while j < len(self.tss):
            d = self.tss[j] - summit
            if d > best_d:
                break
            if d < best_d:
                best, best_d = [self.by_tss[j]], d
            else:
                best.append(self.by_tss[j])
            j += 1
        return best, best_d

    def in_window(self, summit:int, upstream:int, downstream:int) -> list[_GeneRecord]:
        flank = max(upstream, downstream)
        lo = bisect_left(self.starts, summit - self.max_len - flank)
        hi = bisect_right(self.starts, summit + flank)
        hits = []
        for g in self.by_start[lo:hi]:
            if g.strand == "-":
                w_start, w_end = g.start - downstream, g.end + upstream
            else:
                w_start, w_end = g.start - upstream, g.end + downstream
            if w_start <= summit <= w_end:
                hits.append(g)
        return hits


def signed_tss_distance(summit:int, gene:_GeneRecord) -> int:
    """Summit position relative to the TSS, following the gene strand (negative = upstream)."""
    if gene.strand == "-":
        return gene.tss - summit
    return summit - gene.tss


def summit_location(summit:int, gene:_GeneRecord) -> str:
    if gene.start <= summit <= gene.end:
        return "genic"
    upstream = summit < gene.start if gene.strand != "-" else summit > gene.end
    return "upstream" if upstream else "downstream"


def inside_feature(start:int, end:int, gene:_GeneRecord) -> str:
    """
    Location of the whole peak relative to the gene, following the gene strand, as ChIPpeakAnno's insideFeature:
    'upstream' or 'downstream' (no overlap), 'inside' (within the gene), 'overlapStart' (over the TSS),
    'overlapEnd' (over the TES or ending at it) or 'includeFeature' (the gene lies within the peak).
    """
    g_start, g_end = gene.start, gene.end
    if gene.strand == "-":
        # mirrored so that the TSS is the gene start
        start, end, g_start, g_end = -end, -start, -g_end, -g_start
    if end < g_start:
        return "upstream"
    if start > g_end:
        return "downstream"
    if start < g_start:
        return "includeFeature" if end > g_end else "overlapStart"
    return "overlapEnd" if end >= g_end else "inside"


def shortest_distance(start:int, end:int, gene:_GeneRecord) -> int:
    """Shortest distance between an end of the peak and an end of the gene, as ChIPpeakAnno's shortestDistance."""
    return min(abs(start - gene.start), abs(start - gene.end), abs(end - gene.start), abs(end - gene.end))


def _gene_records(annotation:Annotation, coding_only:bool=False, skip_pseudogenes:bool=False, skip_te:bool=False) -> dict[str, list[_GeneRecord]]:
    records = {}
    for ch, genes in annotation.chrs.items():
        for g in genes.values():
            if coding_only and not g.coding:
                continue
            if skip_pseudogenes and g.pseudogene:
                continue
            if skip_te and g.transposable:
                continue
            records.setdefault(ch, []).append(_GeneRecord(g.id, ch, g.start, g.end, g.strand))
    return records


OUTPUT_COLUMNS = [
    "peak_id", "peak_name", "chrom", "peak_start", "peak_end", "summit", "score", "signal", "fold", "p_log10", "q_log10",
    "motif", "gene_id", "gene_chrom", "gene_start", "gene_end", "gene_strand",
    "distance_to_tss", "distance_to_gene", "shortest_distance", "summit_location", "inside_feature", "n_linked_genes",
]


def assign_peaks_to_genes(
    annotation:Annotation,
    peaks:list[Peak]|str,
    mode:Literal["nearest", "window"]="nearest",
    reference:Literal["tss", "gene"]="tss",
    max_distance:int|None=None,
    upstream:int=2000,
    downstream:int=1000,
    one_gene_per_peak:bool=False,
    keep_unassigned:bool=False,
    coding_only:bool=False,
    skip_pseudogenes:bool=False,
    skip_te:bool=False,
    peak_format:str="auto",
    event_width:int=200,
    quiet:bool=False,
) -> pd.DataFrame:
    """
    Link peaks to genes of an annotation. Returns one row per peak-gene link.

    Parameters
    ----------
    annotation: AEGIS Annotation object.
    peaks: list of Peak objects or a path to a GEM events / narrowPeak / BED file.
    mode: 'nearest' (closest gene(s) per peak) or 'window' (all genes whose window
        [TSS - upstream, TES + downstream] contains the summit).
    reference: distance used to pick the nearest gene: 'tss' (summit to TSS) or
        'gene' (peak to gene body, 0 if overlapping).
    max_distance: in 'nearest' mode, drop links farther than this (bp, using the
        chosen reference).
    one_gene_per_peak: in 'nearest' mode, break ties (e.g. a peak inside two
        overlapping genes) keeping the gene with the closest TSS.
    keep_unassigned: also output peaks without any gene, with empty gene columns.
    coding_only, skip_pseudogenes, skip_te: restrict the candidate genes.

    Every sequence name in the peaks must be a sequence of the annotation (one with
    genes or declared in a ##sequence-region header line); otherwise
    ChromosomeMismatchError is raised. Peaks on sequences without genes are left
    unassigned.
    """
    if mode not in ("nearest", "window"):
        raise ValueError(f"Invalid mode '{mode}'. Choose from: nearest, window.")
    if reference not in ("gene", "tss"):
        raise ValueError(f"Invalid reference '{reference}'. Choose from: gene, tss.")
    if upstream < 0 or downstream < 0:
        raise ValueError("Upstream and downstream distances must be >= 0.")
    if max_distance is not None and max_distance < 0:
        raise ValueError("Maximum distance must be >= 0.")

    if isinstance(peaks, str):
        peaks = read_peaks(peaks, peak_format=peak_format, event_width=event_width)

    check_chromosomes((p.chrom for p in peaks), annotation_sequences(annotation))

    genes = _gene_records(annotation, coding_only=coding_only, skip_pseudogenes=skip_pseudogenes, skip_te=skip_te)
    indexes = {ch: _ChromIndex(g) for ch, g in genes.items()}

    if not quiet:
        print(f"Read {len(peaks)} peaks on {len({p.chrom for p in peaks})} sequence(s); {sum(len(g) for g in genes.values())} candidate genes.")

    rows = []
    unassigned = 0
    for p in peaks:
        base = [p.peak_id, p.name, p.chrom, p.start, p.end, p.summit, p.score, p.signal, p.fold, p.p_log10, p.q_log10, p.motif]
        index = indexes.get(p.chrom)
        hits = []
        if index is not None:
            if mode == "window":
                hits = index.in_window(p.summit, upstream, downstream)
            else:
                if reference == "gene":
                    hits, dist = index.nearest_by_gene(p.start, p.end)
                else:
                    hits, dist = index.nearest_by_tss(p.summit)
                if max_distance is not None and dist > max_distance:
                    hits = []
                if one_gene_per_peak and len(hits) > 1:
                    hits = [min(hits, key=lambda g: (abs(signed_tss_distance(p.summit, g)), g.id))]
        if not hits:
            unassigned += 1
            if keep_unassigned:
                rows.append(base + ["", "", pd.NA, pd.NA, "", pd.NA, pd.NA, pd.NA, "", "", 0])
            continue
        hits = sorted(hits, key=lambda g: (abs(signed_tss_distance(p.summit, g)), g.id))
        for g in hits:
            d_gene = max(0, g.start - p.end, p.start - g.end)
            rows.append(base + [g.id, g.ch, g.start, g.end, g.strand,
                                signed_tss_distance(p.summit, g), d_gene, shortest_distance(p.start, p.end, g),
                                summit_location(p.summit, g), inside_feature(p.start, p.end, g), len(hits)])

    df = pd.DataFrame(rows, columns=OUTPUT_COLUMNS)
    for col in ("gene_start", "gene_end", "distance_to_tss", "distance_to_gene", "shortest_distance"):
        df[col] = df[col].astype("Int64")

    if not quiet:
        assigned = len(peaks) - unassigned
        print(f"Assigned {assigned} of {len(peaks)} peaks to {df['gene_id'].replace('', pd.NA).nunique()} genes "
              f"({len(df[df['gene_id'] != ''])} peak-gene links, mode={mode}"
              f"{', reference=' + reference if mode == 'nearest' else f', upstream={upstream}, downstream={downstream}'}).")
    return df


def summarise_target_genes(links:pd.DataFrame) -> pd.DataFrame:
    """One row per target gene: number of peaks and the strongest / closest peak."""
    columns = ["gene_id", "gene_chrom", "gene_start", "gene_end", "gene_strand", "n_peaks", "peak_ids",
               "max_signal", "max_score", "min_abs_distance_to_tss", "summit_locations"]
    df = links[links["gene_id"] != ""].copy()
    if df.empty:
        return pd.DataFrame(columns=columns)
    df["abs_tss"] = df["distance_to_tss"].abs()
    grouped = df.groupby("gene_id", sort=False)
    out = grouped.agg(
        gene_chrom=("gene_chrom", "first"),
        gene_start=("gene_start", "first"),
        gene_end=("gene_end", "first"),
        gene_strand=("gene_strand", "first"),
        n_peaks=("peak_id", "nunique"),
        peak_ids=("peak_id", lambda s: ",".join(dict.fromkeys(s))),
        max_signal=("signal", "max"),
        max_score=("score", "max"),
        min_abs_distance_to_tss=("abs_tss", "min"),
        summit_locations=("summit_location", lambda s: ",".join(sorted(set(s)))),
    ).reset_index()
    out = out.sort_values(["n_peaks", "max_signal", "max_score", "gene_id"], ascending=[False, False, False, True], na_position="last")
    return out[columns].reset_index(drop=True)
