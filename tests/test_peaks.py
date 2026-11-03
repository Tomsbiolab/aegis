"""Tests for peak-to-gene assignment (aegis.peaks and the `aegis peaks` command)."""

import pandas as pd
import pytest
from typer.testing import CliRunner

from aegis.annotation import Annotation
from aegis.cli.peaks import app
from aegis.peaks import (
    assign_peaks_to_genes,
    detect_peak_format,
    read_peaks,
    summarise_target_genes,
)

runner = CliRunner(env={"COLUMNS": "120"})

# chr1: geneA (+) 1001-2000, geneB (-) 5001-6000 overlapping geneC (+) 5501-7000, geneD (+) 20001-21000
# chr2: geneE (-) 101-500
GENES = [
    ("chr1", "geneA", "+", 1001, 2000),
    ("chr1", "geneB", "-", 5001, 6000),
    ("chr1", "geneC", "+", 5501, 7000),
    ("chr1", "geneD", "+", 20001, 21000),
    ("chr2", "geneE", "-", 101, 500),
]


def _gff_text():
    lines = ["##gff-version 3"]
    for ch, g, strand, start, end in GENES:
        lines.append(f"{ch}\ttest\tgene\t{start}\t{end}\t.\t{strand}\t.\tID={g}")
        lines.append(f"{ch}\ttest\tmRNA\t{start}\t{end}\t.\t{strand}\t.\tID={g}.t1;Parent={g}")
        lines.append(f"{ch}\ttest\texon\t{start}\t{end}\t.\t{strand}\t.\tID={g}.t1.exon1;Parent={g}.t1")
        lines.append(f"{ch}\ttest\tCDS\t{start}\t{end}\t.\t{strand}\t0\tID={g}.t1.cds1;Parent={g}.t1")
    return "\n".join(lines) + "\n"


@pytest.fixture
def gff_file(tmp_path):
    path = tmp_path / "toy.gff3"
    path.write_text(_gff_text())
    return str(path)


@pytest.fixture
def annotation(gff_file):
    return Annotation(gff_file, name="toy", quiet=True, coordinate_polishing="none")


def _bed(tmp_path, rows, name="peaks.bed"):
    path = tmp_path / name
    path.write_text("".join("\t".join(map(str, r)) + "\n" for r in rows))
    return str(path)


def _links(df):
    return {(r.peak_id, r.gene_id) for r in df.itertuples()}


# ---------------------------------------------------------------------------
# Reading peaks
# ---------------------------------------------------------------------------

def test_read_bed_converts_to_one_based(tmp_path):
    path = _bed(tmp_path, [("chr1", 2099, 2200, "p1", 50), ("chr1", 2099, 2200, "p1", 50)])
    peaks = read_peaks(path)
    assert detect_peak_format(path) == "bed"
    assert len(peaks) == 1  # duplicated rows (annotated BEDs) are collapsed
    p = peaks[0]
    assert (p.start, p.end, p.summit, p.score, p.peak_id) == (2100, 2200, 2150, 50.0, "p1")


def test_read_narrowpeak_uses_summit_offset(tmp_path):
    path = _bed(tmp_path, [("chr1", 999, 1200, "peakX", 300, ".", 12.5, 4.2, 3.1, 50)], "x.narrowPeak")
    assert detect_peak_format(path) == "narrowpeak"
    p = read_peaks(path)[0]
    assert (p.start, p.end, p.summit, p.signal, p.p_log10, p.q_log10) == (1000, 1200, 1050, 12.5, 4.2, 3.1)


def test_read_gem_narrowpeak_uses_event_position(tmp_path):
    path = _bed(tmp_path, [("chr00", 5064640, 5064841, "00:5064740", 1000, ".", 36399.8, 999.0, 999.0, 100)], "g.narrowPeak")
    p = read_peaks(path)[0]
    assert p.summit == 5064740


def test_read_gem_events_table(tmp_path):
    path = tmp_path / "TF.GEM_events.txt"
    path.write_text(
        "Position\tIP\tControl\tFold\tExpectd\tQ_-lg10\tP_-lg10\tP_poiss\tIPvsEMP\tNoise\tKmerGroup\tMotifId\tKG_score\tStrand\n"
        "00:5064740\t36399.8\t  430.9\t   84.5\t  291.5\t 999.00\t 999.00\t 999.00\t  -0.83\t   0.07\t--------\t-1\t0.00\t*\n"
        "01:4000\t100.0\t 1.0\t 10.0\t 1.0\t 5.00\t 6.00\t 6.00\t 0\t 0\t--------\t-1\t0.00\t*\n"
    )
    assert detect_peak_format(str(path)) == "gem"
    peaks = read_peaks(str(path))
    # same interval as GEM's own BED output (200 bp)
    assert (peaks[0].chrom, peaks[0].start, peaks[0].end, peaks[0].summit) == ("00", 5064641, 5064840, 5064740)
    assert (peaks[0].signal, peaks[0].fold, peaks[0].q_log10) == (36399.8, 84.5, 999.0)
    point = read_peaks(str(path), event_width=0)[1]
    assert (point.start, point.end, point.summit) == (4000, 4000, 4000)


def test_non_unique_names_get_peak_ids(tmp_path):
    path = _bed(tmp_path, [("chr1", 10, 20, "same"), ("chr1", 30, 40, "same")])
    assert [p.peak_id for p in read_peaks(path)] == ["peak_1", "peak_2"]


# ---------------------------------------------------------------------------
# Assignment
# ---------------------------------------------------------------------------

def test_nearest_gene(annotation, tmp_path):
    path = _bed(tmp_path, [
        ("chr1", 2099, 2200, "down_A"),      # 100 bp after geneA
        ("chr1", 5599, 5700, "in_BC"),       # inside overlapping geneB and geneC
        ("chr2", 49, 60, "down_E"),          # downstream of geneE (minus strand, so past its TES)
    ])
    df = assign_peaks_to_genes(annotation, path, reference="gene", quiet=True)
    assert _links(df) == {("down_A", "geneA"), ("in_BC", "geneB"), ("in_BC", "geneC"), ("down_E", "geneE")}

    a = df[df.peak_id == "down_A"].iloc[0]
    assert (a.distance_to_gene, a.distance_to_tss, a.summit_location, a.n_linked_genes) == (100, 1149, "downstream", 1)

    bc = df[df.peak_id == "in_BC"].set_index("gene_id")
    assert list(bc.index) == ["geneC", "geneB"]  # ties sorted by distance to TSS
    assert bc.loc["geneB", "distance_to_tss"] == 6000 - 5650  # minus strand: downstream of TSS is positive
    assert bc.loc["geneC", "distance_to_tss"] == 5650 - 5501
    assert set(bc.distance_to_gene) == {0}
    assert set(bc.summit_location) == {"genic"}

    e = df[df.peak_id == "down_E"].iloc[0]
    assert (e.gene_chrom, e.distance_to_gene, e.summit_location) == ("chr2", 41, "downstream")


def test_reference_tss_vs_gene(annotation, tmp_path):
    path = _bed(tmp_path, [("chr1", 3950, 4050, "p")])  # summit 4000
    by_gene = assign_peaks_to_genes(annotation, path, reference="gene", quiet=True)
    by_tss = assign_peaks_to_genes(annotation, path, reference="tss", quiet=True)
    assert list(by_gene.gene_id) == ["geneB"]  # body 951 bp away
    assert list(by_tss.gene_id) == ["geneC"]   # TSS 1501 bp away (geneB TSS is 2000 bp away)
    assert by_tss.iloc[0].distance_to_tss == -1501


def test_one_gene_per_peak_and_max_distance(annotation, tmp_path):
    path = _bed(tmp_path, [("chr1", 5599, 5700, "in_BC"), ("chr1", 12000, 12100, "far")])
    df = assign_peaks_to_genes(annotation, path, reference="gene", one_gene_per_peak=True, max_distance=5000, quiet=True)
    assert _links(df) == {("in_BC", "geneC")}
    df = assign_peaks_to_genes(annotation, path, reference="gene", one_gene_per_peak=True, max_distance=5000, keep_unassigned=True, quiet=True)
    far = df[df.peak_id == "far"].iloc[0]
    assert far.gene_id == "" and far.n_linked_genes == 0


def test_window_mode(annotation, tmp_path):
    path = _bed(tmp_path, [("chr1", 3950, 4050, "p"), ("chr1", 6950, 7050, "q")])
    df = assign_peaks_to_genes(annotation, path, mode="window", upstream=2000, downstream=1000, quiet=True)
    # p (summit 4000): only geneC's window [3501, 8000] holds it; geneB's (minus) window is [4001, 8000]
    # q (summit 7000): geneB's window [4001, 8000] and inside geneC
    assert _links(df) == {("p", "geneC"), ("q", "geneB"), ("q", "geneC")}


def test_filters(annotation, tmp_path):
    path = _bed(tmp_path, [("chr1", 5599, 5700, "in_BC")])
    annotation.chrs["chr1"]["geneB"].pseudogene = True
    df = assign_peaks_to_genes(annotation, path, reference="gene", skip_pseudogenes=True, quiet=True)
    assert _links(df) == {("in_BC", "geneC")}


def test_summarise_target_genes(annotation, tmp_path):
    path = _bed(tmp_path, [("chr1", 5599, 5700, "a", 10), ("chr1", 6500, 6600, "b", 30), ("chr1", 2099, 2200, "c", 5)])
    genes = summarise_target_genes(assign_peaks_to_genes(annotation, path, reference="gene", quiet=True))
    top = genes.iloc[0]
    assert (top.gene_id, top.n_peaks, top.peak_ids, top.max_score) == ("geneC", 2, "a,b", 30.0)
    assert set(genes.gene_id) == {"geneA", "geneB", "geneC"}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def test_cli_peaks(gff_file, tmp_path):
    peaks = _bed(tmp_path, [("chr1", 2099, 2200, "down_A"), ("chr1", 5599, 5700, "in_BC")], "TF1.bed")
    out = tmp_path / "out"
    # annotation given first: arguments are swapped back
    result = runner.invoke(app, [gff_file, peaks, "-d", str(out), "-q"])
    assert result.exit_code == 0, result.stdout
    links = pd.read_csv(out / "TF1_on_toy_peak_genes.tsv", sep="\t")
    genes = pd.read_csv(out / "TF1_on_toy_target_genes.tsv", sep="\t")
    # default reference is the TSS: the peak inside geneB and geneC goes to geneC (TSS 149 bp away vs 350 bp)
    assert set(links.gene_id) == {"geneA", "geneC"}
    assert set(genes.gene_id) == {"geneA", "geneC"}

    result = runner.invoke(app, [peaks, gff_file, "-r", "gene", "-d", str(out), "-o", "gene.tsv", "-go", "none", "-q"])
    assert result.exit_code == 0, result.stdout
    assert set(pd.read_csv(out / "gene.tsv", sep="\t").gene_id) == {"geneA", "geneB", "geneC"}


def test_cli_peaks_invalid_mode(gff_file, tmp_path):
    peaks = _bed(tmp_path, [("chr1", 2099, 2200, "p")])
    result = runner.invoke(app, [peaks, gff_file, "-m", "closest", "-d", str(tmp_path), "-q"])
    assert result.exit_code != 0
