"""
Utility functions and callbacks for the AEGIS CLI suite.
"""
from __future__ import annotations
from pathlib import Path
from typing import TYPE_CHECKING, List, Sequence, Union, Optional

import typer
from typing_extensions import Annotated

from ..utils.genefunctions import NCBI_GENETIC_CODES, TAXONOMY_ORGANELLE_CODES

if TYPE_CHECKING:
    from ..annotation import Annotation


def split_callback(value: Union[str, Sequence[str], None]) -> List[str]:
    """
    Callback for Typer/Click options that accept comma-separated strings
    (e.g., '--features gene,CDS') or repeated flags (e.g., '-f gene -f CDS').

    Always returns a clean list of trimmed strings.
    """
    if not value:
        return []
    if isinstance(value, (list, tuple, set)):
        result: List[str] = []
        for item in value:
            if isinstance(item, str):
                result.extend([x.strip() for x in item.split(",") if x.strip()])
            elif item is not None:
                result.append(str(item).strip())
        return result
    if isinstance(value, str):
        return [item.strip() for item in value.split(",") if item.strip()]
    return [str(value).strip()]


# ---------------------------------------------------------------------------
# Standard Rich help panel titles used across the AEGIS CLI suite
# ---------------------------------------------------------------------------
EXTRACTION_PANEL = "Feature Extraction Options"
FEATURE_PANEL = "Feature & Biotype Filtering"
FILTER_PANEL = "Filtering Options"
MODEL_PANEL = "Structural Model Sanitisation"
ATTR_PANEL = "Attribute & Metadata Formatting"
CDS_PANEL = "CDS Inference & Reworking"
COORDS_PANEL = "Coordinate & Phase Options"
FORMAT_PANEL = "Format Conversion Options"
RENAME_PANEL = "Feature ID Renaming"
SUBFEATURE_PANEL = "Subfeature & Model Options"
SPLIT_PANEL = "Split Criteria"
SUBSET_PANEL = "Subset Criteria"
SUMMARY_PANEL = "Summary & Comparison Options"
OVERLAP_PANEL = "Overlap Criteria"
MERGE_PANEL = "Overlap & Merging Criteria"
PROMOTER_PANEL = "Promoter & Motif Options"
PEAKS_PANEL = "Peak Assignment Options"
SYMBOLS_PANEL = "Gene Symbol Options"
PRUNE_PANEL = "Pruning Options"
ORTHOLOGY_PANEL = "Orthology Tool Options"
BLAST_PANEL = "BLASTp Options"
GENOME_CLEANING_PANEL = "Genome Cleaning & Renaming"

FASTA_HEADER_PANEL = "Reference FASTA Options"
OUTPUT_HEADER_PANEL = "Output Sequence Header Options"
OUTPUT_FILTER_PANEL = "Result Filtering & Table Formatting"
GENETIC_CODES_PANEL = "Genetic Codes"
ORGANELLE_PANEL = "Organelle Detection"
IO_PANEL = "Input / Output Options"
COLUMNS_PANEL = "Output Columns"
EXEC_PANEL = "Execution & Debugging"


ANNOT_EXTENSIONS = (".gff3.gz", ".gff.gz", ".gtf.gz", ".gff3", ".gff", ".gtf")
FASTA_EXTENSIONS = (".fasta.gz", ".fa.gz", ".fna.gz", ".fas.gz", ".faa.gz", ".fasta", ".fa", ".fna", ".fas", ".faa")
BIO_EXTENSIONS = ANNOT_EXTENSIONS + FASTA_EXTENSIONS


def strip_bio_extension(name: str, extensions: tuple[str, ...] = BIO_EXTENSIONS) -> tuple[str, str]:
    """
    Split a filename or prefix into (stem, ext) only when it ends with a known
    bioinformatics extension (longest match first), preserving dotted prefixes
    such as 'release_v1.2'.
    """
    lower = name.lower()
    for ext in sorted(extensions, key=len, reverse=True):
        if lower.endswith(ext):
            return name[:-len(ext)], name[-len(ext):]
    return name, ""


DEFAULT_NAME_SENTINELS = (
    "{annotation-file}",
    "{genome-file}",
    "{annotation-filename(s)}",
    "{genome-filename(s)}",
)


def derive_bio_name(filepath: str, extensions: tuple[str, ...] = BIO_EXTENSIONS) -> str:
    """
    Derive a clean tag/name from an annotation or genome file path by stripping
    any known bioinformatics extension (including gzipped extensions such as
    '.gff3.gz' or '.fasta.gz'), preserving dotted stems like 'release_v1.2'.
    """
    stem, _ = strip_bio_extension(Path(filepath).name, extensions)
    return stem


def resolve_bio_name(
    filepath: str,
    provided_name: Optional[str] = None,
    *,
    extensions: tuple[str, ...] = BIO_EXTENSIONS,
) -> str:
    """
    Return a single annotation or genome name: the explicit name when provided
    (and not a default placeholder), otherwise derived from ``filepath``.
    """
    cleaned = provided_name.strip() if provided_name else ""
    if cleaned and cleaned not in DEFAULT_NAME_SENTINELS:
        return cleaned
    return derive_bio_name(filepath, extensions)


def _disambiguate_names(
    files: Sequence[str],
    raw_names: Sequence[str],
) -> tuple[List[str], List[tuple[str, str, str]]]:
    """
    Ensure every entry in ``raw_names`` is unique by appending ``_2``, ``_3``, ...
    to subsequent occurrences while avoiding collisions with already-used names
    and with raw names belonging to other entries.
    """
    reserved = set(raw_names)
    used: set[str] = set()
    next_suffix: dict[str, int] = {}
    resolved: List[str] = []
    disambiguated: List[tuple[str, str, str]] = []

    for filepath, base in zip(files, raw_names):
        if base not in used:
            final = base
        else:
            suffix = next_suffix.get(base, 2)
            candidate = f"{base}_{suffix}"
            while candidate in used or candidate in reserved:
                suffix += 1
                candidate = f"{base}_{suffix}"
            next_suffix[base] = suffix + 1
            final = candidate
            disambiguated.append((base, final, filepath))
        used.add(final)
        resolved.append(final)

    return resolved, disambiguated


def resolve_unique_names(
    files: Sequence[str],
    provided_names: Union[str, Sequence[str], None] = None,
    *,
    extensions: tuple[str, ...] = BIO_EXTENSIONS,
    kind: str = "input",
    share_by_file: bool = False,
    strict_count: bool = False,
    quiet: bool = False,
) -> List[str]:
    """
    Resolve and disambiguate names for a sequence of annotation or genome files.

    - Missing entries (or default placeholder sentinels) are derived from each
      filename via ``derive_bio_name``.
    - Duplicate names are disambiguated consistently with ``_2``, ``_3``, ...
      (skipping any name already used or claimed by another entry).
    - When ``share_by_file=True`` (e.g. genomes in ``aegis orthology``), entries
      pointing to the same resolved file path share a single name (raising
      ``typer.BadParameter`` if conflicting explicit names were assigned to the
      same file), while distinct files with colliding names are disambiguated.
    """
    if not files:
        return []

    if isinstance(provided_names, str):
        raw_provided = [n.strip() for n in provided_names.split(",") if n.strip()]
    else:
        raw_provided = [n.strip() for n in provided_names] if provided_names else []

    if raw_provided and all(n in DEFAULT_NAME_SENTINELS or not n for n in raw_provided):
        raw_provided = []

    if raw_provided:
        if strict_count and len(raw_provided) != len(files):
            raise typer.BadParameter(
                f"The provided number of {kind} name(s)/tag(s) ({len(raw_provided)}) "
                f"does not match the number of {kind} file(s) ({len(files)})."
            )
        if len(raw_provided) > len(files):
            raise typer.BadParameter(
                f"Too many {kind} name(s)/tag(s) provided ({len(raw_provided)}) "
                f"for {len(files)} {kind} file(s)."
            )

    raw_names: List[str] = []
    for idx, filepath in enumerate(files):
        if idx < len(raw_provided) and raw_provided[idx] and raw_provided[idx] not in DEFAULT_NAME_SENTINELS:
            raw_names.append(raw_provided[idx])
        else:
            raw_names.append(derive_bio_name(filepath, extensions))

    if share_by_file:
        resolved_paths = [str(Path(f).resolve()) for f in files]
        file_to_names: dict[str, set[str]] = {}
        for rpath, name in zip(resolved_paths, raw_names):
            file_to_names.setdefault(rpath, set()).add(name)
        for rpath, names_set in file_to_names.items():
            if len(names_set) > 1:
                raise typer.BadParameter(
                    f"The {kind} file '{Path(rpath).name}' is mapped to multiple different tags: {sorted(names_set)}. "
                    f"The same {kind} file must use a single consistent tag."
                )

        unique_orig_files: List[str] = []
        unique_rpaths: List[str] = []
        unique_raw_names: List[str] = []
        seen_rpaths: set[str] = set()
        for orig_f, rpath, name in zip(files, resolved_paths, raw_names):
            if rpath not in seen_rpaths:
                seen_rpaths.add(rpath)
                unique_orig_files.append(orig_f)
                unique_rpaths.append(rpath)
                unique_raw_names.append(name)

        unique_resolved, disambiguated = _disambiguate_names(unique_orig_files, unique_raw_names)
        rpath_to_final = dict(zip(unique_rpaths, unique_resolved))
        resolved = [rpath_to_final[rpath] for rpath in resolved_paths]
    else:
        resolved, disambiguated = _disambiguate_names(files, raw_names)

    if disambiguated and not quiet:
        for orig, final, filepath in disambiguated:
            typer.secho(
                f"Warning: Duplicate {kind} name '{orig}' disambiguated to '{final}' (for '{filepath}').",
                fg=typer.colors.YELLOW,
                err=True,
            )

    return resolved


def detect_file_type(filepath: str) -> str:
    """Detect whether a file is FASTA or GFF/GTF annotation."""
    import os
    if not filepath or not os.path.exists(filepath):
        return "unknown"
    fasta_exts = (".fa", ".fasta", ".fna", ".fas", ".fa.gz", ".fasta.gz", ".fna.gz")
    annot_exts = (".gff", ".gff3", ".gtf", ".gff.gz", ".gff3.gz", ".gtf.gz")
    lower = filepath.lower()
    if lower.endswith(fasta_exts):
        return "fasta"
    if lower.endswith(annot_exts):
        return "annotation"
    try:
        from ..utils.misc import open_file
        with open_file(filepath, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    return "fasta"
                if line.startswith("##") or "\t" in line:
                    return "annotation"
                break
    except Exception:
        pass
    return "unknown"


def apply_organelle_chroms(annotations: Sequence[Annotation], mitochondria_chroms: Sequence[str], chloroplast_chroms: Sequence[str]) -> None:
    """
    Gives each loaded annotation the listed organelle chromosomes it has, for commands comparing
    annotations with different contigs. Every listed name must be in at least one annotation or genome.
    """
    known_per_annotation = [set(a.chrs) | (set(a.genome.scaffolds) if a.genome is not None else set()) for a in annotations]
    known = set().union(*known_per_annotation)
    missing = [c for c in (*mitochondria_chroms, *chloroplast_chroms) if c not in known]
    if missing:
        raise typer.BadParameter(f"Organelle chromosome(s) not found in any annotation or genome: {', '.join(missing)}.")
    for a, a_known in zip(annotations, known_per_annotation):
        codes = a.genetic_codes
        codes["mitochondria_chroms"] = tuple(c for c in mitochondria_chroms if c in a_known)
        codes["chloroplast_chroms"] = tuple(c for c in chloroplast_chroms if c in a_known)
        a.set_genetic_codes(**codes)


def warn_organelle_genetic_codes(annotation: Annotation, show_name: bool = False) -> None:
    """
    Warns (on stderr, also with --quiet) which genetic codes translate the annotation's organelle contigs
    (see Annotation.chromosome_compartment): the default taxonomy, plant, mistranslates the mitochondria
    of animals and fungi. show_name prefixes the annotation name, for commands loading several.
    """
    organelles: dict[str, list[str]] = {}
    for chrom in annotation.chrs:
        compartment = annotation.chromosome_compartment(chrom)
        if compartment != "nuclear":
            organelles.setdefault(compartment, []).append(chrom)
    parts = []
    for compartment, label, table in (("mitochondria", "mitochondrial", annotation.mito_table), ("chloroplast", "plastid", annotation.plastid_table)):
        chroms = organelles.get(compartment)
        if chroms:
            shown = ", ".join(chroms[:3]) + (", ..." if len(chroms) > 3 else "")
            name = NCBI_GENETIC_CODES[table]["name"] if table in NCBI_GENETIC_CODES else "custom"
            parts.append(f"proteins on {label} contigs ({shown}) are translated with NCBI genetic code {table} ({name})")
    if parts:
        prefix = f"{annotation.name}: " if show_name else ""
        typer.secho(
            f"Warning: {prefix}{', and '.join(parts)} (taxonomy '{annotation.taxonomy}'). "
            "If this does not fit your organism, set --taxonomy (plant, vertebrate, invertebrate, yeast) or --mito-code/--plastid-code.",
            fg=typer.colors.YELLOW, err=True,
        )


def taxonomy_callback(value: str) -> str:
    if value.lower() not in TAXONOMY_ORGANELLE_CODES:
        raise typer.BadParameter(f"Unknown taxonomy '{value}'. Choose from: {', '.join(TAXONOMY_ORGANELLE_CODES)}.")
    return value.lower()


def genetic_code_callback(value: Optional[int]) -> Optional[int]:
    if value is not None and value not in NCBI_GENETIC_CODES:
        raise typer.BadParameter(f"Unsupported NCBI genetic code {value}. Supported tables: {', '.join(str(t) for t in NCBI_GENETIC_CODES)}.")
    return value


TaxonomyOption = Annotated[str, typer.Option(
    "-tax", "--taxonomy", callback=taxonomy_callback, rich_help_panel=GENETIC_CODES_PANEL,
    help="Taxonomic group preset for organelle genetic codes: 'plant' (default: nuclear=1, mito=1, plastid=11), 'vertebrate' (nuclear=1, mito=2), 'invertebrate' (nuclear=1, mito=5), or 'yeast' (nuclear=1, mito=3). Specific codes can be individually customized with --genetic-code, --mito-code, or --plastid-code."
)]
GeneticCodeOption = Annotated[int, typer.Option(
    "-gc", "--genetic-code", callback=genetic_code_callback, rich_help_panel=GENETIC_CODES_PANEL,
    help="NCBI genetic code table number for nuclear translation (default: 1)."
)]
MitoCodeOption = Annotated[Optional[int], typer.Option(
    "-mc", "--mito-code", callback=genetic_code_callback, rich_help_panel=GENETIC_CODES_PANEL,
    help="NCBI genetic code table number for mitochondrial contigs (overrides --taxonomy default)."
)]
PlastidCodeOption = Annotated[Optional[int], typer.Option(
    "-pc", "--plastid-code", callback=genetic_code_callback, rich_help_panel=GENETIC_CODES_PANEL,
    help="NCBI genetic code table number for plastid/chloroplast contigs (overrides --taxonomy default)."
)]
AutoOrganelleDetectionOption = Annotated[bool, typer.Option(
    "--auto-organelle-detection/--no-auto-organelle-detection", rich_help_panel=ORGANELLE_PANEL,
    help="Detect mitochondrial and plastid contigs from their name (e.g. 'chrM', 'Pt', 'ChrC') or, with a genome, their FASTA description. Organelle contigs are translated with --mito-code/--plastid-code, their genes are renamed as chrMT/chrCP with 'rename --prefix', and they are reported apart from nuclear chromosomes. --mitochondria-chroms/--chloroplast-chroms always apply (default: True)."
)]
MitochondriaChromsOption = Annotated[List[str], typer.Option(
    "--mitochondria-chroms", callback=split_callback, rich_help_panel=ORGANELLE_PANEL,
    help="Explicit list or comma-separated names of mitochondrial chromosomes/scaffolds, recognised as such in addition to (or without) automatic detection."
)]
ChloroplastChromsOption = Annotated[List[str], typer.Option(
    "--chloroplast-chroms", callback=split_callback, rich_help_panel=ORGANELLE_PANEL,
    help="Explicit list or comma-separated names of chloroplast/plastid chromosomes/scaffolds, recognised as such in addition to (or without) automatic detection."
)]
InitiatorMethionineOption = Annotated[bool, typer.Option(
    "--initiator-methionine/--no-initiator-methionine", rich_help_panel=GENETIC_CODES_PANEL,
    help="Translate valid start codons at position 1 to Methionine (M) [default: enabled; matches NCBI RefSeq/UniProt standard]. Use --no-initiator-methionine for literal codon translation."
)]
