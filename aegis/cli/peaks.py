import os
import typer
from typing import Optional
from typing_extensions import Annotated

from ..annotation import Annotation
from ..peaks import assign_peaks_to_genes, summarise_target_genes, read_peaks, PEAK_FORMATS, ChromosomeMismatchError
from .utils import FILTER_PANEL, IO_PANEL, EXEC_PANEL, PEAKS_PANEL

_ANNOTATION_EXTS = (".gff", ".gff3", ".gtf", ".gff.gz", ".gff3.gz", ".gtf.gz")


def _is_annotation_path(path:str) -> bool:
    return path.lower().endswith(_ANNOTATION_EXTS)

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.command()
def main(
    peaks_file: Annotated[Optional[str], typer.Argument(
        help="Path to the peak file: GEM events table (*.GEM_events.txt), narrowPeak (GEM/GPS/MACS2) or BED (or provide via -p/--peaks)."
    )] = None,
    annotation_file: Annotated[Optional[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file (or provide via -a/--annotation)."
    )] = None,

    # 1. Peak assignment options
    mode: Annotated[str, typer.Option(
        "-m", "--mode", help="'nearest': link each peak to its closest gene(s). 'window': link each peak to every gene whose window, from --upstream bp before the TSS to --downstream bp after the TES, contains the peak summit.",
        rich_help_panel=PEAKS_PANEL,
    )] = "nearest",
    reference: Annotated[str, typer.Option(
        "-r", "--reference", help="Distance used to choose the nearest gene: 'tss' (peak summit to TSS, as ChIPpeakAnno's nearestLocation default and HOMER) or 'gene' (peak to gene body, 0 when they overlap, so a peak inside a gene is always linked to it).",
        rich_help_panel=PEAKS_PANEL,
    )] = "tss",
    max_distance: Annotated[Optional[int], typer.Option(
        "-md", "--max-distance", help="In 'nearest' mode, ignore genes farther than this many bp (measured as set by --reference). [default: no limit]",
        rich_help_panel=PEAKS_PANEL,
    )] = None,
    upstream: Annotated[int, typer.Option(
        "-up", "--upstream", help="In 'window' mode, bp upstream of the TSS included in the gene window.",
        rich_help_panel=PEAKS_PANEL,
    )] = 2000,
    downstream: Annotated[int, typer.Option(
        "-down", "--downstream", help="In 'window' mode, bp downstream of the TES included in the gene window.",
        rich_help_panel=PEAKS_PANEL,
    )] = 1000,
    one_gene_per_peak: Annotated[bool, typer.Option(
        "--one-gene-per-peak", help="In 'nearest' mode, keep a single gene per peak when several are equally close (e.g. overlapping genes), choosing the one with the closest TSS.",
        rich_help_panel=PEAKS_PANEL,
    )] = False,
    keep_unassigned: Annotated[bool, typer.Option(
        "--keep-unassigned", help="Also write peaks that were not linked to any gene, with empty gene columns.",
        rich_help_panel=PEAKS_PANEL,
    )] = False,
    peak_format: Annotated[str, typer.Option(
        "-f", "--peak-format", help="Peak file format: 'auto', 'gem' (GEM events table), 'narrowpeak' or 'bed'.",
        rich_help_panel=PEAKS_PANEL,
    )] = "auto",
    event_width: Annotated[int, typer.Option(
        "--event-width", help="Width (bp) given to events of a GEM events table, centred on the event position. GEM's own BED/narrowPeak outputs use 200 bp.",
        rich_help_panel=PEAKS_PANEL,
    )] = 200,

    # 2. Gene filtering
    coding_only: Annotated[bool, typer.Option(
        "--coding-only", "--skip-non-coding", help="Only consider protein-coding genes as targets.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    skip_pseudogenes: Annotated[bool, typer.Option(
        "--skip-pseudogenes", help="Do not consider pseudogenes as targets.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    skip_te: Annotated[bool, typer.Option(
        "--skip-te", "--skip-transposables", help="Do not consider transposable element genes as targets.",
        rich_help_panel=FILTER_PANEL,
    )] = False,

    # 3. Input / Output Options
    peaks_file_opt: Annotated[str, typer.Option(
        "-p", "--peaks", "--peaks-file", help="Path to the peak file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    annotation_file_opt: Annotated[str, typer.Option(
        "-a", "--annotation", "--annotation-file", "--annot", help="Path to input annotation GFF/GTF file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", "--annot-name", help="Annotation version, name or tag. [default: derived from filename]",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    peaks_name: Annotated[str, typer.Option(
        "-pn", "--peaks-name", help="Name or tag of the peak set (e.g. the TF). [default: derived from filename]",
        rich_help_panel=IO_PANEL,
    )] = "{peaks-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Output table with one row per peak-gene link.",
        rich_help_panel=IO_PANEL,
    )] = "{peaks-name}_on_{annotation-name}_peak_genes.tsv",
    genes_output_file: Annotated[str, typer.Option(
        "-go", "--genes-output-file", help="Output table with one row per target gene. Use 'none' to skip it.",
        rich_help_panel=IO_PANEL,
    )] = "{peaks-name}_on_{annotation-name}_target_genes.tsv",
    sep: Annotated[str, typer.Option(
        "-s", "--sep", "--separator", help="Separator for the output files. Tab is default.",
        rich_help_panel=IO_PANEL,
    )] = "\t",

    # 4. Execution & Debugging
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help="Increase terminal reporting verbosity.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
):
    """
    Links ChIP-seq / DAP-seq peaks (GEM output, narrowPeak or BED) to their closest genes, or to every gene within a
    window around them, and writes a peak-gene table plus a target gene table.
    """
    if verbose:
        quiet = False

    peaks_in = peaks_file_opt if peaks_file_opt else peaks_file
    annot_in = annotation_file_opt if annotation_file_opt else annotation_file
    if not peaks_in:
        raise typer.BadParameter("Missing required peak file. Provide as positional argument or via -p/--peaks.")
    if not annot_in:
        raise typer.BadParameter("Missing required annotation file. Provide as positional argument or via -a/--annotation.")

    # Swap if given in reverse order (peak files are tab-separated too, so rely on GFF/GTF extensions)
    if _is_annotation_path(peaks_in) and not _is_annotation_path(annot_in):
        peaks_in, annot_in = annot_in, peaks_in

    mode = mode.lower()
    reference = reference.lower()
    peak_format = peak_format.lower()
    if mode not in ("nearest", "window"):
        raise typer.BadParameter(f"Invalid --mode '{mode}'. Choose from: nearest, window.")
    if reference not in ("gene", "tss"):
        raise typer.BadParameter(f"Invalid --reference '{reference}'. Choose from: gene, tss.")
    if peak_format not in PEAK_FORMATS:
        raise typer.BadParameter(f"Invalid --peak-format '{peak_format}'. Choose from: {', '.join(PEAK_FORMATS)}.")
    if upstream < 0 or downstream < 0:
        raise typer.BadParameter("--upstream and --downstream must be >= 0.")
    if max_distance is not None and max_distance < 0:
        raise typer.BadParameter("--max-distance must be >= 0.")

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annot_in))[0]
    if peaks_name == "{peaks-file}":
        base = os.path.basename(peaks_in)
        if base.endswith(".gz"):
            base = base[:-3]
        for suffix in (".GEM_events.txt", ".GPS_events.txt", ".GEM_events.narrowPeak", ".GPS_events.narrowPeak", ".GEM_events.bed", ".GPS_events.bed"):
            if base.endswith(suffix):
                base = base[:-len(suffix)]
                break
        else:
            base = os.path.splitext(base)[0]
        peaks_name = base

    def resolve(name:str) -> str:
        return name.replace("{peaks-name}", peaks_name).replace("{annotation-name}", annotation_name)

    output_file = resolve(output_file)
    write_genes = genes_output_file.lower() != "none"
    genes_output_file = resolve(genes_output_file)

    try:
        peaks = read_peaks(peaks_in, peak_format=peak_format, event_width=event_width)
    except ValueError as e:
        raise typer.BadParameter(str(e))

    annotation = Annotation(name=annotation_name, annot_file_path=annot_in, quiet=quiet, coordinate_polishing="none", skip_atypical_features=True)

    try:
        links = assign_peaks_to_genes(
            annotation,
            peaks,
            mode=mode,
            reference=reference,
            max_distance=max_distance,
            upstream=upstream,
            downstream=downstream,
            one_gene_per_peak=one_gene_per_peak,
            keep_unassigned=keep_unassigned,
            coding_only=coding_only,
            skip_pseudogenes=skip_pseudogenes,
            skip_te=skip_te,
            quiet=quiet,
        )
    except ChromosomeMismatchError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)

    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, output_file)
    links.to_csv(out_path, sep=sep, index=False)
    if not quiet:
        print(f"Peak-gene links written to {out_path}")

    if write_genes:
        genes_path = os.path.join(output_dir, genes_output_file)
        summarise_target_genes(links).to_csv(genes_path, sep=sep, index=False)
        if not quiet:
            print(f"Target genes written to {genes_path}")


if __name__ == "__main__":
    app()
