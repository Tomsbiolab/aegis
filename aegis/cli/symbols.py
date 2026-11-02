import typer
import os
import warnings
from typing import List, Optional
from typing_extensions import Annotated

from ..annotation import Annotation
from .utils import split_callback, IO_PANEL, EXEC_PANEL, SYMBOLS_PANEL

VALID_TARGETS = ["Name", "Alias", "Symbol"]

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[Optional[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file (or provide via -a/--annotation)."
    )] = None,
    symbols_file: Annotated[Optional[str], typer.Argument(
        help="Input TSV/XLSX file with 'gene-id\tgene-symbol' mapping; with --transcript-names ids may also be transcript ids (or provide via --symbols-file)."
    )] = None,

    # 1. Gene Symbol Options
    symbols_file_opt: Annotated[str, typer.Option(
        "-s", "--symbols-file", "--symbols", help="Path to input TSV/XLSX symbols mapping file. Overrides positional argument if provided.",
        rich_help_panel=SYMBOLS_PANEL,
    )] = "",
    header: Annotated[bool, typer.Option(
        "-H", "--header", help="Indicate the presence of a column header in the input symbols file.",
        rich_help_panel=SYMBOLS_PANEL,
    )] = False,
    sep: Annotated[str, typer.Option(
        "--sep", "--separator", help="Column delimiter for input symbols file (default: tab).",
        rich_help_panel=SYMBOLS_PANEL,
    )] = "\t",
    targets: Annotated[List[str], typer.Option(
        "-t", "--to", "--targets", help=f"Attribute(s) to place the symbols in, separated by commas. Choose from: {VALID_TARGETS}. 'Name' is what genome browsers such as JBrowse label and index; 'Symbol' is written as 'symbol=' unless the annotation already spells it otherwise.",
        callback=split_callback,
        rich_help_panel=SYMBOLS_PANEL,
    )] = ["Name"],
    append: Annotated[bool, typer.Option(
        "--append", help="Add the symbols to the values a gene already has in the target attribute(s) instead of replacing them. Transcripts then keep theirs (unless --transcript-names).",
        rich_help_panel=SYMBOLS_PANEL,
    )] = False,
    clear_unmapped: Annotated[bool, typer.Option(
        "-c", "--clear", help="Also empty the target attribute(s) of genes not in the symbols file, and of their transcripts. By default they are left as read.",
        rich_help_panel=SYMBOLS_PANEL,
    )] = False,
    transcript_names: Annotated[bool, typer.Option(
        "--transcript-names", help="Also set the target attribute(s) of transcripts: ids in the symbols file may be transcript ids, which name that transcript; the other transcripts of a mapped gene get the gene's symbols. Without it, transcripts of mapped genes have the target attribute(s) emptied (as their gene's are replaced), so they do not keep a stale name.",
        rich_help_panel=SYMBOLS_PANEL,
    )] = False,
    on_missing: Annotated[str, typer.Option(
        "--on-missing", help="What to do with ids in the symbols file that are not in the annotation: 'warn', 'raise' (stop with an error) or 'ignore'.",
        rich_help_panel=SYMBOLS_PANEL,
    )] = "warn",

    # 2. Input / Output Options
    annotation_file_opt: Annotated[str, typer.Option(
        "-a", "--annotation", "--annotations", "--annotation-file", help="Path to input annotation GFF/GTF file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", "--annotation-names", help="Annotation version, name or tag.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output directory.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation filename.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name}_symbols.gff3",

    # 3. Execution / Debugging
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help="Enable detailed console output.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
):
    """
    Add gene symbols to an annotation GFF based on an input TSV/Excel mapping of gene (and, with --transcript-names, transcript) IDs to symbols.
    """
    if verbose:
        quiet = False

    annot_in = annotation_file_opt if annotation_file_opt else annotation_file
    symbols_in = symbols_file_opt if symbols_file_opt else symbols_file

    if not annot_in:
        raise typer.BadParameter("Missing required annotation file. Provide as positional argument or via -a/--annotation.")
    if not symbols_in:
        raise typer.BadParameter("Missing required symbols mapping file. Provide as positional argument or via --symbols-file.")

    invalid = [t for t in targets if t.lower() not in {v.lower() for v in VALID_TARGETS}]
    if invalid or not targets:
        raise typer.BadParameter(f"Invalid --to {invalid}. Choose from: {VALID_TARGETS}.")
    if on_missing not in ("warn", "raise", "ignore"):
        raise typer.BadParameter(f"Invalid --on-missing '{on_missing}'. Choose from: warn, raise, ignore.")

    annotation_file = annot_in
    symbols_file = symbols_in

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_dir, exist_ok=True)
    subfolder = False

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet, coordinate_polishing="none", collapse_exons=False, collapse_CDSs=False)

    if "{annotation-name}" in output_file:
        output_file = output_file.replace("{annotation-name}", annotation_name)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", UserWarning)
        try:
            missing = annotation.add_gene_symbols(file_path=symbols_file, targets=targets, mode=("append" if append else "replace"), clear_unmapped=clear_unmapped, header=header, sep=sep, on_missing=on_missing, transcript_names=transcript_names)
        except KeyError as e:
            raise typer.BadParameter(str(e).strip("'\""))
    for w in caught:
        if not quiet:
            print(f"Warning: {w.message}")

    if not quiet:
        print(f"Placed symbols in {', '.join(targets)} ({'appended' if append else 'replaced'}); {len(missing)} ids in {symbols_file} were not found.")

    annotation.export.gff(output_dir=output_dir, filename=output_file, quiet=quiet, subfolder=subfolder)
    
if __name__ == "__main__":
    app()
