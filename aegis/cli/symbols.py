import typer
import os
from typing import Optional
from typing_extensions import Annotated

from ..annotation import Annotation
from .utils import IO_PANEL, EXEC_PANEL, SYMBOLS_PANEL

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[Optional[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file (or provide via -a/--annotation)."
    )] = None,
    symbols_file: Annotated[Optional[str], typer.Argument(
        help="Input TSV/XLSX file with 'gene-id\tgene-symbol' mapping (or provide via --symbols-file)."
    )] = None,

    # 1. Gene Symbol Options
    clear_existing: Annotated[bool, typer.Option(
        "-c", "--clear-existing-symbols", "--clear", help="Clears existing names and symbols from annotation file. Otherwise additional symbols are appended.",
        rich_help_panel=SYMBOLS_PANEL,
    )] = False,
    header: Annotated[bool, typer.Option(
        "-H", "--header", help="Indicate the presence of a column header in the input symbols file.",
        rich_help_panel=SYMBOLS_PANEL,
    )] = False,
    sep: Annotated[str, typer.Option(
        "--sep", "--separator", help="Column delimiter for input symbols file (default: tab).",
        rich_help_panel=SYMBOLS_PANEL,
    )] = "\t",

    # 2. Input / Output Options
    annotation_file_opt: Annotated[str, typer.Option(
        "-a", "--annotation", "--annotation-file", help="Path to input annotation GFF/GTF file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    symbols_file_opt: Annotated[str, typer.Option(
        "--symbols-file", help="Path to input TSV/XLSX symbols mapping file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", help="Annotation version, name or tag.",
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
    Add gene symbols to an annotation GFF based on an input TSV/Excel mapping of gene IDs to symbols.
    """
    if verbose:
        quiet = False

    annot_in = annotation_file_opt if annotation_file_opt else annotation_file
    symbols_in = symbols_file_opt if symbols_file_opt else symbols_file

    if not annot_in:
        raise typer.BadParameter("Missing required annotation file. Provide as positional argument or via -a/--annotation.")
    if not symbols_in:
        raise typer.BadParameter("Missing required symbols mapping file. Provide as positional argument or via --symbols-file.")

    annotation_file = annot_in
    symbols_file = symbols_in

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_dir, exist_ok=True)

    if output_dir == "./aegis_output/":
        subfolder = True
    else:
        subfolder = False

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet, skip_coordinate_polishing=True)

    if "{annotation-name}" in output_file:
        output_file = output_file.replace("{annotation-name}", annotation_name)

    annotation.add_gene_symbols(clear=clear_existing, header=header, sep=sep, file_path=symbols_file)

    annotation.export.gff(output_dir=output_dir, filename=output_file, quiet=quiet, subfolder=subfolder, symbols=True)
    
if __name__ == "__main__":
    app()
