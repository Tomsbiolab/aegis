import typer
import os
import warnings
from typing import List, Optional
from typing_extensions import Annotated

from ..annotation import Annotation
from .utils import split_callback, IO_PANEL, EXEC_PANEL, SYMBOLS_PANEL

VALID_TARGETS = ["Name", "Alias", "Symbol"]
VALID_MATCH_ON = ["id", "name", "alias"]
VALID_OUTPUT_FORMATS = ["auto", "gff3", "gtf"]

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[Optional[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file (or provide via -a/--annotation)."
    )] = None,
    symbols_file: Annotated[Optional[str], typer.Argument(
        help="Input TSV/CSV/XLSX table of ids (first column) and symbols (second column); see --id-column and --to for other layouts (or provide via --symbols-file)."
    )] = None,

    # 1. Gene Symbol Options
    symbols_file_opt: Annotated[str, typer.Option(
        "-s", "--symbols-file", "--symbols", help="Path to input TSV/CSV/XLSX symbols table. Overrides positional argument if provided.",
        rich_help_panel=SYMBOLS_PANEL,
    )] = "",
    header: Annotated[bool, typer.Option(
        "-H", "--header", help="The symbols table has a header row (its column names can then be used in --id-column and --to).",
        rich_help_panel=SYMBOLS_PANEL,
    )] = False,
    sep: Annotated[Optional[str], typer.Option(
        "--sep", "--separator", help="Column separator of the symbols table (default: ',' for .csv files, tab otherwise; not used for .xlsx/.xls).",
        rich_help_panel=SYMBOLS_PANEL,
    )] = None,
    targets: Annotated[List[str], typer.Option(
        "-t", "--to", "--targets", help=f"Attribute(s) to place the symbols in, separated by commas. Choose from: {VALID_TARGETS}. Each can name the table column holding its values as 'Target:column' (a number from 1, or a header name with -H); default column 2. E.g. '--to Name:2,Alias:3'. 'Name' is what genome browsers such as JBrowse label and index; 'Symbol' is written as 'symbol=' unless the annotation already spells it otherwise.",
        callback=split_callback,
        rich_help_panel=SYMBOLS_PANEL,
    )] = ["Name"],
    id_column: Annotated[str, typer.Option(
        "--id-column", help="Table column with the ids: a number from 1, or a header name with -H.",
        rich_help_panel=SYMBOLS_PANEL,
    )] = "1",
    split: Annotated[Optional[str], typer.Option(
        "--split", help="Split each table cell into several values on this separator (e.g. ';' or '|'). By default a cell is one value, commas included.",
        rich_help_panel=SYMBOLS_PANEL,
    )] = None,
    match_on: Annotated[List[str], typer.Option(
        "--match-on", help=f"What the table ids are compared with to find genes, tried in order and separated by commas. Choose from: {VALID_MATCH_ON}. E.g. 'id,alias' also finds genes by their Alias, such as old ids kept by 'aegis rename'. An id that matches several genes by name or alias sets all of them (with a warning).",
        callback=split_callback,
        rich_help_panel=SYMBOLS_PANEL,
    )] = ["id"],
    append: Annotated[bool, typer.Option(
        "--append", help="Add the symbols to the values a gene already has in the target attribute(s) instead of replacing them. Transcripts then keep theirs (unless --transcript-names).",
        rich_help_panel=SYMBOLS_PANEL,
    )] = False,
    clear_unmapped: Annotated[bool, typer.Option(
        "-c", "--clear", help="Also empty the target attribute(s) of genes not in the symbols table, and of their transcripts. By default they are left as read.",
        rich_help_panel=SYMBOLS_PANEL,
    )] = False,
    transcript_names: Annotated[bool, typer.Option(
        "--transcript-names", help="Also set the target attribute(s) of transcripts: table ids may be transcript ids, which name that transcript; the other transcripts of a mapped gene get the gene's symbols. Without it, transcripts of mapped genes have the target attribute(s) emptied (as their gene's are replaced), so they do not keep a stale name.",
        rich_help_panel=SYMBOLS_PANEL,
    )] = False,
    on_missing: Annotated[str, typer.Option(
        "--on-missing", help="What to do with ids in the symbols table that are not in the annotation: 'warn', 'raise' (stop with an error, writing nothing) or 'ignore'.",
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
        "-o", "--output-file", help="Output annotation filename; '{ext}' is replaced by the output format.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name}_symbols.{ext}",
    output_format: Annotated[str, typer.Option(
        "--output-format", help=f"Output format. Choose from: {VALID_OUTPUT_FORMATS}; 'auto' writes the format of the input annotation.",
        rich_help_panel=IO_PANEL,
    )] = "auto",

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
    Add gene symbols to an annotation from a TSV/CSV/Excel table of gene (and, with --transcript-names, transcript) ids and symbols.
    """
    if verbose:
        quiet = False

    annot_in = annotation_file_opt if annotation_file_opt else annotation_file
    symbols_in = symbols_file_opt if symbols_file_opt else symbols_file

    if not annot_in:
        raise typer.BadParameter("Missing required annotation file. Provide as positional argument or via -a/--annotation.")
    if not symbols_in:
        raise typer.BadParameter("Missing required symbols table. Provide as positional argument or via --symbols-file.")

    invalid = [t for t in targets if t.partition(":")[0].strip().lower() not in {v.lower() for v in VALID_TARGETS}]
    if invalid or not targets:
        raise typer.BadParameter(f"Invalid --to {invalid}. Choose from: {VALID_TARGETS}, each optionally as 'Target:column'.")
    invalid = [m for m in match_on if m.lower() not in VALID_MATCH_ON]
    if invalid or not match_on:
        raise typer.BadParameter(f"Invalid --match-on {invalid}. Choose from: {VALID_MATCH_ON}.")
    if on_missing not in ("warn", "raise", "ignore"):
        raise typer.BadParameter(f"Invalid --on-missing '{on_missing}'. Choose from: warn, raise, ignore.")
    output_format = output_format.lower().strip()
    if output_format not in VALID_OUTPUT_FORMATS:
        raise typer.BadParameter(f"Invalid --output-format '{output_format}'. Choose from: {VALID_OUTPUT_FORMATS}.")

    # Read before the annotation, so table errors (e.g. a missing column) show up at once
    try:
        columns = Annotation.read_symbols_table(symbols_in, targets=targets, id_column=id_column, header=header, sep=sep, split=split)
    except ValueError as e:
        raise typer.BadParameter(str(e))

    annotation_file = annot_in

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_dir, exist_ok=True)
    subfolder = False

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet, coordinate_polishing="none", collapse_exons=False, collapse_CDSs=False)

    if output_format == "auto":
        output_format = getattr(annotation, "input_format", "gff3")

    output_file = output_file.replace("{annotation-name}", annotation_name).replace("{ext}", output_format)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", UserWarning)
        try:
            missing = annotation.set_symbol_columns(columns, mode=("append" if append else "replace"), clear_unmapped=clear_unmapped, on_missing=on_missing, transcript_names=transcript_names, match_on=match_on)
        except KeyError as e:
            raise typer.BadParameter(str(e).strip("'\""))
    for w in caught:
        if not quiet:
            print(f"Warning: {w.message}")

    if not quiet:
        placed = "; ".join(f"{', '.join(attributes)} from {len(mapping)} ids" for attributes, mapping in columns)
        print(f"Placed symbols in {placed} ({'appended' if append else 'replaced'}); {len(missing)} ids in {symbols_in} were not found.")

    if output_format == "gtf":
        annotation.export.gtf(output_dir=output_dir, filename=output_file, quiet=quiet, subfolder=subfolder)
    else:
        annotation.export.gff(output_dir=output_dir, filename=output_file, quiet=quiet, subfolder=subfolder)

if __name__ == "__main__":
    app()
