import typer
import os

from typing_extensions import Annotated

from ..annotation import Annotation
from .utils import COLUMNS_PANEL, FILTER_PANEL, IO_PANEL, EXEC_PANEL

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def genes(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    lengths: Annotated[bool, typer.Option(
        "-l", "--lengths", help="Include feature lengths in the output.",
        rich_help_panel=COLUMNS_PANEL,
    )] = False,
    transcript_length: Annotated[bool, typer.Option(
        "--transcript-length", help="Include main transcript lengths instead of gene lengths in the output. Useful for example when 'gene length' is required for count normalisation purposes, where actually the length of the main transcript is what matters.",
        rich_help_panel=COLUMNS_PANEL,
    )] = False,
    coordinates: Annotated[bool, typer.Option(
        "-c", "--coordinates", help="Include feature coordinates in the output.",
        rich_help_panel=COLUMNS_PANEL,
    )] = False,
    chromosomes: Annotated[bool, typer.Option(
        "-ch", "--chromosomes", help="Include chromosome/scaffold information in the output.",
        rich_help_panel=COLUMNS_PANEL,
    )] = False,
    coding_info: Annotated[bool, typer.Option(
        "--coding-info", help="Include whether a gene codes for a protein or not.",
        rich_help_panel=COLUMNS_PANEL,
    )] = False,
    gene_symbols: Annotated[bool, typer.Option(
        "--gene-symbols", help="Whether to include gene symbols in the output.",
        rich_help_panel=COLUMNS_PANEL,
    )] = False,

    skip_coding: Annotated[bool, typer.Option(
        "--skip-coding", "--non-coding-only", help="Whether to skip coding genes.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    skip_non_coding: Annotated[bool, typer.Option(
        "--skip-non-coding", "--coding-only", help="Whether to skip non-coding genes.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    skip_pseudogenes: Annotated[bool, typer.Option(
        "--skip-pseudogenes", help="Whether to skip pseudogenes.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    skip_transposables: Annotated[bool, typer.Option(
        "--skip-transposables", "--skip-te", help="Whether to skip transposable elements.",
        rich_help_panel=FILTER_PANEL,
    )] = False,

    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag. [default: derived from filename]",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output filename, without extension.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name}_genes_list.tsv",
    sep: Annotated[str, typer.Option(
        "-s", "--sep", help="Separator for the output file. Tab is default.",
        rich_help_panel=IO_PANEL,
    )] = "\t",

    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help="Increase terminal reporting verbosity.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
):
    """
    Lists genes from an annotation file and exports them to a TSV/CSV file.
    """
    if verbose:
        quiet = False

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]
    
    if output_file == "{annotation-name}_genes_list.tsv":
        output_file = f"{annotation_name}_genes_list.tsv"

    os.makedirs(output_dir, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet, skip_coordinate_polishing=True)

    annotation.export.gene_list(
        output_dir=output_dir,
        filename=output_file,
        lengths=lengths,
        coordinates=coordinates,
        chromosomes=chromosomes,
        coding_info=coding_info,
        skip_coding=skip_coding,
        skip_non_coding=skip_non_coding,
        sep=sep,
        skip_pseudogenes=skip_pseudogenes,
        skip_transposables=skip_transposables,
        gene_symbols=gene_symbols,
        main_transcript_length_instead_of_gene_length=transcript_length
    )

@app.command()
def transcripts(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    lengths: Annotated[bool, typer.Option(
        "-l", "--lengths", help="Include feature lengths in the output.",
        rich_help_panel=COLUMNS_PANEL,
    )] = False,
    coordinates: Annotated[bool, typer.Option(
        "-c", "--coordinates", help="Include feature coordinates in the output.",
        rich_help_panel=COLUMNS_PANEL,
    )] = False,
    chromosomes: Annotated[bool, typer.Option(
        "-ch", "--chromosomes", help="Include chromosome/scaffold information in the output.",
        rich_help_panel=COLUMNS_PANEL,
    )] = False,
    coding_info: Annotated[bool, typer.Option(
        "--coding-info", help="Include whether a gene codes for a protein or not.",
        rich_help_panel=COLUMNS_PANEL,
    )] = False,
    gene_symbols: Annotated[bool, typer.Option(
        "--gene-symbols", help="Whether to include gene symbols in the output.",
        rich_help_panel=COLUMNS_PANEL,
    )] = False,

    skip_coding: Annotated[bool, typer.Option(
        "--skip-coding", "--non-coding-only", help="Whether to skip coding transcripts.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    skip_non_coding: Annotated[bool, typer.Option(
        "--skip-non-coding", "--coding-only", help="Whether to skip non-coding transcripts.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    skip_pseudogenes: Annotated[bool, typer.Option(
        "--skip-pseudogenes", help="Whether to skip pseudogenes.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    skip_transposables: Annotated[bool, typer.Option(
        "--skip-transposables", "--skip-te", help="Whether to skip transposable elements.",
        rich_help_panel=FILTER_PANEL,
    )] = False,

    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag. [default: derived from filename]",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output filename, without extension.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name}_transcripts_list.tsv",
    sep: Annotated[str, typer.Option(
        "-s", "--sep", help="Separator for the output file. Default is tab.",
        rich_help_panel=IO_PANEL,
    )] = "\t",

    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help="Increase terminal reporting verbosity.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
):
    """
    Lists transcripts from an annotation file and exports them to a TSV/CSV file.
    """
    if verbose:
        quiet = False

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]
    
    if output_file == "{annotation-name}_transcripts_list.tsv":
        output_file = f"{annotation_name}_transcripts_list.tsv"

    os.makedirs(output_dir, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet, skip_coordinate_polishing=True)

    annotation.export.transcript_list(
        output_dir=output_dir,
        filename=output_file,
        lengths=lengths,
        coordinates=coordinates,
        chromosomes=chromosomes,
        coding_info=coding_info,
        skip_coding=skip_coding,
        skip_non_coding=skip_non_coding,
        sep=sep,
        skip_pseudogenes=skip_pseudogenes,
        skip_transposables=skip_transposables,
        gene_symbols=gene_symbols
    )

if __name__ == "__main__":
    app()
