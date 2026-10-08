import typer
import os
from typing import Optional, List
from typing_extensions import Annotated

from ..annotation import Annotation
from .utils import COLUMNS_PANEL, FILTER_PANEL, IO_PANEL, EXEC_PANEL, split_callback

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def genes(
    annotation_file: Annotated[Optional[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file (or provide via -a/--annotation)."
    )] = None,

    # 1. Output Columns
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

    # 2. Filtering Options
    coding_only: Annotated[bool, typer.Option(
        "--coding-only", "--skip-non-coding", help="Whether to include only coding genes.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    non_coding_only: Annotated[bool, typer.Option(
        "--non-coding-only", "--skip-coding", help="Whether to include only non-coding genes.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    biotypes: Annotated[List[str], typer.Option(
        "-b", "-r", "--biotypes", "--rna-classes", help="Filter by transcript biotype (e.g. 'mRNA,lncRNA'). Comma-separated list.",
        callback=split_callback,
        rich_help_panel=FILTER_PANEL,
    )] = [],
    skip_pseudogenes: Annotated[bool, typer.Option(
        "--skip-pseudogenes", help="Whether to skip pseudogenes.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    skip_te: Annotated[bool, typer.Option(
        "--skip-te", "--skip-transposables", help="Whether to skip transposable elements.",
        rich_help_panel=FILTER_PANEL,
    )] = False,

    # 3. Input / Output Options
    annotation_file_opt: Annotated[str, typer.Option(
        "-a", "--annotation", "--annotations", "--annotation-file", "--annot", help="Path to input annotation GFF/GTF file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", "--annotation-names", "--annot-name", help="Annotation version, name or tag. [default: derived from filename]",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output filename, with or without extension.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name}_genes_list.tsv",
    sep: Annotated[str, typer.Option(
        "-s", "--sep", "--separator", help="Separator for the output file. Tab is default.",
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
    Lists genes from an annotation file and exports them to a TSV/CSV file.
    """
    if verbose:
        quiet = False

    annot_in = annotation_file_opt if annotation_file_opt else annotation_file
    if not annot_in:
        raise typer.BadParameter("Missing required annotation file. Provide as positional argument or via -a/--annotation.")
    annotation_file = annot_in

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]
    
    if "{annotation-name}" in output_file:
        output_file = output_file.replace("{annotation-name}", annotation_name)

    os.makedirs(output_dir, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet, skip_coordinate_polishing=True)

    if biotypes:
        annotation.filter_by_rna_class(rna_classes=biotypes, remove_genes_accordingly=True, quiet=quiet)

    skip_coding = non_coding_only
    skip_non_coding = coding_only
    skip_transposables = skip_te

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
        main_transcript_length_instead_of_gene_length=transcript_length,
        quiet=quiet,
    )

@app.command()
def transcripts(
    annotation_file: Annotated[Optional[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file (or provide via -a/--annotation)."
    )] = None,

    # 1. Output Columns
    gene_id: Annotated[bool, typer.Option(
        "--gene-id", "--gene", help="Include parent gene ID in the output.",
        rich_help_panel=COLUMNS_PANEL,
    )] = False,
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

    # 2. Filtering Options
    coding_only: Annotated[bool, typer.Option(
        "--coding-only", "--skip-non-coding", help="Whether to include only coding transcripts.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    non_coding_only: Annotated[bool, typer.Option(
        "--non-coding-only", "--skip-coding", help="Whether to include only non-coding transcripts.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    main_only: Annotated[bool, typer.Option(
        "-m", "--main", "--only-main", "--main-only", help="Only list the primary / main transcript for each gene.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    biotypes: Annotated[List[str], typer.Option(
        "-b", "-r", "--biotypes", "--rna-classes", help="Filter by transcript biotype (e.g. 'mRNA,lncRNA'). Comma-separated list.",
        callback=split_callback,
        rich_help_panel=FILTER_PANEL,
    )] = [],
    skip_pseudogenes: Annotated[bool, typer.Option(
        "--skip-pseudogenes", help="Whether to skip pseudogenes.",
        rich_help_panel=FILTER_PANEL,
    )] = False,
    skip_te: Annotated[bool, typer.Option(
        "--skip-te", "--skip-transposables", help="Whether to skip transposable elements.",
        rich_help_panel=FILTER_PANEL,
    )] = False,

    # 3. Input / Output Options
    annotation_file_opt: Annotated[str, typer.Option(
        "-a", "--annotation", "--annotations", "--annotation-file", "--annot", help="Path to input annotation GFF/GTF file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", "--annotation-names", "--annot-name", help="Annotation version, name or tag. [default: derived from filename]",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output filename, with or without extension.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name}_transcripts_list.tsv",
    sep: Annotated[str, typer.Option(
        "-s", "--sep", "--separator", help="Separator for the output file. Default is tab.",
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
    Lists transcripts from an annotation file and exports them to a TSV/CSV file.
    """
    if verbose:
        quiet = False

    annot_in = annotation_file_opt if annotation_file_opt else annotation_file
    if not annot_in:
        raise typer.BadParameter("Missing required annotation file. Provide as positional argument or via -a/--annotation.")
    annotation_file = annot_in

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]
    
    if "{annotation-name}" in output_file:
        output_file = output_file.replace("{annotation-name}", annotation_name)

    os.makedirs(output_dir, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet, skip_coordinate_polishing=True)

    if biotypes:
        annotation.filter_by_rna_class(rna_classes=biotypes, remove_genes_accordingly=True, quiet=quiet)

    skip_coding = non_coding_only
    skip_non_coding = coding_only
    skip_transposables = skip_te

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
        gene_symbols=gene_symbols,
        only_main=main_only,
        gene_id=gene_id,
        quiet=quiet,
    )

if __name__ == "__main__":
    app()
