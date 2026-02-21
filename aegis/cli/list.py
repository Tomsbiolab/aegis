import typer
import os

from typing_extensions import Annotated

from ..annotation import Annotation

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def genes(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag. [default: derived from filename]"
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder."
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output filename, without extension."
    )] = "{annotation-name}_genes_list.tsv",
    lengths: Annotated[bool, typer.Option(
        "-l", "--lengths", help="Include feature lengths in the output."
    )] = False,
    coordinates: Annotated[bool, typer.Option(
        "-c", "--coordinates", help="Include feature coordinates in the output."
    )] = False,
    chromosomes: Annotated[bool, typer.Option(
        "-ch", "--chromosomes", help="Include chromosome/scaffold information in the output."
    )] = False,
    coding_info: Annotated[bool, typer.Option(
        "--coding-info", help="Include whether a gene codes for a protein or not."
    )] = False,
    skip_coding: Annotated[bool, typer.Option(
        "--skip-coding", help="Whether to skip coding genes."
    )] = False,
    skip_non_coding: Annotated[bool, typer.Option(
        "--skip-non-coding", help="Whether to skip non-coding genes."
    )] = False,
    gene_symbols: Annotated[bool, typer.Option(
        "--gene-symbols", help="Whether to include gene symbols in the output."
    )] = False,
    sep: Annotated[str, typer.Option(
        "-s", "--sep", help="Separator for the output file. Tab is default."
    )] = "\t",
    skip_pseudogenes: Annotated[bool, typer.Option(
        "--skip-pseudogenes", help="Whether to skip pseudogenes."
    )] = False,
    skip_transposables: Annotated[bool, typer.Option(
        "--skip-transposables", help="Whether to skip transposable elements."
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
):
    """
    Lists genes from an annotation file and exports them to a TSV/CSV file.
    """
    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]
    
    if output_file == "{annotation-name}_genes_list.tsv":
        output_file = f"{annotation_name}_genes_list.tsv"

    os.makedirs(output_dir, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet)

    annotation.list_genes(
        custom_path=output_dir,
        output_file=output_file,
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

@app.command()
def transcripts(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag. [default: derived from filename]"
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder."
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output filename, without extension."
    )] = "{annotation-name}_transcripts_list.tsv",
    lengths: Annotated[bool, typer.Option(
        "-l", "--lengths", help="Include feature lengths in the output."
    )] = False,
    coordinates: Annotated[bool, typer.Option(
        "-c", "--coordinates", help="Include feature coordinates in the output."
    )] = False,
    chromosomes: Annotated[bool, typer.Option(
        "-ch", "--chromosomes", help="Include chromosome/scaffold information in the output."
    )] = False,
    coding_info: Annotated[bool, typer.Option(
        "--coding-info", help="Include whether a gene codes for a protein or not."
    )] = False,
    skip_coding: Annotated[bool, typer.Option(
        "--skip-coding", help="Whether to skip coding transcripts."
    )] = False,
    skip_non_coding: Annotated[bool, typer.Option(
        "--skip-non-coding", help="Whether to skip non-coding transcripts."
    )] = False,
    gene_symbols: Annotated[bool, typer.Option(
        "--gene-symbols", help="Whether to include gene symbols in the output."
    )] = False,
    sep: Annotated[str, typer.Option(
        "-s", "--sep", help="Separator for the output file. Default is tab."
    )] = "\t",
    skip_pseudogenes: Annotated[bool, typer.Option(
        "--skip-pseudogenes", help="Whether to skip pseudogenes."
    )] = False,
    skip_transposables: Annotated[bool, typer.Option(
        "--skip-transposables", help="Whether to skip transposable elements."
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
):
    """
    Lists transcripts from an annotation file and exports them to a TSV/CSV file.
    """
    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]
    
    if output_file == "{annotation-name}_transcripts_list.tsv":
        output_file = f"{annotation_name}_transcripts_list.tsv"

    os.makedirs(output_dir, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet)

    annotation.list_transcripts(
        custom_path=output_dir,
        output_file=output_file,
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
