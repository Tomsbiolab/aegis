import typer
import os

from typing_extensions import Annotated

from ..annotation import Annotation

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    symbols_file: Annotated[str, typer.Argument(
        help="Input tsv file with list of 'gene-id\tgene-symbol\n' entries. Excel '.xlsx' files are also permitted."
    )],
    sep: Annotated[str, typer.Option(
        "-s", "--separator", help="Indicate separator if input text file is separated by anything other than tabs."
    )] = "\t",
    header: Annotated[bool, typer.Option(
        "-H", "--header", help="Use this flag to indicate the presence of a header in input symbols_file."
    )] = False,
    clear_existing: Annotated[bool, typer.Option(
        "-c", "--clear-existing-symbols", help="Clears existing names and symbols from annotation file. Otherwise additional symbols are added to the existing ones."
    )] = False,
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output directory."
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation filename, without extension."
    )] = "{annotation-name}_symbols.gff3",
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,    

):
    """
    Add gene symbols to a gff (as 'symbol=' attributes) based on an input tsv with 'gene-id\tgene-symbol\n' entries.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_dir, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet)

    if output_file == "{annotation-name}_symbols.gff3":
        output_file = f"{annotation_name}_symbols.gff3"

    annotation.add_gene_symbols(clear=clear_existing, header=header, sep=sep, file_path=symbols_file)

    annotation.export.gff(custom_path=output_dir, tag=output_file, quiet=quiet)
    
if __name__ == "__main__":
    app()
