import typer
import os

from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    genome_file: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file."
    )],
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-g", "--genome-name", help="Genome assembly version, name or tag."
    )] = "{genome-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder."
    )] = "./aegis_output/",
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
):
    """
    Outputs a summary with descriptive stats.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    if genome_name == "{genome-file}":
        genome_name = os.path.splitext(os.path.basename(genome_file))[0]

    os.makedirs(output_dir, exist_ok=True)

    genome = Genome(name=genome_name, genome_file_path=genome_file, quiet=quiet)
    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, genome=genome, quiet=quiet)

    annotation.stats.update(custom_path=output_dir, export=True, quiet=quiet)


if __name__ == "__main__":
    app()
