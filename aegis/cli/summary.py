import typer
import os
from typing_extensions import Annotated
from aegis.annotation import Annotation
from aegis.genome import Genome

app = typer.Typer(add_completion=False)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    genome_fasta: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file."
    )],
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-g", "--genome-name", help="Genome assembly version, name or tag."
    )] = "{genome-fasta}",
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/"
):
    """
    Outputs a summary with descriptive stats.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_folder, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file)
    genome = Genome(name=genome_name, genome_file_path=genome_fasta)

    annotation.update_stats(custom_path=output_folder, export=True, genome=genome)


if __name__ == "__main__":
    app()
