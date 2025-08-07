import typer
import os
from typing_extensions import Annotated
from aegis.annotation import Annotation, read_file_with_fallback
from aegis.genome import Genome

app = typer.Typer(add_completion=False)
@app.command()
def main(
    genome_file: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file."
    )],
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-gn", "--genome-name", help="Genome assembly version, name or tag."
    )] = "{genome-fasta}",
    remove_scaffolds: Annotated[bool, typer.Option(
        "-rs", "--remove-scaffolds", help="Remove scaffolds."
    )] = False,
    remove_organelles: Annotated[bool, typer.Option(
        "-ro", "--remove-organelles", help="Remove mitochondrial and chloroplast chromosomes."
    )] = False,
    remove_chr00: Annotated[bool, typer.Option(
        "-r0", "--remove_chr00", help="Remove unknown chromosome or chr00."
    )] = False,
    rename_chromosomes: Annotated[str, typer.Option(
        "-rc", "--rename-chromosomes", help="Path to the table with correspondences between original chromosome and new name. Tab sepparated file without header."
    )] = None,
    output_folder: Annotated[str, typer.Option(
        "-o", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/"
):
    """
    Removes scaffolds, chloroplast and mitochondrial chromosomes from genome fasta file.
    """

    if genome_name == "{genome-fasta}" and genome_file != "":
        genome_name = os.path.splitext(os.path.basename(genome_file))[0]
    elif genome_file == "":
        genome_name = "genome"

    os.makedirs(output_folder, exist_ok=True)

    genome = Genome(name = genome_name, genome_file_path = genome_file)

    genome.rename_features_dap(chromosome_dict={'chr01_kk' : 'chr01_jeje'})

    if remove_scaffolds:
        genome.remove_scaffolds(remove_00=remove_chr00, remove_organelles=remove_organelles)

    genome.export(output_folder = output_folder)


if __name__ == "__main__":
    app()
