import typer
import os

from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome

app = typer.Typer(add_completion=False, no_args_is_help=True)
@app.command()
def main(
    genome_file: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file."
    )],
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file. If provided, it will be processed to match the cleaned genome."
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-g", "--genome-name", help="Genome assembly version, name or tag."
    )] = "{genome-fasta}",
    remove_scaffolds: Annotated[bool, typer.Option(
        "--remove-scaffolds", help="Enable the removal of scaffolds and unplaced contigs from the genome."
    )] = False,
    remove_organelles: Annotated[bool, typer.Option(
        "--remove-organelles", help="Remove mitochondrial and chloroplast chromosomes. Only effective if --remove-scaffolds is also enabled."
    )] = False,
    remove_chr00: Annotated[bool, typer.Option(
        "--remove-chr00", help="Remove chromosomes named 'chr00' or similar, often representing unknown chromosomes. Only effective if --remove-scaffolds is also enabled."
    )] = False,
    rename_map: Annotated[str, typer.Option(
        "--rename-map", help="Path to a TSV file for renaming chromosomes. Format: 'old_name<tab>new_name' per line, without a header."
    )] = "",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the directory where output files will be saved."
    )] = "./aegis_output/"
):
    """
    Processes and cleans a genome FASTA file and its corresponding annotation (GFF/GTF).

    This tool can perform several cleaning operations:
    - Rename chromosomes based on a provided map file.
    - Remove scaffolds, unplaced contigs, and organellar DNA.
    """

    if genome_name == "{genome-fasta}" and genome_file != "":
        genome_name = os.path.splitext(os.path.basename(genome_file))[0]
    elif genome_file == "":
        genome_name = "genome"

    if (annotation_name == "{annotation-file}") and (annotation_file != ""):
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    if output_dir == "./aegis_output/":
        subfolder = True
    else:
        subfolder = False
        
    genome = Genome(name = genome_name, genome_file_path = genome_file)

    if annotation_file:
        
        annotation = Annotation(annotation_file, annotation_name, genome=genome)
    
    os.makedirs(output_dir, exist_ok=True)

    if rename_map != "":
        with open(rename_map, encoding='utf-8') as f:
            chromosome_rename_map = {linea.split('\t')[0]: linea.split('\t')[1].strip() for linea in f}

        chromosome_equivalences = genome.rename_features_from_dic(rename_map=chromosome_rename_map)

    if remove_scaffolds:
        genome.remove_scaffolds(remove_00=remove_chr00, remove_organelles=remove_organelles)

    genome.export(output_folder = output_dir)

    if annotation_file and rename_map != "":
        annotation.rename_chromosomes(equivalences=chromosome_equivalences)
        annotation.export.gff(custom_path=output_dir, subfolder=subfolder)


if __name__ == "__main__":
    app()
