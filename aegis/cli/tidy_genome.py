import typer
import os
from typing_extensions import Annotated
from aegis.annotation import Annotation
from aegis.genome import Genome

app = typer.Typer(add_completion=False)
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
        "-gn", "--genome-name", help="Genome assembly version, name or tag."
    )] = "{genome-fasta}",
    remove_scaffolds: Annotated[bool, typer.Option(
        "-rs", "--remove-scaffolds", help="Enable the removal of scaffolds and unplaced contigs from the genome."
    )] = False,
    remove_organelles: Annotated[bool, typer.Option(
        "-ro", "--remove-organelles", help="Remove mitochondrial and chloroplast chromosomes. Only effective if --remove-scaffolds is also enabled."
    )] = False,
    remove_chr00: Annotated[bool, typer.Option(
        "-r0", "--remove_chr00", help="Remove chromosomes named 'chr00' or similar, often representing unknown chromosomes. Only effective if --remove-scaffolds is also enabled."
    )] = False,
    rename_map_path: Annotated[str, typer.Option(
        "-rc", "--rename-map-path", help="Path to a TSV file for renaming chromosomes. Format: 'old_name<tab>new_name' per line, without a header."
    )] = None,
    output_folder: Annotated[str, typer.Option(
        "-o", "--output-folder", help="Path to the directory where output files will be saved."
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

    genome = Genome(name = genome_name, genome_file_path = genome_file)

    if annotation_file:
        
        annotation = Annotation(annotation_file, annotation_name, genome=genome)
    
    os.makedirs(output_folder, exist_ok=True)

    if rename_map_path is not None:
        with open(rename_map_path, encoding='utf-8') as f:
            chromosome_rename_map = {linea.split('\t')[0]: linea.split('\t')[1].strip() for linea in f}

        chromosome_equivalences = genome.rename_features_from_dic(rename_map=chromosome_rename_map)

    if remove_scaffolds:
        genome.remove_scaffolds(remove_00=remove_chr00, remove_organelles=remove_organelles)

    genome.export(output_folder = output_folder)

    if annotation_file:

        annotation.rename_chromosomes(equivalences=chromosome_equivalences)
        annotation.export_gff(custom_path=output_folder)



if __name__ == "__main__":
    app()
