import typer
import os
from typing_extensions import Annotated
from aegis.annotation import Annotation
from aegis.genome import Genome

def split_callback(value:str):
    if value:
        return [item.strip() for item in value.split(",")]
    return []

app = typer.Typer(add_completion=False)
@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    genome_fasta: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file."
    )] = "",
    chr_cap: Annotated[int, typer.Option(
        "-cc", "--chr-cap", help="Add a chromosome cap to generate an annotation gff (and assembly fasta) subset(s)."
    )] = 2,
    chosen_chromosomes: Annotated[str, typer.Option(
    "-c", "--chosen-chromosomes", help="Overrides -cc. Only the chosen chromosomes/scaffolds will be in the resulting annotation gff (and assembly fasta) subset(s)",
    callback=split_callback
    )] = None,
    gene_cap: Annotated[int, typer.Option(
        "-gc", "--gene-cap", help="Add a total gene number cap to reduce size of gff subset."
    )] = 3000,
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-gn", "--genome-name", help="Genome assembly version, name or tag."
    )] = "{genome-fasta}",
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/",
    output_annot_file: Annotated[str, typer.Option(
        "-o", "--output-annot-file", help="Path to the output annotation filename, without extension"
    )] = "{annotation-name}_...",
    output_genome_file: Annotated[str, typer.Option(
        "-o", "--output-genome-file", help="Path to the output genome filename, without extension."
    )] = "{genome-name}_..."
):
    """
    Obtain subsets of a gff, random or directed. A lite version of a gff file and its corresponding genome fasta file can be useful for debugging/trialing tools.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_folder, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file)

    if output_file == "{annotation-name}.{ext}":
        output_file = f"{annotation_name}"

    if genome_name == "{genome-fasta}" and genome_fasta != "":
        genome_name = os.path.splitext(os.path.basename(genome_fasta))[0]
    elif genome_fasta == "":
        genome_name = "genome"

    output_annot_file += ".gff3"
    output_genome_file += ".fasta"

    if genome_fasta:
        g = Genome(genome_name, genome_fasta)

    





    if output_file == "{annotation-name}_...":
        output_file = f"MODIFY_suffix logic"

    output_file += ".gff3"
    annotation.export_gff(custom_path=output_folder, tag=output_file)


if __name__ == "__main__":
    app()
