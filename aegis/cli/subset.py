import typer
import os
import random
import warnings
from typing_extensions import Annotated
from aegis.annotation import Annotation
from aegis.genome import Genome

def split_callback(value:str):
    if value:
        return set([item.strip() for item in value.split(",")])
    return set([])

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
        "-gc", "--gene-cap", help="Add a total gene number cap to reduce size of gff subset. The gene cap will affect scaffolds/chromosomes as uniformly as possible."
    )] = 3000,
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-gn", "--genome-name", help="Genome assembly version, name or tag."
    )] = "{genome-fasta}",
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/subsets/",
    output_annot_file: Annotated[str, typer.Option(
        "-o", "--output-annot-file", help="Path to the output annotation filename, including extension."
    )] = "{annotation-name}_subset.gff3",
    output_genome_file: Annotated[str, typer.Option(
        "-og", "--output-genome-file", help="Path to the output genome filename, including extension."
    )] = "{genome-name}_subset.fasta"
):
    """
    Obtain subsets of a gff, random or directed. A lite version of a gff file and its corresponding genome fasta file can be useful for debugging/trialing tools.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    if genome_name == "{genome-fasta}" and genome_fasta != "":
        genome_name = os.path.splitext(os.path.basename(genome_fasta))[0]
    elif genome_fasta == "":
        genome_name = "genome"

    os.makedirs(output_folder, exist_ok=True)

    if output_annot_file == "{annotation-name}_subset.gff3":
        output_annot_file = f"{annotation_name}_subset.gff3"

    if output_genome_file == "{genome-name}_subset.fasta":
        output_genome_file = f"{genome_name}_subset.fasta"

    if genome_fasta:
        g = Genome(genome_name, genome_fasta)
        a = Annotation(annotation_file, annotation_name, genome=g)

    else:
        a = Annotation(annotation_file, annotation_name)

    if not chosen_chromosomes:
        if genome_fasta:
            common_chromosomes = set(a.chrs).intersection(set(g.scaffolds))
        else:
            common_chromosomes = set(a.chrs)
        
        if chr_cap > len(common_chromosomes):
            chosen_chromosomes = common_chromosomes.copy()
            if genome_fasta:
                warnings.warn(f"Cap value {chr_cap} exceeds the number of available scaffolds/chrosomomes ({len(common_chromosomes)}) common to both genome and annotation files. The subset, in any case, will be based on common chromosomes/scaffolds.", category=UserWarning)
            else:
                warnings.warn(f"Cap value {chr_cap} exceeds the number of available scaffolds/chrosomomes ({len(common_chromosomes)}) in annotation file. The subset, in any case, will be based on common chromosomes/scaffolds.", category=UserWarning)
        else:
            chosen_chromosomes = set(random.sample(list(common_chromosomes), chr_cap))

    a.subset(chosen_features=chosen_chromosomes, gene_cap=gene_cap)
    a.export_gff(custom_path=output_folder, tag=output_annot_file, subfolder=False, skip_atypical_fts=True)

    if genome_fasta:
        g.subset(chosen_features=chosen_chromosomes)
        g.export(output_folder=output_folder, file=output_genome_file)

if __name__ == "__main__":
    app()
