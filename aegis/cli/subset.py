import typer
import os
import random
import warnings

from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome

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
    genome_file: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file."
    )] = "",
    chr_cap: Annotated[int, typer.Option(
        "--chr-cap", help="Add a chromosome cap to generate an annotation gff (and assembly fasta) subset(s)."
    )] = 2,
    chosen_chromosomes: Annotated[str, typer.Option(
        "-c", "--chromosomes", help="Overrides --chr-cap. Only the chosen chromosomes/scaffolds will be in the resulting annotation gff (and assembly fasta) subset(s)",
    callback=split_callback
    )] = None,
    gene_cap: Annotated[int, typer.Option(
        "--gene-cap", help="Add a total gene number cap to reduce size of gff subset. The gene cap will affect scaffolds/chromosomes as uniformly as possible."
    )] = 3000,
    min_genes: Annotated[int, typer.Option(
        "--min-genes", help="Minimum total number of genes in the subset. Overrides --chr-cap if needed. Does not override --chosen-chromosomes if used."
    )] = 1500,
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-g", "--genome-name", help="Genome assembly version, name or tag."
    )] = "{genome-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder."
    )] = "./aegis_output/subsets/",
    output_annot_file: Annotated[str, typer.Option(
        "-oa", "--output-annot-file", help="Path to the output annotation filename, including extension."
    )] = "{annotation-name}_subset.gff3",
    output_genome_file: Annotated[str, typer.Option(
        "-og", "--output-genome-file", help="Path to the output genome filename, including extension."
    )] = "{genome-name}_subset.fasta",
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
):
    """
    Obtain subsets of an annotation file, random or directed. Ramdom subsets prioritise chromosomal features if available. A lite version of a gff file and its corresponding genome fasta file can be useful for debugging/trialing tools.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    if genome_name == "{genome-file}" and genome_file != "":
        genome_name = os.path.splitext(os.path.basename(genome_file))[0]
    elif genome_file == "":
        genome_name = "genome"

    os.makedirs(output_dir, exist_ok=True)

    if output_annot_file == "{annotation-name}_subset.gff3":
        output_annot_file = f"{annotation_name}_subset.gff3"

    if output_genome_file == "{genome-name}_subset.fasta":
        output_genome_file = f"{genome_name}_subset.fasta"

    if genome_file:
        g = Genome(genome_name, genome_file)
        a = Annotation(annotation_file, annotation_name, genome=g, quiet=quiet)
        common_chromosomes = set(a.chrs).intersection(set(g.scaffolds))
        common_actual_chromosomes = common_chromosomes - g.scaffold_names
        common_actual_chromosomes_minus_mt_chl = common_actual_chromosomes - g.accessory_chromosome_names

        if not common_chromosomes:
            raise ValueError(f"There are no common scaffolds/chromosomes between provided annotation and genome file.")
    else:
        a = Annotation(annotation_file, annotation_name, quiet=quiet)
        common_chromosomes = set(a.chrs)

    if not chosen_chromosomes:
        if chr_cap > len(common_chromosomes):
            chosen_chromosomes = common_chromosomes.copy()
            if genome_file:
                warnings.warn(f"Cap value {chr_cap} exceeds the number of available scaffolds/chrosomomes ({len(common_chromosomes)}) common to both genome and annotation files. The subset, in any case, will be based on common chromosomes/scaffolds.", category=UserWarning)
            else:
                warnings.warn(f"Cap value {chr_cap} exceeds the number of available scaffolds/chrosomomes ({len(common_chromosomes)}) in annotation file. The subset, in any case, will be based on common chromosomes/scaffolds.", category=UserWarning)
        elif chr_cap <= len(common_actual_chromosomes_minus_mt_chl):
            chosen_chromosomes = set(random.sample(list(common_actual_chromosomes_minus_mt_chl), chr_cap))
        elif chr_cap <= len(common_actual_chromosomes):
            chosen_chromosomes = set(random.sample(list(common_actual_chromosomes), chr_cap))
        else:
            chosen_chromosomes = set(random.sample(list(common_chromosomes), chr_cap))

        chosen_chromosomes = a.subset(chosen_features=chosen_chromosomes, gene_cap=gene_cap, common_chromosomes=common_chromosomes, min_genes=min_genes)
    
    else:
        # if chosen_chromosomes is selected min_genes parameter is ignored
        chosen_chromosomes = a.subset(chosen_features=chosen_chromosomes, gene_cap=gene_cap, common_chromosomes=common_chromosomes, min_genes=0)

    a.export.gff(custom_path=output_dir, tag=output_annot_file, subfolder=False, skip_atypical_fts=True)

    if genome_file:
        g.subset(chosen_features=chosen_chromosomes)
        g.export(output_folder=output_dir, file=output_genome_file)

if __name__ == "__main__":
    app()
