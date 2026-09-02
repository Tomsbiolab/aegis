import typer
import os
import random
import warnings
import re

from typing import Optional
from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome

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
    genome_file: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file."
    )] = "",
    chr_cap: Annotated[Optional[int], typer.Option(
        "--chr-cap", help="Add a chromosome cap to generate an annotation gff (and assembly fasta) subset(s). Set to 0 to disable."
    )] = 2,
    no_chr_cap: Annotated[bool, typer.Option(
        "--no-chr-cap", help="Do not enforce a chromosome cap. Selects all common chromosomes/scaffolds."
    )] = False,
    chosen_chromosomes: Annotated[str, typer.Option(
        "-c", "--chromosomes", help="Overrides --chr-cap. Only the chosen chromosomes/scaffolds will be in the resulting annotation gff (and assembly fasta) subset(s). Add them and separate them by commas e.g. --chromosomes 'chr1,chr3'.",
        callback=split_callback
    )] = "",
    gene_cap: Annotated[Optional[int], typer.Option(
        "--gene-cap", help="Add a total gene number cap to reduce size of gff subset. Set to 0 to disable. The gene cap will affect scaffolds/chromosomes as uniformly as possible."
    )] = 3000,
    no_gene_cap: Annotated[bool, typer.Option(
        "--no-gene-cap", help="Do not enforce a gene cap. Retains all genes present on the selected chromosomes/scaffolds."
    )] = False,
    min_genes: Annotated[Optional[int], typer.Option(
        "--min-genes", help="Minimum total number of genes in the subset. Overrides --chr-cap if needed. Does not override --chosen-chromosomes if used. Set to 0 to disable."
    )] = 1500,
    no_min_genes: Annotated[bool, typer.Option(
        "--no-min-genes", help="Do not enforce a minimum total number of genes in the subset."
    )] = False,
    seed: Annotated[Optional[int], typer.Option(
        "-s", "--seed", help="Random seed for reproducible subset sampling."
    )] = None,
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
    header_id_tag: Annotated[str, typer.Option(
        "--header-id-tag", help="Extract chromosome/scaffold ID from FASTA header description by tag name (e.g., 'OriSeqID')."
    )] = "",
    header_id_regex: Annotated[str, typer.Option(
        "--header-id-regex", help="Extract chromosome/scaffold ID from FASTA header description using a regex capture group (e.g., 'OriSeqID=(\\S+)')."
    )] = "",
    gwh: Annotated[bool, typer.Option(
        "--gwh", help="Preset for Genome Warehouse (GWH) FASTA files. Automatically extracts original sequence IDs from 'OriSeqID=...' in headers."
    )] = False,
):
    """
    Obtain subsets of an annotation file, random or directed. Ramdom subsets prioritise chromosomal features if available. A lite version of a gff file and its corresponding genome fasta file can be useful for debugging/trialing tools.
    """
    if seed is not None:
        random.seed(seed)

    if no_gene_cap or gene_cap is None or gene_cap <= 0:
        gene_cap = None

    if no_chr_cap or chr_cap is None or chr_cap <= 0:
        chr_cap = None

    if no_min_genes or min_genes is None or min_genes <= 0:
        min_genes = 0

    chosen_chromosomes_set = set(chosen_chromosomes)

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

    fasta_exts = (".fa", ".fasta", ".fna", ".fas", ".fa.gz", ".fasta.gz", ".fna.gz")
    annot_exts = (".gff", ".gff3", ".gtf", ".gff.gz", ".gff3.gz", ".gtf.gz")

    annot_looks_like_fasta = annotation_file.lower().endswith(fasta_exts)
    genome_looks_like_annot = genome_file and genome_file.lower().endswith(annot_exts)

    if annot_looks_like_fasta or genome_looks_like_annot:
        typer.secho("\nWARNING: It looks like the input files might be swapped based on their extensions.", fg=typer.colors.YELLOW)
        typer.echo(f"  - Provided annotation file: {annotation_file}")
        if genome_file:
            typer.echo(f"  - Provided genome file:     {genome_file}")
        
        # Prompt the user
        is_correct = typer.confirm(
            "\nAre you sure the first provided file is an annotation (GFF/GTF) file and the latter an assembly (FASTA) file? (If AEGIS just did not recognise the extension names, go ahead.)"
        )
        if not is_correct:
            typer.secho("Aborting. Please run the command again with the correct file order.", fg=typer.colors.RED)
            raise typer.Abort()

    if genome_file:
        g = Genome(
            name=genome_name,
            genome_file_path=genome_file,
            quiet=quiet,
            header_id_tag=header_id_tag if header_id_tag != "" else None,
            header_id_regex=header_id_regex if header_id_regex != "" else None,
            gwh=gwh,
        )
        a = Annotation(annotation_file, annotation_name, genome=g, quiet=quiet)
        common_chromosomes = set(a.chrs).intersection(set(g.scaffolds))
        common_actual_chromosomes = common_chromosomes - g.scaffold_names
        common_actual_chromosomes_minus_mt_chl = common_actual_chromosomes - g.accessory_chromosome_names

        if not common_chromosomes:
            desc_oriseq_ids = set()
            for s in g.scaffolds.values():
                d = getattr(s, "description", "")
                m = re.search(r"(?:^|[\s;,])OriSeqID=([^\s;,]+)", d)
                if m:
                    desc_oriseq_ids.add(m.group(1))
            if set(a.chrs).intersection(desc_oriseq_ids):
                raise ValueError(
                    f"There are no common scaffolds/chromosomes between provided annotation and genome file. "
                    f"However, annotation chromosomes match 'OriSeqID' in FASTA header descriptions (GWH format). "
                    f"Please run with '--gwh' or '--header-id-tag OriSeqID' to re-header the genome."
                )
            raise ValueError(f"There are no common scaffolds/chromosomes between provided annotation and genome file.")
    else:
        a = Annotation(annotation_file, annotation_name, quiet=quiet)
        common_chromosomes = set(a.chrs)

    if not chosen_chromosomes_set:
        if chr_cap is None or chr_cap >= len(common_chromosomes):
            chosen_chromosomes_set = common_chromosomes.copy()
            if chr_cap is not None and chr_cap > len(common_chromosomes):
                if genome_file:
                    warnings.warn(f"Cap value {chr_cap} exceeds the number of available scaffolds/chrosomomes ({len(chosen_chromosomes_set)}) common to both genome and annotation files. The subset, in any case, will be based on common chromosomes/scaffolds.", category=UserWarning)
                else:
                    warnings.warn(f"Cap value {chr_cap} exceeds the number of available scaffolds/chrosomomes ({len(chosen_chromosomes_set)}) in annotation file. The subset, in any case, will be based on common chromosomes/scaffolds.", category=UserWarning)
        elif genome_file:
            if chr_cap <= len(common_actual_chromosomes_minus_mt_chl):
                chosen_chromosomes_set = set(random.sample(list(common_actual_chromosomes_minus_mt_chl), chr_cap))
            elif chr_cap <= len(common_actual_chromosomes):
                chosen_chromosomes_set = set(random.sample(list(common_actual_chromosomes), chr_cap))
            else:
                chosen_chromosomes_set = set(random.sample(list(common_chromosomes), chr_cap))
        else:
            chosen_chromosomes_set = set(random.sample(list(common_chromosomes), chr_cap))

        chosen_chromosomes_set = a.subset(
            chosen_features=chosen_chromosomes_set,
            gene_cap=gene_cap,
            common_chromosomes=common_chromosomes,
            min_genes=min_genes,
            seed=seed,
            quiet=quiet,
        )
    
    else:
        # if chosen_chromosomes is selected min_genes parameter is ignored
        chosen_chromosomes_set = a.subset(
            chosen_features=chosen_chromosomes_set,
            gene_cap=gene_cap,
            common_chromosomes=common_chromosomes,
            min_genes=0,
            seed=seed,
            quiet=quiet,
        )

    a.export.gff(output_dir=output_dir, filename=output_annot_file, subfolder=False, skip_atypical_fts=True)

    if genome_file:
        g.subset(chosen_features=chosen_chromosomes_set, quiet=quiet)
        g.export(output_dir=output_dir, filename=output_genome_file, quiet=quiet)

if __name__ == "__main__":
    app()
