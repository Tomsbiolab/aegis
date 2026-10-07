import typer
import os

from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome
from .utils import IO_PANEL, EXEC_PANEL, FASTA_HEADER_PANEL, CDS_PANEL

GENOME_CLEANING_PANEL = "Genome Cleaning & Renaming"

app = typer.Typer(add_completion=False, no_args_is_help=True)
@app.command()
def main(
    genome_file: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file."
    )],
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file. If provided, it will be processed to match the cleaned genome."
    )] = "",
    remove_scaffolds: Annotated[bool, typer.Option(
        "--remove-scaffolds", help="Enable the removal of scaffolds and unplaced contigs from the genome.",
        rich_help_panel=GENOME_CLEANING_PANEL,
    )] = False,
    remove_organelles: Annotated[bool, typer.Option(
        "--remove-organelles", help="Remove mitochondrial and chloroplast chromosomes from the genome.",
        rich_help_panel=GENOME_CLEANING_PANEL,
    )] = False,
    remove_chr00: Annotated[bool, typer.Option(
        "--remove-chr00", help="Remove chromosomes named 'chr00' or similar, often representing unplaced/unknown chromosomes.",
        rich_help_panel=GENOME_CLEANING_PANEL,
    )] = False,
    rename_map: Annotated[str, typer.Option(
        "--rename-map", help="Path to a TSV file for renaming chromosomes. Format: 'old_name<tab>new_name' per line, without a header.",
        rich_help_panel=GENOME_CLEANING_PANEL,
    )] = "",

    min_codon_len: Annotated[int, typer.Option(
        "--min-codon-len", help="Minimum codon length required for predicted ORFs (default: 2).",
        rich_help_panel=CDS_PANEL,
    )] = 2,
    recalculate_phases: Annotated[bool, typer.Option(
        "--recalculate-phases", help="Recalculate CDS segment phases based on segment lengths and splicing leftover, preserving 5' initial phase for partial CDSs.",
        rich_help_panel=CDS_PANEL,
    )] = False,
    reset_phases_zero: Annotated[bool, typer.Option(
        "--reset-phases-zero", help="Reset initial CDS phase to 0 and recalculate all downstream segment phases.",
        rich_help_panel=CDS_PANEL,
    )] = False,
    skip_coordinate_polishing: Annotated[bool, typer.Option(
        "--skip-coordinate-polishing", help="Do not mutate feature coordinates when boundaries differ; log discrepancies as warnings instead.",
        rich_help_panel=CDS_PANEL,
    )] = False,

    header_id_tag: Annotated[str, typer.Option(
        "--header-id-tag", help="Extract chromosome/scaffold ID from FASTA header description by tag name (e.g., 'OriSeqID').",
        rich_help_panel=FASTA_HEADER_PANEL,
    )] = "",
    header_id_regex: Annotated[str, typer.Option(
        "--header-id-regex", help="Extract chromosome/scaffold ID from FASTA header description using a regex capture group (e.g., 'OriSeqID=(\\S+)').",
        rich_help_panel=FASTA_HEADER_PANEL,
    )] = "",
    gwh: Annotated[bool, typer.Option(
        "--gwh", help="Preset for Genome Warehouse (GWH) FASTA files. Automatically extracts original sequence IDs from 'OriSeqID=...' in headers.",
        rich_help_panel=FASTA_HEADER_PANEL,
    )] = False,
    keep_description: Annotated[bool, typer.Option(
        "--keep-description/--no-keep-description", help="Preserve full FASTA header descriptions in output genome file.",
        rich_help_panel=FASTA_HEADER_PANEL,
    )] = False,

    genome_name: Annotated[str, typer.Option(
        "-g", "-gn", "--genome-name", help="Genome assembly version, name or tag.",
        rich_help_panel=IO_PANEL,
    )] = "{genome-file}",
    annotation_name: Annotated[str, typer.Option(
        "-a", "-an", "--annotation-name", help="Annotation version, name or tag.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the directory where output files will be saved.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_genome_file: Annotated[str, typer.Option(
        "-og", "-o", "--output-genome-file", "--output-file", help="Path to the output genome filename, with or without extension.",
        rich_help_panel=IO_PANEL,
    )] = "{genome-name}_tidy.fasta",
    output_annot_file: Annotated[str, typer.Option(
        "-oa", "--output-annot-file", help="Path to the output annotation filename, with or without extension.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name}_tidy.gff3",

    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help="Enable detailed console output.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
):
    """
    Processes and cleans a genome FASTA file and its corresponding annotation (GFF/GTF).

    This tool can perform several cleaning operations:
    - Rename chromosomes based on a provided map file.
    - Remove scaffolds, unplaced contigs, and organellar DNA.
    """

    if genome_name == "{genome-file}":
        genome_name = os.path.splitext(os.path.basename(genome_file))[0]

    if (annotation_name == "{annotation-file}") and (annotation_file != ""):
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    if output_dir == "./aegis_output/":
        subfolder = True
    else:
        subfolder = False
        
    g = Genome(
        name=genome_name,
        genome_file_path=genome_file,
        header_id_tag=header_id_tag if header_id_tag != "" else None,
        header_id_regex=header_id_regex if header_id_regex != "" else None,
        gwh=gwh,
    )

    if annotation_file:
        
        a = Annotation(
            annot_file_path=annotation_file,
            name=annotation_name,
            genome=g,
            skip_coordinate_polishing=skip_coordinate_polishing,
            min_codon_len=min_codon_len,
            recalculate_phases=recalculate_phases,
            reset_phases_zero=reset_phases_zero,
        )
    
    os.makedirs(output_dir, exist_ok=True)

    if rename_map != "":
        with open(rename_map, encoding='utf-8') as f:
            chromosome_rename_map = {linea.split('\t')[0]: linea.split('\t')[1].strip() for linea in f}

        chromosome_equivalences = g.rename_features_from_dic(rename_map=chromosome_rename_map)

    if verbose:
        quiet = False

    if remove_scaffolds:
        g.remove_scaffolds(remove_00=remove_chr00)
    elif remove_chr00:
        new_scaffolds = {s_id: sc.copy() for s_id, sc in g.scaffolds.items() if not sc.unknown_chromosome}
        g.scaffolds = new_scaffolds
        g.update()

    if remove_organelles:
        g.remove_organelles()


    if output_genome_file == "{genome-name}_tidy.fasta":
        output_genome_file = f"{genome_name}_tidy.fasta"

    if output_annot_file == "{annotation-name}_tidy.gff3":
        output_annot_file = f"{annotation_name}_matching_genome_tidy.gff3"

    g.export(output_dir = output_dir, filename=output_genome_file, subfolder=subfolder, keep_description=keep_description)

    if annotation_file:
        if rename_map != "":
            a.rename_chromosomes(equivalences=chromosome_equivalences)
        extra_chrs = set(a.chrs.keys()) - set(g.scaffolds.keys())
        if extra_chrs:
            a.remove_chromosomes(extra_chrs)
        a.export.gff(output_dir=output_dir, subfolder=subfolder, filename=output_annot_file)


if __name__ == "__main__":
    app()
