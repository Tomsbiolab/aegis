import typer
import os
import re
import pandas as pd
from typing import Optional
from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome
from .utils import (
    IO_PANEL,
    EXEC_PANEL,
    FASTA_HEADER_PANEL,
    PROMOTER_PANEL,
    detect_file_type,
)

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[Optional[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file (or provide via -a/--annotation)."
    )] = None,
    genome_file: Annotated[Optional[str], typer.Argument(
        help="Path to the input genome FASTA file (or provide via -g/--genome)."
    )] = None,
    genelist: Annotated[Optional[str], typer.Argument(
        help="Input TSV file with list of 'gene-id' entries (or provide via -l/--genelist). Excel '.xlsx' files are also permitted."
    )] = None,
    motif: Annotated[Optional[str], typer.Argument(
        help="DNA motif pattern (plain sequence e.g. 'TATAAA' or regular expression e.g. 'TATA[AT]A'). Can also be provided via -m/--motif."
    )] = None,
    motif_length: Annotated[Optional[int], typer.Argument(
        help="Actual length of motif in bp (optional; auto-deduced for plain sequence motifs)."
    )] = None,

    # 1. Promoter & Motif Options
    promoter_size: Annotated[int, typer.Option(
        "-ps", "--promoter-size", help="Size of the promoter region in base pairs (bp) upstream of TSS or ATG depending on --promoter-type (default: 2000).",
        rich_help_panel=PROMOTER_PANEL,
    )] = 2000,
    promoter_type: Annotated[str, typer.Option(
        "-p", "--promoter-type", help="Reference point for promoter extraction: 'standard' (upstream of TSS), 'upstream_ATG' (upstream of main CDS ATG codon), or 'standard_plus_up_to_ATG' (upstream of TSS plus 5' UTR up to start codon).",
        rich_help_panel=PROMOTER_PANEL,
    )] = "standard",
    motif_length_opt: Annotated[Optional[int], typer.Option(
        "-ml", "--motif-length", help="Actual span/length of the motif in bp. Automatically deduced for plain sequences (e.g. 'TATAAA' -> 6), but must be explicitly specified for regular expressions containing metacharacters (e.g. 'TATA[AT]A' has span 6 bp).",
        rich_help_panel=PROMOTER_PANEL,
    )] = None,

    # 2. Input / Output Options
    annotation_file_opt: Annotated[str, typer.Option(
        "-a", "--annotation", "--annotation-file", help="Path to input annotation GFF/GTF file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    genome_file_opt: Annotated[str, typer.Option(
        "-g", "--genome", "--genome-file", help="Path to input genome FASTA file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    genelist_opt: Annotated[str, typer.Option(
        "-l", "--genelist", "--genelist-file", help="Input TSV/XLSX file with list of 'gene-id' entries. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    motif_opt: Annotated[str, typer.Option(
        "-m", "--motif", help="DNA motif pattern (plain sequence or regular expression). Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", help="Annotation version, name or tag.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-gn", "--genome-name", help="Genome assembly version, name or tag.",
        rich_help_panel=IO_PANEL,
    )] = "{genome-file}",
    query_tag: Annotated[str, typer.Option(
        "--genelist-tag", help="Query gene list tag/name to improve output description.",
        rich_help_panel=IO_PANEL,
    )] = "query_genes",
    motif_tag: Annotated[str, typer.Option(
        "--motif-tag", help="Motif tag/name to improve output description, e.g. '{TF}_{motif_name}'.",
        rich_help_panel=IO_PANEL,
    )] = "query_motif",
    header: Annotated[bool, typer.Option(
        "-H", "--header", help="Indicate the presence of a column header in the input genelist file.",
        rich_help_panel=IO_PANEL,
    )] = False,
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output directory.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",

    # 3. Reference FASTA Options
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

    # 4. Execution / Debugging
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
    Scans a set of 'query' genes to locate all occurrences of a specified DNA motif within their upstream promoter regions.
    """
    if verbose:
        quiet = False

    annot_in = annotation_file_opt if annotation_file_opt else annotation_file
    genome_in = genome_file_opt if genome_file_opt else genome_file
    genelist_in = genelist_opt if genelist_opt else genelist
    motif_in = motif_opt if motif_opt else motif
    final_motif_len = motif_length_opt if motif_length_opt is not None else motif_length

    if not annot_in:
        raise typer.BadParameter("Missing required annotation file. Provide as positional argument or via -a/--annotation.")
    if not genome_in:
        raise typer.BadParameter("Missing required genome file. Provide as positional argument or via -g/--genome.")
    if not genelist_in:
        raise typer.BadParameter("Missing required gene list file. Provide as positional argument or via -l/--genelist.")
    if not motif_in:
        raise typer.BadParameter("Missing required motif pattern. Provide as positional argument or via -m/--motif.")

    # Swap if accidentally provided in reverse order
    if detect_file_type(annot_in) == "fasta" and detect_file_type(genome_in) == "annotation":
        annot_in, genome_in = genome_in, annot_in

    annotation_file = annot_in
    genome_file = genome_in
    genelist = genelist_in

    if final_motif_len is None:
        if re.search(r'[[\](){}+*?|\\]', motif_in):
            raise typer.BadParameter(
                f"Cannot auto-deduce motif length from regular expression '{motif_in}'. "
                f"Please specify the actual span/length of the motif explicitly using -ml/--motif-length (e.g. for 'TATA[AT]A', use -ml 6)."
            )
        else:
            final_motif_len = len(motif_in)

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    if genome_name == "{genome-file}":
        genome_name = os.path.splitext(os.path.basename(genome_file))[0]

    if output_dir == "./aegis_output/":
        subfolder = True
    else:
        subfolder = False

    if header:
        if genelist.endswith(".xlsx"):
            df = pd.read_excel(genelist, skiprows=1, dtype=str)
        else:
            df = pd.read_csv(genelist, skiprows=1, dtype=str)
    else:
        if genelist.endswith(".xlsx"):
            df = pd.read_excel(genelist, dtype=str)
        else:
            df = pd.read_csv(genelist, dtype=str)

    df = df.fillna("")
    genes = df.iloc[:, 0].tolist()
    genes = [gene for gene in genes if gene != ""]

    genome = Genome(
        name=genome_name,
        genome_file_path=genome_file,
        quiet=quiet,
        header_id_tag=header_id_tag if header_id_tag != "" else None,
        header_id_regex=header_id_regex if header_id_regex != "" else None,
        gwh=gwh,
    )
    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, genome=genome, quiet=quiet, skip_coordinate_polishing=True)

    annotation.generate_promoters(promoter_size=promoter_size, promoter_type=promoter_type)

    annotation.motifs.find_and_plot(
        query_genes=genes,
        motif=motif_in,
        motif_length=final_motif_len,
        glistname=query_tag,
        tf_motif_tag=motif_tag,
        output_dir=output_dir,
        subfolder=subfolder,
        quiet=quiet,
    )
    
if __name__ == "__main__":
    app()
