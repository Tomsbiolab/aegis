import typer
import os
import pandas as pd

from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome
from .utils import IO_PANEL, EXEC_PANEL, FASTA_HEADER_PANEL

PROMOTER_PANEL = "Promoter & Motif Options"

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    genome_file: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file."
    )],
    genelist: Annotated[str, typer.Argument(
        help="Input tsv file with list of 'gene-id' entries. Excel '.xlsx' files are also permitted."
    )],
    motif: Annotated[str, typer.Argument(
        help="Regular expression pattern in python describing motif."
    )],
    motif_length: Annotated[int, typer.Argument(
        help="Actual length of motif."
    )],
    # 1. Promoter & Motif Options
    promoter_size: Annotated[int, typer.Option(
        "-ps", "--promoter-size", help="Size of the promoter region in base pairs (bp) upstream of TSS or ATG depending on --promoter-type (default: 2000).",
        rich_help_panel=PROMOTER_PANEL,
    )] = 2000,
    promoter_type: Annotated[str, typer.Option(
        "-p", "--promoter-type", help="Reference point for promoter extraction: 'standard' (upstream of TSS), 'upstream_ATG' (upstream of main CDS ATG codon), or 'standard_plus_up_to_ATG' (upstream of TSS plus 5' UTR up to start codon).",
        rich_help_panel=PROMOTER_PANEL,
    )] = "standard",

    # 2. Reference FASTA Options
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

    # 3. Input / Output Options
    annotation_name: Annotated[str, typer.Option(
        "-a", "-an", "--annotation-name", help="Annotation version, name or tag.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-g", "-gn", "--genome-name", help="Genome assembly version, name or tag.",
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
    Scans a set of “query” genes to locate all occurrences of a specified DNA motif within their upstream promoter regions.
    """
    if verbose:
        quiet = False

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

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

    annotation.motifs.find_and_plot(query_genes=genes, motif=motif, motif_length=motif_length, glistname=query_tag, tf_motif_tag=motif_tag, output_dir=output_dir, subfolder=subfolder, quiet=quiet)
    
if __name__ == "__main__":
    app()
