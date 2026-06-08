import typer
import os
import pandas as pd

from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome


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
    header: Annotated[bool, typer.Option(
        "-H", "--header", help="Use this flag to indicate the presence of a header in input genelist file."
    )] = False,
    promoter_size: Annotated[int, typer.Option(
        "-ps", "--promoter-size", help=f"Only applies if promoter included in '-f'. Promoter size in bp upstream of TSS or ATG depending on '-p'."
    )] = 2000,
    promoter_type: Annotated[str, typer.Option(
        "-p", "--promoter-type", help=f"Only applies if promoter included in '-f'. Defines the reference point for the promoter regions of '-ps' size. 'standard': Generated upstream of the transcript's start site (TSS); 'upstream_ATG': Generated upstream of the main CDS's start codon (ATG). If no CDS, falls back to standard; 'standard_plus_up_to_ATG': Generated upstream of the transcript's start site (TSS) and any gene sequence up to the start codon (ATG) is also added. If no CDS, falls back to standard."
    )] = "standard",
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-g", "--genome-name", help="Genome assembly version, name or tag."
    )] = "{genome-file}",
    query_tag: Annotated[str, typer.Option(
        "--genelist-tag", help="Query gene list tag/name to improve output description."
    )] = "query_genes",
    motif_tag: Annotated[str, typer.Option(
        "--motif-tag", help="Motif tag/name to improve output description, e.g. '{TF}_{motif_name}'."
    )] = "query_motif",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output directory."
    )] = "./aegis_output/",
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,

):
    """
    Scans a set of “query” genes to locate all occurrences of a specified DNA motif within their upstream promoter regions.
    """

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

    genome = Genome(name=genome_name, genome_file_path=genome_file, quiet=quiet)
    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, genome=genome, quiet=quiet)

    annotation.generate_promoters(promoter_size=promoter_size, promoter_type=promoter_type)

    annotation.motifs.find_and_plot(query_genes=genes, motif=motif, motif_length=motif_length, glistname=query_tag, tf_motif_tag=motif_tag, output_dir=output_dir, subfolder=subfolder, quiet=quiet)
    
if __name__ == "__main__":
    app()
