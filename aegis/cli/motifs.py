import typer
import os
from typing_extensions import Annotated
from aegis.annotation import Annotation
from aegis.genome import Genome
import pandas as pd

app = typer.Typer(add_completion=False)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    genome_fasta: Annotated[str, typer.Argument(
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
        "-h", "--header", help="Use this flag to indicate the presence of a header in input symbols_file."
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
    )] = "{genome-fasta}",
    query_tag: Annotated[str, typer.Option(
        "-q", "--genelist_tag", help="Query gene list tag/name to improve output description."
    )] = "query_genes",
    motif_tag: Annotated[str, typer.Option(
        "-m", "--motif_tag", help="Motif tag/name to improve output description, e.g. '{TF}_{motif_name}'."
    )] = "query_motif",
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/",

):
    """
    Scans a set of “query” genes to locate all occurrences of a specified DNA motif within their upstream promoter regions.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_folder, exist_ok=True)

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

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file)
    genome = Genome(name=genome_name, genome_file_path=genome_fasta)

    annotation.generate_promoters(genome=genome, promoter_size=promoter_size, promoter_type=promoter_type, generate_sequence=True)

    annotation.find_motifs(query_genes=genes, motif=motif, motif_length=motif_length, glistname=query_tag, tf_motif_tag=motif_tag, custom_path=output_folder)
    
if __name__ == "__main__":
    app()
