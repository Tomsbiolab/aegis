import typer
import os
from typing_extensions import Annotated
from aegis.genome import Genome
from aegis.annotation import Annotation

features = ["gene", "transcript", "CDS", "protein", "promoter"]

IDs = ["gene", "transcript", "CDS", "feature"]

modes = ["all", "main", "unique", "unique_per_gene"]

promoter_types = ["standard", "upstream_ATG", "standard_plus_up_to_ATG"]

app = typer.Typer(add_completion=False)

def split_callback(value: str):
    if value:
        return [item.strip() for item in value.split(',')]
    return []

@app.command()
def main(
    genome_fasta: Annotated[str, typer.Option(
        "-g", "--genome-fasta", help="Path to the input genome FASTA file."
    )],
    annotation_file: Annotated[str, typer.Option(
        "-a", "--annotation-file", help="Path to the input annotation GFF/GTF file."
    )],
    genome_name: Annotated[str, typer.Option(
        "-gn", "--genome-name", help="Genome assembly version, name or tag."
    )] = "{genome-fasta}",
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-o", "--output-dir", help="Path to the directory where output FASTA files will be saved."
    )] = "./aegis_output/",
    feature_type: Annotated[str, typer.Option(
        "-f", "--feature-type", help=f"One or more feature types to extract, separated by commas. Choose from: {features}.",
        callback=split_callback
    )] = "gene",

    mode: Annotated[str, typer.Option(
        "-m", "--mode", help=f"Select extraction modes separated by commas. Choose from: {features}. 'All' for all transcripts/CDSs/proteins of a gene, 'main' for just the main variant of a gene defined by default as the longest transcript/CDS/protein, 'unique_per_gene' removes redundant CDS/protein sequences per gene, and 'unique' does the same for proteins at the whole output level.",
        callback=split_callback
    )] = "all,main",

    promoter_size: Annotated[int, typer.Option(
        "-ps", "--promoter-size", help=f"Only applies if promoter included in '-f'. Promoter size in bp upstream of TSS or ATG depending on '-p'."
    )] = 2000,

    promoter_type: Annotated[str, typer.Option(
        "-p", "--promoter-type", help=f"Only applies if promoter included in '-f'. Defines the reference point for the promoter regions of '-ps' size. 'standard': Generated upstream of the transcript's start site (TSS); 'upstream_ATG': Generated upstream of the main CDS's start codon (ATG). If no CDS, falls back to standard; 'standard_plus_up_to_ATG': Generated upstream of the transcript's start site (TSS) and any gene sequence up to the start codon (ATG) is also added. If no CDS, falls back to standard."
    )] = "standard",

    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help=f"Whether to include extra details in fasta headers; scaffold/chromosome number, genome co-ordinates, and/or protein tags if applicable."
    )] = False,
    feature_id: Annotated[str, typer.Option(
        "-i", "--feature-id", help=f"Most specific feature ID used in fasta header outputs. E.g. you may want to export transcripts but associated to gene ids directly instead of using the transcript feature IDs. Choose from: {IDs}."
    )] = "feature"
):
    """
    Extract sequences from a genome based on an annotation.

    This command supports multiple output formats and allows selecting
    specific features (e.g. gene, transcript, CDS, protein, promoter).
    
    Use the --mode and --feature-id flags to control sequence filtering
    and ID labeling. Promoter generation supports multiple strategies,
    including upstream of TSS or ATG.
    """

    for f_type in feature_type:
        if f_type not in features:
            raise typer.BadParameter(f"Invalid feature type: {f_type}. Choose from: {features}")

    for m_type in mode:
        if m_type not in modes:
            raise typer.BadParameter(f"Invalid mode: {m_type}. Choose from: {modes}")

    if feature_id not in IDs:
        raise typer.BadParameter(f"Invalid feature ID: {feature_id}. Choose from: {IDs}")
    
    if annotation_name == "filename":
        annotation_name = os.path.splitext(annotation_file)[0]

    if genome_name == "filename":
        genome_name = os.path.splitext(genome_fasta)[0]

    genome = Genome(name=genome_name, genome_file_path=genome_fasta)
    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, genome=genome)

    if "promoter" in feature_type:
        annotation.generate_promoters(genome, promoter_size=promoter_size, promoter_type=promoter_type)

    annotation.generate_sequences(genome)

    if "gene" in feature_type:

        annotation.export_genes(custom_path=output_dir, verbose=verbose)

    if "transcript" in feature_type:

        if "gene" in feature_id:
            used_id = "gene"
        else:
            used_id = "transcript"

        if "all" in mode:
            annotation.export_transcripts(only_main=False, verbose=verbose, custom_path=output_dir, used_id=used_id)
        if "main" in mode or "unique" in mode or "unique_per_gene" in mode:
            annotation.export_transcripts(custom_path=output_dir, verbose=verbose, used_id=used_id)

    if "protein" in feature_type:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        elif "CDS" in feature_id:
            used_id = "CDS"
        else:
            used_id = "protein"

        if "unique_per_gene" in mode:
            annotation.export_proteins(only_main=False, custom_path=output_dir, verbose=verbose, unique_proteins_per_gene=True, used_id=used_id)
        elif "unique" in mode:
            annotation.export_unique_proteins(custom_path=output_dir, verbose=verbose)
        else:
            if "all" in mode:
                annotation.export_proteins(only_main=False, custom_path=output_dir, verbose=verbose, used_id=used_id)
            if "main" in mode or "unique" in mode or "unique_per_gene" in mode:
                annotation.export_proteins(custom_path=output_dir, verbose=verbose, used_id=used_id)

    if "CDS" in feature_type:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        else:
            used_id = "CDS"

        if "all" in mode:
            annotation.export_CDSs(only_main=False, custom_path=output_dir, verbose=verbose, used_id=used_id)
        if "main" in mode or "unique" in mode or "unique_per_gene" in mode:
            annotation.export_CDSs(custom_path=output_dir, verbose=verbose, used_id=used_id)

    if "promoter" in feature_type:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        else:
            used_id = "promoter"

        if "all" in mode:
            annotation.export_promoters(only_main=False, custom_path=output_dir, verbose=verbose, used_id=used_id)
        if "main" in mode or "unique" in mode or "unique_per_gene" in mode:
            annotation.export_promoters(custom_path=output_dir, verbose=verbose, used_id=used_id)

if __name__ == "__main__":
    app()