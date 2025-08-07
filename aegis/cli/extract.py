import typer
import os
from typing_extensions import Annotated
from aegis.genome import Genome
from aegis.annotation import Annotation
from textwrap import dedent

FEATURES = ["gene", "transcript", "CDS", "protein", "promoter"]

VALID_IDS = ["gene", "transcript", "CDS", "feature"]

EXTRACTION_MODES = ["all", "main", "unique", "unique_per_gene"]

PROMOTER_TYPES = ["standard", "upstream_ATG", "standard_plus_up_to_ATG"]

RNA_CLASSES = ["mRNA", "antisense_lncRNA", "antisense_RNA", 
                "miRNA_primary_transcript", "ncRNA", "lncRNA",
                "lnc_RNA", "pseudogenic_tRNA", "rRNA", "snoRNA",
                "snRNA", "tRNA", "pre_miRNA", "tRNA_pseudogene",
                "SRP_RNA", "RNase_MRP_RNA"]

app = typer.Typer(add_completion=False)

def split_callback(value:str):
    if value:
        return [item.strip() for item in value.split(",")]
    return []

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    genome_file: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file."
    )],
    genome_name: Annotated[str, typer.Option(
        "-g", "--genome-name", help="A name or tag for the genome assembly (e.g., 'TAIR10'). [default: a name derived from the genome FASTA filename]"
    )] = "{genome-file}",
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="A name or tag for the annotation version (e.g., 'Araport11'). [default: a name derived from the annotation filename]"
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-o", "--output-dir", help="Path to the directory where output FASTA files will be saved."
    )] = "./aegis_output/features/",
    feature_type: Annotated[str, typer.Option(
        "-f", "--feature-type", help=f"Feature type(s) to extract, as a comma-separated list. Available options: {', '.join(FEATURES)}.",
        callback=split_callback
    )] = "gene",

    mode: Annotated[str, typer.Option(
        "-m", "--mode", help=f"""Extraction mode(s), as a comma-separated list. Controls filtering of features.\n\n

        - 'all': Extract all features (e.g., all transcripts for a gene).\n
        - 'unique_per_gene': Keep one copy of each unique protein/CDS sequence per gene.\n
        - 'main': Extract only the main variant (e.g., the longest transcript).\n
        - 'unique': Keep only one copy of each unique protein sequence across the entire output.""",
        callback=split_callback
    )] = "all,main",
    rna_classes: Annotated[str, typer.Option(
        "-r", "--rna-classes", help=f"Filter transcripts by biotype (e.g., 'mRNA,lncRNA'). Provide a comma-separated list. If empty, all biotypes are included.",
        callback=split_callback
    )] = "",
    promoter_size: Annotated[int, typer.Option(
        "-ps", "--promoter-size", help=f"Size of the promoter region in base pairs (bp). Used only if 'promoter' is a selected feature."
    )] = 2000,

    promoter_type: Annotated[str, typer.Option(
        "-p", "--promoter-type", help="""\
                Defines the reference point for extracting promoter regions. Used only if 'promoter' is selected.\n\n
                
                Options:\n
                - 'standard': Upstream of the transcript's start site (TSS).\n
                - 'upstream_ATG': Upstream of the main CDS's start codon (ATG). Falls back to 'standard' if no CDS is present.\n
                - 'standard_plus_up_to_ATG': The 'standard' promoter plus the 5' UTR (sequence between TSS and ATG). Falls back to 'standard' if no CDS.
                """
    )] = "standard",

    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help=f"Add extra details in fasta headers; scaffold/chromosome number, genome co-ordinates, and/or protein tags if applicable."
    )] = False,
    feature_id: Annotated[str, typer.Option(
        "-i", "--feature-id", help=f"Specifies which feature ID to use in FASTA headers. E.g., use 'gene' to label all outputs (transcripts, proteins) with their parent gene ID. 'feature' uses the most specific ID available. Available: {', '.join(VALID_IDS)}."
    )] = "feature"
):
    """
    Extract sequences from a genome based on an annotation file.

    This command supports multiple output formats and allows selecting
    specific features (e.g. gene, transcript, CDS, protein, promoter).
    
    Use the --mode and --feature-id flags to control sequence filtering
    and ID labeling. Promoter generation supports multiple strategies,
    including upstream of TSS or ATG.
    """

    for f_type in feature_type:
        if f_type not in FEATURES:
            raise typer.BadParameter(f"Invalid feature type: {f_type}. Choose from: {FEATURES}")

    for m_type in mode:
        if m_type not in EXTRACTION_MODES:
            raise typer.BadParameter(f"Invalid mode: {m_type}. Choose from: {EXTRACTION_MODES}")

    if feature_id not in VALID_IDS:
        raise typer.BadParameter(f"Invalid feature ID: {feature_id}. Choose from: {VALID_IDS}")
    
    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    if genome_name == "{genome-file}":
        genome_name = os.path.splitext(os.path.basename(genome_file))[0]

    if promoter_type not in PROMOTER_TYPES:
        raise typer.BadParameter(f"Invalid promoter type: '{promoter_type}'. Choose from: {PROMOTER_TYPES}")

    for rna_class in rna_classes:
        if rna_class not in RNA_CLASSES:
            raise typer.BadParameter(f"Invalid rna class: {rna_class}. Choose from: {RNA_CLASSES}")

    genome = Genome(name=genome_name, genome_file_path=genome_file)
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
            annotation.export_transcripts(only_main=False, verbose=verbose, custom_path=output_dir, used_id=used_id, rna_classes=rna_classes)
        if "main" in mode or "unique" in mode or "unique_per_gene" in mode:
            annotation.export_transcripts(custom_path=output_dir, verbose=verbose, used_id=used_id, rna_classes=rna_classes)

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