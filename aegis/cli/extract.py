import typer
import os

from typing_extensions import Annotated
from textwrap import dedent

from ..genome import Genome
from ..annotation import Annotation

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
        "-d", "--output-dir", help="Path to the directory where output FASTA files will be saved."
    )] = "./aegis_output/features/",
    features: Annotated[str, typer.Option(
        "-f", "--features", help=f"Feature type(s) to extract, as a comma-separated list. Available options: {', '.join(FEATURES)}.",
        callback=split_callback
    )] = "gene",

    mode: Annotated[str, typer.Option(
        "-m", "--mode", help=f"""Extraction mode(s), as a comma-separated list. Controls filtering of features.\n\n

        - 'all': Extract all features (e.g., all transcripts for a gene).\n
        - 'unique_per_gene': Keep one copy of each unique protein/CDS sequence per gene.\n
        - 'main': Extract only the main variant (e.g., the longest transcript).\n
        - 'unique': Keep only one copy of each unique protein/CDS sequence across the entire output.""",
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

    detailed_headers: Annotated[bool, typer.Option(
        "-dh", "--detailed-headers", help=f"Add extra details in fasta headers; scaffold/chromosome number, genome co-ordinates, and/or protein tags if applicable."
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
    no_collapse_exons: Annotated[bool, typer.Option(
        "--no-collapse-exons", help="Do not merge overlapping/adjacent exons."
    )] = False,
    no_collapse_CDSs: Annotated[bool, typer.Option(
        "--no-collapse-CDSs", help="Do not merge overlapping/adjacent CDS segments."
    )] = False,
    feature_id: Annotated[str, typer.Option(
        "--feature-id", help=f"Specifies which feature ID to use in FASTA headers. E.g., use 'gene' to label all outputs (transcripts, proteins) with their parent gene ID. 'feature' uses the most specific ID available. Available: {', '.join(VALID_IDS)}."
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

    collapse_exons = not no_collapse_exons
    collapse_CDSs = not no_collapse_CDSs

    for f_type in features:
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

    genome = Genome(name=genome_name, genome_file_path=genome_file, quiet=quiet)
    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, genome=genome, quiet=quiet, collapse_exons=collapse_exons, collapse_CDSs=collapse_CDSs)

    if "promoter" in features:
        annotation.generate_promoters(genome, promoter_size=promoter_size, promoter_type=promoter_type)

    annotation.generate_sequences(genome, quiet=quiet)

    if "gene" in features:

        annotation.export.genes(custom_path=output_dir, verbose=detailed_headers)

    if "transcript" in features:

        if "gene" in feature_id:
            used_id = "gene"
        else:
            used_id = "transcript"

        if "unique_per_gene" in mode:
            annotation.export.transcripts(only_main=False, verbose=detailed_headers, custom_path=output_dir, used_id=used_id, rna_classes=rna_classes, unique_transcripts_per_gene=True)
        elif "unique" in mode:
            annotation.export.unique_transcripts(custom_path=output_dir, quiet=quiet, rna_classes=rna_classes)
        elif "all" in mode:
            annotation.export.transcripts(only_main=False, verbose=detailed_headers, custom_path=output_dir, used_id=used_id, rna_classes=rna_classes)
        else:
            annotation.export.transcripts(custom_path=output_dir, verbose=detailed_headers, used_id=used_id, rna_classes=rna_classes)

    if "protein" in features:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        elif "CDS" in feature_id:
            used_id = "CDS"
        else:
            used_id = "protein"

        if "unique_per_gene" in mode:
            annotation.export.proteins(only_main=False, custom_path=output_dir, verbose=detailed_headers, unique_proteins_per_gene=True, used_id=used_id)
        elif "unique" in mode:
            annotation.export.unique_proteins(custom_path=output_dir, quiet=quiet)
        elif "all" in mode:
            annotation.export.proteins(only_main=False, custom_path=output_dir, verbose=detailed_headers, used_id=used_id, only_cds_main=False)
        else:
            annotation.export.proteins(custom_path=output_dir, verbose=detailed_headers, used_id=used_id)

    if "CDS" in features:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        else:
            used_id = "CDS"

        if "unique_per_gene" in mode:
            annotation.export.CDSs(only_main=False, custom_path=output_dir, verbose=detailed_headers, used_id=used_id, unique_CDSs_per_gene=True)
        elif "unique" in mode:
            annotation.export.unique_CDSs(custom_path=output_dir, quiet=quiet)
        elif "all" in mode:
            annotation.export.CDSs(only_main=False, custom_path=output_dir, verbose=detailed_headers, used_id=used_id, only_cds_main=False)
        else:
            annotation.export.CDSs(custom_path=output_dir, verbose=detailed_headers, used_id=used_id)

    if "promoter" in features:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        else:
            used_id = "promoter"

        if "all" in mode or "unique_per_gene" in mode or "unique" in mode:
            annotation.export.promoters(only_main=False, custom_path=output_dir, verbose=detailed_headers, used_id=used_id)
        else:
            annotation.export.promoters(custom_path=output_dir, verbose=detailed_headers, used_id=used_id)

if __name__ == "__main__":
    app()