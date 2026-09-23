import typer
import os
from typing import List, Optional
from typing_extensions import Annotated

from ..annotation import Annotation
from .utils import split_callback

RNA_CLASSES = [
    "mRNA", "antisense_lncRNA", "antisense_RNA",
    "miRNA_primary_transcript", "ncRNA", "lncRNA",
    "lnc_RNA", "pseudogenic_tRNA", "rRNA", "snoRNA",
    "snRNA", "tRNA", "pre_miRNA", "tRNA_pseudogene",
    "SRP_RNA", "RNase_MRP_RNA", "Y_RNA", "YRNA",
    "scaRNA", "vault_RNA", "telomerase_RNA", "scRNA",
    "RNase_P_RNA"
]

app = typer.Typer(add_completion=False, no_args_is_help=True)


@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output directory."
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation filename, with or without extension."
    )] = "{annotation-name}_filtered",
    coding_only: Annotated[bool, typer.Option(
        "--coding-only", help="Keep only protein-coding genes and transcripts (removes non-coding genes and non-coding transcripts from mixed genes)."
    )] = False,
    non_coding_only: Annotated[bool, typer.Option(
        "--non-coding-only", help="Keep only non-protein-coding genes and transcripts (removes coding genes and coding transcripts from mixed genes, keeping lncRNAs, etc.)."
    )] = False,
    rna_classes: Annotated[List[str], typer.Option(
        "-r", "--rna-classes", help="Filter transcripts by biotype (e.g. 'mRNA,lncRNA'). Provide a comma-separated list. Transcripts not in the list and genes without remaining transcripts are removed.",
        callback=split_callback
    )] = [],
    skip_pseudogenes: Annotated[bool, typer.Option(
        "--skip-pseudogenes", help="Remove pseudogenes from the annotation."
    )] = False,
    pseudogenes_only: Annotated[bool, typer.Option(
        "--pseudogenes-only", help="Keep only pseudogenes, removing non-pseudogene genes."
    )] = False,
    skip_te: Annotated[bool, typer.Option(
        "--skip-te", "--skip-transposables", help="Remove transposable element genes from the annotation."
    )] = False,
    te_only: Annotated[bool, typer.Option(
        "--te-only", "--transposables-only", help="Keep only transposable element genes, removing non-TE genes."
    )] = False,
    min_cds_size: Annotated[Optional[int], typer.Option(
        "--min-cds-size", help="Remove genes whose main CDS length (bp) is smaller than this threshold."
    )] = None,
    has_symbol: Annotated[bool, typer.Option(
        "--has-symbol", help="Keep only genes that have an assigned gene symbol."
    )] = False,
    main_only: Annotated[bool, typer.Option(
        "-m", "--main", help="Include only the main transcript and main CDS per gene."
    )] = False,
    include_UTRs: Annotated[bool, typer.Option(
        "-u", "--include-UTRs", help="Include UTRs in output GFF."
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
):
    """
    Filter an annotation file based on biotypes, RNA classes, transposable elements, pseudogenes, CDS size, or gene symbols.
    """
    if coding_only and non_coding_only:
        raise typer.BadParameter("Cannot specify both --coding-only and --non-coding-only.")

    if skip_pseudogenes and pseudogenes_only:
        raise typer.BadParameter("Cannot specify both --skip-pseudogenes and --pseudogenes-only.")

    if skip_te and te_only:
        raise typer.BadParameter("Cannot specify both --skip-te/--skip-transposables and --te-only/--transposables-only.")

    if min_cds_size is not None and min_cds_size <= 0:
        raise typer.BadParameter("--min-cds-size must be greater than 0.")

    for rna_class in rna_classes:
        if rna_class not in RNA_CLASSES:
            raise typer.BadParameter(f"Invalid RNA class: '{rna_class}'. Choose from: {RNA_CLASSES}")

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_dir, exist_ok=True)
    subfolder = (output_dir == "./aegis_output/")

    if output_file == "{annotation-name}_filtered":
        output_file = f"{annotation_name}_filtered"

    if not (output_file.endswith(".gff3") or output_file.endswith(".gff")):
        output_file += ".gff3"

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet)

    # 1. Biotype filtering (coding vs non-coding)
    if coding_only:
        annotation.remove_non_coding_genes_and_transcripts(quiet=quiet)
    elif non_coding_only:
        annotation.remove_coding_genes_and_transcripts(quiet=quiet)

    # 2. Specific RNA classes / biotypes
    if rna_classes:
        annotation.filter_by_rna_class(rna_classes=rna_classes, remove_genes_accordingly=True, quiet=quiet)

    # 3. Pseudogene filtering
    if skip_pseudogenes:
        annotation.remove_pseudogene_genes(quiet=quiet)
    elif pseudogenes_only:
        annotation.remove_non_pseudogene_genes(quiet=quiet)

    # 4. Transposable elements filtering
    if skip_te:
        annotation.remove_TE_genes(quiet=quiet)
    elif te_only:
        annotation.remove_non_TE_genes(quiet=quiet)

    # 5. CDS size filtering
    if min_cds_size is not None:
        annotation.remove_genes_with_small_CDSs(CDS_threshold=min_cds_size, quiet=quiet)

    # 6. Gene symbol filtering
    if has_symbol:
        annotation.remove_genes_without_symbols(quiet=quiet)

    annotation.export.gff(
        output_dir=output_dir,
        filename=output_file,
        main_only=main_only,
        UTRs=include_UTRs,
        subfolder=subfolder,
        quiet=quiet,
    )


if __name__ == "__main__":
    app()
