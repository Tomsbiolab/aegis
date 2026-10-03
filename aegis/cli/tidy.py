import typer
import os

from typing import List, Optional
from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome
from ..utils.genefunctions import TAXONOMY_ORGANELLE_CODES, resolve_taxonomy_tables, NCBI_GENETIC_CODES
from .utils import split_callback

RNA_CLASSES = ["mRNA", "antisense_lncRNA", "antisense_RNA", 
                "miRNA_primary_transcript", "ncRNA", "lncRNA",
                "lnc_RNA", "pseudogenic_tRNA", "rRNA", "snoRNA",
                "snRNA", "tRNA", "pre_miRNA", "tRNA_pseudogene",
                "SRP_RNA", "RNase_MRP_RNA"]

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
        "-d", "--output-dir", help="Path to the output folder."
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation filename, without extension."
    )] = "{annotation-name}_tidy.gff3",
    main_only: Annotated[bool, typer.Option(
        "-m", "--main", help="Whether to include only a main transcript and main CDS per gene."
    )] = False,
    include_UTRs: Annotated[bool, typer.Option(
        "-u", "--include-UTRs", help="Include UTRs in output gff, normally these are not required by external tools as they can be deduced from exons and CDS features."
    )] = False,
    just_genes: Annotated[bool, typer.Option(
        "-g", "--just-genes", help="Whether to only include gene level features."
    )] = False,
    features: Annotated[List[str], typer.Option(
        "-f", "--features", help=f"Selects only certain transcripts (e.g., 'mRNA,lncRNA'). Provide a comma-separated list. If empty, all biotypes are included. This option automatically enables 'clean_features'.",
        callback=split_callback
    )] = [],
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
    remove_symbols: Annotated[bool, typer.Option(
        "--remove-symbols", help="Removes symbol attributes from gff output."
    )] = False,
    remove_aliases: Annotated[bool, typer.Option(
        "--remove-aliases", help="Removes alias attributes from gff output."
    )] = False,
    clean_attributes: Annotated[bool, typer.Option(
        "--clean-attributes", help="Removes non-standard attributes from a gff, may help with external tool compatibility issues."
    )] = False,
    clean_features: Annotated[bool, typer.Option(
        "--clean-features", help="Removes non-standard features from a gff, may help with external tool compatibility issues."
    )] = False,
    add_gene_id: Annotated[bool, typer.Option(
        "--add-gene-id", help="Creates a special attribute 'Gene_id' which has the Gene_id as a value but is given to all gene features and subfeatures. This is useful for example when using featureCounts at the exon level but summarising counts at the Gene_id level."
    )] = False,
    symbols_as_description: Annotated[bool, typer.Option(
        "--symbols-as-description", help="Places gene symbols as 'Description=' attributes. Useful for JBrowse(2) display."
    )] = False,
    repeat_exons_utrs: Annotated[bool, typer.Option(
        "--repeat-exons-utrs", help="Creates individual exon/UTR entries with individual parental references for cases where a feature has more than one transcript level parent."
    )] = False,
    unique_cds_entry_ids: Annotated[bool, typer.Option(
        "--unique-cds-entry-ids", help="CDS entries corresponding to a same protein in a gff by default share the same id. However since the default format is incompatible with some external tools, this flag will ensure each CDS entry (line) has a unique id."
    )] = False, 
    for_lifton: Annotated[bool, typer.Option(
        "--for-lifton", help="Ensures output has individual CDS entry ids (-u) as it is required for LifOn compatibility in its current version."
    )] = False,
    genome_file: Annotated[str, typer.Option(
        "--genome-file", help="Path to genome FASTA file. Required when using --rework-all-CDSs or --infer-missing-CDSs."
    )] = "",
    genome_name: Annotated[str, typer.Option(
        "--genome-name", help="A name or tag for the genome assembly."
    )] = "",
    infer_missing_CDSs: Annotated[bool, typer.Option(
        "--infer-missing-CDSs", help="Detects and creates CDSs where missing, without overriding existing CDS annotations. Requires --genome-file."
    )] = False,
    rework_all_CDSs: Annotated[bool, typer.Option(
        "--rework-all-CDSs", help="Recalculates ALL CDSs from the genome sequence, overriding existing ones. More aggressive than --infer-missing-CDSs. Requires --genome-file."
    )] = False,
    fallback_to_trim: Annotated[bool, typer.Option(
        "--fallback-to-trim", help="When recalculating CDSs, fallback to trimming unaligned ends if no high-ratio ORF is found."
    )] = False,
    coding_ratio_threshold: Annotated[float, typer.Option(
        "--coding-ratio-threshold", help="Threshold ratio of coding sequence length to transcript length for rework CDS (default: 0.7)."
    )] = 0.7,
    allow_internal_stops: Annotated[bool, typer.Option(
        "--allow-internal-stops/--no-allow-internal-stops", help="Allow internal stop codons in progressive rework fallback (default: True)."
    )] = True,
    allow_partial: Annotated[bool, typer.Option(
        "--allow-partial/--no-allow-partial", help="Allow partial ORFs without stop codon in progressive rework fallback (default: True)."
    )] = True,
    enforce_start_codon: Annotated[bool, typer.Option(
        "--enforce-start-codon/--no-enforce-start-codon", help="Require start codon (ATG) in initial rework passes (default: True)."
    )] = True,
    orf_choice_mode: Annotated[str, typer.Option(
        "--orf-choice-mode", help="ORF selection criteria: 'longest' or 'earliest' (default: 'longest')."
    )] = "longest",
    skip_coordinate_polishing: Annotated[bool, typer.Option(
        "--skip-coordinate-polishing", help="Do not mutate feature coordinates when boundaries differ; log discrepancies as warnings instead."
    )] = False,
    no_collapse_exons: Annotated[bool, typer.Option(
        "--no-collapse-exons", help="Do not merge overlapping/adjacent exons."
    )] = False,
    no_collapse_CDSs: Annotated[bool, typer.Option(
        "--no-collapse-CDSs", help="Do not merge overlapping/adjacent CDS segments."
    )] = False,
    consider_read_utrs: Annotated[bool, typer.Option(
        "--consider-read-utrs", help="Consider UTRs as read from the GFF rather than inferred."
    )] = False,
    keep_original_subfeature_ids: Annotated[bool, typer.Option(
        "--keep-original-subfeature-ids", help="Keep original subfeature ids for CDS, UTR and exon features. By default, since tidy detects shared exons and UTRs between transcripts of the same gene, it will rename these subfeatures accordingly."
    )] = False,
    standard_features: Annotated[bool, typer.Option(
        "--standard-features", help="Standardises feature names to the most common names, for instance 'transcript' or 'pseudotranscript' just become 'mRNA' for downstream tool compatibility."
    )] = False,
    remove_genes_with_no_transcripts: Annotated[bool, typer.Option(
        "--remove-genes-with-no-transcripts", help="Removes genes with no transcripts."
    )] = False,
    remove_transcripts_with_no_exons: Annotated[bool, typer.Option(
        "--remove-transcripts-with-no-exons", help="Removes transcripts with no exons."
    )] = False,
    keep_missing_transcript_parent_references: Annotated[bool, typer.Option(
        "--keep-missing-transcript-parent-references", help="Keep parental references to missing transcripts. These are removed by default."
    )] = False,
    print_empty_attributes: Annotated[bool, typer.Option(
        "--print-empty-attributes", help="Print empty attributes. These are normally skipped."
    )] = False,
    print_orphaned_features: Annotated[bool, typer.Option(
        "--print-orphaned-features", help="Print orphaned features. Orphaned features are features which are not assigned to any gene, or genes which could not be incorporated into the annotation object. These are normally skipped."
    )] = False,
    adjust_internal_shifts: Annotated[str, typer.Option(
        "--adjust-internal-shifts", help="Frameshift / phase handling mode: 'intra_exon' (default), 'all', or 'none'."
    )] = "intra_exon",
    taxonomy: Annotated[str, typer.Option(
        "-tax", "--taxonomy", help="Taxonomic group preset for organelle genetic codes: 'plant' (default: nuclear=1, mito=1, plastid=11), 'vertebrate' (nuclear=1, mito=2), 'invertebrate' (nuclear=1, mito=5), or 'yeast' (nuclear=1, mito=3). Specific codes can be individually customized with --genetic-code, --mito-code, or --plastid-code."
    )] = "plant",
    genetic_code: Annotated[int, typer.Option(
        "-gc", "--genetic-code", help="NCBI genetic code table number for nuclear genes (default: 1)."
    )] = 1,
    auto_organelle_codes: Annotated[bool, typer.Option(
        "--auto-organelle-codes/--no-auto-organelle-codes", help="Automatically use mitochondrial and plastid genetic codes for organelle contigs (default: True)."
    )] = True,
    mito_code: Annotated[Optional[int], typer.Option(
        "-mc", "--mito-code", help="NCBI genetic code table number for mitochondrial contigs (overrides --taxonomy default)."
    )] = None,
    plastid_code: Annotated[Optional[int], typer.Option(
        "-pc", "--plastid-code", help="NCBI genetic code table number for plastid/chloroplast contigs (overrides --taxonomy default)."
    )] = None,
    mitochondria_chroms: Annotated[List[str], typer.Option(
        "--mitochondria-chroms", help="Explicit list or comma-separated names of mitochondrial chromosomes/scaffolds to translate with --mito-code.",
        callback=split_callback
    )] = [],
    chloroplast_chroms: Annotated[List[str], typer.Option(
        "--chloroplast-chroms", help="Explicit list or comma-separated names of chloroplast/plastid chromosomes/scaffolds to translate with --plastid-code.",
        callback=split_callback
    )] = [],
):
    """
    Cleans and reformats a GFF/GTF file to correct common formatting errors and improve compatibility with other bioinformatics tools.
    
    This script parses an annotation file, allows for extensive filtering and reformatting, and exports a standardized GFF3 file.
    """

    collapse_exons = not(no_collapse_exons)
    collapse_CDSs = not(no_collapse_CDSs)
    remove_missing_transcript_parent_references = not(keep_missing_transcript_parent_references)

    if keep_original_subfeature_ids:
        rename_features = ()
    else:
        rename_features = ("CDS", "UTR", "exon")

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    if output_dir == "./aegis_output/":
        subfolder = True
    else:
        subfolder = False

    for feature in features:
        if feature not in RNA_CLASSES:
            raise typer.BadParameter(f"Invalid feature: {feature}. Choose from: {RNA_CLASSES}")

    if adjust_internal_shifts not in ("intra_exon", "all", "none"):
        raise typer.BadParameter(f"Invalid adjust_internal_shifts: '{adjust_internal_shifts}'. Choose from: 'intra_exon', 'all', 'none'.")

    if orf_choice_mode not in ("longest", "earliest"):
        raise typer.BadParameter(f"Invalid orf_choice_mode: '{orf_choice_mode}'. Choose from: 'longest', 'earliest'.")

    mito_chroms = mitochondria_chroms if len(mitochondria_chroms) > 0 else None
    chloro_chroms = chloroplast_chroms if len(chloroplast_chroms) > 0 else None

    if (rework_all_CDSs or infer_missing_CDSs) and not genome_file:
        raise typer.BadParameter("A genome FASTA file must be provided via --genome-file when using --rework-all-CDSs or --infer-missing-CDSs.")

    if genome_file:
        if not genome_name:
            genome_name = os.path.splitext(os.path.basename(genome_file))[0]
        genome = Genome(name=genome_name, genome_file_path=genome_file, quiet=quiet)
    else:
        genome = None

    os.makedirs(output_dir, exist_ok=True)

    res_gc, res_mito, res_plastid = resolve_taxonomy_tables(
        taxonomy=taxonomy,
        table=genetic_code,
        mito_table=mito_code,
        plastid_table=plastid_code,
    )

    annotation = Annotation(
        name=annotation_name,
        annot_file_path=annotation_file,
        genome=genome,
        rework_all_CDSs=rework_all_CDSs,
        work_out_missing_CDSs=infer_missing_CDSs,
        fallback_to_trim=fallback_to_trim,
        quiet=quiet,
        collapse_exons=collapse_exons,
        collapse_CDSs=collapse_CDSs,
        consider_read_utrs=consider_read_utrs,
        rename_features=rename_features,
        standardise_features=standard_features,
        remove_genes_with_no_transcripts=remove_genes_with_no_transcripts,
        remove_transcripts_with_no_exons=remove_transcripts_with_no_exons,
        remove_missing_transcript_parent_references=remove_missing_transcript_parent_references,
        remove_genes_with_no_transcripts_even_if_pseudogene=remove_genes_with_no_transcripts,
        skip_orphaned_features=not(print_orphaned_features),
        adjust_internal_shifts=adjust_internal_shifts,
        taxonomy=taxonomy,
        table=res_gc,
        auto_organelle_codes=auto_organelle_codes,
        mito_table=res_mito,
        plastid_table=res_plastid,
        mitochondria_chroms=mito_chroms,
        chloroplast_chroms=chloro_chroms,
        skip_coordinate_polishing=skip_coordinate_polishing,
        coding_ratio_threshold=coding_ratio_threshold,
        allow_internal_stops=allow_internal_stops,
        allow_partial=allow_partial,
        enforce_start_codon=enforce_start_codon,
        orf_choice_mode=orf_choice_mode,
    )

    if output_file == "{annotation-name}_tidy.gff3":
        output_file = f"{annotation_name}_tidy.gff3"

    if unique_cds_entry_ids or for_lifton:
        annotation.CDS_to_CDS_segment_ids()
    else:
        annotation.CDS_segment_to_CDS_ids()

    if symbols_as_description:
        remove_symbols = True

    if features:
        annotation.filter_by_rna_class(rna_classes=features, remove_genes_accordingly=True, quiet=quiet)
        clean_features = True

    annotation.export.gff(output_dir=output_dir, filename=output_file, main_only=main_only, UTRs=include_UTRs, just_genes=just_genes, repeat_exons_utrs=repeat_exons_utrs, skip_atypical_fts=clean_features, quiet=quiet, aliases=(not remove_aliases), symbols=(not remove_symbols), symbols_as_description=symbols_as_description, clean_attributes=clean_attributes, featurecountsID=add_gene_id, print_empty_attributes=print_empty_attributes, subfolder=subfolder)

if __name__ == "__main__":
    app()
