import typer
import os

from typing import List
from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome
from .utils import (
    split_callback,
    TaxonomyOption,
    GeneticCodeOption,
    AutoOrganelleCodesOption,
    MitoCodeOption,
    PlastidCodeOption,
    MitochondriaChromsOption,
    ChloroplastChromsOption,
    IO_PANEL,
    EXEC_PANEL,
    CDS_PANEL,
    FASTA_HEADER_PANEL,
)

RNA_CLASSES = ["mRNA", "antisense_lncRNA", "antisense_RNA", 
                "miRNA_primary_transcript", "ncRNA", "lncRNA",
                "lnc_RNA", "pseudogenic_tRNA", "rRNA", "snoRNA",
                "snRNA", "tRNA", "pre_miRNA", "tRNA_pseudogene",
                "SRP_RNA", "RNase_MRP_RNA"]

CLEANING_PANEL = "GFF Cleaning & Compatibility"
FEATURE_PANEL = "Feature Filtering"

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation filename, without extension.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name}_tidy.gff3",
    genome_file: Annotated[str, typer.Option(
        "-g", "--genome", "--genome-file", help="Path to genome FASTA file. Required when using --rework-all-CDSs or --infer-missing-CDSs.",
        rich_help_panel=IO_PANEL,
    )] = "",
    genome_name: Annotated[str, typer.Option(
        "-gn", "--genome-name", help="A name or tag for the genome assembly.",
        rich_help_panel=IO_PANEL,
    )] = "",

    main_only: Annotated[bool, typer.Option(
        "-m", "--main", help="Whether to include only a main transcript and main CDS per gene.",
        rich_help_panel=FEATURE_PANEL,
    )] = False,
    include_UTRs: Annotated[bool, typer.Option(
        "-u", "--include-UTRs", help="Include UTRs in output gff, normally these are not required by external tools as they can be deduced from exons and CDS features.",
        rich_help_panel=FEATURE_PANEL,
    )] = False,
    just_genes: Annotated[bool, typer.Option(
        "--just-genes", help="Whether to only include gene level features.",
        rich_help_panel=FEATURE_PANEL,
    )] = False,
    features: Annotated[List[str], typer.Option(
        "-f", "--features", "--biotypes", "--rna-classes", help=f"Selects only certain transcripts (e.g., 'mRNA,lncRNA'). Provide a comma-separated list. If empty, all biotypes are included. This option automatically enables 'clean_features'.",
        callback=split_callback,
        rich_help_panel=FEATURE_PANEL,
    )] = [],
    remove_genes_with_no_transcripts: Annotated[bool, typer.Option(
        "--remove-genes-with-no-transcripts", help="Removes genes with no transcripts.",
        rich_help_panel=FEATURE_PANEL,
    )] = False,
    remove_transcripts_with_no_exons: Annotated[bool, typer.Option(
        "--remove-transcripts-with-no-exons", help="Removes transcripts with no exons.",
        rich_help_panel=FEATURE_PANEL,
    )] = False,

    clean_features: Annotated[bool, typer.Option(
        "--clean-features", help="Removes non-standard features from a gff, may help with external tool compatibility issues.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    clean_attributes: Annotated[bool, typer.Option(
        "--clean-attributes", help="Removes non-standard attributes from a gff, may help with external tool compatibility issues.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    standard_features: Annotated[bool, typer.Option(
        "--standard-features", help="Standardises feature names to the most common names, for instance 'transcript' or 'pseudotranscript' just become 'mRNA' for downstream tool compatibility.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    for_lifton: Annotated[bool, typer.Option(
        "--for-lifton", help="Ensures output has individual CDS entry ids as required for Lifton/Liftoff compatibility in its current version.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    unique_cds_entry_ids: Annotated[bool, typer.Option(
        "--unique-cds-entry-ids", help="CDS entries corresponding to a same protein in a gff by default share the same id. However since the default format is incompatible with some external tools, this flag will ensure each CDS entry (line) has a unique id.",
        rich_help_panel=CLEANING_PANEL,
    )] = False, 
    repeat_exons_utrs: Annotated[bool, typer.Option(
        "--repeat-exons-utrs", help="Creates individual exon/UTR entries with individual parental references for cases where a feature has more than one transcript level parent.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    keep_original_subfeature_ids: Annotated[bool, typer.Option(
        "--keep-original-subfeature-ids", help="Keep original subfeature ids for CDS, UTR and exon features. By default, since tidy detects shared exons and UTRs between transcripts of the same gene, it will rename these subfeatures accordingly.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    remove_symbols: Annotated[bool, typer.Option(
        "--remove-symbols", help="Removes symbol attributes from gff output.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    remove_aliases: Annotated[bool, typer.Option(
        "--remove-aliases", help="Removes alias attributes from gff output.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    add_gene_id: Annotated[bool, typer.Option(
        "--add-gene-id", help="Creates a special attribute 'Gene_id' which has the Gene_id as a value but is given to all gene features and subfeatures. This is useful for example when using featureCounts at the exon level but summarising counts at the Gene_id level.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    symbols_as_description: Annotated[bool, typer.Option(
        "--symbols-as-description", help="Places gene symbols as 'Description=' attributes. Useful for JBrowse(2) display.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    consider_read_utrs: Annotated[bool, typer.Option(
        "--consider-read-utrs", help="Consider UTRs as read from the GFF rather than inferred.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    keep_missing_transcript_parent_references: Annotated[bool, typer.Option(
        "--keep-missing-transcript-parent-references", help="Keep parental references to missing transcripts. These are removed by default.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    print_empty_attributes: Annotated[bool, typer.Option(
        "--print-empty-attributes", help="Print empty attributes. These are normally skipped.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    print_orphaned_features: Annotated[bool, typer.Option(
        "--print-orphaned-features", help="Print orphaned features. Orphaned features are features which are not assigned to any gene, or genes which could not be incorporated into the annotation object. These are normally skipped.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    no_collapse_exons: Annotated[bool, typer.Option(
        "--no-collapse-exons", help="Do not merge overlapping/adjacent exons.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,
    no_collapse_CDSs: Annotated[bool, typer.Option(
        "--no-collapse-CDSs", help="Do not merge overlapping/adjacent CDS segments.",
        rich_help_panel=CLEANING_PANEL,
    )] = False,

    infer_missing_CDSs: Annotated[bool, typer.Option(
        "--infer-missing-CDSs", help="Detects and creates CDSs where missing, without overriding existing CDS annotations. Requires --genome-file.",
        rich_help_panel=CDS_PANEL,
    )] = False,
    rework_all_CDSs: Annotated[bool, typer.Option(
        "--rework-all-CDSs", help="Recalculates ALL CDSs from the genome sequence, overriding existing ones. More aggressive than --infer-missing-CDSs. Requires --genome-file.",
        rich_help_panel=CDS_PANEL,
    )] = False,
    fallback_to_trim: Annotated[bool, typer.Option(
        "--fallback-to-trim", help="When recalculating CDSs, fallback to trimming unaligned ends if no high-ratio ORF is found.",
        rich_help_panel=CDS_PANEL,
    )] = False,
    coding_ratio_threshold: Annotated[float, typer.Option(
        "--coding-ratio-threshold", help="Threshold ratio of coding sequence length to transcript length for rework CDS (default: 0.7).",
        rich_help_panel=CDS_PANEL,
    )] = 0.7,
    allow_internal_stops: Annotated[bool, typer.Option(
        "--allow-internal-stops/--no-allow-internal-stops", help="Allow internal stop codons in progressive rework fallback (default: True).",
        rich_help_panel=CDS_PANEL,
    )] = True,
    allow_partial: Annotated[bool, typer.Option(
        "--allow-partial/--no-allow-partial", help="Allow partial ORFs without stop codon in progressive rework fallback (default: True).",
        rich_help_panel=CDS_PANEL,
    )] = True,
    enforce_start_codon: Annotated[bool, typer.Option(
        "--enforce-start-codon/--no-enforce-start-codon", help="Require start codon (ATG) in initial rework passes (default: True).",
        rich_help_panel=CDS_PANEL,
    )] = True,
    orf_choice_mode: Annotated[str, typer.Option(
        "--orf-choice-mode", help="ORF selection criteria: 'longest' or 'earliest' (default: 'longest').",
        rich_help_panel=CDS_PANEL,
    )] = "longest",
    min_codon_len: Annotated[int, typer.Option(
        "--min-codon-len", help="Minimum codon length required for predicted ORFs (e.g. 30 or 50 to suppress micro-ORFs, default: 2).",
        rich_help_panel=CDS_PANEL,
    )] = 2,
    adjust_internal_shifts: Annotated[str, typer.Option(
        "--adjust-internal-shifts", help="Frameshift / phase handling mode: 'intra_exon' (default), 'all', or 'none'.",
        rich_help_panel=CDS_PANEL,
    )] = "intra_exon",
    skip_coordinate_polishing: Annotated[bool, typer.Option(
        "--skip-coordinate-polishing", help="Do not mutate feature coordinates when boundaries differ; log discrepancies as warnings instead.",
        rich_help_panel=CDS_PANEL,
    )] = False,
    recalculate_phases: Annotated[bool, typer.Option(
        "--recalculate-phases", help="Recalculate CDS segment phases based on segment lengths and splicing leftover, preserving 5' initial phase for partial CDSs.",
        rich_help_panel=CDS_PANEL,
    )] = False,
    reset_phases_zero: Annotated[bool, typer.Option(
        "--reset-phases-zero", help="Reset initial CDS phase to 0 and recalculate all downstream segment phases.",
        rich_help_panel=CDS_PANEL,
    )] = False,

    taxonomy: TaxonomyOption = "plant",
    genetic_code: GeneticCodeOption = 1,
    auto_organelle_codes: AutoOrganelleCodesOption = True,
    mito_code: MitoCodeOption = None,
    plastid_code: PlastidCodeOption = None,
    mitochondria_chroms: MitochondriaChromsOption = [],
    chloroplast_chroms: ChloroplastChromsOption = [],

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

    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
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


    if (rework_all_CDSs or infer_missing_CDSs) and not genome_file:
        raise typer.BadParameter("A genome FASTA file must be provided via --genome-file when using --rework-all-CDSs or --infer-missing-CDSs.")

    if genome_file:
        if not genome_name:
            genome_name = os.path.splitext(os.path.basename(genome_file))[0]
        genome = Genome(
            name=genome_name,
            genome_file_path=genome_file,
            quiet=quiet,
            header_id_tag=header_id_tag if header_id_tag != "" else None,
            header_id_regex=header_id_regex if header_id_regex != "" else None,
            gwh=gwh,
        )
    else:
        genome = None

    os.makedirs(output_dir, exist_ok=True)


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
        table=genetic_code,
        auto_organelle_codes=auto_organelle_codes,
        mito_table=mito_code,
        plastid_table=plastid_code,
        mitochondria_chroms=mitochondria_chroms,
        chloroplast_chroms=chloroplast_chroms,
        skip_coordinate_polishing=skip_coordinate_polishing,
        coding_ratio_threshold=coding_ratio_threshold,
        allow_internal_stops=allow_internal_stops,
        allow_partial=allow_partial,
        enforce_start_codon=enforce_start_codon,
        orf_choice_mode=orf_choice_mode,
        min_codon_len=min_codon_len,
        recalculate_phases=recalculate_phases,
        reset_phases_zero=reset_phases_zero,
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
