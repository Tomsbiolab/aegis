import typer
import os

from typing_extensions import Annotated

from ..annotation import Annotation

RNA_CLASSES = ["mRNA", "antisense_lncRNA", "antisense_RNA", 
                "miRNA_primary_transcript", "ncRNA", "lncRNA",
                "lnc_RNA", "pseudogenic_tRNA", "rRNA", "snoRNA",
                "snRNA", "tRNA", "pre_miRNA", "tRNA_pseudogene",
                "SRP_RNA", "RNase_MRP_RNA"]

app = typer.Typer(add_completion=False, no_args_is_help=True)

def split_callback(value:str):
    if value:
        return [item.strip() for item in value.split(",")]
    return []

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
    features: Annotated[str, typer.Option(
        "-f", "--features", help=f"Selects only certain transcripts (e.g., 'mRNA,lncRNA'). Provide a comma-separated list. If empty, all biotypes are included. This option automatically enables 'clean_features'.",
        callback=split_callback
    )] = "",
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
    infer_missing_CDSs: Annotated[bool, typer.Option(
        "--infer-missing-CDSs", help="Detects and creates CDSs where missing, without overriding existing CDS annotations."
    )] = False,
    rework_all_CDSs: Annotated[bool, typer.Option(
        "--rework-all-CDSs", help="Recalculates ALL CDSs from the genome sequence, overriding existing ones. More aggressive than --infer-missing-CDSs."
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
    )] = False
    
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

    os.makedirs(output_dir, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, rework_all_CDSs=rework_all_CDSs, work_out_missing_CDSs=infer_missing_CDSs, quiet=quiet, collapse_exons=collapse_exons, collapse_CDSs=collapse_CDSs, consider_read_utrs=consider_read_utrs, rename_features=rename_features, standardise_features=standard_features, remove_genes_with_no_transcripts=remove_genes_with_no_transcripts, remove_transcripts_with_no_exons=remove_transcripts_with_no_exons, remove_missing_transcript_parent_references=remove_missing_transcript_parent_references, remove_genes_with_no_transcripts_even_if_pseudogene=remove_genes_with_no_transcripts, skip_orphaned_features=not(print_orphaned_features))

    if output_file == "{annotation-name}_tidy.gff3":
        output_file = f"{annotation_name}_tidy.gff3"

    if unique_cds_entry_ids or for_lifton:
        annotation.CDS_to_CDS_segment_ids()
    else:
        annotation.CDS_segment_to_CDS_ids()

    if symbols_as_description:
        remove_symbols = True

    if features:
        annotation.filter_by_rna_class(rna_classes=features)
        clean_features = True

    annotation.export.gff(custom_path=output_dir, tag=output_file, main_only=main_only, UTRs=include_UTRs, just_genes=just_genes, repeat_exons_utrs=repeat_exons_utrs, skip_atypical_fts=clean_features, quiet=quiet, aliases=(not remove_aliases), symbols=(not remove_symbols), symbols_as_description=symbols_as_description, clean_attributes=clean_attributes, featurecountsID=add_gene_id, print_empty_attributes=print_empty_attributes, subfolder=subfolder)

if __name__ == "__main__":
    app()
