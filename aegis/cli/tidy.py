import typer
import os
from typing_extensions import Annotated
from aegis.annotation import Annotation

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
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation filename, without extension."
    )] = "{annotation-name}_tidy.gff3",
    main_only: Annotated[bool, typer.Option(
        "-m", "--main", help="Whether to include only a main transcript and main CDS per gene."
    )] = False,
    include_UTRs: Annotated[bool, typer.Option(
        "-ut", "--include_UTRs", help="Include UTRs in output gff, normally these are not required by external tools as they can be deduced from exons and CDS features."
    )] = False,
    just_genes: Annotated[bool, typer.Option(
        "-g", "--just_genes", help="Whether to only include gene level features."
    )] = False,
    remove_symbols: Annotated[bool, typer.Option(
        "-s", "--remove_symbols", help="Removes symbol attributes from gff output."
    )] = False,
    remove_aliases: Annotated[bool, typer.Option(
        "-al", "--remove_aliases", help="Removes alias attributes from gff output."
    )] = False,
    clean_attributes: Annotated[bool, typer.Option(
        "-c", "--clean_attributes", help="Removes non-standard attributes from a gff, may help with external tool compatibility issues."
    )] = False,
    for_featurecounts: Annotated[bool, typer.Option(
        "-f", "--for_featurecounts", help="Creates a special attribute 'featurecounts_id' which has the gene-id as a value but is given to all gene features and subfeatures. This is useful for example when using featureCounts at the exon level but summarising counts at the gene-id level."
    )] = False,
    symbols_as_descriptors: Annotated[bool, typer.Option(
        "-sd", "--symbols_as_descriptors", help="Places gene symbols as 'Description=' attributes. Useful for JBrowse(2) display."
    )] = False,
    repeat_exons_utrs: Annotated[bool, typer.Option(
        "-r", "--repeat_exons_utrs", help="Creates individual exon/UTR entries with individual parental references for cases where a feature has more than one transcript level parent."
    )] = False,
    cds_segment_ids: Annotated[bool, typer.Option(
        "-u", "--unique_cds_entry_ids", help="CDS entries corresponding to a same protein in a gff by default share the same id. However since the default format is incompatible with some external tools, this flag will ensure each CDS entry (line) has a unique id."
    )] = False, 
    for_lifton: Annotated[bool, typer.Option(
        "-l", "--for_lifton", help="Ensures output has individual CDS entry ids (-u) and 1bp features are removed, both of these modifications are required for LifOn compatibility in its current version."
    )] = False,
    clean_features: Annotated[bool, typer.Option(
        "-cf", "--clean_features", help="Removes non-standard features from a gff, may help with external tool compatibility issues."
    )] = False,
    rna_classes: Annotated[str, typer.Option(
        "-rc", "--rna-classes", help=f"Filters out transcripts by biotype (e.g., 'mRNA,lncRNA'). Provide a comma-separated list. If empty, all biotypes are included. This option automatically enables 'clean_features'.",
        callback=split_callback
    )] = ""

):
    """
    Cleans and reformats a GFF/GTF file to correct common formatting errors and improve compatibility with other bioinformatics tools.
    
    This script parses an annotation file, allows for extensive filtering and reformatting, and exports a standardized GFF3 file.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    for rna_class in rna_classes:
        if rna_class not in RNA_CLASSES:
            raise typer.BadParameter(f"Invalid rna class: {rna_class}. Choose from: {RNA_CLASSES}")

    os.makedirs(output_folder, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file)

    if output_file == "{annotation-name}_tidy.gff3":
        output_file = f"{annotation_name}_tidy.gff3"

    if cds_segment_ids or for_lifton:
        annotation.CDS_to_CDS_segment_ids()
    else:
        annotation.CDS_segment_to_CDS_ids()

    if symbols_as_descriptors:
        remove_symbols = True

    if rna_classes:
        annotation.filter_by_rna_class(rna_classes=rna_classes)
        clean_features = True

    annotation.update_attributes(clean=clean_attributes, featurecountsID=for_featurecounts, symbols=(not remove_symbols), symbols_as_descriptors=symbols_as_descriptors, aliases=(not remove_aliases))

    annotation.export_gff(custom_path=output_folder, tag=output_file, main_only=main_only, UTRs=include_UTRs, just_genes=just_genes, repeat_exons_utrs=repeat_exons_utrs, no_1bp_features=for_lifton, skip_atypical_fts=clean_features)

if __name__ == "__main__":
    app()
