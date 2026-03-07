import typer
import os

from typing_extensions import Annotated

from ..annotation import Annotation

app = typer.Typer(add_completion=False)


def split_callback(value:str) -> list[str]:
    if value:
        return [item.strip() for item in value.split(",")]
    return []

VALID_FEATURES: list[str] = ["gene", "transcript", "CDS", "exon", "UTR"]

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
        "-o", "--output-file", help="Path to the output annotation file."
    )] = "{annotation-name}_renamed.gff3",
    rename_features: Annotated[str, typer.Option(
        "-f", "--features", help=f"Choose what feature levels will have ids renamed, separated by commas. Choose from: {VALID_FEATURES}.",
        callback=split_callback
    )] = "transcript,CDS,exon,UTR",
    keep_ids_with_gene_id_contained: Annotated[bool, typer.Option(
        "--rename-minimal", help="Only rename a gene subfeature id if it does not include the parental 'gene_id' base. I.e. leave features such as 'gene_id_t001' untouched but rename 't001' as it does not contain the parental gene_id."
    )] = False,
    keep_numbering: Annotated[bool, typer.Option(
        "--keep-numbering", help="Try to retain original gene subfeature id numbering. I.e. rename a transcript id from 'gene_id_t004' to 'gene_id_T.4' without losing the original transcript number."
    )] = False,
    unique_cds_entry_ids: Annotated[bool, typer.Option(
        "--unique-cds-entry-ids", help="CDS entries corresponding to a same protein in a gff by default share the same id. However since the default format is incompatible with some external tools, this flag will ensure each CDS entry (line) has a unique id."
    )] = False, 
    prefix: Annotated[str, typer.Option(
        "--prefix", help="Choose a new gene id prefix to rename the whole annotation. e.g. swich from 'VIT...' to 'Vitvi...'. Together with other options such as -sf, -sp, -se, and -gd the general feature id structure can be designed: i.e. '{prefix}{chromosome/scaffold}g{gene_count:0{gene_num_digits}d}{separator}{suffix}'. Gene subfeatures will be renamed on the basis of the configured parental gene-id."
    )] = "",
    suffix: Annotated[str, typer.Option(
        "--suffix", help="Choose a new gene id suffix to include within the new gene id structure; See -b option."
    )] = "",
    spacer: Annotated[int, typer.Option(
        "--spacer", help="Gene id number jump between one gene-id and the next, i.e. a spacer of 10 would result in '{prefix}{chromosome/scaffold}g00010' followed by '{prefix}{chromosome/scaffold}g00020'; See -b option."
    )] = 10,
    sep: Annotated[str, typer.Option(
        "--separator", help="Choose a new gene id suffix to include within the new gene id structure. e.g. '{prefix}{chromosome/scaffold}g{gene_count:0{gene_num_digits}d}.{suffix}' instead of '{prefix}{chromosome/scaffold}g{gene_count:0{gene_num_digits}d}_{suffix}'; See -b option."
    )] = "_",
    g_id_digits: Annotated[int, typer.Option(
        "--gene-id-digits", help="Choose the number of digits to use for the gene id number. e.g. '{prefix}{chromosome/scaffold}g00010_{suffix}' would be the default first gene number in a particular chromosome or scaffold; See -b option."
    )] = 5,
    t_id_digits: Annotated[int, typer.Option(
        "--transcript-id-digits", help="Choose the number of digits to use for the transcript id number suffix. With the default three digits and the '_' separator the first transcript of every gene would have the '_t001' suffix."
    )] = 3,
    strip_gene_tag: Annotated[bool, typer.Option(
        "--remove-gene-tag", help="Some gffs have a flanking literal 'gene' tag. Use this flag to remove it. e.g. 'gene-Solyc00g174340' would become just 'Solyc00g174340'"
    )] = False,
    remove_point_suffix: Annotated[bool, typer.Option(
        "--remove-point-suffix", help="Some gene id formats carry an '.annotation-version' suffix that is in some cases not welcome. Use this flag to remove it: e.g. 'Solyc00g174340.2' would become just 'Solyc00g174340'."
    )] = False,
    gene_id_correspondences: Annotated[bool, typer.Option(
        "--gene-id-correspondences", help="Whether to produce a tsv file with correspondences between old and renamed gene ids '{annotation-name}_renamed_correspondences.tsv'."
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
    no_collapse_exons: Annotated[bool, typer.Option(
        "--no-collapse-exons", help="Do not merge overlapping/adjacent exons."
    )] = False,
    no_collapse_CDSs: Annotated[bool, typer.Option(
        "--no-collapse-CDSs", help="Do not merge overlapping/adjacent CDS segments."
    )] = False
):
    """
    Rename feature ids of an annotation file.
    """
    collapse_exons = not no_collapse_exons
    collapse_CDSs = not no_collapse_CDSs

    for feature in rename_features:
        if feature not in VALID_FEATURES:
            raise typer.BadParameter(f"Invalid feature level: {feature}. Choose from: {VALID_FEATURES}")
        
    if (remove_point_suffix or strip_gene_tag) and "gene" not in rename_features:
        typer.echo(f"'gene' was not included in features={rename_features} but --remove_point_suffix or --strip_gene_tag flags were used, therefore, 'gene' was added to the list of modified features.", err=True)
        rename_features.append("gene")

    if rename_features == []:
        raise typer.BadParameter(f"No features were chosen to rename their ids. Select from: {VALID_FEATURES}.")

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_dir, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet, collapse_exons=collapse_exons, collapse_CDSs=collapse_CDSs)

    if output_file == "{annotation-name}_renamed.gff3":
        output_file = f"{annotation_name}_renamed.gff3"

    annotation.rename_ids(custom_path=output_dir, features=rename_features, keep_ids_with_gene_id_contained=keep_ids_with_gene_id_contained, remove_point_suffix=remove_point_suffix, strip_gene_tag=strip_gene_tag, keep_subfeature_numbers=keep_numbering, cds_segment_ids=unique_cds_entry_ids, prefix=prefix, suffix=suffix, spacer=spacer, sep=sep, g_id_digits=g_id_digits, t_id_digits=t_id_digits, correspondences=gene_id_correspondences, collapse_exons=collapse_exons, collapse_CDSs=collapse_CDSs, quiet=quiet)
    annotation.export.gff(custom_path=output_dir, tag=output_file, quiet=quiet)

if __name__ == "__main__":
    app()
