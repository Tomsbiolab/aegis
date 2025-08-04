import typer
import os
from typing_extensions import Annotated
from aegis.annotation import Annotation

app = typer.Typer(add_completion=False)


def split_callback(value:str):
    if value:
        return [item.strip() for item in value.split(",")]
    return []

features: list = ["gene", "transcript", "CDS", "exon", "UTR"]

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
        "-o", "--output-file", help="Path to the output annotation file."
    )] = "{annotation-name}_renamed.gff3",
    feature_type: Annotated[str, typer.Option(
        "-f", "--feature_type", help=f"Choose what feature levels will have ids renamed, separated by commas. Choose from: {features}.",
        callback=split_callback
    )] = "transcript,CDS,exon,UTR",
    keep_ids_with_gene_id_contained: Annotated[bool, typer.Option(
        "-k", "--rename-minimal", help="Only rename a gene subfeature id if it does not include the parental 'gene_id' base. I.e. leave features such as 'gene_id_t001' untouched but rename 't001' as it does not contain the parental gene_id."
    )] = False,
    keep_numbering: Annotated[bool, typer.Option(
        "-n", "--keep-numbering", help="Try to retain original gene subfeature id numbering. I.e. rename a transcript id from 'gene_id_t004' to 'gene_id_T.4' without losing the original transcript number."
    )] = False,
    cds_segment_ids: Annotated[bool, typer.Option(
        "-u", "--unique_cds_entry_ids", help="CDS entries corresponding to a same protein in a gff by default share the same id. However since the default format is incompatible with some external tools, this flag will ensure each CDS entry (line) has a unique id."
    )] = False,
    prefix: Annotated[str, typer.Option(
        "-b", "--gene-id-prefix", help="Choose a new gene id prefix to rename the whole annotation. e.g. swich from 'VIT...' to 'Vitvi...'. Together with other options such as -sf, -sp, -se, and -gd the general feature id structure can be designed: i.e. '{prefix}{chromosome/scaffold}g{gene_count:0{gene_num_digits}d}{separator}{suffix}'. Gene subfeatures will be renamed on the basis of the configured parental gene-id."
    )] = "",
    suffix: Annotated[str, typer.Option(
        "-sf", "--suffix", help="Choose a new gene id suffix to include within the new gene id structure; See -b option."
    )] = "",
    spacer: Annotated[int, typer.Option(
        "-sp", "--spacer", help="Gene id number jump between one gene-id and the next, i.e. a spacer of 10 would result in '{prefix}{chromosome/scaffold}g00010' followed by '{prefix}{chromosome/scaffold}g00020'; See -b option."
    )] = 10,
    sep: Annotated[str, typer.Option(
        "-se", "--separator", help="Choose a new gene id suffix to include within the new gene id structure. e.g. '{prefix}{chromosome/scaffold}g{gene_count:0{gene_num_digits}d}.{suffix}' instead of '{prefix}{chromosome/scaffold}g{gene_count:0{gene_num_digits}d}_{suffix}'; See -b option."
    )] = "_",
    g_id_digits: Annotated[int, typer.Option(
        "-gd", "--gene-id-digits", help="Choose the number of digits to use for the gene id number. e.g. '{prefix}{chromosome/scaffold}g00010_{suffix}' would be the default first gene number in a particular chromosome or scaffold; See -b option."
    )] = 5,
    t_id_digits: Annotated[int, typer.Option(
        "-td", "--transcript-id-digits", help="Choose the number of digits to use for the transcript id number suffix. With the default three digits and the '_' separator the first transcript of every gene would have the '_t001' suffix."
    )] = 3,
    strip_gene_tag: Annotated[int, typer.Option(
        "-rt", "--remove-gene-tag", help="Some gffs have a flanking literal 'gene' tag. Use this flag to remove it. e.g. 'gene-Solyc00g174340' would become just 'Solyc00g174340'"
    )] = True,
    remove_point_suffix: Annotated[int, typer.Option(
        "-rp", "--remove-point-suffix", help="Some gene id formats carry an '.annotation-version' suffix that is in some cases not welcome. Use this flag to remove it: e.g. 'Solyc00g174340.2' would become just 'Solyc00g174340'."
    )] = True,
    gene_id_correspondences: Annotated[bool, typer.Option(
        "-c", "--gene-id-correspondences", help="Whether to produce a tsv file with correspondences between old and renamed gene ids '{annotation-name}_renamed_correspondences.tsv'."
    )] = False,
):
    """
    Rename feature ids of an annotation file.
    """

    if feature_type not in features:
        raise typer.BadParameter(f"Invalid feature leve: {feature_type}. Choose from: {features}")

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_folder, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file)

    if output_file == "{annotation-name}_renamed.gff3":
        output_file = f"{annotation_name}_renamed.gff3"

    annotation.rename_ids(custom_path=output_folder, features=features, keep_ids_with_gene_id_contained=keep_ids_with_gene_id_contained, remove_point_suffix=remove_point_suffix, strip_gene_tag=strip_gene_tag, keep_subfeature_numbers=keep_numbering, cds_segment_ids=cds_segment_ids, prefix=prefix, suffix=suffix, spacer=spacer, sep=sep, g_id_digits=g_id_digits, t_id_digits=t_id_digits, correspondences=gene_id_correspondences)

    annotation.export_gff(custom_path=output_folder, tag=output_file)

if __name__ == "__main__":
    app()
