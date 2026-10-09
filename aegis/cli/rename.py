import typer
import os
from typing import Optional
from typing_extensions import Annotated

from ..annotation import Annotation
from .utils import split_callback, RENAME_PANEL, SUBFEATURE_PANEL, IO_PANEL, EXEC_PANEL

app = typer.Typer(add_completion=False, no_args_is_help=True)

VALID_FEATURES: list[str] = ["gene", "transcript", "CDS", "exon", "UTR"]

@app.command()
def main(
    annotation_file: Annotated[Optional[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file (or provide via -a/--annotation)."
    )] = None,

    # 1. Feature ID Renaming
    rename_features: Annotated[list[str], typer.Option(
        "-f", "--features", "--rename-features", help=f"Choose what feature levels will have ids renamed, separated by commas. Choose from: {VALID_FEATURES}.",
        callback=split_callback,
        rich_help_panel=RENAME_PANEL,
    )] = ["transcript", "CDS", "exon", "UTR"],
    prefix: Annotated[str, typer.Option(
        "-p", "--prefix", help="Choose a new gene id prefix to rename the whole annotation. e.g. switch from 'VIT...' to 'Vitvi...'. Together with other options such as --suffix, --spacer, --separator, and --gene-id-digits the general feature id structure can be designed: i.e. '{prefix}{chromosome/scaffold}g{gene_count:0{gene_num_digits}d}{separator}{suffix}'. Gene subfeatures will be renamed on the basis of the configured parental gene-id.",
        rich_help_panel=RENAME_PANEL,
    )] = "",
    suffix: Annotated[str, typer.Option(
        "--suffix", help="Choose a new gene id suffix to include within the new gene id structure; See --prefix option.",
        rich_help_panel=RENAME_PANEL,
    )] = "",
    spacer: Annotated[int, typer.Option(
        "--spacer", help="Gene id number jump between one gene-id and the next, i.e. a spacer of 10 would result in '{prefix}{chromosome/scaffold}g00010' followed by '{prefix}{chromosome/scaffold}g00020'; See --prefix option.",
        rich_help_panel=RENAME_PANEL,
    )] = 10,
    sep: Annotated[str, typer.Option(
        "--sep", "--separator", help="Choose a new gene id separator to include within the new gene id structure. e.g. '{prefix}{chromosome/scaffold}g{gene_count:0{gene_num_digits}d}.{suffix}' instead of '{prefix}{chromosome/scaffold}g{gene_count:0{gene_num_digits}d}_{suffix}'; See --prefix option.",
        rich_help_panel=RENAME_PANEL,
    )] = "_",
    g_id_digits: Annotated[int, typer.Option(
        "--gene-id-digits", help="Choose the number of digits to use for the gene id number. e.g. '{prefix}{chromosome/scaffold}g00010_{suffix}' would be the default first gene number in a particular chromosome or scaffold; See --prefix option.",
        rich_help_panel=RENAME_PANEL,
    )] = 5,
    t_id_digits: Annotated[int, typer.Option(
        "--transcript-id-digits", help="Choose the number of digits to use for the transcript id number suffix. With the default three digits and the '_' separator the first transcript of every gene would have the '_t001' suffix.",
        rich_help_panel=RENAME_PANEL,
    )] = 3,
    strip_gene_tag: Annotated[bool, typer.Option(
        "--strip-gene-tags", "--strip-gene-tag", "--remove-gene-tag", "--remove-gene-tags", help="Some GFFs have a flanking literal 'gene-' or 'gene:' prefix. Use this flag to remove it. e.g. 'gene-Solyc00g174340' becomes 'Solyc00g174340'.",
        rich_help_panel=RENAME_PANEL,
    )] = False,
    remove_point_suffix: Annotated[bool, typer.Option(
        "--remove-point-suffix", help="Some gene id formats carry an '.annotation-version' suffix that is in some cases not welcome. Use this flag to remove it: e.g. 'Solyc00g174340.2' becomes 'Solyc00g174340'.",
        rich_help_panel=RENAME_PANEL,
    )] = False,
    gene_id_correspondences: Annotated[bool, typer.Option(
        "--gene-id-correspondences", help="Whether to produce a TSV file with correspondences between old and renamed gene ids '{annotation-name}_renamed_correspondences.tsv'.",
        rich_help_panel=RENAME_PANEL,
    )] = False,

    # 2. Subfeature & Model Options
    keep_existing_ids_if_derived_from_base_id: Annotated[bool, typer.Option(
        "--rename-minimal", help="Only rename a gene subfeature id if it does not include the parental 'gene_id' base. I.e. leave features such as 'gene_id_t001' untouched but rename 't001' as it does not contain the parental gene_id.",
        rich_help_panel=SUBFEATURE_PANEL,
    )] = False,
    keep_numbering: Annotated[bool, typer.Option(
        "--keep-numbering", help="Try to retain original gene subfeature id numbering. I.e. rename a transcript id from 'gene_id_t004' to 'gene_id_T.4' without losing the original transcript number.",
        rich_help_panel=SUBFEATURE_PANEL,
    )] = False,
    unique_cds_entry_ids: Annotated[bool, typer.Option(
        "--unique-cds-entry-ids", help="CDS entries corresponding to a same protein in a gff by default share the same id. However since the default format is incompatible with some external tools, this flag will ensure each CDS entry (line) has a unique id.",
        rich_help_panel=SUBFEATURE_PANEL,
    )] = False,
    no_collapse_exons: Annotated[bool, typer.Option(
        "--no-collapse-exons", help="Do not merge overlapping/adjacent exons.",
        rich_help_panel=SUBFEATURE_PANEL,
    )] = False,
    no_collapse_cds: Annotated[bool, typer.Option(
        "--no-collapse-cds", help="Do not merge overlapping/adjacent CDS segments.",
        rich_help_panel=SUBFEATURE_PANEL,
    )] = False,

    # 3. Input / Output Options
    annotation_file_opt: Annotated[str, typer.Option(
        "-a", "--annotation", "--annotations", "--annotation-file", "--annot", help="Path to input annotation GFF/GTF file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", "--annotation-names", help="Annotation version, name or tag.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output directory.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation file.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name}_renamed.gff3",

    # 4. Execution & Debugging
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help="Increase terminal reporting verbosity.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
):
    """
    Rename feature ids of an annotation file.
    """

    if verbose:
        quiet = False

    annot_in = annotation_file_opt if annotation_file_opt else annotation_file
    if not annot_in:
        raise typer.BadParameter("Missing required annotation file. Provide as positional argument or via -a/--annotation.")
    annotation_file = annot_in

    collapse_exons = not no_collapse_exons
    collapse_CDSs = not no_collapse_cds

    for feature in rename_features:
        if feature not in VALID_FEATURES:
            raise typer.BadParameter(f"Invalid feature level: {feature}. Choose from: {VALID_FEATURES}")

    subfolder = False
        
    if (remove_point_suffix or strip_gene_tag) and "gene" not in rename_features:
        typer.echo(f"'gene' was not included in features={rename_features} but --remove_point_suffix or --strip_gene_tag flags were used, therefore, 'gene' was added to the list of modified features.", err=True)
        rename_features.append("gene") #type: ignore

    if rename_features == []:
        raise typer.BadParameter(f"No features were chosen to rename their ids. Select from: {VALID_FEATURES}.")

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_dir, exist_ok=True)

    annotation = Annotation(
        name=annotation_name,
        annot_file_path=annotation_file,
        quiet=quiet,
        collapse_exons=collapse_exons,
        collapse_CDSs=collapse_CDSs,
        skip_coordinate_polishing=True,
    )

    if "{annotation-name}" in output_file:
        output_file = output_file.replace("{annotation-name}", annotation_name)
    name_only = os.path.splitext(output_file)[0]
    output_file_correspondences = f"{name_only}_correspondences.tsv"

    annotation.rename_ids(correspondences_output_dir=output_dir, correspondences_filename=output_file_correspondences, features=tuple(rename_features), keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id, remove_point_suffix=remove_point_suffix, strip_gene_tag=strip_gene_tag, keep_subfeature_numbers=keep_numbering, cds_segment_ids=unique_cds_entry_ids, prefix=prefix, suffix=suffix, spacer=spacer, sep=sep, g_id_digits=g_id_digits, t_id_digits=t_id_digits, correspondences=gene_id_correspondences, quiet=quiet)
    annotation.export.gff(output_dir=output_dir, filename=output_file, quiet=quiet, subfolder=subfolder)

if __name__ == "__main__":
    app()
