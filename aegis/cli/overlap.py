import typer
import os
import warnings

from typing import List, Optional
from typing_extensions import Annotated

from ..annotation import Annotation
from ..utils.genefunctions import export_group_equivalences
from .utils import split_callback, IO_PANEL, EXEC_PANEL, OVERLAP_PANEL

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_files: Annotated[List[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file(s) associated to the same genome assembly. Input only one to measure gene overlaps within a single annotation, input several to compare between annotation files."
    )] = [],

    # 1. Overlap Criteria
    overlap_threshold: Annotated[int, typer.Option(
        "-ot", "--overlap-threshold", help="Select the required overlap threshold to report a gene-id pair match (default: 6). Increase for more stringent comparisons, or decrease for more extensive reporting.",
        rich_help_panel=OVERLAP_PANEL,
    )] = 6,
    include_nas: Annotated[bool, typer.Option(
        "-na", "--include-nas", help="Whether to include unmapped / non-overlapping gene IDs in the output table.",
        rich_help_panel=OVERLAP_PANEL,
    )] = False,
    simple: Annotated[bool, typer.Option(
        "-s", "--simple", help="Whether to remove percentage overlap details at different feature levels for a simplified output table.",
        rich_help_panel=OVERLAP_PANEL,
    )] = False,

    # 2. Input / Output Options
    annotation_files_opt: Annotated[List[str], typer.Option(
        "-a", "--annotation", "--annotations", "--annotation-files", "--annotation-file", help="Path to input annotation GFF/GTF file(s). Overrides positional arguments if provided.",
        callback=split_callback,
        rich_help_panel=IO_PANEL,
    )] = [],
    annotation_names: Annotated[List[str], typer.Option(
        "-an", "--annotation-names", "--annotation-name", help="Annotation versions, names or tags. Provide in the same order as annotation files, separated by commas.",
        callback=split_callback,
        rich_help_panel=IO_PANEL,
    )] = ["{annotation-filename(s)}"],
    reference_annotation: Annotated[Optional[str], typer.Option(
        "-r", "--reference-annotation", help="Select a single annotation (by name or filename) to use as reference. Only matches to/from this annotation are reported.",
        rich_help_panel=IO_PANEL,
    )] = None,
    original_annotation_files: Annotated[List[str], typer.Option(
        "--original-annotation-files", help="Optional original annotation files before coordinate transfer / liftover to evaluate conservation of synteny.",
        callback=split_callback,
        rich_help_panel=IO_PANEL,
    )] = [],
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output directory.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_filetag: Annotated[str, typer.Option(
        "-o", "--output-file", "--output-filetag", "--filetag", help="Output file prefix/name.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name(s)}",

    # 3. Execution / Debugging
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help="Verbose logging, useful if encountering a problem or error.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
):
    """
    Calculates degree of gene overlaps between annotations associated to the same assembly and results in a gene-id equivalence table. If only one annotation file is provided as input, gene overlaps within the same annotation will be measured.
    """
    if verbose:
        quiet = False
    detailed_output = not simple

    annot_files = list(annotation_files_opt) if annotation_files_opt else list(annotation_files)
    if not annot_files:
        raise typer.BadParameter("At least one annotation file must be provided. Provide as positional argument or via -a/--annotations.")

    if len(annot_files) > 1 and annot_files[-1].lower() in ("true", "false"):
        typer.echo(
            "⚠️  Detected extra value 'true' or 'false' at the end of positional arguments.\n"
            "👉 Did you mean to use the '--include-nas' or '--simple' flags? Use them like this: '-na' or '-s' (no 'true' needed).",
            err=True,
        )
        raise typer.Exit(code=1)

    os.makedirs(output_dir, exist_ok=True)

    if annotation_names == ["{annotation-filename(s)}"]:
        annotation_names = []
        for annotation_file in annot_files:
            annotation_names.append(os.path.splitext(os.path.basename(annotation_file))[0])

    if len(annot_files) != len(annotation_names):
        raise typer.BadParameter(f"The provided number of annotation name(s)/tag(s) do not match the number of annotation file(s).")

    if len(annotation_names) != len(set(annotation_names)):
        raise typer.BadParameter("Avoid repeated annotation tag(s)/name(s).")
    
    if len(annot_files) != len(set(annot_files)):
        raise typer.BadParameter("Avoid repeated annotation filename(s).")

    if original_annotation_files != []:
        synteny = True
        if len(annot_files) != len(original_annotation_files):
            raise typer.BadParameter(f"The provided number of original annotation files do not match the number of annotation file(s).")
    else:
        synteny = False
        original_annotation_files = ["NA"] * len(annot_files)
    
    if reference_annotation is not None and reference_annotation != "None":
        if reference_annotation not in annot_files and reference_annotation not in annotation_names:
            raise typer.BadParameter(f"The provided reference-annotation = {reference_annotation} is not present neither in annotation-files ({annot_files}) nor annotation-names ({annotation_names}).")
        
    if len(annot_files) == 1:
        if original_annotation_files[0] != "NA":
            warnings.warn(f"Note that the provided original annotation file {original_annotation_files[0]} will not be used as synteny analysis is not implemented when evaluating gene overlaps within a single annotation = {annotation_names[0]}.", category=UserWarning)

    annotations = []

    for n, annotation_file in enumerate(annot_files):
        if original_annotation_files[n].lower() != "na":
            original_annotation = Annotation(name=f"{annotation_names[n]}_original", annot_file_path=original_annotation_files[n], quiet=quiet, skip_coordinate_polishing=True)
            annotations.append(Annotation(name=annotation_names[n], annot_file_path=annotation_file, original_annotation=original_annotation, quiet=quiet, skip_coordinate_polishing=True))
        else:
            annotations.append(Annotation(name=annotation_names[n], annot_file_path=annotation_file, quiet=quiet, skip_coordinate_polishing=True))

        if annotation_names[n] == reference_annotation or annotation_file == reference_annotation:
            annotations[n].target = True

    if len(annot_files) == 1:
        if output_filetag == "{annotation-name(s)}":
            output_file = annotations[0].name
        else:
            output_file = output_filetag
            
        output_file += f"_self_overlaps_t{overlap_threshold}.csv"

        annotations[0].overlaps.detect()
        _ = annotations[0].overlaps.export(output_dir=output_dir, filename=output_file, verbose=detailed_output, overlap_threshold=overlap_threshold, export_self=True, save_csv=True, NAs=include_nas, quiet=quiet)

    elif len(annot_files) > 1:
        if output_filetag == "{annotation-name(s)}":
            output_filetag = ""
            
        export_group_equivalences(annotations, output_folder=output_dir, verbose=detailed_output, synteny=synteny, group_tag=output_filetag, overlap_threshold=overlap_threshold, include_NAs=include_nas, output_also_single_files=False, quiet=quiet)


if __name__ == "__main__":
    app()
