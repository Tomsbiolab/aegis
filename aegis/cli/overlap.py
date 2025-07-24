import typer
import os
from typing import List
from typing_extensions import Annotated
from aegis.annotation import Annotation
from aegis.genefunctions import export_group_equivalences

app = typer.Typer(add_completion=False)

@app.command()
def main(
    annotation_files: Annotated[List[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file(s) associated to the same genome assembly. Input only one to measure gene overlaps within a single annotation, input several to compare between annotation files."
    )],
    annotation_names: Annotated[List[str], typer.Option(
        "-a", "--annotation-names", help="Annotation versions, names or tags. Provide them in the same number and order as the corresponding annotation files."
    )] = ["{annotation-filename(s)}"],
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/",
    overlap_threshold: Annotated[int, typer.Option(
        "-o", "--overlap_threshold", help="Select the required overlap threshold to report a gene-id pair match. The default value of 6 is expected to result in a valid set of id equivalences between annotation files. Increase it for more stringent comparisons, or decrease it for more extensive reporting of overlaps."
    )] = 6,
    include_NAs: Annotated[bool, typer.Option(
        "-n", "--include_NAs", help="Whether to include NAs in output file, i.e. whether gene ids without overlaps are listed or not."
    )] = False,
    simple: Annotated[bool, typer.Option(
        "-s", "--simple", help="Whether to remove percentage overlap details at different feature levels for a more simple output table."
    )] = False,
    original_annotation_files: Annotated[List[str], typer.Option(
        "-t", "--original-annotation-files", help="Should some of the annotations be a result of a liftover or coordinate transfer, you can optionally provide a list of the original files before the transfer. If at least 2 annotation files are being compared, conservation of synteny will be calculated wherever possible based on gene order before/after transfer. These original annotation files must be in the same number and order as the corresponding annotation files. Use NA as a placemarker for annotation files without an original annotation file. e.g. '-t original_file_1 NA original_file_3'"
    )] = [],
):
    """
    Calculates degree of gene overlaps between annotations associated to the same assembly and results in a gene-id equivalence table. If only one annotation file is provided as input, gene overlaps within the same annotation will be measured.
    """

    verbose = not simple

    if len(annotation_files) > 1 and annotation_files[-1].lower() in ("true", "false"):
        typer.echo(
            "⚠️  Detected extra value 'true' or 'false' at the end of positional arguments.\n"
            "👉 Did you mean to use the '--include_NAs' or '--simple' flags? Use them like this: '-n' or '-s' (no 'true' needed).",
            err=True,
        )
        raise typer.Exit(code=1)

    os.makedirs(output_folder, exist_ok=True)

    if annotation_names == ["{annotation-filename(s)}"]:
        annotation_names = []
        for annotation_file in annotation_files:
            annotation_names.append(os.path.splitext(os.path.basename(annotation_file))[0])

    if len(annotation_files) != len(annotation_names):
        raise ValueError(f"The provided number of annotation name(s)/tag(s) do not match the number of annotation file(s).")

    if original_annotation_files != []:
        if len(annotation_files) != len(original_annotation_files):
            raise ValueError(f"The provided number of original annotation files do not match the number of annotation file(s).")
        
    else:
        for annotation_file in annotation_files:
            original_annotation_files.append("NA")
    
    annotations = []

    for n, annotation_file in enumerate(annotation_files):

        if original_annotation_files[n].lower() != "na":
            original_annotation = Annotation(name=f"{annotation_names[n]}_original", annot_file_path=original_annotation_files[n])
            annotations.append(Annotation(name=annotation_names[n], annot_file_path=annotation_file, original_annotation=original_annotation))
        else:
            annotations.append(Annotation(name=annotation_names[n], annot_file_path=annotation_file))

    if len(annotation_files) == 1:

        output_file = f"{annotations[0].name}_self_overlaps_t{overlap_threshold}.tsv"

        annotations[0].detect_gene_overlaps()

        annotations[0].export_equivalences(custom_path=output_folder, output_file=output_file, verbose=verbose, overlap_threshold=overlap_threshold, export_self=True, export_csv=True, return_df=False, NAs=include_NAs)

    elif len(annotation_files) > 1:

        if len(annotation_files) == 2:

            output_tag = f"{annotations[0].name}_{annotations[1].name}"

        else:
            output_tag = f"{annotations[0].name}_..._{annotations[-1].name}"

        export_group_equivalences(annotations, custom_path=output_folder, verbose=verbose, multidirectional=True, group_tag=output_tag, overlap_threshold=overlap_threshold)
    else:
        raise ValueError(f"No annotation-files provided.")


if __name__ == "__main__":
    app()
