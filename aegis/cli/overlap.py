import typer
import os
import warnings
from typing import List
from typing_extensions import Annotated
from aegis.annotation import Annotation
from aegis.genefunctions import export_group_equivalences

app = typer.Typer(add_completion=False)

def split_callback(value:str):
    if value:
        return [item.strip() for item in value.split(",")]
    return []

@app.command()
def main(
    annotation_files: Annotated[List[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file(s) associated to the same genome assembly. Input only one to measure gene overlaps within a single annotation, input several to compare between annotation files."
    )],
    annotation_names: Annotated[str, typer.Option(
        "-a", "--annotation-names", help="Annotation versions, names or tags. Provide them in the same number and order as the corresponding annotation files, separated by commas. e.g. name1,name2",
        callback=split_callback
    )] = "{annotation-filename(s)}",
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/",
    output_filetag: Annotated[str, typer.Option(
        "-f", "--output-filetag", help="Optional output filetag prefix to prevent auto-naming based on annotation names, specially useful when comparing several annotations."
    )] = "{annotation-name(s)}",
    overlap_threshold: Annotated[int, typer.Option(
        "-o", "--overlap-threshold", help="Select the required overlap threshold to report a gene-id pair match. The default value of 6 is expected to result in a valid set of id equivalences between annotation files. Increase it for more stringent comparisons, or decrease it for more extensive reporting of overlaps."
    )] = 6,
    include_NAs: Annotated[bool, typer.Option(
        "-n", "--include-NAs", help="Whether to include NAs in output file, i.e. whether gene ids without overlaps are listed or not."
    )] = False,
    simple: Annotated[bool, typer.Option(
        "-s", "--simple", help="Whether to remove percentage overlap details at different feature levels for a more simple output table."
    )] = False,
    original_annotation_files: Annotated[str, typer.Option(
        "-t", "--original-annotation-files", help="Should some of the annotations be a result of a liftover or coordinate transfer, you can optionally provide a list of the original files before the transfer, separated by commas. If at least 2 annotation files are being compared, conservation of synteny will be calculated wherever possible based on gene order before/after transfer. These original annotation files must be in the same number and order as the corresponding annotation files. Use NA as a placemarker for annotation files without an original annotation file. e.g. '-t original_file_1,NA,original_file_3'",
        callback=split_callback
    )] = "",
    reference_annotation: Annotated[str, typer.Option(
        "-r", "--reference-annotation", help="Select a single annotation, by providing its name/tag or filename, to use as a reference. Only matches to and from this annotation will be reported. Otherwise matches are reported between all annotations."
    )] = "None",
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

    if annotation_names == "{annotation-filename(s)}":
        annotation_names = []
        for annotation_file in annotation_files:
            annotation_names.append(os.path.splitext(os.path.basename(annotation_file))[0])

    if len(annotation_files) != len(annotation_names):
        raise ValueError(f"The provided number of annotation name(s)/tag(s) do not match the number of annotation file(s).")

    if len(annotation_names) != len(set(annotation_names)):
        raise ValueError("Avoid repeated annotation tag(s)/name(s).")
    
    if len(annotation_files) != len(set(annotation_files)):
        raise ValueError("Avoid repeated annotation filename(s).")

    if original_annotation_files != []:
        synteny = True
        if len(annotation_files) != len(original_annotation_files):
            raise ValueError(f"The provided number of original annotation files do not match the number of annotation file(s).")
        
    else:
        synteny = False
        original_annotation_files = ["NA"] * len(annotation_files)
    
    if reference_annotation != "None":
        if reference_annotation not in annotation_files and reference_annotation not in annotation_names:
            raise ValueError(f"The provided reference-annotation = {reference_annotation} is not present neither in annotation-files ({annotation_files}) nor annotation-names ({annotation_names}).")
        
    if len(annotation_files) == 1:
        if original_annotation_files[0] != "NA":
            warnings.warn(f"Note that he provided original annotation file {original_annotation_files[0]} will not be used as synteny analysis is not implemented when evaluating gene overlaps within a single annotation = {annotation_names[0]}.", category=UserWarning)

    annotations = []

    for n, annotation_file in enumerate(annotation_files):

        if original_annotation_files[n].lower() != "na":
            original_annotation = Annotation(name=f"{annotation_names[n]}_original", annot_file_path=original_annotation_files[n])
            annotations.append(Annotation(name=annotation_names[n], annot_file_path=annotation_file, original_annotation=original_annotation))
        else:
            annotations.append(Annotation(name=annotation_names[n], annot_file_path=annotation_file))

        if annotation_names[n] == reference_annotation or annotation_file == reference_annotation:
            annotations[n].target = True

    if len(annotation_files) == 1:

        if output_filetag == "{annotation-name(s)}":
            output_file = annotations[0].name
        else:
            output_file = output_filetag
            
        output_file += f"_self_overlaps_t{overlap_threshold}.csv"

        annotations[0].detect_gene_overlaps()

        annotations[0].export_equivalences(custom_path=output_folder, output_file=output_file, verbose=verbose, overlap_threshold=overlap_threshold, export_self=True, export_csv=True, return_df=False, NAs=include_NAs)

    elif len(annotation_files) > 1:

        if output_filetag == "{annotation-name(s)}":
            output_filetag = ""
            
        export_group_equivalences(annotations, output_folder=output_folder, verbose=verbose, synteny=synteny, group_tag=output_filetag, overlap_threshold=overlap_threshold, include_NAs=include_NAs, output_also_single_files=False)

    else:
        raise ValueError(f"No annotation-files provided.")


if __name__ == "__main__":
    app()
