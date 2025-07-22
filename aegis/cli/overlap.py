import typer
import os
from typing_extensions import Annotated
from aegis.annotation import Annotation

app = typer.Typer(add_completion=False)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Option(
        "-a", "--annotation-file", help="Path to the input annotation GFF/GTF file."
    )],
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output TSV file."
    )] = "{annotation-name}_self_overlaps.tsv",
):
    """
    Finds overlapping features in an annotation file and exports them to a table.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(annotation_file)[0]

    if output_file == "{annotation-name}_self_overlaps.tsv":
        output_file = f"{annotation_name}_self_overlaps.tsv"

    os.system(f"mkdir -p {output_folder}")

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file)

    annotation.detect_gene_overlaps()

    annotation.export_equivalences(custom_path=output_folder, output_file=output_file, export_self=True, export_csv=True, return_df=False, NAs=False)

if __name__ == "__main__":
    app()
