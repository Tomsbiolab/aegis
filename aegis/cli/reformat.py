import typer
import os
from typing_extensions import Annotated
from aegis.annotation import Annotation, detect_file_format, read_file_with_fallback

app = typer.Typer(add_completion=False)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    input_format: Annotated[str, typer.Option(
        "-m", "--input_format", help="GTF/GFF format is automatically detected. Choose GTF or GFF to override."
    )] = "Auto Detect",
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation filename, without extension."
    )] = "{annotation-name}.{ext}"
):
    """
    Convert between GFF and GTF formats.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_folder, exist_ok=True)

    encoding = read_file_with_fallback(annotation_file)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file)

    input_format = input_format.lower()

    if input_format == "auto detect":
        input_format = detect_file_format(annotation_file, encoding=encoding)
    elif input_format.lower() == "gff":
        input_format = "gff3"

    if output_file == "{annotation-name}.{ext}":
        output_file = f"{annotation_name}"

    if input_format == "gff3":
        output_file += ".gtf"
        annotation.export_gtf(custom_path=output_folder, tag=output_file, UTRs=True)

    elif input_format == "gtf":
        output_file += ".gff3"
        annotation.export_gff(custom_path=output_folder, tag=output_file, UTRs=True)

if __name__ == "__main__":
    app()
