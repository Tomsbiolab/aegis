import typer
import os

from typing_extensions import Annotated

from ..annotation import Annotation, detect_file_format, read_file_with_fallback

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    input_format: Annotated[str, typer.Option(
        "-f", "--input-format", help="GTF/GFF format is automatically detected. Choose GTF or GFF to override."
    )] = "Auto Detect",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder."
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation filename, without extension."
    )] = "{annotation-name}.{ext}",
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
):
    """
    Convert between GFF and GTF formats.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_dir, exist_ok=True)

    if output_dir == "./aegis_output/":
        subfolder = True
    else:
        subfolder = False

    encoding = read_file_with_fallback(annotation_file)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet)

    input_format = input_format.lower()

    if input_format == "auto detect":
        input_format = detect_file_format(annotation_file, encoding=encoding)
    elif input_format.lower() == "gff" or input_format.lower() == "gff3":
        input_format = "gff3"
    elif input_format.lower() == "gtf":
        input_format = "gtf"
    else:
        raise ValueError(f"Invalid input format: {input_format}. Choose 'Auto Detect', 'GFF', 'GFF3' or 'GTF'.")

    if output_file == "{annotation-name}.{ext}":
        output_file = f"{annotation_name}"
        if input_format == "gff3":
            output_file += ".gtf"
        elif input_format == "gtf":
            output_file += ".gff3"

    if input_format == "gff3":
        annotation.export.gtf(custom_path=output_dir, tag=output_file, UTRs=True, quiet=quiet, subfolder=subfolder)
    elif input_format == "gtf":
        annotation.export.gff(custom_path=output_dir, tag=output_file, UTRs=True, quiet=quiet, subfolder=subfolder)

if __name__ == "__main__":
    app()
