import typer
import os

from typing_extensions import Annotated

from ..annotation import Annotation, detect_file_format, read_file_with_fallback
from .utils import IO_PANEL, EXEC_PANEL

FORMAT_PANEL = "Format Options"
COORDS_PANEL = "Coordinate & Phase Options"

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    input_format: Annotated[str, typer.Option(
        "-f", "--input-format", "--format", help="GTF/GFF format is automatically detected. Choose GTF or GFF to override.",
        rich_help_panel=FORMAT_PANEL,
    )] = "Auto Detect",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation filename, without extension.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name}.{ext}",
    strict_gtf_2_2: Annotated[bool, typer.Option(
        "--strict-gtf-2-2", help="Export strict GTF 2.2 format without top-level gene/transcript lines and with 5UTR/3UTR features.",
        rich_help_panel=FORMAT_PANEL,
    )] = False,
    polish_coordinates: Annotated[bool, typer.Option(
        "--polish-coordinates/--skip-coordinate-polishing", help="Mutate feature coordinates when boundaries differ (default: False, preserves original coordinates).",
        rich_help_panel=COORDS_PANEL,
    )] = False,
    recalculate_phases: Annotated[bool, typer.Option(
        "--recalculate-phases", help="Recalculate CDS segment phases based on segment lengths and splicing leftover, preserving 5' initial phase for partial CDSs.",
        rich_help_panel=COORDS_PANEL,
    )] = False,
    reset_phases_zero: Annotated[bool, typer.Option(
        "--reset-phases-zero", help="Reset initial CDS phase to 0 and recalculate all downstream segment phases.",
        rich_help_panel=COORDS_PANEL,
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
):
    """
    Convert between GFF and GTF formats.
    """

    skip_coordinate_polishing = not polish_coordinates

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_dir, exist_ok=True)

    if output_dir == "./aegis_output/":
        subfolder = True
    else:
        subfolder = False

    encoding = read_file_with_fallback(annotation_file)

    annotation = Annotation(
        name=annotation_name,
        annot_file_path=annotation_file,
        quiet=quiet,
        skip_coordinate_polishing=skip_coordinate_polishing,
        recalculate_phases=recalculate_phases,
        reset_phases_zero=reset_phases_zero,
    )

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
        annotation.export.gtf(output_dir=output_dir, filename=output_file, UTRs=True, quiet=quiet, subfolder=subfolder, strict_gtf_2_2=strict_gtf_2_2)
    elif input_format == "gtf":
        annotation.export.gff(output_dir=output_dir, filename=output_file, UTRs=True, quiet=quiet, subfolder=subfolder)

if __name__ == "__main__":
    app()
