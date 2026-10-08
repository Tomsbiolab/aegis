import typer
import os
from typing import Optional
from typing_extensions import Annotated

from ..annotation import Annotation, detect_file_format, read_file_with_fallback
from .utils import IO_PANEL, EXEC_PANEL, FORMAT_PANEL, COORDS_PANEL

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file (or provide via -a/--annotation)."
    )] = "",

    # 1. Format Options
    input_format: Annotated[str, typer.Option(
        "-f", "--input-format", "--format", help="Input format: 'Auto Detect', 'GFF', 'GFF3', or 'GTF'. Automatically detected by default.",
        rich_help_panel=FORMAT_PANEL,
    )] = "Auto Detect",
    to_format: Annotated[Optional[str], typer.Option(
        "-t", "--to", help="Explicit target output format: 'gff3' or 'gtf'. Optional; if omitted, format is automatically inverted based on input (GFF to GTF or GTF to GFF).",
        rich_help_panel=FORMAT_PANEL,
    )] = None,
    strict_gtf_2_2: Annotated[bool, typer.Option(
        "--strict-gtf-2-2", help="Export strict GTF 2.2 format without top-level gene/transcript lines and with 5UTR/3UTR features.",
        rich_help_panel=FORMAT_PANEL,
    )] = False,
    strip_utrs: Annotated[bool, typer.Option(
        "--strip-utrs", "--no-utrs", help="Strip UTRs from exported output (by default, UTRs are retained).",
        rich_help_panel=FORMAT_PANEL,
    )] = False,

    # 2. Coordinate & Phase Options
    polish_coordinates: Annotated[bool, typer.Option(
        "--polish-coordinates/--skip-coordinate-polishing", help="Mutate feature coordinates when boundaries differ [default: disabled; preserves original coordinates].",
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

    # 3. Input / Output Options
    annotation_file_opt: Annotated[str, typer.Option(
        "-a", "--annotation", "--annotations", "--annotation-file", "--annot", help="Path to input annotation GFF/GTF file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", "--annotation-names", "--annot-name", help="Annotation version, name or tag.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation filename, with or without extension.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name}.{ext}",

    # 4. Execution / Debugging
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help="Enable detailed console output.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
):
    """
    Convert between GFF and GTF formats.
    """
    if verbose:
        quiet = False

    annot_in = annotation_file_opt if annotation_file_opt else annotation_file
    if not annot_in:
        raise typer.BadParameter("Missing required annotation file. Provide as positional argument or via -a/--annotation.")
    annotation_file = annot_in

    skip_coordinate_polishing = not polish_coordinates
    export_utrs = not strip_utrs

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_dir, exist_ok=True)
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

    norm_in_format = input_format.lower().strip()
    if norm_in_format == "auto detect":
        detected_in_format = detect_file_format(annotation_file, encoding=encoding)
    elif norm_in_format in ("gff", "gff3"):
        detected_in_format = "gff3"
    elif norm_in_format == "gtf":
        detected_in_format = "gtf"
    else:
        raise typer.BadParameter(f"Invalid input format: {input_format}. Choose 'Auto Detect', 'GFF', 'GFF3' or 'GTF'.")

    # Determine target format
    if to_format is not None:
        norm_to = to_format.lower().strip()
        if norm_to in ("gff", "gff3"):
            target_format = "gff3"
        elif norm_to == "gtf":
            target_format = "gtf"
        else:
            raise typer.BadParameter(f"Invalid target format: '{to_format}'. Choose 'gff3' or 'gtf'.")
    else:
        # Invert based on detected input format
        target_format = "gtf" if detected_in_format == "gff3" else "gff3"

    if "{annotation-name}" in output_file:
        output_file = output_file.replace("{annotation-name}", annotation_name)
    if "{ext}" in output_file:
        output_file = output_file.replace("{ext}", target_format)
    elif not output_file.endswith(f".{target_format}") and not output_file.endswith(".gff") and not output_file.endswith(".gtf") and not output_file.endswith(".gff3"):
        output_file += f".{target_format}"

    if target_format == "gtf":
        annotation.export.gtf(output_dir=output_dir, filename=output_file, UTRs=export_utrs, quiet=quiet, subfolder=subfolder, strict_gtf_2_2=strict_gtf_2_2)
    elif target_format == "gff3":
        annotation.export.gff(output_dir=output_dir, filename=output_file, UTRs=export_utrs, quiet=quiet, subfolder=subfolder)

if __name__ == "__main__":
    app()
