import typer
import os

from typing_extensions import Annotated

from ..annotation import Annotation, read_file_with_fallback
from .utils import IO_PANEL, EXEC_PANEL

features = ["gene", "transcript"]
PRUNE_PANEL = "Pruning Options"

app = typer.Typer(add_completion=False, no_args_is_help=True)
@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    target_ids: Annotated[str, typer.Argument(
        help="Input file with list of IDs (one per line) OR a comma-separated list of IDs to remove."
    )],

    # 1. Pruning Options
    feature_type: Annotated[str, typer.Option(
        "-f", "--feature-type", "--feature", help=f"Feature level to be removed based on input IDs. Choose from {features}.",
        rich_help_panel=PRUNE_PANEL,
    )] = "gene",

    # 2. Input / Output Options
    annotation_name: Annotated[str, typer.Option(
        "-a", "-an", "--annotation-name", help="Annotation version, name or tag.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output directory.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation filename, with or without extension.",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-name}_pruned",

    # 3. Execution / Debugging
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
    Remove a list of IDs (from a file or comma-separated argument) from an annotation.
    """
    if verbose:
        quiet = False

    if feature_type not in features:
        raise typer.BadParameter(f"Invalid feature level: {feature_type}. Choose from: {features}")

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_dir, exist_ok=True)
    subfolder = (output_dir == "./aegis_output/")

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, quiet=quiet, skip_coordinate_polishing=True)

    input_ids = set()
    if os.path.isfile(target_ids):
        encoding = read_file_with_fallback(target_ids)
        with open(target_ids, encoding=encoding) as f_in:
            for line in f_in:
                if line.startswith("#"):
                    continue
                clean_id = line.strip()
                if clean_id:
                    input_ids.add(clean_id)
    else:
        for item in target_ids.split(","):
            clean_id = item.strip()
            if clean_id:
                input_ids.add(clean_id)

    if not input_ids:
        raise typer.BadParameter("No valid IDs provided to prune.")

    if feature_type == "gene":
        annotation.remove_genes(to_remove=input_ids, override_rescue=True, quiet=quiet)
    else:
        annotation.remove_transcripts(to_remove=input_ids, remove_genes_accordingly=True, quiet=quiet)

    if not (output_file.endswith(".gff3") or output_file.endswith(".gff")):
        output_file += ".gff3"

    annotation.export.gff(output_dir=output_dir, filename=output_file, quiet=quiet, subfolder=subfolder)


if __name__ == "__main__":
    app()
