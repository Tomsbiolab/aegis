import typer
import os
from typing import Optional
from typing_extensions import Annotated

from ..annotation import Annotation, read_file_with_fallback
from .utils import IO_PANEL, EXEC_PANEL, PRUNE_PANEL

features = ["gene", "transcript"]

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[Optional[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file (or provide via -a/--annotation)."
    )] = None,
    target_ids: Annotated[Optional[str], typer.Argument(
        help="Input file with list of IDs (one per line) OR a comma-separated list of IDs (or provide via -i/--ids)."
    )] = None,

    # 1. Pruning Options
    feature_type: Annotated[str, typer.Option(
        "-f", "--feature-type", "--feature", help=f"Feature level to prune based on input IDs. Choose from {features}.",
        rich_help_panel=PRUNE_PANEL,
    )] = "gene",
    keep: Annotated[bool, typer.Option(
        "-k", "--keep", "--whitelist", "--invert", help="Invert selection: retain only the specified IDs and remove all others.",
        rich_help_panel=PRUNE_PANEL,
    )] = False,

    # 2. Input / Output Options
    annotation_file_opt: Annotated[str, typer.Option(
        "-a", "--annotation", "--annotations", "--annotation-file", "--annot", help="Path to input annotation GFF/GTF file. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    target_ids_opt: Annotated[str, typer.Option(
        "-i", "--ids", "--target-ids", help="Input file with list of IDs OR comma-separated list of IDs. Overrides positional argument if provided.",
        rich_help_panel=IO_PANEL,
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-an", "--annotation-name", "--annotation-names", "--annot-name", help="Annotation version, name or tag.",
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
    Remove (or retain with --keep) a list of IDs (from a file or comma-separated argument) in an annotation.
    """
    if verbose:
        quiet = False

    annot_in = annotation_file_opt if annotation_file_opt else annotation_file
    ids_in = target_ids_opt if target_ids_opt else target_ids

    if not annot_in:
        raise typer.BadParameter("Missing required annotation file. Provide as positional argument or via -a/--annotation.")
    if not ids_in:
        raise typer.BadParameter("Missing required target IDs. Provide as positional argument or via -i/--ids.")

    annotation_file = annot_in
    target_ids = ids_in

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

    if keep:
        if feature_type == "gene":
            to_remove = set(annotation.genes.keys()) - input_ids
            annotation.remove_genes(to_remove=to_remove, override_rescue=True, quiet=quiet)
        else:
            to_remove = set(annotation.transcripts.keys()) - input_ids
            annotation.remove_transcripts(to_remove=to_remove, remove_genes_accordingly=True, quiet=quiet)
    else:
        if feature_type == "gene":
            annotation.remove_genes(to_remove=input_ids, override_rescue=True, quiet=quiet)
        else:
            annotation.remove_transcripts(to_remove=input_ids, remove_genes_accordingly=True, quiet=quiet)

    if "{annotation-name}" in output_file:
        output_file = output_file.replace("{annotation-name}", annotation_name)

    if not (output_file.endswith(".gff3") or output_file.endswith(".gff")):
        output_file += ".gff3"

    annotation.export.gff(output_dir=output_dir, filename=output_file, quiet=quiet, subfolder=subfolder)


if __name__ == "__main__":
    app()
