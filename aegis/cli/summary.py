import typer
import os

from typing_extensions import Annotated

from ..annotation import Annotation
from ..genome import Genome

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    genome_file: Annotated[str, typer.Argument(
        help="Path to the input genome FASTA file."
    )] = "",
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    genome_name: Annotated[str, typer.Option(
        "-g", "--genome-name", help="Genome assembly version, name or tag."
    )] = "{genome-file}",
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder."
    )] = "./aegis_output/stats/",
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
    header_id_tag: Annotated[str, typer.Option(
        "--header-id-tag", help="Extract chromosome/scaffold ID from FASTA header description by tag name (e.g., 'OriSeqID')."
    )] = "",
    header_id_regex: Annotated[str, typer.Option(
        "--header-id-regex", help="Extract chromosome/scaffold ID from FASTA header description using a regex capture group (e.g., 'OriSeqID=(\\S+)')."
    )] = "",
    gwh: Annotated[bool, typer.Option(
        "--gwh", help="Preset for Genome Warehouse (GWH) FASTA files. Automatically extracts original sequence IDs from 'OriSeqID=...' in headers."
    )] = False,
):
    """
    Outputs a summary with descriptive stats.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    if genome_name == "{genome-file}":
        genome_name = os.path.splitext(os.path.basename(genome_file))[0]

    os.makedirs(output_dir, exist_ok=True)

    if genome_file:
        genome = Genome(
            name=genome_name,
            genome_file_path=genome_file,
            quiet=quiet,
            header_id_tag=header_id_tag if header_id_tag != "" else None,
            header_id_regex=header_id_regex if header_id_regex != "" else None,
            gwh=gwh,
        )
        annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, genome=genome, quiet=quiet)
    else:
        annotation = Annotation(name=annotation_name, annot_file_path=annotation_file, genome=None, quiet=quiet)

    annotation.stats.update(output_dir=output_dir, export=True, quiet=quiet)


if __name__ == "__main__":
    app()
