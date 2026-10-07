import typer
import os
from typing import List
from typing_extensions import Annotated

from ..annotation import Annotation
from .utils import MERGE_PANEL, SUBFEATURE_PANEL, IO_PANEL, EXEC_PANEL, split_callback

app = typer.Typer(add_completion=False, no_args_is_help=True)

@app.command()
def main(
    annotation_files: Annotated[List[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file(s) associated to the same genome assembly. If overlaps are restricted, priority is given in the order files are passed."
    )] = [],

    # 1. Overlap & Merging Criteria
    max_gene_overlap: Annotated[float, typer.Option(
        "-go", "--max-gene-overlap", "--gene-overlap",
        help="Maximum allowed gene overlap percentage (0-100) with prioritized annotations. Genes exceeding this are excluded. Default: 100 (allows any degree of gene overlap).",
        rich_help_panel=MERGE_PANEL,
    )] = 100,
    max_exon_overlap: Annotated[float, typer.Option(
        "-eo", "--max-exon-overlap", "--exon-overlap",
        help="Maximum allowed exon overlap percentage (0-100) with prioritized annotations. Genes exceeding this are excluded. Default: 100 (allows any degree of exon overlap).",
        rich_help_panel=MERGE_PANEL,
    )] = 100,
    max_cds_overlap: Annotated[float, typer.Option(
        "-co", "--max-cds-overlap", "--cds-overlap",
        help="Maximum allowed CDS overlap percentage (0-100) with prioritized annotations. Genes exceeding this are excluded. Default: 100 (allows any degree of CDS overlap).",
        rich_help_panel=MERGE_PANEL,
    )] = 100,

    # 2. Subfeature & Model Options
    skip_renaming: Annotated[bool, typer.Option(
        "--skip-renaming",
        help="Skip renaming of subgene features (transcript, CDS, exon, UTRs).",
        rich_help_panel=SUBFEATURE_PANEL,
    )] = False,
    no_collapse_exons: Annotated[bool, typer.Option(
        "--no-collapse-exons", help="Do not merge overlapping/adjacent exons.",
        rich_help_panel=SUBFEATURE_PANEL,
    )] = False,
    no_collapse_CDSs: Annotated[bool, typer.Option(
        "--no-collapse-CDSs", help="Do not merge overlapping/adjacent CDS segments.",
        rich_help_panel=SUBFEATURE_PANEL,
    )] = False,

    # 3. Input / Output Options
    annotation_files_opt: Annotated[List[str], typer.Option(
        "-a", "--annotations", "--annotation-files",
        help="Path to input annotation GFF/GTF file(s). Overrides positional arguments if provided.",
        callback=split_callback,
        rich_help_panel=IO_PANEL,
    )] = [],
    annotation_names: Annotated[List[str], typer.Option(
        "-an", "--annotation-names",
        help="Optional names or tags for the annotations, separated by commas. [default: derived from filenames]",
        callback=split_callback,
        rich_help_panel=IO_PANEL,
    )] = [],
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output directory.",
        rich_help_panel=IO_PANEL,
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output filename (e.g., 'merged.gff3').",
        rich_help_panel=IO_PANEL,
    )] = "{annotation-names}.gff3",

    # 4. Execution & Debugging
    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help="Increase terminal reporting verbosity.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum.",
        rich_help_panel=EXEC_PANEL,
    )] = False,
):
    """
    Merge two or more GFF3 annotation files associated with the same genome assembly.
    """
    if verbose:
        quiet = False

    annot_files = list(annotation_files_opt) if annotation_files_opt else list(annotation_files)
    if not annot_files or len(annot_files) < 2:
        raise typer.BadParameter("At least two annotation files must be provided to merge. Provide as positional arguments or via -a/--annotations.")

    collapse_exons = not no_collapse_exons
    collapse_CDSs = not no_collapse_CDSs

    if output_dir == "./aegis_output/":
        subfolder = True
    else:
        subfolder = False

    if output_file == "{annotation-names}.gff3":
        output_file = None  # type: ignore

    if skip_renaming:
        features = []
    else:
        features = ["transcript", "CDS", "exon", "UTR"]

    name_1 = annotation_names[0] if len(annotation_names) > 0 else os.path.splitext(os.path.basename(annot_files[0]))[0]
    if not quiet:
        print(f"Loading base annotation from {annot_files[0]}...")
    base_annotation = Annotation(
        name=name_1,
        annot_file_path=annot_files[0],
        quiet=quiet,
        collapse_exons=collapse_exons,
        collapse_CDSs=collapse_CDSs,
        skip_coordinate_polishing=True,
    )

    for i, annotation_file in enumerate(annot_files[1:], start=1):
        name_i = annotation_names[i] if len(annotation_names) > i else os.path.splitext(os.path.basename(annotation_file))[0]
        if not quiet:
            print(f"Loading annotation to merge from {annotation_file}...")
        merge_annotation = Annotation(
            name=name_i,
            annot_file_path=annotation_file,
            quiet=quiet,
            collapse_exons=collapse_exons,
            collapse_CDSs=collapse_CDSs,
            skip_coordinate_polishing=True,
        )

        if not quiet:
            print("Merging annotations...")
        base_annotation.merge(
            other=merge_annotation,
            max_gene_overlap=max_gene_overlap,
            max_exon_overlap=max_exon_overlap,
            max_cds_overlap=max_cds_overlap,
            features_to_rename=tuple(features),
        )

    if not quiet:
        print(f"Writing merged annotation to {output_dir}...")
    base_annotation.export.gff(output_dir=output_dir, quiet=quiet, subfolder=subfolder, filename=output_file)

    if not quiet:
        print("Done.")

if __name__ == "__main__":
    app()
