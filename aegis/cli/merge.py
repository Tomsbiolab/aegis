import typer

from typing import List
from typing_extensions import Annotated

from ..annotation import Annotation

app = typer.Typer(add_completion=False)

def split_callback(value:str):
    if value:
        return [item.strip() for item in value.split(",")]
    return []

features: list = ["gene", "transcript", "CDS", "exon", "UTR"]

@app.command()
def main(
    annotation_files: Annotated[List[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file(s) associated to the same genome assembly. If gene and exon overlaps are chosen, priority will be given to the gene models of the first annotations and the subsequent ones."
    )],
    output_dir: Annotated[str, typer.Option(
        "-d", "--output-dir", help="Path to the output folder."
    )] = "./aegis_output/",
    gene_threshold: Annotated[float, typer.Option(
        "-g", "--gene-threshold",
        help="Gene overlap threshold in percentage (0-100). A gene will not be added if the maximum overlap of its exons with those of the prioritized annotations exceeds this value. The default 100 disables this check."
    )] = 100,
    exon_threshold: Annotated[float, typer.Option(
        "-e", "--exon-threshold",
        help="Exon overlap threshold in percentage (0-100). A gene will not be added if the maximum overlap of its exons with those of the prioritized annotations exceeds this value. The default 100 disables this check."
    )] = 100,
    skip_renaming: Annotated[bool, typer.Option(
        "--skip-renaming",
        help=f"Skip renaming of subgene features: transcript,CDS,exon,UTRs.",
    )] = False,
    quiet: Annotated[bool, typer.Option(
        "-q", "--quiet", help="Keeps terminal reporting to a minimum."
    )] = False,
    no_collapse_exons: Annotated[bool, typer.Option(
        "--no-collapse-exons", help="Do not merge overlapping/adjacent exons."
    )] = False
):
    """
    Merge two or more GFF3 annotation files.
    """

    collapse_exons = not no_collapse_exons

    if skip_renaming:
        features = []
    else:
        features = ["transcript", "CDS", "exon", "UTR"]

    print(f"Loading base annotation from {annotation_files[0]}...")
    base_annotation = Annotation(annot_file_path=annotation_files[0], quiet=quiet, collapse_exons=collapse_exons)

    for annotation_file in annotation_files[1:]:

        print(f"Loading annotation to merge from {annotation_file}...")
        merge_annotation = Annotation(annot_file_path=annotation_file, quiet=quiet, collapse_exons=collapse_exons)

        print("Merging annotations...")
        base_annotation.merge(
            other=merge_annotation,
            gene_overlap_threshold=gene_threshold,
            exon_overlap_threshold=exon_threshold,
            features_to_rename=features
        )

    print(f"Writing merged annotation to {output_dir}...")
    base_annotation.export.gff(output_dir, quiet=quiet)

    print("Done.")

if __name__ == "__main__":
    app()
