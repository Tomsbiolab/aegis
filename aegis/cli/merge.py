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
    max_gene_overlap: Annotated[float, typer.Option(
        "-g", "--max-gene-overlap",
        help="Max gene overlap percentage. Maximum allowed gene overlap (0-100) with prioritised annotations. Genes exceeding this are excluded. Default: 100 (allows any degree of gene overlap)."
    )] = 100,
    max_exon_overlap: Annotated[float, typer.Option(
        "-e", "--max-exon-overlap",
        help="Max exon overlap percentage. Maximum allowed exon overlap (0-100) with prioritised annotations. Genes exceeding this are excluded. Default: 100 (allows any degree of exon overlap)."
    )] = 100,
    max_cds_overlap: Annotated[float, typer.Option(
        "-c", "--max-cds-overlap",
        help="Max CDS overlap percentage. Maximum allowed CDS overlap (0-100) with prioritised annotations. Genes exceeding this are excluded. Default: 100 (allows any degree of CDS overlap)."
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
    )] = False,
    no_collapse_CDSs: Annotated[bool, typer.Option(
        "--no-collapse-CDSs", help="Do not merge overlapping/adjacent CDS segments."
    )] = False
):
    """
    Merge two or more GFF3 annotation files.
    """

    collapse_exons = not no_collapse_exons
    collapse_CDSs = not no_collapse_CDSs

    if output_dir == "./aegis_output/":
        subfolder = True
    else:
        subfolder = False

    if skip_renaming:
        features = []
    else:
        features = ["transcript", "CDS", "exon", "UTR"]

    print(f"Loading base annotation from {annotation_files[0]}...")
    base_annotation = Annotation(annot_file_path=annotation_files[0], quiet=quiet, collapse_exons=collapse_exons, collapse_CDSs=collapse_CDSs)

    for annotation_file in annotation_files[1:]:

        print(f"Loading annotation to merge from {annotation_file}...")
        merge_annotation = Annotation(annot_file_path=annotation_file, quiet=quiet, collapse_exons=collapse_exons, collapse_CDSs=collapse_CDSs)

        print("Merging annotations...")
        base_annotation.merge(
            other=merge_annotation,
            max_gene_overlap=max_gene_overlap,
            max_exon_overlap=max_exon_overlap,
            max_cds_overlap=max_cds_overlap,
            features_to_rename=features
        )

    print(f"Writing merged annotation to {output_dir}...")
    base_annotation.export.gff(output_dir, quiet=quiet, subfolder=subfolder)

    print("Done.")

if __name__ == "__main__":
    app()
