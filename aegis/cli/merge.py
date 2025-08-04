import typer
from typing import List
from typing_extensions import Annotated
from aegis.annotation import Annotation

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
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/",
    gene_threshold: Annotated[float, typer.Option(
        "-g", "--gene-threshold",
        help="Gene overlap threshold in percentage (0-100). A gene will not be added if the maximum overlap of its exons with those of the prioritized annotations exceeds this value. The default 100 disables this check."
    )] = 100,
    exon_threshold: Annotated[float, typer.Option(
        "-e", "--exon-threshold",
        help="Exon overlap threshold in percentage (0-100). A gene will not be added if the maximum overlap of its exons with those of the prioritized annotations exceeds this value. The default 100 disables this check."
    )] = 100,
    skip_feature_renaming: Annotated[bool, typer.Option(
        "-f", "--feature_type",
        help=f"Skip feature renaming of transcript,CDS,exon,UTR.",
    )] = False,
):
    """
    Merge two GFF3 annotation files.
    """

    if skip_feature_renaming:
        features = []
    else:
        features = ["transcript", "CDS", "exon", "UTR"]

    print(f"Loading base annotation from {annotation_files[0]}...")
    base_annotation = Annotation(annot_file_path=annotation_files[0])

    for annotation_file in annotation_files[1:]:

        print(f"Loading annotation to merge from {annotation_file}...")
        merge_annotation = Annotation(annot_file_path=annotation_file)

        print("Merging annotations...")
        base_annotation.merge(
            other=merge_annotation,
            gene_overlap_threshold=gene_threshold,
            exon_overlap_threshold=exon_threshold,
            features_to_rename=features
        )

    print(f"Writing merged annotation to {output_folder}...")
    base_annotation.export_gff(output_folder)

    print("Done.")

if __name__ == "__main__":
    app()
