import typer
from typing_extensions import Annotated
from aegis.genome import Genome
from aegis.annotation import Annotation

app = typer.Typer(help="Find overlapping features in an annotation file and export them to a table.")

@app.command()
def main(
    annotation_gff: Annotated[str, typer.Option(
        "-a", "--annotation-gff", help="Path to the input annotation GFF3 file."
    )],
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output TSV file."
    )],
):
    """
    Finds overlapping features in an annotation file and exports them to a table.
    """
    genome = Genome(name="genome", genome_file_path=None)
    annotation = Annotation(name="annotation", annot_file_path=annotation_gff, genome=genome)

    annotation.self_overlap_genes()

    annotation.export_equivalences(output_file=output_file, export_self=True, export_csv=True, return_df=False)


if __name__ == "__main__":
    app()
