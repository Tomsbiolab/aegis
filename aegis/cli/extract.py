import typer
from typing_extensions import Annotated, Literal
from typing import List
from aegis.genome import Genome
from aegis.annotation import Annotation

Feature = Literal["transcript", "protein", "CDS", "gene", "promoter"]


app = typer.Typer(help="Extract sequences from a genome based on an annotation.")

@app.command()
def main(
    genome_fasta: Annotated[str, typer.Option(
        "-g", "--genome-fasta", help="Path to the input genome FASTA file."
    )],
    annotation_gff: Annotated[str, typer.Option(
        "-a", "--annotation-gff", help="Path to the input annotation GFF3 file."
    )],
    output_dir: Annotated[str, typer.Option(
        "-o", "--output-dir", help="Path to the directory where output FASTA files will be saved."
    )],
    feature_type: Annotated[List[Feature], typer.Option(
        "-t", "--feature-type", help="One or more feature types to extract. Choose from: transcript, protein.",
        case_sensitive=False
    )]
):
    genome = Genome(name="genome", genome_file_path=genome_fasta)
    annotation = Annotation(name="annotation", annot_file_path=annotation_gff, genome=genome)

    annotation.generate_sequences(genome)

    if 'transcript' in feature_type:
        annotation.export_transcripts(custom_path=output_dir)

    if 'protein' in feature_type:
        annotation.export_proteins(custom_path=output_dir)

if __name__ == "__main__":
    app()