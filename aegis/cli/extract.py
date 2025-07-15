import typer
from typing_extensions import Annotated, Literal
from typing import List
from aegis.genome import Genome
from aegis.annotation import Annotation
import sys
import re

Feature = Literal["gene", "transcript", "CDS", "protein", "promoter"]

IDs = Literal["gene", "transcript", "CDS", "feature"]

Mode = Literal["all", "main", "unique", "unique_per_gene"]

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
    feature_type: Annotated[List[str], typer.Option(
        "-f", "--feature-type", help=f"One or more feature types to extract. Choose from: {Feature}."
    )],
    mode: Annotated[List[str], typer.Option(
        "-m", "--mode", help=f"Extract main or all features, or both. Choose from: {Mode}."
    )],
    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help=f"Verbose."
    )] = False,
    feature_id: Annotated[str, typer.Option(
        "-i", "--feature-id", help=f"Most specific feature ID used in fasta header outputs. E.g. you may want to export transcripts but associated to gene ids directly instead of using the transcript feature IDs. Choose from: {IDs}."
    )] = "feature"
):
    genome = Genome(name="genome", genome_file_path=genome_fasta)
    annotation = Annotation(name="annotation", annot_file_path=annotation_gff, genome=genome)

    annotation.generate_sequences(genome)

    valid_features = [f.value for f in Feature.__args__]
    for f_type in feature_type:
        if f_type not in valid_features:
            raise typer.BadParameter(f"Invalid feature type: {f_type}. Choose from: {valid_features}")

    valid_modes = [m.value for m in Mode.__args__]
    for m_type in mode:
        if m_type not in valid_modes:
            raise typer.BadParameter(f"Invalid mode: {m_type}. Choose from: {valid_modes}")

    valid_ids = [i.value for i in IDs.__args__]
    if feature_id not in valid_ids:
        raise typer.BadParameter(f"Invalid feature ID: {feature_id}. Choose from: {valid_ids}")

    if "transcript" in feature_type:

        if "gene" in feature_id:
            used_id = "gene"
        else:
            used_id = "transcript"

        if "all" in mode:
            annotation.export_transcripts(only_main=False, verbose=verbose, custom_path=output_dir, used_id=used_id)
        if "main" in mode or "unique" in mode or "unique_per_gene" in mode:
            annotation.export_transcripts(custom_path=output_dir, verbose=verbose, used_id=used_id)

    if "protein" in feature_type:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        elif "CDS" in feature_id:
            used_id = "CDS"
        else:
            used_id = "protein"

        if "unique_per_gene" in mode:
            annotation.export_proteins(only_main=False, custom_path=output_dir, verbose=verbose, unique_proteins_per_gene=True, used_id=used_id)
        elif "unique" in mode:
            annotation.export_unique_proteins(custom_path=output_dir, verbose=verbose)
        else:
            if "all" in mode:
                annotation.export_proteins(only_main=False, custom_path=output_dir, verbose=verbose, used_id=used_id)
            if "main" in mode or "unique" in mode or "unique_per_gene" in mode:
                annotation.export_proteins(custom_path=output_dir, verbose=verbose, used_id=used_id)

    if "CDS" in feature_type:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        else:
            used_id = "CDS"

        if "all" in mode:
            annotation.export_CDSs(only_main=False, custom_path=output_dir, verbose=verbose, used_id=used_id)
        if "main" in mode or "unique" in mode or "unique_per_gene" in mode:
            annotation.export_CDSs(custom_path=output_dir, verbose=verbose, used_id=used_id)

    if "promoter" in feature_type:

        if "gene" in feature_id:
            used_id = "gene"
        elif "transcript" in feature_id:
            used_id = "transcript"
        else:
            used_id = "promoter"

        if "all" in mode:
            annotation.export_promoters(only_main=False, custom_path=output_dir, verbose=verbose, used_id=used_id)
        if "main" in mode or "unique" in mode or "unique_per_gene" in mode:
            annotation.export_promoters(custom_path=output_dir, verbose=verbose, used_id=used_id)

if __name__ == "__main__":
    app()