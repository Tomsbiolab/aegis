import typer
import os
from typing_extensions import Annotated
from aegis.annotation import Annotation, read_file_with_fallback

features = ["gene", "transcript"]

app = typer.Typer(add_completion=False)
@app.command()
def main(
    annotation_file: Annotated[str, typer.Argument(
        help="Path to the input annotation GFF/GTF file."
    )],
    input_id_file: Annotated[str, typer.Argument(
        help="Input file with list of ids, one per line."
    )],
    annotation_name: Annotated[str, typer.Option(
        "-a", "--annotation-name", help="Annotation version, name or tag."
    )] = "{annotation-file}",
    feature_type: Annotated[str, typer.Option(
        "-f", "--feature-type", help=f"Identify feature level to be removed, based on input ids. Choose from {features}."
    )] = "gene",
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/",
    output_file: Annotated[str, typer.Option(
        "-o", "--output-file", help="Path to the output annotation filename, without extension."
    )] = "{annotation-name}_pruned"
):
    """
    Remove a list of ids, from a file, from the current annotation.
    """

    if annotation_name == "{annotation-file}":
        annotation_name = os.path.splitext(os.path.basename(annotation_file))[0]

    os.makedirs(output_folder, exist_ok=True)

    annotation = Annotation(name=annotation_name, annot_file_path=annotation_file)

    if output_file == "{annotation-name}.{ext}":
        output_file = f"{annotation_name}"

    if feature_type not in features:
        raise typer.BadParameter(f"Invalid feature level: {feature_type}. Choose from: {features}")

    if output_file == "{annotation-name}_pruned":
        output_file = f"{annotation_name}_pruned"

    input_ids = set()
   
    encoding = read_file_with_fallback(input_id_file)

    f_in = open(input_id_file, encoding=encoding)
    for line in f_in:
        if line.startswith("#"):
            continue
        id = line.strip()
        input_ids.add(id)
    f_in.close()

    if feature_type == "gene":
        annotation.remove_genes(input_ids)
    
    else:
        annotation.remove_transcripts(input_ids)

    output_file += ".gff3"
    annotation.export_gff(custom_path=output_folder, tag=output_file)


if __name__ == "__main__":
    app()
