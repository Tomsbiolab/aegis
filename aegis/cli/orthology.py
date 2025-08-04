import typer
import os
import shutil
from pathlib import Path
from typing import List
from typing_extensions import Annotated
from aegis.annotation import Annotation
from aegis.genome import Genome
from aegis.equivalence import Simple_annotation, perform_pairwise_comparison, run_command

app = typer.Typer(add_completion=False)

def split_callback(value:str):
    if value:
        return [item.strip() for item in value.split(",")]
    return []

@app.command()
def main(
    annotation_files: Annotated[List[str], typer.Argument(
        help="Path to the input annotation GFF/GTF file(s) associated to the same genome assembly. Input only one to measure gene overlaps within a single annotation, input several to compare between annotation files."
    )],
    genome_files: Annotated[str, typer.Option(
        "-g", "--genome-fastas", help="Genome assemblies corresponding to annotation files. Provide them in the same number and order, separated by commas. e.g. genomefile1,genomefile2,genomefile3,genomefile4",
        callback=split_callback
    )],
    annotation_names: Annotated[str, typer.Option(
        "-a", "--annotation-names", help="Annotation versions, names or tags. Provide them in the same number and order as the corresponding annotation files, separated by commas. e.g. name1,name2,name3,name4",
        callback=split_callback
    )] = "{annotation-filename(s)}",
    output_folder: Annotated[str, typer.Option(
        "-d", "--output-folder", help="Path to the output folder."
    )] = "./aegis_output/",
    group_names: Annotated[str, typer.Option(
        "-gn", "--group-names", help="Optional grouping of input annotations, into species for example. Use NA as a placemarker for annotation files without a group label. e.g. '-g group1,NA,group1,group2'",
        callback=split_callback
    )] = "",
    original_annotation_files: Annotated[str, typer.Option(
        "-ot", "--original-annotation-files", help="Should some of the annotations be a result of a liftover or coordinate transfer, you can optionally provide a list of the original files before the transfer, separated by commas. If at least 2 annotation files are being compared, conservation of synteny will be calculated wherever possible based on gene order before/after transfer. These original annotation files must be in the same number and order as the corresponding annotation files. Use NA as a placemarker for annotation files without an original annotation file. e.g. '-t original_file_1,NA,original_file_3'",
        callback=split_callback
    )] = "",
    reference_annotation: Annotated[str, typer.Option(
        "-r", "--reference-annotation", help="Select a single annotation, by providing its name/tag or filename, to use as a reference. Only matches to and from this annotation will be reported. Otherwise matches are reported between all annotations."
    )] = "None",
    include_single_blasts: Annotated[bool, typer.Option(
        "-b", "--include_single_blasts", help="Decide whether to report unidirectional (i.e. just fw or rv) blasts in the orthologue summary."
    )] = False,
    threads: Annotated[int, typer.Option(
        "-t", "--threads", help="Number of threads."
    )] = 5,
    skip_rbhs: Annotated[bool, typer.Option(
        "-rb", "--skip_RBHs", help="Decide whether to skip RBHs which are not RBBHs, these are reported by default in the orthologue summary."
    )] = False,
    verbose: Annotated[bool, typer.Option(
        "-v", "--verbose", help="Whether to include more details in the orthologue summary."
    )] = False
):
    """
    Provides a set of orthologues relationships leveraging external and internal tools such as Litfoff + AEGIS overlaps, LiftOn + AEGIS overlaps, MCScan, orthofinder and BLAST. Wherever relevant, tools are run reciprocally for an extra confidence mark in orthologous relationships.

    Script to perform gene correspondence analysis between multiple genomes.

    This script automates a bioinformatics pipeline that uses several tools
    (Liftoff, Lifton, DIAMOND, JCVI) to find orthologous genes between all possible
    pairs from a given list of genomes. For each pair, it performs:

    1.  Annotation Liftover: Transfers annotations between the two genomes based
        on sequence homology (using Liftoff and Lifton).
    2.  Reciprocal Best Hit (RBH): Uses DIAMOND (a fast alternative to BLAST) to
        find the best reciprocal homologs at the protein level.
    3.  Synteny and Collinearity: Uses the JCVI toolkit to identify conserved gene
        blocks (synteny) and find orthologs within that context.

    After all pairwise comparisons, it collects all unique proteomes and runs
    OrthoFinder to infer orthogroups across all species.

    The script is designed to be scalable, processing any number of genomes provided
    in the configuration section. All results for each pairwise comparison are stored
    in a separate, clearly named directory.

    All external tool commands are executed via Docker to ensure a reproducible
    environment.
    """

    if len(annotation_files) < 2:
        raise typer.BadParameter(f"At least 2 annotation-files must be provided.")
    
    if len(annotation_files) > 1 and annotation_files[-1].lower() in ("true", "false"):
        typer.echo(
            "⚠️  Detected extra value 'true' or 'false' at the end of positional arguments.\n"
            "👉 Did you mean to use the '--include_NAs' or '--simple' flags? Use them like this: '-n' or '-s' (no 'true' needed).",
            err=True,
        )
        raise typer.Exit(code=1)

    if annotation_names == "{annotation-filename(s)}":
        annotation_names = []
        for annotation_file in annotation_files:
            annotation_names.append(os.path.splitext(os.path.basename(annotation_file))[0])

    if len(annotation_files) != len(annotation_names):
        raise ValueError(f"The provided number of annotation name(s)/tag(s) do not match the number of annotation file(s).")

    if len(annotation_names) != len(set(annotation_names)):
        raise ValueError("Avoid repeated annotation tag(s)/name(s).")
    
    if len(annotation_files) != len(set(annotation_files)):
        raise ValueError("Avoid repeated annotation filename(s).")
    
    if len(genome_files) != len(set(genome_files)):
        raise ValueError("Avoid repeated genome assemblies. If looking to compare annotation versions associated to the same genome assembly, 'aegis-overlap' may be more appropriate.")

    if original_annotation_files != []:
        synteny = True
        if len(annotation_files) != len(original_annotation_files):
            raise ValueError(f"The provided number of original annotation files do not match the number of annotation file(s).")
        
    else:
        synteny = False
        original_annotation_files = ["NA"] * len(annotation_files)
    
    if group_names != "":
        if len(annotation_files) != len(group_names):
            raise ValueError(f"The provided number of groups do not match the number of annotation file(s).")
        
    else:
        group_names = ["NA"] * len(annotation_files)

    if reference_annotation != "None":
        if reference_annotation not in annotation_files and reference_annotation not in annotation_names:
            raise ValueError(f"The provided reference-annotation = {reference_annotation} is not present neither in annotation-files ({annotation_files}) nor annotation-names ({annotation_names}).")
        

    skip_unidirectional_blasts = not (include_single_blasts)

    if skip_rbhs and not skip_unidirectional_blasts:
        raise ValueError(f"Do not include single blasts if rbhs are to be skipped as these provide higher support for orthology.")

    genomes = [ Genome(name=f"{annotation_names[n]}_genome", genome_file_path=g) for n, g in enumerate(genome_files) ]

    annotations = []

    for n, annotation_file in enumerate(annotation_files):

        if original_annotation_files[n].lower() != "na":
            original_annotation = Annotation(name=f"{annotation_names[n]}_original", genome=genome_files[n], annot_file_path=original_annotation_files[n])
            annotations.append(Annotation(name=annotation_names[n], annot_file_path=annotation_file, original_annotation=original_annotation))
        else:
            annotations.append(Annotation(name=annotation_names[n], annot_file_path=annotation_file))

        if annotation_names[n] == reference_annotation or annotation_file == reference_annotation:
            annotations[n].target = True


    output_folder = Path(output_folder).resolve() / "orthologues"
    output_folder.mkdir(parents=True, exist_ok=True)
    output_folder = str(output_folder) + "/"

    RESULTS_DIRECTORY = Path(f"{output_folder}temp/")

    all_protein_files = set()

    for n1, a1 in enumerate(annotations):

        for n2, a2 in enumerate(annotations):

            if n1 == n2:
                continue
            # Create a dedicated directory for the results of this pair
            pair_directory = RESULTS_DIRECTORY / f"{a1.name}_vs_{a2.name}"
            
            ### CAMBIO ###
            # Capturar las rutas de los ficheros de proteínas devueltos por la función.
            protein_path_1, protein_path_2 = perform_pairwise_comparison(
                annot1=a1,
                annot2=a2,
                genome1=genomes[n1],
                genome2=genomes[n2],
                working_directory=pair_directory,
                num_threads=threads
            )
            # Añadir las rutas al conjunto
            all_protein_files.add(protein_path_1)
            all_protein_files.add(protein_path_2)

    print("\n--- ALL PAIRWISE ANALYSES COMPLETE ---")
    
    ### CAMBIO ###
    # --- 3. RUN ORTHOFINDER ON ALL PROTEOMES ---
    print("\n[STARTING ORTHOFINDER ANALYSIS]")
    
    # 3.1 Crear el directorio de entrada para OrthoFinder
    orthofinder_input_dir = RESULTS_DIRECTORY / "orthofinder_input"
    orthofinder_input_dir.mkdir(parents=True, exist_ok=True)
    print(f"  Created OrthoFinder input directory: {orthofinder_input_dir}")

    # 3.2 Copiar todos los ficheros de proteínas al directorio de entrada
    print("  Copying protein files to OrthoFinder input directory...")
    for protein_file in all_protein_files:
        if protein_file.exists():
            print(f"    - Copying {protein_file.name}")
            shutil.copy(protein_file, orthofinder_input_dir)
        else:
            print(f"    - WARNING: Protein file not found, skipping: {protein_file}")

    # 3.3 Ejecutar OrthoFinder
    print("\n  Running OrthoFinder... (this can take a very long time)")
    orthofinder_cmd = [
        "orthofinder",
        "-f", str(orthofinder_input_dir),
        "-t", str(threads) # Usar los mismos hilos que en los otros pasos
    ]
    
    # OrthoFinder se ejecuta en el directorio de resultados principal.
    # Sus resultados se crearán en una subcarpeta dentro de `orthofinder_input`.
    run_command(RESULTS_DIRECTORY, orthofinder_cmd)
    
    print("\n--- ORTHOFINDER ANALYSIS COMPLETE ---")


    
    print("\n--- ALL ANALYSES COMPLETE ---")
    print(f"All results can be found in subdirectories within: {RESULTS_DIRECTORY}")
    

    simple_annotations = []

    for n, a in enumerate(annotations):
        simple_annotations.append(Simple_annotation(a.name, a, group_names[n]))

    del annotations

    extra_tag = ""
    if verbose:
        extra_tag = "_verbose"

    protein_file_tag = "_proteins_g_id_main"

    for n1, a1 in enumerate(simple_annotations):

        for n2, a2 in enumerate(simple_annotations):

            if n1 == n2:
                continue
            
            if reference_annotation != "":
                if not a1.target and not a2.target:
                    continue

            mcscan_folder = RESULTS_DIRECTORY / f"{a1.name}_vs_{a2.name}"
            orthofinder_folder = RESULTS_DIRECTORY
            overlaps_path = f"MODIFY/"
            blasts_path = f"MODIFY/"

            a1.add_mcscan_equivalences(f"{mcscan_folder}{a1.name}.{a2.name}.anchors", "0", a2.name, group_names[n2])
            a1.add_mcscan_equivalences(f"{mcscan_folder}{a1.name}.{a2.name}.last.filtered", "0", a2.name, group_names[n2])
            a1.add_orthofinder_equivalences(f"{orthofinder_folder}Orthologues/{a1.name}{protein_file_tag}.tsv", a2.name, group_names[n2])
            
            a1.add_reciprocal_overlap_equivalences(overlaps_path, a1.name, a2.name, group_names[n2])

            a1.add_reciprocal_overlap_equivalences(overlaps_path, a1.name, a2.name, group_names[n2], liftoff=False)

            a1.add_blast_equivalences(f"{blasts_path}", a1.name, a2.name, group_names[n2], skip_rbhs=skip_rbhs, skip_unidirectional_blasts=skip_unidirectional_blasts)

        output_file = f"{equivalences_path}{a1.name}_equivalences{extra_tag}.tsv"

        output_file_filtered = f"{equivalences_path}{a1.name}_equivalences_filtered{extra_tag}.tsv"
        output_file_filtered_just_rbbhs_and_rbhs = f"{equivalences_path}{a1.name}_equivalences_filtered_just_rbbhs_and_rbhs{extra_tag}.tsv"
        output_file_filtered_just_rbbhs = f"{equivalences_path}{a1.name}_equivalences_filtered_just_rbbhs{extra_tag}.tsv"

        a1.export_summary_equivalences(output_file, verbose=verbose)

        if skip_rbhs and skip_unidirectional_blasts:
            a1.export_summary_equivalences(output_file_filtered_just_rbbhs, filtered=True, simple_rbh_blasts=False, unidirectional_blasts=False, verbose=verbose)

        elif skip_unidirectional_blasts:
            a1.export_summary_equivalences(output_file_filtered_just_rbbhs_and_rbhs, filtered=True, unidirectional_blasts=False, coverage_threshold=30, identity_threshold=30, verbose=verbose)

        else:
            a1.export_summary_equivalences(output_file_filtered, filtered=True, coverage_threshold=30, identity_threshold=30, verbose=verbose)

if __name__ == "__main__":
    app()
