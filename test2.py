# -*- coding: utf-8 -*-
"""
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
import os
import re
import shutil
import subprocess
from pathlib import Path
from itertools import combinations
from dataclasses import dataclass

# You might need to install aegis: pip install aegis-scaffold
from aegis.genome import Genome
from aegis.annotation import Annotation

# --- DATA STRUCTURE FOR GENOME INPUT ---
@dataclass
class GenomeData:
    """A simple class to hold paths and a name for a single species."""
    name: str
    genome_path: Path
    gff_path: Path

# --- HELPER FUNCTIONS ---

def run_command(working_directory: Path, command: list):
    """
    Executes a generic command inside a Docker container.

    Args:
        working_directory (Path): The working directory for the command.
        command (list): The command and its arguments as a list of strings.

    Raises:
        subprocess.CalledProcessError: If the command fails.
    """
    try:
        print(f"    Executing command: {' '.join(command)}")
        result = subprocess.run(command, check=True, capture_output=True, text=True, cwd=working_directory)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {' '.join(command)}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        raise

def prepare_gff_for_jcvi(working_directory: Path, original_gff_path: Path, name: str) -> Path:
    """
    Copies and renames a GFF3 file to match the base name, as required by JCVI.

    Args:
        working_directory (Path): The directory where the copy will be created.
        original_gff_path (Path): Path to the original GFF file.
        name (str): The base name to be used for the new file.

    Returns:
        Path: The path to the new, renamed GFF file.
    """
    renamed_gff_path = working_directory / f"{name}.gff3"
    print(f"    Copying and renaming GFF for JCVI: '{renamed_gff_path.name}'")
    shutil.copyfile(original_gff_path, renamed_gff_path)
    return renamed_gff_path

# --- CORE ANALYSIS FUNCTION FOR A PAIR OF GENOMES ---

def perform_pairwise_comparison(species1: GenomeData, species2: GenomeData, working_directory: Path, num_threads: int) -> tuple[Path, Path]:
    """
    Runs the full analysis pipeline for a single pair of genomes.

    Args:
        species1 (GenomeData): Data for the first species.
        species2 (GenomeData): Data for the second species.
        working_directory (Path): The directory to store results for this pair.
        num_threads (int): Number of threads for multicore tools.
    
    ### CAMBIO ###
    Returns:
        tuple[Path, Path]: A tuple containing the paths to the generated protein
                           FASTA files for species1 and species2.
    """
    print(f"\n[STARTING ANALYSIS] Comparing {species1.name} with {species2.name}")
    print(f"Results will be stored in: {working_directory}")
    working_directory.mkdir(exist_ok=True)

    # --- 1. ANNOTATION LIFTOVER (Liftoff and Lifton) ---
    print("\n[STEP 1.1] Running Liftoff to map annotations...")

    # Direction: Species 1 -> Species 2
    liftoff_1_to_2_gff = working_directory / f"liftoff_{species1.name}_to_{species2.name}.gff"
    liftoff_cmd_1 = [
        "liftoff", str(species2.genome_path), str(species1.genome_path),
        "-g", str(species1.gff_path), "-o", str(liftoff_1_to_2_gff)
    ]
    # run_command(working_directory, liftoff_cmd_1)

    # Direction: Species 2 -> Species 1
    liftoff_2_to_1_gff = working_directory / f"liftoff_{species2.name}_to_{species1.name}.gff"
    liftoff_cmd_2 = [
        "liftoff", str(species1.genome_path), str(species2.genome_path),
        "-g", str(species2.gff_path), "-o", str(liftoff_2_to_1_gff)
    ]
    # run_command(working_directory, liftoff_cmd_2)

    print("\n[STEP 1.2] Running Lifton to map annotations...")

    # Direction: Species 1 -> Species 2
    lifton_1_to_2_gff = working_directory / f"lifton_{species1.name}_to_{species2.name}.gff3"
    lifton_cmd_1 = [
        "lifton", "-g", str(species1.gff_path), "-o", str(lifton_1_to_2_gff),
        "-copies", str(species2.genome_path), str(species1.genome_path)
    ]
    # run_command(working_directory, lifton_cmd_1)

    # Direction: Species 2 -> Species 1
    lifton_2_to_1_gff = working_directory / f"lifton_{species2.name}_to_{species1.name}.gff3"
    lifton_cmd_2 = [
        "lifton", "-g", str(species2.gff_path), "-o", str(lifton_2_to_1_gff),
        "-copies", str(species1.genome_path), str(species2.genome_path)
    ]
    # run_command(working_directory, lifton_cmd_2)

    # --- 2. SEQUENCE EXTRACTION (Proteins and CDS) ---
    print("\n[STEP 2] Extracting protein and CDS sequences...")
    
    # Process species 1
    print(f"  Processing {species1.name}...")
    genome1 = Genome(species1.name, str(species1.genome_path))
    annot1 = Annotation(f"{species1.name}_annot", str(species1.gff_path), genome1)
    annot1.generate_sequences(genome1)
    annot1.export_proteins(only_main=True, custom_path=str(working_directory), used_id="gene", verbose=False)
    annot1.export_CDSs(only_main=True, custom_path=str(working_directory), used_id="gene", verbose=False)

    # Process species 2
    print(f"  Processing {species2.name}...")
    genome2 = Genome(species2.name, str(species2.genome_path))
    annot2 = Annotation(f"{species2.name}_annot", str(species2.gff_path), genome2)
    annot2.generate_sequences(genome2)
    annot2.export_proteins(only_main=True, custom_path=str(working_directory), used_id="gene", verbose=False)
    annot2.export_CDSs(only_main=True, custom_path=str(working_directory), used_id="gene", verbose=False)

    features_dir = working_directory / "features"
    protein_fasta_1 = features_dir / f"{species1.name}_annot_proteins_g_id_main.fasta"
    protein_fasta_2 = features_dir / f"{species2.name}_annot_proteins_g_id_main.fasta"
    cds_fasta_1 = features_dir / f"{species1.name}_annot_CDSs_g_id_main.fasta"
    cds_fasta_2 = features_dir / f"{species2.name}_annot_CDSs_g_id_main.fasta"

    # --- 3. RECIPROCAL BEST HIT (RBH) HOMOLOGY SEARCH (DIAMOND) ---
    print("\n[STEP 3] Running reciprocal homology analysis with DIAMOND...")
    
    diamond_db_1 = working_directory / f"{species1.name}_diamond_db"
    diamond_db_2 = working_directory / f"{species2.name}_diamond_db"
    diamond_result_1_to_2 = working_directory / f"diamond_{species1.name}_to_{species2.name}.txt"
    diamond_result_2_to_1 = working_directory / f"diamond_{species2.name}_to_{species1.name}.txt"

    print(f"  Creating DIAMOND database for {species1.name}...")
    makedb_cmd_1 = [
        "diamond", "makedb", "-p", str(num_threads), "--in", str(protein_fasta_1), "-d", str(diamond_db_1)
    ]
    run_command(working_directory, makedb_cmd_1)

    print(f"  Creating DIAMOND database for {species2.name}...")
    makedb_cmd_2 = [
        "diamond", "makedb", "-p", str(num_threads), "--in", str(protein_fasta_2), "-d", str(diamond_db_2)
    ]
    run_command(working_directory, makedb_cmd_2)

    print(f"  Running DIAMOND search ({species1.name} -> {species2.name})...")
    blastp_cmd_1 = [
        "diamond", "blastp", "--threads", str(num_threads), "--db", str(diamond_db_2), "--ultra-sensitive", 
        "--out", str(diamond_result_1_to_2), "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", 
        "qlen", "slen", "length", "bitscore", "evalue", "--query", str(protein_fasta_1), 
        "--max-target-seqs", "1", "--evalue", "0.00001", "--max-hsps", "1"
    ]
    run_command(working_directory, blastp_cmd_1)

    print(f"  Running DIAMOND search ({species2.name} -> {species1.name})...")
    blastp_cmd_2 = [
        "diamond", "blastp", "--threads", str(num_threads), "--db", str(diamond_db_1), "--ultra-sensitive",
        "--out", str(diamond_result_2_to_1), "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", 
        "qlen", "slen", "length", "bitscore", "evalue", "--query", str(protein_fasta_2),
        "--max-target-seqs", "1", "--evalue", "0.00001", "--max-hsps", "1"
    ]
    run_command(working_directory, blastp_cmd_2)

    # --- 4. SYNTENY-BASED ORTHOLOG SEARCH (JCVI) ---
    print("\n[STEP 4] Preparing files for JCVI synteny analysis...")

    # JCVI requires specific file naming conventions
    jcvi_name_1 = f"{species1.name}_jcvi"
    jcvi_name_2 = f"{species2.name}_jcvi"
    
    cleaned_cds_1 = working_directory / f"{jcvi_name_1}.cds"
    cleaned_cds_2 = working_directory / f"{jcvi_name_2}.cds"
    bed_file_1 = working_directory / f"{jcvi_name_1}.bed"
    bed_file_2 = working_directory / f"{jcvi_name_2}.bed"
    
    # 4.1 Format CDS files for JCVI
    jcvi_format_cmd_1 = ["python", "-m", "jcvi.formats.fasta", "format", str(cds_fasta_1), str(cleaned_cds_1)]
    run_command(working_directory, jcvi_format_cmd_1)
    
    jcvi_format_cmd_2 = ["python", "-m", "jcvi.formats.fasta", "format", str(cds_fasta_2), str(cleaned_cds_2)]
    run_command(working_directory, jcvi_format_cmd_2)

    # 4.2 Convert GFF3 to BED format
    gff_for_jcvi_1 = prepare_gff_for_jcvi(working_directory, species1.gff_path, jcvi_name_1)
    gff_for_jcvi_2 = prepare_gff_for_jcvi(working_directory, species2.gff_path, jcvi_name_2)

    gff_to_bed_cmd_1 = [
        "python", "-m", "jcvi.formats.gff", "bed", "--type=mRNA",
        "--key=Parent", "--primary_only", str(gff_for_jcvi_1.name), "-o", str(bed_file_1)
    ]
    run_command(working_directory, gff_to_bed_cmd_1)

    gff_to_bed_cmd_2 = [
        "python", "-m", "jcvi.formats.gff", "bed", "--type=mRNA",
        "--key=Parent", "--primary_only", str(gff_for_jcvi_2.name), "-o", str(bed_file_2)
    ]
    run_command(working_directory, gff_to_bed_cmd_2)
    
    # 4.3. Run JCVI ortholog analysis
    print("\n[STEP 4.3] Running JCVI ortholog analysis (this may take a while)...")
    
    # This command finds reciprocal best hits and creates the .anchors file
    jcvi_ortho_cmd = [
        "python", "-m", "jcvi.compara.catalog", "ortholog",
        jcvi_name_1, jcvi_name_2, "--no_strip_names"
    ]
    run_command(working_directory, jcvi_ortho_cmd)
    
    print(f"\n--- PAIRWISE ANALYSIS COMPLETE: {species1.name} vs {species2.name} ---")

    ### CAMBIO ###
    # Devuelve las rutas a los ficheros de proteínas generados.
    return protein_fasta_1, protein_fasta_2

# --- MAIN FUNCTION ---

def main():
    """
    Main function to orchestrate the gene correspondence analysis for all
    pairs of genomes.
    """
    # --- CONFIGURATION ---
    
    # Directory to store all analysis results. A subdirectory will be created
    # for each pairwise comparison.
    RESULTS_DIRECTORY = Path('/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/pairwise_analysis_results')

    # Tool parameters
    NUM_THREADS = 5

    # List of all genomes and their annotations to be compared.
    # Add as many GenomeData objects as needed.
    GENOME_DATA_LIST = [
        GenomeData(
            name="T2T_Ref",
            genome_path=Path('/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/test_genomes/T2T_ref_dapfit_organellefree_chr00.fasta'),
            gff_path=Path('/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/test_genomes/5.1_on_T2T_ref_aegis_fcounts_dapfit_for_lifton.gff3')
        ),
        GenomeData(
            name="CBDRx",
            genome_path=Path('/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/test_genomes/CBDRx_organellefree.fasta'),
            gff_path=Path('/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/test_genomes/CBDRx_on_CBDRx_aegis_fcounts_for_lifton.gff3')
        ),
        # Example: Add a third genome for comparison
        # GenomeData(
        #     name="Another_Species",
        #     genome_path=Path('/path/to/another/genome.fasta'),
        #     gff_path=Path('/path/to/another/annotation.gff3')
        # ),
    ]

    # --- 1. ENVIRONMENT SETUP ---
    print(f"Creating the main results directory: {RESULTS_DIRECTORY}")
    RESULTS_DIRECTORY.mkdir(exist_ok=True)
    
    # --- 2. GENERATE PAIRS AND RUN ANALYSIS ---
    if len(GENOME_DATA_LIST) < 2:
        print("Error: At least two genomes are required for comparison.")
        return

    # Create all unique 2-genome combinations
    genome_pairs = combinations(GENOME_DATA_LIST, 2)
    
    ### CAMBIO ###
    # Conjunto para almacenar las rutas a todos los ficheros de proteomas únicos.
    # Se usa un conjunto para evitar duplicados, aunque los ficheros se generen
    # en directorios diferentes para cada par. Al copiarlos, se sobreescribirán
    # y quedará solo uno por especie.
    all_protein_files = set()

    for species1, species2 in genome_pairs:
        # Create a dedicated directory for the results of this pair
        pair_directory = RESULTS_DIRECTORY / f"{species1.name}_vs_{species2.name}"
        
        ### CAMBIO ###
        # Capturar las rutas de los ficheros de proteínas devueltos por la función.
        protein_path_1, protein_path_2 = perform_pairwise_comparison(
            species1=species1,
            species2=species2,
            working_directory=pair_directory,
            num_threads=NUM_THREADS
        )
        # Añadir las rutas al conjunto
        all_protein_files.add(protein_path_1)
        all_protein_files.add(protein_path_2)


    print("\n--- ALL PAIRWISE ANALYSES COMPLETE ---")
    
    ### CAMBIO ###
    # --- 3. RUN ORTHOFINDER ON ALL PROTEOMES ---
    print("\n[STARTING ORTHOFINDER ANALYSIS]")
    
    try:
        # 3.1 Crear el directorio de entrada para OrthoFinder
        orthofinder_input_dir = RESULTS_DIRECTORY / "orthofinder_input"
        orthofinder_input_dir.mkdir(exist_ok=True)
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
            "-t", str(NUM_THREADS) # Usar los mismos hilos que en los otros pasos
        ]
        
        # OrthoFinder se ejecuta en el directorio de resultados principal.
        # Sus resultados se crearán en una subcarpeta dentro de `orthofinder_input`.
        run_command(RESULTS_DIRECTORY, orthofinder_cmd)
        
        print("\n--- ORTHOFINDER ANALYSIS COMPLETE ---")

    except Exception as e:
        print("\n!!!!!! FAILED ORTHOFINDER ANALYSIS !!!!!!")
        print(f"Encountered an error: {e}")
    
    print("\n--- ALL ANALYSES COMPLETE ---")
    print(f"All results can be found in subdirectories within: {RESULTS_DIRECTORY}")


if __name__ == "__main__":
    main()