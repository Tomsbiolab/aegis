# -*- coding: utf-8 -*-
"""
Script to perform gene correspondence analysis between two genomes.

This script automates a bioinformatics pipeline that uses several tools
(Liftoff, Lifton, DIAMOND, JCVI) to find orthologous genes between a "source"
and a "target" genome using three main methods:
1.  Annotation Liftover: Transfers annotations from one genome to another based
    on sequence homology (using Liftoff and Lifton).
2.  Reciprocal Best Hit (RBH): Uses DIAMOND (a fast alternative to BLAST) to
    find the best reciprocal homologs at the protein level.
3.  Synteny and Collinearity: Uses the JCVI toolkit to identify conserved gene
    blocks (synteny) and find orthologs within that context.

All external tool commands are executed via Docker to ensure a reproducible
environment.
"""
import os
import re
import shutil
import subprocess
from pathlib import Path

from aegis.genome import Genome
from aegis.annotation import Annotation

# --- HELPER FUNCTIONS ---

def run_docker_command(working_directory: Path, command: list):
    """
    Executes a generic command inside a Docker container.

    Args:
        working_directory (Path): The working directory to be mounted in Docker.
        command (list): The command and its arguments as a list of strings.

    Raises:
        subprocess.CalledProcessError: If the Docker command fails.
    """
    try:
        # Use capture_output=True to get stdout/stderr for better error reporting
        result = subprocess.run(command, check=True, capture_output=True, text=True)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error executing Docker command: {' '.join(command)}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        raise

def prepare_gff_for_jcvi(working_directory: Path, original_gff_path: Path, name: str) -> Path:
    """
    Copies and renames a GFF3 file to match the file name, as required by JCVI.

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

# --- MAIN FUNCTION ---

def main():
    """
    Main function that orchestrates the gene correspondence analysis.
    """
    # --- CONFIGURATION ---
    # Define all input paths and parameters here.
    
    # Paths to input genomes and annotations
    SOURCE_GENOME_PATH = Path('/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/test_genomes/T2T_ref_dapfit_organellefree_chr00.fasta')
    TARGET_GENOME_PATH = Path('/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/test_genomes/CBDRx_organellefree.fasta')
    SOURCE_GFF_PATH = Path('/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/test_genomes/5.1_on_T2T_ref_aegis_fcounts_dapfit_for_lifton.gff3')
    TARGET_GFF_PATH = Path('/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/test_genomes/CBDRx_on_CBDRx_aegis_fcounts_for_lifton.gff3')

    # Directory to store all intermediate and final results
    WORKING_DIRECTORY = Path('/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/tmp')

    # Tool parameters
    NUM_THREADS = 5

    # Docker image names
    LIFTOFF_IMG = "quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0"
    LIFTON_IMG = "antsanpaj/lifton"
    JCVI_IMG = "jcvi-image2"  # Assumes this image exists locally

    # --- 1. ENVIRONMENT SETUP ---
    print(f"Creating the working directory: {WORKING_DIRECTORY}")
    WORKING_DIRECTORY.mkdir(exist_ok=True)
    
    # --- 2. ANNOTATION LIFTOVER (Liftoff and Lifton) ---
    # Transfer annotations from one genome to another in both directions.
    print("\n[STEP 2.1] Running Liftoff to map annotations...")

    # Direction: Source -> Target
    liftoff_forward_gff = WORKING_DIRECTORY / "liftoff_fwd.gff"
    liftoff_fwd_cmd = [
        "docker", "run", "--rm", "-v", "/media:/media", "-w", str(WORKING_DIRECTORY),
        LIFTOFF_IMG, "liftoff", str(TARGET_GENOME_PATH), str(SOURCE_GENOME_PATH),
        "-g", str(SOURCE_GFF_PATH), "-o", str(liftoff_forward_gff)
    ]
    run_docker_command(WORKING_DIRECTORY, liftoff_fwd_cmd)

    # Direction: Target -> Source
    liftoff_reverse_gff = WORKING_DIRECTORY / "liftoff_rev.gff"
    liftoff_rev_cmd = [
        "docker", "run", "--rm", "-v", "/media:/media", "-w", str(WORKING_DIRECTORY),
        LIFTOFF_IMG, "liftoff", str(SOURCE_GENOME_PATH), str(TARGET_GENOME_PATH),
        "-g", str(TARGET_GFF_PATH), "-o", str(liftoff_reverse_gff)
    ]
    run_docker_command(WORKING_DIRECTORY, liftoff_rev_cmd)

    print("\n[STEP 2.2] Running Lifton to map annotations...")
    
    # Direction: Source -> Target
    lifton_forward_gff = WORKING_DIRECTORY / "lifton_fwd.gff3"
    lifton_fwd_cmd = [
        "docker", "run", "--rm", "-v", "/media:/media", "-w", str(WORKING_DIRECTORY),
        LIFTON_IMG, "lifton", "-g", str(SOURCE_GFF_PATH), "-o", str(lifton_forward_gff),
        "-copies", str(TARGET_GENOME_PATH), str(SOURCE_GENOME_PATH)
    ]
    run_docker_command(WORKING_DIRECTORY, lifton_fwd_cmd)

    # Direction: Target -> Source
    lifton_reverse_gff = WORKING_DIRECTORY / "lifton_rev.gff3"
    lifton_rev_cmd = [
        "docker", "run", "--rm", "-v", "/media:/media", "-w", str(WORKING_DIRECTORY),
        LIFTON_IMG, "lifton", "-g", str(TARGET_GFF_PATH), "-o", str(lifton_reverse_gff),
        "-copies", str(SOURCE_GENOME_PATH), str(TARGET_GENOME_PATH)
    ]
    run_docker_command(WORKING_DIRECTORY, lifton_rev_cmd)

    # --- 3. SEQUENCE EXTRACTION (Proteins and CDS) ---
    # Use AEGIS library to extract sequences needed for homology analyses
    # (proteins for DIAMOND, CDS for JCVI).
    print("\n[STEP 3] Extracting protein and CDS sequences...")
    
    # Process source genome and annotation
    print("  Processing SOURCE species...")
    source_genome = Genome('source_genome', str(SOURCE_GENOME_PATH))
    source_annot = Annotation('source_annot', str(SOURCE_GFF_PATH), source_genome)
    source_annot.generate_sequences(source_genome)
    source_annot.export_proteins(only_main=True, custom_path=str(WORKING_DIRECTORY), used_id="gene", verbose=False)
    source_annot.export_CDSs(only_main=True, custom_path=str(WORKING_DIRECTORY), used_id="gene", verbose=False)

    # Process target genome and annotation
    print("  Processing TARGET species...")
    target_genome = Genome('target_genome', str(TARGET_GENOME_PATH))
    target_annot = Annotation('target_annot', str(TARGET_GFF_PATH), target_genome)
    target_annot.generate_sequences(target_genome)
    target_annot.export_proteins(only_main=True, custom_path=str(WORKING_DIRECTORY), used_id="gene", verbose=False)
    target_annot.export_CDSs(only_main=True, custom_path=str(WORKING_DIRECTORY), used_id="gene", verbose=False)

    # Paths to the sequence files generated by 'aegis'
    features_dir = WORKING_DIRECTORY / "features"
    source_protein_fasta = features_dir / "source_annot_proteins_g_id_main.fasta"
    target_protein_fasta = features_dir / "target_annot_proteins_g_id_main.fasta"
    source_cds_fasta = features_dir / "source_annot_CDSs_g_id_main.fasta"
    target_cds_fasta = features_dir / "target_annot_CDSs_g_id_main.fasta"
    
    # --- 4. RECIPROCAL BEST HIT (RBH) HOMOLOGY SEARCH (DIAMOND) ---
    print("\n[STEP 4] Running reciprocal homology analysis with DIAMOND...")
    
    source_diamond_db = WORKING_DIRECTORY / "source_diamond_db"
    diamond_forward_result = WORKING_DIRECTORY / "diamond_fwd.txt"

    target_diamond_db = WORKING_DIRECTORY / "target_diamond_db"
    diamond_reverse_result = WORKING_DIRECTORY / "diamond_rev.txt"

    print("  Creating DIAMOND database for the source species...")
    makedb_cmd_1 = [
        "diamond", "makedb", "-p", str(NUM_THREADS),
        "--in", str(source_protein_fasta), "-d", str(source_diamond_db)
    ]
    subprocess.run(makedb_cmd_1, check=True)

    print("  Creating DIAMOND database for the target species...")
    makedb_cmd_2 = [
        "diamond", "makedb", "-p", str(NUM_THREADS),
        "--in", str(target_protein_fasta), "-d", str(target_diamond_db)
    ]
    subprocess.run(makedb_cmd_2, check=True)

    print("  Running DIAMOND search (Source -> Target)...")
    blastp_cmd_1 = [
        "diamond", "blastp", "--threads", str(NUM_THREADS),
        "--db", str(target_diamond_db), "--ultra-sensitive", "--out", str(diamond_forward_result),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "qlen", "slen", "length", "bitscore", "evalue",
        "--query", str(source_protein_fasta), "--max-target-seqs", "1",
        "--evalue", "0.00001", "--max-hsps", "1"
    ]
    subprocess.run(blastp_cmd_1, check=True)

    print("  Running DIAMOND search (Target -> Source)...")
    blastp_cmd_2 = [
        "diamond", "blastp", "--threads", str(NUM_THREADS),
        "--db", str(source_diamond_db), "--ultra-sensitive", "--out", str(diamond_reverse_result),
        "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", "qlen", "slen", "length", "bitscore", "evalue",
        "--query", str(target_protein_fasta), "--max-target-seqs", "1",
        "--evalue", "0.00001", "--max-hsps", "1"
    ]
    subprocess.run(blastp_cmd_2, check=True)
    
    # --- 5. SYNTENY-BASED ORTHOLOG SEARCH (JCVI) ---
    print("\n[STEP 5] Preparing files for JCVI synteny analysis...")

    # JCVI requires filenames to match the file names
    source_name = source_cds_fasta.stem
    target_name = target_cds_fasta.stem
    
    # Paths for files required by JCVI
    cleaned_source_cds = WORKING_DIRECTORY / f"{source_name}.cds"
    cleaned_target_cds = WORKING_DIRECTORY / f"{target_name}.cds"
    source_bed = WORKING_DIRECTORY / f"{source_name}.bed"
    target_bed = WORKING_DIRECTORY / f"{target_name}.bed"
    
    temp_source_cds = WORKING_DIRECTORY / "temp_source.cds"
    temp_target_cds = WORKING_DIRECTORY / "temp_target.cds"
    shutil.copyfile(source_cds_fasta, temp_source_cds)
    shutil.copyfile(target_cds_fasta, temp_target_cds)
    
    jcvi_format_cmd_1 = [
        "docker", "run", "-i", "--rm", "-v", "/media:/media", "-w", str(WORKING_DIRECTORY),
        JCVI_IMG, "python", "-m", "jcvi.formats.fasta", "format",
        str(temp_source_cds.name), str(cleaned_source_cds)
    ]
    run_docker_command(WORKING_DIRECTORY, jcvi_format_cmd_1)
    
    jcvi_format_cmd_2 = [
        "docker", "run", "-i", "--rm", "-v", "/media:/media", "-w", str(WORKING_DIRECTORY),
        JCVI_IMG, "python", "-m", "jcvi.formats.fasta", "format",
        str(temp_target_cds.name), str(cleaned_target_cds)
    ]
    run_docker_command(WORKING_DIRECTORY, jcvi_format_cmd_2)

    # Clean up temporary files
    temp_source_cds.unlink()
    temp_target_cds.unlink()

    # 5.2. Convert GFF3 to BED format
    source_gff_for_jcvi = prepare_gff_for_jcvi(WORKING_DIRECTORY, SOURCE_GFF_PATH, source_name)
    target_gff_for_jcvi = prepare_gff_for_jcvi(WORKING_DIRECTORY, TARGET_GFF_PATH, target_name)

    gff_to_bed_cmd_1 = [
        "docker", "run", "-i", "--rm", "-v", "/media:/media", "-w", str(WORKING_DIRECTORY),
        JCVI_IMG, "python", "-m", "jcvi.formats.gff", "bed", "--type=mRNA",
        "--key=Parent", "--primary_only", str(source_gff_for_jcvi.name), "-o", str(source_bed)
    ]
    run_docker_command(WORKING_DIRECTORY, gff_to_bed_cmd_1)

    gff_to_bed_cmd_2 = [
        "docker", "run", "-i", "--rm", "-v", "/media:/media", "-w", str(WORKING_DIRECTORY),
        JCVI_IMG, "python", "-m", "jcvi.formats.gff", "bed", "--type=mRNA",
        "--key=Parent", "--primary_only", str(target_gff_for_jcvi.name), "-o", str(target_bed)
    ]
    run_docker_command(WORKING_DIRECTORY, gff_to_bed_cmd_2)
    
    # 5.3. Run JCVI ortholog analysis
    print("\n[STEP 5.3] Running JCVI ortholog analysis (this may take a while)...")
    
    # The `jcvi.compara.catalog ortholog` tool performs a reciprocal search and
    # generates an anchors file (`.anchors`) listing the ortholog pairs.
    jcvi_ortho_cmd_1 = [
        "docker", "run", "-i", "--rm", "-v", "/media:/media", "-w", str(WORKING_DIRECTORY),
        JCVI_IMG, "python", "-m", "jcvi.compara.catalog", "ortholog",
        source_name, target_name, "--no_strip_names"
    ]
    run_docker_command(WORKING_DIRECTORY, jcvi_ortho_cmd_1)
    
    jcvi_ortho_cmd_2 = [
        "docker", "run", "-i", "--rm", "-v", "/media:/media", "-w", str(WORKING_DIRECTORY),
        JCVI_IMG, "python", "-m", "jcvi.compara.catalog", "ortholog",
        target_name, source_name, "--no_strip_names"
    ]
    run_docker_command(WORKING_DIRECTORY, jcvi_ortho_cmd_2)

    print("\n--- ANALYSIS COMPLETE ---")
    print(f"All results can be found in: {WORKING_DIRECTORY}")


if __name__ == "__main__":
    main()