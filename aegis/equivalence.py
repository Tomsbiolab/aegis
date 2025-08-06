import pandas as pd
import re
import time
import os
import subprocess

from pathlib import Path


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


def perform_pairwise_comparison(annot1: object, annot2: object, genome1: object, genome2: object, working_directory: Path, num_threads: int) -> tuple[Path, Path]:

    annot1.update_attributes(clean=True, symbols=False, symbols_as_descriptors=False)
    annot1.export_gff(custom_path=str(working_directory), tag=f"{annot1.name}_clean.gff3", subfolder=False)

    annot2.update_attributes(clean=True, symbols=False, symbols_as_descriptors=False)
    annot2.export_gff(custom_path=str(working_directory), tag=f"{annot2.name}_clean.gff3", subfolder=False)

    annot1_lifton = annot1.copy()
    annot1_lifton.CDS_to_CDS_segment_ids()
    annot1_lifton.update_attributes(clean=True, symbols=False, symbols_as_descriptors=False)
    annot1_lifton.export_gff(custom_path=str(working_directory), tag=f"{annot1.name}_lifton_clean.gff3", subfolder=False, no_1bp_features=True)
    del annot1_lifton

    annot2_lifton = annot2.copy()
    annot2_lifton.CDS_to_CDS_segment_ids()
    annot2_lifton.update_attributes(clean=True, symbols=False, symbols_as_descriptors=False)
    annot2_lifton.export_gff(custom_path=str(working_directory), tag=f"{annot2.name}_lifton_clean.gff3", subfolder=False, no_1bp_features=True)
    del annot2_lifton

    print(f"\n[STARTING ANALYSIS] Comparing {annot1.name} with {annot2.name}")
    print(f"Results will be stored in: {working_directory}")
    working_directory.mkdir(parents=True, exist_ok=True)

    # --- 1. ANNOTATION LIFTOVER (Liftoff and Lifton) ---
    print("\n[STEP 1.1] Running Liftoff to map annotations...")

    # Direction: Species 1 -> Species 2
    liftoff_1_to_2_gff = working_directory / f"liftoff_{annot1.name}_to_{annot2.name}.gff"
    liftoff_cmd_1 = [
        "liftoff", str(genome2.file), str(genome1.file),
        "-g", f"{annot1.name}_clean.gff3", "-o", str(liftoff_1_to_2_gff)
    ]
    run_command(working_directory, liftoff_cmd_1)

    # Direction: Species 2 -> Species 1
    liftoff_2_to_1_gff = working_directory / f"liftoff_{annot2.name}_to_{annot1.name}.gff"
    liftoff_cmd_2 = [
        "liftoff", str(genome1.file), str(genome2.file),
        "-g", f"{annot2.name}_clean.gff3", "-o", str(liftoff_2_to_1_gff)
    ]
    run_command(working_directory, liftoff_cmd_2)

    print("\n[STEP 1.2] Running Lifton to map annotations...")

    # Direction: Species 1 -> Species 2
    lifton_1_to_2_gff = working_directory / f"lifton_{annot1.name}_to_{annot2.name}.gff3"
    lifton_cmd_1 = [
        "lifton", "-g", f"{annot1.name}_lifton_clean.gff3", "-o", str(lifton_1_to_2_gff),
        "-copies", str(genome2.file), str(genome1.file)
    ]
    run_command(working_directory, lifton_cmd_1)

    # Direction: Species 2 -> Species 1
    lifton_2_to_1_gff = working_directory / f"lifton_{annot2.name}_to_{annot1.name}.gff3"
    lifton_cmd_2 = [
        "lifton", "-g", f"{annot2.name}_lifton_clean.gff3", "-o", str(lifton_2_to_1_gff),
        "-copies", str(genome1.file), str(genome2.file)
    ]
    run_command(working_directory, lifton_cmd_2)

    # --- 2. SEQUENCE EXTRACTION (Proteins and CDS) ---
    print("\n[STEP 2] Extracting protein and CDS sequences...")
    
    # Process species 1
    print(f"  Processing {annot1.name}...")
    annot1.generate_sequences(genome1)
    annot1.export_proteins(only_main=True, custom_path=str(working_directory), used_id="gene", verbose=False)
    annot1.export_CDSs(only_main=True, custom_path=str(working_directory), used_id="gene", verbose=False)

    # Process species 2
    print(f"  Processing {annot2.name}...")
    annot2.generate_sequences(genome2)
    annot2.export_proteins(only_main=True, custom_path=str(working_directory), used_id="gene", verbose=False)
    annot2.export_CDSs(only_main=True, custom_path=str(working_directory), used_id="gene", verbose=False)

    features_dir = working_directory / "features"
    protein_fasta_1 = features_dir / f"{annot1.name}_proteins_g_id_main.fasta"
    protein_fasta_2 = features_dir / f"{annot2.name}_proteins_g_id_main.fasta"
    cds_fasta_1 = features_dir / f"{annot1.name}_CDSs_g_id_main.fasta"
    cds_fasta_2 = features_dir / f"{annot2.name}_CDSs_g_id_main.fasta"

    # --- 3. RECIPROCAL BEST HIT (RBH) HOMOLOGY SEARCH (DIAMOND) ---
    print("\n[STEP 3] Running reciprocal homology analysis with DIAMOND...")
    
    diamond_db_1 = working_directory / f"{annot1.name}_diamond_db"
    diamond_db_2 = working_directory / f"{annot2.name}_diamond_db"
    diamond_result_1_to_2 = working_directory / f"diamond_{annot1.name}_to_{annot2.name}.txt"
    diamond_result_2_to_1 = working_directory / f"diamond_{annot2.name}_to_{annot1.name}.txt"

    print(f"  Creating DIAMOND database for {annot1.name}...")
    makedb_cmd_1 = [
        "diamond", "makedb", "-p", str(num_threads), "--in", str(protein_fasta_1), "-d", str(diamond_db_1)
    ]
    run_command(working_directory, makedb_cmd_1)

    print(f"  Creating DIAMOND database for {annot2.name}...")
    makedb_cmd_2 = [
        "diamond", "makedb", "-p", str(num_threads), "--in", str(protein_fasta_2), "-d", str(diamond_db_2)
    ]
    run_command(working_directory, makedb_cmd_2)

    print(f"  Running DIAMOND search ({annot1.name} -> {annot2.name})...")
    blastp_cmd_1 = [
        "diamond", "blastp", "--threads", str(num_threads), "--db", str(diamond_db_2), "--ultra-sensitive", 
        "--out", str(diamond_result_1_to_2), "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", 
        "qlen", "slen", "length", "bitscore", "evalue", "--query", str(protein_fasta_1), 
        "--max-target-seqs", "1", "--evalue", "0.00001", "--max-hsps", "1"
    ]
    run_command(working_directory, blastp_cmd_1)

    print(f"  Running DIAMOND search ({annot2.name} -> {annot1.name})...")
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
    jcvi_name_1 = f"{annot1.name}_clean"
    jcvi_name_2 = f"{annot2.name}_clean"
    
    cleaned_cds_1 = working_directory / f"{jcvi_name_1}.cds"
    cleaned_cds_2 = working_directory / f"{jcvi_name_2}.cds"
    bed_file_1 = working_directory / f"{jcvi_name_1}.bed"
    bed_file_2 = working_directory / f"{jcvi_name_2}.bed"
    
    # 4.1 Format CDS files for JCVI
    jcvi_format_cmd_1 = ["python", "-m", "jcvi.formats.fasta", "format", str(cds_fasta_1), str(cleaned_cds_1)]
    run_command(working_directory, jcvi_format_cmd_1)
    
    jcvi_format_cmd_2 = ["python", "-m", "jcvi.formats.fasta", "format", str(cds_fasta_2), str(cleaned_cds_2)]
    run_command(working_directory, jcvi_format_cmd_2)



    #os.chdir(working_directory)

    gff_to_bed_cmd_1 = [
        "python", "-m", "jcvi.formats.gff", "bed", "--type=mRNA",
        "--key=Parent", "--primary_only", f"{annot1.name}_clean.gff3", "-o", str(bed_file_1)
    ]

    print(gff_to_bed_cmd_1)
    run_command(working_directory, gff_to_bed_cmd_1)

    gff_to_bed_cmd_2 = [
        "python", "-m", "jcvi.formats.gff", "bed", "--type=mRNA",
        "--key=Parent", "--primary_only", f"{annot2.name}_clean.gff3", "-o", str(bed_file_2)
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
    
    print(f"\n--- PAIRWISE ANALYSIS COMPLETE: {annot1.name} vs {annot2.name} ---")

    ### CAMBIO ###
    # Devuelve las rutas a los ficheros de proteínas generados.
    return protein_fasta_1, protein_fasta_2

def parse_evalue(e):
    e = e.strip()
    try:
        if e.startswith('>'):
            return float(round_evalue(e[1:]))  # handle cases like '>10'
        elif e.lower() in ('na', 'nan', ''):
            return float('nan')  # handle missing values
        else:
            return float(round_evalue(e))
            
    except ValueError:
        print(f"Error: Could not parse E-value: {e}")
        return float('nan')  # or raise error, depending on use case

def round_evalue(e):
    e = float(e)
    e = f"{e:.2e}"

    return e

class Equivalence():

    preferred_type_order = ["rec_liftoff_aegis", "rec_lifton_aegis", "fwd_liftoff_aegis", "rev_liftoff_aegis", "rev_lifton_aegis", "rev_lifton_aegis", "mcscan_anchors", "mcscan_last_filtered", "rbbh", "rbh", "orthofinder", "fwd_blastp", "rev_blastp", "fwd_blast", "rev_blast"]
    reliability_order = ["vvvtop_reliable", "vvtop_reliable", "vtop_reliable", "top_reliable", "vvvvv_reliable", "vvvv_reliable", "vvv_reliable", "vv_reliable", "v_reliable", "reliable", "NA"]

    def __init__(self, id_, type_, target_annotation, species, score:str="", evalue:str=None, reliability:str="NA"):
        self.id = id_
        self.type = type_
        self.species = species
        self.score = score
        self.target_ = target_annotation
        self.reliability = reliability

        self.aegis_score = None
        self.aegis_synteny = None
        self.aegis_copies = None
        self.match_num = None

        self.aegis_score_rev = None
        self.aegis_synteny_rev = None
        self.aegis_copies_rev = None
        self.match_num_rev = None

        self.coverage = None
        self.identity = None

        self.coverage_rev = None
        self.identity_rev = None

        if evalue != None:
            self.evalue = parse_evalue(evalue)
        else:
            self.evalue = None


        if "confidence" not in self.score:

            if "blast" in self.type or "rbh" in self.type or "rbbh" in self.type:
                for score_component in self.score.split(", "):
                    if score_component.startswith("coverage="):
                        self.identity = round(float(score_component.split("=")[-1].split("/")[0]), 1)
                    elif score_component.startswith("identity="):
                        self.coverage = round(float(score_component.split("=")[-1].split("/")[0]), 1)

                if "rbh" in self.type or "rbbh" in self.type:
                    for score_component in self.score.split(", "):
                        if score_component.startswith("coverage="):
                            self.identity_rev = round(float(score_component.split("=")[-1].split("/")[1]), 1)
                        elif score_component.startswith("identity="):
                            self.coverage_rev = round(float(score_component.split("=")[-1].split("/")[1]), 1)
                
            elif "aegis" in self.type:

                scores = self.score.split("/")

                self.aegis_score = int(scores[0].split(" (")[0])

                if "copies" in scores[0]:
                    self.aegis_copies = True
                else:
                    self.aegis_copies = False

                if "multiple" in scores[0]:
                    self.match_num = int(scores[0].split("multiple ")[-1].split(")")[0])
                else:
                    self.match_num = 1            

                if "synteny" in scores[0]:
                    self.aegis_synteny = True
                else:
                    self.aegis_synteny = False

                if "rec_" in self.type:

                    self.aegis_score_rev = int(scores[1].split(" (")[0])

                    if "copies" in scores[1]:
                        self.aegis_copies_rev = True
                    else:
                        self.aegis_copies_rev = False

                    if "multiple" in scores[1]:
                        self.match_num_rev = int(scores[1].split("multiple ")[-1].split(")")[0])
                    else:
                        self.match_num_rev = 1            

                    if "synteny" in scores[1]:
                        self.aegis_synteny_rev = True
                    else:
                        self.aegis_synteny_rev = False

        else:
            for master_score_component in self.type.split("], "):
                if master_score_component.endswith("]"):
                    master_score_component = master_score_component[:-1]

                if "blast" in master_score_component or "rbh" in master_score_component or "rbbh" in master_score_component:
                    for score_component in master_score_component.split("[")[1].split(", "):
                        if score_component.startswith("coverage="):
                            self.identity = round(float(score_component.split("=")[-1].split("/")[0]), 1)
                        elif score_component.startswith("identity="):
                            self.coverage = round(float(score_component.split("=")[-1].split("/")[0]), 1)

                    if "rbh" in master_score_component or "rbbh" in master_score_component:
                        for score_component in master_score_component.split("[")[1].split(", "):
                            if score_component.startswith("coverage="):
                                self.identity_rev = round(float(score_component.split("=")[-1].split("/")[1]), 1)
                            elif score_component.startswith("identity="):
                                self.coverage_rev = round(float(score_component.split("=")[-1].split("/")[1]), 1)

                elif "aegis" in master_score_component:

                    if "liftoff" in self.type and "lifton" in master_score_component:
                        continue

                    scores = master_score_component.split("[")[1].split("/")

                    self.aegis_score = int(scores[0].split(" (")[0])

                    if "copies" in scores[0]:
                        self.aegis_copies = True
                    else:
                        self.aegis_copies = False

                    if "multiple" in scores[0]:
                        self.match_num = int(scores[0].split("multiple ")[-1].split(")")[0])
                    else:
                        self.match_num = 1            

                    if "synteny" in scores[0]:
                        self.aegis_synteny = True
                    else:
                        self.aegis_synteny = False

                    if "rec_" in master_score_component:

                        self.aegis_score_rev = int(scores[1].split(" (")[0])

                        if "copies" in scores[1]:
                            self.aegis_copies_rev = True
                        else:
                            self.aegis_copies_rev = False

                        if "multiple" in scores[1]:
                            self.match_num_rev = int(scores[1].split("multiple ")[-1].split(")")[0])
                        else:
                            self.match_num_rev = 1            

                        if "synteny" in scores[1]:
                            self.aegis_synteny_rev = True
                        else:
                            self.aegis_synteny_rev = False

    def _rank(self):
        if self.reliability == "NA":
            target_rank = self.target_annotation.lower()
            type_rank = self.preferred_type_order.index(self.type) if self.type in self.preferred_type_order else float('inf')
            evalue = self.evalue if isinstance(self.evalue, float) else float('inf')
            aegis_score = self.aegis_score if self.aegis_score != None else 100
            aegis_score_rev = self.aegis_score_rev if self.aegis_score_rev != None else 100

            if self.aegis_synteny == None:
                synteny = 3
            elif self.aegis_synteny:
                synteny = 1
            else:
                synteny = 2

            if self.aegis_synteny_rev == None:
                synteny_rev = 3
            elif self.aegis_synteny_rev:
                synteny_rev = 1
            else:
                synteny_rev = 2

            if self.aegis_copies == None:
                copies = 3
            elif self.aegis_copies:
                copies = 1
            else:
                copies = 2

            if self.aegis_copies_rev == None:
                copies_rev = 3
            elif self.aegis_copies_rev:
                copies_rev = 1
            else:
                copies_rev = 2

            return (target_rank, type_rank, evalue, -aegis_score, -aegis_score_rev, synteny, synteny_rev, copies, copies_rev)
        else:
            target_rank = self.target_annotation.lower()
            reliability_rank = self.reliability_order.index(self.reliability) if self.reliability in self.reliability_order else float('inf')
            evalue = self.evalue if isinstance(self.evalue, float) else float('inf')
            aegis_score = self.aegis_score if self.aegis_score != None else 100
            aegis_score_rev = self.aegis_score_rev if self.aegis_score_rev != None else 100

            if self.aegis_synteny == None:
                synteny = 3
            elif self.aegis_synteny:
                synteny = 1
            else:
                synteny = 2

            if self.aegis_synteny_rev == None:
                synteny_rev = 3
            elif self.aegis_synteny_rev:
                synteny_rev = 1
            else:
                synteny_rev = 2

            if self.aegis_copies == None:
                copies = 3
            elif self.aegis_copies:
                copies = 1
            else:
                copies = 2

            if self.aegis_copies_rev == None:
                copies_rev = 3
            elif self.aegis_copies_rev:
                copies_rev = 1
            else:
                copies_rev = 2
            
            return (target_rank, reliability_rank, evalue, -aegis_score, -aegis_score_rev, synteny, synteny_rev, copies, copies_rev)

    def __lt__(self, other):
        return self._rank() < other._rank()

    def __eq__(self, other):
        return self._rank() == other._rank()

    def __str__(self):
        return f"{self.id}\t{self.type}\t{self.target_annotation}\t{self.species}\t{self.score}"

    def __repr__(self):
        return f"{self.id}\t{self.type}\t{self.target_annotation}\t{self.species}\t{self.score}"

    def verbose(self):
        return f"{self.id}\t{self.type}\t{self.target_annotation}\t{self.species}\t{self.score}\t{self.reliability}"

class Simple_gene():
    def __init__(self, id):
        self.id = id
        self.equivalences = []
        self.filtered_equivalences = []

    def filter_equivalences(self, simple_rbh_blasts:bool=True, unidirectional_blasts:bool=True, replace:bool=True, identity_threshold=0, coverage_threshold=0, evalue_threshold=float('inf')):
        """
        For now thresholds only apply to rbhs and unidirectional_blasts.
        """

        if self.filtered_equivalences == [] or replace:

            self.filtered_equivalences = []

            equivalent_genes_temp_d = {}

            for equivalence in self.equivalences:
                if equivalence.target_annotation not in equivalent_genes_temp_d:
                    equivalent_genes_temp_d[equivalence.target_annotation] = []
            
            for equivalence in self.equivalences:
                if equivalence.id not in equivalent_genes_temp_d[equivalence.target_annotation]:
                    equivalent_genes_temp_d[equivalence.target_annotation].append(equivalence.id)

            for target_annotation, equivalent_genes in equivalent_genes_temp_d.items():
                for equivalent_gene in equivalent_genes:

                    evidences = []
                    liftoff_score = 0
                    lifton_score = 0

                    lowest_evalue = float('inf')
                    found_evalue = False
                    species = ""

                    for equivalence in self.equivalences:
                        if equivalence.id == equivalent_gene:
                            species = equivalence.species
                            if "aegis" in equivalence.type:
                                if "rec_" in equivalence.type:
                                    score = max(equivalence.aegis_score, equivalence.aegis_score_rev)
                                else:
                                    score = equivalence.aegis_score
                                if "lifton" in equivalence.type:
                                    if score > lifton_score:
                                        lifton_score = score
                                else:
                                    if score > liftoff_score:
                                        liftoff_score = score
                            elif "rbbh" in equivalence.type or "rbh" in equivalence.type or "blast" in equivalence.type:

                                if not unidirectional_blasts and "blast" in equivalence.type:
                                    continue

                                if not simple_rbh_blasts and "rbh" in equivalence.type:
                                    continue

                                if "blast" in equivalence.type or "rbh" in equivalence.type:
                                    if equivalence.identity < identity_threshold:
                                        continue
                                    if equivalence.coverage < coverage_threshold:
                                        continue
                                    if equivalence.evalue > evalue_threshold:
                                        continue

                                if equivalence.evalue < lowest_evalue:
                                    lowest_evalue = equivalence.evalue
                                found_evalue = True

                            evidences.append(f"{equivalence.type} [{equivalence.score}]")

                    summary_evidence = ", ".join(evidences)

                    if (liftoff_score >= 9 or lifton_score >= 9) and "rbbh" in summary_evidence and "mcscan_anchors" in summary_evidence:
                        reliability = "vvvtop_reliable"
                    elif liftoff_score >= 9 or lifton_score >= 9:
                        reliability = "vvtop_reliable"

                    elif liftoff_score >= 7 or lifton_score >= 7 and "rbbh" in summary_evidence and "mcscan_anchors" in summary_evidence:
                        reliability = "vtop_reliable"
                    elif liftoff_score >= 7 or lifton_score >= 7:
                        reliability = "top_reliable"

                    elif liftoff_score >= 6 or lifton_score >= 6 and "rbbh" in summary_evidence and "mcscan_anchors" in summary_evidence:
                        reliability = "vvvvv_reliable"
                    elif "rbbh" in summary_evidence or (liftoff_score >= 6 or lifton_score >= 6 and "rbh" in summary_evidence and "mcscan_anchors" in summary_evidence):
                        reliability = "vvvv_reliable"
                    elif liftoff_score >= 6 or lifton_score >= 6 and "rbh" in summary_evidence in summary_evidence:
                        reliability = "vvv_reliable"

                    elif "mcscan_anchors" in summary_evidence and "rbh" in summary_evidence:
                        reliability = "vv_reliable"

                    elif "mcscan_anchors" in summary_evidence or "rbh" in summary_evidence:
                        reliability = "v_reliable"                   

                    else:
                        reliability = "reliable"

                    if found_evalue:
                        evalue = str(lowest_evalue)
                    else:
                        evalue = "NA"

                    if "top_reliable" in reliability:
                        score = "high_confidence"
                    elif "vv_reliable" in reliability:
                        score = "medium_confidence"
                    else:
                        score = "lower_confidence"

                    if evidences != []:

                        self.filtered_equivalences.append(Equivalence(equivalent_gene, summary_evidence, target_annotation, species, score, evalue, reliability))


class Simple_annotation():
    def __init__(self, name, annotation_object:object, species:str):
        self.name = name
        self.genes = {}
        self.species = species
        self.target = annotation_object.target
        for gene in annotation_object.all_gene_ids:
            if gene not in self.genes:
                self.genes[gene] = Simple_gene(gene)
            else:
                print(f"Warning: repeated gene id: {gene}")

        self.added_equivalences = {}

    def export_summary_equivalences(self, output_file, filtered:bool=False, simple_rbh_blasts:bool=True, unidirectional_blasts:bool=True, replace:bool=True, identity_threshold=0, coverage_threshold=0, evalue_threshold=float('inf'), verbose:bool=True):
        
        start = time.time()
        
        out = ["gene\tequivalence\ttype\ttarget_annotation\ttarget_species\tscore"]
        if verbose:
            out[0] += "\treliability"

        for gene in self.genes.values():
            gene.equivalences = sorted(gene.equivalences)
        
        if filtered:
            self.filter_equivalences(simple_rbh_blasts=simple_rbh_blasts, unidirectional_blasts=unidirectional_blasts, replace=replace, identity_threshold=identity_threshold, coverage_threshold=coverage_threshold, evalue_threshold=evalue_threshold)

            for gene in self.genes.values():
                gene.filtered_equivalences = sorted(gene.filtered_equivalences)
                for equivalence in gene.filtered_equivalences:
                    if verbose:
                        out.append(f"{gene.id}\t{equivalence.verbose()}")
                    else:
                        out.append(f"{gene.id}\t{equivalence}")
        else:
            for gene in self.genes.values():
                for equivalence in gene.equivalences:
                    if verbose:
                        out.append(f"{gene.id}\t{equivalence.verbose()}")
                    else:
                        out.append(f"{gene.id}\t{equivalence}")

        out = "\n".join(out)

        f_out = open(output_file, "w", encoding="utf-8")
        f_out.write(out)
        f_out.close()

        end = time.time()
        lapse = end - start
        print(f"Exporting all filtered={filtered} equivalences for {self.name} took {round(lapse/60, 2)} minutes\n")

    def add_mcscan_equivalences(self, file, key_col, target_annotation, species):
        """
        Adds MCScan equivalences from ".anchors" or ".last.filtered"
        """

        if os.path.isfile(file):

            df = pd.read_csv(file, comment='#', header=None, sep="\t", dtype=str)

            if key_col == "0":
                target_col = "1"
            else:
                target_col = "0"

            if file.endswith(".filtered"):
                equivalence_type = "mcscan_last_filtered"
                df.columns = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11"]
                score_col = "11"
            else:
                equivalence_type = "mcscan_anchors"
                df.columns = ["0", "1", "2"]
                score_col = "2"

            for index, row in df.iterrows():
                gene_query = row[key_col]
                gene_target = row[target_col]
                score = row[score_col]
                self.genes[gene_query].equivalences.append(Equivalence(gene_target, equivalence_type, target_annotation, species, score))


    def add_orthofinder_equivalences(self, file, target_annotation, species):

        if os.path.isfile(file):
            df = pd.read_csv(file, sep="\t", dtype=str)
            for index, row in df.iterrows():
                gene_queries = row[file.split("/")[-1].split(".tsv")[0]].split(", ")
                gene_targets = row["Orthologs"].split(", ")
                score = row["Orthogroup"]
                for gene_query in gene_queries:
                    for gene_target in gene_targets:
                        self.genes[gene_query].equivalences.append(Equivalence(gene_target, "orthofinder", target_annotation, species, score))

    def add_blast_equivalences(self, blast_folder, query_annotation, target_annotation, species, skip_rbhs:bool=False, skip_unidirectional_blasts:bool=False, proteins:bool=True):

        if not skip_unidirectional_blasts and skip_rbhs:
            print("Warning: It makes no sense to skip RBHs when unidirectional blasts are not being skipped.")

        else:

            start = time.time()

            query_col = "query_x"
            target_col = "subject_x"
            query_col_evalue = "E-value_x"
            target_col_evalue = "E-value_y"
            query_col_normbit = "norm_bitscore_x"
            target_col_normbit = "norm_bitscore_y"
            query_col_coverage = "coverage_x"
            target_col_coverage = "coverage_y"
            query_col_identity = "identity_x"
            target_col_identity = "identity_y"

            rbbh_file = f"{blast_folder}/rbh/rbbh_{query_annotation}_against_{target_annotation}.csv"
            rbh_file = f"{blast_folder}/rbh/rbh_{query_annotation}_against_{target_annotation}.csv"
            if not os.path.isfile(rbbh_file):
                query_col = "subject_x"
                target_col = "query_x"
                query_col_evalue = "E-value_y"
                target_col_evalue = "E-value_x"
                query_col_normbit = "norm_bitscore_y"
                target_col_normbit = "norm_bitscore_x"
                query_col_coverage = "coverage_y"
                target_col_coverage = "coverage_x"
                query_col_identity = "identity_y"
                target_col_identity = "identity_x"
                rbbh_file = f"{blast_folder}/rbh/rbbh_{target_annotation}_against_{query_annotation}.csv"
                rbh_file = f"{blast_folder}/rbh/rbh_{target_annotation}_against_{query_annotation}.csv"

            fwd_file = f"{blast_folder}fwd_{query_annotation}_against_{target_annotation}.csv"
            if not os.path.isfile(fwd_file):
                fwd_file = f"{blast_folder}rev_{query_annotation}_against_{target_annotation}.csv"

            rev_file = f"{blast_folder}rev_{target_annotation}_against_{query_annotation}.csv"
            if not os.path.isfile(rev_file):
                rev_file = f"{blast_folder}fwd_{target_annotation}_against_{query_annotation}.csv"

            rbbh_hits = set()

            target_annotation = clean_annotation_tag(target_annotation)

            temp_start = time.time()

            if os.path.isfile(rbbh_file):
                df = pd.read_csv(rbbh_file, sep="\t", dtype=str)
                for index, row in df.iterrows():
                    gene_query = row[query_col]
                    gene_target = row[target_col]
                    evalue = row[query_col_evalue]
                    score = f"evalue={row[query_col_evalue]}/{row[target_col_evalue]}"
                    score += f", norm_bitscore={round(float(row[query_col_normbit]), 1)}/{round(float(row[target_col_normbit]), 1)}"
                    score += f", coverage={round(float(row[query_col_coverage]), 1)}/{round(float(row[target_col_coverage]), 1)}"
                    score += f", identity={round(float(row[query_col_identity]), 1)}/{round(float(row[target_col_identity]), 1)}"
                    pair = sorted([gene_query, gene_target])
                    pair = "-".join(pair)
                    rbbh_hits.add(pair)
                    self.genes[gene_query].equivalences.append(Equivalence(gene_target, "rbbh", target_annotation, species, score, evalue))

            end = time.time()
            lapse = end - temp_start
            print(f"Adding RBBH equivalences to {self.name} '{query_annotation} vs {target_annotation}' took {round(lapse/60, 2)} minutes")

            if not skip_rbhs:

                temp_start = time.time()

                rbh_hits = set()

                if os.path.isfile(rbh_file):
                    df = pd.read_csv(rbh_file, sep="\t", dtype=str)
                    for index, row in df.iterrows():
                        gene_query = row[query_col]
                        gene_target = row[target_col]
                        evalue = row[query_col_evalue]
                        score = f"evalue={row[query_col_evalue]}/{row[target_col_evalue]}"
                        score += f", norm_bitscore={round(float(row[query_col_normbit]), 1)}/{round(float(row[target_col_normbit]), 1)}"
                        score += f", coverage={round(float(row[query_col_coverage]), 1)}/{round(float(row[target_col_coverage]), 1)}"
                        score += f", identity={round(float(row[query_col_identity]), 1)}/{round(float(row[target_col_identity]), 1)}"
                        pair = sorted([gene_query, gene_target])
                        pair = "-".join(pair)
                        if pair not in rbbh_hits:
                            rbh_hits.add(pair)
                            self.genes[gene_query].equivalences.append(Equivalence(gene_target, "rbh", target_annotation, species, score, evalue))

                end = time.time()
                lapse = end - temp_start
                print(f"Adding RBH equivalences to {self.name} '{query_annotation} vs {target_annotation}' took {round(lapse/60, 2)} minutes")

                if not skip_unidirectional_blasts:

                    temp_start = time.time()

                    if proteins:
                        equivalence_type = "fwd_blastp"
                    else:
                        equivalence_type = "fwd_blast"

                    if os.path.isfile(fwd_file):
                        df = pd.read_csv(fwd_file, sep="\t", dtype=str, header=None)
                        df.columns = ["0", "1", "2", "3", "4", "5", "6", "7", "8"]
                        for index, row in df.iterrows():
                            gene_query = row["0"]
                            gene_target = row["1"]
                            norm_bitscore = round((float(row["7"]) / float(row["4"])), 1)
                            evalue = row['8']
                            score = f"evalue={evalue}"
                            score += f", norm_bitscore={norm_bitscore}"
                            score += f", coverage={round(float(row['3']), 1)}"
                            score += f", identity={round(float(row['2']), 1)}"
                            pair = sorted([gene_query, gene_target])
                            pair = "-".join(pair)
                            if pair not in rbh_hits and pair not in rbbh_hits:
                                self.genes[gene_query].equivalences.append(Equivalence(gene_target, equivalence_type, target_annotation, species, score, evalue))

                    if proteins:
                        equivalence_type = "rev_blastp"
                    else:
                        equivalence_type = "rev_blast"

                    if os.path.isfile(rev_file):
                        df = pd.read_csv(rev_file, sep="\t", dtype=str, header=None)
                        df.columns = ["0", "1", "2", "3", "4", "5", "6", "7", "8"]
                        for index, row in df.iterrows():
                            gene_query = row["1"]
                            gene_target = row["0"]
                            norm_bitscore = round((float(row["7"]) / float(row["4"])), 1)
                            evalue = row['8']
                            score = f"evalue={evalue}"
                            score += f", norm_bitscore={norm_bitscore}"
                            score += f", coverage={round(float(row['3']), 1)}"
                            score += f", identity={round(float(row['2']), 1)}"
                            pair = sorted([gene_query, gene_target])
                            pair = "-".join(pair)
                            if pair not in rbh_hits and pair not in rbbh_hits:
                                self.genes[gene_query].equivalences.append(Equivalence(gene_target, equivalence_type, target_annotation, species, score, evalue))

                    end = time.time()
                    lapse = end - temp_start
                    print(f"Adding fwd/rev blast equivalences to {self.name} '{query_annotation} vs {target_annotation}' took {round(lapse/60, 2)} minutes")

            end = time.time()
            lapse = end - start
            
            if skip_unidirectional_blasts and not skip_rbhs:
                print(f"Adding RBBH, and RBH equivalences to {self.name} '{query_annotation} vs {target_annotation}' took {round(lapse/60, 2)} minutes")
            else:
                print(f"Adding RBBH, RBH, and fwd/rev blast equivalences to {self.name} '{query_annotation} vs {target_annotation}' took {round(lapse/60, 2)} minutes")

    def add_reciprocal_overlap_equivalences(self, folder, query_tag, target_tag, species, liftoff:bool=True):

        program = "Liftoff"
        if not liftoff:
            program = "Lifton"

        fwd_file = f"{folder}{query_tag}_equivalences_{program}_on_{target_tag}.tsv"
        rev_file = f"{folder}{target_tag}_equivalences_{program}_on_{query_tag}.tsv"

        if os.path.isfile(fwd_file) and os.path.isfile(rev_file):
            start = time.time()

            go_ahead = True

            fwd_df = pd.read_csv(fwd_file, sep="\t", encoding="utf-8", dtype=str, na_filter=False)

            duplicates = fwd_df[fwd_df.duplicated(subset=["query_id", "target_id"], keep=False)]

            if not duplicates.empty:
                print(f"Error: Duplicate query and target id pairs for fwd {program}:")
                print(duplicates)
                go_ahead = False

            rev_df = pd.read_csv(rev_file, sep="\t", encoding="utf-8", dtype=str, na_filter=False)

            duplicates = rev_df[rev_df.duplicated(subset=["query_id", "target_id"], keep=False)]

            if not duplicates.empty:
                print(f"Error: Duplicate query and target id pairs for rev {program}:")
                print(duplicates)
                go_ahead = False

            if go_ahead:

                fwd_df = fwd_df[fwd_df["query_id"] != "NA"]
                rev_df = rev_df[rev_df["query_id"] != "NA"]

                fwd_df = fwd_df[fwd_df["target_id"] != "NA"]
                rev_df = rev_df[rev_df["target_id"] != "NA"]

                fwd_df["query_origin"] = fwd_df["query_origin"].apply(clean_annotation_tag)

                rev_df["target_origin"] = rev_df["target_origin"].apply(clean_annotation_tag)

                fwd_df = fwd_df[fwd_df["query_origin"] == self.name]
                rev_df = rev_df[rev_df["target_origin"] == self.name]

                if liftoff:
                    equivalence_suffix = "liftoff_aegis"
                else:
                    equivalence_suffix = "lifton_aegis"

                fwd_hits = {}
            
                for index, row in fwd_df.iterrows():

                    gene_query = row["query_id"]
                    gene_target = row["target_id"]

                    copies = False
                    if re.search(r"_\d{1,4}$", gene_query):
                        copies = True
                        gene_query = gene_query.split("_")[:-1]
                        gene_query = "_".join(gene_query)

                    if re.search(r"_\d{1,4}$", gene_target):
                        copies = True
                        gene_target = gene_target.split("_")[:-1]
                        gene_target = "_".join(gene_target)

                    pair = f"{gene_query}---{gene_target}"

                    synteny = True
                    if row["query_synteny_conserved"] == "False" or row["target_synteny_conserved"] == "False":
                        synteny = False
                    score = row["overlap_score"]

                    if synteny and copies:
                        score = f"{score} (synteny and copies)"
                    elif synteny:
                        score = f"{score} (synteny)"
                    elif copies:
                        score = f"{score} (copies)"

                    if pair not in fwd_hits:
                        fwd_hits[pair] = score
                    else:

                        old_score = int(fwd_hits[pair].split(" (")[0])
                        new_score = int(score.split(" (")[0])

                        if old_score > new_score:
                            new_score = old_score

                        if "multiple " in fwd_hits[pair]:
                            match_num = int(fwd_hits[pair].split("multiple ")[-1].split(")")[0])
                            match_num += 1
                        else:
                            match_num = 2

                        synteny = False
                        if "synteny" in fwd_hits[pair] or "synteny" in score:
                            synteny = True

                        copies = False
                        if "copies" in fwd_hits[pair] or "copies" in score:
                            copies = True
                        
                        if copies and synteny:
                            fwd_hits[pair] = f"{new_score} (synteny and copies) (multiple {match_num})"
                        elif copies:
                            fwd_hits[pair] = f"{new_score} (copies) (multiple {match_num})"
                        else:
                            raise ValueError(f"{pair} with score {score} should not have been already added to fwd_hits dictionary since at least one extra match should have been found through copies. The previous pair entry was already added with score {fwd_hits[pair]}. Query={query_tag} and Target={target_tag}.")

                rev_hits = {}
            
                for index, row in rev_df.iterrows():

                    gene_query = row["target_id"]
                    gene_target = row["query_id"]

                    copies = False
                    if re.search(r"_\d{1,4}$", gene_query):
                        copies = True
                        gene_query = gene_query.split("_")[:-1]
                        gene_query = "_".join(gene_query)

                    if re.search(r"_\d{1,4}$", gene_target):
                        copies = True
                        gene_target = gene_target.split("_")[:-1]
                        gene_target = "_".join(gene_target)

                    pair = f"{gene_query}---{gene_target}"

                    synteny = True
                    if row["target_synteny_conserved"] == "False" or row["query_synteny_conserved"] == "False":
                        synteny = False
                    score = row["overlap_score"]

                    if synteny and copies:
                        score = f"{score} (synteny and copies)"
                    elif synteny:
                        score = f"{score} (synteny)"
                    elif copies:
                        score = f"{score} (copies)"

                    if pair not in rev_hits:
                        rev_hits[pair] = score
                    else:

                        old_score = int(rev_hits[pair].split(" (")[0])
                        new_score = int(score.split(" (")[0])

                        if old_score > new_score:
                            new_score = old_score

                        if "multiple " in rev_hits[pair]:
                            match_num = int(rev_hits[pair].split("multiple ")[-1].split(")")[0])
                            match_num += 1
                        else:
                            match_num = 2

                        synteny = False
                        if "synteny" in rev_hits[pair] or "synteny" in score:
                            synteny = True

                        copies = False
                        if "copies" in rev_hits[pair] or "copies" in score:
                            copies = True
                        
                        if copies and synteny:
                            rev_hits[pair] = f"{new_score} (synteny and copies) (multiple {match_num})"
                        elif copies:
                            rev_hits[pair] = f"{new_score} (copies) (multiple {match_num})"
                        else:
                            raise ValueError(f"{pair} with score {score} should not have been already added to rev_hits dictionary since at least one extra match should have been found through copies. The previous pair entry was already added with score {rev_hits[pair]}. Query={query_tag} and Target={target_tag}.")

                for fwd_pair in fwd_hits:

                    if fwd_pair in rev_hits:
                        equivalence_type = f"rec_{equivalence_suffix}"
                        score = f"{fwd_hits[fwd_pair]}/{rev_hits[fwd_pair]}"
                    else:
                        equivalence_type = f"fwd_{equivalence_suffix}"
                        score = f"{fwd_hits[fwd_pair]}"

                    gene_query = fwd_pair.split("---")[0]
                    gene_target = fwd_pair.split("---")[1]

                    self.genes[gene_query].equivalences.append(Equivalence(gene_target, equivalence_type, target_tag, species, score))


                for rev_pair in rev_hits:
                    if rev_pair not in fwd_hits:
                        equivalence_type = f"rev_{equivalence_suffix}"
                        score = f"{rev_hits[rev_pair]}"

                        gene_query = rev_pair.split("---")[0]
                        gene_target = rev_pair.split("---")[1]
                        
                        self.genes[gene_query].equivalences.append(Equivalence(gene_target, equivalence_type, target_tag, species, score))

                end = time.time()
                lapse = end - start
                print(f"Adding {program} overlap equivalences for tags = [{query_tag}, {target_tag}] to {self.name} took {round(lapse/60, 2)} minutes")

        else:
            print(f"Warning: {fwd_file} or {rev_file} is missing.")


    def filter_equivalences(self, simple_rbh_blasts:bool=True, unidirectional_blasts:bool=True, replace:bool=True, identity_threshold=0, coverage_threshold=0, evalue_threshold=float("inf")):

        for gene in self.genes.values():
            gene.filter_equivalences(simple_rbh_blasts=simple_rbh_blasts, unidirectional_blasts=unidirectional_blasts, replace=replace, identity_threshold=identity_threshold, coverage_threshold=coverage_threshold, evalue_threshold=evalue_threshold)


def clean_annotation_tag(annotation_tag):
    annotation_tag = annotation_tag.replace("Lifton", "")

    annotation_tag = annotation_tag.split("_from_")[0]

    annotation_tag = annotation_tag.split("_on_")[0]
    
    return annotation_tag