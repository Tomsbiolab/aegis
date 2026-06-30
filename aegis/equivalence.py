from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .genome import Genome

import pandas as pd
import time
import os
import shutil
import subprocess

from collections import defaultdict
from pathlib import Path

from .annotation import Annotation
from .utils.misc import run_command
from .utils.evalue import parse_evalue

def pairwise_orthology(annot1_name:str, annot2_name: str, annot1_file: str, annot2_file: str, genome1: Genome, genome2: Genome, working_directory: Path, num_threads: int, types: str, evalue:float=0.00001, coverage:float=30, max_hsps:int=1, copies:bool=True, synteny:bool=False, skip_liftoff:bool=False, skip_lifton:bool=False, skip_mcscan:bool=False, skip_blasts:bool=False, pairwise_orthofinder:bool=False, skip_orthofinder:bool=False, quiet:bool=True):

    print(f"\n\tRunning pairwise orthology for {annot1_name} and {annot2_name}...")

    annot1: Annotation
    annot2: Annotation

    annot1 = Annotation(annot_file_path=annot1_file, genome=genome1, name=annot1_name, quiet=quiet, define_synteny=synteny)
    annot2 = Annotation(annot_file_path=annot2_file, genome=genome2, name=annot2_name, quiet=quiet, define_synteny=synteny)

    liftless_overlaps_dir = working_directory / "liftless_overlaps"
    liftoff_dir = working_directory / "liftoff"
    lifton_dir = working_directory / "lifton"
    protein_dir = working_directory / "proteins"
    diamond_dir = working_directory / "diamond"
    mcscan_dir = working_directory / "mcscan"
    orthofinder_dir = protein_dir / f"orthofinder_{annot1.name}__to__{annot2.name}"

    liftoff_pair_dir = liftoff_dir / f"{annot1.name}__to__{annot2.name}"
    lifton_pair_dir = lifton_dir / f"{annot1.name}__to__{annot2.name}"

    if genome1.file == genome2.file:
        skip_mcscan = True
        skip_orthofinder = True
        skip_lifton = True
        skip_liftoff = True
        liftless_overlaps = True
    else:
        liftless_overlaps = False

    if not skip_liftoff:

        liftoff_pair_dir.mkdir(parents=True, exist_ok=True)

        liftoff_overlaps_tsv = liftoff_dir / f"liftoff__{annot1.name}__to__{annot2.name}_overlaps.tsv"
        liftoff_gff = liftoff_pair_dir / f"liftoff__{annot1.name}__to__{annot2.name}.gff"

        if liftoff_overlaps_tsv.exists():
            if not quiet:
                print(f"\n\tExisting Liftoff overlap output found for {annot1.name} on {annot2.name}. Skipping Liftoff and overlaps calculation.")
        else:

            if liftoff_gff.exists():
                if not quiet:
                    print(f"\n\tExisting Liftoff output found for {annot1.name} on {annot2.name}. Reusing liftoff output for aegis overlap calculation.")
            else:
                if not quiet:
                    print(f"\n\tRunning Liftoff to map annotations from {annot1.name} on {annot2.name}")

                shared_gff_path = Path(working_directory) / "gffs" / f"{annot1.name}.gff3"
                local_gff_path = liftoff_pair_dir / f"{annot1.name}.gff3"
                shutil.copy(shared_gff_path, local_gff_path)

                liftoff_cmd = [
                    "liftoff", str(genome2.file), str(genome1.file),
                    "-g", str(local_gff_path), "-o", str(liftoff_gff), "-flank",  "0.1", "-f", types, "-p", str(num_threads)
                ]
                if copies:
                    liftoff_cmd.append("-copies")

                try:
                    run_command(liftoff_pair_dir , liftoff_cmd)
                except subprocess.CalledProcessError:
                    print(f"\n\t⚠️  Liftoff failed to map annotations from {annot1.name} on {annot2.name}.")
                    print("\tThis can happen with highly divergent genomes. Skipping Liftoff mapping...")

                to_remove = liftoff_pair_dir  / "intermediate_files"

                if os.path.exists(str(to_remove)):
                    shutil.rmtree(str(to_remove))

                unmapped_file = f"{str(liftoff_pair_dir)}/unmapped_features.txt"
                if os.path.isfile(unmapped_file):
                    os.remove(unmapped_file)

                if local_gff_path.exists():
                    os.remove(local_gff_path)


                if not os.path.isfile(liftoff_gff):
                    with open(liftoff_gff, "w") as f:
                        f.write("##gff-version 3\n")

            if not quiet:
                print(f"\t\tRunning aegis overlap on liftoff result.")

            if synteny:
                a_liftoff = Annotation(annot_file_path=str(liftoff_gff), name=f"liftoff_{annot1.name}", genome=genome2, original_annotation=annot1, quiet=quiet, define_synteny=synteny)
            else:
                a_liftoff = Annotation(annot_file_path=str(liftoff_gff), name=f"liftoff_{annot1.name}", genome=genome2, quiet=quiet, define_synteny=synteny)

            a_liftoff.overlaps.detect(annot2, quiet=quiet, clear=True)

            _ = a_liftoff.overlaps.export(output_dir=str(liftoff_dir), filename=f"liftoff__{annot1.name}__to__{annot2.name}_overlaps.tsv", verbose=True, save_csv=True, NAs=False, quiet=quiet, synteny=synteny, copies_info=True)

            del a_liftoff

    elif liftless_overlaps:

        liftless_overlaps_tsv = liftless_overlaps_dir / f"{annot1.name}__to__{annot2.name}_overlaps.tsv"

        if liftless_overlaps_tsv.exists():
            if not quiet:
                print(f"\n\tExisting aegis overlap output found for {annot1.name} on {annot2.name}. Skipping.")
        else:
            annot1.overlaps.detect(annot2, quiet=quiet, clear=True)
            _ = annot1.overlaps.export(output_dir=str(liftless_overlaps_dir), filename=f"{annot1.name}__to__{annot2.name}_overlaps.tsv", verbose=True, save_csv=True, NAs=False, quiet=quiet, synteny=False, copies_info=False)

    if not skip_lifton:

        lifton_pair_dir.mkdir(parents=True, exist_ok=True)

        lifton_overlaps_tsv = lifton_dir / f"lifton__{annot1.name}__to__{annot2.name}_overlaps.tsv"
        lifton_gff = lifton_pair_dir / f"lifton__{annot1.name}__to__{annot2.name}.gff3"

        if lifton_overlaps_tsv.exists():
            if not quiet:
                print(f"\n\tExisting LiftOn overlap output found for {annot1.name} on {annot2.name}. Skipping LiftOn and overlaps calculation.")
        else:
            if lifton_gff.exists():
                if not quiet:
                    print(f"\n\tExisting LiftOn output found for {annot1.name} on {annot2.name}. Reusing lifton output for aegis overlap calculation.")
            else:
                if not quiet:
                    print(f"\n\tRunning LiftOn to map annotations from {annot1.name} on {annot2.name}")
                
                shared_lifton_gff_path = Path(working_directory) / "gffs" / f"{annot1.name}__for__lifton.gff3"
                local_lifton_gff_path = lifton_pair_dir / f"{annot1.name}__for__lifton.gff3"
                shutil.copy(shared_lifton_gff_path, local_lifton_gff_path)

                lifton_cmd = [
                    "lifton", "-g", str(local_lifton_gff_path), "-o", str(lifton_gff),
                    "-flank",  "0.1", "-f", types, "-t", str(num_threads)
                ]
                if copies:
                    lifton_cmd.append("-copies")

                lifton_cmd.append(str(genome2.file))
                lifton_cmd.append(str(genome1.file))

                try:
                    run_command(lifton_pair_dir, lifton_cmd)
                except subprocess.CalledProcessError:
                    print(f"\n\t⚠️  LiftOn failed to map annotations from {annot1.name} on {annot2.name}.")
                    print("\tThis can happen with highly divergent genomes. Skipping LiftOn mapping...")

                to_remove = lifton_pair_dir / "lifton_output"

                if os.path.exists(str(to_remove)):
                    shutil.rmtree(str(to_remove))

                if local_lifton_gff_path.exists():
                    os.remove(local_lifton_gff_path)


                if not os.path.isfile(lifton_gff):
                    with open(lifton_gff, "w") as f:
                        f.write("##gff-version 3\n")

            if not quiet:
                print(f"\t\tRunning aegis overlap on lifton result.")

            if synteny:
                a_lifton = Annotation(annot_file_path=str(lifton_gff), name=f"lifton_{annot1.name}", genome=genome2, original_annotation=annot1, quiet=quiet, define_synteny=synteny)
            else:
                a_lifton = Annotation(annot_file_path=str(lifton_gff), name=f"lifton_{annot1.name}", genome=genome2, quiet=quiet, define_synteny=synteny)

            a_lifton.overlaps.detect(annot2, quiet=quiet, clear=True)

            _ = a_lifton.overlaps.export(output_dir=str(lifton_dir), filename=f"lifton__{annot1.name}__to__{annot2.name}_overlaps.tsv", verbose=True, save_csv=True, NAs=False, quiet=quiet, synteny=synteny, copies_info=True)

            del a_lifton

    protein_fasta = protein_dir / f"{annot1.name}_proteins_g_id_main.fasta"

    if not skip_blasts:

        diamond_result = diamond_dir / f"single_{annot1.name}__to__{annot2.name}.txt"
        diamond_result_best = diamond_dir / f"single_best_{annot1.name}__to__{annot2.name}.txt"
        diamond_db = diamond_dir / f"{annot2.name}_diamond_db"

        if diamond_result.exists():
            if not quiet:
                print(f"\n\tExisting standard DIAMOND results found for {annot1.name} on {annot2.name}. Skipping standard search.")
        else:
            if not quiet:
                print(f"\n\tRunning standard DIAMOND search ({annot1.name} on {annot2.name})")
            blastp_cmd = [
                "diamond", "blastp", "--threads", str(num_threads), "--db", str(diamond_db), "--ultra-sensitive", 
                "--out", str(diamond_result), "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", 
                "qlen", "slen", "length", "bitscore", "evalue", "--query", str(protein_fasta), 
                "--evalue", str(evalue), "--max-hsps", str(max_hsps), "--query-cover", str(coverage)
            ]
            run_command(diamond_dir, blastp_cmd)

        if diamond_result_best.exists():
            if not quiet:
                print(f"\n\tExisting 'best' DIAMOND results found for {annot1.name} on {annot2.name}. Skipping 'best' search.")
        else:
            if not quiet:
                print(f"\n\tRunning 'best' DIAMOND search ({annot1.name} on {annot2.name})")
            blastp_cmd = [
                "diamond", "blastp", "--threads", str(num_threads), "--db", str(diamond_db), "--ultra-sensitive", 
                "--out", str(diamond_result_best), "--outfmt", "6", "qseqid", "sseqid", "pident", "qcovhsp", 
                "qlen", "slen", "length", "bitscore", "evalue", "--query", str(protein_fasta), 
                "--max-target-seqs", "1", "--evalue", str(evalue), "--max-hsps", str(max_hsps)
            ]
            run_command(diamond_dir, blastp_cmd)
    
    if not skip_mcscan:
        if not quiet:
            print(f"\n\tRunning JCVI ortholog analysis between {annot1.name} and {annot2.name}")

        mcscan_name1 = annot1.name.replace(".", "_")
        mcscan_name2 = annot2.name.replace(".", "_")

        anchors_file = mcscan_dir / f"{mcscan_name1}.{mcscan_name2}.anchors"
        filtered_file = mcscan_dir / f"{mcscan_name1}.{mcscan_name2}.last.filtered"

        if anchors_file.exists() and filtered_file.exists():
            if not quiet:
                print(f"\n\tExisting MCscan output files found for {annot1.name} and {annot2.name}. Skipping MCscan step.")
        else:
            if not quiet:
                print(f"\n\tRunning JCVI ortholog analysis between {annot1.name} and {annot2.name}")

            pair_tmp_dir = mcscan_dir / f"tmp_{mcscan_name1}_v_{mcscan_name2}"
            pair_tmp_dir.mkdir(parents=True, exist_ok=True)

            try:
                for name in [mcscan_name1, mcscan_name2]:
                    for ext in [".cds", ".bed"]:
                        src = mcscan_dir / f"{name}{ext}"
                        dst = pair_tmp_dir / f"{name}{ext}"
                        if src.exists() and not dst.exists():
                            try:
                                dst.symlink_to(src)
                            except OSError:
                                shutil.copy(src, dst)

                jcvi_ortho_cmd = [
                    "python", "-m", "jcvi.compara.catalog", "ortholog",
                    mcscan_name1, mcscan_name2, "--no_strip_names", "--cpus", str(num_threads)
                ]
                
                run_command(pair_tmp_dir, jcvi_ortho_cmd)

                for ext in [".anchors", ".last", ".last.filtered"]:
                    generated_file = pair_tmp_dir / f"{mcscan_name1}.{mcscan_name2}{ext}"
                    target_file = mcscan_dir / f"{mcscan_name1}.{mcscan_name2}{ext}"
                    if generated_file.exists():
                        shutil.copy(generated_file, target_file)

            except subprocess.CalledProcessError as e:
                error_msg = f"{e.stdout} {e.stderr}"
                if "0 anchor was found" in error_msg:
                    print(f"\n\t⚠️  JCVI found 0 syntenic anchors between {annot1.name} and {annot2.name}.")
                    print("\tThis is a normal biological result. Skipping anchor generation and continuing pipeline...")
                else:
                    raise
            finally:
                if pair_tmp_dir.exists():
                    shutil.rmtree(pair_tmp_dir)

    if pairwise_orthofinder and not skip_orthofinder:


        orthofinder_results_dir = orthofinder_dir / "orthofinder"
        existing_results = list(orthofinder_results_dir.glob("Results*")) if orthofinder_results_dir.exists() else []

        if existing_results:
            if not quiet:
                print(f"\n\tExisting pairwise OrthoFinder results found for {annot1.name} and {annot2.name}. Skipping.")
        else:
            if not quiet:
                print(f"\nRunning OrthoFinder between {annot1.name} and {annot2.name}")

            orthofinder_dir.mkdir(parents=True, exist_ok=True)

            protein_fasta_1 = protein_dir / f"{annot1.name}_proteins_g_id_main.fasta"
            protein_fasta_2 = protein_dir / f"{annot2.name}_proteins_g_id_main.fasta"

            shutil.copy(protein_fasta_1, orthofinder_dir / protein_fasta_1.name)
            if protein_fasta_1 != protein_fasta_2:
                shutil.copy(protein_fasta_2, orthofinder_dir / protein_fasta_2.name)

            orthofinder_cmd = [
                "orthofinder",
                "-f", str(orthofinder_dir),
                "-t", str(num_threads),
                "-a", str(num_threads),
                "-o", f"{str(orthofinder_dir)}/orthofinder/"
            ]
            run_command(orthofinder_dir, orthofinder_cmd)

class Equivalence():

    __slots__ = (
        'id', 'type', 'species', 'score', 'target_annotation', 'reliability',
        'aegis_score', 'aegis_synteny', 'aegis_copies', 'match_num',
        'aegis_score_rev', 'aegis_synteny_rev', 'aegis_copies_rev', 'match_num_rev',
        'coverage', 'identity', 'coverage_rev', 'identity_rev', 'evalue'
    )

    preferred_type_order = ["aegis_overlap", "rec_liftoff_aegis", "rec_lifton_aegis", "fwd_liftoff_aegis", "rev_liftoff_aegis", "fwd_lifton_aegis", "rev_lifton_aegis", "mcscan_anchors", "mcscan_last_filtered", "rbbh", "rbh", "orthofinder", "fwd_blastp", "rev_blastp", "fwd_blast", "rev_blast"]
    reliability_order = ["vvvtop_reliable", "vvtop_reliable", "vtop_reliable", "top_reliable", "vvvvv_reliable", "vvvv_reliable", "vvv_reliable", "vv_reliable", "v_reliable", "reliable", "NA"]

    def __init__(self, id_, type_, target_annotation, species, score:str="", evalue:str|None=None, reliability:str="NA"):
        self.id = id_
        self.type = type_
        self.species = species
        self.score = score
        self.target_annotation = target_annotation
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
        return f"{self.id}\t{self.type}\t{self.score}\t{self.target_annotation}\t{self.species}"

    def __repr__(self):
        return f"{self.id}\t{self.type}\t{self.score}\t{self.target_annotation}\t{self.species}"

    def verbose(self):
        return f"{self.id}\t{self.type}\t{self.score}\t{self.target_annotation}\t{self.species}\t{self.reliability}"

class Simple_gene():
    __slots__ = ('id', 'equivalences', 'filtered_equivalences')
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

            grouped_eqs = defaultdict(list)

            for equivalence in self.equivalences:
                grouped_eqs[(equivalence.target_annotation, equivalence.id)].append(equivalence)

            for (target_annotation, equivalent_gene), eqs in grouped_eqs.items():
                evidences = []
                liftoff_score = 0
                lifton_score = 0

                lowest_evalue = float('inf')
                found_evalue = False
                species = ""

                for equivalence in eqs:
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
    def __init__(self, name, annotation_object:Annotation, species:str):
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

    def export_summary_equivalences(self, output_file, filtered:bool=False, simple_rbh_blasts:bool=True, unidirectional_blasts:bool=True, replace:bool=True, identity_threshold=0, coverage_threshold=0, evalue_threshold=float('inf'), verbose:bool=True, quiet:bool=False, export_csv:bool=True, return_df:bool=False):
        
        start = time.time()
        
        out = ["gene_id_A\tgene_id_B\tscore\tsummary_score\tannotation_A\tannotation_B\tspecies_A\tspecies_B"]
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
                        out.append(f"{gene.id}\t{equivalence.id}\t{equivalence.type}\t{equivalence.score}\t{self.name}\t{equivalence.target_annotation}\t{self.species}\t{equivalence.species}\t{equivalence.reliability}")
                    else:
                        out.append(f"{gene.id}\t{equivalence.id}\t{equivalence.type}\t{equivalence.score}\t{self.name}\t{equivalence.target_annotation}\t{self.species}\t{equivalence.species}")
        else:
            for gene in self.genes.values():
                for equivalence in gene.equivalences:
                    if verbose:
                        out.append(f"{gene.id}\t{equivalence.id}\t{equivalence.type}\t{equivalence.score}\t{self.name}\t{equivalence.target_annotation}\t{self.species}\t{equivalence.species}\t{equivalence.reliability}")
                    else:
                        out.append(f"{gene.id}\t{equivalence.id}\t{equivalence.type}\t{equivalence.score}\t{self.name}\t{equivalence.target_annotation}\t{self.species}\t{equivalence.species}")

        if return_df:
            columns = out[0].split('\t')
            data = [line.split('\t') for line in out[1:]]
            df = pd.DataFrame(data, columns=columns, dtype=str)
        
        if export_csv:
            out = "\n".join(out)
            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()

        end = time.time()
        lapse = end - start
        if not quiet:
            print(f"Exporting all filtered={filtered} equivalences for {self.name} took {round(lapse/60, 2)} minutes\n")

        if return_df:
            return df
        


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

            for gene_query, gene_target, score in zip(df[key_col], df[target_col], df[score_col]):
                self.genes[gene_query].equivalences.append(
                    Equivalence(gene_target, equivalence_type, target_annotation, species, score)
                )


    def add_orthofinder_equivalences(self, file, target_annotation, species):

        if os.path.isfile(file):
            df = pd.read_csv(file, sep="\t", dtype=str)
            df.columns = ["orthogroup", "query", "target"]
            for queries_str, targets_str, score in zip(df["query"], df["target"], df["orthogroup"]):
                gene_queries = queries_str.split(", ")
                gene_targets = targets_str.split(", ")
                for gene_query in gene_queries:
                    for gene_target in gene_targets:
                        self.genes[gene_query].equivalences.append(
                            Equivalence(gene_target, "orthofinder", target_annotation, species, score)
                        )

    def add_blast_equivalences(self, blast_folder, query_annotation, target_annotation, species, skip_rbhs:bool=False, skip_unidirectional_blasts:bool=False, proteins:bool=True, quiet:bool=False):

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

            rbbh_file = f"{blast_folder}/rbbh_{query_annotation}__to__{target_annotation}.txt"
            rbh_file = f"{blast_folder}/rbh_{query_annotation}__to__{target_annotation}.txt"
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

                rbbh_file = f"{blast_folder}/rbbh_{target_annotation}__to__{query_annotation}.txt"
                rbh_file = f"{blast_folder}/rbh_{target_annotation}__to__{query_annotation}.txt"

            fwd_file = f"{blast_folder}/single_{query_annotation}__to__{target_annotation}.txt"
            rev_file = f"{blast_folder}/single_{target_annotation}__to__{query_annotation}.txt"

            rbbh_hits = set()

            temp_start = time.time()

            if os.path.isfile(rbbh_file):
                df = pd.read_csv(rbbh_file, sep="\t", dtype=str)
                for gene_query, gene_target, q_ev, t_ev, q_nb, t_nb, q_cov, t_cov, q_id, t_id in zip(
                    df[query_col], df[target_col], df[query_col_evalue], df[target_col_evalue],
                    df[query_col_normbit], df[target_col_normbit], df[query_col_coverage], df[target_col_coverage],
                    df[query_col_identity], df[target_col_identity]
                ):
                    score = f"evalue={q_ev}/{t_ev}"
                    score += f", norm_bitscore={round(float(q_nb), 1)}/{round(float(t_nb), 1)}"
                    score += f", coverage={round(float(q_cov), 1)}/{round(float(t_cov), 1)}"
                    score += f", identity={round(float(q_id), 1)}/{round(float(t_id), 1)}"
                    pair = sorted([gene_query, gene_target])
                    pair = "-".join(pair)
                    rbbh_hits.add(pair)
                    self.genes[gene_query].equivalences.append(
                        Equivalence(gene_target, "rbbh", target_annotation, species, score, q_ev)
                    )

            end = time.time()
            lapse = end - temp_start
            if not quiet:
                print(f"Adding RBBH equivalences to {self.name} '{query_annotation} vs {target_annotation}' took {round(lapse/60, 2)} minutes")

            if not skip_rbhs:

                temp_start = time.time()

                rbh_hits = set()

                if os.path.isfile(rbh_file):
                    df = pd.read_csv(rbh_file, sep="\t", dtype=str)
                    for gene_query, gene_target, q_ev, t_ev, q_nb, t_nb, q_cov, t_cov, q_id, t_id in zip(
                        df[query_col], df[target_col], df[query_col_evalue], df[target_col_evalue],
                        df[query_col_normbit], df[target_col_normbit], df[query_col_coverage], df[target_col_coverage],
                        df[query_col_identity], df[target_col_identity]
                    ):
                        score = f"evalue={q_ev}/{t_ev}"
                        score += f", norm_bitscore={round(float(q_nb), 1)}/{round(float(t_nb), 1)}"
                        score += f", coverage={round(float(q_cov), 1)}/{round(float(t_cov), 1)}"
                        score += f", identity={round(float(q_id), 1)}/{round(float(t_id), 1)}"
                        pair = sorted([gene_query, gene_target])
                        pair = "-".join(pair)
                        if pair not in rbbh_hits:
                            rbh_hits.add(pair)
                            self.genes[gene_query].equivalences.append(
                                Equivalence(gene_target, "rbh", target_annotation, species, score, q_ev)
                            )

                end = time.time()
                lapse = end - temp_start
                if not quiet:
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
                        for q, s, pident, qcov, qlen, _, _, bitscore, evalue in zip(
                            df["0"], df["1"], df["2"], df["3"], df["4"], df["5"], df["6"], df["7"], df["8"]
                        ):
                            norm_bitscore = round((float(bitscore) / float(qlen)), 1)
                            score = f"evalue={evalue}"
                            score += f", norm_bitscore={norm_bitscore}"
                            score += f", coverage={round(float(qcov), 1)}"
                            score += f", identity={round(float(pident), 1)}"
                            pair = sorted([q, s])
                            pair = "-".join(pair)
                            if pair not in rbh_hits and pair not in rbbh_hits:
                                self.genes[q].equivalences.append(Equivalence(s, equivalence_type, target_annotation, species, score, evalue))

                    if proteins:
                        equivalence_type = "rev_blastp"
                    else:
                        equivalence_type = "rev_blast"

                    if os.path.isfile(rev_file):
                        df = pd.read_csv(rev_file, sep="\t", dtype=str, header=None)
                        df.columns = ["0", "1", "2", "3", "4", "5", "6", "7", "8"]
                        for q, s, pident, qcov, qlen, _, _, bitscore, evalue in zip(
                            df["0"], df["1"], df["2"], df["3"], df["4"], df["5"], df["6"], df["7"], df["8"]
                        ):
                            gene_query = s
                            gene_target = q
                            norm_bitscore = round((float(bitscore) / float(qlen)), 1)
                            score = f"evalue={evalue}"
                            score += f", norm_bitscore={norm_bitscore}"
                            score += f", coverage={round(float(qcov), 1)}"
                            score += f", identity={round(float(pident), 1)}"
                            pair = sorted([gene_query, gene_target])
                            pair = "-".join(pair)
                            if pair not in rbh_hits and pair not in rbbh_hits:
                                self.genes[gene_query].equivalences.append(Equivalence(gene_target, equivalence_type, target_annotation, species, score, evalue))

                    end = time.time()
                    lapse = end - temp_start
                    if not quiet:
                        print(f"Adding fwd/rev blast equivalences to {self.name} '{query_annotation} vs {target_annotation}' took {round(lapse/60, 2)} minutes")

            end = time.time()
            lapse = end - start
            
            if not quiet:
                if skip_unidirectional_blasts and not skip_rbhs:
                    print(f"Adding RBBH, and RBH equivalences to {self.name} '{query_annotation} vs {target_annotation}' took {round(lapse/60, 2)} minutes")
                else:
                    print(f"Adding RBBH, RBH, and fwd/rev blast equivalences to {self.name} '{query_annotation} vs {target_annotation}' took {round(lapse/60, 2)} minutes")

    def add_reciprocal_overlap_equivalences(self, folder, query_tag, target_tag, species, liftoff:bool=True, quiet:bool=False, synteny_present:bool=False):

        program = "liftoff" if liftoff else "lifton"
        fwd_file = f"{folder}/{program}__{query_tag}__to__{target_tag}_overlaps.tsv"
        rev_file = f"{folder}/{program}__{target_tag}__to__{query_tag}_overlaps.tsv"

        if os.path.isfile(fwd_file) and os.path.isfile(rev_file):
            start = time.time()

            fwd_df = pd.read_csv(fwd_file, sep="\t", encoding="utf-8", dtype=str, na_filter=False)
            if not fwd_df[fwd_df.duplicated(subset=["gene_id_A", "gene_id_B"], keep=False)].empty:
                raise ValueError(f"Duplicate query {query_tag} and target {target_tag} id pairs for fwd {program}.")

            rev_df = pd.read_csv(rev_file, sep="\t", encoding="utf-8", dtype=str, na_filter=False)
            if not rev_df[rev_df.duplicated(subset=["gene_id_A", "gene_id_B"], keep=False)].empty:
                raise ValueError(f"Duplicate query {query_tag} and target {target_tag} id pairs for rev {program}.")

            fwd_df = fwd_df[(fwd_df["gene_id_A"] != "NA") & (fwd_df["gene_id_B"] != "NA")]
            rev_df = rev_df[(rev_df["gene_id_A"] != "NA") & (rev_df["gene_id_B"] != "NA")]

            fwd_df = fwd_df[fwd_df["gene_id_A_origin"] == f"{program}_{self.name}"]
            rev_df = rev_df[rev_df["gene_id_B_origin"] == self.name]

            equivalence_suffix = "liftoff_aegis" if liftoff else "lifton_aegis"

            fwd_hits = {}
            fwd_synteny_col = fwd_df["gene_id_A_synteny_conserved"] if (synteny_present and "gene_id_A_synteny_conserved" in fwd_df.columns) else [None] * len(fwd_df)
            
            for g_q, g_t, q_cp, t_cp, synt_val, score_val in zip(
                fwd_df["gene_id_A"], fwd_df["gene_id_B"], fwd_df["gene_id_A_copy"], fwd_df["gene_id_B_copy"], fwd_synteny_col, fwd_df["overlap_score"]
            ):
                copies = (q_cp == "True" or t_cp == "True")
                gene_query = "_".join(g_q.split("_")[:-1]) if q_cp == "True" else g_q
                gene_target = "_".join(g_t.split("_")[:-1]) if t_cp == "True" else g_t

                pair = f"{gene_query}---{gene_target}"
                synteny = (synt_val == "True") if synt_val is not None else False

                modifiers = []
                if synteny: modifiers.append("synteny")
                if copies: modifiers.append("copies")
                modifier_str = f" ({' and '.join(modifiers)})" if modifiers else ""
                score = f"{score_val}{modifier_str}"

                if pair not in fwd_hits:
                    fwd_hits[pair] = score
                else:
                    old_score = int(fwd_hits[pair].split(" (")[0])
                    new_score = max(old_score, int(score_val))
                    
                    match_num = 2
                    if "multiple " in fwd_hits[pair]:
                        match_num = int(fwd_hits[pair].split("multiple ")[-1].split(")")[0]) + 1

                    has_copies = "copies" in fwd_hits[pair] or copies
                    has_synteny = "synteny" in fwd_hits[pair] or synteny

                    mods = []
                    if has_synteny: mods.append("synteny")
                    if has_copies: mods.append("copies")
                    mod_str = f" ({' and '.join(mods)})" if mods else ""
                    fwd_hits[pair] = f"{new_score}{mod_str} (multiple {match_num})"

            rev_hits = {}
            rev_synteny_col = rev_df["gene_id_B_synteny_conserved"] if (synteny_present and "gene_id_B_synteny_conserved" in rev_df.columns) else [None] * len(rev_df)
            
            for g_q, g_t, q_cp, t_cp, synt_val, score_val in zip(
                rev_df["gene_id_B"], rev_df["gene_id_A"], rev_df["gene_id_B_copy"], rev_df["gene_id_A_copy"], rev_synteny_col, rev_df["overlap_score"]
            ):
                copies = (q_cp == "True" or t_cp == "True")
                gene_query = "_".join(g_q.split("_")[:-1]) if q_cp == "True" else g_q
                gene_target = "_".join(g_t.split("_")[:-1]) if t_cp == "True" else g_t

                pair = f"{gene_query}---{gene_target}"
                synteny = (synt_val == "True") if synt_val is not None else False

                modifiers = []
                if synteny: modifiers.append("synteny")
                if copies: modifiers.append("copies")
                modifier_str = f" ({' and '.join(modifiers)})" if modifiers else ""
                score = f"{score_val}{modifier_str}"

                if pair not in rev_hits:
                    rev_hits[pair] = score
                else:
                    old_score = int(rev_hits[pair].split(" (")[0])
                    new_score = max(old_score, int(score_val))
                    
                    match_num = 2
                    if "multiple " in rev_hits[pair]:
                        match_num = int(rev_hits[pair].split("multiple ")[-1].split(")")[0]) + 1

                    has_copies = "copies" in rev_hits[pair] or copies
                    has_synteny = "synteny" in rev_hits[pair] or synteny

                    mods = []
                    if has_synteny: mods.append("synteny")
                    if has_copies: mods.append("copies")
                    mod_str = f" ({' and '.join(mods)})" if mods else ""
                    rev_hits[pair] = f"{new_score}{mod_str} (multiple {match_num})"

            for fwd_pair, f_score in fwd_hits.items():
                gene_query, gene_target = fwd_pair.split("---")
                if fwd_pair in rev_hits:
                    equivalence_type = f"rec_{equivalence_suffix}"
                    score = f"{f_score}/{rev_hits[fwd_pair]}"
                else:
                    equivalence_type = f"fwd_{equivalence_suffix}"
                    score = f_score

                self.genes[gene_query].equivalences.append(Equivalence(gene_target, equivalence_type, target_tag, species, score))

            for rev_pair, r_score in rev_hits.items():
                if rev_pair not in fwd_hits:
                    equivalence_type = f"rev_{equivalence_suffix}"
                    gene_query, gene_target = rev_pair.split("---")
                    self.genes[gene_query].equivalences.append(Equivalence(gene_target, equivalence_type, target_tag, species, r_score))

            end = time.time()
            lapse = end - start
            if not quiet:
                print(f"Adding {program} overlap equivalences for tags = [{query_tag}, {target_tag}] to {self.name} took {round(lapse/60, 2)} minutes")
        else:
            print(f"Warning: {fwd_file} or {rev_file} is missing.")

    def add_liftless_overlap_equivalences(self, folder, query_tag, target_tag, species, quiet:bool=False):

        fwd_file = f"{folder}/{query_tag}__to__{target_tag}_overlaps.tsv"

        if os.path.isfile(fwd_file):
            start = time.time()

            fwd_df = pd.read_csv(fwd_file, sep="\t", encoding="utf-8", dtype=str, na_filter=False)

            duplicates = fwd_df[fwd_df.duplicated(subset=["gene_id_A", "gene_id_B"], keep=False)]

            if not duplicates.empty:
                raise ValueError(f"Duplicate query {query_tag} and target {target_tag} id pairs for aegis_overlap:\n{duplicates}")

            fwd_df = fwd_df[fwd_df["gene_id_A"] != "NA"]
            fwd_df = fwd_df[fwd_df["gene_id_B"] != "NA"]

            fwd_df = fwd_df[fwd_df["gene_id_A_origin"] == self.name]

            equivalence_suffix = "aegis_overlap"

            hits = {}
        
            for gene_query, gene_target, score in zip(fwd_df["gene_id_A"], fwd_df["gene_id_B"], fwd_df["overlap_score"]):
                pair = f"{gene_query}---{gene_target}"
                if pair not in hits:
                    hits[pair] = score
                else:
                    raise ValueError(f"There should not be repeated {query_tag} and {target_tag} id pairs within '{query_tag}__to__{target_tag}_overlaps.tsv' file. The pair: {pair} was already found.")

            for pair, score in hits.items():
                gene_query, gene_target = pair.split("---")
                self.genes[gene_query].equivalences.append(Equivalence(gene_target, equivalence_suffix, target_tag, species, score))

            end = time.time()
            lapse = end - start
            if not quiet:
                print(f"Adding aegis overlap equivalences for tags = [{query_tag}, {target_tag}] to {self.name} took {round(lapse/60, 2)} minutes")

        else:
            print(f"Warning: {fwd_file} is missing.")


    def filter_equivalences(self, simple_rbh_blasts:bool=True, unidirectional_blasts:bool=True, replace:bool=True, identity_threshold=0, coverage_threshold=0, evalue_threshold=float("inf")):

        for gene in self.genes.values():
            gene.filter_equivalences(simple_rbh_blasts=simple_rbh_blasts, unidirectional_blasts=unidirectional_blasts, replace=replace, identity_threshold=identity_threshold, coverage_threshold=coverage_threshold, evalue_threshold=evalue_threshold)
