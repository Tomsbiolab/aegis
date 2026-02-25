
from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..annotation import Annotation

import time
import sys
import networkx as nx
import pandas as pd

from pathlib import Path
from tqdm import tqdm

from ..hits import OverlapHit

class AnnotationOverlaps:
    """
    Component for handling overlaps methods for the Annotation class.
    Accessed via 'annotation_object.overlaps'.
    """
    def __init__(self, annotation:Annotation):
        self._annot = annotation
        self.self = []
        self.other = []
    
    def detect(self, other:Annotation|None=None, sort_processes:int=1, clear=True, quiet:bool=True):
        """
        Detecting gene overlaps within the same annotation object or between
        annotation objects, provided they refer to the same genome.
        """
        start_time = time.time()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False

        if not self._annot.sorted:
            self._annot.sort_genes(processes=sort_processes)
        
        if other != None:
            if not other.sorted:
                other.sort_genes(processes=sort_processes)

        if clear:
            self.clear()
            if other != None:
                other.overlaps.clear()

        if other != None:

            if self._annot.genome == other.genome:

                if self._annot.genome == None:
                    if not quiet:
                        print(f"Note: Make sure that both annotations that are being compared are associated to the same genome version. Otherwise the resulting coordinate overlaps will not be correct.")
                
                if other.name in self.other or self._annot.name in other.overlaps.other:
                    print(f"Overlaps between {self._annot.id} and {other.id} "
                           "annotations have already been detected, please "
                           "run 'self.overlaps.clear()' and/or other.overlaps.clear() if you want to "
                           "recalculate them")
                else:
                    start_time = time.time()
                    progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mWorking out overlaps between {self._annot.id} and {other.id} annotations:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

                    if other.name not in self.other:
                        self.other.append(other.name)
                    if self._annot.name not in other.overlaps.other:
                        other.overlaps.other.append(self._annot.name)
                    for chr, genes in self._annot.chrs.items():
                        for g1 in genes.values():
                            progress_bar.update(1)
                            found_overlap = False
                            if chr in other.chrs:
                                for g2 in other.chrs[chr].values():
                                    overlapping, overlap_bp = g1.overlap(g2)
                                    if overlapping:
                                        found_overlap = True
                                    elif found_overlap and g1.end < g2.start:
                                        break
                                    else:
                                        continue
                                    if g1.strand == g2.strand:
                                        gene_orientation = True
                                    else:
                                        gene_orientation = False

                                    gene_query_percent = (overlap_bp / g1.size) * 100                  
                                    gene_target_percent = (overlap_bp / g2.size) * 100

                                    target_exons = False
                                    query_exons = False
                                    best_exon_overlap = 0
                                    exon_orientation = False
                                    overlapping = False
                                    for t1 in g1.transcripts.values():
                                        if t1.exons != []:
                                            query_exons = True
                                        for t2 in g2.transcripts.values():
                                            if t2.exons != []:
                                                target_exons = True
                                            if t1.strand == t2.strand:
                                                exon_orientation = True
                                            overlap_exon_temp = 0
                                            for e1 in t1.exons:
                                                for e2 in t2.exons:
                                                    overlap_temp, overlap_bp = e1.overlap(e2)
                                                    if overlap_temp:
                                                        overlap_exon_temp += overlap_bp
                                                        overlapping = True
                                            if overlap_exon_temp > best_exon_overlap:
                                                best_exon_overlap = overlap_exon_temp
                                                exon_query_size = t1.size
                                                exon_target_size = t2.size
                                    
                                    if target_exons and query_exons:
                                        exons_in_both = True
                                        if gene_orientation != exon_orientation:
                                            print(f"Warning: {self._annot.id} query and {other.id} target have discrepancies in the orientation of gene and exons. Genes: {g1.id} and {g2.id}.")
                                        if overlapping:
                                            exon_query_percent = (best_exon_overlap / exon_query_size) * 100
                                            exon_target_percent = (best_exon_overlap / exon_target_size) * 100
                                        else:
                                            exon_query_percent = 0
                                            exon_target_percent = 0
                                    else:
                                        exons_in_both = False
                                        exon_orientation = None
                                        exon_query_percent = None
                                        exon_target_percent = None


                                    target_CDS = False
                                    query_CDS = False
                                    best_CDS_overlap = 0
                                    CDS_orientation = False
                                    overlapping = False
                                    for t1 in g1.transcripts.values():
                                        if t1.CDSs != {}:
                                            query_CDS = True
                                        for t2 in g2.transcripts.values():
                                            if t2.CDSs != {}:
                                                target_CDS = True
                                            if t1.strand == t2.strand:
                                                CDS_orientation = True
                                            for CDS1 in t1.CDSs.values():
                                                for CDS2 in t2.CDSs.values():
                                                    overlap_CDS_temp = 0
                                                    for c1 in CDS1.CDS_segments:
                                                        for c2 in CDS2.CDS_segments:
                                                            overlap_temp, overlap_bp = c1.overlap(c2)
                                                            if not overlap_temp:
                                                                continue
                                                            overlap_CDS_temp += overlap_bp
                                                            overlapping = True
                                                    if overlap_CDS_temp > best_CDS_overlap:
                                                        best_CDS_overlap = overlap_CDS_temp
                                                        CDS_query_size = CDS1.size
                                                        CDS_target_size = CDS2.size
                                                
                                    if target_CDS and query_CDS:
                                        CDSs_in_both = True
                                        if gene_orientation != CDS_orientation:
                                            print(f"Warning: {self._annot.id} query and {other.id} target have discrepancies in the orientation of gene and CDS. Genes: {g1.id} and {g2.id}.")
                                        if overlapping:
                                            CDS_query_percent = (best_CDS_overlap / CDS_query_size) * 100
                                            CDS_target_percent = (best_CDS_overlap / CDS_target_size) * 100
                                        else:
                                            CDS_query_percent = 0
                                            CDS_target_percent = 0
                                    else:
                                        CDSs_in_both = False
                                        CDS_orientation = None
                                        CDS_query_percent = None
                                        CDS_target_percent = None


                                    protein_query_percent = None
                                    protein_target_percent = None
                                    if CDS_query_percent != None and CDS_query_percent != 0:
                                        if CDS_orientation:
                                            target_protein = False
                                            query_protein = False
                                            best_protein_overlap = 0
                                            overlapping = False
                                            for t1 in g1.transcripts.values():
                                                if t1.CDSs != {}:
                                                    query_protein = True
                                                for t2 in g2.transcripts.values():
                                                    if t2.CDSs != {}:
                                                        target_protein = True
                                                    for CDS1 in t1.CDSs.values():
                                                        for CDS2 in t2.CDSs.values():
                                                            overlap_protein_temp = 0
                                                            for c1 in CDS1.CDS_segments:
                                                                for c2 in CDS2.CDS_segments:
                                                                    overlap_temp, overlap_bp = c1.overlap(c2)
                                                                    if not overlap_temp:
                                                                        continue
                                                                    if c1.frame != c2.frame:
                                                                        continue
                                                                    
                                                                    overlap_protein_temp += overlap_bp
                                                                    overlapping = True
                                                            if overlap_protein_temp > best_protein_overlap:
                                                                best_protein_overlap = overlap_protein_temp
                                                                protein_query_size = CDS1.size
                                                                protein_target_size = CDS2.size
                                                        
                                            if target_protein and query_protein:
                                                if overlapping:
                                                    protein_query_percent = (best_protein_overlap / protein_query_size) * 100
                                                    protein_target_percent = (best_protein_overlap / protein_target_size) * 100
                                                else:
                                                    protein_query_percent = 0
                                                    protein_target_percent = 0

                                    g1.overlaps["other"].append(OverlapHit(g2.id, 
                                                                    other.name,
                                                                    gene_orientation,
                                                                    gene_query_percent,
                                                                    gene_target_percent,
                                                                    exons_in_both,
                                                                    exon_query_percent,
                                                                    exon_target_percent,
                                                                    CDSs_in_both,
                                                                    CDS_query_percent,
                                                                    CDS_target_percent,
                                                                    protein_query_percent,
                                                                    protein_target_percent,
                                                                    g2.conserved_synteny,
                                                                    g2.extra_copy))

                                    g2.overlaps["other"].append(OverlapHit(g1.id,
                                                                    self._annot.name,
                                                                    gene_orientation,
                                                                    gene_target_percent,
                                                                    gene_query_percent,
                                                                    exons_in_both,
                                                                    exon_target_percent,
                                                                    exon_query_percent,
                                                                    CDSs_in_both,
                                                                    CDS_target_percent,
                                                                    CDS_query_percent,
                                                                    protein_query_percent,
                                                                    protein_target_percent,
                                                                    g1.conserved_synteny,
                                                                    g1.extra_copy))
                    self._annot.add_aliases()
                    other.add_aliases()
                    now = time.time()
                    lapse = now - start_time
                    progress_bar.close()
                    if not quiet:
                        print(f"\nDetecting overlaps between {other.id} and {self._annot.id} annotations took {round(lapse/60, 1)} minutes")
            else:
                print(f"Did not generate overlaps between {other.id} and {self._annot.id} annotations as they are associated to different genomes")

        else:
            if self.self != []:
                print("There are already detected 'self' gene overlaps, please run 'self.clear_overlaps()' if you want to recalculate them")
            else:
                progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                            bar_format=(
                f'\033[1;91mWorking out overlaps within {self._annot.id} annotation:\033[0m '
                '{percentage:3.0f}%|'
                f'\033[1;91m{{bar}}\033[0m| '
                '{n}/{total} [{elapsed}<{remaining}]'))

                # making sure self overlaps are not added twice
                start_time = time.time()
                self.self = set(self.self)
                for chr, genes in self._annot.chrs.items():
                    gl = list(genes.keys())[1:]
                    for g1 in genes.values():
                        progress_bar.update(1)
                        found_overlap = False
                        for gl_id in gl:
                            g2 = genes[gl_id]
                            if g1.id == g2.id:
                                continue
                            overlapping, overlap_bp = g1.overlap(g2)

                            if overlapping:
                                self.self.add(g1.id)
                                self.self.add(g2.id)
                                found_overlap = True

                            elif found_overlap and g1.end < g2.start:
                                break
                            else:
                                continue

                            if g1.strand == g2.strand:
                                gene_orientation = True
                            else:
                                gene_orientation = False

                            gene_query_percent = (overlap_bp / g1.size) * 100                  
                            gene_target_percent = (overlap_bp / g2.size) * 100

                            target_exons = False
                            query_exons = False
                            best_exon_overlap = 0
                            exon_orientation = False
                            overlapping = False
                            for t1 in g1.transcripts.values():
                                if t1.exons != []:
                                    query_exons = True
                                for t2 in g2.transcripts.values():
                                    if t2.exons != []:
                                        target_exons = True
                                    if t1.strand == t2.strand:
                                        exon_orientation = True
                                    overlap_exon_temp = 0
                                    for e1 in t1.exons:
                                        for e2 in t2.exons:
                                            overlap_temp, overlap_bp = e1.overlap(e2)
                                            if overlap_temp:
                                                overlap_exon_temp += overlap_bp
                                                overlapping = True
                                    if overlap_exon_temp > best_exon_overlap:
                                        best_exon_overlap = overlap_exon_temp
                                        exon_query_size = t1.size
                                        exon_target_size = t2.size
                            
                            if target_exons and query_exons:
                                exons_in_both = True
                                if gene_orientation != exon_orientation:
                                    print(f"Warning: {self._annot.id} query and target have discrepancies in the orientation of gene and exons. Genes: {g1.id} and {g2.id}")
                                if overlapping:
                                    exon_query_percent = (best_exon_overlap / exon_query_size) * 100
                                    exon_target_percent = (best_exon_overlap / exon_target_size) * 100
                                else:
                                    exon_query_percent = 0
                                    exon_target_percent = 0
                            else:
                                exons_in_both = False
                                exon_orientation = None
                                exon_query_percent = None
                                exon_target_percent = None

                            target_CDS = False
                            query_CDS = False
                            best_CDS_overlap = 0
                            CDS_orientation = False
                            overlapping = False
                            for t1 in g1.transcripts.values():
                                if t1.CDSs != {}:
                                    query_CDS = True
                                for t2 in g2.transcripts.values():
                                    if t2.CDSs != {}:
                                        target_CDS = True
                                    if t1.strand == t2.strand:
                                        CDS_orientation = True
                                    for CDS1 in t1.CDSs.values():
                                        for CDS2 in t2.CDSs.values():
                                            overlap_CDS_temp = 0
                                            for c1 in CDS1.CDS_segments:
                                                for c2 in CDS2.CDS_segments:
                                                    overlap_temp, overlap_bp = c1.overlap(c2)
                                                    if not overlap_temp:
                                                        continue
                                                    overlap_CDS_temp += overlap_bp
                                                    overlapping = True
                                            if overlap_CDS_temp > best_CDS_overlap:
                                                best_CDS_overlap = overlap_CDS_temp
                                                CDS_query_size = CDS1.size
                                                CDS_target_size = CDS2.size
                                        
                            if target_CDS and query_CDS:
                                CDSs_in_both = True
                                if gene_orientation != CDS_orientation:
                                    print(f"Error: {self._annot.id} query and target have discrepancies in the orientation of gene and CDS. Genes: {g1.id} and {g2.id}. DO NOT CONTINUE! -> fix the problem!")
                                if overlapping:
                                    CDS_query_percent = (best_CDS_overlap / CDS_query_size) * 100
                                    CDS_target_percent = (best_CDS_overlap / CDS_target_size) * 100
                                else:
                                    CDS_query_percent = 0
                                    CDS_target_percent = 0
                            else:
                                CDSs_in_both = False
                                CDS_orientation = None
                                CDS_query_percent = None
                                CDS_target_percent = None

                            if CDSs_in_both:
                                protein_query_percent = 0
                                protein_target_percent = 0
                                if CDS_query_percent != None and CDS_query_percent != 0:
                                    if CDS_orientation:
                                        target_protein = False
                                        query_protein = False
                                        best_protein_overlap = 0
                                        overlapping = False
                                        for t1 in g1.transcripts.values():
                                            if t1.CDSs != {}:
                                                query_protein = True
                                            for t2 in g2.transcripts.values():
                                                if t2.CDSs != {}:
                                                    target_protein = True
                                                for CDS1 in t1.CDSs.values():
                                                    for CDS2 in t2.CDSs.values():
                                                        overlap_protein_temp = 0
                                                        for c1 in CDS1.CDS_segments:
                                                            for c2 in CDS2.CDS_segments:
                                                                overlap_temp, overlap_bp = c1.overlap(c2)
                                                                if not overlap_temp:
                                                                    continue
                                                                if c1.frame != c2.frame:
                                                                    continue
                                                                
                                                                overlap_protein_temp += overlap_bp
                                                                overlapping = True
                                                        if overlap_protein_temp > best_protein_overlap:
                                                            best_protein_overlap = overlap_protein_temp
                                                            protein_query_size = CDS1.size
                                                            protein_target_size = CDS2.size
                                                    
                                        if target_protein and query_protein:
                                            if overlapping:
                                                protein_query_percent = (best_protein_overlap / protein_query_size) * 100
                                                protein_target_percent = (best_protein_overlap / protein_target_size) * 100
                            else:

                                protein_query_percent = None
                                protein_target_percent = None

                            g1.overlaps["self"].append(OverlapHit(g2.id, self._annot.name,
                                                                gene_orientation,
                                                                gene_query_percent,
                                                                gene_target_percent,
                                                                exons_in_both,
                                                                exon_query_percent,
                                                                exon_target_percent,
                                                                CDSs_in_both,
                                                                CDS_query_percent,
                                                                CDS_target_percent,
                                                                protein_query_percent,
                                                                protein_target_percent,
                                                                g2.conserved_synteny,
                                                                g2.extra_copy))

                            g2.overlaps["self"].append(OverlapHit(g1.id, self._annot.name,
                                                                gene_orientation,
                                                                gene_target_percent,
                                                                gene_query_percent,
                                                                exons_in_both,
                                                                exon_target_percent,
                                                                exon_query_percent,
                                                                CDSs_in_both,
                                                                CDS_target_percent,
                                                                CDS_query_percent,
                                                                protein_query_percent,
                                                                protein_target_percent,
                                                                g1.conserved_synteny,
                                                                g1.extra_copy))
                        try:
                            gl.remove(g1.id)
                        except:
                            pass

                self.self = list(self.self)
                progress_bar.close()
                now = time.time()
                lapse = now - start_time
                if not quiet:
                    print(f"\nDetecting gene overlaps within the {self._annot.id} annotation took {round(lapse/60, 1)} minutes\n")
                    print(f"\nThere are {len(self.self)} genes overlapping with other genes in {self._annot.id} annotation\n")
                self.add_qualitative_info()

    def as_networks(self, self_mode:bool=True):
        self.networks = {}
        for chr, genes in self._annot.chrs.items():
            G = nx.Graph()
            for g in genes.values():
                if self_mode:
                    overlaps = g.overlaps["self"]
                else:
                    overlaps = g.overlaps["other"]
                for o in overlaps:
                    G.add_edge(g.id, o.id)
            self.networks[chr] = list(nx.connected_components(G))


    def clear(self, keep_self=False, keep_other=False):
        if not keep_self:
            self.self = []
            for genes in self._annot.chrs.values():
                for g in genes.values():
                    g.overlaps["self"] = []
        if not keep_other:
            self.other = []
            for genes in self._annot.chrs.values():
                for g in genes.values():
                    g.overlaps["other"] = []

    def add_qualitative_info(self, quiet:bool=True):
        """
        Number of unique full segment overlaps between genes including all transcript variants.
        """
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        total = len(self._annot.all_gene_ids.keys())

        progress_bar = tqdm(total=total, disable=disable,
                                bar_format=(
                    f'\033[1;95mAdding qualitative info to {self._annot.id} overlaps:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;95m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self._annot.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                for o in g.overlaps["self"]:
                    if o.exon_query_percent > 0:

                        g_exons, g_CDSs, g_UTRs = set(), set(), set()
                        o_exons, o_CDSs, o_UTRs = set(), set(), set()

                        for t in g.transcripts.values():
                            if t.main:
                                g_exons.update(set([(e.start, e.end) for e in t.exons]))
                                for c in t.CDSs.values():
                                    if c.main:
                                        g_CDSs.update(set([(cs.start, cs.end, cs.frame) for cs in c.CDS_segments]))
                                        g_UTRs.update(set([(u.start, u.end) for u in c.UTRs]))
                                
                        chrom = self._annot.all_gene_ids[o.id]
                        for t in self._annot.chrs[chrom][o.id].transcripts.values():
                            if t.main:
                                o_exons.update(set([(e.start, e.end) for e in t.exons]))
                                for c in t.CDSs.values():
                                    if c.main:
                                        o_CDSs.update(set([(cs.start, cs.end, cs.frame) for cs in c.CDS_segments]))
                                        o_UTRs.update(set([(u.start, u.end) for u in c.UTRs]))

                        for e1 in g_exons:
                            for e2 in o_exons:
                                if e1 == e2:
                                    o.full_exon_overlaps += 1

                        for c1 in g_CDSs:
                            for c2 in o_CDSs:
                                if c1[0] == c2[0] and c1[1] == c2[1]:
                                    o.full_CDS_overlaps += 1
                                    if c1[2] == c2[2]:
                                        o.full_protein_overlaps += 1

                        for u1 in g_UTRs:
                            for u2 in o_UTRs:
                                if u1 == u2:
                                    o.full_UTR_overlaps += 1

        progress_bar.close()

    def clear_with_selected_CDSs(self):
        for genes in self._annot.chrs.values():
            for g in genes.values():
                g.overlap_with_selected_CDS = False   

    def clear_with_selected_exons(self):
        for genes in self._annot.chrs.values():
            for g in genes.values():
                g.overlap_with_selected_exon = False

    def export(self, custom_path: str = "", overlap_threshold: int = 6, verbose: bool = True, synteny: bool = False, NAs: bool = True, export_csv: bool = False, export_self: bool = False, output_file: str = "", quiet: bool = False, copies_info: bool = False) -> pd.DataFrame:
        start_time = time.time()
        if export_self:
            export = "self"
            export_tag = "self_"
        else:
            export = "other"
            export_tag = ""

        tag = f"{export_tag}{self._annot.id}{self._annot.feature_suffix}_overlap_t{overlap_threshold}"
        

        if custom_path:
            export_folder = Path(custom_path)
        else:
            export_folder = Path(custom_path or self._annot.path) / "overlaps"
        export_folder.mkdir(parents=True, exist_ok=True)
        export_folder = str(export_folder) + "/"

        correct_order = ["gene_id_A", "gene_id_B", "gene_id_A_synteny_conserved", "gene_id_B_synteny_conserved", "same_strand", "min_gene_percent", "min_exon_percent", "min_CDS_percent", "gene_id_A_origin", "gene_id_B_origin", "overlap_score", "gene_id_A_copy", "gene_id_B_copy"]

        rows = []

        for genes in self._annot.chrs.values():
            for g in genes.values():
                for name, hits in g.overlaps.items():
                    if name == export:
                        for hit in hits:
                            if hit.score < overlap_threshold:
                                continue

                            row = {
                                "gene_id_A": g.id,
                                "gene_id_B": hit.id,
                                "gene_id_A_origin": self._annot.name,
                                "gene_id_B_origin": hit.origin,
                                "overlap_score": hit.score,
                            }

                            if synteny:
                                row["gene_id_A_synteny_conserved"] = g.conserved_synteny if self._annot.liftover else pd.NA
                                row["gene_id_B_synteny_conserved"] = hit.target_synteny_conserved

                            if verbose:
                                row["same_strand"] = hit.orientation
                                row["min_gene_percent"] = hit.min_gene_percent
                                row["min_exon_percent"] = hit.min_exon_percent
                                row["min_CDS_percent"] = hit.min_CDS_percent

                            if copies_info:
                                row["gene_id_A_copy"] = g.extra_copy
                                row["gene_id_B_copy"] = hit.extra_copy

                            rows.append(row)

        eq_df = pd.DataFrame(rows)


        if eq_df.empty:

            if export == "self":
                print(f"\nNo {self._annot.id} self overlaps were detected.")
            else:
                print(f"\nNo {self._annot.id} overlaps to the following annotation(s) '{self.other}' were detected.")
            cols = [
                "gene_id_A",
                "gene_id_B",
                "gene_id_A_origin",
                "gene_id_B_origin",
                "overlap_score",
            ]

            if synteny:
                cols += [
                    "gene_id_A_synteny_conserved",
                    "gene_id_B_synteny_conserved",
                ]

            if verbose:
                cols += [
                    "same_strand",
                    "min_gene_percent",
                    "min_exon_percent",
                    "min_CDS_percent",
                ]

            if copies_info:
                cols += [
                    "gene_id_A_copy",
                    "gene_id_B_copy",
                ]
            rows = []
            eq_df = pd.DataFrame(rows, columns=cols)



        if export == "self":
            eq_df["sorted_id_pair"] = eq_df.apply(lambda row: tuple(sorted([row["gene_id_A"], row["gene_id_B"]])), axis=1)
            eq_df = eq_df.drop_duplicates(subset="sorted_id_pair").drop(columns="sorted_id_pair")
            eq_df.drop(inplace=True, columns=["gene_id_A_origin", "gene_id_B_origin"])

        else:
            eq_df["sorted_id_pair"] = eq_df.apply(lambda row: tuple(sorted([f'{row["gene_id_A"]}_{row["gene_id_A_origin"]}', f'{row["gene_id_B"]}_{row["gene_id_B_origin"]}'])), axis=1)
            eq_df = eq_df.drop_duplicates(subset="sorted_id_pair").drop(columns="sorted_id_pair")

        if NAs:
            tag += "_gene_id_A_NAs"
            if export == "self":
                overlapping_genes = set(pd.concat([eq_df["gene_id_A"], eq_df["gene_id_B"]]).dropna())
            else:
                overlapping_genes = set(eq_df["gene_id_A"].dropna())

            na_rows = []

            if not copies_info:

                if export == "self":

                    for genes in self._annot.chrs.values():
                        for g in genes.values():
                            if g.id not in overlapping_genes:
                                na_rows.append({
                                    "gene_id_A": g.id,
                                    "overlap_score": 0,
                                    "gene_id_A_synteny_conserved": g.conserved_synteny if self._annot.liftover else pd.NA
                                })

                else:

                    for genes in self._annot.chrs.values():
                        for g in genes.values():
                            if g.id not in overlapping_genes:
                                na_rows.append({
                                    "gene_id_A": g.id,
                                    "gene_id_A_origin": self._annot.name,
                                    "overlap_score": 0,
                                    "gene_id_A_synteny_conserved": g.conserved_synteny if self._annot.liftover else pd.NA
                                })

                if synteny:
                    for g_id in self._annot.unmapped:
                        na_rows.append({
                            "gene_id_A": g_id,
                            "gene_id_A_origin": self._annot.name
                        })

            else:

                if export == "self":

                    for genes in self._annot.chrs.values():
                        for g in genes.values():
                            if g.id not in overlapping_genes:
                                na_rows.append({
                                    "gene_id_A": g.id,
                                    "gene_id_A_copy": g.extra_copy,
                                    "overlap_score": 0,
                                    "gene_id_A_synteny_conserved": g.conserved_synteny if self._annot.liftover else pd.NA
                                })

                else:

                    for genes in self._annot.chrs.values():
                        for g in genes.values():
                            if g.id not in overlapping_genes:
                                na_rows.append({
                                    "gene_id_A": g.id,
                                    "gene_id_A_origin": self._annot.name,
                                    "overlap_score": 0,
                                    "gene_id_A_copy": g.extra_copy,
                                    "gene_id_A_synteny_conserved": g.conserved_synteny if self._annot.liftover else pd.NA
                                })

                if synteny:
                    for g_id in self._annot.unmapped:
                        na_rows.append({
                            "gene_id_A": g_id,
                            "gene_id_A_origin": self._annot.name
                        })          

            # Combine with the original df
            if na_rows:
                eq_df = pd.concat([eq_df, pd.DataFrame(na_rows)], ignore_index=True)

        eq_df = eq_df[[col for col in correct_order if col in eq_df.columns]]

        if synteny:
            column_sort_order = ["gene_id_A_origin", "gene_id_B_origin", "overlap_score", "gene_id_A_synteny_conserved", "gene_id_B_synteny_conserved", "gene_id_A", "gene_id_B"]
            ascending = [True, True, False, False, False, True, True]
        elif export == "self":
            column_sort_order = ["overlap_score", "gene_id_A", "gene_id_B"]
            ascending = [False, True, True]
        else:
            column_sort_order = ["gene_id_A_origin", "gene_id_B_origin", "overlap_score", "gene_id_A", "gene_id_B"]
            ascending = [True, True, False, True, True]            

        eq_df.sort_values(by=column_sort_order, ascending=ascending, inplace=True)
        eq_df.reset_index(drop=True, inplace=True)
        
        if export_csv:
            if output_file:
                export_path = f"{export_folder}{output_file}"
            else:
                export_path = f"{export_folder}{tag}.csv"

            eq_df.to_csv(export_path, sep="\t", index=False, na_rep="NA")

            now = time.time()
            lapse = now - start_time
            if not quiet:
                if export == "self":
                    print(f"\nExporting {self._annot.id} self overlaps took {round(lapse/60, 1)} minutes")
                else:
                    print(f"\nExporting {self._annot.id} overlaps to the following annotation(s) '{self.other}' took {round(lapse/60, 1)} minutes")
        
        return eq_df