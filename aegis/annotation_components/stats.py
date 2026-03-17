from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..annotation import Annotation

import pandas as pd

from pathlib import Path
from statistics import mean

from ..utils.plots import barplot, pie_chart
from ..subfeatures import Intron

class AnnotationStats:
    """
    Component for handling statistical methods and metric calculations for the Annotation class.
    Accessed via 'annotation_object.stats'.
    """
    data: dict
    _annot: Annotation

    def __init__(self, annotation:Annotation):
        self._annot = annotation
        self.data = {}

    def calculate_transcript_masking(self):
        for genes in self._annot.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    t.calculate_masking()
                    for c in t.CDSs.values():
                        c.calculate_masking()
        self._annot.update()
    
    def calculate_gc_content(self):
        for genes in self._annot.chrs.values():
            for g in genes.values():
                g.calculate_gc_content()
                for t in g.transcripts.values():
                    t.calculate_gc_content()
                    for c in t.CDSs.values():
                        c.calculate_gc_content()
    def gene_count(self):
        gene_objects = 0
        unique_gene_ids_in_overlaps = set()
        for genes in self._annot.chrs.values():
            gene_objects += len(genes)
            for g in genes.values():
                if g.overlaps is not None:
                    for o in g.overlaps["self"]:
                        unique_gene_ids_in_overlaps.add(o.id)
        print(f"There are {gene_objects} gene objects and {len(self._annot.all_gene_ids)} genes in all gene ids and {len(unique_gene_ids_in_overlaps)} ids contained in self overlaps.")

    def update(self, custom_path:str="", export:bool=False, max_x:int|None=None, quiet:bool=True):

        self._annot.update_features(quiet=quiet)
        self._annot.generate_introns()
        if not quiet:
            print(f"\nUpdating stats for {self._annot.id}")

        if not self._annot.contains_protein_sequences:
            if self._annot.genome is not None:
                self._annot.generate_proteins(readthrough="both")

        self.calculate_gc_content()

        if export:
            export_folder = Path(custom_path or self._annot.path) / "out_stats"
            export_folder.mkdir(parents=True, exist_ok=True)
            export_folder = str(export_folder) + "/"

        to_tally = ["coding_genes", "noncoding_genes", "CDSs_without_stop", "CDSs_with_stop"]

        self.data = {"mean_transcripts" : [], "mean_exons" : [], "mean_exon_size" : [], "mean_gene_size" : [], "mean_intron_size" : [], "mean_CDS_size" : [], "mean_UTR_size" : [], "mean_transcript_size" : [], "mean_five_prime_UTR_size" : [], "mean_three_prime_UTR_size" : []}

        if self._annot.genome is not None:
            self.data["mean_protein_size"] = []

        for key in to_tally:
            self.data[key] = []

        for ft, value in self._annot.features.items():
            self.data[ft] = value

        self.data["five_prime_UTRs"] = 0
        self.data["three_prime_UTRs"] = 0

        for genes in self._annot.chrs.values():
            for g in genes.values():
                self.data["mean_transcripts"].append(len(g.transcripts))
                self.data["mean_gene_size"].append(g.size)
                for t in g.transcripts.values():
                    if t.main:
                        if t.introns:
                            for i in t.introns:
                                self.data["mean_intron_size"].append(i.size)

                        if t.coding:
                            for c in t.CDSs.values():
                                if c.main:
                                    self.data["mean_CDS_size"].append(c.size)
                                    if hasattr(c, "UTRs"):
                                        utr_5 = False
                                        utr_3 = False
                                        u_size = 0
                                        u5_size = 0
                                        u3_size = 0
                                        for u in c.UTRs:
                                            if u.prime == "3'":
                                                utr_3 = True
                                                u3_size += u.size
                                            else:
                                                utr_5 = True
                                                u5_size += u.size
                                            u_size += u.size
                                        if utr_5:
                                            self.data["five_prime_UTRs"] += 1
                                        if utr_3:
                                            self.data["three_prime_UTRs"] += 1
                                        self.data["mean_UTR_size"].append(u_size)
                                        self.data["mean_five_prime_UTR_size"].append(u5_size)
                                        self.data["mean_three_prime_UTR_size"].append(u3_size)
                                    if self._annot.genome is not None and c.protein is not None:
                                        self.data["mean_protein_size"].append(c.protein.size)
                            self.data["coding_genes"].append(g.id)
                        else:
                            self.data["noncoding_genes"].append(g.id)

                        for e in t.exons:
                            self.data["mean_exon_size"].append(e.size)

                        self.data["mean_exons"].append(len(t.exons))
                        
                        self.data["mean_transcript_size"].append(t.size)

        # anything with mean will be also plotted as distribution plots:
        if export:
            for key in self.data:
                if "mean" in key:
                    tag = key.split("mean_")[1]
                    if tag[-1] != "s":
                        tag += "s"
                    barplot(self.data[key], export_folder, f"{self._annot.id}{self._annot.feature_suffix}_{tag}", f"Distribution of {self._annot.id} {tag}", max_x)   

        if self._annot.genome is not None:
            self.data["CDSs_without_stop"] = []
            self.data["CDSs_with_stop"] = []
            self.data["intron_composition"] = set()
            intron_stats = {}
            for b in Intron.canonical_seqs:
                intron_stats[f"intron-exon boundary: {b}"] = 0
            intron_stats["other_intron_seqs"] = 0
            for genes in self._annot.chrs.values():
                for g in genes.values():                    
                    for t in g.transcripts.values():
                        if t.main:
                            if t.coding:
                                for c in t.CDSs.values():
                                    if c.main:
                                        if c.protein != None:
                                            if not c.protein.end_stop:
                                                self.data["CDSs_without_stop"].append(g.id)
                                            else:
                                                self.data["CDSs_with_stop"].append(g.id)
                            if t.introns:
                                for i in t.introns:
                                    self.data["intron_composition"].add(i.boundary)
                                    if i.canonical:
                                        intron_stats[f"intron-exon boundary: {i.boundary}"] += 1
                                    else:
                                        intron_stats["other_intron_seqs"] += 1
            self.data["intron_composition"] = list(self.data["intron_composition"])

            if export:
                labels = list(intron_stats.keys())
                mod_labels = []
                for l in labels:
                    if "other" in l:
                        mod_labels.append("other")
                    else:
                        mod_labels.append(l.split(": ")[1])

                values = list(intron_stats.values())
                pie_chart(mod_labels, values, export_folder, f"{self._annot.id}{self._annot.feature_suffix}_intron_composition", f"Intron composition of {self._annot.id} annotation")
            self.data.update(intron_stats)

        for key in self.data:
            if "mean" in key:
                if len(self.data[key]) > 0:
                    self.data[key] = mean(self.data[key])
                else:
                    self.data[key] = 0
            #tallying
            elif key in to_tally:
                self.data[key] = len(self.data[key])

        if export:
            if self._annot.genome is not None:
                out_file = f"{self._annot.id}{self._annot.feature_suffix}_full_stats.csv"
            else:
                out_file = f"{self._annot.id}{self._annot.feature_suffix}_basic_stats.csv"

            warning_df = pd.DataFrame({k: pd.Series(sorted(v)) for k, v in self._annot.warnings.items()}).fillna('')
            error_df = pd.DataFrame({k: pd.Series(sorted(v)) for k, v in self._annot.errors.items()}).fillna('')

            warning_df.to_csv(f"{export_folder}{self._annot.id}{self._annot.feature_suffix}_warnings.csv", sep="\t", index=False)
            error_df.to_csv(f"{export_folder}{self._annot.id}{self._annot.feature_suffix}_errors.csv", sep="\t", index=False)

            f_out = open(f"{export_folder}{out_file}", "w", encoding="utf-8")
            f_out.write("")
            f_out.close()

            f_out = open(f"{export_folder}{out_file}", "a")
            x = -1
            for key, value in self.data.items():
                x += 1
                value_temp = value
                if isinstance(value, list):
                    if len(value) > 200:
                        value_temp = value[:199]
                        value_temp.append("...")
                if x == 0:
                    f_out.write(f"{key}\t{value_temp}")
                else:
                    f_out.write(f"\n{key}\t{value_temp}")
            f_out.close()
        if not quiet:
            print(f"\nUpdated stats for {self._annot.id}")
