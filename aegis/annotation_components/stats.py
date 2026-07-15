from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..annotation import Annotation

import pandas as pd
import warnings

from pathlib import Path
from statistics import mean

from ..utils.plots import barplot, pie_chart
from ..subfeatures import Intron
from .base import AnnotationComponent

KEY_DESCRIPTIONS = {
    # Mean metrics
    "mean_transcripts":           "Mean transcripts per gene",
    "mean_exons":                 "Mean exons per transcript",
    "mean_exon_size":             "Mean exon size (bp)",
    "mean_gene_size":             "Mean gene size (bp)",
    "mean_intron_size":           "Mean intron size (bp)",
    "mean_CDS_size":              "Mean CDS size (bp)",
    "mean_UTR_size":              "Mean UTR size (bp)",
    "mean_transcript_size":       "Mean transcript size (bp)",
    "mean_five_prime_UTR_size":   "Mean 5' UTR size (bp)",
    "mean_three_prime_UTR_size":  "Mean 3' UTR size (bp)",
    "mean_protein_size":          "Mean protein size (aa)",
    # Tally counts
    "coding_genes":               "Number of coding genes",
    "noncoding_genes":            "Number of non-coding genes",
    "CDSs_without_stop":          "Number of CDS without stop codon",
    "CDSs_with_stop":             "Number of CDS with stop codon",
    # Feature counts (generic ones present in most annotations)
    "gene":                       "Number of genes",
    "mRNA":                       "Number of mRNAs",
    "exon":                       "Number of exons",
    "CDS":                        "Number of CDS",
    "lnc_RNA":                    "Number of lncRNAs",
    "antisense_RNA":              "Number of antisense RNAs",
    "antisense_lncRNA":           "Number of antisense lncRNAs",
    "tRNA":                       "Number of tRNAs",
    "ncRNA":                      "Number of ncRNAs",
    "snoRNA":                     "Number of snoRNAs",
    "snRNA":                      "Number of snRNAs",
    "miRNA":                      "Number of miRNAs",
    "miRNA_primary_transcript":   "Number of miRNA primary transcripts",
    "pseudogenic_tRNA":           "Number of pseudogenic tRNAs",
    # UTR presence
    "five_prime_UTRs":            "Number of transcripts with 5' UTR",
    "three_prime_UTRs":           "Number of transcripts with 3' UTR",
    # Intron composition
    "intron_composition":         "Intron boundary types observed",
    "other_intron_seqs":          "Number of non-canonical introns",
    # Other stats
    "single_exon_genes":                "Number of single-exon genes",
    "mrnas_with_utr_both_sides":        "Number of mRNAs with UTR on both sides",
    "mrnas_with_at_least_one_utr":      "Number of mRNAs with at least one UTR",
    "introns_in_5utr":                  "Number of introns in 5' UTR",
    "introns_in_3utr":                  "Number of introns in 3' UTR",
    # Total lengths
    "total_length_gene":                "Total length of genes (bp)",
    "total_length_mRNA":                "Total spliced mRNA length (bp)",
    "total_length_CDS":                 "Total CDS length (bp)",
    "total_length_exon":                "Total exon length (bp)",
    "total_length_intron":              "Total intron length (bp)",
    "total_length_5utr":                "Total 5' UTR length (bp)",
    "total_length_3utr":                "Total 3' UTR length (bp)",
    # Longest features
    "longest_gene":                     "Longest gene",
    "longest_mRNA":                     "Longest mRNA (spliced)",
    "longest_exon":                     "Longest exon (bp)",
    "longest_intron":                   "Longest intron (bp)",
    "longest_CDS":                      "Longest CDS (bp)",
    "longest_CDS_segment":              "Longest CDS segment (bp)",
    # Shortest features
    "shortest_gene":                    "Shortest gene",
    "shortest_mRNA":                    "Shortest mRNA (spliced)",
    "shortest_exon":                    "Shortest exon (bp)",
    "shortest_intron":                  "Shortest intron (bp)",
    "shortest_CDS":                     "Shortest CDS (bp)",
    "shortest_CDS_segment":             "Shortest CDS segment (bp)",
}


class AnnotationStats(AnnotationComponent):
    """
    Component for handling statistical methods and metric calculations for the Annotation class.
    Accessed via 'annotation_object.stats'.
    """
    data: dict

    def __init__(self, annot:Annotation):
        super().__init__(annot)
        self.data = {}

    def calculate_transcript_masking(self):
        if self._annot.genome is not None:
            for genes in self._annot.chrs.values():
                for g in genes.values():
                    g.quality.calculate_masking()
                    for t in g.transcripts.values():
                        t.quality.calculate_masking()
                        for c in t.CDSs.values():
                            c.quality.calculate_masking()
        self._annot.update()
    
    def calculate_gc_content(self):
        if self._annot.genome is not None:
            for genes in self._annot.chrs.values():
                for g in genes.values():
                    g.quality.calculate_gc_content()
                    for t in g.transcripts.values():
                        t.quality.calculate_gc_content()
                        for c in t.CDSs.values():
                            c.quality.calculate_gc_content()
                            
    def gene_count(self):
        gene_objects = 0
        unique_gene_ids_in_overlaps = set()
        for genes in self._annot.chrs.values():
            gene_objects += len(genes)
            for g in genes.values():
                for o in g.overlaps["self"]:
                    unique_gene_ids_in_overlaps.add(o.id)
        print(f"There are {gene_objects} gene objects and {len(self._annot.all_gene_ids)} genes in all gene ids and {len(unique_gene_ids_in_overlaps)} ids contained in self overlaps.")

    def update(self, output_dir: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "stats", export:bool=False, max_x:int|None=None, quiet:bool=True,
        #deprecated_arguments
        custom_path:str=""):

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        self._annot.update_features(quiet=quiet)
        self._annot.generate_introns()
        if not quiet:
            print(f"\nUpdating stats for {self._annot.id}")

        if not self._annot.contains_protein_sequences:
            if self._annot.genome is not None:
                self._annot.generate_proteins()

        self.calculate_gc_content()

        if export:
            export_folder = self._resolve_output_dir(output_dir=output_dir, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder)
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

        # Other counters
        self.data["single_exon_genes"] = 0
        self.data["mrnas_with_utr_both_sides"] = 0
        self.data["mrnas_with_at_least_one_utr"] = 0
        self.data["introns_in_5utr"] = 0
        self.data["introns_in_3utr"] = 0
        # Total lengths
        self.data["total_length_gene"] = 0
        self.data["total_length_mRNA"] = 0
        self.data["total_length_CDS"] = 0
        self.data["total_length_exon"] = 0
        self.data["total_length_intron"] = 0
        self.data["total_length_5utr"] = 0
        self.data["total_length_3utr"] = 0
        # Longest / shortest trackers (bp)
        _longest_gene       = 0
        _shortest_gene      = None
        _longest_mRNA       = 0
        _shortest_mRNA      = None
        _longest_exon       = 0
        _shortest_exon      = None
        _longest_intron     = 0
        _shortest_intron    = None
        _longest_CDS        = 0
        _shortest_CDS       = None
        _longest_CDS_seg    = 0
        _shortest_CDS_seg   = None

        for genes in self._annot.chrs.values():
            for g in genes.values():
                self.data["mean_transcripts"].append(len(g.transcripts))
                self.data["mean_gene_size"].append(g.size)
                # Total and longest/shortest gene
                self.data["total_length_gene"] += g.size
                if g.size > _longest_gene:
                    _longest_gene = g.size
                if _shortest_gene is None or g.size < _shortest_gene:
                    _shortest_gene = g.size

                for t in g.transcripts.values():
                    if t.main:
                        if t.introns:
                            for i in t.introns:
                                self.data["mean_intron_size"].append(i.size)
                                self.data["total_length_intron"] += i.size
                                if i.size > _longest_intron:
                                    _longest_intron = i.size
                                if _shortest_intron is None or i.size < _shortest_intron:
                                    _shortest_intron = i.size

                        if t.coding:
                            for c in t.CDSs.values():
                                if c.main:
                                    self.data["mean_CDS_size"].append(c.size)
                                    self.data["total_length_CDS"] += c.size
                                    if c.size > _longest_CDS:
                                        _longest_CDS = c.size
                                    if _shortest_CDS is None or c.size < _shortest_CDS:
                                        _shortest_CDS = c.size
                                    for seg in c.CDS_segments:
                                        if seg.size > _longest_CDS_seg:
                                            _longest_CDS_seg = seg.size
                                        if _shortest_CDS_seg is None or seg.size < _shortest_CDS_seg:
                                            _shortest_CDS_seg = seg.size
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
                                                self.data["total_length_3utr"] += u.size
                                            else:
                                                utr_5 = True
                                                u5_size += u.size
                                                self.data["total_length_5utr"] += u.size
                                            u_size += u.size
                                        if utr_5:
                                            self.data["five_prime_UTRs"] += 1
                                        if utr_3:
                                            self.data["three_prime_UTRs"] += 1
                                        if utr_5 and utr_3:
                                            self.data["mrnas_with_utr_both_sides"] += 1
                                        if utr_5 or utr_3:
                                            self.data["mrnas_with_at_least_one_utr"] += 1
                                        self.data["mean_UTR_size"].append(u_size)
                                        self.data["mean_five_prime_UTR_size"].append(u5_size)
                                        self.data["mean_three_prime_UTR_size"].append(u3_size)
                                    # UTR introns: introns not inside the CDS
                                    if t.introns:
                                        cds_start = c.CDS_segments[0].start
                                        cds_end   = c.CDS_segments[-1].end
                                        for i in t.introns:
                                            if not i.intra_coding:
                                                if t.strand == "+":
                                                    if i.end < cds_start:
                                                        self.data["introns_in_5utr"] += 1
                                                    elif i.start > cds_end:
                                                        self.data["introns_in_3utr"] += 1
                                                elif t.strand == "-":
                                                    if i.start > cds_end:
                                                        self.data["introns_in_5utr"] += 1
                                                    elif i.end < cds_start:
                                                        self.data["introns_in_3utr"] += 1
                                    if self._annot.genome is not None and c.protein is not None:
                                        self.data["mean_protein_size"].append(c.protein.size)
                            self.data["coding_genes"].append(g.id)
                        else:
                            self.data["noncoding_genes"].append(g.id)

                        for e in t.exons:
                            self.data["mean_exon_size"].append(e.size)
                            self.data["total_length_exon"] += e.size
                            if e.size > _longest_exon:
                                _longest_exon = e.size
                            if _shortest_exon is None or e.size < _shortest_exon:
                                _shortest_exon = e.size

                        n_exons = len(t.exons)
                        self.data["mean_exons"].append(n_exons)
                        if n_exons == 1:
                            self.data["single_exon_genes"] += 1

                        t_size = t.size
                        self.data["mean_transcript_size"].append(t_size)
                        self.data["total_length_mRNA"] += t_size
                        if t_size > _longest_mRNA:
                            _longest_mRNA = t_size
                        if _shortest_mRNA is None or t_size < _shortest_mRNA:
                            _shortest_mRNA = t_size

        # Finalise longest/shortest as plain bp values
        self.data["longest_gene"]         = _longest_gene
        self.data["shortest_gene"]        = _shortest_gene if _shortest_gene is not None else 0
        self.data["longest_mRNA"]         = _longest_mRNA
        self.data["shortest_mRNA"]        = _shortest_mRNA if _shortest_mRNA is not None else 0
        self.data["longest_exon"]         = _longest_exon
        self.data["shortest_exon"]        = _shortest_exon if _shortest_exon is not None else 0
        self.data["longest_intron"]       = _longest_intron
        self.data["shortest_intron"]      = _shortest_intron if _shortest_intron is not None else 0
        self.data["longest_CDS"]          = _longest_CDS
        self.data["shortest_CDS"]         = _shortest_CDS if _shortest_CDS is not None else 0
        self.data["longest_CDS_segment"]  = _longest_CDS_seg
        self.data["shortest_CDS_segment"] = _shortest_CDS_seg if _shortest_CDS_seg is not None else 0

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
                # Build a human-readable label for the key
                if key in KEY_DESCRIPTIONS:
                    label = KEY_DESCRIPTIONS[key]
                elif key.startswith("intron-exon boundary:"):
                    boundary = key.split(": ", 1)[1]
                    label = f"Intron-exon boundary {boundary} count"
                else:
                    # Generic fallback: replace underscores with spaces and capitalise
                    label = key.replace("_", " ").capitalize()
                if x == 0:
                    f_out.write(f"{label}\t{value_temp}")
                else:
                    f_out.write(f"\n{label}\t{value_temp}")
            f_out.close()
        if not quiet:
            print(f"\nUpdated stats for {self._annot.id}")
