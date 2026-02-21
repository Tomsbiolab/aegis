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
    def __init__(self, annotation):
        self._annot = annotation

    def calculate_transcript_masking(self, hard_masked_genome:object):
        for genes in self._annot.chrs.values():
            for g in genes.values():
                g.generate_hard_sequence(hard_masked_genome)
                for t in g.transcripts.values():
                    t.generate_hard_sequence(hard_masked_genome)
                    t.calculate_masking()
                    t.clear_sequence(just_hard=True)
                    for c in t.CDSs.values():
                        c.generate_hard_sequence(hard_masked_genome)
                        c.calculate_masking()
                        c.clear_sequence(just_hard=True)
        self._annot.update()
    
    def calculate_gc_content(self):
        for genes in self._annot.chrs.values():
            for g in genes.values():
                g.calculate_gc_content()
                for t in g.transcripts.values():
                    t.calculate_gc_content()
                    for c in t.CDSs.values():
                        c.calculate_gc_content()

    def update_stats(self, custom_path:str="", export:bool=False, genome:object=None, max_x:int=None, quiet:bool=True):
        if not quiet:
            print(f"\nUpdating stats for {self._annot.id}")
        if not self._annot.generated_all_sequences or not self._annot.contains_protein_sequences:
            if genome != None:
                self._annot.generate_sequences(genome, quiet=quiet)
                self._annot.clear_sequences(keep_proteins=True)

        self.calculate_gc_content()

        if export:
            export_folder = Path(custom_path or self._annot.path) / "out_stats"
            export_folder.mkdir(parents=True, exist_ok=True)
            export_folder = str(export_folder) + "/"

        to_tally = ["coding_genes", "noncoding_genes", "CDSs_without_stop", "CDSs_with_stop"]

        self._annot.stats = {"mean_transcripts" : [], "mean_exons" : [], "mean_exon_size" : [], "mean_gene_size" : [], "mean_intron_size" : [], "mean_CDS_size" : [], "mean_UTR_size" : [], "mean_transcript_size" : [], "mean_five_prime_UTR_size" : [], "mean_three_prime_UTR_size" : []}

        if genome != None:
            self._annot.stats["mean_protein_size"] = []

        for key in to_tally:
            self._annot.stats[key] = []

        for ft, value in self._annot.features.items():
            self._annot.stats[ft] = value

        self._annot.stats["five_prime_UTRs"] = 0
        self._annot.stats["three_prime_UTRs"] = 0

        for genes in self._annot.chrs.values():
            for g in genes.values():
                self._annot.stats["mean_transcripts"].append(len(g.transcripts))
                self._annot.stats["mean_gene_size"].append(g.size)
                for t in g.transcripts.values():
                    if t.main:
                        if hasattr(t, "introns"):
                            for i in t.introns:
                                self._annot.stats["mean_intron_size"].append(i.size)

                        if t.coding:
                            for c in t.CDSs.values():
                                if c.main:
                                    self._annot.stats["mean_CDS_size"].append(c.size)
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
                                            self._annot.stats["five_prime_UTRs"] += 1
                                        if utr_3:
                                            self._annot.stats["three_prime_UTRs"] += 1
                                        self._annot.stats["mean_UTR_size"].append(u_size)
                                        self._annot.stats["mean_five_prime_UTR_size"].append(u5_size)
                                        self._annot.stats["mean_three_prime_UTR_size"].append(u3_size)
                                    if genome != None:
                                        self._annot.stats["mean_protein_size"].append(c.protein.size)
                            self._annot.stats["coding_genes"].append(g.id)
                        else:
                            self._annot.stats["noncoding_genes"].append(g.id)

                        for e in t.exons:
                            self._annot.stats["mean_exon_size"].append(e.size)

                        self._annot.stats["mean_exons"].append(len(t.exons))
                        
                        self._annot.stats["mean_transcript_size"].append(t.size)

        # anything with mean will be also plotted as distribution plots:
        if export:
            for key in self._annot.stats:
                if "mean" in key:
                    tag = key.split("mean_")[1]
                    if tag[-1] != "s":
                        tag += "s"
                    barplot(self._annot.stats[key], export_folder, f"{self._annot.id}{self._annot.feature_suffix}_{tag}", f"Distribution of {self._annot.id} {tag}", max_x)   

        if genome != None:
            self._annot.stats["CDSs_without_stop"] = []
            self._annot.stats["CDSs_with_stop"] = []
            self._annot.stats["intron_composition"] = set()
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
                                                self._annot.stats["CDSs_without_stop"].append(g.id)
                                            else:
                                                self._annot.stats["CDSs_with_stop"].append(g.id)
                                            
                            for i in t.introns:
                                self._annot.stats["intron_composition"].add(i.boundary)
                                if i.canonical:
                                    intron_stats[f"intron-exon boundary: {i.boundary}"] += 1
                                else:
                                    intron_stats["other_intron_seqs"] += 1
            self._annot.stats["intron_composition"] = list(self._annot.stats["intron_composition"])

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
            self._annot.stats.update(intron_stats)

        for key in self._annot.stats:
            if "mean" in key:
                if len(self._annot.stats[key]) > 0:
                    self._annot.stats[key] = mean(self._annot.stats[key])
                else:
                    self._annot.stats[key] = 0
            #tallying
            elif key in to_tally:
                self._annot.stats[key] = len(self._annot.stats[key])

        if export:
            if genome != None:
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
            for key, value in self._annot.stats.items():
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
