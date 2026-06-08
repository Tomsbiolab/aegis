from __future__ import annotations

import matplotlib.pyplot as plt
from scipy.stats import fisher_exact
import pandas as pd

from ..utils.genefunctions import reverse_complement
from ..utils.misc import find_all_occurrences
from .base import AnnotationComponent
from ..utils.misc import start_progress_bar

class AnnotationMotifs (AnnotationComponent):
    """
    Component for handling motif methods for the Annotation class.
    Accessed via 'annotation_object.motifs'.
    """

    def find_and_plot(self, query_genes:list[str], motif:str, motif_length:int, glistname:str, tf_motif_tag:str, backlist:list[str]=[], backlistname:str="", filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, subfolder_name: str = "motifs", subfolder: bool = False, use_annot_dir: bool = False, quiet:bool=False):

        if subfolder_name != "motifs":
            subfolder = True

        bin_division = 30
        bins_genome_division = 30

        if backlist != [] and not backlistname:
            backlistname = "custom_background"

        if motif_length < 4 or len(motif) < 4:
            raise ValueError(f"Chosen motif={motif} is too short (len={motif_length}) for promoter search.")

        random_ids = self._annot.return_random_gene_ids(len(query_genes), to_avoid=query_genes)
        if backlist == []:
            total = (len(query_genes) * 2) + len(self._annot.all_gene_ids.keys())
        else:
            total = (len(query_genes) * 2) + len(backlist)

        progress_bar = start_progress_bar(total=total, description=f"Scanning {glistname} genes for {tf_motif_tag} ({motif})", colour="94", quiet=quiet)

        extra_suffixes = [tf_motif_tag]

        motif = motif.upper()
        promoter_length = 0

        for genes in self._annot.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    if t.main:
                        if t.promoter:
                            if t.promoter.seq != "":
                                promoter_length = t.promoter.size
                                break
                break
            break
                    
        # critical thing to understand here is that we are looking at how a motif
        # is oriented with regards to the TSS
        towards_occurrences = []
        against_occurrences = []
        interest_proportion = 0
        avg_motifs_interest = []
        for id in query_genes:
            progress_bar.update(1)
            ch = self._annot.all_gene_ids[id]
            g = self._annot.chrs[ch][id]
            for t in g.transcripts.values():
                if t.main:
                    p = t.promoter
                    occurrences_t = find_all_occurrences(motif, p.seq)
                    occurrences_a = find_all_occurrences(motif, reverse_complement(p.seq)) #type: ignore
                    occurrence_count_total = len(occurrences_t) + len(occurrences_a)
                    if occurrence_count_total != 0:
                        avg_motifs_interest.append(occurrence_count_total)
                    if occurrences_t != [] or occurrences_a != []:
                        interest_proportion += 1
                    for m in occurrences_t:
                        towards_occurrences.append((m[0], m[1], p.start + m[0], p.start + m[1], m[2], "same", g.id, str(g.names), g.strand, g.ch))
                    for m in occurrences_a:
                        against_occurrences.append((m[0], m[1], p.end - m[0], p.end - m[1], m[2], "different", g.id, str(g.names), g.strand, g.ch))


        midpoints = []
        for o in towards_occurrences:
            midpoint = (o[0] + o[1]) // 2
            midpoint -= promoter_length
            midpoints.append(midpoint)

        for o in against_occurrences:
            midpoint = (o[0] + o[1]) // 2
            midpoint += 1
            midpoint = (midpoint * -1)
            midpoints.append(midpoint)


        all_occurrences = towards_occurrences + against_occurrences
        df = pd.DataFrame(all_occurrences, columns=["start", "end", "genomic_start", "genomic_end", "sequence", "orientation with respect to gene", "gene_id", "gene_name", "gene_strand", "chromosome"])
        if backlist == []:

            backlist_suffixes = extra_suffixes.copy()
            backlist_suffixes.append(glistname)

            backlist_file_path = self._resolve_output_path(output_dir=output_dir, subfolder_name=subfolder_name, subfolder=subfolder, use_annot_dir=use_annot_dir, filepath=filepath, extra_suffixes=backlist_suffixes, filename=filename, extension=".csv")

            df.to_csv(str(backlist_file_path), sep="\t")
        
        interest_count = len(midpoints)
        plt.hist(midpoints, bins=(promoter_length//bin_division), color='skyblue', edgecolor='skyblue')
        plt.grid(False)
        plt.title(f"{self._annot.id} {tf_motif_tag} {glistname} histogram\npromoters with motif: ({interest_proportion}/{len(query_genes)})")
        plt.xlabel("promoter position")
        plt.ylabel(f"motif occurrence count (total: {interest_count})")
        plt.grid(True)
        
        if backlist == []:

            backlist_file_path = self._resolve_output_path(output_dir=output_dir, subfolder_name=subfolder_name, subfolder=subfolder, use_annot_dir=use_annot_dir, filepath=filepath, extra_suffixes=backlist_suffixes, filename=filename, extension=".pdf")
            
            plt.savefig(str(backlist_file_path))

        plt.close()

        # critical thing to understand here is that we are looking at how a motif
        # is oriented with regards to the TSS
        towards_occurrences = []
        against_occurrences = []
        avg_motifs_random = []
        random_proportion = 0
        for id in random_ids:
            progress_bar.update(1)
            ch = self._annot.all_gene_ids[id]
            g = self._annot.chrs[ch][id]
            for t in g.transcripts.values():
                if t.main:
                    p = t.promoter
                    occurrences_t = find_all_occurrences(motif, p.seq)
                    occurrences_a = find_all_occurrences(motif, reverse_complement(p.seq)) #type: ignore
                    occurrence_count_total = len(occurrences_t) + len(occurrences_a)
                    if occurrence_count_total != 0:
                        avg_motifs_random.append(occurrence_count_total)
                    if occurrences_t != [] or occurrences_a != []:
                        random_proportion += 1
                    for m in occurrences_t:
                        towards_occurrences.append((m[0], m[1], p.start + m[0], p.start + m[1], m[2], "same", g.id, str(g.names), g.strand, g.ch))
                    for m in occurrences_a:
                        against_occurrences.append((m[0], m[1], p.end - m[0], p.end - m[1], m[2], "different", g.id, str(g.names), g.strand, g.ch))

        midpoints = []
        for o in towards_occurrences:
            midpoint = (o[0] + o[1]) // 2
            midpoint -= promoter_length
            midpoints.append(midpoint)

        for o in against_occurrences:
            midpoint = (o[0] + o[1]) // 2
            midpoint += 1
            midpoint = (midpoint * -1)
            midpoints.append(midpoint)


        all_occurrences = towards_occurrences + against_occurrences
        df = pd.DataFrame(all_occurrences, columns=["start", "end", "genomic_start", "genomic_end", "sequence", "orientation with respect to gene", "gene_id", "gene_name", "gene_strand", "chromosome"])
        if backlist == []:

            backlist_suffixes.append("random")

            backlist_file_path = self._resolve_output_path(output_dir=output_dir, subfolder_name=subfolder_name, subfolder=subfolder, use_annot_dir=use_annot_dir, filepath=filepath, extra_suffixes=backlist_suffixes, filename=filename, extension=".csv")

            df.to_csv(str(backlist_file_path), sep="\t")

        random_count = len(midpoints)
        plt.hist(midpoints, bins=(promoter_length//bin_division), color='grey', edgecolor='grey')
        plt.grid(False)
        plt.title(f"{self._annot.id} {tf_motif_tag} {glistname} random histogram\npromoters with motif: ({random_proportion}/{len(query_genes)})")
        plt.xlabel("promoter position")
        plt.ylabel(f"motif occurrence count (total: {random_count})")
        plt.grid(True)
        if backlist == []:
            backlist_file_path = self._resolve_output_path(output_dir=output_dir, subfolder_name=subfolder_name, subfolder=subfolder, use_annot_dir=use_annot_dir, filepath=filepath, extra_suffixes=backlist_suffixes, filename=filename, extension=".pdf")

            plt.savefig(str(backlist_file_path))
        plt.close()

        # critical thing to understand here is that we are looking at how a motif
        # is oriented with regards to the TSS
        towards_occurrences = []
        against_occurrences = []
        avg_motifs_genomic = []
        genomic_proportion = 0
        if backlist == []:
            for id in self._annot.all_gene_ids:
                progress_bar.update(1)
                ch = self._annot.all_gene_ids[id]
                g = self._annot.chrs[ch][id]
                for t in g.transcripts.values():
                    if t.main:
                        p = t.promoter
                        occurrences_t = find_all_occurrences(motif, p.seq)
                        occurrences_a = find_all_occurrences(motif, reverse_complement(p.seq)) #type: ignore
                        occurrence_count_total = len(occurrences_t) + len(occurrences_a)
                        if occurrence_count_total != 0:
                            avg_motifs_genomic.append(occurrence_count_total)
                        if occurrences_t != [] or occurrences_a != []:
                            genomic_proportion += 1
                        for m in occurrences_t:
                            towards_occurrences.append((m[0], m[1], p.start + m[0], p.start + m[1], m[2], "same", g.id, str(g.names), g.strand, g.ch))
                        for m in occurrences_a:
                            against_occurrences.append((m[0], m[1], p.end - m[0], p.end - m[1], m[2], "different", g.id, str(g.names), g.strand, g.ch))

            midpoints = []
            for o in towards_occurrences:
                midpoint = (o[0] + o[1]) // 2
                midpoint -= promoter_length
                midpoints.append(midpoint)

            for o in against_occurrences:
                midpoint = (o[0] + o[1]) // 2
                midpoint += 1
                midpoint = (midpoint * -1)
                midpoints.append(midpoint)

            genomic_count = len(midpoints)
            plt.hist(midpoints, bins=(promoter_length//bins_genome_division), color='grey', edgecolor='grey')
            plt.grid(False)
            plt.title(f"{self._annot.id} {tf_motif_tag} full genome histogram\npromoters with motif: ({genomic_proportion}/{len(self._annot.all_gene_ids.keys())})")
            plt.xlabel("promoter position")
            plt.ylabel(f"motif occurrence count (total: {genomic_count})")
            plt.grid(True)

            whole_genome_suffixes = extra_suffixes.copy()
            whole_genome_suffixes.append("whole_genome")

            whole_genome_file_path = self._resolve_output_path(output_dir=output_dir, subfolder_name=subfolder_name, subfolder=subfolder, use_annot_dir=use_annot_dir, filepath=filepath, extra_suffixes=whole_genome_suffixes, filename=filename, extension=".pdf")

            plt.savefig(str(whole_genome_file_path))
            plt.close()
        
        else:
            for id in backlist:
                progress_bar.update(1)
                ch = self._annot.all_gene_ids[id]
                g = self._annot.chrs[ch][id]
                for t in g.transcripts.values():
                    if t.main:
                        p = t.promoter
                        occurrences_t = find_all_occurrences(motif, p.seq)
                        occurrences_a = find_all_occurrences(motif, reverse_complement(p.seq)) #type: ignore
                        occurrence_count_total = len(occurrences_t) + len(occurrences_a)
                        if occurrence_count_total != 0:
                            avg_motifs_genomic.append(occurrence_count_total)
                        if occurrences_t != [] or occurrences_a != []:
                            genomic_proportion += 1
                        for m in occurrences_t:
                            towards_occurrences.append((m[0], m[1], p.start + m[0], p.start + m[1], m[2], "same", g.id, str(g.names), g.strand, g.ch))
                        for m in occurrences_a:
                            against_occurrences.append((m[0], m[1], p.end - m[0], p.end - m[1], m[2], "different", g.id, str(g.names), g.strand, g.ch))

            midpoints = []
            for o in towards_occurrences:
                midpoint = (o[0] + o[1]) // 2
                midpoint -= promoter_length
                midpoints.append(midpoint)

            for o in against_occurrences:
                midpoint = (o[0] + o[1]) // 2
                midpoint += 1
                midpoint = (midpoint * -1)
                midpoints.append(midpoint)

            genomic_count = len(midpoints)
            plt.hist(midpoints, bins=(promoter_length//bins_genome_division), color='grey', edgecolor='grey')
            plt.grid(False)
            plt.title(f"{self._annot.id} {tf_motif_tag} {backlistname} as background histogram\npromoters with motif: ({genomic_proportion}/{len(backlist)})")
            plt.xlabel("promoter position")
            plt.ylabel(f"motif occurrence count (total: {genomic_count})")
            plt.grid(True)

            background_suffixes = extra_suffixes.copy()
            background_suffixes.append(f"{backlistname}_as_background")            

            back_file_path = self._resolve_output_path(output_dir=output_dir, subfolder_name=subfolder_name, subfolder=subfolder, use_annot_dir=use_annot_dir, filepath=filepath, extra_suffixes=background_suffixes, filename=filename, extension=".pdf")

            plt.savefig(str(back_file_path))
            plt.close()

        progress_bar.close()

        # Counts of non-motif occurrences
        interest_non_count = (len(query_genes) * (int(promoter_length/motif_length)) * 2) - interest_count
        if backlist == []:
            genomic_non_count = (len(self._annot.all_gene_ids.keys()) * (int(promoter_length/motif_length)) * 2) - genomic_count
        else:
            genomic_non_count = (len(backlist) * (int(promoter_length/motif_length)) * 2) - genomic_count
        interest_non_proportion = len(query_genes) - interest_proportion
        if backlist == []:
            genomic_non_proportion = len(self._annot.all_gene_ids.keys()) - genomic_proportion
        else:
            genomic_non_proportion = len(backlist) - genomic_proportion

        odds_ratio_occurrences, p_value_occurrences = fisher_exact([[interest_count, genomic_count], 
                                                                    [interest_non_count, genomic_non_count]])

        odds_ratio_proportion, p_value_proportion = fisher_exact([[interest_proportion, genomic_proportion], 
                                                                  [interest_non_proportion, genomic_non_proportion]])

        promoter_percentage_interest = (interest_proportion / len(query_genes)) * 100
        if backlist == []:
            promoter_percentage_genome = (genomic_proportion / len(self._annot.all_gene_ids.keys())) * 100
        else:
            promoter_percentage_genome = (genomic_proportion / len(backlist)) * 100

        output_file = self._resolve_output_path(output_dir=output_dir, subfolder_name=subfolder_name, subfolder=subfolder, use_annot_dir=use_annot_dir, filepath=filepath, filename=filename, extension=".csv")

        return interest_count, genomic_count, p_value_occurrences, odds_ratio_occurrences, promoter_percentage_interest, promoter_percentage_genome, p_value_proportion, odds_ratio_proportion, avg_motifs_interest, avg_motifs_genomic, output_file
