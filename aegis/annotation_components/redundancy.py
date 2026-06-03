from __future__ import annotations
from typing import TYPE_CHECKING

import networkx as nx
import sys
from tqdm import tqdm

if TYPE_CHECKING:
    from ..annotation import Annotation
    from ..genome import Genome

class AnnotationRedundancy:
    """
    Component for handling redundancy removal for the Annotation class.
    Accessed via 'annotation_object.redundancy'.
    """
    _annot: Annotation
    def __init__(self, annotation: Annotation):
        self._annot = annotation

    def mark_noisy_genes(self, protein_size:int=50, intron_size:int=100000, remove_noncoding:bool=True, remove_masked:bool=True, quiet:bool=False):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self._annot.id} noisy genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self._annot.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if remove_masked:
                    if g.quality.masked_fraction == 1:
                        g.quality.remove = True
                        continue
                if remove_noncoding:
                    if not g.coding:
                        g.quality.remove = True
                        continue
                for t in g.transcripts.values():
                    if t.main:
                        if t.quality.masked_fraction == 1:
                            g.quality.remove = True
                            break
                        for c in t.CDSs.values():
                            if c.main:
                                if c.size < (protein_size * 3) or c.quality.masked_fraction == 1:
                                    g.quality.remove = True
                                    break
                        if t.introns:
                            for i in t.introns:
                                if i.size > intron_size:
                                    g.quality.remove = True
                                    break
                        break

        progress_bar.close()

    def mark_transcriptomic_supported_genes(self, quiet:bool=False):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self._annot.id} transcriptomic supported genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self._annot.chrs.values():
            for g in genes.values():
                g.quality.transcriptomic_evidence = False
                progress_bar.update(1)
                if "stringtie" in g.source or "psiclass" in g.source:
                    g.quality.transcriptomic_evidence = True
                else:
                    for o in g.overlaps["self"]:
                        if o.score >= 5:
                            if o.CDS_query_percent > 30:
                                if "stringtie" in self._annot.chrs[g.ch][o.id].source or "psiclass" in self._annot.chrs[g.ch][o.id].source:
                                    g.quality.transcriptomic_evidence = True
        progress_bar.close()

    def mark_abinitio_supported_genes(self, reliable_sources:list=["augustus", "genemark"], quiet:bool=False):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self._annot.id} abinitio supported genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self._annot.chrs.values():
            for g in genes.values():
                g.quality.abinitio_evidence = False
                progress_bar.update(1)
                if g.source in reliable_sources:
                    g.quality.abinitio_evidence = True
                else:
                    for o in g.overlaps["self"]:
                        if o.score >= 5:
                            if o.CDS_query_percent > 30:
                                if self._annot.chrs[g.ch][o.id].source in reliable_sources:
                                    g.quality.abinitio_evidence = True
        progress_bar.close()

    def add_reliable_CDS_evidence_score(self, quiet:bool=False):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self._annot.id} reliable CDS evidences:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

        for genes in self._annot.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if not g.quality.remove:
                    for o in g.overlaps["self"]:
                        if not self._annot.chrs[g.ch][o.id].quality.remove:
                            if o.score == 11:
                                g.quality.reliable_score += 1

        for genes in self._annot.chrs.values():
            for g in genes.values():
                if g.quality.reliable_score == 0:
                    g.quality.remove = True
        progress_bar.close()

    def mark_reliable_CDS_evidences(self, unreliable_sources:list=["GlimmerHMM", "geneid_v1.4"], quiet:bool=False):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self._annot.id} reliable CDS evidences:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

        for genes in self._annot.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if not g.quality.remove:
                    for o in g.overlaps["self"]:
                        if not self._annot.chrs[g.ch][o.id].quality.remove:
                            if (self._annot.chrs[g.ch][o.id].source not in unreliable_sources) or (g.source not in unreliable_sources):
                                if o.score == 11:
                                    g.quality.reliable = True
                                    if not self._annot.chrs[g.ch][o.id].quality.reliable:
                                        self._annot.chrs[g.ch][o.id].quality.remove = True
                                    self._annot.chrs[g.ch][o.id].quality.reliable = True

        for genes in self._annot.chrs.values():
            for g in genes.values():
                if not g.quality.reliable:
                    g.quality.remove = True
        progress_bar.close()

    def find_best_gene_model(self, source_priority:list, just_with_reliables:bool=True, quiet:bool=False):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self._annot.id} noisy genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        
        if just_with_reliables:
            for genes in self._annot.chrs.values():
                for g in genes.values():
                    progress_bar.update(1)
                    if not g.quality.remove:
                        for o in g.overlaps["self"]:
                            if o.score >= 5 or (o.score == 1 and o.antiscore >= 5):
                                if not self._annot.chrs[g.ch][o.id].quality.remove and not g.quality.remove:
                                    if g.quality.reliable_score > self._annot.chrs[g.ch][o.id].quality.reliable_score:
                                        self._annot.chrs[g.ch][o.id].quality.remove = True
                                    elif self._annot.chrs[g.ch][o.id].quality.reliable_score > g.quality.reliable_score:
                                        g.quality.remove = True
                                    else:
                                        query_best = g.compare_protein_blast_hits(self._annot.chrs[g.ch][o.id], source_priority)
                                        if query_best:
                                            self._annot.chrs[g.ch][o.id].quality.remove = True
                                        else:
                                            g.quality.remove = True

        else:
            for genes in self._annot.chrs.values():
                for g in genes.values():
                    progress_bar.update(1)
                    # mark unrescuable genes in different strand to a reliable gene and overlapping exon with cds of the reliable gene
                    if g.quality.remove and not g.quality.overlap_reliable and g.quality.transcriptomic_evidence and not g.quality.unrescuable:
                        for o in g.overlaps["self"]:
                            if not self._annot.chrs[g.ch][o.id].quality.remove:
                                if o.score == 1:
                                    for t3 in g.transcripts.values():
                                        for t4 in self._annot.chrs[g.ch][o.id].transcripts.values():
                                            for e1 in t3.exons:
                                                for c1 in t4.CDSs.values():
                                                    for cs1 in c1.CDS_segments:
                                                        overlapping, _ = e1.overlap(cs1)
                                                        if overlapping:
                                                            g.quality.unrescuable = True
                                                            g.quality.rescue = False
                                                            break

                            if self._annot.chrs[g.ch][o.id].quality.remove and not self._annot.chrs[g.ch][o.id].quality.overlap_reliable and self._annot.chrs[g.ch][o.id].quality.transcriptomic_evidence and not self._annot.chrs[g.ch][o.id].quality.unrescuable and not g.quality.unrescuable:
                                if o.score >= 5 or (o.score == 1 and o.antiscore >= 5):
                                    query_best = g.compare_protein_blast_hits(self._annot.chrs[g.ch][o.id], source_priority)
                                    if query_best:
                                        g.quality.rescue = True
                                        self._annot.chrs[g.ch][o.id].quality.rescue = False
                                        self._annot.chrs[g.ch][o.id].quality.remove = True
                                        self._annot.chrs[g.ch][o.id].quality.unrescuable = True
                                    else:
                                        self._annot.chrs[g.ch][o.id].quality.rescue = True
                                        g.quality.rescue = False
                                        g.quality.remove = True
                                        g.quality.unrescuable = True

        progress_bar.close()

    def mark_overlap_with_reliable_genes(self, quiet:bool=False):
                # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self._annot.id} overlap with reliable genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self._annot.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                for o in g.overlaps["self"]:
                    if o.score >= 5 or (o.score == 1 and o.antiscore >= 5):
                        if not self._annot.chrs[g.ch][o.id].quality.remove:
                            g.quality.overlap_reliable = True
        progress_bar.close()

    def add_better_ab_initio_models_as_alternative_transcripts(self, source_priority, reliable_sources:list=["augustus", "genemark"], quiet:bool=False):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys())*3, disable=disable,
                                bar_format=(
                    f'\033[1;91mSelecting {self._annot.id} alternative transcripts:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

        for genes in self._annot.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                g_coding_ratio = 0
                full_UTR_exons = 0
                for t in g.transcripts.values():
                    if t.main:
                        g_coding_ratio = t.coding_ratio
                        for c in t.CDSs.values():
                            if c.main:
                                full_UTR_exons = c.full_UTR_exons

                if (not g.quality.remove or g.quality.rescue) and g_coding_ratio < 0.7 and full_UTR_exons > 0:
                    for o in g.overlaps["self"]:
                        if o.score >= 5:
                            if self._annot.chrs[g.ch][o.id].quality.remove and not self._annot.chrs[g.ch][o.id].quality.rescue and self._annot.chrs[g.ch][o.id].source in reliable_sources:
                                query_best = g.compare_protein_blast_hits(self._annot.chrs[g.ch][o.id], source_priority)
                                if not query_best:
                                    g.alternative_transcript_rescue.append(o.id)
        
        # cluster genes with alternative transcript rescue
        gene_groups = []

        for genes in self._annot.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if (not g.quality.remove or g.quality.rescue) and g.alternative_transcript_rescue != []:
                    for g2 in genes.values():
                        if g.id == g2.id:
                            continue
                        if (not g2.quality.remove or g2.quality.rescue) and g2.alternative_transcript_rescue != []:
                            if bool(set(g.alternative_transcript_rescue).intersection(set(g2.alternative_transcript_rescue))):
                                temp_group = set(g.alternative_transcript_rescue) | set(g2.alternative_transcript_rescue)
                                temp_group.add(g.id)
                                temp_group.add(g2.id)
                                found = False
                                index = 0
                                for n, group in enumerate(gene_groups):
                                    if bool(temp_group.intersection(group)):
                                        found = True
                                        index = n
                                        break
                                if found:
                                    gene_groups[index] = gene_groups[index] | temp_group
                                else:
                                    gene_groups.append(temp_group)


        merge_genes = set()

        for gene_set in gene_groups:
            merge_genes = merge_genes | gene_set
            best_reliable = ""
            best_unreliable = ""
            for g_id in gene_set:
                chrom = self._annot.all_gene_ids[g_id]
                if not self._annot.chrs[chrom][g_id].quality.remove or self._annot.chrs[chrom][g_id].quality.rescue:
                    best_reliable = g_id
                else:
                    best_unreliable = g_id
            
            for g_id in gene_set:
                chrom = self._annot.all_gene_ids[g_id]
                if not self._annot.chrs[chrom][g_id].quality.remove or self._annot.chrs[chrom][g_id].quality.rescue:
                    query_best = self._annot.chrs[chrom][g_id].compare_protein_blast_hits(self._annot.chrs[chrom][best_reliable], source_priority)
                    if query_best:
                        best_reliable = g_id
                else:
                    query_best = self._annot.chrs[chrom][g_id].compare_protein_blast_hits(self._annot.chrs[chrom][best_unreliable], source_priority)
                    if query_best:
                        best_unreliable = g_id

            for g_id in gene_set:
                chrom = self._annot.all_gene_ids[g_id]

                if g_id == best_reliable:
                    continue
                elif g_id == best_unreliable:
                    for t in self._annot.chrs[chrom][best_unreliable].transcripts.values():
                        if t.main:
                            t_copy = t.copy()
                            t_copy.id = "alternative_transcript"
                            t_copy.parents = [g_id]
                            t_copy.symbols = []
                            t_copy.names = []
                            t_copy.synonyms = []
                    self._annot.chrs[chrom][best_reliable].transcripts["alternative_transcript"] = t_copy.copy()
                    del t_copy
                else:
                    self._annot.chrs[chrom][g_id].quality.rescue = False
                    self._annot.chrs[chrom][g_id].quality.remove = True

        for genes in self._annot.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if g.id not in merge_genes and g.alternative_transcript_rescue != []:
                    best = g.alternative_transcript_rescue[0]
                    for alt in g.alternative_transcript_rescue:
                        if alt == best:
                            continue
                        query_best = g.compare_protein_blast_hits(self._annot.chrs[g.ch][alt], source_priority)
                        if not query_best:
                            best = alt

                    for t in self._annot.chrs[g.ch][best].transcripts.values():
                        if t.main:
                            t_copy = t.copy()
                            t_copy.id = "alternative_transcript"
                            t_copy.parents = [g.id]
                            t_copy.symbols = []
                            t_copy.names = []
                            t_copy.synonyms = []
                    g.transcripts[t_copy.id] = t_copy.copy()
                    del t_copy
        progress_bar.close()
        self._annot.update(rename_features=["transcript", "CDS", "exon", "UTR"])

    def rescue_longer_same_frame_CDS(self, reliable_sources:list[str]=["augustus", "genemark"], quiet:bool=False):

        for genes in self._annot.chrs.values():
            for g in genes.values():
                if not g.quality.remove or g.quality.rescue:

                    main_CDS_size = 0

                    for t in g.transcripts.values():
                        if t.main:
                            for c in t.CDSs.values():
                                if c.main:
                                    main_CDS_size = c.size

                    posible_alternative_transcripts = []

                    for o in g.overlaps["self"]:
                        if o.CDSs_in_both:

                            overlap_main_CDS_size = 0

                            for t in self._annot.chrs[g.ch][o.id].transcripts.values():
                                if t.main:
                                    for c in t.CDSs.values():
                                        if c.main:
                                            overlap_main_CDS_size = c.size

                            if (self._annot.chrs[g.ch][o.id].source in reliable_sources) and ((o.full_protein_overlaps >= 2) or (o.protein_query_percent >= 90)) and (overlap_main_CDS_size > main_CDS_size):

                                for t in self._annot.chrs[g.ch][o.id].transcripts.values():
                                    if t.main:
                                        t_copy = t.copy()
                                        t_copy.id = "alternative_transcript2"
                                        t_copy.parents = [g.id]
                                        t_copy.symbols = []
                                        t_copy.names = []
                                        t_copy.synonyms = []

                                posible_alternative_transcripts.append(t_copy.copy())
                                del t_copy

                    best_candidate_size = 0
                    best_candidate = None

                    if posible_alternative_transcripts:

                        for t_candidate in posible_alternative_transcripts:

                            if t_candidate.size > best_candidate_size:

                                best_candidate = t_candidate
                                best_candidate_size = t_candidate.size

                        if best_candidate:
                            g.transcripts[best_candidate.id] = best_candidate.copy()

        self._annot.update(rename_features=["transcript", "CDS", "exon", "UTR"], quiet=quiet)

    def remove_CDS_overlaps(self, source_priority, blast:bool=False, anti:bool=True):
        for genes in self._annot.chrs.values():
            for g in genes.values():
                if not g.quality.remove or g.quality.rescue:
                    for o in g.overlaps["self"]:
                        overlapping_CDS = False
                        if anti:
                            if o.score >= 5 or (o.score == 1 and o.antiscore >= 5):
                                overlapping_CDS = True
                        else:
                            if o.score >= 5:
                                overlapping_CDS = True

                        if overlapping_CDS:
                            if (not self._annot.chrs[g.ch][o.id].quality.remove or self._annot.chrs[g.ch][o.id].quality.rescue) and (not g.quality.remove or g.quality.rescue):
                                if blast:
                                    query_best = g.compare_protein_blast_hits(self._annot.chrs[g.ch][o.id], source_priority)
                                else:
                                    query_best = g.longer_CDS(self._annot.chrs[g.ch][o.id])
                                if query_best:
                                    self._annot.chrs[g.ch][o.id].quality.remove = True
                                    self._annot.chrs[g.ch][o.id].quality.rescue = False
                                else:
                                    g.quality.remove = True
                                    g.quality.rescue = False

    def mark_intron_nesting(self, ignore_removed:bool=True):
        for genes in self._annot.chrs.values():
            for g in genes.values():
                for o in g.overlaps["self"]:
                    if o.score < 5:
                        chrom = self._annot.all_gene_ids[o.id]
                        if ignore_removed:
                            if self._annot.chrs[chrom][o.id].quality.remove and not self._annot.chrs[chrom][o.id].quality.rescue:
                                continue
                        if (o.exon_query_percent == 0 and o.exon_target_percent == 0) and (g.start > self._annot.chrs[chrom][o.id].start or g.end < self._annot.chrs[chrom][o.id].end):
                            g.quality.intron_nested = True
                            if self._annot.chrs[chrom][o.id].start < g.start and self._annot.chrs[chrom][o.id].end > g.end:
                                g.quality.intron_nested_fully_contained = True

                            target_cds = self._annot.chrs[chrom][o.id].get_main_CDS_range()
                            query_cds = g.get_main_CDS_range()


                            if target_cds and query_cds:
                                c_start_target, c_end_target = target_cds
                                c_start_query, c_end_query = query_cds
                                
                                if c_start_query > c_end_target or c_end_query < c_start_target:
                                    # UTR intron nested means that a main CDS of a gene finishes and starts outside of the overlaped gene's CDS region
                                    g.quality.UTR_intron_nested = True

                            for t in self._annot.chrs[chrom][o.id].transcripts.values():
                                if t.main:
                                    if t.introns:
                                        for i in t.introns:
                                            if i.start < g.start and i.end > g.end:
                                                g.quality.intron_nested_single = True
                                                break
                                        break

    def remove_fully_intron_nested_genes(self):
        for genes in self._annot.chrs.values():
            for g in genes.values():
                if not g.quality.remove or g.quality.rescue:
                    if g.quality.intron_nested and not g.quality.UTR_intron_nested:
                        for t in g.transcripts.values():
                            if t.main:
                                for c in t.CDSs.values():
                                    if c.main:
                                        if c.size <= 450 or t.coding_ratio < 0.4:
                                            g.quality.remove = True
                                            g.quality.rescue = False

    def mark_overlap_with_other_selected_CDSs(self, quiet:bool=False):
        self._annot.overlaps.clear_with_selected_CDSs()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self._annot.id} overlap with other selected CDSs:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))        
        for genes in self._annot.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if not g.quality.remove or g.quality.rescue:
                    for o in g.overlaps["self"]:
                        if o.score < 5 and (self._annot.chrs[g.ch][o.id].quality.rescue or not self._annot.chrs[g.ch][o.id].quality.remove):
                            for t in g.transcripts.values():
                                for t2 in self._annot.chrs[g.ch][o.id].transcripts.values():
                                    for e in t.exons:
                                        for c in t2.CDSs.values():
                                            for cs in c.CDS_segments:
                                                overlapping, _ = e.overlap(cs)
                                                if overlapping:
                                                    g.quality.overlap_with_selected_CDS = True
                                                    break
                    for o in g.overlaps["self"]:
                        if o.score == 11:
                            for o2 in self._annot.chrs[g.ch][o.id].overlaps["self"]:
                                if o2.score < 5 and (not self._annot.chrs[g.ch][o2.id].quality.remove or self._annot.chrs[g.ch][o2.id].quality.rescue) and o2.id != g.id:
                                    for t in self._annot.chrs[g.ch][o.id].transcripts.values():
                                        for t2 in self._annot.chrs[g.ch][o2.id].transcripts.values():
                                            for e in t.exons:
                                                for c in t2.CDSs.values():
                                                    for cs in c.CDS_segments:
                                                        overlapping, _ = e.overlap(cs)
                                                        if overlapping:
                                                            self._annot.chrs[g.ch][o.id].quality.overlap_with_selected_CDS = True
                                                            break
        progress_bar.close()

    def mark_overlap_with_other_selected_exons(self, quiet:bool=False):
        self._annot.overlaps.clear_with_selected_exons()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mMark {self._annot.id} overlap with other selected exons:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))   
        for genes in self._annot.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if not g.quality.remove or g.quality.rescue:
                    for o in g.overlaps["self"]:
                        if o.score < 5 and (self._annot.chrs[g.ch][o.id].quality.rescue or not self._annot.chrs[g.ch][o.id].quality.remove):
                            for t in g.transcripts.values():
                                for t2 in self._annot.chrs[g.ch][o.id].transcripts.values():
                                    for e1 in t.exons:
                                        for e2 in t2.exons:
                                            overlapping, _ = e1.overlap(e2)
                                            if overlapping:
                                                g.quality.overlap_with_selected_exon = True
                                                break
                    for o in g.overlaps["self"]:
                        if o.score == 11:
                            for o2 in self._annot.chrs[g.ch][o.id].overlaps["self"]:
                                if o2.score < 5 and (not self._annot.chrs[g.ch][o2.id].quality.remove or self._annot.chrs[g.ch][o2.id].quality.rescue) and o2.id != g.id:
                                    for t in self._annot.chrs[g.ch][o.id].transcripts.values():
                                        for t2 in self._annot.chrs[g.ch][o2.id].transcripts.values():
                                            for e1 in t.exons:
                                                for e2 in t2.exons:
                                                    overlapping, _ = e1.overlap(e2)
                                                    if overlapping:
                                                        self._annot.chrs[g.ch][o.id].quality.overlap_with_selected_exon = True
                                                        break
        progress_bar.close()        

    def select_best_possible_non_overlapping_UTR(self, exon=False, quiet:bool=False):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mSelect {self._annot.id} best possible non overlapping UTR:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        
        if not exon:
            for genes in self._annot.chrs.values():
                for g in genes.values():
                    g.update()
                    progress_bar.update(1)
                    if not g.quality.overlap_with_selected_CDS:
                        sizes = [g.size]
                    else:
                        sizes = []
                    if not g.quality.remove or g.quality.rescue:
                        if len(g.transcripts) < 2:
                            for o in g.overlaps["self"]:
                                if o.score == 11:
                                    if not self._annot.chrs[g.ch][o.id].quality.overlap_with_selected_CDS:
                                        sizes.append(self._annot.chrs[g.ch][o.id].size)

                            rescue_transcripts = {}
                            if sizes != []:
                                max_gene_size = max(sizes)
                                for o in g.overlaps["self"]:
                                    if o.score == 11:
                                        if not self._annot.chrs[g.ch][o.id].quality.overlap_with_selected_CDS and self._annot.chrs[g.ch][o.id].size == max_gene_size and rescue_transcripts == {}:
                                            for t in self._annot.chrs[g.ch][o.id].transcripts.values():
                                                rescue_transcripts[t.id] = t.copy()

                            else:
                                sizes = [g.size]
                                for o in g.overlaps["self"]:
                                    if o.score == 11:
                                        sizes.append(self._annot.chrs[g.ch][o.id].size)
                                min_gene_size = min(sizes)
                                for o in g.overlaps["self"]:
                                    if o.score == 11:
                                        if self._annot.chrs[g.ch][o.id].size == min_gene_size and rescue_transcripts == {}:
                                            for t in self._annot.chrs[g.ch][o.id].transcripts.values():
                                                rescue_transcripts[t.id] = t.copy()

                            if rescue_transcripts != {}:
                                g.transcripts = rescue_transcripts.copy()
                                g.update()
        else:
            for genes in self._annot.chrs.values():
                for g in genes.values():
                    g.update()
                    progress_bar.update(1)
                    if not g.quality.overlap_with_selected_exon:
                        sizes = [g.size]
                    else:
                        sizes = []
                    if not g.quality.remove or g.quality.rescue:
                        if len(g.transcripts) < 2:
                            for o in g.overlaps["self"]:
                                if o.score == 11:
                                    if not self._annot.chrs[g.ch][o.id].quality.overlap_with_selected_exon:
                                        sizes.append(self._annot.chrs[g.ch][o.id].size)

                            rescue_transcripts = {}
                            if sizes != []:
                                max_gene_size = max(sizes)
                                for o in g.overlaps["self"]:
                                    if o.score == 11:
                                        if not self._annot.chrs[g.ch][o.id].quality.overlap_with_selected_exon and self._annot.chrs[g.ch][o.id].size == max_gene_size and rescue_transcripts == {}:
                                            for t in self._annot.chrs[g.ch][o.id].transcripts.values():
                                                rescue_transcripts[t.id] = t.copy()

                            else:
                                sizes = [g.size]
                                for o in g.overlaps["self"]:
                                    if o.score == 11:
                                        sizes.append(self._annot.chrs[g.ch][o.id].size)
                                min_gene_size = min(sizes)
                                for o in g.overlaps["self"]:
                                    if o.score == 11:
                                        if self._annot.chrs[g.ch][o.id].size == min_gene_size and rescue_transcripts == {}:
                                            for t in self._annot.chrs[g.ch][o.id].transcripts.values():
                                                rescue_transcripts[t.id] = t.copy()

                            if rescue_transcripts != {}:
                                g.transcripts = rescue_transcripts.copy()      
                                g.update()      

        progress_bar.close()
        self._annot.rename_ids(quiet=quiet)
        self._annot.remove_duplicate_transcripts(quiet=quiet)
        self._annot.update(rename_features=["transcript", "CDS", "exon", "UTR"], quiet=quiet)

    def find_best_gene_model_nested_overlaps(self, source_priority, blast=False):
        """
        For genes fully contained in other genes which have exon overlap choose best blast hit.
        """
        for genes in self._annot.chrs.values():
            for g in genes.values():
                if not g.quality.remove or g.quality.rescue:
                    for o in g.overlaps["self"]:
                        if (not self._annot.chrs[g.ch][o.id].quality.remove or self._annot.chrs[g.ch][o.id].quality.rescue) and (not g.quality.remove or g.quality.rescue):
                            if o.gene_query_percent >= 100 or o.gene_target_percent >= 100:
                                if o.exon_query_percent > 0:
                                    if blast:
                                        query_best = g.compare_protein_blast_hits(self._annot.chrs[g.ch][o.id], source_priority)
                                    else:
                                        query_best = g.longer_CDS(self._annot.chrs[g.ch][o.id])
                                    if query_best:
                                        self._annot.chrs[g.ch][o.id].quality.remove = True
                                        self._annot.chrs[g.ch][o.id].quality.rescue = False
                                    else:
                                        g.quality.remove = True
                                        g.quality.rescue = False

    def find_best_gene_model_exon_num_overlaps(self, source_priority, blast:bool=False, exon_num:int=2):
        """
        For genes with more than X exons exactly the same.
        """
        for genes in self._annot.chrs.values():
            for g in genes.values():
                if not g.quality.remove or g.quality.rescue:
                    for o in g.overlaps["self"]:
                        if (not self._annot.chrs[g.ch][o.id].quality.remove or self._annot.chrs[g.ch][o.id].quality.rescue) and (not g.quality.remove or g.quality.rescue):
                            if o.full_exon_overlaps >= exon_num:
                                if blast:
                                    query_best = g.compare_protein_blast_hits(self._annot.chrs[g.ch][o.id], source_priority)
                                else:
                                    query_best = g.longer_CDS(self._annot.chrs[g.ch][o.id])
                                if query_best:
                                    self._annot.chrs[g.ch][o.id].quality.remove = True
                                    self._annot.chrs[g.ch][o.id].quality.rescue = False
                                else:
                                    g.quality.remove = True
                                    g.quality.rescue = False

    def remove_exon_overlaps(self, source_priority, blast:bool=False):
        for genes in self._annot.chrs.values():
            for g in genes.values():
                if not g.quality.remove or g.quality.rescue:
                    for o in g.overlaps["self"]:
                        if (not self._annot.chrs[g.ch][o.id].quality.remove or self._annot.chrs[g.ch][o.id].quality.rescue) and (not g.quality.remove or g.quality.rescue):
                            if o.exon_query_percent > 0 or o.exon_target_percent > 0:
                                if blast:
                                    query_best = g.compare_protein_blast_hits(self._annot.chrs[g.ch][o.id], source_priority)
                                else:
                                    query_best = g.longer_CDS(self._annot.chrs[g.ch][o.id])
                                if query_best:
                                    self._annot.chrs[g.ch][o.id].quality.remove = True
                                    self._annot.chrs[g.ch][o.id].quality.rescue = False
                                else:
                                    g.quality.remove = True
                                    g.quality.rescue = False

    def remove_UTRs_from_exon_overlaps(self):

        for genes in self._annot.chrs.values():
            for g in genes.values():
                if not g.quality.remove or g.quality.rescue:
                    for o in g.overlaps["self"]:
                        if (not self._annot.chrs[g.ch][o.id].quality.remove or self._annot.chrs[g.ch][o.id].quality.rescue) and (not g.quality.remove or g.quality.rescue):
                            if o.exon_query_percent > 0 or o.exon_target_percent > 0:
                                self._annot.chrs[g.ch][o.id].clear_UTRs()
                                g.clear_UTRs()

        self._annot.sorted = False

    def _take_snapshot(self):
        snapshot = {}
        for chrom, genes in self._annot.chrs.items():
            for g_id, g in genes.items():
                snapshot[g_id] = {
                    attr: getattr(g.quality, attr, None) for attr in [
                        'remove', 'rescue', 'reliable', 'reliable_score', 'overlap_reliable', 'unrescuable', 
                        'intron_nested', 'intron_nested_fully_contained', 'intron_nested_single', 'UTR_intron_nested', 
                        'overlap_with_selected_CDS', 'overlap_with_selected_exon', 'transcriptomic_evidence', 'abinitio_evidence'
                    ]
                }
                snapshot[g_id]['transcripts'] = tuple(sorted(g.transcripts.keys()))
        return snapshot

    def _print_changes(self, step_name, prev_snapshot, quiet=False):
        current_snapshot = self._take_snapshot()
        if quiet:
            return current_snapshot
        
        # Find deleted genes
        deleted_genes = set(prev_snapshot.keys()) - set(current_snapshot.keys())
        # Find new genes
        new_genes = set(current_snapshot.keys()) - set(prev_snapshot.keys())
        # Find changed genes
        changed_genes = {}
        for g_id in set(prev_snapshot.keys()) & set(current_snapshot.keys()):
            prev_vals = prev_snapshot[g_id]
            curr_vals = current_snapshot[g_id]
            diffs = {}
            for attr, val in curr_vals.items():
                if val != prev_vals[attr]:
                    diffs[attr] = (prev_vals[attr], val)
            if diffs:
                changed_genes[g_id] = diffs

        if deleted_genes or new_genes or changed_genes:
            print(f"--- [Filter Step: {step_name}] Changes ---")
            if new_genes:
                print(f"  New genes: {', '.join(sorted(new_genes))}")
            if deleted_genes:
                print(f"  Removed genes: {', '.join(sorted(deleted_genes))}")
            if changed_genes:
                for g_id, diffs in sorted(changed_genes.items()):
                    diff_str = ", ".join(f"{attr}: {prev} -> {curr}" for attr, (prev, curr) in diffs.items())
                    print(f"  Gene {g_id}: {diff_str}")
            print("------------------------------------------")
        return current_snapshot

    def filter(self, source_priority:list, quiet:bool=False):
        if not quiet:
            print("=== Initial BLAST values ===")
            for chrom, genes in self._annot.chrs.items():
                for g_id, g in sorted(genes.items()):
                    gene_blasts = []
                    for b in g.quality.blast_hits:
                        gene_blasts.append(f"{b.source}:{b.score}:{b.evalue}")
                    
                    prot_blasts = []
                    for t in g.transcripts.values():
                        for c in t.CDSs.values():
                            if c.protein:
                                for b in c.protein.blast_hits:
                                    prot_blasts.append(f"{b.source}:{b.score}:{b.evalue}")
                    
                    blast_info = []
                    if gene_blasts:
                        blast_info.append(f"gene_blasts: {sorted(gene_blasts)}")
                    if prot_blasts:
                        blast_info.append(f"protein_blasts: {sorted(prot_blasts)}")
                    if blast_info:
                        print(f"  Gene {g_id}: {'; '.join(blast_info)}")
            print("============================")

        snapshot = self._take_snapshot()

        self._annot.remove_duplicate_transcripts(quiet=quiet)
        snapshot = self._print_changes("remove_duplicate_transcripts", snapshot, quiet=quiet)

        self._annot.make_alternative_transcripts_into_genes(quiet=quiet)
        snapshot = self._print_changes("make_alternative_transcripts_into_genes", snapshot, quiet=quiet)

        self._annot.overlaps.detect(quiet=quiet)
        snapshot = self._print_changes("overlaps.detect", snapshot, quiet=quiet)

        self._annot.stats.calculate_transcript_masking()
        snapshot = self._print_changes("calculate_transcript_masking", snapshot, quiet=quiet)

        self.mark_noisy_genes(quiet=quiet)
        snapshot = self._print_changes("mark_noisy_genes", snapshot, quiet=quiet)

        self._annot.remove_genes(quiet=quiet)
        snapshot = self._print_changes("remove_genes", snapshot, quiet=quiet)

        self.mark_transcriptomic_supported_genes(quiet=quiet)
        snapshot = self._print_changes("mark_transcriptomic_supported_genes", snapshot, quiet=quiet)

        self.mark_abinitio_supported_genes(quiet=quiet)
        snapshot = self._print_changes("mark_abinitio_supported_genes", snapshot, quiet=quiet)

        self.add_reliable_CDS_evidence_score(quiet=quiet)
        snapshot = self._print_changes("add_reliable_CDS_evidence_score", snapshot, quiet=quiet)

        self.find_best_gene_model(source_priority, quiet=quiet)
        snapshot = self._print_changes("find_best_gene_model (with reliables)", snapshot, quiet=quiet)

        self.mark_overlap_with_reliable_genes(quiet=quiet)
        snapshot = self._print_changes("mark_overlap_with_reliable_genes", snapshot, quiet=quiet)

        self.find_best_gene_model(source_priority, just_with_reliables=False, quiet=quiet)
        snapshot = self._print_changes("find_best_gene_model (without reliables)", snapshot, quiet=quiet)

        self.add_better_ab_initio_models_as_alternative_transcripts(source_priority, reliable_sources=["augustus", "Liftoff", "genemark"], quiet=quiet)
        snapshot = self._print_changes("add_better_ab_initio_models_as_alternative_transcripts", snapshot, quiet=quiet)

        self._annot.overlaps.detect(quiet=quiet)
        snapshot = self._print_changes("overlaps.detect", snapshot, quiet=quiet)

        self.rescue_longer_same_frame_CDS(quiet=quiet)
        snapshot = self._print_changes("rescue_longer_same_frame_CDS", snapshot, quiet=quiet)

        self._annot.overlaps.detect(quiet=quiet)
        snapshot = self._print_changes("overlaps.detect", snapshot, quiet=quiet)

        self.remove_CDS_overlaps(source_priority)
        snapshot = self._print_changes("remove_CDS_overlaps", snapshot, quiet=quiet)

        self.mark_intron_nesting()
        snapshot = self._print_changes("mark_intron_nesting", snapshot, quiet=quiet)

        self.remove_fully_intron_nested_genes()
        snapshot = self._print_changes("remove_fully_intron_nested_genes", snapshot, quiet=quiet)

        self._annot.make_alternative_transcripts_into_genes(quiet=quiet)
        snapshot = self._print_changes("make_alternative_transcripts_into_genes", snapshot, quiet=quiet)

        self._annot.overlaps.detect(quiet=quiet)
        snapshot = self._print_changes("overlaps.detect", snapshot, quiet=quiet)

        self.mark_overlap_with_other_selected_CDSs(quiet=quiet)
        snapshot = self._print_changes("mark_overlap_with_other_selected_CDSs", snapshot, quiet=quiet)

        self.mark_overlap_with_other_selected_exons(quiet=quiet)
        snapshot = self._print_changes("mark_overlap_with_other_selected_exons", snapshot, quiet=quiet)

        self.select_best_possible_non_overlapping_UTR(quiet=quiet)
        snapshot = self._print_changes("select_best_possible_non_overlapping_UTR", snapshot, quiet=quiet)

        self._annot.overlaps.detect(quiet=quiet)
        snapshot = self._print_changes("overlaps.detect", snapshot, quiet=quiet)

        self.mark_overlap_with_other_selected_CDSs(quiet=quiet)
        snapshot = self._print_changes("mark_overlap_with_other_selected_CDSs", snapshot, quiet=quiet)

        self.mark_overlap_with_other_selected_exons(quiet=quiet)
        snapshot = self._print_changes("mark_overlap_with_other_selected_exons", snapshot, quiet=quiet)

        self.select_best_possible_non_overlapping_UTR(exon=True, quiet=quiet)

        self._annot.remove_genes(quiet=quiet)
        snapshot = self._print_changes("remove_genes", snapshot, quiet=quiet)

        self._annot.overlaps.detect(quiet=quiet)
        snapshot = self._print_changes("overlaps.detect", snapshot, quiet=quiet)

        self._annot.make_alternative_genes_into_transcripts(quiet=quiet)
        snapshot = self._print_changes("make_alternative_genes_into_transcripts", snapshot, quiet=quiet)
        self._annot.overlaps.detect(quiet=quiet)
        snapshot = self._print_changes("overlaps.detect", snapshot, quiet=quiet)

        self.find_best_gene_model_nested_overlaps(source_priority)
        snapshot = self._print_changes("find_best_gene_model_nested_overlaps", snapshot, quiet=quiet)

        self.find_best_gene_model_exon_num_overlaps(source_priority)
        snapshot = self._print_changes("find_best_gene_model_exon_num_overlaps", snapshot, quiet=quiet)

        self.remove_exon_overlaps(source_priority)
        snapshot = self._print_changes("remove_exon_overlaps", snapshot, quiet=quiet)

        self.remove_UTRs_from_exon_overlaps()
        snapshot = self._print_changes("remove_UTRs_from_exon_overlaps", snapshot, quiet=quiet)

        self._annot.remove_genes(quiet=quiet)
        snapshot = self._print_changes("remove_genes", snapshot, quiet=quiet)

        self._annot.update(rename_features=["gene", "transcript", "CDS", "exon", "UTR"], quiet=quiet)

        snapshot = self._print_changes("update (rename features)", snapshot, quiet=quiet)

        self._annot.overlaps.detect(quiet=quiet)
        snapshot = self._print_changes("overlaps.detect", snapshot, quiet=quiet)

    def remove_alternative(self):
        nodes = self._annot.overlaps.networks[chr][0].nodes()
        print("Nodes in the graph:")
        for node in nodes:
            print(node)

        # Find articulation points (connector nodes)
        connector_nodes = list(nx.articulation_points(self._annot.overlaps.networks[chr][0]))

        # Remove connector nodes from the graph
        for node in connector_nodes:
            self._annot.overlaps.networks[chr][0].remove_node(node)

