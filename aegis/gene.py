from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .genome import Genome
    from .transcript import Transcript
    from .hits import OverlapHit

from .feature import Feature
from .subfeatures import Exon

class Gene(Feature):

    __slots__ = (
        'pseudogene', 'transposable', 'transcripts', 'previous_gene',
        'next_gene', 'synteny_order', 'old_previous_gene', 'old_next_gene',
        'old_synteny_order', 'conserved_synteny', 'coding', 'noncoding',
        'overlaps', 'remove', 'rescue', 'reliable', 'reliable_score',
        'overlap_reliable', 'unrescuable', 'overlap_with_selected_CDS',
        'overlap_with_selected_exon', 'alternative_transcript_rescue',
        'intron_nested', 'intron_nested_fully_contained', 'intron_nested_single',
        'UTR_intron_nested', 'transcriptomic_evidence', 'abinitio_evidence',
        'base_id', 'original_base_id'
    )

    transcripts:dict[str, Transcript]
    synteny_order:int|None
    old_synteny_order:int|None
    previous_gene:str|None|bool
    next_gene:str|None|bool
    old_previous_gene:str|None|bool
    old_next_gene:str|None|bool
    conserved_synteny:bool|None
    alternative_transcript_rescue:list
    overlaps:dict[str, list[OverlapHit]]
    
    def __init__(self, pseudogene:bool, transposable:bool, feature_id:str, 
                 ch:str, source:str, feature:str, strand:str,
                 start:int, end:int, score:str, phase:str, attributes:str|list|dict):
        super().__init__(feature_id, ch, source, feature, strand, start, end,
                         score, phase, attributes)
        self.pseudogene = pseudogene
        self.transposable = transposable
        # transcripts will be added as {"transcript_id" : transcript_object}
        self.transcripts = {}
        # order within the chromosome
        self.previous_gene = None
        self.next_gene = None
        self.synteny_order = None
        self.old_previous_gene = None
        self.old_next_gene = None
        self.old_synteny_order = None
        self.conserved_synteny = None
        # A gene can have both coding and noncoding transcripts
        self.coding = False
        self.noncoding = False

        self.overlaps = {"self" : [], "other" : []}

        self.remove = False
        self.rescue = False
        self.reliable = False
        self.reliable_score = 0
        self.overlap_reliable = False
        self.unrescuable = False

        self.overlap_with_selected_CDS = False
        self.overlap_with_selected_exon = False
        self.alternative_transcript_rescue = []

        self.intron_nested = False
        self.intron_nested_fully_contained = False
        self.intron_nested_single = False
        self.UTR_intron_nested = False

        self.transcriptomic_evidence = False
        self.abinitio_evidence = False


        self.obtain_base_id(original=True)


    def update(self):
        self.update_size()
        self.sort_transcripts()
        self.coding = False
        self.noncoding = False
        for t in self.transcripts.values():
            if t.coding:
                self.coding = True
            else:
                self.noncoding = True
            if self.strand == ".":
                if t.strand != ".":
                    self.strand = t.strand

        if self.coding:
            best_id = None
            best_cds_size = -1
            best_exon_size = -1
            for t in self.transcripts.values():
                t.main = False
                cds_size = 0
                if t.CDSs != {}:
                    t.determine_main_CDS()
                    for c in t.CDSs.values():
                        if c.main:
                            cds_size = c.size
                exon_size = t.size
                # Select by largest CDS, break ties by largest exon span
                if (cds_size > best_cds_size or (cds_size == best_cds_size and exon_size > best_exon_size)):
                    best_cds_size = cds_size
                    best_exon_size = exon_size
                    best_id = t.id
            if best_id:
                self.transcripts[best_id].main = True

        elif self.noncoding:
            best_id = None
            best_exon_size = -1
            for t in self.transcripts.values():
                t.main = False
                if t.size > best_exon_size:
                    best_exon_size = t.size
                    best_id = t.id
            if best_id:
                self.transcripts[best_id].main = True
        else:
            print(f"Error: gene {self.id} has no transcripts annotated")

        self.homogenise_exon_scores()

    def obtain_base_id(self, original:bool=False):

        if self.id.endswith("_gene"):
            self.base_id = self.id[:-5]
        elif self.id.endswith("gene"):
            self.base_id = self.id[:-4]
        elif self.id.startswith("gene:"):
            self.base_id = self.id[5:]
        elif self.id.startswith("gene-"):
            self.base_id = self.id[5:]
        elif self.id.startswith("gene"):
            self.base_id = self.id[4:]
        else:
            self.base_id = self.id

        if original:
            self.original_base_id = self.base_id

    def rename(self, count:int, sep:str="_", digits:int=5, prefix:str="", suffix:str="", base_id_as_id:bool=False, remove_point_suffix:bool=False):

        if remove_point_suffix:
            if "." in self.id:
                self.id = self.id.split(".")[0]
                if self.original_id != self.id:
                    self.renamed = True
                    self.obtain_base_id()

        if base_id_as_id:
            if self.base_id != self.id:
                self.id = self.base_id
                if self.original_id != self.id:
                    self.renamed = True

        if prefix:

            if suffix:
                self.id = f"{prefix}{self.ch}g{count:0{digits}d}{sep}{suffix}"
            else:
                self.id = f"{prefix}{self.ch}g{count:0{digits}d}"
    
            if self.original_id != self.id:
                self.renamed = True
                self.obtain_base_id()
            
        if self.renamed:
            self.update_numbering()

    def sort_transcripts(self):
        sorted_transcripts = sorted(self.transcripts.values())
        self.transcripts = {t.id: t for t in sorted_transcripts}

    def homogenise_exon_scores(self):
        all_exons = {}
        for t in self.transcripts.values():
            for e in t.exons:
                tag = (e.start, e.end, e.strand)
                if tag not in all_exons:
                    all_exons[tag] = e.score
                elif e.score != ".":
                    if all_exons[tag] == ".":
                        all_exons[tag] = e.score
                    elif float(e.score) > float(all_exons[tag]):
                        all_exons[tag] = e.score

        for t in self.transcripts.values():
            for e in t.exons:
                e.score = all_exons[(e.start, e.end, e.strand)]

    def clear_UTRs(self):
        for t in self.transcripts.values():
            t.clear_UTRs()
        self.update()

    def combine_transcripts(self, genome:Genome, low_memory:bool=True, respect_non_coding:bool=False, quiet:bool=False):
        """
        Useful for RNA-Seq read counting for transcript variants as "one" gene.
        """
        temp_fts = []
        for t in self.transcripts.values():
            for e in t.exons:
                temp_fts.append(e.copy())
        self.transcripts = {}
        if temp_fts != []:
            still_overlapping_fts = True
            while still_overlapping_fts:
                features_to_remove = []
                features_to_add = []
                for i, tempft1 in enumerate(temp_fts):
                    for j, tempft2 in enumerate(temp_fts):
                        if i != j:
                            interval1 = tempft1.size
                            interval2 = tempft2.size
                            small = min(tempft1.start, tempft1.end, tempft2.start, tempft2.end)
                            large = max(tempft1.start, tempft1.end, tempft2.start, tempft2.end)
                            overlap_bp = ((interval1 + interval2) - ((large - small) + 1))
                            # >= is crucial to combine contiguous features
                            if overlap_bp >= 0:
                                temp = Exon("combined", self.ch, self.source, "exon", self.strand, small, large, self.score, ".", "")
                                add = True
                                # this is to avoid adding a same overlap twice
                                for f in features_to_add:
                                    if temp.equal_coordinates(f):
                                        add = False
                                        break
                                if add:
                                    features_to_add.append(temp)
                                if tempft1 not in features_to_remove:
                                    features_to_remove.append(tempft1)
                                if tempft2 not in features_to_remove:
                                    features_to_remove.append(tempft2)
                for sub_to_add in features_to_add:
                    temp_fts.append(sub_to_add)
                for sub_to_rem in features_to_remove:
                    temp_fts.remove(sub_to_rem)
                if features_to_remove == [] and features_to_add == []:
                    still_overlapping_fts = False
            temp_fts.sort()

            temp_coding_feature = "mRNA"
            if not self.coding:
                temp_coding_feature = "lncRNA"
            self.transcripts[f"{self.id}_t001"] = Transcript(f"{self.id}_t001", self.ch, self.source, temp_coding_feature, self.strand, temp_fts[0].start, temp_fts[-1].end, self.score, ".", f"ID={self.id}_t001;Parent={self.id}")
            self.transcripts[f"{self.id}_t001"].exons = temp_fts.copy()
            counter = 0
            for e in self.transcripts[f"{self.id}_t001"].exons:
                counter += 1
                e.feature = "exon"
                e.id = f"{self.id}_generated_exon_{counter}"
                e.attributes = [f"ID={e.id}", f"Parent={self.id}_t001"]
                e.parents = [f"{self.id}_t001"]

        for t in self.transcripts.values():
            t.update(consider_polycistronic=False, consider_read_utrs=False, quiet=quiet)
            if respect_non_coding:
                if not self.coding:
                    continue
            t.generate_sequence(genome, low_memory)
            t.generate_best_protein(genome)
            t.generate_CDSs_based_on_ORF(low_memory)
            for c in t.CDSs.values():
                c.generate_sequence(genome, low_memory)
            if t.coding_ratio < 0.80:
                t.generate_sequence(genome, low_memory)
                t.generate_best_protein(genome, must_have_stop=False)
                t.generate_CDSs_based_on_ORF(low_memory)
                for c in t.CDSs.values():
                    c.generate_sequence(genome, low_memory)
            t.update(consider_polycistronic=False, consider_read_utrs=False, quiet=quiet)

    def longer_CDS(self, other:Gene):
        for t1 in self.transcripts.values():
            if t1.main:
                for c1 in t1.CDSs.values():
                    if c1.main:
                        for t2 in other.transcripts.values():
                            if t2.main:
                                for c2 in t2.CDSs.values():
                                    if c2.main:
                                        if c1.longer(c2):
                                            return True
                                        else:
                                            return False
                                        
    def compare_protein_blast_hits(self, other:Gene, source_priority:list):
        """
        Method required to deal with the fact that a gene may have several blast hits due to several proteins...
        """
        self_proteins = []
        for t in self.transcripts.values():
            for c in t.CDSs.values():
                self_proteins.append(c.protein)

        other_proteins = []
        for t in other.transcripts.values():
            for c in t.CDSs.values():
                other_proteins.append(c.protein)

        best_self_protein = None
        if self_proteins == []:
            print(f"Warning: {self.id} gene has no associated proteins.")
        elif len(self_proteins) == 1:
            best_self_protein = self_proteins[0]
        else:
            best_self_protein = self_proteins[0]
            for n, p in enumerate(self_proteins):
                if n > 0:
                    current_best = best_self_protein.compare_blast_hits(p, source_priority)
                    if not current_best:
                        best_self_protein = p

        best_other_protein = None
        if other_proteins == []:
            print(f"Warning: {other.id} gene has no associated proteins.")
        elif len(other_proteins) == 1:
            best_other_protein = other_proteins[0]
        else:
            best_other_protein = other_proteins[0]
            for n, p in enumerate(other_proteins):
                if n > 0:
                    current_best = best_other_protein.compare_blast_hits(p, source_priority)
                    if not current_best:
                        best_other_protein = p
        
        if best_self_protein != None and best_other_protein != None:
            query_best = best_self_protein.compare_blast_hits(best_other_protein, source_priority)

            return query_best
        
        elif best_other_protein == None and best_other_protein == None:
            print(f"Warning: {self.id} or {other.id} genes have no associated proteins.")
            
            return None
        
        elif best_other_protein == None:

            return True

        elif best_self_protein == None:
            
            return False

    def get_main_CDS_range(self):
        for t in self.transcripts.values():
            if t.main:
                for c in t.CDSs.values():
                    if c.main:
                        return c.start, c.end
                break

        return None

    def __str__(self):
        if self.symbols != []:
            return f"{self.id}: {self.symbols}"
        elif self.names != []:
            return f"{self.id}: {self.names}"
        else:
            return f"{self.id}"
