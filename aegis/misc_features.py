from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .hits import BlastHit

import copy

from .utils.genefunctions import flexible_translate
from .feature import Feature

class Protein():

    __slots__ = ("id", "ch", "summary_tag", "readthrough", "_blast_hits", "start", "end_stop", "early_stop", "nucleotide_surplus", "gaps", "seq", "coding_start", "coding_end", "partial", "truncated", "size")

    _blast_hits: list[BlastHit] | None
    size: int

    def __init__(self, prot_id:str, nucleotides:str, chrom:str, readthrough:str="both"):
        self.id = prot_id
        self.ch = chrom
        self.summary_tag = ""
        self.readthrough = readthrough
        self._blast_hits = None
        # readthrough can be start, end, both or none 
        self.start, self.end_stop, self.early_stop, self.nucleotide_surplus, self.gaps, self.seq, self.coding_start, self.coding_end = flexible_translate(nucleotides, readthrough=readthrough)
        if self.start == "late" or self.start == "absent" or self.end_stop == False or self.nucleotide_surplus or self.gaps:
            self.partial = True
            self.summary_tag += "partial"
        else:
            self.partial = False
        if self.early_stop and self.end_stop:
            self.truncated = True
            if self.summary_tag == "":
                self.summary_tag += "truncated"
            else:
                self.summary_tag += "_truncated"
        else:
            self.truncated = False
        self.size = len(self.seq)

    def copy(self):
        return copy.deepcopy(self)

    def compare_blast_hits(self, other:Protein, source_priority:list):
        compared = False
        query_best = True
        while not compared:
            for s in source_priority:
                query_evalue = float(2)
                query_bitscore = float(-1)
                target_evalue = float(2)
                target_bitscore = float(-1)
                
                for b in self.blast_hits:
                    if b.source == s:
                        if b.evalue < query_evalue:
                            query_evalue = b.evalue
                        if b.score > query_bitscore:
                            query_bitscore = b.score

                for b in other.blast_hits:
                    if b.source == s:
                        if b.evalue < target_evalue:
                            target_evalue = b.evalue
                        if b.score > target_bitscore:
                            target_bitscore = b.score

                if query_evalue < target_evalue:
                    compared = True
                    break
                elif query_evalue > target_evalue:
                    query_best = False
                    compared = True
                    break
                elif query_bitscore > target_bitscore:
                    compared = True
                    break
                elif query_bitscore < target_bitscore:
                    query_best = False
                    compared = True
                    break
                elif s == source_priority[-1]:
                    compared = True

        return query_best

    @property
    def blast_hits(self):
        if self._blast_hits is None:
            self._blast_hits = []
        return self._blast_hits

class Promoter(Feature):
    __slots__ = ('type',)
    def __init__(self, promoter_type, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)
        self.type = promoter_type
