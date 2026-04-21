from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .hits import BlastHit

import copy

from .feature import Feature

class Protein():

    __slots__ = ("id", "ch", "readthrough", "_blast_hits", "nucleotide_surplus", "seq", "partial", "truncated", "start", "end")

    _blast_hits: list[BlastHit] | None
    start: int
    end: int
    ch: str
    readthrough: str
    partial: bool
    truncated: bool
    seq: str
    nucleotide_surplus: bool

    def __init__(self, prot_id:str, sequence:str, chrom:str, start:int, end:int, nucleotide_surplus:bool, readthrough:str):
        self.id = prot_id
        self.ch = chrom
        self.start = start
        self.end = end
        self.readthrough = readthrough
        self._blast_hits = None

        self.seq = sequence
        self.nucleotide_surplus = nucleotide_surplus

        if self.ATG_start == False or self.end_stop == False or self.nucleotide_surplus or self.gaps:
            self.partial = True
        else:
            self.partial = False

        if self.early_stop and self.end_stop:
            self.truncated = True
        else:
            self.truncated = False

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

    @property
    def gaps(self):
        if "-" in self.seq:
            return True
        return False
    
    @property
    def size(self) -> int:
        return (self.end - self.start) + 1
    
    @property
    def ATG_start(self) -> bool:
        return self.seq[0] == "M"

    @property
    def end_stop(self) -> bool:
        return self.seq[-1] == "*"

    @property
    def early_stop(self) -> bool:
        return "*" in self.seq[:-1]

    @property
    def ATG_late(self) -> bool:
        return "M" in self.seq[1:]

    @property
    def summary_tag(self) -> str:
        summary_tag = []
        if self.partial:
            summary_tag.append("partial")
        if self.truncated:
            summary_tag.append("truncated")
        return "_".join(summary_tag)


class Promoter(Feature):
    __slots__ = ('type',)
    def __init__(self, promoter_type, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)
        self.type = promoter_type
