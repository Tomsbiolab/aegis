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

    def get_blast_hits_key(self, source_priority:list) -> tuple:
        best = {}
        for b in self.blast_hits:
            if b.source not in best:
                best[b.source] = {"evalue": b.evalue, "score": b.score}
            else:
                if b.evalue < best[b.source]["evalue"]:
                    best[b.source]["evalue"] = b.evalue
                if b.score > best[b.source]["score"]:
                    best[b.source]["score"] = b.score
                    
        key = []
        for s in source_priority:
            hit = best.get(s, {"evalue": float('inf'), "score": float('-inf')})
            key.append(hit["evalue"])
            key.append(-hit["score"])
            
        return tuple(key)

    def compare_blast_hits(self, other:Protein, source_priority:list) -> bool:
        if not source_priority:
            return True
        return self.get_blast_hits_key(source_priority) <= other.get_blast_hits_key(source_priority)

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
