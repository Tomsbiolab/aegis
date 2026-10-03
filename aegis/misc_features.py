from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .hits import BlastHit

import copy

from .feature import Feature

class Protein():

    __slots__ = ("id", "ch", "readthrough", "_blast_hits", "nucleotide_surplus", "seq", "nuc_seq", "partial", "truncated", "start", "end", "_segments")

    _blast_hits: list[BlastHit] | None
    _segments: tuple[tuple[int, int], ...] | None
    start: int
    end: int
    ch: str
    readthrough: str
    partial: bool
    truncated: bool
    seq: str
    nuc_seq: str
    nucleotide_surplus: bool

    def __init__(self, prot_id:str, sequence:str, chrom:str, start:int, end:int, nucleotide_surplus:bool, readthrough:str, nuc_seq:str="", segments:tuple[tuple[int, int], ...]|None=None):
        self.id = prot_id
        self.ch = chrom
        self.start = start
        self.end = end
        self.readthrough = readthrough
        self._blast_hits = None
        self._segments = segments

        self.seq = sequence
        self.nuc_seq = nuc_seq
        self.nucleotide_surplus = nucleotide_surplus

        if self.ATG_start == False or self.end_stop == False or self.nucleotide_surplus or self.gaps:
            self.partial = True
        else:
            self.partial = False

        if self.early_stop:
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
        if "X" in self.seq:
            return True
        return False
    
    def __len__(self) -> int:
        return len(self.seq)

    @property
    def size(self) -> int:
        """Protein size in amino acids (excluding terminal stop codon)."""
        return len(self.seq.rstrip("*"))

    @property
    def genomic_span(self) -> int:
        """Genomic chromosome span (end - start + 1), including any introns."""
        return (self.end - self.start) + 1
    
    @property
    def ATG_start(self) -> bool:
        if not self.seq:
            return False
        if self.seq.startswith("M"):
            return True
        if self.nuc_seq and len(self.nuc_seq) >= 3:
            st = self.nuc_seq[:3].upper().replace("U", "T")
            if st in ("ATG", "GTG", "TTG", "ATA", "ATT", "ATC"):
                return True
        return False

    @property
    def end_stop(self) -> bool:
        return self.seq.endswith("*") if self.seq else False

    @property
    def early_stop(self) -> bool:
        return "*" in self.seq[:-1]

    @property
    def ATG_late(self) -> bool:
        return "M" in self.seq[1:]

    @property
    def segments(self) -> tuple[tuple[int, int], ...]:
        """Genomic coordinate intervals for the translated protein codons."""
        if self._segments is not None:
            return self._segments
        return ((self.start, self.end),)

    @property
    def partial_5prime(self) -> bool:
        """Whether the protein is partial at the 5' end (lacks ATG start)."""
        return not self.ATG_start

    @property
    def partial_3prime(self) -> bool:
        """Whether the protein is partial at the 3' end (lacks stop codon or has nucleotide surplus)."""
        return not self.end_stop or self.nucleotide_surplus

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
