from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .hits import BlastHit

import copy

from .feature import Feature
from .utils.genefunctions import get_start_codons, get_alt_start_codons

INITIATOR_METHIONINE_MODES = ("canonical", "all", "none")

class Protein():

    __slots__ = ("id", "ch", "readthrough", "_blast_hits", "seq", "nuc_seq", "start", "end", "_segments", "table", "trimmed_5p", "trimmed_3p", "frameshifts", "start_status")

    _blast_hits: list[BlastHit] | None
    _segments: tuple[tuple[int, int], ...] | None
    start: int
    end: int
    ch: str
    readthrough: str
    seq: str
    nuc_seq: str
    table: int | str | dict[str, str]
    trimmed_5p: int
    trimmed_3p: int
    frameshifts: int
    start_status: str

    def __init__(self, prot_id:str, sequence:str, chrom:str, start:int, end:int, readthrough:str, nuc_seq:str="", segments:tuple[tuple[int, int], ...]|None=None, table:int|str|dict[str, str]=1, trimmed_5p:int=0, trimmed_3p:int=0, frameshifts:int=0, initiator_methionine:str="canonical"):
        """
        trimmed_5p: CDS bases skipped before the first translated codon (nonzero
            initial phase or 5' trimming); the first codon is then not a start codon.
        trimmed_3p: CDS bases left untranslated after the last complete codon.
        frameshifts: segment junctions where bases were dropped to follow annotated
            phases (internal frameshifts corrected during translation).
        initiator_methionine: which start codons are translated as M when they
            initiate the protein (NCBI convention): "canonical" (curated starts of
            the table), "all" (also alternative NCBI initiators) or "none".
        """
        if initiator_methionine not in INITIATOR_METHIONINE_MODES:
            raise ValueError(f"initiator_methionine must be one of {INITIATOR_METHIONINE_MODES}, got '{initiator_methionine}'.")
        self.id = prot_id
        self.table = table
        self.trimmed_5p = trimmed_5p
        self.trimmed_3p = trimmed_3p
        self.frameshifts = frameshifts
        self.ch = chrom
        self.start = start
        self.end = end
        self.readthrough = readthrough
        self._blast_hits = None
        self._segments = segments

        self.seq = sequence
        self.nuc_seq = nuc_seq

        self.start_status = self._grade_start()
        if self.start_status == "canonical" and initiator_methionine != "none" or self.start_status == "alternative" and initiator_methionine == "all":
            self.seq = "M" + self.seq[1:]

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
    def ambiguous_residues(self) -> int:
        """Number of residues translated from ambiguous codons (X)."""
        return self.seq.count("X")

    @property
    def gaps(self) -> bool:
        return "X" in self.seq

    @property
    def nucleotide_surplus(self) -> bool:
        """Whether CDS bases were left untranslated at either end (see trimmed_5p / trimmed_3p)."""
        return bool(self.trimmed_5p or self.trimmed_3p)
    
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
    
    def _grade_start(self) -> str:
        """
        "canonical" if the first codon is a curated start codon of the table, "alternative"
        if it is another NCBI initiator of the table, otherwise "none" (also when 5' bases
        of the CDS were skipped). Without nucleotide sequence, an initial M is canonical.
        """
        if not self.seq or self.trimmed_5p:
            return "none"
        if self.nuc_seq and len(self.nuc_seq) >= 3:
            codon = self.nuc_seq[:3].upper().replace("U", "T")
            if codon in get_start_codons(self.table):
                return "canonical"
            if codon in get_alt_start_codons(self.table):
                return "alternative"
            return "none"
        return "canonical" if self.seq.startswith("M") else "none"

    @property
    def ATG_start(self) -> bool:
        """Whether the protein begins with a canonical or alternative start codon of its table (5' complete)."""
        return self.start_status != "none"

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
        """Whether the protein is partial at the 5' end (does not begin with a start codon of its table)."""
        return not self.ATG_start

    @property
    def partial_3prime(self) -> bool:
        """Whether the protein is partial at the 3' end (does not end with a stop codon)."""
        return not self.end_stop

    @property
    def partial(self) -> bool:
        return self.partial_5prime or self.partial_3prime

    @property
    def truncated(self) -> bool:
        """Whether the protein has premature (internal) stop codons."""
        return self.early_stop

    @property
    def summary_tag(self) -> str:
        summary_tag = []
        if self.partial:
            summary_tag.append("partial")
        if self.truncated:
            summary_tag.append("truncated")
        if self.start_status == "alternative":
            summary_tag.append("alt_start")
        if self.gaps:
            summary_tag.append("ambiguous")
        if self.frameshifts:
            summary_tag.append("frameshift")
        return "_".join(summary_tag)


class Promoter(Feature):
    __slots__ = ('type',)
    def __init__(self, promoter_type, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)
        self.type = promoter_type
