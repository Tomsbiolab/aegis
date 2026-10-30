from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .hits import BlastHit

import copy

from .feature import Feature
from .utils.genefunctions import get_start_codons, get_alt_start_codons, TRANSL_EXCEPT_NAMES

class Protein():

    __slots__ = ("id", "ch", "readthrough", "_blast_hits", "seq", "nuc_seq", "start", "end", "_segments", "table", "trimmed_5p", "trimmed_3p", "frameshifts", "exception", "transl_except", "start_status", "polyA_bases")

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
    frameshifts: tuple[int, ...]
    exception: str|None
    transl_except: tuple[tuple[int, int, str], ...]
    start_status: str

    def __init__(self, prot_id:str, sequence:str, chrom:str, start:int, end:int, readthrough:str, nuc_seq:str="", segments:tuple[tuple[int, int], ...]|None=None, table:int|str|dict[str, str]=1, trimmed_5p:int=0, trimmed_3p:int=0, frameshifts:tuple[int, ...]=(), exception:str|None=None, transl_except:tuple[tuple[int, int, str], ...]=(), polyA_bases:int=0, initiator_methionine:bool=True, alt_start_methionine:bool=True):
        """
        trimmed_5p: CDS bases skipped before the first translated codon (nonzero
            initial phase or 5' trimming); the first codon is then not a start codon.
        trimmed_3p: CDS bases left untranslated after the last complete codon.
        frameshifts: net shift (in bases) of each translated segment junction in transcript
            order: negative when bases are read twice (annotated overlap), positive when bases
            are skipped (annotated gap, or dropped to follow annotated phases).
        exception: GFF exception of the CDS (e.g. "ribosomal slippage"), if any.
        transl_except: residues recoded from the CDS transl_except, as (start, end, residue)
            1-based positions within nuc_seq; they take precedence over the initiator methionine.
        polyA_bases: A's the polyA tail adds to complete the stop codon (transl_except aa:TERM): they end
            nuc_seq but have no genomic position, so nuc_seq is the segments' genomic bases plus these A's.
        initiator_methionine: whether valid start codons at position 1 are translated
            as Methionine (M) [default: True, matching NCBI RefSeq/UniProt convention].
        alt_start_methionine: whether alternative start codons (e.g. CTG in table 1/11)
            are also translated as M when initiator_methionine is True [default: True].
        """
        self.id = prot_id
        self.table = table
        self.trimmed_5p = trimmed_5p
        self.trimmed_3p = trimmed_3p
        self.frameshifts = frameshifts
        self.exception = exception
        self.transl_except = transl_except
        self.polyA_bases = polyA_bases
        self.ch = chrom
        self.start = start
        self.end = end
        self.readthrough = readthrough
        self._blast_hits = None
        self._segments = segments

        self.seq = sequence
        self.nuc_seq = nuc_seq

        self.start_status = self._grade_start()
        if initiator_methionine and not any(start == 1 for start, _, _ in transl_except):
            if self.start_status == "canonical" or (self.start_status == "alternative" and alt_start_methionine):
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

    def __len__(self) -> int:
        return len(self.seq)

    @property
    def size(self) -> int:
        """Protein size in amino acids, excluding the terminal stop (internal stops are residues, see internal_stops)."""
        return len(self.seq) - 1 if self.end_stop else len(self.seq)

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
    def internal_stops(self) -> int:
        """Number of stop codons before the last residue."""
        return self.seq[:-1].count("*")

    @property
    def tags(self) -> dict[str, str]:
        """
        Protein tags in output order, the single source for detailed FASTA headers and summaries.
        Quality tags only appear when they apply; readthrough and table are always present.
        """
        tags = {}
        if self.partial:
            tags["partial"] = ",".join(end for end, is_partial in (("5'", self.partial_5prime), ("3'", self.partial_3prime)) if is_partial)
        if self.start_status == "alternative":
            tags["alt_start"] = self.nuc_seq[:3].upper()
        if self.internal_stops:
            tags["internal_stops"] = str(self.internal_stops)
        if self.ambiguous_residues:
            tags["ambiguous_residues"] = str(self.ambiguous_residues)
        if self.frameshifts:
            tags["frameshifts"] = ",".join(f"{shift:+d}" for shift in self.frameshifts)
        if self.exception:
            tags["exception"] = self.exception
        if self.transl_except:
            tags["transl_except"] = ",".join(f"(pos:{start}..{end},aa:{TRANSL_EXCEPT_NAMES[residue]})" for start, end, residue in self.transl_except)
        if self.polyA_bases:
            tags["polyA_bases"] = str(self.polyA_bases)
        tags["readthrough"] = self.readthrough
        tags["table"] = "custom" if isinstance(self.table, dict) else str(self.table)
        return tags



class Promoter(Feature):
    __slots__ = ('type',)
    def __init__(self, promoter_type, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)
        self.type = promoter_type
