from __future__ import annotations

from typing import TYPE_CHECKING, Literal

from .feature import Feature
from .misc_features import Protein
from .utils.genefunctions import trim_surplus, map_relative_to_genomic, translate, find_ORFs, choose_orf, reverse_complement

class CDS(Feature):

    __slots__ = ('main', 'CDS_segments', 'full_UTR_exons', 'protein', 'UTRs')

    protein: Protein|None
    CDS_segments: list[Feature]
    UTRs: list[UTR]
    full_UTR_exons: int
    main:bool

    def __init__(self, CDS_segments:list, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)    
        self.main = False
        self.CDS_segments = CDS_segments
        self.full_UTR_exons = 0
        self.protein = None
        self.update()

    def update(self):
        self.update_phase()
        self.update_frame()

    @property
    def size(self):
        size = 0
        for segment in self.CDS_segments:
            size += segment.size
        return size

    def update_phase(self):
        leftover = 0
        if self.strand == "+":
            for cs in self.CDS_segments:
                if leftover == 0:
                    cs.phase = 0
                else:
                    cs.phase = 3 - leftover
                leftover = (cs.size - cs.phase) % 3
        elif self.strand == "-":
            for cs in reversed(self.CDS_segments):
                if leftover == 0:
                    cs.phase = 0
                else:
                    cs.phase = 3 - leftover
                leftover = (cs.size - cs.phase) % 3    

    def update_frame(self):
        if self.strand == "+":
            for cs in self.CDS_segments:
                frame = (cs.start + cs.phase) % 3 #type: ignore
                if frame == 0:
                    frame = 3
                cs.frame = frame

        if self.strand == "-":
            for cs in reversed(self.CDS_segments):
                frame = (cs.start + cs.phase) % 3 #type: ignore
                if frame == 0:
                    frame = 3
                frame = 7 - frame
                cs.frame = frame

    def rename(self, base_id:str, base_gene_id:str, count:int, sep:str="_", digits:int=3, keep_numbering:bool=False, keep_existing_ids_if_derived_from_base_id:bool=False, cds_segment_ids:bool=False):

        rename = False
        rename_cs = False

        if keep_existing_ids_if_derived_from_base_id:
            if base_gene_id not in self.id:
                rename = True
            for cs in self.CDS_segments:
                if base_gene_id not in cs.id:
                    rename_cs = True
        else:
            rename = True
            rename_cs = True

        if rename:
            self.renamed = True

            if keep_numbering and self.id_number != None:
                self.id = f"{base_id}{sep}CDS{self.id_number:0{digits}d}"
            else:
                if self.main:
                    self.id = f"{base_id}{sep}CDS{1:0{digits}d}"
                else:
                    self.id = f"{base_id}{sep}CDS{count:0{digits}d}"

            if self.original_id != self.id:
                self.renamed = True
                self.update_numbering()
                if self.protein is not None:
                    self.protein.id = f"{self.id}.prot"

        cs_count = 0
        for cs in self.CDS_segments:
            cs_count += 1

            if rename_cs:

                if keep_numbering and cs.id_number != None:
                    cs.id = f"{base_id}{sep}CDS{cs.id_number}"
                elif cds_segment_ids:
                    cs.id = f"{base_id}{sep}CDS{self.id_number}{sep}{cs_count}"
                else:
                    cs.id = self.id

                if cs.original_id != cs.id:
                    self.renamed = True
                    cs.update_numbering()

    def clear_UTRs(self):
        self.full_UTR_exons = 0
        self.UTRs = []

    @property
    def seq(self) -> str:
        if not self._ACTIVE_GENOME:
            raise ValueError("No genome loaded and you are trying to access the sequence. Load your genome together with your annotation.")
        else:
            cds_seq = ""

            if self.strand == "-":
                for cs in reversed(self.CDS_segments):
                    cds_seq += cs.seq
            else:
                for cs in self.CDS_segments:
                    cds_seq += cs.seq

            return cds_seq

    @property
    def hard_seq(self) -> str:
        if not self._ACTIVE_HARD_GENOME:
            raise ValueError("No hard masked genome loaded and you are trying to access the hard masked sequence. Load your hard masked genome together with your annotation.")
        else:
            cds_seq = ""

            if self.strand == "-":
                for cs in reversed(self.CDS_segments):
                    cds_seq += cs.hard_seq
            else:
                for cs in self.CDS_segments:
                    cds_seq += cs.hard_seq

            return cds_seq

    @property
    def seqs(self) -> list[str]:
        if not self._ACTIVE_GENOME:
            raise ValueError("No genome loaded and you are trying to access the sequence. Load your genome together with your annotation.")
        else:
            cds_seqs = ["", ""]
            for cs in self.CDS_segments:
                cds_seqs[0] += cs.seq

            for cs in reversed(self.CDS_segments):
                cds_seqs[1] += reverse_complement(cs.seq)

            return cds_seqs

    @property
    def hard_seqs(self) -> list[str]:
        if not self._ACTIVE_HARD_GENOME:
            raise ValueError("No hard masked genome loaded and you are trying to access the hard masked sequence. Load your hard masked genome together with your annotation.")
        else:
            cds_seqs = ["", ""]
            for cs in self.CDS_segments:
                cds_seqs[0] += cs.hard_seq

            for cs in reversed(self.CDS_segments):
                cds_seqs[1] += reverse_complement(cs.hard_seq)

            return cds_seqs

    @property
    def five_prime_UTR_seq(self) -> str:
        five_prime_UTR_seq = ""
        if self.strand == "+":
            for u in self.UTRs:
                if u.prime == "5'":
                    five_prime_UTR_seq += u.seq # type: ignore
        elif self.strand == "-":
            for u in reversed(self.UTRs):
                if u.prime == "5'":
                    five_prime_UTR_seq += u.seq # type: ignore
        return five_prime_UTR_seq
    
    @property
    def three_prime_UTR_seq(self) -> str:
        three_prime_UTR_seq = ""
        if self.strand == "+":
            for u in self.UTRs:
                if u.prime == "3'":
                    three_prime_UTR_seq += u.seq # type: ignore
        elif self.strand == "-":
            for u in reversed(self.UTRs):
                if u.prime == "3'":
                    three_prime_UTR_seq += u.seq # type: ignore
        return three_prime_UTR_seq

    def generate_protein(self, mode: Literal["start", "end", "orf", "orf_or_end"] = "end", max_nucleotide_trim: int | None = None, tolerated_stops: int = 0, orf_choice_mode: Literal["longest", "earliest"]="longest", must_have_stop: bool = False, enforce_start_codon: bool = True, min_codon_len: int = 2, start_codons: tuple[str, ...] = ("ATG",), stop_codons: tuple[str, ...] = ("TAA", "TAG", "TGA"), correct_CDS:bool=False, quiet:bool=True):

        if self.strand == ".":
            seq_fw, seq_rv = self.seqs
            fw_orf = choose_orf(find_ORFs(seq_fw, min_codon_len=min_codon_len, enforce_start_codon=enforce_start_codon, must_have_stop=must_have_stop, tolerated_stops=tolerated_stops, start_codons=start_codons, stop_codons=stop_codons), mode=orf_choice_mode)
            rv_orf = choose_orf(find_ORFs(seq_rv, min_codon_len=min_codon_len, enforce_start_codon=enforce_start_codon, must_have_stop=must_have_stop, tolerated_stops=tolerated_stops, start_codons=start_codons, stop_codons=stop_codons), mode=orf_choice_mode)

            if len(fw_orf[0]) >= len(rv_orf[0]):
                self.strand = "+"
                for cs in self.CDS_segments:
                    cs.strand = "+"
            else:
                self.strand = "-"
                for cs in self.CDS_segments:
                    cs.strand = "-"

        coding_seq, nucleotide_surplus, relative_coding_start, relative_coding_end = trim_surplus(self.seq, mode=mode, max_nucleotide_trim=max_nucleotide_trim, orf_choice_mode=orf_choice_mode, must_have_stop=must_have_stop, tolerated_stops=tolerated_stops, enforce_start_codon=enforce_start_codon, start_codons=start_codons, stop_codons=stop_codons, min_codon_len=min_codon_len)

        if relative_coding_end != 0:

            protein_seq = translate(coding_seq)

            corrected_segments = map_relative_to_genomic(segments=self.CDS_segments, rel_start=relative_coding_start, rel_end=relative_coding_end, strand=self.strand)

            protein_start = corrected_segments[0][0]
            protein_end = corrected_segments[-1][1]

            if correct_CDS:

                new_CDS_segments = []
                if self.parents:
                    new_parents = self.parents[:]
                else:
                    new_parents = []

                for start, end in corrected_segments:
                    new_CDS_segments.append(Feature(feature_id=self.id, ch=self.ch, start=start, end=end, strand=self.strand, parents=new_parents, source=self.source, score=self.score, feature=self.feature))

                self.CDS_segments = new_CDS_segments

                self.start = protein_start
                self.end = protein_end

            self.protein = Protein(prot_id=f"{self.id}.prot", sequence=protein_seq, chrom=self.ch, start=protein_start, end=protein_end, nucleotide_surplus=nucleotide_surplus, readthrough=mode)

            if not quiet and nucleotide_surplus:
                print(f"{self.id} has a nucleotide surplus when translating to protein, the annotated CDS might be incorrect.")
        elif not quiet:
            print(f"{self.id} CDS could not be translated to a protein with mode={mode}")

    def clear_protein(self):
        self.protein = None

    def equal_segments(self, other:CDS):
        self.CDS_segments.sort()
        other.CDS_segments.sort()
        same = True
        if len(self.CDS_segments) == len(other.CDS_segments):
            for n, segment in enumerate(self.CDS_segments):
                if not segment.equal_sequence(other.CDS_segments[n]):
                    same = False
        else:
            same = False
        
        return same

    @property
    def relative_coding_start(self):
        """ Returns python index of first protein nucleotide within the CDS sequence string, or 0 if no protein was generated yet."""
        if self.protein:
            return (self.protein.start - self.start)
        else:
            return 0

    @property
    def relative_coding_end(self):
        """ Returns python index of last protein nucleotide within the CDS sequence string, or the last CDS nucleotide index if no protein was generated yet."""
        if self.protein:
            return (self.protein.end - self.start)
        else:
            return self.end - self.start

class Exon(Feature):
    __slots__ = ()

    def __init__(self, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)

class UTR(Feature):

    __slots__ = ('prime',)
    def __init__(self, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)
        self.prime = "3'"

class Intron(Feature):
    __slots__ = ('intra_coding',)
    canonical_seqs = ["GT-AG", "GC-AG", "AT-AC"]

    def __init__(self, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)
        self.intra_coding = False

    @property
    def boundary(self):
        return f"{self.splice_site_donor}-{self.splice_site_acceptor}"
    
    @property
    def splice_site_donor(self):
        return self.seq[0:2]
    
    @property
    def splice_site_acceptor(self):
        return self.seq[-2:]

    @property
    def canonical(self):
        if self.boundary in Intron.canonical_seqs:
            return True
        else:
            return False
