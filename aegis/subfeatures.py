from __future__ import annotations

import warnings
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
        if self.CDS_segments:
            if len(self.CDS_segments) > 1:
                self.CDS_segments.sort()
            self.start = self.CDS_segments[0].start
            self.end = self.CDS_segments[-1].end
        self.full_UTR_exons = 0
        self.protein = None
        self.update()

    def update(self):
        if self.CDS_segments:
            if len(self.CDS_segments) > 1:
                self.CDS_segments.sort()
            self.start = self.CDS_segments[0].start
            self.end = self.CDS_segments[-1].end
        self.update_phase()
        self.update_frame()

    @property
    def size(self):
        size = 0
        for segment in self.CDS_segments:
            size += segment.size
        return size

    def update_phase(self, override: bool = False, full_override: bool = False):
        if not self.CDS_segments:
            return

        if len(self.CDS_segments) > 1:
            self.CDS_segments.sort()

        if self.strand == "-":
            initial_seg = self.CDS_segments[-1]
            working_segs = list(reversed(self.CDS_segments))
        else:
            initial_seg = self.CDS_segments[0]
            working_segs = self.CDS_segments

        # If not overriding and all segments already have valid phases, preserve them
        if not override and not full_override and all(cs.phase in (0, 1, 2) for cs in self.CDS_segments):
            self.phase = initial_seg.phase
            return

        if full_override:
            initial_phase = 0
        else:
            initial_phase = initial_seg.phase if initial_seg.phase in (0, 1, 2) else 0

        leftover = 0
        for i, cs in enumerate(working_segs):
            if i == 0:
                cs.phase = initial_phase
            else:
                if leftover == 0:
                    cs.phase = 0
                else:
                    cs.phase = 3 - leftover
            leftover = (cs.size - cs.phase) % 3

        self.phase = initial_phase

    def update_frame(self):
        if not self.CDS_segments:
            return

        if len(self.CDS_segments) > 1:
            self.CDS_segments.sort()

        if self.strand == "+":
            for cs in self.CDS_segments:
                if cs.phase is not None:
                    frame = (cs.start + cs.phase) % 3
                    if frame == 0:
                        frame = 3
                    cs.frame = frame

        if self.strand == "-":
            for cs in reversed(self.CDS_segments):
                if cs.phase is not None:
                    frame = (cs.end - cs.phase) % 3
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
            if not self.CDS_segments:
                return ["", ""]
            if self.ch not in self._ACTIVE_GENOME.scaffolds:
                warnings.warn(
                    f"CDS '{self.id}' is on contig '{self.ch}', which is missing from genome '{self._ACTIVE_GENOME.name}'. Returning empty sequences.",
                    category=UserWarning,
                    stacklevel=2,
                )
                return ["", ""]
            scf_seq = self._ACTIVE_GENOME.scaffolds[self.ch].seq
            scf_len = len(scf_seq)
            if len(self.CDS_segments) > 1:
                self.CDS_segments.sort()
            fw_seq = "".join(scf_seq[max(0, cs.start - 1):min(scf_len, cs.end)] for cs in self.CDS_segments)
            return [fw_seq, reverse_complement(fw_seq)]

    @property
    def hard_seqs(self) -> list[str]:
        if not self._ACTIVE_HARD_GENOME:
            raise ValueError("No hard masked genome loaded and you are trying to access the hard masked sequence. Load your hard masked genome together with your annotation.")
        else:
            if not self.CDS_segments:
                return ["", ""]
            if self.ch not in self._ACTIVE_HARD_GENOME.scaffolds:
                warnings.warn(
                    f"CDS '{self.id}' is on contig '{self.ch}', which is missing from hard-masked genome '{self._ACTIVE_HARD_GENOME.name}'. Returning empty sequences.",
                    category=UserWarning,
                    stacklevel=2,
                )
                return ["", ""]
            scf_seq = self._ACTIVE_HARD_GENOME.scaffolds[self.ch].seq
            scf_len = len(scf_seq)
            if len(self.CDS_segments) > 1:
                self.CDS_segments.sort()
            fw_seq = "".join(scf_seq[max(0, cs.start - 1):min(scf_len, cs.end)] for cs in self.CDS_segments)
            return [fw_seq, reverse_complement(fw_seq)]

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

    def generate_protein(
        self,
        mode: Literal["start", "end", "orf", "orf_or_end", "orf_or_start"] = "end",
        max_nucleotide_trim: int | None = None,
        tolerated_stops: int | None = 0,
        orf_choice_mode: Literal["longest", "earliest"]="longest",
        must_have_stop: bool = False,
        enforce_start_codon: bool = True,
        min_codon_len: int = 2,
        start_codons: tuple[str, ...] | None = None,
        stop_codons: tuple[str, ...] | None = None,
        correct_CDS: bool = False,
        always_resolve_strand: bool = True,
        ignore_ambiguous_strands: bool = False,
        quiet: bool = True,
        adjust_internal_shifts: Literal["intra_exon", "all", "none"] | bool = "intra_exon",
        table: int | str | dict[str, str] = 1,
    ):
        if table is None:
            table = 1

        if len(self.CDS_segments) > 1:
            self.CDS_segments.sort()

        if (self.strand == "." or self.strand == "?") and not ignore_ambiguous_strands:
            seq_fw, seq_rv = self.seqs
            fw_orf = choose_orf(find_ORFs(seq_fw, min_codon_len=min_codon_len, enforce_start_codon=enforce_start_codon, must_have_stop=must_have_stop, tolerated_stops=tolerated_stops, start_codons=start_codons, stop_codons=stop_codons, table=table), mode=orf_choice_mode)
            rv_orf = choose_orf(find_ORFs(seq_rv, min_codon_len=min_codon_len, enforce_start_codon=enforce_start_codon, must_have_stop=must_have_stop, tolerated_stops=tolerated_stops, start_codons=start_codons, stop_codons=stop_codons, table=table), mode=orf_choice_mode)

            has_orf = len(fw_orf[0]) > 0 or len(rv_orf[0]) > 0
            if has_orf or always_resolve_strand:
                new_strand = "+" if len(fw_orf[0]) >= len(rv_orf[0]) else "-"
                strand_changed = (new_strand != self.strand)
                self.strand = new_strand
                for cs in self.CDS_segments:
                    cs.strand = new_strand
                if strand_changed:
                    self.update_phase(override=True)
                    self.update_frame()
                else:
                    self.update()

        cds_phase = self.phase if self.phase in (1, 2) else 0

        # Determine internal shift adjustment mode: "intra_exon" (default), "all", or "none"
        if adjust_internal_shifts is True or adjust_internal_shifts == "intra_exon":
            adj_mode = "intra_exon"
        elif adjust_internal_shifts == "all":
            adj_mode = "all"
        else:
            adj_mode = "none"

        working_segs = self.CDS_segments if self.strand != "-" else list(reversed(self.CDS_segments))
        has_internal_shift = False
        shifts_to_adjust = False
        if len(working_segs) > 1 and all(cs.phase in (0, 1, 2) for cs in working_segs):
            prev_lo = (working_segs[0].size - (working_segs[0].phase or 0)) % 3
            prev_seg = working_segs[0]
            for cs in working_segs[1:]:
                if (prev_lo + (cs.phase or 0)) % 3 != 0:
                    has_internal_shift = True
                    is_contig = (cs.start <= prev_seg.end + 2) if self.strand != "-" else (prev_seg.start <= cs.end + 2)
                    if adj_mode == "all" or (adj_mode == "intra_exon" and is_contig):
                        shifts_to_adjust = True
                prev_lo = (cs.size - (cs.phase or 0)) % 3
                prev_seg = cs

        if not quiet and has_internal_shift:
            if shifts_to_adjust:
                print(f"Warning: {self.id} has an internal phase shift at segment boundaries; junction bases were adjusted.")
            else:
                print(f"Warning: {self.id} has a phase mismatch across introns; mature mRNA spliced continuously without adjusting junction bases.")

        if shifts_to_adjust and mode in ("start", "end"):
            coding_parts = []
            segment_intervals = []
            pending_leftover = ""
            n_segs = len(working_segs)
            for i, cs in enumerate(working_segs):
                p = cs.phase or 0
                s = cs.seq
                l = len(s)
                lo = (l - p) % 3
                if i == 0:
                    trim_5p = p
                    coding_parts.append(s[p : l - lo])
                    pending_leftover = s[l - lo :]
                else:
                    leading = s[:p]
                    if len(pending_leftover) + len(leading) == 3:
                        trim_5p = 0
                        coding_parts.append(pending_leftover + leading)
                    else:
                        trim_5p = p
                    coding_parts.append(s[p : l - lo])
                    pending_leftover = s[l - lo :]

                if i == n_segs - 1:
                    trim_3p = lo
                else:
                    next_cs = working_segs[i + 1]
                    next_p = next_cs.phase or 0
                    if lo + next_p == 3:
                        trim_3p = 0
                    else:
                        trim_3p = lo

                if self.strand != "-":
                    g_start = cs.start + trim_5p
                    g_end = cs.end - trim_3p
                else:
                    g_end = cs.end - trim_5p
                    g_start = cs.start + trim_3p

                if g_start <= g_end:
                    segment_intervals.append((int(g_start), int(g_end)))

            coding_seq = "".join(coding_parts)
            nucleotide_surplus = True
            relative_coding_start = cds_phase
            relative_coding_end = relative_coding_start + len(coding_seq) - 1

            if self.strand == "-":
                corrected_segments = list(reversed(segment_intervals))
            else:
                corrected_segments = segment_intervals
        else:
            coding_seq, nucleotide_surplus, relative_coding_start, relative_coding_end = trim_surplus(
                self.seq, 
                mode=mode, 
                max_nucleotide_trim=max_nucleotide_trim, 
                orf_choice_mode=orf_choice_mode, 
                must_have_stop=must_have_stop, 
                tolerated_stops=tolerated_stops, 
                enforce_start_codon=enforce_start_codon, 
                start_codons=start_codons, 
                stop_codons=stop_codons, 
                min_codon_len=min_codon_len,
                phase=cds_phase,
                table=table,
            )
            corrected_segments = (
                map_relative_to_genomic(segments=self.CDS_segments, rel_start=relative_coding_start, rel_end=relative_coding_end, strand=self.strand)
                if coding_seq and len(coding_seq) >= 3 and relative_coding_end >= relative_coding_start
                else []
            )

        if coding_seq and len(coding_seq) >= 3 and relative_coding_end >= relative_coding_start:

            protein_seq = translate(coding_seq, table=table)

            if corrected_segments:
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
                    self.update()

                self.protein = Protein(prot_id=f"{self.id}.prot", sequence=protein_seq, chrom=self.ch, start=protein_start, end=protein_end, nucleotide_surplus=nucleotide_surplus, readthrough=mode, nuc_seq=coding_seq, segments=tuple(corrected_segments))

                if not quiet and nucleotide_surplus:
                    print(f"{self.id} has a nucleotide surplus when translating to protein, the annotated CDS might be incorrect.")
            else:
                self.protein = None
                if not quiet:
                    print(f"{self.id} CDS could not be mapped to genomic coordinates with mode={mode}")
        else:
            self.protein = None
            if not quiet:
                if nucleotide_surplus:
                    print(f"{self.id} has a nucleotide surplus and could not be translated to a protein with mode={mode}")
                else:
                    print(f"{self.id} CDS could not be translated to a protein with mode={mode}")

    def clear_protein(self):
        self.protein = None

    def equal_segments(self, other:CDS):
        if len(self.CDS_segments) > 1:
            self.CDS_segments.sort()
        if len(other.CDS_segments) > 1:
            other.CDS_segments.sort()
        same = True
        if len(self.CDS_segments) == len(other.CDS_segments):
            for n, segment in enumerate(self.CDS_segments):
                if not segment.equal_sequence(other.CDS_segments[n]):
                    same = False
        else:
            same = False
        
        return same

    def _calculate_relative_coding_coords(self) -> tuple[int, int]:
        """Calculates 0-based slice indices within self.seq corresponding to the protein."""
        if not self.protein or not self.CDS_segments:
            return 0, max(0, self.size - 1)

        prot_start = self.protein.start
        prot_end = self.protein.end

        if len(self.CDS_segments) > 1:
            self.CDS_segments.sort()
        working_segs = self.CDS_segments if self.strand != "-" else reversed(self.CDS_segments)

        rel_start = None
        rel_end = None
        offset = 0

        for cs in working_segs:
            if self.strand != "-":
                if rel_start is None and cs.start <= prot_start <= cs.end:
                    rel_start = offset + (prot_start - cs.start)
                if rel_end is None and cs.start <= prot_end <= cs.end:
                    rel_end = offset + (prot_end - cs.start)
            else:
                if rel_start is None and cs.start <= prot_end <= cs.end:
                    rel_start = offset + (cs.end - prot_end)
                if rel_end is None and cs.start <= prot_start <= cs.end:
                    rel_end = offset + (cs.end - prot_start)
            offset += cs.size

        final_start = rel_start if rel_start is not None else 0
        final_end = rel_end if rel_end is not None else max(0, self.size - 1)
        return final_start, final_end

    @property
    def relative_coding_start(self) -> int:
        """ Returns python index of first protein nucleotide within the CDS sequence string, or 0 if no protein was generated yet."""
        return self._calculate_relative_coding_coords()[0]

    @property
    def relative_coding_end(self) -> int:
        """ Returns python index of last protein nucleotide within the CDS sequence string, or the last CDS nucleotide index if no protein was generated yet."""
        return self._calculate_relative_coding_coords()[1]

class Exon(Feature):
    __slots__ = ()

    def __init__(self, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)

class UTR(Feature):

    __slots__ = ('prime',)
    def __init__(self, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)
        feat_lower = str(feature).lower()
        if "5" in feat_lower or "five" in feat_lower:
            self.prime = "5'"
        else:
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
        return self.seq[0:2].upper()
    
    @property
    def splice_site_acceptor(self):
        return self.seq[-2:].upper()

    @property
    def canonical(self):
        if self.boundary in Intron.canonical_seqs:
            return True
        else:
            return False
