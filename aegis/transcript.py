from __future__ import annotations
import warnings
from typing import Literal

from .feature import Feature
from .subfeatures import Exon, Intron, CDS, UTR
from .misc_features import Promoter
from .utils.genefunctions import reverse_complement
from .utils.misc import unambiguous_numbers, covered_length

class Transcript(Feature):

    __slots__ = (
        'exons', 'CDSs', 'temp_CDSs', 'temp_UTRs', 'read_UTRs_CDS', 'main',
        'miRNAs', 'renamed_exons', 'renamed_utrs',
        'promoter', 'introns', 'collapsed_exons', 'collapsed_CDS_segments', 'generated_exons'
    )

    CDSs: dict[str, CDS]
    exons: list[Exon]
    introns: list[Intron]| None
    UTRs: list[UTR]
    temp_CDSs: list[CDS|Feature]|None
    temp_UTRs: list[UTR]|None
    read_UTRs_CDS: tuple[int, int]|None
    promoter: Promoter | None

    def __init__(self, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)
        self.exons = []
        self.CDSs = {}
        self.temp_CDSs = []
        self.temp_UTRs = []
        # Span of the CDS whose read UTRs the transcript keeps instead of deriving them (see keep_read_UTRs)
        self.read_UTRs_CDS = None
        self.main = False
        self.miRNAs = []
        self.renamed_exons = False
        self.renamed_utrs = False
        self.introns = None
        self.promoter = None

        self.collapsed_exons = False
        self.collapsed_CDS_segments = False
        self.generated_exons = False

    def update(self, quiet:bool=False):
        if self.exons == []:
            # Without exons the UTRs read from the GFF are kept, and the exons are built from them and the CDS
            self.generate_CDSs(quiet=quiet)
            self.generate_exons()
            if len(self.exons) > 1:
                self.exons.sort()
        else:
            if len(self.exons) > 1:
                self.exons.sort()
            self.generate_CDSs(quiet=quiet)

        for i, c in enumerate(self.CDSs.values()):
            if len(c.CDS_segments) > 1:
                c.CDS_segments.sort()
            if i == 0:
                c_start = c.start
                c_end = c.end
            else:
                if c.start < c_start:
                    c_start = c.start
                if c.end > c_end:
                    c_end = c.end

        if self.CDSs:
            if self.strand == "+":
                for e in self.exons:
                    if e.end >= c_start and e.start <= c_end:
                        e.coding = True
                    
            elif self.strand == "-":
                for e in self.exons:
                    if e.start <= c_end and e.end >= c_start:
                        e.coding = True

    def rename(self, base_id:str, count:int, sep:str="_", digits:int=3, keep_numbering:bool=False, keep_existing_ids_if_derived_from_base_id:bool=False):

        rename = False

        if keep_existing_ids_if_derived_from_base_id:
            if base_id not in self.id:
                rename = True
        else:
            rename = True

        if rename:

            if keep_numbering and self.id_number != None:
                self.id = f"{base_id}{sep}t{self.id_number:0{digits}d}"
            else:
                if self.main:
                    self.id = f"{base_id}{sep}t{1:0{digits}d}"
                else:
                    self.id = f"{base_id}{sep}t{count:0{digits}d}"

            if self.original_id != self.id:
                self.renamed = True
                self.update_numbering()

    def rename_exons(self, base_id:str, sep:str="_", digits:int=3, keep_numbering:bool=False, keep_existing_ids_if_derived_from_base_id:bool=False):

        rename = False

        if keep_existing_ids_if_derived_from_base_id:
            for e in self.exons:
                if base_id not in e.id:
                    rename = True
        else:
            rename = True

        if rename:
            if len(self.exons) > 1:
                self.exons.sort()
            keep_numbering = keep_numbering and unambiguous_numbers(((e.start, e.end), e.id_number) for e in self.exons)

            if self.strand == "+" or self.strand == ".":
                for x, e in enumerate(self.exons):
                    if keep_numbering and e.id_number != None:
                        e.id = f"{base_id}{sep}e{e.id_number:0{digits}d}"
                    else:
                        e.id = f"{base_id}{sep}e{x+1:0{digits}d}"

                    if e.original_id != e.id:
                        self.renamed_exons = True
                        e.update_numbering()
                
            elif self.strand == "-":
                for x, e in enumerate(reversed(self.exons)):
                    if keep_numbering and e.id_number != None:
                        e.id = f"{base_id}{sep}e{e.id_number:0{digits}d}"
                    else:
                        e.id = f"{base_id}{sep}e{x+1:0{digits}d}"

                    if e.original_id != e.id:
                        self.renamed_exons = True
                        e.update_numbering()

    def rename_utrs(self, base_id:str, sep:str="_", digits:int=3, keep_numbering:bool=False, keep_existing_ids_if_derived_from_base_id:bool=False):

        rename = False

        if keep_existing_ids_if_derived_from_base_id:
            for c in self.CDSs.values():
                for u in c.UTRs:
                    if base_id not in u.id:
                        rename = True
        else:
            rename = True

        if rename:
            keep_numbering = keep_numbering and unambiguous_numbers(((u.start, u.end), u.id_number) for c in self.CDSs.values() for u in c.UTRs)

            if self.strand == "+" or self.strand == ".":
                for c in self.CDSs.values():
                    for x, u in enumerate(c.UTRs):
                        if keep_numbering and u.id_number != None:
                            u.id = f"{base_id}{sep}u{u.id_number:0{digits}d}"
                        else:
                            u.id = f"{base_id}{sep}u{x+1:0{digits}d}"
                        u.provisional = False

                        if u.original_id != u.id:
                            self.renamed_utrs = True
                            u.update_numbering()
                
            elif self.strand == "-":
                for c in reversed(self.CDSs.values()):
                    for x, u in enumerate(reversed(c.UTRs)):
                        if keep_numbering and u.id_number != None:
                            u.id = f"{base_id}{sep}u{u.id_number:0{digits}d}"
                        else:
                            u.id = f"{base_id}{sep}u{x+1:0{digits}d}"
                        u.provisional = False

                        if u.original_id != u.id:
                            self.renamed_utrs = True
                            u.update_numbering()

    def collapse_exons(self) -> list[Exon]:
        """
        Merges overlapping or directly adjacent exons into single exons. A merged exon keeps the id
        and attributes of its 5'-most exon; exons that are not merged are kept as they are.
        Returns the merged exons.
        """
        merged_exons = []
        if len(self.exons) > 1:

            self.exons.sort()
            groups = [[self.exons[0]]]
            cur_end = self.exons[0].end
            for e in self.exons[1:]:
                if e.start <= cur_end + 1:
                    groups[-1].append(e)
                    cur_end = max(cur_end, e.end)
                else:
                    groups.append([e])
                    cur_end = e.end

            if len(groups) < len(self.exons):
                merged = []
                for group in groups:
                    if len(group) == 1:
                        merged.append(group[0])
                        continue
                    rep = group[0] if self.strand != "-" else max(group, key=lambda e: e.end)
                    exon = Exon(rep.id, rep.ch, rep.source, rep.feature, rep.strand, group[0].start, max(e.end for e in group), rep.score, [self.id])
                    if rep.misc_attributes:
                        exon.misc_attributes = rep.misc_attributes[:]
                    merged.append(exon)
                    merged_exons.append(exon)
                self.collapsed_exons = True
                self.exons = merged
                self.update()
        return merged_exons

    def collapse_CDS_segments(self):
        """
        Merges CDS segments that are rows of one continuous coding stretch: duplicated or nested rows,
        overlaps of 3 bp or more, and adjacent or 1-2 bp overlapping rows whose phases continue the
        reading frame across the merge. Overlaps of 1-2 bp whose phases instead read the shared bases
        twice (programmed -1/-2 frameshifts, e.g. PEG10), or that have no phases, are kept apart.
        A merged segment takes the phase of its 5'-most row, where it starts. Phases describe each row
        on its own (GFF3), so all other segments keep their annotated phases.
        """
        if self.temp_CDSs and not self.CDSs:
            self.generate_CDSs()

        parents = [self.id]
        for cds in self.CDSs.values():

            if len(cds.CDS_segments) > 1:
                cds.CDS_segments.sort()
                groups = [[cds.CDS_segments[0]]]
                cur_end = cds.CDS_segments[0].end
                for seg in cds.CDS_segments[1:]:
                    if seg.start <= cur_end + 1:
                        # Adjacent segments and 1-2 bp overlaps only collapse when merging keeps the reading frame:
                        # the leftover bases of the upstream segment, minus the overlap, must complete a codon
                        # with the downstream segment's phase. Longer overlaps always collapse
                        prev_seg = groups[-1][-1]
                        overlap = cur_end - seg.start + 1
                        phase_compatible = True
                        if overlap <= 2:
                            if prev_seg.phase is not None and seg.phase is not None:
                                if cds.strand != "-":
                                    prev_lo = (prev_seg.size - prev_seg.phase) % 3
                                    phase_compatible = (prev_lo - overlap + seg.phase) % 3 == 0
                                else:
                                    seg_lo = (seg.size - seg.phase) % 3
                                    phase_compatible = (seg_lo - overlap + prev_seg.phase) % 3 == 0
                            else:
                                # Without phases a 1-2 bp overlap may be a programmed frameshift, so only adjacent rows merge
                                phase_compatible = overlap == 0

                        if phase_compatible:
                            cur_end = max(cur_end, seg.end)
                            groups[-1].append(seg)
                            continue

                    groups.append([seg])
                    cur_end = seg.end

                if len(groups) < len(cds.CDS_segments):
                    merged = []
                    for group in groups:
                        # A merged row starts where its 5'-most row starts, so it takes that row's phase
                        rep = group[0] if cds.strand != "-" else max(group, key=lambda s: s.end)
                        merged.append(Feature(cds.id, rep.ch, rep.source, "CDS", rep.strand, group[0].start, max(s.end for s in group), rep.score, parents, phase=rep.phase))
                        if rep.misc_attributes:
                            merged[-1].misc_attributes = rep.misc_attributes[:]
                    cds.CDS_segments = merged
                    self.collapsed_CDS_segments = True
                    cds.update()

        if self.collapsed_CDS_segments:
            self.update()

    def clear_UTRs(self):
        for c in self.CDSs.values():
            c.clear_UTRs()
        self.temp_UTRs = None
        self.read_UTRs_CDS = None

    def keep_read_UTRs(self):
        """
        Makes the transcript keep the UTRs read from the GFF instead of deriving them from its exons, for as
        long as its CDS keeps the span it was read with and the UTRs fit the exons and CDS (see read_UTRs_fit).
        Called before the CDS is generated from its rows.
        """
        if self.temp_UTRs and self.temp_CDSs:
            self.read_UTRs_CDS = (min(cs.start for cs in self.temp_CDSs), max(cs.end for cs in self.temp_CDSs))

    def read_UTRs_fit(self) -> bool:
        """
        Whether the read UTRs the transcript keeps still fit it: its CDS has the span the UTRs were read with,
        and each UTR is on the CDS strand, outside the CDS and within one exon. Parts of a UTR beyond the
        first or last exon are allowed, as coordinate polishing extends that exon (as for CDS rows).
        """
        if self.read_UTRs_CDS is None:
            return False
        span_start = min(e.start for e in self.exons)
        span_end = max(e.end for e in self.exons)
        for c in self.CDSs.values():
            if (c.start, c.end) != self.read_UTRs_CDS:
                return False
            for u in c.UTRs:
                if u.strand != c.strand or (u.start <= c.end and c.start <= u.end):
                    return False
                start, end = max(u.start, span_start), min(u.end, span_end)
                if start <= end and not any(e.start <= start and end <= e.end for e in self.exons):
                    return False
        return True

    def generate_promoter(self, promoter_size:int, ch_size:int, promoter_type:str = "standard"):
        """
        promoter_type (str): Defines the reference point for the promoter.
            - standard (default): Promoter based on 'promoter_size' is generated upstream of the transcript's start site (TSS)
            - upstream_ATG : Promoter based on 'promoter_size' is generated upstream of the main CDS's start codon (ATG). If no CDS, falls back to standard.
            - standard_plus_up_to_ATG : Promoter based on 'promoter_size' is generated upstream of the transcript's start site (TSS) and any gene sequence up to the start codon (ATG) is also added. If no CDS, falls back to standard.
        """

        valid_options = ["standard", "upstream_ATG", "standard_plus_up_to_ATG"]

        if promoter_type not in valid_options:
            raise ValueError(f"promoter_type={promoter_type} selected is not a valid option. Choose from: {valid_options}.")

        prom_id = self.id + "_promoter"
        if self.strand == "+":
            temp_start = self.start - promoter_size
            temp_end = self.start - 1
            if promoter_type != "standard":
                if self.CDSs != {}:
                    for c in self.CDSs.values():
                        if c.main:
                            temp_end = c.start - 1
                            if promoter_type == "standard_plus_up_to_ATG":
                                temp_start = self.start - promoter_size
                            else:
                                temp_start = c.start - promoter_size
                else:
                    promoter_type = "standard"
        
        elif self.strand == "-":
            temp_start = self.end + 1
            temp_end = self.end + promoter_size
            if promoter_type != "standard":
                if self.CDSs != {}:
                    for c in self.CDSs.values():
                        if c.main:
                            temp_start = c.end + 1
                            if promoter_type == "standard_plus_up_to_ATG":
                                temp_end = self.end + promoter_size
                            else:
                                temp_end = c.end + promoter_size
                else:
                    promoter_type = "standard"
        else:
            print(f"Warning: {self.id} does not have a strand.")
            promoter_type = "none"
        
        if promoter_type != "none":
            # Making sure genes at the end of chromosomes or contigs have a defined
            # promoter, even if it is smaller than the requested.
            if temp_start < 1:
                temp_start = 1
            # making sure that if gene starts at base 1, no promoter seq is given
            if temp_end < 1:
                temp_start = 1
                temp_end = 0

            # similar thing at the other contig/chr end
            if temp_end > ch_size:
                temp_end = ch_size
            if temp_start > ch_size:
                temp_start = 1
                temp_end = 0

            self.promoter = Promoter(promoter_type, prom_id, self.ch, self.source, self.feature, self.strand, temp_start, temp_end, self.score, [self.id])

    def generate_best_protein(
        self,
        mode: Literal["orf", "orf_or_end", "orf_or_start"] = "orf",
        start_codons: tuple[str, ...] | None = None,
        stop_codons: tuple[str, ...] | None = None,
        min_codon_len: int = 2,
        enforce_start_codon: bool = True,
        must_have_stop: bool = True,
        tolerated_stops: int | None = 0,
        orf_choice_mode: Literal["longest", "earliest"] = "longest",
        max_nucleotide_trim: int | None = None,
        always_resolve_strand: bool = True,
        ignore_ambiguous_strands: bool = False,
        quiet: bool = True,
        table: int | str | dict[str, str] = 1,
        initiator_methionine: bool = True,
        alt_start_methionine: bool = True,
    ):
        """
        Extracts the best protein/ORF across the spliced exons of the transcript,
        updates the CDS and UTRs, and updates the transcript strand if ambiguous.
        """
        self.CDSs = {}
        self.clear_UTRs()
        self.temp_CDSs = []

        cds_id = f"{self.id}_CDS1"
        for e in self.exons:
            self.temp_CDSs.append(Feature(feature_id=cds_id, ch=self.ch, source=self.source, feature="CDS", strand=self.strand, start=e.start, end=e.end, score=self.score, parents=[self.id]))

        self.generate_CDSs(quiet=quiet)

        candidate_cds = self.CDSs.get(cds_id)
        if candidate_cds is not None:
            candidate_cds.generate_protein(
                mode=mode,
                max_nucleotide_trim=max_nucleotide_trim,
                tolerated_stops=tolerated_stops,
                orf_choice_mode=orf_choice_mode,
                must_have_stop=must_have_stop,
                enforce_start_codon=enforce_start_codon,
                min_codon_len=min_codon_len,
                start_codons=start_codons,
                stop_codons=stop_codons,
                correct_CDS=True,
                always_resolve_strand=always_resolve_strand,
                ignore_ambiguous_strands=ignore_ambiguous_strands,
                quiet=quiet,
                table=table,
                initiator_methionine=initiator_methionine,
                alt_start_methionine=alt_start_methionine,
            )

            if candidate_cds.protein is not None:
                self.strand = candidate_cds.strand
                for e in self.exons:
                    e.strand = self.strand
                self.update(quiet=quiet)
            else:
                if cds_id in self.CDSs:
                    del self.CDSs[cds_id]
                self.coding = False
                self.clear_UTRs()
                for e in self.exons:
                    e.coding = False
                self.update(quiet=quiet)

    def keep_unchanged_subfeatures(self, old_CDSs:list[CDS], old_UTRs:list[UTR], quiet:bool=False):
        """
        After generate_best_protein, keeps what the rework did not change: a replaced CDS with the strand,
        segments and phases of the new one is kept instead of it (with its id and attributes, taking the new
        protein), and
        UTRs with the coordinates of replaced ones are kept (see generate_UTRs). Otherwise the CDS is named
        <transcript>_CDS001. New UTRs keep provisional ids until Annotation.update names them.
        """
        if not self.CDSs:
            return []
        c = next(iter(self.CDSs.values()))
        model = (c.strand, [(cs.start, cs.end, cs.phase) for cs in c.CDS_segments])
        unchanged = next((old for old in old_CDSs if (old.strand, [(cs.start, cs.end, cs.phase) for cs in old.CDS_segments]) == model), None)
        if unchanged is not None:
            unchanged.protein = c.protein
            if unchanged.protein is not None:
                unchanged.protein.id = f"{unchanged.id}.prot"
            c = unchanged
        else:
            c.rename(base_id=self.id, base_gene_id="", count=1)
        self.CDSs = {c.id: c}
        c.UTRs = [u for u in old_UTRs if u.strand == c.strand]
        self.update(quiet=quiet)

    def almost_equal(self, other:Transcript):
        almost_equal = True
        if len(self.exons) != len(other.exons):
            almost_equal = False
        else:
            if len(self.exons) > 1:
                self.exons.sort()
            if len(other.exons) > 1:
                other.exons.sort()
            for n, exon in enumerate(self.exons):
                if exon.start != other.exons[n].start or exon.end != other.exons[n].end:
                    almost_equal = False
                    break
        return almost_equal

    def _detect_polycistronic(self) -> list[str]:
        """
        Evaluates raw parsed temp_CDSs for polycistronic structures, overlapping segments,
        and formatting inconsistencies prior to CDS generation.
        Returns a list of warning flags (empty if standard monocistronic):
          - 'overlapping_CDS_segments': Physical coordinate overlap between adjacent CDS segments.
          - 'possible_polycistronic_transcript': Structural evidence of multiple CDSs (grouped IDs or overlapping distinct IDs).
            Exemplified by the canonical GFF3 specification (canonical_gene_example.gff3),
            where mRNA00003 defines two distinct logical CDSs (cds00003 and cds00004) under the same transcript.
          - 'independent_CDS_segment_ids': Distinct IDs per segment without polycistronic evidence (GFF row-ID style).
        """
        if not self.temp_CDSs or len(self.temp_CDSs) <= 1:
            return []

        self.temp_CDSs.sort()
        segs = self.temp_CDSs
        warnings = []

        auto_prefix = f"{self.id}_CDS_"
        all_auto_ids = all(getattr(s, "id", "").startswith(auto_prefix) for s in segs)

        has_overlap = any(
            segs[i].start <= segs[i - 1].end
            for i in range(1, len(segs))
        )

        if has_overlap:
            warnings.append("overlapping_CDS_segments")

        unique_ids = {s.id for s in segs if getattr(s, "id", "")}
        num_ids = len(unique_ids)
        num_segs = len(segs)

        # Polycistronic evidence:
        # 1. Grouped IDs: e.g. seg1 & seg2 share ID_A, seg3 & seg4 share ID_B (1 < num_ids < num_segs)
        # 2. Overlapping segments with distinct IDs (e.g. dual-coding / overlapping ORFs)
        is_polycistronic = (num_ids > 1 and num_ids < num_segs) or (has_overlap and num_ids > 1)

        if is_polycistronic:
            warnings.append("possible_polycistronic_transcript")
        elif num_ids > 1 and not all_auto_ids:
            warnings.append("independent_CDS_segment_ids")

        return warnings

    def generate_CDSs(self, quiet:bool=False):
        """
        Creates a CDS object based on the CDS segments of a transcript, and gives it its UTRs (see place_UTRs).
        AEGIS enforces 1 logical CDS per transcript. Any polycistronic or
        multi-CDS structures (e.g. theoretical multi-CDS annotations as in
        canonical_gene_example.gff3) are flagged via _detect_polycistronic()
        for warning logs in Annotation.warnings.
        """
        if self.temp_CDSs:
            if len(self.temp_CDSs) > 1:
                self.temp_CDSs.sort()

            parents = [self.id]
            self.coding = True

            if self.strand == "+":
                temp_id = self.temp_CDSs[0].id
            else:
                temp_id = self.temp_CDSs[-1].id

            self.CDSs[temp_id] = CDS(
                self.temp_CDSs, temp_id,
                self.temp_CDSs[0].ch,
                self.temp_CDSs[0].source,
                self.temp_CDSs[0].feature,
                self.temp_CDSs[0].strand,
                self.temp_CDSs[0].start,
                self.temp_CDSs[-1].end,
                self.temp_CDSs[0].score,
                parents
            )

            self.temp_CDSs = None

        self.determine_main_CDS()
        self.place_UTRs()

        for c in self.CDSs.values():
            c.update()

    def place_UTRs(self):
        """
        Gives each CDS its UTRs. The UTRs read from the GFF are used while the transcript keeps them (see
        keep_read_UTRs) and they fit it, or when it has no exons to derive UTRs from (its exons are then built
        from the CDS and UTRs); otherwise the UTRs are derived from the exons, and read UTRs are no longer kept.
        """
        if self.temp_UTRs and (self.read_UTRs_CDS is not None or not self.exons):
            self.assign_UTRs()
        if not self.exons or self.read_UTRs_fit():
            self.update_UTRs()
        else:
            self.read_UTRs_CDS = None
            self.generate_UTRs()

    def determine_main_CDS(self):
        if self.CDSs != {}:
            x = 0
            main = ""
            for c in self.CDSs.values():
                c.main = False
            for c_id, c in self.CDSs.items():
                x += 1
                if x == 1:
                    main = c_id
                    size = c.size
                    start = c.start
                else:
                    if c.size > size:
                        main = c_id
                        size = c.size
                        start = c.start
                    elif c.size == size:
                        if c.start < start:
                            main = c_id
                            size = c.size
                            start = c.start
            if main:
                self.CDSs[main].main = True

    def assign_UTRs(self):
        if self.temp_UTRs:
            self.temp_UTRs.sort()
            for c in self.CDSs.values():
                c.UTRs = self.temp_UTRs[:]
            self.temp_UTRs = None
            self.update_UTRs()

    def generate_UTRs(self):
        """
        Derives the UTRs of each CDS from the exons. A derived UTR with the coordinates of a current
        UTR, or of one read from the GFF, is that UTR (keeping its id and attributes); only new UTRs
        get provisional ids (<CDS id>_u<n>, see UTR.provisional).
        """
        previous = {(u.start, u.end): u for u in (self.temp_UTRs or [])}
        self.temp_UTRs = None
        for c in self.CDSs.values():
            previous.update({(u.start, u.end): u for u in c.UTRs})
            c.UTRs = []
            parents = [self.id]
            c.update_span()
            for exon in self.exons:
                if c.strand != exon.strand:
                    continue
                if exon.end <= c.end and exon.start >= c.start:
                    continue
                if exon.end < c.start:
                    c.UTRs.append(UTR("", exon.ch, exon.source, "UTR",
                                      exon.strand, exon.start, exon.end,
                                      exon.score, parents))
                elif exon.start < c.start:
                    c.UTRs.append(UTR("", exon.ch, exon.source, "UTR",
                                      exon.strand, exon.start, c.start-1,
                                      exon.score, parents))
                if exon.start > c.end:
                    c.UTRs.append(UTR("", exon.ch, exon.source, "UTR",
                                      exon.strand, exon.start, exon.end,
                                      exon.score, parents))
                elif exon.end > c.end:
                    c.UTRs.append(UTR("", exon.ch, exon.source, "UTR",
                                      exon.strand, c.end+1, exon.end,
                                      exon.score, parents))
            c.UTRs = sorted(previous.pop((u.start, u.end), u) for u in c.UTRs)
            # New UTRs are numbered in transcript order, skipping the ids of kept ones
            used_ids = {u.id for u in c.UTRs if u.id}
            n = 0
            for u in (c.UTRs if c.strand != "-" else reversed(c.UTRs)):
                u.parents = parents[:]
                if not u.id:
                    n += 1
                    while f"{c.id}_u{n}" in used_ids:
                        n += 1
                    u.id = f"{c.id}_u{n}"
                    u.provisional = True
        self.update_UTRs()

    def update_UTRs(self):
        for c in self.CDSs.values():
            if c.UTRs:
                c.UTRs.sort()
                if c.strand == "+":
                    for u in c.UTRs:
                        if u.start < c.start:
                            u.prime = "5'"
                            u.feature = "five_prime_UTR"
                        else:
                            u.prime = "3'"
                            u.feature = "three_prime_UTR"
                elif c.strand == "-":
                    for u in c.UTRs:
                        if u.end > c.end:
                            u.prime = "5'"
                            u.feature = "five_prime_UTR"
                        else:
                            u.prime = "3'"
                            u.feature = "three_prime_UTR"

                c.full_UTR_exons = len(self.exons) - len(c.CDS_segments)

    def generate_exons(self):
        """
        Creates exon features from CDSs and UTRs or directly from the transcript
        """
        temp_fts = []
        for c in self.CDSs.values():
            if c.main:
                for cs in c.CDS_segments:
                    temp_fts.append(cs)
                for u in c.UTRs:
                    temp_fts.append(u)

        parents = [self.id]
        # Exons reconstructed from CDS/UTRs
        if temp_fts != []:
            self.exons = [ Exon("temp", ft.ch, ft.source, "exon", ft.strand, ft.start, ft.end, ft.score, parents) for ft in temp_fts ]
            self.collapse_exons()

        # Exons rebuilt from the transcript
        else:
            self.exons = [Exon(f"temp", self.ch, self.source, "exon", self.strand, self.start, self.end, self.score, parents)]

        self.rename_exons(base_id=self.id)
        self.generated_exons = True

    def generate_introns(self):
        if len(self.exons) > 1:
            self.exons.sort()
        parents = [self.id]
        # Track the furthest exon end so far, as with overlapping/nested exons the
        # previous exon in (start, end) order may not reach as far as earlier ones
        intervals = []
        covered_end = None
        for exon in self.exons:
            if covered_end is not None and exon.start > covered_end + 1:
                intervals.append((covered_end + 1, exon.start - 1))
            if covered_end is None or exon.end > covered_end:
                covered_end = exon.end
        # Introns are kept in coordinate order but numbered 5' to 3', as exons are
        total = len(intervals)
        self.introns = []
        for n, (start, end) in enumerate(intervals):
            number = total - n if self.strand == "-" else n + 1
            self.introns.append(Intron(f"{self.id}_intron_{number}", self.ch,
                                       self.source, "intron", self.strand,
                                       start, end, self.score, parents))
        if self.strand == "+":
            for i in self.introns:
                for c in self.CDSs.values():
                    if c.main:
                        if i.start > c.start and i.end < c.end:
                            i.intra_coding = True
        elif self.strand == "-":
            for i in self.introns:
                for c in self.CDSs.values():
                    if c.main:
                        if i.end < c.end and i.start > c.start:
                            i.intra_coding = True

    def clear_promoter(self):
        self.promoter = None

    @property
    def seq(self) -> str:
        if not self._ACTIVE_GENOME:
            raise ValueError("No genome loaded and you are trying to access the sequence. Load your genome together with your annotation.")
        else:
            transcript_seq = ""
            if self.strand == "-":
                for exon in reversed(self.exons):
                    transcript_seq += exon.seq
            else:
                for exon in self.exons:
                    transcript_seq += exon.seq

            return transcript_seq

    @property
    def hard_seq(self) -> str:
        if not self._ACTIVE_HARD_GENOME:
            raise ValueError("No hard masked genome loaded and you are trying to access the hard masked sequence. Load your hard masked genome together with your annotation.")
        else:
            transcript_seq = ""
            if self.strand == "-":
                for exon in reversed(self.exons):
                    transcript_seq += exon.hard_seq
            else:
                for exon in self.exons:
                    transcript_seq += exon.hard_seq
            return transcript_seq

    @property
    def seqs(self) -> list[str]:
        if not self._ACTIVE_GENOME:
            raise ValueError("No genome loaded and you are trying to access the sequence. Load your genome together with your annotation.")
        else:
            if not self.exons:
                return ["", ""]
            if self.ch not in self._ACTIVE_GENOME.scaffolds:
                warnings.warn(
                    f"Transcript '{self.id}' is on contig '{self.ch}', which is missing from genome '{self._ACTIVE_GENOME.name}'. Returning empty sequences.",
                    category=UserWarning,
                    stacklevel=2,
                )
                return ["", ""]
            scf_seq = self._ACTIVE_GENOME.scaffolds[self.ch].seq
            scf_len = len(scf_seq)
            if len(self.exons) > 1:
                self.exons.sort()
            fw_seq = "".join(scf_seq[max(0, exon.start - 1):min(scf_len, exon.end)] for exon in self.exons)
            return [fw_seq, reverse_complement(fw_seq)]

    @property
    def hard_seqs(self) -> list[str]:
        if not self._ACTIVE_HARD_GENOME:
            raise ValueError("No hard masked genome loaded and you are trying to access the hard masked sequence. Load your hard masked genome together with your annotation.")
        else:
            if not self.exons:
                return ["", ""]
            if self.ch not in self._ACTIVE_HARD_GENOME.scaffolds:
                warnings.warn(
                    f"Transcript '{self.id}' is on contig '{self.ch}', which is missing from hard-masked genome '{self._ACTIVE_HARD_GENOME.name}'. Returning empty sequences.",
                    category=UserWarning,
                    stacklevel=2,
                )
                return ["", ""]
            scf_seq = self._ACTIVE_HARD_GENOME.scaffolds[self.ch].seq
            scf_len = len(scf_seq)
            if len(self.exons) > 1:
                self.exons.sort()
            fw_seq = "".join(scf_seq[max(0, exon.start - 1):min(scf_len, exon.end)] for exon in self.exons)
            return [fw_seq, reverse_complement(fw_seq)]
    
    @property
    def size(self) -> int:
        """Mature transcript length: bases covered by its exons, each counted once (genomic_span includes introns)."""
        return covered_length((e.start, e.end) for e in self.exons)

    @property
    def CDS_size(self):
        size = 0
        for c in self.CDSs.values():
            if c.main:
                size = c.size
                break
        return size

    @property
    def coding_ratio(self):
        if self.CDS_size != 0:
            return round((self.CDS_size / self.size), 2)
        else:
            return 0