from __future__ import annotations
from typing import Literal

from .feature import Feature
from .subfeatures import Exon, Intron, CDS, UTR
from .misc_features import Promoter
from .utils.genefunctions import reverse_complement

class Transcript(Feature):

    __slots__ = (
        'exons', 'CDSs', 'temp_CDSs', 'temp_UTRs', 'main',
        'miRNAs', 'renamed_exons', 'renamed_utrs', 'polycistronic',
        'promoter', 'introns', 'collapsed_exons', 'collapsed_CDS_segments', 'generated_exons'
    )

    CDSs: dict[str, CDS]
    exons: list[Exon]
    introns: list[Intron]| None
    UTRs: list[UTR]
    temp_CDSs: list[CDS|Feature]|None
    temp_UTRs: list[UTR]|None
    promoter: Promoter | None

    def __init__(self, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)
        self.exons = []
        self.CDSs = {}
        self.temp_CDSs = []
        self.temp_UTRs = []
        self.main = False
        self.miRNAs = []
        self.renamed_exons = False
        self.renamed_utrs = False
        self.polycistronic = "no"
        self.introns = None
        self.promoter = None

        self.collapsed_exons = False
        self.collapsed_CDS_segments = False
        self.generated_exons = False

    def update(self, quiet:bool=False, consider_read_utrs:bool=False, consider_polycistronic:bool=False):
        if self.exons == []:
            self.generate_CDSs(quiet=quiet, consider_read_utrs=True, consider_polycistronic=consider_polycistronic)
            self.generate_exons()
            if len(self.exons) > 1:
                self.exons.sort()
        else:
            if len(self.exons) > 1:
                self.exons.sort()
            self.generate_CDSs(quiet=quiet, consider_read_utrs=consider_read_utrs, consider_polycistronic=consider_polycistronic)

        if self.exons:
            if len(self.exons) > 1:
                self.exons.sort()
            self.start = self.exons[0].start
            self.end = self.exons[-1].end

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

            if self.strand == "+" or self.strand == ".":
                for c in self.CDSs.values():
                    for x, u in enumerate(c.UTRs):
                        if keep_numbering and u.id_number != None:
                            u.id = f"{base_id}{sep}u{u.id_number:0{digits}d}"
                        else:
                            u.id = f"{base_id}{sep}u{x+1:0{digits}d}"

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

                        if u.original_id != u.id:
                            self.renamed_utrs = True
                            u.update_numbering()

    def collapse_exons(self):
        """
        Merges overlapping or directly adjacent exons into single exons.
        """
        
        if len(self.exons) > 1:

            self.exons.sort()
            merged = []
            cur_start = self.exons[0].start
            cur_end = self.exons[0].end
            parents = [self.id]
            for x, e in enumerate(self.exons[1:]):
                if e.start <= cur_end + 1:
                    if e.end > cur_end:
                        
                        cur_end = e.end
                else:
                    merged.append(Exon("combined", self.exons[x].ch, self.exons[x].source, "exon", self.exons[x].strand, cur_start, cur_end, self.exons[x].score, parents))
                    cur_start = e.start
                    cur_end = e.end

            merged.append(Exon("combined", self.exons[-1].ch, self.exons[-1].source, "exon", self.exons[-1].strand, cur_start, cur_end, self.exons[-1].score, parents))

            if len(merged) < len(self.exons):
                self.collapsed_exons = True
                self.exons = merged
                self.update()

    def collapse_CDS_segments(self):
        """
        Merges overlapping or directly adjacent CDS segments into single segments
        """
        parents = [self.id]
        for cds in self.CDSs.values():

            if len(cds.CDS_segments) > 1:
                cds.CDS_segments.sort()
                merged = []
                cur_start = cds.CDS_segments[0].start
                cur_end = cds.CDS_segments[0].end
                group_segs = [cds.CDS_segments[0]]
                for seg in cds.CDS_segments[1:]:
                    if seg.start <= cur_end + 1:
                        # Overlapping segments (seg.start <= cur_end) always collapse
                        # Directly adjacent segments (seg.start == cur_end + 1) only collapse if phase-compatible
                        prev_seg = group_segs[-1]
                        phase_compatible = True
                        is_adjacent = seg.start == cur_end + 1
                        is_short_overlap = seg.start <= cur_end and (cur_end - seg.start + 1) in (1, 2)
                        if (is_adjacent or is_short_overlap) and prev_seg.phase is not None and seg.phase is not None:
                            overlap = (cur_end - seg.start + 1) if is_short_overlap else 0
                            if cds.strand != "-":
                                prev_lo = (prev_seg.size - prev_seg.phase) % 3
                                phase_compatible = (prev_lo - overlap + seg.phase) % 3 == 0
                            else:
                                seg_lo = (seg.size - seg.phase) % 3
                                phase_compatible = (seg_lo - overlap + prev_seg.phase) % 3 == 0

                        if phase_compatible:
                            if seg.end > cur_end:
                                cur_end = seg.end
                            group_segs.append(seg)
                            continue

                    rep = group_segs[0] if cds.strand != "-" else group_segs[-1]
                    merged.append(Feature(cds.id, rep.ch, rep.source, "CDS", rep.strand, cur_start, cur_end, rep.score, parents, phase=rep.phase))
                    cur_start = seg.start
                    cur_end = seg.end
                    group_segs = [seg]

                rep = group_segs[0] if cds.strand != "-" else group_segs[-1]
                merged.append(Feature(cds.id, rep.ch, rep.source, "CDS", rep.strand, cur_start, cur_end, rep.score, parents, phase=rep.phase))

                if len(merged) < len(cds.CDS_segments):
                    cds.CDS_segments = merged
                    self.collapsed_CDS_segments = True
                    cds.update_phase(override=True)
                    cds.update_frame()
                    cds.update()
                    
        if self.collapsed_CDS_segments:
            self.update()

    def clear_UTRs(self):
        for c in self.CDSs.values():
            c.clear_UTRs()
        self.temp_UTRs = None

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

    def generate_CDSs(self, quiet:bool=False, consider_polycistronic:bool=False, consider_read_utrs:bool=False):
        """
        Creates CDS objects are based on the CDS segments of a transcript.
        Assumptions (when considering polycistronic transcripts):
        1 - if more than 1 CDS segment shares an ID that means that we can rely
            on the IDs to generate the different CDSs in a polycistronic gene
        2 - polycistronic genes with non-overlapping CDSs are unlikely to exist
            and if they do exist they are probably annotated as different ID
            transcripts
        """

        if self.temp_CDSs:
            if len(self.temp_CDSs) > 1:
                self.temp_CDSs.sort()

            parents = [self.id]
        
            self.coding = True
            if consider_polycistronic:
                grant_ids = True
                for segment in self.temp_CDSs:
                    if segment.id != "":
                        grant_ids = False

                if grant_ids:
                    for n, _ in enumerate(self.temp_CDSs):
                        self.temp_CDSs[n].id = f"{self.id}_CDS1"

                more_than_1_CDS = False
                more_than_1_segment_with_same_ID = False
                more_than_1_segment_with_different_ID = False
                self.temp_CDSs.sort()
                # more than 1 CDS is determined by overlaps
                if len(self.temp_CDSs) > 1:
                    seen_ids = {self.temp_CDSs[0].id}
                    for sn in range(1, len(self.temp_CDSs)):
                        prev = self.temp_CDSs[sn - 1]
                        curr = self.temp_CDSs[sn]
                        if curr.start <= prev.end:
                            more_than_1_CDS = True
                        seg_id = curr.id
                        if seg_id in seen_ids:
                            more_than_1_segment_with_same_ID = True
                        else:
                            more_than_1_segment_with_different_ID = True
                            seen_ids.add(seg_id)
                    if len(seen_ids) > 1:
                        more_than_1_segment_with_different_ID = True
                    if len(seen_ids) < len(self.temp_CDSs):
                        more_than_1_segment_with_same_ID = True
                if not more_than_1_CDS and self.temp_CDSs != []:
                    if self.strand == "+":
                        temp_id = self.temp_CDSs[0].id
                    else:
                        temp_id = self.temp_CDSs[-1].id
                    self.CDSs[temp_id] = CDS(self.temp_CDSs, temp_id, 
                                            self.temp_CDSs[0].ch, 
                                            self.temp_CDSs[0].source, 
                                            self.temp_CDSs[0].feature,
                                            self.temp_CDSs[0].strand, 
                                            self.temp_CDSs[0].start,
                                            self.temp_CDSs[-1].end,
                                            self.temp_CDSs[0].score,
                                            parents)
                    if more_than_1_segment_with_same_ID and more_than_1_segment_with_different_ID:
                        if not quiet:
                            print(f"Warning: Transcript {self.id} may be "
                                "polycistronic although CDS segments were all "
                                "combined into the same CDS since the most likely "
                                "scenario is that some mistake has been made in the"
                                " gff, please check")   
                        self.polycistronic = "maybe"
                elif more_than_1_CDS and more_than_1_segment_with_different_ID:
                    CDS_temp = {}
                    for c in self.temp_CDSs:
                        if c.id not in CDS_temp:
                            CDS_temp[c.id] = [c]
                        else:
                            CDS_temp[c.id].append(c)
                    for c_id, segments in CDS_temp.items():
                        self.CDSs[c_id] = CDS(segments, c_id, segments[0].ch,
                                        segments[0].source, segments[0].feature,
                                        segments[0].strand, segments[0].start,
                                        segments[-1].end, segments[0].score,
                                        parents)
                    if not quiet:
                        print(f"Warning: Transcript {self.id} is likely to be "
                            "polycistronic since CDS segments overlap and they "
                            "have different IDs, the CDS segments have been "
                            "separated into their corresponding CDS ids, however, "
                            "please check that it truly is a polycistronic gene "
                            "and not a gff mistake")
                    self.polycistronic = "yes" 
                elif more_than_1_CDS:
                    if not quiet:
                        print(f"Error: Transcript {self.id} is likely to have a "
                            "problem in the annotation of CDS segments (it could "
                            "also be a consequence of liftoff) as the segments "
                            "overlap but they share the same id, please fix the gff.")
                    self.polycistronic = "maybe"
                

            else:
                if self.strand == "+":
                    temp_id = self.temp_CDSs[0].id
                else:
                    temp_id = self.temp_CDSs[-1].id
                self.CDSs[temp_id] = CDS(self.temp_CDSs, temp_id, 
                                        self.temp_CDSs[0].ch, 
                                        self.temp_CDSs[0].source, 
                                        self.temp_CDSs[0].feature,
                                        self.temp_CDSs[0].strand, 
                                        self.temp_CDSs[0].start,
                                        self.temp_CDSs[-1].end,
                                        self.temp_CDSs[0].score,
                                        parents)

            self.temp_CDSs = None

        self.determine_main_CDS()

        if not consider_read_utrs:
            self.generate_UTRs()
        else:
            # consider_read_utrs is True
            if not self.temp_UTRs or self.polycistronic == "yes":
                self.generate_UTRs()
            else:
                self.assign_UTRs()

        for c in self.CDSs.values():
            c.update()

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
                c.UTRs = self.temp_UTRs
            self.temp_UTRs = None

    def generate_UTRs(self):
        if self.temp_UTRs:
            self.temp_UTRs = None
        for c in self.CDSs.values():
            c.UTRs = []
            parents = [self.id]
            for exon in self.exons:
                if c.strand != exon.strand:
                    continue
                if exon.end <= c.CDS_segments[-1].end and exon.start >= c.CDS_segments[0].start:
                    continue
                if exon.end < c.CDS_segments[0].start:
                    c.UTRs.append(UTR("", exon.ch, exon.source, "UTR",
                                      exon.strand, exon.start, exon.end,
                                      exon.score, parents))
                elif exon.start < c.CDS_segments[0].start:
                    c.UTRs.append(UTR("", exon.ch, exon.source, "UTR",
                                      exon.strand, exon.start, c.CDS_segments[0].start-1,
                                      exon.score, parents))
                if exon.start > c.CDS_segments[-1].end:
                    c.UTRs.append(UTR("", exon.ch, exon.source, "UTR",
                                      exon.strand, exon.start, exon.end,
                                      exon.score, parents))
                elif exon.end > c.CDS_segments[-1].end:
                    c.UTRs.append(UTR("", exon.ch, exon.source, "UTR",
                                      exon.strand, c.CDS_segments[-1].end+1, exon.end,
                                      exon.score, parents))
            c.UTRs.sort()
            if c.strand == "+" or c.strand == ".":
                for n, u in enumerate(c.UTRs):
                    u.id = f"{c.id}_u{n+1}"
                    u.parents = parents
            elif c.strand == "-":
                counter = len(c.UTRs)
                for n, u in enumerate(c.UTRs):
                    u.id = f"{c.id}_u{counter}"      
                    u.parents = parents 
                    counter -= 1
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
                            u.feature = "three_prime_UTR"
                elif c.strand == "-":
                    for u in c.UTRs:
                        if u.end > c.end:
                            u.prime = "5'"
                            u.feature = "five_prime_UTR"
                        else:
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
        self.introns = []
        counter = 0
        parents = [self.id]
        for n, exon in enumerate(self.exons):
            counter += 1
            if n == (len(self.exons) - 1):
                continue
            self.introns.append(Intron(f"{self.id}_intron_{counter}", self.ch,
                                       self.source, "intron", self.strand,
                                       exon.end + 1, self.exons[n+1].start - 1,
                                       self.score, parents))
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
            transcript_seqs = ["", ""]
            for exon in self.exons:
                transcript_seqs[0] += exon.seq

            for exon in reversed(self.exons):
                transcript_seqs[1] += reverse_complement(exon.seq)

            return transcript_seqs

    @property
    def hard_seqs(self) -> list[str]:
        if not self._ACTIVE_HARD_GENOME:
            raise ValueError("No hard masked genome loaded and you are trying to access the hard masked sequence. Load your hard masked genome together with your annotation.")
        else:
            transcript_seqs = ["", ""]
            for exon in self.exons:
                transcript_seqs[0] += exon.hard_seq

            for exon in reversed(self.exons):
                transcript_seqs[1] += reverse_complement(exon.hard_seq)

            return transcript_seqs
    
    @property
    def size(self):
        size = 0
        for exon in self.exons:
            size += exon.size
        return size

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