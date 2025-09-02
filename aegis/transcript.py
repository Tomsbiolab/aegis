from .feature import Feature
from .subfeatures import Exon, Intron, CDS, UTR
from .misc_features import Promoter
from .genefunctions import find_ORFs, longest_ORF, translate, overlap


class Transcript(Feature):
    def __init__(self, feature_id:str, ch:str, source:str, 
                 feature:str, strand:str, start:int, end:int, score:str, 
                 phase:str, attributes:str):
        super().__init__(feature_id, ch, source, feature, strand, start, end,
                         score, phase, attributes)
        self.exons = []
        self.CDSs = {}
        self.temp_CDSs = []
        self.temp_UTRs = []
        self.coding = False
        self.main = False
        self.miRNAs = []
        self.overlaps = {"self" : [], "other" : []}
        self.renamed_exons = False
        self.renamed_utrs = False
    
    def update_size(self):
        self.size = 0
        for exon in self.exons:
            self.size += exon.size

    def update(self, quiet:bool=False):
        if self.feature == "lncRNA" or not self.coding:
            self.CDSs = {}
            self.coding = False

        self.generate_CDSs(quiet=quiet)
        if self.exons == []:
            self.generate_exons()
        self.exons.sort()
        self.generate_UTRs()
        self.update_UTRs()
        self.update_size()
        self.generate_introns()
        self.exon_update()
        c_start = None
        c_end = None
        for i, c in enumerate(self.CDSs.values()):
            if i == 0:
                c_start = c.start
                c_end = c.end
            else:
                if c.start < c_start:
                    c_start = c.start
                if c.end > c_end:
                    c_end = c.end
        if c_start != None:
            if self.strand == "+":
                for e in self.exons:
                    if e.end > c_start and e.start < c_end:
                        e.coding = True
                    
            elif self.strand == "-":
                for e in self.exons:
                    if e.start < c_end and e.end > c_start:
                        e.coding = True

    def rename(self, base_id, count, sep:str="_", digits:int=3, keep_numbering:bool=False, keep_ids_with_base_id_contained:bool=False):

        rename = False

        if keep_ids_with_base_id_contained:
            if base_id not in self.id:
                rename = True
        else:
            rename = True

        if rename:

            if keep_numbering and self.id_number != None:
                if self.main:
                    self.id = f"{base_id}{sep}t{self.id_number:0{digits}d}"
                else:
                    self.id = f"{base_id}{sep}t{self.id_number:0{digits}d}"
            else:
                if self.main:
                    self.id = f"{base_id}{sep}t{1:0{digits}d}"
                else:
                    self.id = f"{base_id}{sep}t{count:0{digits}d}"

            if self.original_id != self.id:
                self.renamed = True
                self.update_numbering()

    def rename_exons(self, count, base_id, sep:str="_", digits:int=3, keep_numbering:bool=False, keep_ids_with_base_id_contained:bool=False, rev:bool=False):

        rename = False

        if keep_ids_with_base_id_contained:
            for e in self.exons:
                if base_id not in e.id:
                    rename = True
        else:
            rename = True

        if rename:
            for x, e in enumerate(self.exons):
                if rev and x != 0:
                    count -= 1
                else:
                    count += 1
                if keep_numbering and e.id_number != None:
                    e.id = f"{base_id}{sep}e{e.id_number:0{digits}d}"
                else:
                    e.id = f"{base_id}{sep}e{count:0{digits}d}"

                if e.original_id != e.id:
                    self.renamed_exons = True
                    e.update_numbering()

    def rename_utrs(self, count, base_id, sep:str="_", digits:int=3, keep_numbering:bool=False, keep_ids_with_base_id_contained:bool=False, rev:bool=False):

        rename = False

        if keep_ids_with_base_id_contained:
            for c in self.CDSs.values():
                for u in c.UTRs:
                    if base_id not in u.id:
                        rename = True
        else:
            rename = True

        if rename:
            for x, c in enumerate(self.CDSs.values()):
                for j, u in enumerate(c.UTRs):
                    if rev and x != 0 and j != 0:
                        count -= 1
                    else:
                        count += 1
                    if keep_numbering and u.id_number != None:
                        u.id = f"{base_id}{sep}u{u.id_number:0{digits}d}"
                    else:
                        u.id = f"{base_id}{sep}u{count:0{digits}d}"

                    if u.original_id != u.id:
                        self.renamed_utrs = True
                        u.update_numbering()

    def exon_update(self):
        CDS_size = 0
        c_start = None
        c_end = None
        for c in self.CDSs.values():
            if c.main:
                CDS_size = c.size
        if CDS_size != 0:
            self.coding_ratio = round((CDS_size / self.size), 2)
        else:
            self.coding_ratio = 0
        if c_start != None:
            if self.strand == "+":
                for e in self.exons:
                    if e.end > c_start and e.start < c_end:
                        e.coding = True
            elif self.strand == "-":
                for e in self.exons:
                    if e.start < c_end and e.end > c_start:
                        e.coding = True                

    def clear_UTRs(self):
        for c in self.CDSs.values():
            c.clear_UTRs()
        self.temp_UTRs = []
        self.exons = []
        self.update()

    def generate_promoter(self, promoter_size, ch_size, promoter_type:str = "standard"):
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


            self.promoter = Promoter(promoter_type, prom_id, self.ch, self.source,
                                    self.feature, self.strand, temp_start,
                                    temp_end, self.score, ".",
                                    self.attributes)

    def generate_best_protein(self, genome:object, must_have_stop:bool=True):
        if (self.strand == "+") or (self.strand == "-"):
            self.protein_start, self.protein_end_stop, self.protein_early_stop, self.protein_nucleotide_surplus, self.protein_gaps, self.protein_seq, self.coding_start, self.coding_end = translate(self.seq, "none", must_have_stop=must_have_stop)
        elif self.strand == ".":
            plus_orfs = find_ORFs(self.seqs[0], must_have_stop)
            neg_orfs = find_ORFs(self.seqs[1], must_have_stop)
            plus_long_orf, _, _ = longest_ORF(plus_orfs)
            neg_long_orf, _, _ = longest_ORF(neg_orfs)

            if plus_long_orf != "" or neg_long_orf != "":
                if len(plus_long_orf) >= len(neg_long_orf):
                    self.strand = "+"
                    for e in self.exons:
                        e.strand = "+"
                        e.generate_sequence(genome)
                    self.generate_sequence(genome)
                    self.generate_best_protein(genome, must_have_stop)

                else:
                    self.strand = "-"
                    for e in self.exons:
                        e.strand = "-"
                        e.generate_sequence(genome)
                    self.generate_sequence(genome)
                    self.generate_best_protein(genome, must_have_stop)

    def generate_CDSs_based_on_ORF(self, low_memory:bool=True):
        if not hasattr(self, "temp_CDSs"):
            self.temp_CDSs = []
        if not hasattr(self, "temp_UTRs"):
            self.temp_UTRs = []
        if self.temp_CDSs == []:
            if self.protein_seq != "":
                start_exon = ""
                end_exon = ""
                surplus_start = ""
                surplus_end = ""
                if self.strand == "+":
                    temp_size = 0
                    for index, e in enumerate(self.exons):
                        temp_size += e.size
                        if temp_size > self.coding_start:
                            surplus_start =  self.coding_start - (temp_size-e.size)
                            start_exon = index
                            break
                    temp_size = 0
                    for index, e in enumerate(self.exons):
                        temp_size += e.size
                        if temp_size > self.coding_end:
                            surplus_end = self.coding_end - (temp_size-e.size)
                            end_exon = index
                            break
                    
                    for index, e in enumerate(self.exons):

                        if (index == start_exon) and (index == end_exon):
                            self.temp_CDSs.append(Feature(f"{self.id}_CDS1", e.ch, e.source, "CDS", e.strand, 
                                                        e.start+surplus_start, e.start+surplus_end,
                                                        e.score, e.phase, 
                                                        f"ID={self.id}_CDS1;Parent={self.id}"))
                        elif index == start_exon:
                            self.temp_CDSs.append(Feature(f"{self.id}_CDS1", e.ch, e.source, "CDS", e.strand, 
                                                        e.start+surplus_start, e.end, e.score, e.phase, 
                                                        f"ID={self.id}_CDS1;Parent={self.id}"))
                        elif (index > start_exon) and (index < end_exon):
                            self.temp_CDSs.append(Feature(f"{self.id}_CDS1", e.ch, e.source, "CDS", e.strand, e.start, e.end, e.score, e.phase, 
                                                        f"ID={self.id}_CDS1;Parent={self.id}"))
                        elif index == end_exon:
                            self.temp_CDSs.append(Feature(f"{self.id}_CDS1", e.ch, e.source, "CDS", e.strand, 
                                                        e.start, e.start+surplus_end, e.score, e.phase, 
                                                        f"ID={self.id}_CDS1;Parent={self.id}"))

                elif self.strand == "-":
                    temp_size = 0
                    for index, e in enumerate(reversed(self.exons)):
                        temp_size += e.size
                        if temp_size > self.coding_start:
                            surplus_start = self.coding_start - (temp_size-e.size)
                            start_exon = index
                            break
                    temp_size = 0
                    for index, e in enumerate(reversed(self.exons)):
                        temp_size += e.size
                        if temp_size > self.coding_end:
                            surplus_end = self.coding_end - (temp_size-e.size)
                            end_exon = index
                            break
                    
                    for index, e in enumerate(reversed(self.exons)):

                        if (index == start_exon) and (index == end_exon):
                            self.temp_CDSs.append(Feature(f"{self.id}_CDS1", e.ch, e.source, "CDS", e.strand, 
                                                        e.end-surplus_end, e.end-surplus_start,
                                                        e.score, e.phase, 
                                                        f"ID={self.id}_CDS1;Parent={self.id}"))
                        elif index == start_exon:
                            self.temp_CDSs.append(Feature(f"{self.id}_CDS1", e.ch, e.source, "CDS", e.strand, 
                                                        e.start, e.end-surplus_start, e.score, e.phase, 
                                                        f"ID={self.id}_CDS1;Parent={self.id}"))
                        elif (index > start_exon) and (index < end_exon):
                            self.temp_CDSs.append(Feature(f"{self.id}_CDS1", e.ch, e.source, "CDS", e.strand, e.start, e.end, e.score, e.phase, 
                                                        f"ID={self.id}_CDS1;Parent={self.id}"))
                        elif index == end_exon:
                            self.temp_CDSs.append(Feature(f"{self.id}_CDS1", e.ch, e.source, "CDS", e.strand, 
                                                        e.end-surplus_end, e.end, e.score, e.phase, 
                                                        f"ID={self.id}_CDS1;Parent={self.id}"))

                elif self.strand == ".":
                    pass


                self.generate_CDSs()
                if low_memory:
                    self.clear_sequence()
        else:
            print(f"CDS segments exist already for {self.id}")

    def almost_equal(self, other):
        almost_equal = True
        if len(self.exons) != len(other.exons):
            almost_equal = False
        else:
            for n, exon in enumerate(self.exons):
                if exon.start != other.exons[n].start or exon.end != other.exons[n].end:
                    almost_equal = False
                    break
        return almost_equal

    def generate_CDSs(self, quiet:bool=False):
        """
        Creates CDS objects are based on the CDS segments of a transcript.
        Assumptions:
        1 - if more than 1 CDS segment shares an ID that means that we can rely
            on the IDs to generate the different CDSs in a polycistronic gene
        2 - polycistronic genes with non-overlapping CDSs are unlikely to exist
            and if they do exist they are probably annotated as different ID
            transcripts
        """
        # This is one of three options = ["yes", "no", "maybe"]
        self.polycistronic = "no"
        if hasattr(self, "temp_CDSs") and hasattr(self, "temp_UTRs"):
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
            if len(self.temp_CDSs) > 0:
                self.coding = True
            self.temp_CDSs.sort()
            # more than 1 CDS is determined by overlaps
            if len(self.temp_CDSs) > 1:
                for sn1, segment1 in enumerate(self.temp_CDSs):
                    for sn2, segment2 in enumerate(self.temp_CDSs):
                        if sn1 == sn2:
                            continue
                        # not interested in overlap_bp
                        overlapping, _ = overlap(segment1, segment2)
                        if overlapping:
                            more_than_1_CDS = True
                        if segment1.id == segment2.id:
                            more_than_1_segment_with_same_ID = True
                        else:
                            more_than_1_segment_with_different_ID = True
            if not more_than_1_CDS and self.temp_CDSs != []:
                # CDS object with segments as a property
                if self.strand == "+":
                    temp_id = self.temp_CDSs[0].id
                else:
                    temp_id = self.temp_CDSs[-1].id
                #these two lines resulted in repeated CDS ids
                #pattern = r"(?i)(CDS|cds)(\d+)$"
                #temp_id = sub(pattern, r"\g<1>1", temp_id)
                self.CDSs[temp_id] = CDS(self.temp_CDSs.copy(), temp_id, 
                                         self.temp_CDSs[0].ch, 
                                         self.temp_CDSs[0].source, 
                                         self.temp_CDSs[0].feature,
                                         self.temp_CDSs[0].strand, 
                                         self.temp_CDSs[0].start,
                                         self.temp_CDSs[-1].end,
                                         self.temp_CDSs[0].score,
                                         ".", self.temp_CDSs[0].attributes)
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
                    self.CDSs[c_id] = CDS(segments.copy(), c_id, segments[0].ch,
                                       segments[0].source, segments[0].feature,
                                       segments[0].strand, segments[0].start,
                                       segments[-1].end, segments[0].score,
                                       ".", segments[0].attributes)
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
                        "overlap but they share the same id, please fix the "
                        "gff.")
                self.polycistronic = "maybe"
            del self.temp_CDSs
        self.determine_main_CDS()
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
        if hasattr(self, "temp_UTRs"):
            if self.temp_UTRs == []:
                self.generate_UTRs()
            else:
                self.temp_UTRs.sort()
                if self.polycistronic == "no" or self.polycistronic == "maybe":
                    if len(self.CDSs) == 1:
                        for c in self.CDSs.values():
                            c.UTRs = self.temp_UTRs.copy()
                        del self.temp_UTRs
                    elif len(self.CDSs) > 1:
                        print(f"Error in code: More than 1 CDS object for {self.id} transcript")
                elif self.exons == []:
                    print("Error: Unable to assign UTRs as GFF has no exons "
                        f"and several CDSs exist for this {self.id} transcript")
                elif self.polycistronic == "yes":
                    self.generate_UTRs()
    
    def generate_UTRs(self):
        for c in self.CDSs.values():
            c.UTRs = []
            for exon in self.exons:
                if c.strand != exon.strand:
                    continue
                if exon.end <= c.CDS_segments[-1].end and exon.start >= c.CDS_segments[0].start:
                    pass
                if exon.end < c.CDS_segments[0].start:
                    c.UTRs.append(UTR("", exon.ch, exon.source, "UTR",
                                      exon.strand, exon.start, exon.end,
                                      exon.score, ".", ""))
                elif exon.start < c.CDS_segments[0].start:
                    c.UTRs.append(UTR("", exon.ch, exon.source, "UTR",
                                      exon.strand, exon.start, c.CDS_segments[0].start-1,
                                      exon.score, ".", ""))
                if exon.start > c.CDS_segments[-1].end:
                    c.UTRs.append(UTR("", exon.ch, exon.source, "UTR",
                                      exon.strand, exon.start, exon.end,
                                      exon.score, ".", ""))
                elif exon.end > c.CDS_segments[-1].end:
                    c.UTRs.append(UTR("", exon.ch, exon.source, "UTR",
                                      exon.strand, c.CDS_segments[-1].end+1, exon.end,
                                      exon.score, ".", ""))
            c.UTRs.sort()
            if c.strand == "+":
                for n, u in enumerate(c.UTRs):
                    u.id = f"{c.id}_u{n+1}"
                    u.attributes = f"ID={u.id};Parent={self.id}"
                    u.parents = [self.id]
            elif c.strand == "-":
                counter = len(c.UTRs)
                for n, u in enumerate(c.UTRs):
                    u.id = f"{c.id}_u{counter}"
                    u.attributes = f"ID={u.id};Parent={self.id}"      
                    u.parents = [self.id]  
                    counter -= 1
        self.update_UTRs()

    def update_UTRs(self):
        for c in self.CDSs.values():
            if hasattr(c, "UTRs"):
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
        # WARNING! double check function
        temp_fts = []
        for c in self.CDSs.values():
            if c.main:
                for cs in c.CDS_segments:
                    temp_fts.append(cs.copy())
                for u in c.UTRs:
                    temp_fts.append(u.copy())
        # Exons reconstructed from CDS/UTRs
        if temp_fts != []:
            still_overlapping_fts = True
            while still_overlapping_fts:
                features_to_remove = []
                features_to_add = []
                for i, tempft1 in enumerate(temp_fts):
                    for j, tempft2 in enumerate(temp_fts):
                        if i != j:
                            interval1 = tempft1.size
                            interval2 = tempft2.size
                            small = min(tempft1.start, tempft1.end, 
                                        tempft2.start, tempft2.end)
                            large = max(tempft1.start, tempft1.end,
                                        tempft2.start, tempft2.end)
                            overlap_bp = ((interval1 + interval2)
                                            - ((large - small) + 1))
                            # >= is crucial to combine contiguous features
                            if overlap_bp >= 0:
                                temp = Exon("combined", self.ch, self.source, 
                                            "exon", self.strand, small, large, 
                                            self.score, ".", "")
                                add = True
                                # this is to avoid adding a same overlap twice
                                for f in features_to_add:
                                    if temp == f:
                                        add = False
                                        break
                                if add:
                                    features_to_add.append(temp)
                                if tempft1 not in features_to_remove:
                                    features_to_remove.append(tempft1)
                                if tempft2 not in features_to_remove:
                                    features_to_remove.append(tempft2)
                for sub_to_add in features_to_add:
                    temp_fts.append(sub_to_add)
                for sub_to_rem in features_to_remove:
                    if sub_to_rem in temp_fts:
                        if sub_to_rem in temp_fts:
                            temp_fts.remove(sub_to_rem)
                if (features_to_remove == []
                    and features_to_add == []):
                    still_overlapping_fts = False
            
            self.exons = temp_fts.copy()
            self.exons.sort()

            if self.strand == "+":
                for n, e in enumerate(self.exons):
                    e.feature = "exon"
                    e.id = f"{self.id}_e{n+1}"
                    e.attributes = f"ID={e.id};Parent={self.id}"
                    e.misc_attributes = ""
                    e.parents = [self.id]
                    e.phase = "."
            elif self.strand == "-":
                counter = len(self.exons)
                for n, e in enumerate(self.exons):
                    e.id = f"{self.id}_e{counter}"
                    e.attributes = f"ID={e.id};Parent={self.id}"      
                    e.parents = [self.id]  
                    counter -= 1

        # Exons rebuilt from the transcript
        else:
            self.exons = [Exon(f"{self.id}_e1", self.ch, self.source, "exon", self.strand, self.start, self.end, self.score, ".", "")]
            self.exons[0].attributes = (f"ID={self.exons[0].id};Parent={self.id}")
            self.exons[0].parents = [self.id]

    def generate_introns(self):
        self.introns = []
        counter = 0
        for n, exon in enumerate(self.exons):
            counter += 1
            if n == (len(self.exons) - 1):
                continue
            self.introns.append(Intron(f"{self.id}_intron_{counter}", self.ch,
                                       self.source, "intron", self.strand,
                                       exon.end + 1, self.exons[n+1].start - 1,
                                       self.score, ".", f"ID={self.id}_intron_{counter};Parent={self.id}"))
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

    def generate_sequence(self, genome:object, low_memory:bool=False):
        for exon in self.exons:
            exon.generate_sequence(genome)
        if not low_memory:
            for intron in self.introns:
                intron.generate_sequence(genome)
            if hasattr(self, "promoter"):
                self.promoter.generate_sequence(genome)
        self.seq = ""
        if self.strand == "+":
            for segment in self.exons:
                self.seq += segment.seq
        elif self.strand == "-":
            for segment in reversed(self.exons):
                self.seq += segment.seq
        elif self.strand == ".":
            self.seqs = ["", ""]
            for segment in self.exons:
                self.seqs[0] += segment.seqs[0]
            for segment in reversed(self.exons):
                self.seqs[1] += segment.seqs[1]
            
    def generate_hard_sequence(self, hard_masked_genome:object, low_memory:bool=False):
        for exon in self.exons:
            exon.generate_hard_sequence(hard_masked_genome)
        if not low_memory:
            for intron in self.introns:
                intron.generate_hard_sequence(hard_masked_genome)
            if hasattr(self, "promoter"):
                self.promoter.generate_hard_sequence(hard_masked_genome)
        self.hard_seq = ""
        if self.strand == "+":
            for segment in self.exons:
                self.hard_seq += segment.hard_seq
        elif self.strand == "-":
            for segment in reversed(self.exons):
                self.hard_seq += segment.hard_seq
        elif self.strand == ".":
            self.hard_seqs = ["", ""]
            for segment in self.exons:
                self.hard_seqs[0] += segment.hard_seqs[0]
            for segment in reversed(self.exons):
                self.hard_seqs[1] += segment.hard_seqs[1]

    def clear_sequence(self, just_hard=False):
        self.hard_seq = ""
        if hasattr(self, "promoter"):
            self.promoter.hard_seq = ""
        for exon in self.exons:
            exon.clear_sequence(just_hard=just_hard)
        if hasattr(self, "introns"):
            for intron in self.introns:
                intron.clear_sequence(just_hard=just_hard)        

        if not just_hard:
            self.seq = ""
            self.protein_seq = ""
            if hasattr(self, "promoter"):
                self.promoter.seq = ""
