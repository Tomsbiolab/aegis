from .feature import Feature
from .misc_features import Protein
from .genefunctions import reverse_complement

class CDS(Feature):
    def __init__(self, CDS_segments:list, feature_id:str, 
                 ch:str, source:str, feature:str, strand:str, start:int, 
                 end:int, score:str, phase:str, attributes:str):
        super().__init__(feature_id, ch, source, feature, strand, start, end,
                         score, phase, attributes)    
        self.main = False
        self.CDS_segments = CDS_segments
        #self.kozak = ""
        #self.codon_usage = {}
        self.five_prime_UTR_seq = ""
        self.three_prime_UTR_seq = ""
        self.size = 0
        self.full_UTR_exons = 0
        self.frame = "."
        self.protein = None
        for segment in self.CDS_segments:
            self.size += segment.size

    def update(self):
        self.update_size()
        self.update_phase()
        self.update_frame()

    def update_size(self):
        self.size = 0
        for segment in self.CDS_segments:
            self.size += segment.size

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
                frame = (cs.start + cs.phase) % 3
                if frame == 0:
                    frame = 3
                cs.frame = frame

        if self.strand == "-":
            for cs in reversed(self.CDS_segments):
                frame = (cs.start + cs.phase) % 3
                if frame == 0:
                    frame = 3
                frame = 7 - frame
                cs.frame = frame


    def rename(self, base_id, base_gene_id, count, sep:str="_", digits:int=3, keep_numbering:bool=False, keep_ids_with_base_id_contained:bool=False, cds_segment_ids:bool=False):

        rename = False
        rename_cs = False

        if keep_ids_with_base_id_contained:
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
                if self.main:
                    self.id = f"{base_id}{sep}CDS{self.id_number:0{digits}d}"
                else:
                    self.id = f"{base_id}{sep}CDS{self.id_number:0{digits}d}"
            else:
                if self.main:
                    self.id = f"{base_id}{sep}CDS{1:0{digits}d}"
                else:
                    self.id = f"{base_id}{sep}CDS{count:0{digits}d}"

            if self.original_id != self.id:
                self.renamed = True
                self.update_numbering()
                self.generate_protein()

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
        self.five_prime_UTR_seq = ""
        self.three_prime_UTR_seq = ""
        self.full_UTR_exons = 0
        del self.UTRs

    def generate_sequence(self, genome:object, low_memory:bool=False):
        self.seq = ""
        for segment in self.CDS_segments:
            segment.generate_sequence(genome)
        if self.strand == "+":
            for segment in self.CDS_segments:
                self.seq += segment.seq
        elif self.strand == "-":
            for segment in reversed(self.CDS_segments):
                self.seq += segment.seq
        elif self.strand == ".":
            self.seqs = ["", ""]
            for segment in self.CDS_segments:
                self.seqs[0] += segment.seqs[0]
            for segment in reversed(self.CDS_segments):
                self.seqs[1] += segment.seqs[1]

        if low_memory:
            for segment in self.CDS_segments:
                segment.clear_sequence()
        else:
            if hasattr(self, "UTRs"):
                for u in self.UTRs:
                    u.generate_sequence(genome)
                if self.strand == "+":
                    for u in self.UTRs:
                        if u.prime == "5'":
                            self.five_prime_UTR_seq += u.seq
                        else:
                            self.three_prime_UTR_seq += u.seq
                elif self.strand == "-":
                    for u in reversed(self.UTRs):
                        if u.prime == "5'":
                            self.five_prime_UTR_seq += u.seq
                        else:
                            self.three_prime_UTR_seq += u.seq
        self.generate_protein(low_memory=low_memory)

    def generate_hard_sequence(self, hard_masked_genome:object, low_memory:bool=False):
        self.hard_seq = ""
        for segment in self.CDS_segments:
            segment.generate_hard_sequence(hard_masked_genome)
        if self.strand == "+":
            for segment in self.CDS_segments:
                self.hard_seq += segment.hard_seq
        elif self.strand == "-":
            for segment in reversed(self.CDS_segments):
                self.hard_seq += segment.hard_seq
        elif self.strand == ".":
            self.hard_seqs = ["", ""]
            for segment in self.CDS_segments:
                self.hard_seqs[0] += segment.hard_seqs[0]
            for segment in reversed(self.CDS_segments):
                self.hard_seqs[1] += segment.hard_seqs[1]

        if low_memory:
            for segment in self.CDS_segments:
                segment.clear_sequence(just_hard=False)
        else:
            if hasattr(self, "UTRs"):
                for u in self.UTRs:
                    u.generate_hard_sequence(hard_masked_genome)

    def clear_sequence(self, just_hard=False, keep_proteins:bool=False):
        self.hard_seq = ""
        for segment in self.CDS_segments:
            segment.clear_sequence(just_hard=just_hard)
        if hasattr(self, "UTRs"):
            for u in self.UTRs:    
                u.clear_sequence(just_hard=just_hard)
        if not just_hard:
            self.seq = ""
            self.five_prime_UTR_seq = ""
            self.three_prime_UTR_seq = ""
            if not keep_proteins:
                self.protein = None

    def generate_protein(self, readthrough:str="both", low_memory:bool=False):
        self.protein = Protein(f"{self.id}.prot", self.seq, self.ch, readthrough)
        if low_memory:
            self.seq = ""

    def equal_segments(self, other):
        self.CDS_segments.sort()
        other.CDS_segments.sort()
        same = True
        if len(self.CDS_segments) == len(other.CDS_segments):
            for n, segment in enumerate(self.CDS_segments):
                if not segment.almost_equal(other.CDS_segments[n]):
                    same = False
        else:
            same = False
        
        return same

class Exon(Feature):
    def __init__(self, feature_id:str, ch:str, source:str, feature:str,
                 strand:str, start:int, end:int, score:str, phase:str, 
                 attributes:str):
        super().__init__(feature_id, ch, source, feature, strand, start, end,
                         score, phase, attributes)
        self.coding = False

class UTR(Feature):
    def __init__(self, feature_id:str, ch:str, source:str, feature:str,
                 strand:str, start:int, end:int, score:str, phase:str, 
                 attributes:str):
        super().__init__(feature_id, ch, source, feature, strand, start, end,
                         score, phase, attributes)
        self.prime = "3'"

class Intron(Feature):
    canonical_seqs = ["GT-AG", "GC-AG", "AT-AC"]
    def __init__(self, feature_id:str, ch:str, source:str, feature:str,
                 strand:str, start:int, end:int, score:str, phase:str, 
                 attributes:str):
        super().__init__(feature_id, ch, source, feature, strand, start, end,
                         score, phase, attributes)
        self.intra_coding = False
        self.boundary = ""
        self.canonical = False
        self.splice_site_donor = ""
        self.splice_site_acceptor = ""

    def generate_sequence(self, genome):
        if self.strand == "+":
            self.seq = genome.scaffolds[self.ch].seq[self.start-1:self.end]
        elif self.strand == "-":
            self.seq = reverse_complement(genome.scaffolds[self.ch].seq[self.start-1:self.end])
        elif self.strand == ".":
            self.seqs = (genome.scaffolds[self.ch].seq[self.start-1:self.end], reverse_complement(genome.scaffolds[self.ch].seq[self.start-1:self.end]))
        self.splice_site_donor = self.seq[0:2]
        self.splice_site_acceptor = self.seq[-2:]
        self.boundary = f"{self.splice_site_donor}-{self.splice_site_acceptor}"
        if self.boundary in Intron.canonical_seqs:
            self.canonical = True
        else:
            self.canonical = False

    def clear_sequence(self, just_hard=False):
        self.hard_seq = ""
        if not just_hard:
            self.seq = ""