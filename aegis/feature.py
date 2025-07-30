from .genefunctions import count_occurrences, reverse_complement
import copy
import re

class Feature():
    """
    Parent class including basic gff feature properties, that can then be
    inherited by Gene, Transcript etc. -> this facilitates code reading and
    editing. Or used on its own for CDS segments.
    """
    # These attributes cannot be mistaken by misc attributes or any other
    attributes_to_ignore_when_reading_gff = ["id", "reliable_score", "remove", "rescue", "blasts", "gene_masked_fraction", "transcript_masked_fraction", "cds_masked_fraction", "gene_gc_content", "transcript_gc_content", "cds_gc_content", "intron_nested", "intron_nested_fully_contained", "intron_nested_single", "intron_utr_nested", "pseudogene", "transposable", "alternative_transcript_rescue", "cds_orientated_overlaps"]
    def __init__(self, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, phase:str, attributes:str):
        
        self.id = feature_id
        self.original_id = feature_id
        self.ch = ch
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = phase
        self.frame = "."
        self.attributes = attributes
        self.gtf_attributes = ""
        self.size = (self.end - self.start) + 1
        self.seq = ""
        self.hard_seq = ""
        self.parents = []
        self.names = []
        self.symbols = []
        self.descriptors = []
        self.processes = []
        self.synonyms = []
        self.gc_content = 0
        self.aliases = []
        self.blast_hits = []
        self.renamed = False

        self.id_number = None
        self.original_id_number = None

        attributes_l = attributes.split(";")
        if attributes_l == [""]:
            attributes_l = []
        misc_attributes = []

        for a in attributes_l:
            a_label = a.split("=")[0].lower()
            if len(a.split("=")) > 1:
                a_value = a.split("=")[1]
            else:
                a_value = ""
            if "parent" == a_label:
                parents = a_value.split(",")
                for p in parents:
                    p = p.strip()
                    if p not in self.parents and p != "":
                        self.parents.append(p)
            elif "alias" == a_label:
                aliases = a_value.split(",")
                for alias in aliases:
                    if alias not in self.aliases and alias != "":
                        alias = alias.strip()
                        self.aliases.append(alias)
            elif a_label in Feature.attributes_to_ignore_when_reading_gff:
                continue
            elif "symbol" == a_label:
                symbols = a_value.split(",")
                for symbol in symbols:
                    if symbol not in self.symbols and symbol != "":
                        symbol = symbol.strip()
                        self.symbols.append(symbol)
            elif "name" in a_label:
                names = a_value.split(",")
                for name in names:
                    if name not in self.names and name != "":
                        name = name.strip()
                        self.names.append(name)
            elif a.strip() != "":
                misc_attributes.append(a.strip())
        self.update_numbering(original=True)

        self.misc_attributes = ";".join(misc_attributes)
        self.masked_fraction = 0

    def update_numbering(self, original:bool=False):

        match = re.search(r'(\d+)$', self.id)
        if match:
            if original:
                self.original_id_number = int(match.group(1))
            self.id_number = int(match.group(1))

    def update_size(self):
        self.size = (self.end - self.start) + 1

    def calculate_masking(self):
        self.masked_fraction = round(((count_occurrences(self.hard_seq, "X") + (count_occurrences(self.hard_seq, "N"))) / self.size), 2)

    def generate_sequence(self, genome:object):
        if self.start != "NA" and self.end != "NA":
            if self.strand == "+":
                self.seq = genome.scaffolds[self.ch].seq[self.start-1:self.end]
            elif self.strand == "-":
                self.seq = reverse_complement(genome.scaffolds[self.ch].seq[self.start-1:self.end])
            elif self.strand == ".":
                self.seqs = (genome.scaffolds[self.ch].seq[self.start-1:self.end], reverse_complement(genome.scaffolds[self.ch].seq[self.start-1:self.end]))

    def clear_sequence(self, just_hard=False):
        self.hard_seq = ""
        if hasattr(self, "hard_seqs"):
            del self.hard_seqs
        if not just_hard:
            self.seq = ""
            if hasattr(self, "seqs"):
                del self.seqs

    def generate_hard_sequence(self, hard_masked_genome:object):
        if self.start != "NA" and self.end != "NA":
            if self.strand == "+":
                self.hard_seq = hard_masked_genome.scaffolds[self.ch].seq[self.start-1:self.end]
            elif self.strand == "-":
                self.hard_seq = reverse_complement(hard_masked_genome.scaffolds[self.ch].seq[self.start-1:self.end])
            elif self.strand == ".":
                self.hard_seqs = (hard_masked_genome.scaffolds[self.ch].seq[self.start-1:self.end], reverse_complement(hard_masked_genome.scaffolds[self.ch].seq[self.start-1:self.end]))
            
    def calculate_gc_content(self):
        if self.seq != "":
            gc_count = self.seq.count('G') + self.seq.count('C')
            self.gc_content = round((gc_count / self.size), 2)

    def print_gff(self):
        return(f"{self.ch}\t{self.source}\t{self.feature}\t{self.start}\t"
               f"{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t"
               f"{self.attributes}\n")

    def print_gtf(self):

        return(f"{self.ch}\t{self.source}\t{self.feature}\t{self.start}\t"
               f"{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t"
               f"{self.gtf_attributes}\n")
    
    def copy(self):
        return copy.deepcopy(self)
    
    def __str__(self):
        return str(self.id)

    def almost_equal(self, other):
        """
        Id or attributes may or may not be equal, it still returns True
        """
        return (self.start == other.start and self.end == other.end
                and self.ch == other.ch
                and self.strand == other.strand and self.score == other.score 
                and self.phase == other.phase)
    
    def equal_coordinates(self, other):
        return self.start == other.start and self.end == other.end and self.ch == other.ch
    
    def __lt__(self, other):
        return (self.start < other.start) or (self.start == other.start and self.end < other.end)
    
    def __le__(self, other):
        return (self.start < other.start) or (self.start == other.start and self.end <= other.end)
    
    def __eq__(self, other):
        """
        Is a feature is exactly the same as each other
        """
        return (self.start == other.start and self.end == other.end
                and self.ch == other.ch and self.id == other.id
                and self.source == other.source
                and self.strand == other.strand and self.score == other.score 
                and self.phase == other.phase and
                self.attributes == other.attributes)

    def longer(self, other:object):
        if self.seq != "" and other.seq != "":
            if len(self.seq) >= len(other.seq):
                return True
            else:
                return False
        else:
            print(f"Error: Either {self.id} or {other.id} sequences are empty!")

    def compare_blast_hits(self, other:object, source_priority:list):
        compared = False
        query_best = True
        while not compared:
            for s in source_priority:
                query_evalue = float(2)
                query_bitscore = float(-1)
                target_evalue = float(2)
                target_bitscore = float(-1)
                
                for b in self.blast_hits:
                    if b.source == s:
                        if b.evalue < query_evalue:
                            query_evalue = b.evalue
                        if b.score > query_bitscore:
                            query_bitscore = b.score

                for b in other.blast_hits:
                    if b.source == s:
                        if b.evalue < target_evalue:
                            target_evalue = b.evalue
                        if b.score > target_bitscore:
                            target_bitscore = b.score

                if query_evalue < target_evalue:
                    compared = True
                    break
                elif query_evalue > target_evalue:
                    query_best = False
                    compared = True
                    break
                elif query_bitscore > target_bitscore:
                    compared = True
                    break
                elif query_bitscore < target_bitscore:
                    query_best = False
                    compared = True
                    break
                elif s == source_priority[-1]:
                    if hasattr(self, "coding") and hasattr(other, "coding"):
                        if other.coding and not self.coding:
                            query_best = False
                    elif hasattr(other, "coding"):
                        query_best = False
                    compared = True

        return query_best 