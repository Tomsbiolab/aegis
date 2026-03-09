from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .genome import Genome

import copy
import re

from .utils.genefunctions import reverse_complement
from .utils.misc import count_occurrences

from functools import total_ordering

@total_ordering
class Feature():
    """
    Parent class including basic gff feature properties, that can then be
    inherited by Gene, Transcript etc. -> this facilitates code reading and
    editing. Or used on its own for CDS segments.
    """

    __slots__ = (
        'id', 'original_id', 'ch', 'source', 'feature', 'start', 'end',
        'score', 'strand', 'phase', 'frame', 'gtf_attributes', 'gene_id',
        'size', 'seq', 'hard_seq', 'seqs', 'hard_seqs', 'parents', 'names',
        'symbols', 'descriptors', 'processes', 'synonyms', 'gc_content',
        'aliases', 'blast_hits', 'renamed', 'id_number', 'original_id_number',
        'misc_attributes', 'extra_copy', 'masked_fraction', 'coding'
    )

    gtf_attributes: list[str]
    size: int
    start: int
    end: int
    parents:list[str]
    misc_attributes:list[str]
    seqs: list[str]
    hard_seqs: list[str]
    id_number: int|None
    original_id_number: int|None
    gene_id: str|None
    temp_attributes: list[str]

    # These attributes cannot be mistaken by misc attributes or any other
    attributes_to_ignore_when_reading_gff:set = {"id", "parent", "reliable_score", "remove", "rescue", "blasts", "gene_masked_fraction", "transcript_masked_fraction", "cds_masked_fraction", "gene_gc_content", "transcript_gc_content", "cds_gc_content", "intron_nested", "intron_nested_fully_contained", "intron_nested_single", "intron_utr_nested", "pseudogene", "transposable", "alternative_transcript_rescue", "cds_orientated_overlaps", "gene_id"}
    def __init__(self, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, phase:str, parents:list[str]=[], attributes:dict={}):
        
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
        self.coding = False
        self.gene_id = None
        self.parents = parents[:]
       
        self.gtf_attributes = []
        self.size = (self.end - self.start) + 1
        self.seq = ""
        self.hard_seq = ""
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

        self.misc_attributes = []

        self.extra_copy = False

        for key, value in attributes.items():

            if key == "Alias":
                self.aliases = value[:]

            elif key == "Name":
                self.names = value[:]
        
            elif key == "Symbol":
                self.symbols = value[:]

            elif "extra_copy_number" in key:
                value = int(value.strip())
                if value > 0:
                    self.extra_copy = True
                self.misc_attributes.append(f"{key}={value}")
            elif key.lower in Feature.attributes_to_ignore_when_reading_gff:
                continue
            else:
                self.misc_attributes.append(f"{key}={value}")

        self.update_numbering(original=True)

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

    def generate_sequence(self, genome:Genome):
        if self.start != "NA" and self.end != "NA":
            if self.strand == "+":
                self.seq = genome.scaffolds[self.ch].seq[self.start-1:self.end]
            elif self.strand == "-":
                self.seq = reverse_complement(genome.scaffolds[self.ch].seq[self.start-1:self.end])
            elif self.strand == ".":
                self.seqs = [genome.scaffolds[self.ch].seq[self.start-1:self.end], reverse_complement(genome.scaffolds[self.ch].seq[self.start-1:self.end])]

    def clear_sequence(self, just_hard=False):
        self.hard_seq = ""
        if hasattr(self, "hard_seqs"):
            del self.hard_seqs
        if not just_hard:
            self.seq = ""
            if hasattr(self, "seqs"):
                del self.seqs

    def generate_hard_sequence(self, hard_masked_genome:Genome):
        if self.start != "NA" and self.end != "NA":
            if self.strand == "+":
                self.hard_seq = hard_masked_genome.scaffolds[self.ch].seq[self.start-1:self.end]
            elif self.strand == "-":
                self.hard_seq = reverse_complement(hard_masked_genome.scaffolds[self.ch].seq[self.start-1:self.end])
            elif self.strand == ".":
                self.hard_seqs = [hard_masked_genome.scaffolds[self.ch].seq[self.start-1:self.end], reverse_complement(hard_masked_genome.scaffolds[self.ch].seq[self.start-1:self.end])]
            
    def calculate_gc_content(self):
        if self.seq != "":
            gc_count = self.seq.count('G') + self.seq.count('C')
            self.gc_content = round((gc_count / self.size), 2)

    def print_gff(self, clean:bool=False, names:bool=False, symbols:bool=False, aliases:bool=False, symbols_as_description:bool=False, featurecountsID:bool=False):

        parent_string = ",".join(self.parents)
        temp_attributes = [f"ID={self.id}", f"Parent={parent_string}"]
        if featurecountsID:
            temp_attributes.append(f"Gene_id={self.gene_id}")

        if symbols:
            symbol_string = ",".join(self.symbols)
            temp_attributes.append(f"Symbol={symbol_string}")

        if symbols_as_description:
            symbol_string = ",".join(self.symbols)
            temp_attributes.append(f"Description={symbol_string}")

        if names:
            name_string = ",".join(self.names)
            temp_attributes.append(f"Name={name_string}")

        if aliases:
            alias_string = ",".join(self.aliases)
            temp_attributes.append(f"Alias={alias_string}")

        if not clean:
            temp_attributes.extend(self.misc_attributes)

        attribute_string = ";".join(temp_attributes)

        return(f"{self.ch}\t{self.source}\t{self.feature}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t{attribute_string}\n")

    def print_gtf(self):
        attribute_string = "; ".join(self.gtf_attributes)
        return(f"{self.ch}\t{self.source}\t{self.feature}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t{attribute_string}\n")
    
    def copy(self):
        return copy.deepcopy(self)
    
    def __str__(self):
        return str(self.id)

    def equal_sequence(self, other):
        return (self.start == other.start and self.end == other.end and self.ch == other.ch and self.strand == other.strand)
    
    def equal_coordinates(self, other):
        return self.start == other.start and self.end == other.end and self.ch == other.ch
    
    def __lt__(self, other):
        return (self.start < other.start) or (self.start == other.start and self.end < other.end)
    
    def __le__(self, other):
        return (self.start < other.start) or (self.start == other.start and self.end <= other.end)

    def longer(self, other:Feature):
        if self.seq != "" and other.seq != "":
            if len(self.seq) >= len(other.seq):
                return True
            else:
                return False
        else:
            print(f"Error: Either {self.id} or {other.id} sequences are empty!")

    def identical(self, other:Feature) -> bool:
        """
        Checks if all attributes in __slots__ are identical.
        Two attributes are identical if:
        1. Both are uninitialized.
        2. Both are initialized and have the same value.
        """
        if not isinstance(other, type(self)):
            return False

        for attr in self.__slots__:
            self_has = hasattr(self, attr)
            other_has = hasattr(other, attr)

            if self_has != other_has:
                return False

            if self_has:
                if getattr(self, attr) != getattr(other, attr):
                    return False
        return True

    def overlap(self, other:Feature) -> tuple[bool, int]:
        """
        Make sure the features have the same .ch (or chromosome/scaffold) before calling this function
        """
        overlapping = False

        overlap_start = max(self.start, other.start)
        overlap_end = min(self.end, other.end)
        
        # Calculate overlap length
        overlap_bp = (overlap_end - overlap_start) + 1

        # checking only overlapping features
        if overlap_bp > 0:
            overlapping = True

        return overlapping, overlap_bp

    def compare_blast_hits(self, other:Feature, source_priority:list):
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