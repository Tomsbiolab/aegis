from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .genome import Genome

import copy
import re

from .utils.genefunctions import reverse_complement
from .other_components import FeatureQuality

from functools import total_ordering

@total_ordering
class Feature():
    """
    Parent class including basic gff feature properties, that can then be
    inherited by Gene, Transcript etc. -> this facilitates code reading and
    editing. Or used on its own for CDS segments.
    """

    __slots__ = ('id', 'original_id', 'ch', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'frame', 'gtf_attributes', 'gene_id', 'size', 'parents', 'names', 'symbols', 'descriptors', 'processes', 'synonyms', 'aliases', 'renamed', 'id_number', 'original_id_number', 'misc_attributes', 'extra_copy', 'coding', '_quality')
    
    _ACTIVE_GENOME: Genome|None = None
    _ACTIVE_HARD_GENOME: Genome | None = None

    _ID_NUMBER_RE = re.compile(r'(\d+)$')

    gtf_attributes: list[str]|None
    size: int
    start: int
    end: int
    parents:list[str]|None
    misc_attributes:None|list[str]
    id_number: int|None
    original_id_number: int|None
    gene_id: str|None
    temp_attributes: list[str]
    descriptors:list|None
    synonyms:list|None
    aliases:list|None
    names:list|None
    symbols:list|None
    frame:None|int
    phase:None|int
    strand:str
    feature_id:str
    original_id:str
    coding:bool
    _quality: FeatureQuality | None

    # These attributes cannot be mistaken by misc attributes or any other
    attributes_to_ignore_when_reading_gff:set = {"id", "parent", "reliable_score", "remove", "rescue", "blasts", "gene_masked_fraction", "transcript_masked_fraction", "cds_masked_fraction", "gene_gc_content", "transcript_gc_content", "cds_gc_content", "intron_nested", "intron_nested_fully_contained", "intron_nested_single", "intron_utr_nested", "pseudogene", "transposable", "alternative_transcript_rescue", "cds_orientated_overlaps", "gene_id"}
    def __init__(self, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        
        self.id = feature_id
        self.original_id = feature_id
        self.ch = ch
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.phase = None
        self.frame = None
        self.coding = False
        self.gene_id = None
        if parents:
            self.parents = parents[:]
        else:
            self.parents = None
       
        self.gtf_attributes = None
        self.size = (self.end - self.start) + 1
        self.names = None
        self.symbols = None
        self.descriptors = None
        self.synonyms = None
        
        self.aliases = None
        
        self.renamed = False
        self.id_number = None
        self.original_id_number = None

        self.misc_attributes = None

        self.extra_copy = False

        self._quality = None

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
                if self.misc_attributes is None:
                    self.misc_attributes = []
                self.misc_attributes.append(f"{key}={value}")
            elif key.lower() in Feature.attributes_to_ignore_when_reading_gff:
                continue
            else:
                if self.misc_attributes is None:
                    self.misc_attributes = []
                self.misc_attributes.append(f"{key}={value}")

        self.update_numbering(original=True)

    def update_numbering(self, original:bool=False):

        match = Feature._ID_NUMBER_RE.search(self.id)
        if match:
            if original:
                self.original_id_number = int(match.group(1))
            self.id_number = int(match.group(1))

    def update_size(self):
        self.size = (self.end - self.start) + 1

    def print_gff(self, clean:bool=False, names:bool=False, symbols:bool=False, aliases:bool=False, symbols_as_description:bool=False, featurecountsID:bool=False, print_empty_attributes:bool=False):

        temp_attributes = [f"ID={self.id}"]

        if self.parents:
            parent_string = ",".join(self.parents)
            if parent_string or print_empty_attributes:
                temp_attributes.append(f"Parent={parent_string}")

        if featurecountsID and self.gene_id is not None:
            temp_attributes.append(f"Gene_id={self.gene_id}")

        if symbols and self.symbols is not None:
            symbol_string = ",".join(self.symbols)
            if symbol_string or print_empty_attributes:
                temp_attributes.append(f"Symbol={symbol_string}")

        if symbols_as_description and self.symbols is not None:
            symbol_string = ",".join(self.symbols)
            if symbol_string or print_empty_attributes:
                temp_attributes.append(f"Description={symbol_string}")

        if names and self.names is not None:
            name_string = ",".join(self.names)
            if name_string or print_empty_attributes:
                temp_attributes.append(f"Name={name_string}")

        if aliases and self.aliases is not None:
            alias_string = ",".join(self.aliases)
            if alias_string or print_empty_attributes:
                temp_attributes.append(f"Alias={alias_string}")

        if not clean and self.misc_attributes is not None:
            temp_attributes.extend(self.misc_attributes)

        attribute_string = ";".join(temp_attributes)
        phase = self.phase if self.phase is not None else "."

        return(f"{self.ch}\t{self.source}\t{self.feature}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{phase}\t{attribute_string}\n")

    def print_gtf(self):
        attribute_string = "; ".join(self.gtf_attributes) # type: ignore
        phase = self.phase if self.phase is not None else "."
        return(f"{self.ch}\t{self.source}\t{self.feature}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{phase}\t{attribute_string}\n")
    
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

        if self.size > other.size:
            return True
        else:
            return False

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

        elif overlap_bp < 0:
            overlap_bp = 0

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
                
                for b in self.quality.blast_hits:
                    if b.source == s:
                        if b.evalue < query_evalue:
                            query_evalue = b.evalue
                        if b.score > query_bitscore:
                            query_bitscore = b.score

                for b in other.quality.blast_hits:
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

    @property
    def quality(self) -> FeatureQuality:
        if self._quality is None:
            self._quality = FeatureQuality(self)
        return self._quality

    @property
    def seq(self) -> str|None:
        if not self._ACTIVE_GENOME:
            raise ValueError("No genome loaded and you are trying to access the sequence. Load your genome together with your annotation.")
        else:
            if self.strand == "+":
                return self._ACTIVE_GENOME.scaffolds[self.ch].seq[self.start-1:self.end]
            elif self.strand == "-":
                return reverse_complement(self._ACTIVE_GENOME.scaffolds[self.ch].seq[self.start-1:self.end])

    @property
    def hard_seq(self) -> str|None:
        if not self._ACTIVE_HARD_GENOME:
            raise ValueError("No hard masked genome loaded and you are trying to access the hard masked sequence. Load your hard masked genome together with your annotation.")
        else:
            if self.strand == "+":
                return self._ACTIVE_HARD_GENOME.scaffolds[self.ch].seq[self.start-1:self.end]
            elif self.strand == "-":
                return reverse_complement(self._ACTIVE_HARD_GENOME.scaffolds[self.ch].seq[self.start-1:self.end])

    @property
    def seqs(self) -> list[str]|None:
        if not self._ACTIVE_GENOME:
            raise ValueError("No genome loaded and you are trying to access the sequence. Load your genome together with your annotation.")
        else:
            raw = self._ACTIVE_GENOME.scaffolds[self.ch].seq[self.start-1:self.end]
            return [raw, reverse_complement(raw)]

    @property
    def hard_seqs(self) -> list[str]|None:
        if not self._ACTIVE_HARD_GENOME:
            raise ValueError("No hard masked genome loaded and you are trying to access the hard masked sequence. Load your hard masked genome together with your annotation.")
        else:
            raw = self._ACTIVE_HARD_GENOME.scaffolds[self.ch].seq[self.start-1:self.end]
            return [raw, reverse_complement(raw)]