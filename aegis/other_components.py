from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .gene import Gene
    from .feature import Feature
    from .transcript import Transcript
    from .subfeatures import CDS
    from .hits import BlastHit

import copy

from .utils.misc import count_occurrences

class FeatureAttributes():
    __slots__ = ('names', 'symbols', 'aliases', 'descriptors', 'synonyms', 'misc')

    names: list[str] | None
    symbols: list[str] | None
    aliases: list[str] | None
    descriptors: list[str] | None
    synonyms: list[str] | None
    misc: list[str] | None

    def __init__(self, names=None, symbols=None, aliases=None, descriptors=None, synonyms=None, misc=None):
        self.names = names
        self.symbols = symbols
        self.aliases = aliases
        self.descriptors = descriptors
        self.synonyms = synonyms
        self.misc = misc

    def __eq__(self, other):
        if not isinstance(other, FeatureAttributes):
            return False
        for slot in self.__slots__:
            if getattr(self, slot) != getattr(other, slot):
                return False
        return True

    def copy(self):
        return copy.deepcopy(self)

class GeneSynteny():
    __slots__ = ('previous', 'next', 'order', 'old_previous', 'old_next', 'old_order', 'liftover_conserved')

    previous: str | None
    next: str | None
    order: int | None
    old_previous: str | None
    old_next: str | None
    old_order: int | None
    liftover_conserved: bool
            
    def __init__(self):
        self.previous = None
        self.next = None
        self.order = None
        self.old_previous = None
        self.old_next = None
        self.old_order = None
        self.liftover_conserved = False

class FeatureQuality():
    __slots__ = ('_feature', 'remove', 'rescue', 'reliable', 'reliable_score', 'overlap_reliable', 'unrescuable', 'overlap_with_selected_CDS', 'overlap_with_selected_exon', 'intron_nested', 'intron_nested_fully_contained', 'intron_nested_single', 'UTR_intron_nested', 'transcriptomic_evidence', 'abinitio_evidence', 'gc_content', 'masked_fraction', '_blast_hits', '_alternative_transcript_rescue', 'potential_alt_transcript_group')

    _feature:Feature|Gene|Transcript|CDS
    remove:bool
    rescue:bool
    reliable:bool
    reliable_score:int
    overlap_reliable:bool
    unrescuable:bool
    overlap_with_selected_CDS:bool
    overlap_with_selected_exon:bool
    _alternative_transcript_rescue:list[str]|None
    intron_nested:bool
    intron_nested_fully_contained:bool
    intron_nested_single:bool
    UTR_intron_nested:bool
    transcriptomic_evidence:bool
    abinitio_evidence:bool
    gc_content:float
    masked_fraction:float
    potential_alt_transcript_group:int|None
    _blast_hits:list[BlastHit]|None

    def __init__(self, feature:Gene|Feature|Transcript|CDS):

        self._feature = feature
        self._blast_hits = None
        self._alternative_transcript_rescue = None

        self.gc_content = 0
        self.masked_fraction = 0

        self.remove = False
        self.rescue = False
        self.reliable = False
        self.reliable_score = 0
        self.overlap_reliable = False
        self.unrescuable = False

        self.overlap_with_selected_CDS = False
        self.overlap_with_selected_exon = False
        self.intron_nested = False
        self.intron_nested_fully_contained = False
        self.intron_nested_single = False
        self.UTR_intron_nested = False

        self.transcriptomic_evidence = False
        self.abinitio_evidence = False

        self.potential_alt_transcript_group = None

    def calculate_masking(self):
        self.masked_fraction = round(((count_occurrences(self._feature.hard_seq, "X") + (count_occurrences(self._feature.hard_seq, "N"))) / self._feature.size), 2)

    def calculate_gc_content(self):
        if self._feature.seq:
            gc_count = self._feature.seq.count('G') + self._feature.seq.count('C')
            self.gc_content = round((gc_count / self._feature.size), 1)

    def get_attributes(self) -> list[str]:
        temp_attributes = []

        temp_attributes.append(f"reliable_score={self.reliable_score}")
        temp_attributes.append(f"remove={self.remove}")
        temp_attributes.append(f"rescue={self.rescue}")
        blasts = []
        if self.blast_hits is not None:
            for b in self.blast_hits:
                blasts.append(f"{b.source}_{b.score}")
        blasts = ",".join(blasts)

        if blasts:
            temp_attributes.append(f"blasts={blasts}")
        
        alternative_transcript_rescue = ",".join(self.alternative_transcript_rescue)

        if alternative_transcript_rescue:
            temp_attributes.append(f"alternative_transcript_rescue={alternative_transcript_rescue}")

        if isinstance(self._feature, Gene):

            overlaps = []
            for o in self._feature.overlaps["self"]:
                if o.score >= 5:
                    overlaps.append(o.id)
            overlaps = ",".join(overlaps)
            if overlaps:
                temp_attributes.append(f"CDS_orientated_overlaps={overlaps}")

            gene_GC_content = self.gc_content
            gene_masked_fraction = self.masked_fraction

            transcript_masked_fraction = 0
            CDS_masked_fraction = 0

            transcript_GC_content = 0
            CDS_GC_content = 0

            for t in self._feature.transcripts.values():
                if t.main:
                    transcript_masked_fraction = t.quality.masked_fraction
                    transcript_GC_content = t.quality.gc_content
                    for c in t.CDSs.values():
                        if c.main:
                            CDS_masked_fraction = c.quality.masked_fraction
                            CDS_GC_content = c.quality.gc_content

            temp_attributes.append(f"gene_masked_fraction={gene_masked_fraction}")
            temp_attributes.append(f"transcript_masked_fraction={transcript_masked_fraction}")
            temp_attributes.append(f"CDS_masked_fraction={CDS_masked_fraction}")
            temp_attributes.append(f"gene_GC_content={gene_GC_content}")
            temp_attributes.append(f"transcript_GC_content={transcript_GC_content}")
            temp_attributes.append(f"CDS_GC_content={CDS_GC_content}")

        elif isinstance(self._feature, Transcript):
            transcript_GC_content = self.gc_content
            transcript_masked_fraction = self.masked_fraction
            CDS_masked_fraction = 0
            CDS_GC_content = 0

            transcript_masked_fraction = self._feature.quality.masked_fraction
            transcript_GC_content = self._feature.quality.gc_content
            for c in self._feature.CDSs.values():
                if c.main:
                    CDS_masked_fraction = c.quality.masked_fraction
                    CDS_GC_content = c.quality.gc_content

            temp_attributes.append(f"transcript_masked_fraction={transcript_masked_fraction}")
            temp_attributes.append(f"CDS_masked_fraction={CDS_masked_fraction}")
            temp_attributes.append(f"transcript_GC_content={transcript_GC_content}")
            temp_attributes.append(f"CDS_GC_content={CDS_GC_content}")

        elif isinstance(self._feature, CDS):
            CDS_GC_content = self.gc_content
            CDS_masked_fraction = self.masked_fraction

            if self._feature.main:
                CDS_masked_fraction = self._feature.quality.masked_fraction
                CDS_GC_content = self._feature.quality.gc_content

            temp_attributes.append(f"CDS_masked_fraction={CDS_masked_fraction}")
            temp_attributes.append(f"CDS_GC_content={CDS_GC_content}")

        temp_attributes.append(f"intron_nested={self.intron_nested}")
        temp_attributes.append(f"intron_nested_fully_contained={self.intron_nested_fully_contained}")
        temp_attributes.append(f"intron_nested_single={self.intron_nested_single}")
        temp_attributes.append(f"intron_UTR_nested={self.UTR_intron_nested}")

        return temp_attributes

    @property
    def blast_hits(self):
        if self._blast_hits is None:
            self._blast_hits = []
        return self._blast_hits

    @property
    def alternative_transcript_rescue(self):
        if self._alternative_transcript_rescue is None:
            self._alternative_transcript_rescue = []
        return self._alternative_transcript_rescue
