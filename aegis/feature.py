from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .genome import Genome

import copy
import re
import warnings

from .utils.genefunctions import reverse_complement
from .utils.gtf_gff import format_gff3_values, DEFAULT_WRITE_KEYS, SOURCE_KEYS, EXTRA_ATTRIBUTES, INT_NAME, INT_SYMBOL, INT_ALIAS, INT_PARENT
from .other_components import FeatureQuality, FeatureAttributes
from .conf import default_features

from functools import total_ordering

# Name, Alias and symbol, in the order they are written unless the input had them elsewhere
IDENTITY_ATTRIBUTES = (INT_NAME, INT_ALIAS, INT_SYMBOL)
_IDENTITY_RANK = {a: r for r, a in enumerate(IDENTITY_ATTRIBUTES)}

# Shared copies of FeatureAttributes.order tuples: a file tends to repeat a few layouts on every row
_ORDER_CACHE: dict[tuple, tuple] = {}

def _is_default_identity_order(order) -> bool:
    """Whether identity attributes were read before every misc attribute and in IDENTITY_ATTRIBUTES order."""
    return all(n == 0 for _, n in order) and [a for a, _ in order] == sorted((a for a, _ in order), key=_IDENTITY_RANK.__getitem__)

def merge_identity_and_misc(identity: list[tuple[str, str]], misc: list[str] | None, order: tuple | None) -> list[str]:
    """
    Interleaves identity attribute strings ((attribute, 'key=value') pairs, in default order) with
    misc attribute strings as read (see FeatureAttributes.order). Identity attributes the input did
    not have go first; positions past the end of misc (e.g. misc changed since) go last.
    """
    misc = misc or []
    if not order:
        return [s for _, s in identity] + misc
    position = {a: (n, r) for r, (a, n) in enumerate(order)}
    out = [s for a, s in identity if a not in position]
    placed = sorted(((position[a], s) for a, s in identity if a in position), key=lambda x: x[0])
    j = 0
    for i, m in enumerate(misc):
        while j < len(placed) and placed[j][0][0] <= i:
            out.append(placed[j][1])
            j += 1
        out.append(m)
    out.extend(s for _, s in placed[j:])
    return out

@total_ordering
class Feature():
    """
    Parent class including basic gff feature properties, that can then be
    inherited by Gene, Transcript etc. -> this facilitates code reading and
    editing. Or used on its own for CDS segments.
    """

    __slots__ = ('id', 'original_id', 'ch', 'source', 'feature', 'start', 'end', 'score', 'strand', 'phase', 'frame', 'gtf_attributes', 'gene_id', 'parents', 'renamed', 'id_number', 'original_id_number', 'extra_copy', 'coding', '_quality', '_attributes')
    
    _ACTIVE_GENOME: Genome | None = None
    _ACTIVE_HARD_GENOME: Genome | None = None

    # Trailing number of an id (used to keep numbering when renaming); subfeature classes only
    # accept it after their type name, so e.g. the version in 'cds-NP_065184.2' is not read as one
    _ID_NUMBER_RE = re.compile(r'(\d+)$')

    id: str
    gtf_attributes: list[str]|None
    start: int
    end: int
    parents:list[str]|None
    id_number: int|None
    original_id_number: int|None
    gene_id: str|None
    temp_attributes: list[str]
    frame:None|int
    phase:None|int
    strand:str
    feature_id:str
    original_id:str
    coding:bool
    _quality: FeatureQuality | None
    _attributes: FeatureAttributes | None
    ch:str

    # Structural keys and AEGIS's own processing state (written by AEGIS and recomputed rather than read back)
    aegis_state_attributes:set = {"id", "parent", "reliable_score", "remove", "rescue", "blasts", "gene_masked_fraction", "transcript_masked_fraction", "cds_masked_fraction", "gene_gc_content", "transcript_gc_content", "cds_gc_content", "intron_nested", "intron_nested_fully_contained", "intron_nested_single", "intron_utr_nested", "alternative_transcript_rescue", "cds_orientated_overlaps", "gene_id"}
    # Model status written by other tools (NCBI, Liftoff) that AEGIS derives from its own translation; kept, it would go stale once a model changes
    derived_state_attributes:set = {"partial", "start_range", "end_range", "missing_start_codon", "missing_stop_codon", "inframe_stop_codon", "valid_orf", "valid_orfs"}
    attributes_to_ignore_when_reading_gff:set = aegis_state_attributes | derived_state_attributes
    # CDS attributes that change or explain translation; kept on CDS rows even in clean GFF output
    translation_attributes:set = {"transl_table", "transl_except", "exception"}
    def __init__(self, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}, phase:int|None=None):
        
        self.id = feature_id
        self.original_id = feature_id
        self.ch = ch
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        if strand == "+" or strand == "-":
            self.strand = strand
        else:
            self.strand = "."
        self.phase = phase
        self.frame = None
        self.coding = False
        self.gene_id = None
        if parents:
            self.parents = parents[:]
        else:
            self.parents = None
       
        self.gtf_attributes = None
        
        self.renamed = False
        self.id_number = None
        self.original_id_number = None

        self._attributes = None

        self.extra_copy = False

        self._quality = None

        # (identity attribute, number of misc attributes read before it), in column-9 order
        identity_order = []

        for key, value in attributes.items():

            if key in IDENTITY_ATTRIBUTES:
                if self._attributes is None:
                    self._attributes = FeatureAttributes()
                if key == INT_NAME:
                    self._attributes.names = value[:]
                elif key == INT_ALIAS:
                    self._attributes.aliases = value[:]
                else:
                    self._attributes.symbols = value[:]
                identity_order.append((key, len(self._attributes.misc or ())))

            elif key == SOURCE_KEYS:
                if self._attributes is None:
                    self._attributes = FeatureAttributes()
                self._attributes.keys = dict(value)

            elif key.startswith(EXTRA_ATTRIBUTES):
                if value.partition("=")[0].lower() not in Feature.attributes_to_ignore_when_reading_gff:
                    if self._attributes is None:
                        self._attributes = FeatureAttributes()
                    if self._attributes.misc is None:
                        self._attributes.misc = []
                    self._attributes.misc.append(value)

            elif "extra_copy_number" in key:
                value = int(value.strip())
                if value > 0:
                    self.extra_copy = True
                if self._attributes is None:
                    self._attributes = FeatureAttributes()
                    self._attributes.misc = []
                elif self._attributes.misc is None:
                    self._attributes.misc = []
                self._attributes.misc.append(f"{key}={value}")
            elif key.lower() in Feature.attributes_to_ignore_when_reading_gff:
                continue
            else:
                if self._attributes is None:
                    self._attributes = FeatureAttributes()
                    self._attributes.misc = []
                elif self._attributes.misc is None:
                    self._attributes.misc = []
                self._attributes.misc.append(f"{key}={value}")

        # Only kept when it is not the order identity attributes are written in anyway
        if identity_order and not _is_default_identity_order(identity_order):
            self.attribute_order = tuple(identity_order)

        self.update_numbering(original=True)

    def update_numbering(self, original:bool=False):

        match = self._ID_NUMBER_RE.search(self.id)
        if match:
            if original:
                self.original_id_number = int(match.group(1))
            self.id_number = int(match.group(1))

    @property
    def size(self) -> int:
        return (self.end - self.start) + 1

    @property
    def genomic_span(self) -> int:
        """Genomic chromosome span (end - start + 1), including any introns."""
        return (self.end - self.start) + 1

    @property
    def names(self) -> list[str] | None:
        return self._attributes.names if self._attributes else None

    @names.setter
    def names(self, value: list[str] | None):
        if self._attributes is None:
            self._attributes = FeatureAttributes()
        self._attributes.names = value

    @property
    def symbols(self) -> list[str] | None:
        return self._attributes.symbols if self._attributes else None

    @symbols.setter
    def symbols(self, value: list[str] | None):
        if self._attributes is None:
            self._attributes = FeatureAttributes()
        self._attributes.symbols = value

    @property
    def aliases(self) -> list[str] | None:
        return self._attributes.aliases if self._attributes else None

    @aliases.setter
    def aliases(self, value: list[str] | None):
        if self._attributes is None:
            self._attributes = FeatureAttributes()
        self._attributes.aliases = value

    @property
    def descriptors(self) -> list[str] | None:
        return self._attributes.descriptors if self._attributes else None

    @descriptors.setter
    def descriptors(self, value: list[str] | None):
        if self._attributes is None:
            self._attributes = FeatureAttributes()
        self._attributes.descriptors = value

    @property
    def attribute_keys(self) -> dict[str, str] | None:
        return self._attributes.keys if self._attributes else None

    @attribute_keys.setter
    def attribute_keys(self, value: dict[str, str] | None):
        if self._attributes is None:
            self._attributes = FeatureAttributes()
        self._attributes.keys = value

    def gff_key(self, attribute: str) -> str:
        """GFF key written for attribute ("Name", "Alias" or "Symbol"): as read, else the default."""
        keys = self._attributes.keys if self._attributes else None
        if keys and attribute in keys:
            return keys[attribute]
        return DEFAULT_WRITE_KEYS[attribute]

    def set_gff_key(self, attribute: str, key: str | None):
        """Sets the GFF key spelling for attribute; None or the default spelling clears it."""
        keys = self.attribute_keys
        if key is None or key == DEFAULT_WRITE_KEYS[attribute]:
            if keys:
                keys.pop(attribute, None)
                if not keys:
                    self._attributes.keys = None # type: ignore
            return
        if keys is None:
            keys = {}
            self.attribute_keys = keys
        keys[attribute] = key

    @property
    def misc_attributes(self) -> list[str] | None:
        return self._attributes.misc if self._attributes else None

    @misc_attributes.setter
    def misc_attributes(self, value: list[str] | None):
        if self._attributes is None:
            self._attributes = FeatureAttributes()
        self._attributes.misc = value

    @property
    def attribute_order(self) -> tuple | None:
        return self._attributes.order if self._attributes else None

    @attribute_order.setter
    def attribute_order(self, value: tuple | None):
        if self._attributes is None:
            self._attributes = FeatureAttributes()
        self._attributes.order = _ORDER_CACHE.setdefault(value, value) if value else None

    def place_attribute_like(self, attribute: str, like: str):
        """
        Writes attribute (one the input did not have) right after like, where like was read in
        column 9, e.g. an Alias made from a Name. Without a recorded position for like (default
        layout) nothing changes: attribute is written in default order.
        """
        order = self.attribute_order
        if not order or any(a == attribute for a, _ in order):
            return
        for i, (a, n) in enumerate(order):
            if a == like:
                self.attribute_order = order[:i + 1] + ((attribute, n),) + order[i + 1:]
                return

    def identity_gff_attributes(self, names:bool=True, symbols:bool=True, aliases:bool=True, symbols_as_description:bool=False, print_empty_attributes:bool=False, misc:list[str]|None=None, keep_order:bool=True, missing_as_empty:bool=False) -> list[str]:
        """
        Name, Alias and symbol attributes, percent-encoded and keyed as read (see gff_key), merged
        with the misc attribute strings given in the column-9 order they were read in. The flags
        only allow an attribute: it is written when the feature has values, or always with
        print_empty_attributes. symbols_as_description adds 'description=' with symbol_values(), where
        the attribute they come from was read.
        keep_order=False writes them first, in IDENTITY_ATTRIBUTES order. missing_as_empty writes
        attributes the feature does not have as empty too (with print_empty_attributes).
        """
        pairs = []
        empty = [] if missing_as_empty else None
        for allowed, attribute, values in ((names, INT_NAME, self.names), (aliases, INT_ALIAS, self.aliases), (symbols, INT_SYMBOL, self.symbols)):
            if values is None:
                values = empty
            if allowed and values is not None:
                value_string = format_gff3_values(values)
                if value_string or print_empty_attributes:
                    pairs.append((attribute, f"{self.gff_key(attribute)}={value_string}"))
        if symbols_as_description:
            # What the feature is called (see symbol_values), placed where that attribute was read
            label = self.non_id_names()
            description, source = (label, INT_NAME) if label else (self.symbols or [], INT_SYMBOL)
            value_string = format_gff3_values(description)
            if value_string or (print_empty_attributes and missing_as_empty):
                pairs.append((source, f"description={value_string}"))
        return merge_identity_and_misc(pairs, misc, self.attribute_order if keep_order else None)

    def non_id_names(self) -> list[str]:
        """Names that are not, ignoring case, an ID the feature has or had (e.g. Araport 'ID=AT4G00005;Name=AT4G00005')."""
        ids = {self.id.casefold(), self.original_id.casefold()}
        return [n for n in (self.names or ()) if n.casefold() not in ids]

    def symbol_values(self) -> list[str]:
        """
        What the feature is called: its non_id_names() (Name is where 'aegis symbols' writes and what
        genome browsers label), else its Symbol values. Empty if it only has an ID.
        """
        return self.non_id_names() or self.symbols or []

    def has_misc_attribute(self, keys:set[str]) -> bool:
        """Whether a kept attribute has one of these (lowercase) keys."""
        return any(a.partition("=")[0].lower() in keys for a in (self.misc_attributes or ()))

    def print_gff(self, clean:bool=False, names:bool=False, symbols:bool=False, aliases:bool=False, symbols_as_description:bool=False, featurecountsID:bool=False, print_empty_attributes:bool=False):

        temp_attributes = [f"ID={self.id}"]

        if self.parents:
            parent_string = ",".join(self.parents)
            if parent_string or print_empty_attributes:
                # 'Derives_from' kept for the atypical features that use it (e.g. Araport11 protein rows)
                temp_attributes.append(f"{self.gff_key(INT_PARENT)}={parent_string}")

        if featurecountsID and self.gene_id is not None:
            temp_attributes.append(f"Gene_id={self.gene_id}")

        misc = self.misc_attributes
        if misc is not None and clean:
            misc = [a for a in misc if a.partition("=")[0] in Feature.translation_attributes] if self.feature in default_features["CDS"] else None

        temp_attributes.extend(self.identity_gff_attributes(names=names, symbols=symbols, aliases=aliases, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes, misc=misc, keep_order=not clean))

        attribute_string = ";".join(temp_attributes)
        phase = self.phase if (self.feature in default_features["CDS"] and self.phase is not None) else "."

        return(f"{self.ch}\t{self.source}\t{self.feature}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{phase}\t{attribute_string}\n")

    def print_gtf(self):
        attribute_string = "; ".join(self.gtf_attributes) # type: ignore
        phase = self.phase if (self.feature in default_features["CDS"] and self.phase is not None) else "."
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

    def _genome_slice(self, genome, genome_label: str = "genome") -> str | None:
        """
        Plus-strand sequence of the feature in genome, or None if its contig is missing.
        On circular sequences (organelles, or Is_circular=true in the GFF region line)
        coordinates past the end are read across the origin, as NCBI writes them;
        otherwise out-of-bounds coordinates are truncated.
        """
        scaffold = genome.scaffolds.get(self.ch)
        if scaffold is None:
            warnings.warn(
                f"Feature '{self.id}' is on contig '{self.ch}', which is missing from {genome_label} '{genome.name}'. Returning empty sequence.",
                category=UserWarning,
                stacklevel=3,
            )
            return None
        scf_seq = scaffold.seq
        scf_len = len(scf_seq)
        start = self.start
        end = self.end
        circular = getattr(scaffold, "organelle", False) or getattr(scaffold, "circular", False)
        if circular and scf_len and start is not None and end is not None and start >= 1 and end > scf_len:
            offset = ((start - 1) // scf_len) * scf_len
            start -= offset
            end -= offset
            if end <= scf_len:
                return scf_seq[start - 1:end]
            if end - scf_len < start:
                return scf_seq[start - 1:] + scf_seq[:end - scf_len]
        start_clipped = (start is not None and start < 1)
        end_clipped = (end is not None and end > scf_len)
        if start_clipped or end_clipped:
            warnings.warn(
                f"Feature '{self.id}' coordinates ({self.start}..{self.end}) exceed bounds for contig '{self.ch}' (length {scf_len}) in {genome_label} '{genome.name}'. Sequence will be truncated.",
                category=UserWarning,
                stacklevel=3,
            )
        start = max(0, start - 1) if start is not None else 0
        end = min(scf_len, end) if end is not None else scf_len
        if start >= end:
            return ""
        return scf_seq[start:end]

    @property
    def seq(self) -> str:
        if not self._ACTIVE_GENOME:
            raise ValueError("No genome loaded and you are trying to access the sequence. Load your genome together with your annotation.")
        raw = self._genome_slice(self._ACTIVE_GENOME)
        if not raw:
            return ""
        return reverse_complement(raw) if self.strand == "-" else raw

    @property
    def hard_seq(self) -> str:
        if not self._ACTIVE_HARD_GENOME:
            raise ValueError("No hard masked genome loaded and you are trying to access the hard masked sequence. Load your hard masked genome together with your annotation.")
        raw = self._genome_slice(self._ACTIVE_HARD_GENOME, "hard-masked genome")
        if not raw:
            return ""
        return reverse_complement(raw) if self.strand == "-" else raw

    @property
    def seqs(self) -> list[str]:
        if not self._ACTIVE_GENOME:
            raise ValueError("No genome loaded and you are trying to access the sequence. Load your genome together with your annotation.")
        raw = self._genome_slice(self._ACTIVE_GENOME)
        if not raw:
            return ["", ""]
        return [raw, reverse_complement(raw)]

    @property
    def hard_seqs(self) -> list[str]:
        if not self._ACTIVE_HARD_GENOME:
            raise ValueError("No hard masked genome loaded and you are trying to access the hard masked sequence. Load your hard masked genome together with your annotation.")
        raw = self._genome_slice(self._ACTIVE_HARD_GENOME, "hard-masked genome")
        if not raw:
            return ["", ""]
        return [raw, reverse_complement(raw)]
