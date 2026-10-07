"""
Created on Tue Dec 27 15:20:59 2022

Module with an array of genomic functions.

@authors: David Navarro, Antonio Santiago
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Literal, Optional, Union

if TYPE_CHECKING:
    from ..annotation import Annotation
    from ..gene import Gene
    from ..subfeatures import Feature

import pandas as pd
import time
import warnings
import itertools
import hashlib

from collections import defaultdict
from pathlib import Path

# nucleotides
_STR_FROM = "ACGTRYSWKMBDHVNXUacgtryswkmbdhvnxu-"
_STR_TO   = "TGCAYRSWMKVHDBNXAtgcayrswmkvhdbnxa-"
_BYTES_COMP_TABLE = bytes.maketrans(_STR_FROM.encode(), _STR_TO.encode())

def reverse_complement(in_seq: str) -> str:
    """Return the reverse complement of a nucleotide sequence (handles IUPAC, case, and RNA)."""
    return in_seq.encode('ascii').translate(_BYTES_COMP_TABLE)[::-1].decode('ascii')


def sequence_hash(in_seq: str) -> str:
    """Return a fast, deterministic MD5 checksum of an upper-cased nucleotide sequence."""
    return hashlib.md5(in_seq.upper().encode('ascii', errors='ignore')).hexdigest()

iupac_dna_nucleotides = {
    "W": ["A", "T"],
    "S": ["C", "G"],
    "M": ["A", "C"],
    "K": ["G", "T"],
    "R": ["A", "G"],
    "Y": ["C", "T"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"],
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"]
}

codon_dict = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C', 'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G', "TAA": "*", "TAG": "*", "TGA": "*"}

# Codon assignments and "stops" match NCBI's gc.prt (v4.6). NCBI's initiator codons are split
# in two: "starts" is a curated, strict subset (e.g. table 1 is ATG only) used for ORF
# prediction, since rare initiators would extend predicted ORFs upstream; "alt_starts" holds
# NCBI's remaining initiators (e.g. CTG/TTG in table 1), which only grade the start codon of
# proteins as alternative instead of missing. Pass start_codons=/stop_codons= explicitly to
# find_ORFs and the protein generation methods to use other sets. Tables 27-31 are not
# included: their stop codons are context-dependent and also encode amino acids.
NCBI_GENETIC_CODES: dict[int, dict] = {
    1: {
        "name": "Standard",
        "diff": {},
        "stops": ("TAA", "TAG", "TGA"),
        "starts": ("ATG",),
        "alt_starts": ("CTG", "TTG"),
    },
    2: {
        "name": "Vertebrate Mitochondrial",
        "diff": {"AGA": "*", "AGG": "*", "ATA": "M", "TGA": "W"},
        "stops": ("TAA", "TAG", "AGA", "AGG"),
        "starts": ("ATG", "ATA", "ATT", "ATC", "GTG"),
        "alt_starts": (),
    },
    3: {
        "name": "Yeast Mitochondrial",
        "diff": {"ATA": "M", "CTT": "T", "CTC": "T", "CTA": "T", "CTG": "T", "TGA": "W"},
        "stops": ("TAA", "TAG"),
        "starts": ("ATG", "ATA", "GTG"),
        "alt_starts": (),
    },
    4: {
        "name": "Mold, Protozoan, and Coelenterate Mitochondrial and Mycoplasma/Spiroplasma",
        "diff": {"TGA": "W"},
        "stops": ("TAA", "TAG"),
        "starts": ("ATG", "ATA", "ATT", "ATC", "GTG", "TTG"),
        "alt_starts": ("CTG", "TTA"),
    },
    5: {
        "name": "Invertebrate Mitochondrial",
        "diff": {"AGA": "S", "AGG": "S", "ATA": "M", "TGA": "W"},
        "stops": ("TAA", "TAG"),
        "starts": ("ATG", "ATA", "ATT", "ATC", "GTG", "TTG"),
        "alt_starts": (),
    },
    6: {
        "name": "Ciliate, Dasycladacean and Hexamita Nuclear",
        "diff": {"TAA": "Q", "TAG": "Q"},
        "stops": ("TGA",),
        "starts": ("ATG",),
        "alt_starts": (),
    },
    9: {
        "name": "Echinoderm and Flatworm Mitochondrial",
        "diff": {"AAA": "N", "AGA": "S", "AGG": "S", "TGA": "W"},
        "stops": ("TAA", "TAG"),
        "starts": ("ATG", "GTG"),
        "alt_starts": (),
    },
    10: {
        "name": "Euplotid Nuclear",
        "diff": {"TGA": "C"},
        "stops": ("TAA", "TAG"),
        "starts": ("ATG",),
        "alt_starts": (),
    },
    11: {
        "name": "Bacterial, Archaeal and Plant Plastid",
        "diff": {},
        "stops": ("TAA", "TAG", "TGA"),
        "starts": ("ATG", "GTG", "TTG"),
        "alt_starts": ("ATA", "ATC", "ATT", "CTG"),
    },
    12: {
        "name": "Alternative Yeast Nuclear",
        "diff": {"CTG": "S"},
        "stops": ("TAA", "TAG", "TGA"),
        "starts": ("ATG", "CTG"),
        "alt_starts": (),
    },
    13: {
        "name": "Ascidian Mitochondrial",
        "diff": {"AGA": "G", "AGG": "G", "ATA": "M", "TGA": "W"},
        "stops": ("TAA", "TAG"),
        "starts": ("ATG", "ATA", "GTG", "TTG"),
        "alt_starts": (),
    },
    14: {
        "name": "Alternative Flatworm Mitochondrial",
        "diff": {"AAA": "N", "AGA": "S", "AGG": "S", "TGA": "W", "TAA": "Y"},
        "stops": ("TAG",),
        "starts": ("ATG",),
        "alt_starts": (),
    },
    15: {
        "name": "Blepharisma Macronuclear",
        "diff": {"TAG": "Q"},
        "stops": ("TAA", "TGA"),
        "starts": ("ATG",),
        "alt_starts": (),
    },
    16: {
        "name": "Chlorophycean Mitochondrial",
        "diff": {"TAG": "L"},
        "stops": ("TAA", "TGA"),
        "starts": ("ATG",),
        "alt_starts": (),
    },
    21: {
        "name": "Trematode Mitochondrial",
        "diff": {"AAA": "N", "AGA": "S", "AGG": "S", "ATA": "M", "TGA": "W"},
        "stops": ("TAA", "TAG"),
        "starts": ("ATG", "GTG"),
        "alt_starts": (),
    },
    22: {
        "name": "Scenedesmus obliquus Mitochondrial",
        "diff": {"TCA": "*", "TAG": "L"},
        "stops": ("TAA", "TCA", "TGA"),
        "starts": ("ATG",),
        "alt_starts": (),
    },
    23: {
        "name": "Thraustochytrium Mitochondrial",
        "diff": {"TTA": "*"},
        "stops": ("TAA", "TAG", "TGA", "TTA"),
        "starts": ("ATG", "GTG", "ATT"),
        "alt_starts": (),
    },
    24: {
        "name": "Rhabdopleuridae Mitochondrial",
        "diff": {"AGA": "S", "AGG": "K", "TGA": "W"},
        "stops": ("TAA", "TAG"),
        "starts": ("ATG", "GTG"),
        "alt_starts": ("CTG", "TTG"),
    },
    25: {
        "name": "Candidate Division SR1 and Gracilibacteria",
        "diff": {"TGA": "G"},
        "stops": ("TAA", "TAG"),
        "starts": ("ATG", "GTG", "TTG"),
        "alt_starts": (),
    },
    26: {
        "name": "Pachysolen tannophilus Nuclear",
        "diff": {"CTG": "A"},
        "stops": ("TAA", "TAG", "TGA"),
        "starts": ("ATG",),
        "alt_starts": ("CTG",),
    },
    32: {
        "name": "Balanophoraceae Plastid",
        "diff": {"TAG": "W"},
        "stops": ("TAA", "TGA"),
        "starts": ("ATG", "GTG", "TTG"),
        "alt_starts": ("ATA", "ATC", "ATT", "CTG"),
    },
    33: {
        "name": "Cephalodiscidae Mitochondrial",
        "diff": {"TAA": "Y", "TGA": "W", "AGA": "S", "AGG": "K"},
        "stops": ("TAG",),
        "starts": ("ATG", "GTG"),
        "alt_starts": ("CTG", "TTG"),
    },
}

TAXONOMY_ORGANELLE_CODES: dict[str, dict[str, int]] = {
    "plant": {
        "table": 1,
        "mito_table": 1,      # Plant mitochondria use Standard Code (NCBI Table 1)
        "plastid_table": 11,  # Plant plastid / Chloroplast uses NCBI Table 11
    },
    "vertebrate": {
        "table": 1,
        "mito_table": 2,      # Vertebrate Mitochondrial (NCBI Table 2)
        "plastid_table": 11,
    },
    "invertebrate": {
        "table": 1,
        "mito_table": 5,      # Invertebrate Mitochondrial (NCBI Table 5)
        "plastid_table": 11,
    },
    "yeast": {
        "table": 1,
        "mito_table": 3,      # Yeast Mitochondrial (NCBI Table 3)
        "plastid_table": 11,
    },
}

def resolve_taxonomy_tables(
    taxonomy: str = "plant",
    table: int | str | None = None,
    mito_table: int | str | None = None,
    plastid_table: int | str | None = None,
) -> tuple[int | str, int | str, int | str]:
    """
    Resolves nuclear, mitochondrial, and plastid translation tables based on a taxonomy preset,
    allowing explicit overrides for individual table types.
    """
    tax_key = str(taxonomy).lower() if taxonomy else "plant"
    if tax_key not in TAXONOMY_ORGANELLE_CODES:
        raise ValueError(
            f"Unknown taxonomy preset: '{taxonomy}'. Available presets: {list(TAXONOMY_ORGANELLE_CODES.keys())}"
        )
    preset = TAXONOMY_ORGANELLE_CODES[tax_key]
    res_table = table if table is not None else preset["table"]
    res_mito = mito_table if mito_table is not None else preset["mito_table"]
    res_plastid = plastid_table if plastid_table is not None else preset["plastid_table"]
    return res_table, res_mito, res_plastid

def _build_extended_codon_table(base_dict: dict[str, str]) -> tuple[dict[str, str], defaultdict]:
    """Generates an IUPAC-extended codon table considering all IUPAC ambiguous nucleotide combinations."""
    extended = {}
    all_iupac = list(iupac_dna_nucleotides.keys())
    for c1, c2, c3 in itertools.product(all_iupac, repeat=3):
        ambiguous = c1 + c2 + c3
        possible_aas = set()
        for b1, b2, b3 in itertools.product(iupac_dna_nucleotides[c1], iupac_dna_nucleotides[c2], iupac_dna_nucleotides[c3]):
            possible_aas.add(base_dict[b1 + b2 + b3])
        if len(possible_aas) == 1:
            extended[ambiguous] = possible_aas.pop()
        else:
            extended[ambiguous] = "X"
    byte_codon = {tuple(k.encode('ascii')): v for k, v in extended.items()}
    b_dict = defaultdict(lambda: "X", byte_codon)
    return extended, b_dict

extended_codon_dict, byte_dict = _build_extended_codon_table(codon_dict)
byte_codon_dict = {tuple(k.encode('ascii')): v for k, v in extended_codon_dict.items()}

_GENETIC_CODE_CACHE: dict[int, tuple[dict[str, str], dict[str, str], defaultdict, tuple[str, ...], tuple[str, ...]]] = {
    1: (codon_dict, extended_codon_dict, byte_dict, ("TAA", "TAG", "TGA"), ("ATG",))
}

def _resolve_table_id(table: str) -> int:
    """Resolves an NCBI genetic code given as a numeric string or a table name to its ID."""
    if table.isdigit():
        return int(table)
    for t_id, info in NCBI_GENETIC_CODES.items():
        if info["name"].lower() == table.lower():
            return t_id
    raise ValueError(f"Unknown genetic code name: '{table}'. Available IDs: {list(NCBI_GENETIC_CODES.keys())}")

def get_genetic_code_tables(table: int | str | dict[str, str] = 1) -> tuple[dict[str, str], dict[str, str], defaultdict, tuple[str, ...], tuple[str, ...]]:
    """
    Returns (base_codon_dict, extended_codon_dict, byte_dict, default_stops, default_starts)
    for a given NCBI genetic code table ID, name, or custom dictionary.
    All tables are fully IUPAC-extended for ambiguous base calling.
    """
    if table is None:
        table = 1
    if isinstance(table, str):
        table = _resolve_table_id(table)

    if isinstance(table, int):
        if table in _GENETIC_CODE_CACHE:
            return _GENETIC_CODE_CACHE[table]
        if table not in NCBI_GENETIC_CODES:
            raise ValueError(f"Unsupported NCBI genetic code ID: {table}. Supported tables: {list(NCBI_GENETIC_CODES.keys())}")

        info = NCBI_GENETIC_CODES[table]
        base = codon_dict.copy()
        base.update(info["diff"])
        ext, b_dict = _build_extended_codon_table(base)
        _GENETIC_CODE_CACHE[table] = (base, ext, b_dict, info["stops"], info["starts"])
        return _GENETIC_CODE_CACHE[table]

    elif isinstance(table, dict):
        base = codon_dict.copy()
        base.update(table)
        stops = tuple(sorted([k for k, v in base.items() if v == "*"]))
        starts = ("ATG",)
        ext, b_dict = _build_extended_codon_table(base)
        return base, ext, b_dict, stops, starts
    else:
        raise TypeError(f"table must be an int, str, or dict, got {type(table).__name__}")

def get_start_codons(table: int | str | dict[str, str] = 1) -> tuple[str, ...]:
    """Curated (strict) start codons of an NCBI genetic code table (ID or name). Custom dictionary tables only use ATG."""
    if isinstance(table, dict):
        return ("ATG",)
    return get_genetic_code_tables(table)[4]

def get_alt_start_codons(table: int | str | dict[str, str] = 1) -> tuple[str, ...]:
    """NCBI initiator codons of a table beyond its curated start codons. Custom dictionary tables have none."""
    if isinstance(table, dict):
        return ()
    if table is None:
        table = 1
    if isinstance(table, str):
        table = _resolve_table_id(table)
    if table not in NCBI_GENETIC_CODES:
        raise ValueError(f"Unsupported NCBI genetic code ID: {table}. Supported tables: {list(NCBI_GENETIC_CODES.keys())}")
    return NCBI_GENETIC_CODES[table]["alt_starts"]

def translate(seq: str, table: int | str | dict[str, str] = 1) -> str:
    """
    Translates a DNA or RNA sequence to Amino Acids.
    Handles all ambiguous IUPAC bases according to the specified genetic code table.
    Case-insensitive.
    Assumes input is a multiple of 3.
    """
    if table == 1 or table == "1":
        b_dict = byte_dict
    else:
        _, _, b_dict, _, _ = get_genetic_code_tables(table)

    normalized_seq = seq.upper().replace("U", "T")
    it = iter(normalized_seq.encode('ascii'))
    return "".join(map(b_dict.__getitem__, zip(it, it, it)))

def map_relative_to_genomic(segments:list[Feature], rel_start:int, rel_end:int, strand:str):

    working_segments = segments if strand != "-" else reversed(segments)
    
    output_segments = []
    current_offset = 0
    
    for seg in working_segments:
        if current_offset > rel_end:
            break
            
        seg_len = seg.end - seg.start + 1
        
        seg_rel_start = current_offset
        seg_rel_end = current_offset + seg_len - 1
        
        overlap_start = max(rel_start, seg_rel_start)
        overlap_end = min(rel_end, seg_rel_end)
        
        if overlap_start <= overlap_end:
            dist_5p = overlap_start - seg_rel_start
            overlap_len = overlap_end - overlap_start + 1
            
            if strand != "-":
                g_start = seg.start + dist_5p
                g_end = g_start + overlap_len - 1
            else:
                g_end = seg.end - dist_5p
                g_start = g_end - overlap_len + 1
            
            output_segments.append((int(g_start), int(g_end)))
            
        current_offset += seg_len

    if strand == "-":
        output_segments.reverse()

    return output_segments

def find_ORFs(
    in_seq: str,
    must_have_stop: bool = True,
    tolerated_stops: Union[int, float, None] = 0,
    min_codon_len: int = 2,
    enforce_start_codon: bool = True,
    start_codons: tuple[str, ...] | None = None,
    stop_codons: tuple[str, ...] | None = None,
    table: int | str | dict[str, str] = 1,
) -> list[tuple[str, int, int]]:

    if tolerated_stops is None or tolerated_stops < 0:
        tolerated_stops = float('inf')

    if stop_codons is None or start_codons is None:
        _, _, _, def_stops, def_starts = get_genetic_code_tables(table)
        if stop_codons is None:
            stop_codons = def_stops
        if start_codons is None:
            start_codons = def_starts

    stop_set = frozenset(s.upper().replace("U", "T") for s in (stop_codons if isinstance(stop_codons, (set, frozenset, tuple, list)) else [stop_codons]))
    start_codons = tuple(s.upper().replace("U", "T") for s in start_codons)

    seq_len = len(in_seq)
    seq_upper = in_seq.upper().replace("U", "T")
    min_seq_len = min_codon_len * 3

    f0, f1, f2 = [], [], []
    appends = (f0.append, f1.append, f2.append)

    starts_set = set()
    limit_start = seq_len - 3

    if enforce_start_codon:
        for st_codon in start_codons:
            i = seq_upper.find(st_codon)
            while i != -1:
                if i <= limit_start:
                    starts_set.add(i)
                i = seq_upper.find(st_codon, i + 1)
    else:
        for init_idx in range(min(3, seq_len - 2)):
            starts_set.add(init_idx)
            
        for stop in stop_set:
            idx = seq_upper.find(stop)
            while idx != -1:
                if idx + 3 <= limit_start:
                    starts_set.add(idx + 3)
                idx = seq_upper.find(stop, idx + 1)
                
    starts = sorted(starts_set)
    limit_stop = seq_len - 2

    for i in starts:
        if not enforce_start_codon and seq_upper[i:i+3] in stop_set:
            continue
            
        append_func = appends[i % 3]
        stops_seen = 0
        last_stop_idx = -1
        
        for j in range(i + 3, limit_stop, 3):
            end_idx = j + 3
            codon = seq_upper[j:end_idx]
            
            if codon in stop_set:
                stops_seen += 1
                last_stop_idx = end_idx
                
                if stops_seen > tolerated_stops:
                    if end_idx - i >= min_seq_len:
                        append_func((in_seq[i:end_idx], i, end_idx - 1))
                    break
        else:
            if not must_have_stop:
                frame_end = i + ((seq_len - i) // 3) * 3
                if frame_end - i >= min_seq_len:
                    append_func((in_seq[i:frame_end], i, frame_end - 1))
            else:
                if stops_seen > 0 and last_stop_idx != -1:
                    if last_stop_idx - i >= min_seq_len:
                        append_func((in_seq[i:last_stop_idx], i, last_stop_idx - 1))

    return f0 + f1 + f2

def choose_orf(orfs: list[tuple[str, int, int]], mode: Literal["longest", "earliest"]="longest") -> tuple[str, int, int]:
    if mode == "longest":
        return max(orfs, key=lambda x: (len(x[0]), -x[1]), default=("", 0, -1))
        
    elif mode == "earliest":
        return max(orfs, key=lambda x: (-x[1], len(x[0])), default=("", 0, -1))
        
    else:
        raise ValueError(f"Invalid mode: '{mode}'. Expected 'longest' or 'earliest'.")

def trim_surplus(
    in_seq: str,
    mode: Literal["start", "end", "orf", "orf_or_end", "orf_or_start"] = "orf_or_end",
    max_nucleotide_trim: int | None = None,
    tolerated_stops: int | None = 0,
    orf_choice_mode: Literal["longest", "earliest"]="longest",
    must_have_stop: bool = True,
    enforce_start_codon: bool = True,
    start_codons: tuple[str, ...] | None = None,
    stop_codons: tuple[str, ...] | None = None,
    min_codon_len: int = 2,
    phase: int = 0,
    table: int | str | dict[str, str] = 1,
) -> tuple[str, bool, int, int]:
    """
    Trims surplus nucleotides to ensure sequence length is a multiple of 3, or extracts an ORF.
    
    in_seq: The input nucleotide sequence.
    mode: Trimming strategy. 
        - "start": Trims from the 5' end.
        - "end": Trims from the 3' end.
        - "orf": Extracts the best ORF, must return an ORF or nothing.
        - "orf_or_end": Extracts best ORF. Falls back to 3' trimming if criteria fail.
        - "orf_or_start": Extracts best ORF. Falls back to 5' trimming if criteria fail.
        - if "tolerated_stops" is negative or None, an infinite number will be tolerated (full readthrough mode)
    max_nucleotide_trim: Maximum allowed nucleotides to trim when using ORF modes.
    phase: Number of bases to skip at the 5' end (reading frame phase: 0, 1, or 2).
    table: NCBI genetic code table ID (default: 1 for Standard) or custom dictionary.
    """

    phase_offset = phase if phase in (1, 2) else 0
    available_len = max(0, len(in_seq) - phase_offset)
    surplus = available_len % 3
    nucleotide_surplus = (surplus != 0 or phase_offset != 0)
    coding_start = phase_offset
    coding_end = len(in_seq) - 1

    if mode == "end":
        if surplus:
            coding_end -= surplus
        out_seq = in_seq[coding_start:coding_end + 1] if coding_end >= coding_start else ""
    
    elif mode == "start":
        coding_start += surplus
        out_seq = in_seq[coding_start:coding_end + 1] if coding_end >= coding_start else ""

    elif mode in ("orf", "orf_or_end", "orf_or_start"):
        search_seq = in_seq[phase_offset:] if phase_offset else in_seq
        orfs = find_ORFs(
            search_seq,
            tolerated_stops=tolerated_stops,
            must_have_stop=must_have_stop,
            enforce_start_codon=enforce_start_codon,
            start_codons=start_codons,
            stop_codons=stop_codons,
            min_codon_len=min_codon_len,
            table=table,
        )
        orf, orf_start, orf_end = choose_orf(orfs, mode=orf_choice_mode)

        if orf and (max_nucleotide_trim is None or (len(search_seq) - len(orf)) <= max_nucleotide_trim):
            out_seq = orf
            coding_start = orf_start + phase_offset
            coding_end = orf_end + phase_offset
            nucleotide_surplus = False
        else:
            if mode == "orf_or_end":
                if surplus:
                    coding_end -= surplus
                out_seq = in_seq[coding_start:coding_end + 1] if coding_end >= coding_start else ""

            elif mode == "orf_or_start":
                coding_start += surplus
                out_seq = in_seq[coding_start:coding_end + 1] if coding_end >= coding_start else ""

            else: # mode == "orf"
                out_seq = ""
                nucleotide_surplus = False
                coding_start = 0
                coding_end = -1

    else:
        raise ValueError(f"Invalid mode: {mode}")

    return out_seq, nucleotide_surplus, coding_start, coding_end

def sort_and_update_genes(chrom:str, genes_dict:dict[str, Gene]) -> tuple[str, dict[str, Gene]]:
    genes = sorted(genes_dict.values())
    sorted_genes = {g.id: g for g in genes} 
    return chrom, sorted_genes

def export_group_equivalences(annotations:list[Annotation], output_folder:str|Path, group_tag:str="", synteny:bool=False, overlap_threshold:int=6, verbose:bool=True, clear_overlaps:bool=False, include_NAs:bool=False, output_also_single_files:bool=False, quiet:bool=False):
    """
    This generates equivalences between a set of annotation objects, whether only reporting equivalences to a particular target or between all annotations.
    """

    start = time.time()

    if synteny:
        column_sort_order = ["gene_id_A_origin", "gene_id_B_origin", "overlap_score", "gene_id_A_synteny_conserved", "gene_id_B_synteny_conserved", "gene_id_A", "gene_id_B"]
        ascending = [True, True, False, False, False, True, True]
    else:
        column_sort_order = ["gene_id_A_origin", "gene_id_B_origin", "overlap_score", "gene_id_A", "gene_id_B"]
        ascending = [True, True, False, True, True]

    genome_none = False
    genome_name = ""
    for a in annotations:
        if a.genome == None:
            genome_none = True
        else:
            genome_name = a.genome

    if genome_none:
        warnings.warn("Please verify that all annotations are associated to the same genome version/assembly, this could not be checked based on annotation files alone.", category=UserWarning)

    if genome_name != "":
        for a in annotations:
            if a.genome != None:
                if a.genome != genome_name:
                    raise ValueError("The provided annotations are not based on the same genome version/assembly. Please review input.")
                
    if len(annotations) < 2:
        raise ValueError(f"Not enough annotations ({annotations}) have been provided to export group equivalences.")

    if clear_overlaps:
        for a in annotations:
            a.overlaps.clear()

    export_folder = Path(output_folder) / "overlaps"
    export_folder.mkdir(parents=True, exist_ok=True)
    export_folder = str(export_folder) + "/"

    reference = ""

    for a in annotations:
        if a.target:
            reference = a.name
            break

    processed_pairs = set()
    
    for a1 in annotations:
        for a2 in annotations:
            if a1.name == a2.name:
                continue
                
            pair = tuple(sorted([a1.name, a2.name]))
            if pair in processed_pairs:
                continue

            if reference != "":
                if not a1.target and not a2.target:
                    continue

            a1.overlaps.detect(a2, clear=False)

            processed_pairs.add(pair)

    all_genes = {}
    unmapped_genes = {}

    if include_NAs:
        for a in annotations:
            all_genes[a.name] = set(a.all_gene_ids.keys())
            unmapped_genes[a.name] = set(a.unmapped)

    for x, a in enumerate(annotations):

        if group_tag and (reference or len(annotations) == 2):
            prefix = group_tag

        elif reference and len(annotations) == 2:
            for o in annotations:
                if not o.target:
                    other = o.name
                    break
            prefix = f"{reference}_{other}"
        
        elif len(annotations) == 2:
            prefix = f"{annotations[0].name}_{annotations[1].name}"
        
        else:
            prefix = a.name

        if genome_name:
            single_tag = f"{prefix}_on_{genome_name}_overlaps_t{overlap_threshold}.csv"
        else:
            single_tag = f"{prefix}_overlaps_t{overlap_threshold}.csv"

        if reference:
            if a.name != reference:
                continue

        elif len(annotations) == 2:
            if x != 0:
                continue

        single_df = a.overlaps.export(overlap_threshold=overlap_threshold, synteny=synteny, verbose=verbose, NAs=False, save_csv=False)

        if len(annotations) > 2:

            if x == 0:
                eq_df = single_df.copy()
            else:
                eq_df = pd.concat([eq_df, single_df])

        if len(annotations) == 2 or output_also_single_files:

            if include_NAs:
                na_rows = []

                for a_name, genes in all_genes.items():

                    temp_df = single_df[single_df["gene_id_A_origin"] == a_name]
                    present = set(temp_df["gene_id_A"].dropna())
                    temp_df = single_df[single_df["gene_id_B_origin"] == a_name]
                    present = present | set(temp_df["gene_id_B"].dropna())

                    if a_name == a.name:
                        for g in genes:
                            if g not in present:
                                na_rows.append({
                                    "gene_id_A": g,
                                    "gene_id_A_origin": a_name,
                                    "overlap_score": 0
                                })

                    else:
                        for g in genes:
                            if g not in present:
                                na_rows.append({
                                    "gene_id_B": g,
                                    "gene_id_B_origin": a_name,
                                    "overlap_score": 0
                                })

                if synteny:
                    for a_name, unmapped in unmapped_genes.items():

                        if a_name == a.name:

                            for g_id in unmapped:
                                na_rows.append({
                                    "gene_id_A": g_id,
                                    "gene_id_A_origin": a_name
                                })
                        else:
                            for g_id in unmapped:
                                na_rows.append({
                                    "gene_id_B": g_id,
                                    "gene_id_B_origin": a_name
                                })

                if na_rows:
                    single_df = pd.concat([single_df, pd.DataFrame(na_rows)], ignore_index=True)

            single_df.sort_values(by=column_sort_order, ascending=ascending, inplace=True)
            single_df.reset_index(drop=True, inplace=True)
            single_df.to_csv(f"{export_folder}{single_tag}", sep="\t", index=False, na_rep="NA")

    if len(annotations) > 2:
        if group_tag:
            prefix = group_tag
        else:
            prefix = f"{annotations[0].name}...{annotations[-1].name}"

        if genome_name:
            tag = f"{prefix}_on_{genome_name}_overlaps_t{overlap_threshold}.csv"
        else:
            tag = f"{prefix}_overlaps_t{overlap_threshold}.csv"

        if include_NAs:

            na_rows = []

            for a_name, genes in all_genes.items():

                temp_df = eq_df[eq_df["gene_id_A_origin"] == a_name]
                present = set(temp_df["gene_id_A"].dropna())
                temp_df = eq_df[eq_df["gene_id_B_origin"] == a_name]
                present = present | set(temp_df["gene_id_B"].dropna())

                for g in genes:
                    if g not in present:
                        na_rows.append({
                            "gene_id_A": g,
                            "gene_id_A_origin": a_name,
                            "overlap_score": 0
                        })

            if synteny:
                for a_name, unmapped in unmapped_genes.items():
                    for g_id in unmapped:
                        na_rows.append({
                            "gene_id_A": g_id,
                            "gene_id_A_origin": a_name
                        })

            if na_rows:
                eq_df = pd.concat([eq_df, pd.DataFrame(na_rows)], ignore_index=True)

        eq_df.sort_values(by=column_sort_order, ascending=ascending, inplace=True)
        eq_df.reset_index(drop=True, inplace=True)
        eq_df.to_csv(f"{export_folder}{tag}", sep="\t", index=False, na_rep="NA")

    now = time.time()
    lapse = now - start
    if not quiet:
        print(f"\nGenerating overlaps for annotations = '{annotations}' took {round(lapse/60, 1)} minutes\n")

        