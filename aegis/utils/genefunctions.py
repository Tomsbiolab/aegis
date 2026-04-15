"""
Created on Tue Dec 27 15:20:59 2022

Module with an array of genomic functions.

@authors: David Navarro, Antonio Santiago
"""
from __future__ import annotations
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from ..annotation import Annotation
    from ..gene import Gene

import pandas as pd
import time
import warnings
import itertools

from collections import defaultdict
from pathlib import Path

# nucleotides
_STR_FROM = "ACGTRYSWKMBDHVNX-"
_STR_TO   = "TGCAYRSWMKVHDBNX-"
_BYTES_COMP_TABLE = bytes.maketrans(_STR_FROM.encode(), _STR_TO.encode())

def reverse_complement(in_seq: str) -> str:
    return in_seq.encode('ascii').translate(_BYTES_COMP_TABLE)[::-1].decode('ascii')

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


extended_codon_dict = {}
all_iupac_chars = list(iupac_dna_nucleotides.keys())

for c1, c2, c3 in itertools.product(all_iupac_chars, repeat=3):
    ambiguous_codon = c1 + c2 + c3
    
    possible_aas = set()
    for b1, b2, b3 in itertools.product(iupac_dna_nucleotides[c1], iupac_dna_nucleotides[c2], iupac_dna_nucleotides[c3]):
        standard_codon = b1 + b2 + b3
        possible_aas.add(codon_dict[standard_codon])
    
    if len(possible_aas) == 1:
        extended_codon_dict[ambiguous_codon] = possible_aas.pop()
    else:
        extended_codon_dict[ambiguous_codon] = "X"

byte_codon_dict = {tuple(k.encode('ascii')): v for k, v in extended_codon_dict.items()}
byte_dict = defaultdict(lambda: "X", byte_codon_dict)

def translate(seq: str) -> str:
    """
    Translates an uppercase DNA sequence to Amino Acids.
    Handles all ambiguous IUPAC bases.
    Assumes input is a multiple of 3
    """

    it = iter(seq.encode('ascii'))

    return "".join(map(byte_dict.__getitem__, zip(it, it, it)))

def find_ORFs(in_seq: str, must_have_stop: bool = True, readthrough_stop: bool = False, min_codon_len: int = 2, start_codon: str = "ATG", stop_codons=["TAA", "TAG", "TGA"]) -> list[tuple[str, int, int]]:
    
    stop_set = frozenset(stop_codons) if not isinstance(stop_codons, (set, frozenset)) else stop_codons
    seq_len = len(in_seq)
    min_seq_len = min_codon_len * 3
    
    f0, f1, f2 = [], [], []
    appends = (f0.append, f1.append, f2.append)

    starts =[]
    starts_append = starts.append
    i = in_seq.find(start_codon)
    limit_start = seq_len - 3
    
    while i != -1:
        if i <= limit_start:
            starts_append(i)
        i = in_seq.find(start_codon, i + 1)

    append_start = (not must_have_stop) and (3 >= min_seq_len)
    limit_stop = seq_len - 2

    if must_have_stop and not readthrough_stop:
        for i in starts:
            append_func = appends[i % 3]
            for j in range(i + 3, limit_stop, 3):
                if in_seq[j:j+3] in stop_set:
                    if j + 3 - i >= min_seq_len:
                        append_func((in_seq[i:j+3], i, j + 2))
                    break
                    
    elif must_have_stop and readthrough_stop:
        for i in starts:
            append_func = appends[i % 3]
            for j in range(i + 3, limit_stop, 3):
                if in_seq[j:j+3] in stop_set:
                    if j + 3 - i >= min_seq_len:
                        append_func((in_seq[i:j+3], i, j + 2))
                        
    elif not must_have_stop and not readthrough_stop:
        for i in starts:
            append_func = appends[i % 3]
            if append_start:
                append_func((in_seq[i:i+3], i, i + 2))
            for j in range(i + 3, limit_stop, 3):
                if j + 3 - i >= min_seq_len:
                    append_func((in_seq[i:j+3], i, j + 2))
                if in_seq[j:j+3] in stop_set:
                    break
                    
    else:
        for i in starts:
            append_func = appends[i % 3]
            if append_start:
                append_func((in_seq[i:i+3], i, i + 2))
            for j in range(i + 3, limit_stop, 3):
                if j + 3 - i >= min_seq_len:
                    append_func((in_seq[i:j+3], i, j + 2))

    return f0 + f1 + f2

def longest_ORF(orfs: list[tuple[str, int, int]]) -> tuple[str, int, int]:
    return max(orfs, key=lambda x: len(x[0]), default=("", 0, 0))

def trim_surplus(in_seq: str, mode: Literal["start", "end", "orf", "orf_or_end"] = "orf_or_end", max_nucleotide_trim: int | None = 6, readthrough_stop: bool = True) -> tuple[str, bool]:
    """
    Trims surplus nucleotides to ensure sequence length is a multiple of 3, or extracts an ORF.
    
    in_seq: The input nucleotide sequence.
    mode: Trimming strategy. 
        - "start": Trims from the 5' end.
        - "end": Trims from the 3' end.
        - "orf": Extracts the longest ORF. Falls back to original sequence if criteria fail.
        - "orf_or_end": Extracts longest ORF. Falls back to 3' trimming if criteria fail.
    max_nucleotide_trim: Maximum allowed nucleotides to trim when using ORF modes.
    readthrough_stop: Passed to find_ORFs.
    """

    surplus = len(in_seq) % 3
    nucleotide_surplus = surplus != 0

    if mode == "end":
        out_seq = in_seq[:-surplus] if surplus else in_seq
    
    elif mode == "start":
        out_seq = in_seq[surplus:]

    elif mode in ("orf", "orf_or_end"):
        orfs = find_ORFs(in_seq, readthrough_stop=readthrough_stop)
        orf, _, _ = longest_ORF(orfs)

        if orf and (max_nucleotide_trim is None or (len(in_seq) - len(orf)) <= max_nucleotide_trim):
            out_seq = orf
        else:
            if mode == "orf_or_end":
                out_seq = in_seq[:-surplus] if surplus else in_seq
            else: # mode == "orf"
                out_seq = in_seq
                
    else:
        raise ValueError(f"Invalid mode: {mode}")

    return out_seq, nucleotide_surplus

def flexible_translate(in_seq:str, codon_dict:dict[str, str]=codon_dict, readthrough:str="both", must_have_stop:bool=True):
    """
    Translates a nucleotide sequence into a protein sequence.
    
    in_seq: The input nucleotide sequence.
    readthrough: Translation strategy. 
        - "both": Translates the entire sequence as is.
        - "start": Translates from the first codon even if it is not an ATG.
        - "end": Translates from an ATG and readsthrough stop codons.
    must_have_stop: Whether the protein must have a stop codon.
    codon_table: The codon table to use for translation, standard genetic code by default.
    """

    out_seq = ""
    in_seq = in_seq.upper()
    start = "present"
    early_stop = False
    end_stop = True
    gaps = False
    coding_start = False
    coding_end = False
 

    if readthrough == "both" or readthrough == "start" or readthrough == "end":
        in_seq, nucleotide_surplus = trim_surplus(in_seq)
    else:
        nucleotide_surplus = False

    ambiguous_letters = ["B", "D", "H", "K", "M", "N", "R", "S", "V", "W", "Y", "X"]
    # for masked genomes
    ambiguous_letters.append("X")

    # default option where CDS sequence is read as it has been annotated
    if readthrough == "both":
        for x in range(0, len(in_seq), 3):
            temp_codon = in_seq[x] + in_seq[x+1] + in_seq[x+2]
            if any(letter in ambiguous_letters for letter in temp_codon):
                amino = "X"
                gaps = True
            else:
                amino = codon_dict[temp_codon]
                out_seq += amino
            if x == 0 and amino != "M":
                start = "absent"
                if "ATG" in in_seq:
                    start = "late"
            elif x == (len(in_seq) - 3) and amino != "*":
                end_stop = False
            if amino == "*" and x < (len(in_seq) - 3):
                early_stop = True

    # the protein will stop being read after the first
    # stop codon, but the first ATG will not be searched for, i.e. the first
    # codon is readthrough no matter what
    elif readthrough == "start":
        for x in range(0, len(in_seq), 3):
            temp_codon = in_seq[x] + in_seq[x+1] + in_seq[x+2]
            if any(letter in ambiguous_letters for letter in temp_codon):
                amino = "X"
                gaps = True
            else:
                amino = codon_dict[temp_codon]
                out_seq += amino
            if x == 0 and amino != "M":
                start = "absent"
                if "ATG" in in_seq:
                    start = "late"
            elif x == (len(in_seq) - 3) and amino != "*":
                end_stop = False
            if amino == "*" and x < (len(in_seq) - 3):
                early_stop = True
                break

    # this begins reading from the first ATG till the end of the sequence
    elif readthrough == "end":
        index = in_seq.find("ATG")
        if index == -1:
            # ATG codon not found and hence output is empty, however end_stop
            # and early_stop, or gaps are determined
            out_seq = ""
            for x in range(0, len(in_seq), 3):
                temp_codon = in_seq[x] + in_seq[x+1] + in_seq[x+2]
                if any(letter in ambiguous_letters for letter in temp_codon):
                    amino = "X"
                    gaps = True
                if x == (len(in_seq) - 3) and amino != "*":
                    end_stop = False
                if amino == "*" and x < (len(in_seq) - 3):
                    early_stop = True
        else:
            # Trim nucleotides preceding ATG codon
            return in_seq[index:]

    # longest orf is read as a protein, this is with readthrough == "none"
    else:
        orfs = find_ORFs(in_seq, must_have_stop=must_have_stop)
        if len(orfs) > 0:
            longest, coding_start, coding_end = longest_ORF(orfs)

            for x in range(0, len(longest), 3):
                temp_codon = longest[x] + longest[x+1] + longest[x+2]
                if any(letter in ambiguous_letters for letter in temp_codon):
                    amino = "X"
                    gaps = True
                else:
                    amino = codon_dict[temp_codon]
                    out_seq += amino
                if x == 0 and amino != "M":
                    start = "absent"
                    if "ATG" in longest:
                        start = "late"
                elif x == (len(in_seq) - 3) and amino != "*":
                    end_stop = False
                if amino == "*" and x < (len(in_seq) - 3):
                    early_stop = True
                    break
        else:
            start = "absent"
            end_stop = False
            out_seq = ""

    return start, end_stop, early_stop, nucleotide_surplus, gaps, out_seq, coding_start, coding_end

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

        single_df = a.overlaps.export(overlap_threshold=overlap_threshold, synteny=synteny, verbose=verbose, NAs=False)

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

        