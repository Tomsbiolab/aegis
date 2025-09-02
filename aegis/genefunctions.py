#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 15:20:59 2022

Module with an array of genomic functions.

@authors: David Navarro, Antonio Santiago
"""

from pathlib import Path
from collections import Counter
import pandas as pd
from Bio.Data import CodonTable
from Bio.Seq import Seq
import time
import re
import warnings

def count_occurrences(string, char):
    return Counter(string)[char]


def find_all_occurrences(pattern, text):
    matches = []
    for match in re.finditer(pattern, text):
        matches.append((match.start(), match.end(), match.group()))

    return matches


def reverse_complement(in_seq:str):
    in_seq = Seq(in_seq)
    out_seq = str(in_seq.reverse_complement())
        
    return out_seq


def find_ORFs(in_seq:str, must_have_stop=True, readthrough_stop=False):
    orfs = []
    start_codon = "ATG"
    stop_codons = ["TAA", "TAG", "TGA"]
    for frame in range(3):
        for i in range(frame, len(in_seq)-2, 3):
            codon = in_seq[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(in_seq), 3):
                    codon2 = in_seq[j:j+3]
                    orf = in_seq[i:j+3], i, j + 2
                    if must_have_stop:
                        if codon2 in stop_codons:
                            if len(orf[0]) % 3 == 0:
                                orfs.append(orf)
                    else:
                        if len(orf[0]) % 3 == 0:
                            orfs.append(orf)
                    if codon2 in stop_codons:
                        if not readthrough_stop:
                            break
    return orfs

def longest_ORF(orfs:list):
    longest = ("", 0, 0) #dumb ORF
    for orf in orfs:
        if len(orf[0]) > len(longest[0]):
            longest = orf

    return longest

def trim_surplus(in_seq:str):
    """
    Function that trims surplus nucleotides in the event of a sequence not
    being a multiple of three. The trimming is orientated by the presence of
    start and end codons provided they are very close to the start and end
    of input sequences, otherwise just the end of the sequence is trimmed.
    """
    nucleotide_surplus = False
    surplus = len(in_seq) % 3
    if surplus != 0:
        nucleotide_surplus = True
        orfs = find_ORFs(in_seq, readthrough_stop=True)
        orf, _, _ = longest_ORF(orfs)
        excess = len(in_seq) - len(orf)
        # this aims to trim at most 2 complete or incomplete codons: one from
        # the beginning of the sequence and one from the end
        if excess < 6:
            out_seq = orf
        else:
            out_seq = in_seq[:-surplus]
    else:
        out_seq = in_seq

    return out_seq, nucleotide_surplus

def translate(in_seq:str, readthrough:str="both", must_have_stop:bool=True,
              # standard genetic code
              codon_table:CodonTable=CodonTable.unambiguous_dna_by_id[1]):
    # translating a protein
    nucleotide_surplus = 0
    out_seq = ""
    in_seq = in_seq.upper()
    start = "present"
    early_stop = False
    end_stop = True
    gaps = False
    coding_start = False
    coding_end = False
    
    codon_dict = {}
    for codon, aa in codon_table.forward_table.items():
        codon_dict[codon] = aa

    codon_dict["TAA"] = "*"
    codon_dict["TAG"] = "*"
    codon_dict["TGA"] = "*"

    if readthrough == "both" or readthrough == "start" or readthrough == "end":
        in_seq, nucleotide_surplus = trim_surplus(in_seq)

    ambiguous_letters = ["B", "D", "H", "K", "M", "N", "R", "S", "V", "W", "Y"]
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

def overlap(feat1, feat2):
    overlapping = False
    interval1 = feat1.size
    interval2 = feat2.size
    small = min(feat1.start, feat1.end, feat2.start, feat2.end)
    large = max(feat1.start, feat1.end, feat2.start, feat2.end)
    overlap_bp = (interval1 + interval2) - ((large - small) + 1)
    # checking only overlapping features
    if overlap_bp > 0:
        overlapping = True

    return overlapping, overlap_bp

def export_for_dapseq(annotation, genome, chromosome_dictionary:dict={}, genome_out_folder:str="", gff_out_folder:str="", tag:str="_for_dap.gff3", skip_atypical_fts:bool=True, main_only:bool=False, UTRs:bool=False, exclude_non_coding:bool=False):
    equivalences = genome.rename_features_dap(custom_path=genome_out_folder, return_equivalences=True, export=True, chromosome_dictionary=chromosome_dictionary)
    annotation.rename_chromosomes(equivalences)
    annotation.export_gff(output_folder=gff_out_folder, tag=tag, skip_atypical_fts=skip_atypical_fts, main_only=main_only, UTRs=UTRs, exclude_non_coding=exclude_non_coding)

def export_group_equivalences(annotations:list, output_folder, group_tag:str="", synteny:bool=False, overlap_threshold:int=6, verbose:bool=True, clear_overlaps=False, include_NAs=False, output_also_single_files=False, quiet:bool=False):
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
            genome_name = a.genome.name

    if genome_none:
        warnings.warn("Please verify that all annotations are associated to the same genome version/assembly, this could not be checked based on annotation files alone.", category=UserWarning)

    if genome_name != "":
        for a in annotations:
            if a.genome != None:
                if a.genome.name != genome_name:
                    raise ValueError("The provided annotations are not based on the same genome version/assembly. Please review input.")
                
    if len(annotations) < 2:
        raise ValueError(f"Not enough annotations ({annotations}) have been provided to export group equivalences.")

    if clear_overlaps:
        for a in annotations:
            a.clear_overlaps()

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

            a1.detect_gene_overlaps(a2, clear=False)

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

        single_df = a.export_equivalences(overlap_threshold=overlap_threshold, synteny=synteny, verbose=verbose, NAs=False)

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

        