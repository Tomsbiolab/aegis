#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 27 15:20:59 2022

Module with an array of genomic functions.

@authors: David Navarro, Antonio Santiago
"""

from collections import Counter
from os import system
import pandas as pd
from Bio.Data import CodonTable
from Bio.Seq import Seq
import time
import re

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

def find_ORFs_old(in_seq:str, stop_codon_within_orf=False):
    """
    This function identifies and returns all open reading frames (ORFs) 
    in a given nucleotide sequence. An ORF is by default only considered to
    exist if it begins with a start codon and ends with a stop codon and has
    no stop codons inbetween. An option is included to get ORFs even when they
    have stop codons within the sequence.
    """
    orfs = []
    start_codon = "ATG"
    stop_codons = ['TAA', 'TAG', 'TGA']
    for frame in range(3):
        for i in range(frame, len(in_seq)-2, 3):
            codon = in_seq[i:i+3]
            if codon == start_codon:
                for j in range(i+3, len(in_seq), 3):
                    codon = in_seq[j:j+3]
                    if codon in stop_codons:
                        if not stop_codon_within_orf:
                            break
                    else:
                        orf = in_seq[i:j+3]
                        if len(orf) % 3 == 0:
                            orfs.append((orf, i, j + 2))
                            
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
        orfs = find_ORFs_old(in_seq, stop_codon_within_orf=True)
        orf, coding_start, coding_end = longest_ORF(orfs)
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

def export_group_equivalences(annotations:list, custom_path:str="", multidirectional:bool=False, group_tag:str="all", stringent:bool=True, verbose:bool=True, clear_overlaps=True):
    """
    This generates equivalences between a set of annotation objects, whether only reporting equivalences to a particular target or between all annotations.
    """
    
    columns = ["query_origin", "target_origin", "overlap_score", "query_synteny_conserved", "target_synteny_conserved", "query_id", "target_id"]
    ascending = [True, True, False, False, False, True, True]

    # Force clear_overlaps before doing any overlaps is on by default to be safe,
    # but turning it off saves time if you know what you are doing

    annotations_to_overlap = []
    for a in annotations:
        if a.to_overlap:
            annotations_to_overlap.append(a.copy())

    del annotations

    genome = annotations_to_overlap[0].genome
    same_genome = True
    none_genome = False
    different_annotations = True
    annotation_ids = []

    for a in annotations_to_overlap:
        if a.id not in annotation_ids:
            annotation_ids.append(a.id)
        else:
            different_annotations = False

        if a.genome == None or genome == None:
            none_genome = True
            continue
        if a.genome != genome:
            same_genome = False

    if clear_overlaps:
        for a in annotations_to_overlap:
            a.clear_overlaps()


    if same_genome:

        if none_genome:
            print(f"Warning: Please make sure all submitted annotations are associated to the same genome version, this could not be verified as at least one of the submitted annotations has a 'None' genome.")

        if not multidirectional:
            print("Non multidirectional or targeted mode is not working yet.")

            # count = 0
            # for a in annotations_to_overlap:
            #     if a.target:
            #         count += 1
            # if count == 1:
            #     print(f"\nRunning 'export_group_equivalences' with 'stringent={stringent}' and 'verbose={verbose}' in target mode\n")
            #     for a1 in annotations_to_overlap:
            #         if a1.target:
            #             for a2 in annotations_to_overlap:
            #                 if not a2.to_overlap:
            #                     continue
            #                 if a1.name != a2.name:
            #                     a1.detect_gene_overlaps(a2, clear=False)
            #     for a1 in annotations_to_overlap:
            #         if a1.target:
            #             tag = f"{a1.id}_equivalences"
            #             overlapped_annotations = a1.overlapped_annotations
            #             if stringent:
            #                 tag += "_filtered"
            #             if not verbose:
            #                 tag += "_simple"
            #             if custom_path == "":
            #                 export_path = a1.path + "overlaps/"
            #             else:
            #                 export_path = custom_path
            #             if export_path[-1] != "/":
            #                 export_path += "/"
            #             system(f"mkdir -p {export_path}")
            #             eq_df = a1.export_equivalences(custom_path, stringent, verbose)
            #             start = time.time()
            #             overlapping_genes = eq_df["query_id"].to_list()
            #             for a2 in annotations_to_overlap:
            #                 if not a2.to_overlap:
            #                     continue
            #                 if a1.name != a2.name:
            #                     # adding query genes with no equivalence:
            #                     for genes in a2.chrs.values():
            #                         for g in genes.values():
            #                             if g.id not in overlapping_genes:
            #                                 index = len(eq_df.index)
            #                                 eq_df.loc[index] = pd.NA
            #                                 eq_df.loc[index, "query_id"] = g.id
            #                                 if a2.liftoff:
            #                                     eq_df.loc[index, "query_origin"] = a2.name
            #                                 else:
            #                                     eq_df.loc[index, "query_origin"] = a2.name
            #                                 eq_df.loc[index, "overlap_score"] = 0
            #                     # Only for liftoff annotations
            #                     for g_id in a2.unmapped:
            #                         if g_id not in overlapping_genes:
            #                             index = len(eq_df.index)
            #                             eq_df.loc[index] = pd.NA
            #                             eq_df.loc[index, "query_id"] = g_id
            #                             eq_df.loc[index, "query_origin"] = a2.name
            #     # There should not be any duplicates but just in case
            #     eq_df.drop_duplicates(inplace=True)
            #     eq_df.sort_values(by=columns, ascending=ascending, inplace=True)
            #     eq_df.reset_index(drop=True, inplace=True)
            #     eq_df.to_csv(f"{export_path}{tag}.csv", sep="\t", index=False, na_rep="NA")
            #     now = time.time()
            #     lapse = now - start
            #     print(f"Appending query NAs from {overlapped_annotations} to the exported {a1.id} equivalences took {round(lapse/60, 1)} minutes\n")

            # else:
            #     print(f"\nERROR: Make one annotation the target instead of {count} or set it to multidirectional mode if you want to obtain equivalences between all annotations\n")

        else:
            if custom_path == "":
                print(f"\nERROR: please define an output folder to save the combined csv for {group_tag} annotations\n")
            else:
                print(f"\nRunning 'export_group_equivalences' with 'stringent={stringent}' and 'verbose={verbose}' in multidirectional mode for {group_tag} annotations\n")
                
                start = time.time()
                
                processed_pairs = []
                processed_annotations = []

                for a1 in annotations_to_overlap:
                    processed_annotations.append(a1.id)
                    genome = a1.genome
                    for a2 in annotations_to_overlap:
                        if a1.id != a2.id:
                            if f"{a1.id}-pair-{a2.id}" in processed_pairs:
                                continue
                            if f"{a2.id}-pair-{a1.id}" in processed_pairs:
                                continue
                            a1.detect_gene_overlaps(a2, clear=False)
                            processed_pairs.append(f"{a1.id}-pair-{a2.id}")
                            processed_pairs.append(f"{a2.id}-pair-{a1.id}")

                if genome != None:
                    tag = f"{group_tag}_on_{genome}_equivalences"
                else:
                    tag = f"{group_tag}_equivalences"

                
                if stringent:
                    tag += "_filtered"
                if not verbose:
                    tag += "_simple"
                
                export_path = custom_path
                    
                if export_path[-1] != "/":
                    export_path += "/"

                system(f"mkdir -p {export_path}")

                all_genes = {}

                for x, a1 in enumerate(annotations_to_overlap):
                    single_tag = f"{a1.id}_equivalences"
                    overlapped_annotations = a1.overlapped_annotations
                    if stringent:
                        single_tag += "_filtered"
                    if not verbose:
                        single_tag += "_simple"

                    single_df = a1.export_equivalences(custom_path, stringent, verbose)

                    all_genes[a1.id] = [single_df["query_id"].to_list(), a1.name]

                    # this removes unmatched genes only between query and target annotations
                    # unmapped genes during liftoff (overlap_score == NA) will never have a match even between target annotations
                    single_df = single_df[single_df["overlap_score"] != 0]

                    if x == 0:
                        eq_df = single_df.copy()
                    else:
                        eq_df = pd.concat([eq_df, single_df])

                for a_id in all_genes:
                    genes = all_genes[a_id][0]
                    name = all_genes[a_id][1]

                    temp_df = eq_df[eq_df["query_origin"] == name]
                    present_genes = temp_df["query_id"].to_list()

                    for g in genes:
                        if g not in present_genes:
                            index = len(eq_df.index)
                            eq_df.loc[index] = pd.NA
                            eq_df.loc[index, "query_id"] = g
                            eq_df.loc[index, "query_origin"] = name

                # There should not be any duplicates but just in case
                eq_df.drop_duplicates(inplace=True)
                eq_df.sort_values(by=columns, ascending=ascending, inplace=True)
                eq_df.reset_index(drop=True, inplace=True)
                eq_df.to_csv(f"{export_path}{tag}.csv", sep="\t", index=False, na_rep="NA")

                now = time.time()
                lapse = now - start
                print(f"\nGenerating multidirectional equivalences for all {processed_annotations} took {round(lapse/60, 1)} minutes\n")
    else:
        print("The annotations submitted do not share a common genome. Please try again.")