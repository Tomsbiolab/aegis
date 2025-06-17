#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on June 04 2025

@author: David Navarro
"""

import os

"""Main Program"""
# INPUT
root = "../../genomes_and_annotation/"

root = f"{os.path.abspath(root)}/"

output_folder = root + "aegis_output/annotation_transfers/"
feature_folder = output_folder + "liftoff_feature_types/"
unmapped_folder = output_folder + "unmapped/"
os.system(f"mkdir -p {unmapped_folder}")
os.system(f"mkdir -p {output_folder}original/")

gff_files = ["cannabis/aegis_output/gffs/CBDRx_on_CBDRx_aegis_fcounts.gff3", "cannabis/aegis_output/gffs/PurpleKush_on_PurpleKush_aegis_fcounts.gff3", "cannabis/aegis_output/gffs/Finola_on_Finola_aegis_fcounts.gff3", "cannabis/aegis_output/gffs/JamaicanLion_on_JamaicanLion_aegis_fcounts.gff3", "grapevine/PN40024/aegis_output/gffs/5.1_on_T2T_ref_aegis_fcounts_dapfit.gff3"]

lifton_gff_files = ["cannabis/aegis_output/gffs/CBDRx_on_CBDRx_aegis_fcounts_for_lifton.gff3", "cannabis/aegis_output/gffs/PurpleKush_on_PurpleKush_aegis_fcounts_for_lifton.gff3", "cannabis/aegis_output/gffs/Finola_on_Finola_aegis_fcounts_for_lifton.gff3", "cannabis/aegis_output/gffs/JamaicanLion_on_JamaicanLion_aegis_fcounts_for_lifton.gff3", "grapevine/PN40024/aegis_output/gffs/5.1_on_T2T_ref_aegis_fcounts_dapfit_for_lifton.gff3"]

genome_files = ["cannabis/aegis_output/genomes/CBDRx_organellefree.fasta", "cannabis/aegis_output/genomes/PurpleKush_organellefree.fasta", "cannabis/aegis_output/genomes/Finola_organellefree.fasta", "cannabis/aegis_output/genomes/JamaicanLion_organellefree.fasta", "grapevine/PN40024/aegis_output/genomes/T2T_ref_dapfit_organellefree_chr00.fasta"]

tags = ["CBDRx", "PurpleKush", "Finola", "JamaicanLion", "5.1"]

if len(gff_files) != len(genome_files) or len(genome_files) != len(tags) or len(tags) != len(lifton_gff_files):
    raise ValueError(f"Input gffs, genomes, or tags do not match.")

indispensable_members = ["5.1"]

liftoff = True
lifton = True

"""
Exploring liftoff options:
    
typical command:

liftoff -g input_gff3 -o output_gff3 -f types_file target_genome_fasta original_genome_fasta
        
but there are extra options
"""

# change to empty string if you want to use the native envitonment's lifton installation
singularity_prefix = ""

options = ["", "-copies", "-copies -flank 0.1"]
liftoff_tags = ["default_liftoff", "copies_liftoff", "copies_flank_liftoff"]
lifton_tags = ["default_lifton", "copies_lifton", "copies_flank_lifton"]

options = ["-copies -flank 0.1"]
liftoff_tags = ["copies_flank_liftoff"]
lifton_tags = ["copies_flank_lifton"]

cores = 8

ext = ".gff3"

for n, file in enumerate(gff_files):

    for n2, file2 in enumerate(gff_files):
        if file == file2:
            continue

        if indispensable_members != []:
            indispensable_indexes = []

            for member in indispensable_members:
                indispensable_indexes.append(tags.index(member))

            # this excludes certain liftoff/liftoff combinations
            if n not in indispensable_indexes and n2 not in indispensable_indexes:
                continue

        if lifton and liftoff:
            print(f"Processing liftoffs and liftons for {tags[n]} to {tags[n2]}")
        elif lifton:
            print(f"Processing liftons for {tags[n]} to {tags[n2]}")
        else:
            print(f"Processing liftoffs for {tags[n]} to {tags[n2]}")

        for o, opt in enumerate(options):

            if liftoff:

                liftoff_tag = f"{tags[n]}_{liftoff_tags[o]}_{tags[n2]}"
                liftoff_file = f"{liftoff_tag}_raw{ext}"

                liftoff_file_path = f"{output_folder}original/{liftoff_file}"

                unmapped_file = f"{unmapped_folder}{liftoff_tag}_unmapped.txt"

                intermediate_path = f"{output_folder}intermediate_files"

                if os.path.isdir(intermediate_path):
                    os.system(f"rm -r {intermediate_path}")

                if not os.path.isfile(liftoff_file_path):
                    print(f"Running a {liftoff_tags[o]} of {tags[n]} against {tags[n2]} target genome\n")
                    
                    command = (f"liftoff {opt} -g {root}{file} -o {liftoff_file_path} -u {unmapped_file} -dir {intermediate_path} -f {feature_folder}{tags[n]}_types.txt {root}{genome_files[n2]} {root}{genome_files[n]}")
                
                    os.system(command)

                    if os.path.isdir(intermediate_path):
                        os.system(f"rm -r {intermediate_path}")

                if not os.path.isfile(f"{output_folder}{liftoff_tag}.gff3"):
                    f = open(f"{liftoff_file_path}", "r", encoding="utf-8")
                    header = []
                    out = []
                    for line in f:
                        if not line:
                            continue
                        if line.startswith("#"):
                            if line == "###\n":
                                out.append(line)
                            elif line.startswith("##gff") or line.startswith("# Liftoff"):
                                if line not in header:
                                    header.append(line)
                            continue
                        
                        temp = line.strip().split("\t")
                        if len(temp) <= 2:
                            continue

                        out.append(line)

                    f.close()

                    out = header + out

                    out = "".join(out)

                    f = open(f"{output_folder}{liftoff_tag}.gff3", "w", encoding="utf-8")
                    f.write(out)
                    f.close()

            if lifton:

                lifton_tag = f"{tags[n]}_{lifton_tags[o]}_{tags[n2]}"

                lifton_file = f"{lifton_tag}{ext}"

                lifton_file_path = f"{output_folder}{lifton_file}"

                unmapped_file = f"{unmapped_folder}{lifton_tag}_unmapped.txt"

                lifton_output_path = f"{output_folder}lifton_output/"

                if os.path.isdir(lifton_output_path):
                    os.system(f"rm -r {lifton_output_path}")

                if not os.path.isfile(lifton_file_path):

                    print(f"Running a {lifton_tags[o]} of {tags[n]} against {tags[n2]} target genome\n")
                    
                    command = (f"{singularity_prefix}lifton {opt} -t {cores} -g {root}{lifton_gff_files[n]} -o {lifton_file_path} -u {unmapped_file} -f {feature_folder}{tags[n]}_types.txt {root}{genome_files[n2]} {root}{genome_files[n]}")
                
                    os.system(command)

                    if os.path.isdir(lifton_output_path):
                        os.system(f"rm -r {lifton_output_path}")