#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on May 19 2025

Reading annotation objects, calculating and exporting overlaps.

@author: David Navarro
"""

import time
import os
from modules.tools import pickle_load
from modules.annotation import Annotation

from config.cannabis import annotation_files, annotation_transfer_files, root
from config.PN40024 import annotation_files as annotation_files_2
from config.PN40024 import annotation_transfer_files as annotation_transfer_files_2

start = time.time()

aegis_output = f"{root}aegis_output/"
pickle_paths = [f"{root}cannabis/aegis_output/pickles/", f"{root}grapevine/PN40024/aegis_output/pickles/"]
overlaps_path = f"{aegis_output}overlaps/"

simple_annotation_tags = ["CBDRx", "PurpleKush", "JamaicanLion", "Finola", "5.1"]
genome_tags = ["CBDRx", "PurpleKush", "JamaicanLion", "Finola", "T2T_ref"]

if len(simple_annotation_tags) != len(genome_tags):
    raise ValueError(f"Input annotation or genome tags do not match.")

indispensable_members = ["5.1"]

os.system(f"mkdir -p {overlaps_path}")

lifton = True

annotation_files.update(annotation_transfer_files)
annotation_files_2.update(annotation_transfer_files_2)

annotations = {}

# Completing TAGS:
lifton_annotation_tags = []
liftoff_annotation_tags = []
annotation_tags = []

for n, tag in enumerate(simple_annotation_tags):
    for n2, tag2 in enumerate(simple_annotation_tags):
        if n == n2:
            continue

        lifton_annotation_tags.append(f"{tag}Lifton_from_{genome_tags[n]}_on_{genome_tags[n2]}")

for n, tag in enumerate(simple_annotation_tags):
    for n2, tag2 in enumerate(simple_annotation_tags):
        if n == n2:
            continue

        liftoff_annotation_tags.append(f"{tag}_from_{genome_tags[n]}_on_{genome_tags[n2]}")

for n, tag in enumerate(simple_annotation_tags):
    annotation_tags.append(f"{tag}_on_{genome_tags[n]}")

if lifton:
    program = "Lifton"
else:
    program = "Liftoff"

for tag in annotation_tags:
    if tag not in annotation_files and tag not in annotation_files_2:
        print(f"Warning: tag {tag} chosen is missing in config file.")

for tag in lifton_annotation_tags:
    if tag not in annotation_files and tag not in annotation_files_2:
        print(f"Warning: tag {tag} chosen is missing in config file.")

for tag in liftoff_annotation_tags:
    if tag not in annotation_files and tag not in annotation_files_2:
        print(f"Warning: tag {tag} chosen is missing in config file.")

tags_to_process = annotation_tags.copy()

if lifton:
    tags_to_process += lifton_annotation_tags
else:
    tags_to_process += liftoff_annotation_tags

for n, tag in enumerate(tags_to_process):
    location1 = f"{pickle_paths[0]}{tag}_annotation.pkl"
    location2 = f"{pickle_paths[1]}{tag}_annotation.pkl"

    if os.path.isfile(location1):
        print(f"Loading {tag} annotation object\n")
        annotations[tag] = pickle_load(location1)
    elif os.path.isfile(location2):
        print(f"Loading {tag} annotation object\n")
        annotations[tag] = pickle_load(location2)
    else:
        print(f"Warning: file {location1} or {location2} does not exist yet.")

for n, simple_tag in enumerate(simple_annotation_tags):

    for n2, simple_tag2 in enumerate(simple_annotation_tags):

        if simple_tag == simple_tag2:
            continue

        if indispensable_members != []:

            # this excludes certain overlap combinations
            if simple_tag not in indispensable_members and simple_tag2 not in indispensable_members:
                continue

        if lifton:
            query_tag = f"{simple_tag2}Lifton_from_{genome_tags[n2]}_on_{genome_tags[n]}"

            target_tag = f"{simple_tag}_on_{genome_tags[n]}"
        else:

            query_tag = f"{simple_tag2}_from_{genome_tags[n2]}_on_{genome_tags[n]}"

            target_tag = f"{simple_tag}_on_{genome_tags[n]}"

        print(f"Detecting and exporting overlaps for query='{query_tag}' vs target='{target_tag}'")

        annotations[query_tag].detect_gene_overlaps(annotations[target_tag])
        annotations[query_tag].export_equivalences(custom_path=overlaps_path, stringent=False, return_df=False, NAs=False, export_csv=True, output_file=f"{simple_tag2}_equivalences_{program}_on_{simple_tag}.tsv")

now = time.time()
lapse = now - start
print(f"\nGenerating and exporting all pairwise cannabis annotation overlaps took {round((lapse/60)/60, 1)} hours\n")