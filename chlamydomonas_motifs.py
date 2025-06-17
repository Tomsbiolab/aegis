#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creating all annotation objects for a single species.

@author: David Navarro
"""

# INPUT: change for other species
from config.chlamydomonas import genome_files, annotation_files, annotation_transfer_files, features_path, pickle_path, species

from datetime import datetime
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objs as go
import time
import os
import pandas as pd
from modules.tools import pickle_load, pickle_save
from modules.geneclasses import Genome, Annotation
from scipy.stats import false_discovery_control
from statistics import mean

jaspar_motifs = "../../chlamviz/motifs/jaspar_motifs.csv"
jaspar_motifs_other = "../../chlamviz/motifs/jaspar_motifs_other.csv"
early_up_file = "../../chlamviz/motifs/early_up.txt"
early_down_file = "../../chlamviz/motifs/early_down.txt"
control_10_file = "../../chlamviz/motifs/control_10.txt"
control_30_file = "../../chlamviz/motifs/control_30.txt"
out_stats_main = "../../genomes_and_annotation/algae/chlamydomonas_reinhardtii/6.1_annotation/tf_motifs/"

os.system(f"mkdir -p {out_stats_main}")

start = time.time()

genomes = {}

for k, v in genome_files.items():
    genomes[k] = Genome(k, v)

k = "Chlre_6.1_on_Chlre_6.0"
annot_name, genome = k.split("_on_")
location = f"{pickle_path}{k}_annotation.pkl"
if os.path.isfile(location):
    print(f"Loading {k} annotation object\n")
    a = pickle_load(location)
else:
    print(f"{k} annotation does not exist")
    a = Annotation(annot_name, annotation_files[k], genomes[genome])
    pickle_save(location, a)

a.add_gene_names_and_descriptors("../../genomes_and_annotation/algae/chlamydomonas_reinhardtii/6.1_annotation/names/names.csv", sep=",", just_gene_names=True)

early_up = []
f_in = open(early_up_file, encoding="utf-8")
for line in f_in:
    line = line.strip()
    early_up.append(line)
f_in.close()

early_down = []
f_in = open(early_down_file, encoding="utf-8")
for line in f_in:
    line = line.strip()
    early_down.append(line)
f_in.close()

control_10 = []
f_in = open(control_10_file, encoding="utf-8")
for line in f_in:
    line = line.strip()
    control_10.append(line)
f_in.close()

control_30 = []
f_in = open(control_30_file, encoding="utf-8")
for line in f_in:
    line = line.strip()
    control_30.append(line)
f_in.close()

motif_files = [jaspar_motifs, jaspar_motifs_other]
motif_tags = ["original_list", "other_list"]
promoter_sizes = [2000, 1000, 500]
promoter_tags = ["2kb", "1kb", "0.5kb"]
backgrounds = ["genomic", "control"]
query_list_files = [early_up, early_down, control_10, control_30]
query_list_tags = ["early_up", "early_down", "control_10", "control_30"]
control_background_index = 3

for h, motif_file in enumerate(motif_files):
    for y, promoter_size in enumerate(promoter_sizes):
        a.generate_promoters(genomes[genome], promoter_size, from_ATG=False, count_prom_size_from_TSS=True)
        for background in backgrounds:
            out_stats_path = f"{out_stats_main}{motif_tags[h]}_{promoter_tags[y]}/"
            mega_d = {}
            for n, l in enumerate(query_list_files):
                if background != "genomic" and n == control_background_index:
                    continue
                f_in = open(motif_file, encoding="utf-8")
                x = -1

                motifs_list = []
                odds_p_values = []

                z = -1
                for line in f_in:
                    x += 1
                    if x > 0:
                        z += 1
                        line = line.strip().split("\t")
                        tf = line[0]
                        motif = line[1]
                        motfif_length = int(line[2])
                        motifid = line[3]
                        if background == "genomic":
                            odds_p_values.append(a.find_motifs(l, motif, motfif_length, query_list_tags[n], motifid, tf, custom_path=out_stats_path))
                        else:
                            odds_p_values.append(a.find_motifs(l, motif, motfif_length, query_list_tags[n], motifid, tf, query_list_files[control_background_index], query_list_tags[control_background_index], custom_path=out_stats_path))
                        if motifid:
                            motifs_list.append(f"{tf}_{motifid}_{query_list_tags[n]}")
                            if f"{tf}_{motifid}" not in mega_d:
                                mega_d[f"{tf}_{motifid}"] = {query_list_tags[n]:(odds_p_values[z][8], odds_p_values[z][9], odds_p_values[z][10])}
                            else:
                                mega_d[f"{tf}_{motifid}"][query_list_tags[n]] = (odds_p_values[z][8], odds_p_values[z][9], odds_p_values[z][10])
                        else:
                            motifs_list.append(f"{tf}_{query_list_tags[n]}")
                            if f"{tf}" not in mega_d:
                                mega_d[f"{tf}"] = {query_list_tags[n]:(odds_p_values[z][8], odds_p_values[z][9], odds_p_values[z][10])}
                            else:
                                mega_d[f"{tf}"][query_list_tags[n]] = (odds_p_values[z][8], odds_p_values[z][9], odds_p_values[z][10])

                f_in.close()

                p_values_occurrences = []
                p_values_proportion = []
                for o in odds_p_values:
                    p_values_occurrences.append(o[2])
                    p_values_proportion.append(o[6])

                p_values_occurrences_adjusted = false_discovery_control(p_values_occurrences)
                p_values_occurrences_adjusted = list(p_values_occurrences_adjusted)
                p_values_proportion_adjusted = false_discovery_control(p_values_proportion)
                p_values_proportion_adjusted = list(p_values_proportion_adjusted)

                d = {}

                for j, m in enumerate(motifs_list):
                    if odds_p_values[j][8] != []:
                        mean_interest = mean(odds_p_values[j][8])
                    else:
                        mean_interest = 0

                    if odds_p_values[j][9] != []:
                        mean_background = mean(odds_p_values[j][9])
                    else:
                        mean_background = 0
                    d[m] = [odds_p_values[j][0], int(mean_interest), odds_p_values[j][1], int(mean_background), odds_p_values[j][4], odds_p_values[j][5], odds_p_values[j][6], p_values_proportion_adjusted[j], odds_p_values[j][7]]
                
                columns = ["occurrences in list", "avg motif num in list promoters with motif", "occurrences in background", "avg motif num in background promoters with motif", "percentage promoter in list", "percentage promoter in background distribution", "pvalue proportion", "pvalue proportion adjusted", "odds ratio proportion"]

                df = pd.DataFrame.from_dict(d, orient="index", columns=columns)

                if background == "genomic":
                    df.to_csv(f"{out_stats_path}{query_list_tags[n]}_summary_stats_genomic_background.csv", sep="\t", encoding="utf-8")
                else:
                    df.to_csv(f"{out_stats_path}{query_list_tags[n]}_summary_stats_control_background.csv", sep="\t", encoding="utf-8")

            if background == "genomic":
                for m in mega_d:
                    boxdata = []
                    x_axis = []
                    x_count = []
                    x = -1
                    for l in mega_d[m]:
                        x += 1
                        if x == 0:
                            x_count.append(x)
                            x_axis.append("background")
                            if mega_d[m][l][1] != []:
                                boxdata.append(mega_d[m][l][1])
                            else:
                                boxdata.append([0])
                            outfile = mega_d[m][l][2]
                        x_count.append(x+1)
                        x_axis.append(l)
                        if mega_d[m][l][0] != []:
                            boxdata.append((mega_d[m][l][0]))
                        else:
                            boxdata.append([0])


                    plt.clf()
                    # Create boxplot
                    sns.boxplot(data=boxdata)

                    # Adding labels
                    plt.xlabel("Selected Gene Lists")
                    plt.ylabel("Number of occurrences in promoters with motif")
                    plt.title(f"{m} motif frequency distribution")
                    plt.xticks(x_count, x_axis)

                    # Save the plot as a PDF file
                    #timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
                    plt.savefig(f"{outfile}_frequency_distribution.pdf")

                    # # Create layout
                    # layout = go.Layout(
                    #     title='Weight Boxplot',
                    #     xaxis=dict(title='Samples'),
                    #     yaxis=dict(title='Weight')
                    # )

                    # # Create figure
                    # fig = go.Figure(data=traces, layout=layout)

                    # # Save the plot as HTML file
                    # fig.write_html(f"{outfile}_frequency_distribution.html")

                    # plt.boxplot(boxdata, whiskerprops={'linewidth': 0})

                    # plt.xticks(x_count, x_axis)
                    # plt.grid(False)
                    # plt.savefig(f"{outfile}_{m}_frequency_distribution.png")


now = time.time()
lapse = now - start
print(f"\nCreating and updating {species} annotation objects as well as "
      f"finding motifs took {round((lapse/60)/60, 1)} hours\n")