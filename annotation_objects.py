#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Creating all annotation objects for a single species.

@author: David Navarro
"""

import argparse
import time
import os
import importlib
import gc

from config.paths import all_species
from modules.tools import pickle_load, pickle_save
from modules.annotation import Annotation
from modules.genome import Genome

def str_to_bool(value):
    """Convert a string to a boolean."""
    if value.lower() in ('true', '1', 'yes', "t"):
        return True
    elif value.lower() in ('false', '0', 'no', "f"):
        return False
    else:
        raise argparse.ArgumentTypeError(f"Invalid boolean value: '{value}'. Use 'True' or 'False'.")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--species', nargs='+', type=str, help='A list of config species tags (e.g. tomato mulberry). By default all species will be processed. Please specify here just the correct species if you are going to use the --annotations parameter', required=False, default=[])
parser.add_argument('--annotations', nargs='+', type=str, help='Use a list of tags if not all of the species annotaions are to be processed (e.g. v3_on_12Xv2 v4.3_on_40X_ref).', required=False, default=[])
parser.add_argument('--combine_transcripts', type=str_to_bool, help='Generate combined version of annotation files where alternative transcripts are merged into a "meta" transcript. This may be useful for exon-level read counting.', required=False, default=False)
parser.add_argument('--aegis', type=str_to_bool, help='Generate aegis versions of annotation files. This includes renaming ids.', required=False, default=True)
parser.add_argument('--clean', type=str_to_bool, help='Remove miscelaneous tags in annotation attributes.', required=False, default=True)
parser.add_argument('--export_features', type=str_to_bool, help='Export feature sequences in fasta format.', required=False, default=True)
parser.add_argument('--featurecounts', type=str_to_bool, help='Add featurecounts ID.', required=False, default=True)
parser.add_argument('--export_gffs', type=str_to_bool, help='Export gffs.', required=False, default=True)
parser.add_argument('--utrs', type=str_to_bool, help='Include UTRs in output gffs.', required=False, default=True)
parser.add_argument('--skip_atypical_features', type=str_to_bool, help='Exclude atypical features in output gffs.', required=False, default=True)
parser.add_argument('--export_lengths', type=str_to_bool, help='Export gene lengths.', required=False, default=True)
parser.add_argument('--export_coordinates', type=str_to_bool, help='Export gene coordinates.', required=False, default=True)
parser.add_argument('--export_genome_feature_sizes', type=str_to_bool, help='Export genome feature sizes.', required=False, default=True)
parser.add_argument('--export_stats', type=str_to_bool, help='Export annotation stats.', required=False, default=True)
parser.add_argument('--prepare_for_dap', type=str_to_bool, help='Prepare and export genomes and annotation files for dap.', required=False, default=True)
parser.add_argument('--promoters', type=str_to_bool, help='Whether promoter sequences are generated, --export_features must be True for export.', required=False, default=False)
parser.add_argument('--promoter_size', type=int, help='Promoter size.', required=False, default=5000)
parser.add_argument('--export_genomes', type=int, help='Export genomes.', required=False, default=True)
parser.add_argument('--id_separator', type=str, help='Separator to use for subfeatures when renaming ids e.g. geneid_t001 uses "_" separator.', required=False, default="_")
parser.add_argument('--remove_scaffolds', type=int, help='Remove scaffolds including chr00 or similar.', required=False, default=True)
parser.add_argument('--remove_organelles', type=int, help='Remove organelles including ChrC, ChrM or similar.', required=False, default=True)
parser.add_argument('--only_process_genomes', type=str_to_bool, help='Process only genomes.', required=False, default=False)
parser.add_argument('--rename_chromosomes_based_on_config', type=str_to_bool, help='Rename chromosomes based on config.', required=False, default=True)
parser.add_argument('--sort_processes', type=int, help='Number of processes for sort genes function.', required=False, default=2)
parser.add_argument('--most_specific_id_level_for_feature_fastas', type=str, choices=["gene", "transcript", "CDS", "protein", "promoter"], help='The most specific id which is used in exported fastas.', required=False, default="promoter")
parser.add_argument('--add_coordinates_to_feature_fastas', type=str_to_bool, help='Include coordinates and other details in fasta headers.', required=False, default=False)
parser.add_argument('--process_annotations', type=str_to_bool, help='Change to "False" if only pickles are sought for. Overrides many other parameters if set to "False".', required=False, default=True)


args = parser.parse_args()

species_to_process = args.species
aegis = args.aegis
combine = args.combine_transcripts
dap = args.prepare_for_dap
features = args.export_features
gffs = args.export_gffs
utrs = args.utrs
skip_atypical = args.skip_atypical_features
lengths = args.export_lengths
coordinates = args.export_coordinates
genome_sizes = args.export_genome_feature_sizes
featurecounts = args.featurecounts
clean = args.clean
objects = args.annotations
promoters = args.promoters
promoter_size = int(args.promoter_size)
export_genomes = args.export_genomes
export_stats = args.export_stats
remove_scaffolds = args.remove_scaffolds
remove_organelles = args.remove_organelles
only_genomes = args.only_process_genomes
rename_chromosomes_based_on_config = args.rename_chromosomes_based_on_config
sort_processes = args.sort_processes
feature_id_level = args.most_specific_id_level_for_feature_fastas
fasta_header_details = args.add_coordinates_to_feature_fastas
id_separator = args.id_separator
process_annotations = args.process_annotations

max_x = 1000

parameters=f"aegis={aegis}, combine_transcripts={combine}, prepare_for_dap={dap}, export_features={features}, export_gffs={gffs}, utrs={utrs}, skip_atypical_features={skip_atypical}, export_lengths={lengths}, export_coordinates={coordinates}, export_genome_feature_sizes={genome_sizes}, featurecounts={featurecounts}, clean={clean}, promoters={promoters}, promoter_size={promoter_size}, export_genomes={export_genomes}, remove_scaffolds={remove_scaffolds}, remove_organelles={remove_organelles}, only_genomes={only_genomes}, rename_chromosomes_based_on_config={rename_chromosomes_based_on_config}, sort_processes={sort_processes}, feature_id_level={feature_id_level}, fasta_header_details={fasta_header_details}, id_separator={id_separator}, process_annotations={process_annotations}"

def import_config_module(module_name):
    module = importlib.import_module(f"config.{module_name}")
    globals().update(vars(module))

def process_annotation(annot, genomes, dap, coordinates_path, gff_path, features_path, clean, featurecounts, lengths, coordinates, features, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, stats, max_x, rename_chromosomes_based_on_config:bool, feature_id_level, fasta_header_details, id_separator, conf_genomes:dict={}, dap_genomes:dict={}):
    
    print(f"Entering process_annotation for {annot.id}")

    a = annot.copy()

    temp_genome = genomes[a.genome]

    if a.genome in conf_genomes and rename_chromosomes_based_on_config:
        temp_conf_genome = conf_genomes[a.genome]
    else:
        temp_conf_genome = temp_genome

    if a.genome in dap_genomes and dap:
        temp_dap_genome = dap_genomes[a.genome]
    else:
        temp_dap_genome = temp_conf_genome

    if dap and not a.dapfit:
        a.rename_chromosomes(temp_dap_genome.equivalences, dap=True)
    elif rename_chromosomes_based_on_config:
        a.rename_chromosomes(temp_conf_genome.equivalences)

    if aegis:
        if a.name == "4.0_4.2_merge":
            a.rename_ids(remove_point_suffix=True, sep=id_separator)
        else:
            a.rename_ids(sep=id_separator)

    a.update_attributes(clean=clean, featurecountsID=featurecounts)

    if features:
        if promoters:
            a.generate_promoters(temp_dap_genome, promoter_size)

        a.generate_sequences(temp_dap_genome)

        a.export_all_features("both", promoters, fasta_header_details, features_path, feature_id_level)

    if lengths:
        a.export_lengths(custom_path=coordinates_path)
        a.export_lengths(custom_path=coordinates_path, just_genes=False)
    
    if coordinates:
        a.export_coordinates(custom_path=coordinates_path)

    if stats:
        a.update_stats(custom_path=stats_path, export=stats, genome=temp_dap_genome, max_x=max_x)

    a.clear_sequences()

    if gffs:
        a.export_gff(custom_path=gff_path, skip_atypical_fts=skip_atypical, UTRs=utrs)
        if dap:
            a.export_gff(custom_path=gff_path, just_genes=True)

    if combine and gffs:
        b = a.copy()
        b.combine_transcripts(temp_dap_genome, respect_non_coding=True)
        b.update_attributes(featurecountsID=featurecounts)
        b.export_gff(custom_path=gff_path, skip_atypical_fts=skip_atypical, UTRs=utrs)
        if dap:
            b.export_gff(custom_path=gff_path, just_genes=True)
    
    a.CDS_to_CDS_segment_ids()

    if gffs:
        a.export_gff(custom_path=gff_path, skip_atypical_fts=skip_atypical, UTRs=utrs, no_1bp_features=True)

    print("Ending process_annotation")


if species_to_process == []:
    species_to_process = all_species

if objects == []:
    print(f"Working on the following species: {species_to_process}")
else:
    print(f"Working on the following species: {species_to_process}")
    print(f"Working on the following annotations: {objects}")

print(f"The chosen parameters are the following: {parameters}")

for config_species in species_to_process:

    print(f"Starting to work with {config_species}.")

    gc.collect()
    import_config_module(config_species)

    start = time.time()

    aegis_output = f"{path}aegis_output/"
    pickle_path = f"{aegis_output}pickles/"
    features_path = f"{aegis_output}features/"
    gff_path = f"{aegis_output}gffs/"
    genome_path = f"{aegis_output}genomes/"
    coordinates_path = f"{aegis_output}coordinates/"
    stats_path = f"{aegis_output}stats/"

    os.system(f"mkdir -p {features_path}")
    os.system(f"mkdir -p {pickle_path}")
    os.system(f"mkdir -p {gff_path}")
    os.system(f"mkdir -p {genome_path}")
    os.system(f"mkdir -p {coordinates_path}")
    os.system(f"mkdir -p {stats_path}")

    genomes = {}
    organelleless_genomes = {}
    scaffoldless_genomes = {}
    scaffoldless_organelleless_genomes = {}

    conf_genomes = {}
    conf_organelleless_genomes = {}
    conf_scaffoldless_genomes = {}
    conf_scaffoldless_organelleless_genomes = {}

    dap_genomes = {}
    dap_organelleless_genomes = {}
    dap_scaffoldless_genomes = {}
    dap_scaffoldless_organelleless_genomes = {}

    used_genomes = []

    for k in objects:
        g = k.split("_on_")[-1]
        used_genomes.append(g)

    for k, v in genome_files.items():

        if used_genomes != []:
            if k not in used_genomes:
                continue

        genomes[k] = Genome(k, v, chromosome_dict=genome_dictionary[k], rename_chromosomes=False)
        if export_genomes:
            genomes[k].export(output_folder=genome_path)
        if genome_sizes:
            genomes[k].export_feature_sizes(custom_path=coordinates_path)
        if remove_scaffolds:
            if genomes[k].non_chromosomal_scaffolds:
                scaffoldless_genomes[k] = genomes[k].copy()
                scaffoldless_genomes[k].remove_scaffolds(output_folder=genome_path, export=export_genomes, chromosome_dict=genome_dictionary[k])
        if remove_organelles:
            if genomes[k].organelles:
                organelleless_genomes[k] = genomes[k].copy()
                organelleless_genomes[k].remove_organelles(export=export_genomes, output_folder=genome_path)
        if remove_scaffolds and remove_organelles:
            if genomes[k].non_chromosomal_scaffolds and genomes[k].organelles:
                scaffoldless_organelleless_genomes[k] = genomes[k].copy()
                scaffoldless_organelleless_genomes[k].remove_scaffolds(output_folder=genome_path, export=export_genomes, chromosome_dict=genome_dictionary[k], remove_organelles=remove_organelles)

        if rename_chromosomes_based_on_config and genomes[k].chromosome_dict != {}:
            conf_genomes[k] = Genome(k, v, chromosome_dict=genome_dictionary[k], rename_chromosomes=rename_chromosomes_based_on_config)
            if export_genomes:
                conf_genomes[k].export(output_folder=genome_path)
            if genome_sizes:
                conf_genomes[k].export_feature_sizes(custom_path=coordinates_path)
            if remove_scaffolds:
                if conf_genomes[k].non_chromosomal_scaffolds:
                    conf_scaffoldless_genomes[k] = conf_genomes[k].copy()
                    conf_scaffoldless_genomes[k].remove_scaffolds(output_folder=genome_path, export=export_genomes, chromosome_dict=genome_dictionary[k])
            if remove_organelles:
                if conf_genomes[k].organelles:
                    conf_organelleless_genomes[k] = conf_genomes[k].copy()
                    conf_organelleless_genomes[k].remove_organelles(export=export_genomes, output_folder=genome_path)
            if remove_scaffolds and remove_organelles:
                if conf_genomes[k].non_chromosomal_scaffolds and conf_genomes[k].organelles:
                    conf_scaffoldless_organelleless_genomes[k] = conf_genomes[k].copy()
                    conf_scaffoldless_organelleless_genomes[k].remove_scaffolds(output_folder=genome_path, export=export_genomes, chromosome_dict=genome_dictionary[k], remove_organelles=remove_organelles)

        if dap and not genomes[k].dapfit:
            dap_genomes[k] = Genome(k, v, chromosome_dict=genome_dictionary[k], rename_chromosomes=rename_chromosomes_based_on_config)
            if not dap_genomes[k].dapfit:
                dap_genomes[k].rename_features_dap(output_folder=genome_path, return_equivalences=False, export=export_genomes, chromosome_dict=genome_dictionary[k])

                if genome_sizes:
                    dap_genomes[k].export_feature_sizes(custom_path=coordinates_path)
                if remove_scaffolds:
                    if dap_genomes[k].non_chromosomal_scaffolds:
                        dap_scaffoldless_genomes[k] = dap_genomes[k].copy()
                        dap_scaffoldless_genomes[k].remove_scaffolds(output_folder=genome_path, export=export_genomes, chromosome_dict=genome_dictionary[k])
                if remove_organelles:
                    if dap_genomes[k].organelles:
                        dap_organelleless_genomes[k] = dap_genomes[k].copy()
                        dap_organelleless_genomes[k].remove_organelles(export=export_genomes, output_folder=genome_path)
                if remove_scaffolds and remove_organelles:
                    if dap_genomes[k].non_chromosomal_scaffolds and dap_genomes[k].organelles:
                        dap_scaffoldless_organelleless_genomes[k] = dap_genomes[k].copy()
                        dap_scaffoldless_organelleless_genomes[k].remove_scaffolds(output_folder=genome_path, export=export_genomes, chromosome_dict=genome_dictionary[k], remove_organelles=remove_organelles)

    if not only_genomes:

        for o in objects:
            if o not in annotation_files and o not in annotation_transfer_files:
                print(f"Warning: Annotation object {o} appears to be missing in config files.")

        annotations = {}

        for k, v in annotation_files.items():
            location = f"{pickle_path}{k}_annotation.pkl"
            if os.path.isfile(location):

                if objects != []:
                    if k in objects:
                        print(f"\nLoading {k} annotation object\n")
                        annotations[k] = pickle_load(location)
                else:
                    print(f"\nLoading {k} annotation object\n")
                    annotations[k] = pickle_load(location)
            else:
                annot_name, genome = k.split("_on_")
                if objects != []:
                    if k in objects:
                        if os.path.isfile(v):
                            annotations[k] = Annotation(annot_name, v, genomes[genome], define_synteny=True, sort_processes=sort_processes)
                            pickle_save(location, annotations[k])
                        else:
                            print(f"Warning: {v} file does not exist.")

                else:
                    if os.path.isfile(v):
                        annotations[k] = Annotation(annot_name, v, genomes[genome], define_synteny=True, sort_processes=sort_processes)
                        pickle_save(location, annotations[k])
                    else:
                        print(f"Warning: {v} file does not exist.")

        if process_annotations:

            for k, a in annotations.items():
                if objects != []:
                    if k in objects:
                        if rename_chromosomes_based_on_config:
                            process_annotation(a, genomes, False, coordinates_path, gff_path, features_path, clean, featurecounts, lengths, coordinates, features, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, export_stats, max_x, False, feature_id_level, fasta_header_details, id_separator)
                        process_annotation(a, genomes, False, coordinates_path, gff_path, features_path, clean, featurecounts, lengths, coordinates, features, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, export_stats, max_x, rename_chromosomes_based_on_config, feature_id_level, fasta_header_details, id_separator, conf_genomes)
                        if dap and not a.dapfit:
                            process_annotation(a, genomes, dap, coordinates_path, gff_path, features_path, clean, featurecounts, False, coordinates, False, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, export_stats, max_x, rename_chromosomes_based_on_config, feature_id_level, fasta_header_details, id_separator, conf_genomes, dap_genomes)

                else:
                    if rename_chromosomes_based_on_config:
                        process_annotation(a, genomes, False, coordinates_path, gff_path, features_path, clean, featurecounts, lengths, coordinates, features, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, export_stats, max_x, False, feature_id_level, fasta_header_details, id_separator)
                    process_annotation(a, genomes, False, coordinates_path, gff_path, features_path, clean, featurecounts, lengths, coordinates, features, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, export_stats, max_x, rename_chromosomes_based_on_config, feature_id_level, fasta_header_details, id_separator, conf_genomes)
                    if dap and not a.dapfit:
                        process_annotation(a, genomes, dap, coordinates_path, gff_path, features_path, clean, featurecounts, False, coordinates, False, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, export_stats, max_x, rename_chromosomes_based_on_config, feature_id_level, fasta_header_details, id_separator, conf_genomes, dap_genomes)


        transferred_annotations = {}

        for k, v in annotation_transfer_files.items():
            annot_name, origin = k.split("_on_")[0].split("_from_")
            origin_annot_name = annot_name
            origin_annot_name = origin_annot_name.replace("Lifton", "")
            genome = k.split("_on_")[1]
            tag = f"{annot_name}_from_{origin}_on_{genome}"
            location = f"{pickle_path}{tag}_annotation.pkl"
            if os.path.isfile(location):
                if objects != []:
                    if tag in objects:
                        print(f"\nLoading {tag} annotation object\n")
                        transferred_annotations[tag] = pickle_load(location)
                else:
                    print(f"\nLoading {tag} annotation object\n")
                    transferred_annotations[tag] = pickle_load(location)
            else:
                if objects != []:
                    if tag in objects:
                        if os.path.isfile(v):
                            transferred_annotations[tag] = Annotation(f"{annot_name}_from_{origin}", v, genomes[genome], annotations[f"{origin_annot_name}_on_{origin}"], define_synteny=True, sort_processes=sort_processes)
                            pickle_save(location, transferred_annotations[tag])       
                        else:
                            print(f"Warning: {v} file does not exist.")
                else:
                    if os.path.isfile(v):
                        transferred_annotations[tag] = Annotation(f"{annot_name}_from_{origin}", v, genomes[genome], annotations[f"{origin_annot_name}_on_{origin}"], define_synteny=True, sort_processes=sort_processes)
                        pickle_save(location, transferred_annotations[tag])
                    else:
                        print(f"Warning: {v} file does not exist.")

        if process_annotations:

            for k, a in transferred_annotations.items():
                if objects != []:
                    if k in objects:
                        if rename_chromosomes_based_on_config:
                            process_annotation(a, genomes, False, coordinates_path, gff_path, features_path, clean, featurecounts, lengths, coordinates, features, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, export_stats, max_x, False, feature_id_level, fasta_header_details, id_separator)
                        process_annotation(a, genomes, False, coordinates_path, gff_path, features_path, clean, featurecounts, lengths, coordinates, features, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, export_stats, max_x, rename_chromosomes_based_on_config, feature_id_level, fasta_header_details, id_separator, conf_genomes)
                        if dap and not a.dapfit:
                            process_annotation(a, genomes, dap, coordinates_path, gff_path, features_path, clean, featurecounts, False, coordinates, False, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, export_stats, max_x, rename_chromosomes_based_on_config, feature_id_level, fasta_header_details, id_separator, conf_genomes, dap_genomes)

                else:
                    if rename_chromosomes_based_on_config:
                        process_annotation(a, genomes, False, coordinates_path, gff_path, features_path, clean, featurecounts, lengths, coordinates, features, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, export_stats, max_x, False, feature_id_level, fasta_header_details, id_separator)
                    process_annotation(a, genomes, False, coordinates_path, gff_path, features_path, clean, featurecounts, lengths, coordinates, features, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, export_stats, max_x, rename_chromosomes_based_on_config, feature_id_level, fasta_header_details, id_separator, conf_genomes)
                    if dap and not a.dapfit:
                        process_annotation(a, genomes, dap, coordinates_path, gff_path, features_path, clean, featurecounts, False, coordinates, False, promoters, promoter_size, aegis, gffs, skip_atypical, utrs, combine, stats_path, export_stats, max_x, rename_chromosomes_based_on_config, feature_id_level, fasta_header_details, id_separator, conf_genomes, dap_genomes)

    else:
        if objects != []:
            print(f"\nWarning: only genomes mode was selected so annotations objects selected: '{objects}' will not be processed.")
        else:
            print(f"\nNote: only genomes mode was selected so annotations selected so no annotation will be processed")

    now = time.time()
    lapse = now - start

    print(f"\nCreating and updating {species} annotation objects with annotations='{objects}' took {round((lapse/60)/60, 1)} hours\n")
    
    gc.collect()