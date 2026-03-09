"""
Created on Fri Oct 7 15:09:46 2022

Module defining several genomic classes.

@authors: David Navarro, Antonio Santiago
"""

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .genome import Genome

import sys
import copy
import random
import time
import pandas as pd
import os
import warnings
import math
import gc

from tqdm import tqdm
from multiprocessing import Pool
from pathlib import Path

from .feature import Feature
from .gene import Gene
from .transcript import Transcript
from .subfeatures import Exon, UTR
from .hits import BlastHit
from .utils.genefunctions import sort_and_update_genes
from .utils.misc import read_file_with_fallback
from .utils.gtf_gff import parse_gff_line, convert_gtf_to_gff3, detect_file_format
from .annotation_components.stats import AnnotationStats
from .annotation_components.export import AnnotationExport
from .annotation_components.motifs import AnnotationMotifs
from .annotation_components.overlaps import AnnotationOverlaps
from .annotation_components.redundancy import AnnotationRedundancy
from .conf import default_noncoding_transcripts, default_features_r


class Annotation():

    stats: AnnotationStats
    export: AnnotationExport
    motifs: AnnotationMotifs
    genome: str | None
    all_gene_ids:dict[str, str]
    # the tuple values consist of transcript id keys consist of (chromosome, gene_id)
    all_transcript_ids:dict[str, tuple[str, str]]
    _gene_info:dict[str, set[str]]
    _transcript_info:dict[str, set[str]]

    chrs:dict[str, dict[str, Gene]]
    
    bar_colors = ["31", "32", "33", "33", "33", "34"]

    def __init__(self, annot_file_path:str, name:str|None=None, genome:Genome|None=None, original_annotation:Annotation|None=None, target:bool=False, to_overlap:bool=True, rework_all_CDSs:bool=False, work_out_missing_CDSs:bool=False, chosen_chromosomes:tuple[str, ...]|None=None, chosen_coordinates:tuple[int, int]|None=None, sort_processes:int=1, define_synteny=False, rename_features:list=[], keep_existing_ids_if_derived_from_base_id:bool=False, quiet:bool=False, consider_polycistronic:bool=False, consider_read_utrs:bool=False, infer_genes_from_transcripts:bool=True, infer_genes_from_subfeatures:bool=True, skip_orphaned_features:bool=True, skip_atypical_features:bool=True, incorporate_and_rename_repeated_ids:bool=True, collapse_exons:bool=True, collapse_CDSs:bool=True, standardise_features:bool=False):
        
        start_time = time.time()

        pid = os.getpid()

        self.file = str(Path(annot_file_path).resolve())
        self.path = str(Path(annot_file_path).resolve().parent) + "/" 

        if name is None:
            self.name = Path(annot_file_path).stem
        else:
            self.name = name

        if chosen_chromosomes is not None and len(chosen_chromosomes) > 0:
            if len(chosen_chromosomes) > 1:
                self.name = f"{self.name}_{chosen_chromosomes[0]}-{chosen_chromosomes[-1]}"
            else:
                self.name = f"{self.name}_{chosen_chromosomes[0]}"

        self.gff_header = []
        self.target = target
        self.to_overlap = to_overlap
        self.liftover = False
        self.merged = False
        self.sorted = False
        
        if genome is not None:
            self.genome = genome.name
            self.id = f"{self.name}_on_{self.genome}"
        else:
            self.genome = None
            self.id = self.name

        if not quiet:
            print(f"\nProcessing {self.id} annotation object\n")
        self.excluded_chromosomes = []
        self.features = {}
        # genes will be added as {"ch":{"gene_id" : gene_object}}
        self.chrs = {}
        # we save here chr as the value corresponding to each gene id
        self.all_gene_ids = {}
        # we save here chr, gene id as the tuple corresponding to each transcript id
        self.all_transcript_ids = {}
        self.all_protein_ids = {}
        self.protein_equivalences = {}
        self.unmapped = []
        # Here we insert any feature which is of an unknown category
        self.atypical_features = []
        # Here we insert any recognisable feature that was impossible to fit into the current structure
        self.orphaned_features = []

        self.feature_suffix = ""
        self.suffix = ""

        self.coding_removed = False
        if "_minus_coding" in annot_file_path:
            self.coding_removed = True

        self.non_coding_removed = False
        if "_minus_non_coding" in annot_file_path:
            self.non_coding_removed = True

        self.transposable_removed = False
        if "_minus_TE" in annot_file_path:
            self.transposable_removed = True

        self.non_transposable_removed = False
        if "_minus_non_TE" in annot_file_path:
            self.non_transposable_removed = True

        self.featurecounts = False
        if "_fcounts" in annot_file_path:
            self.featurecounts = True

        self.combined = False
        if "_combined" in annot_file_path:
            self.combined = True

        self.clean = False
        if "_clean" in annot_file_path:
            self.clean = True

        self.dapmod = False
        if "_dapmod" in annot_file_path:
            self.dapmod = True

        self.confrenamed = False
        if "_confrenamed" in annot_file_path:
            self.confrenamed = True

        self.full_renamed_ids = False
        if "_full_renamed_ids" in annot_file_path:
            self.confrenamed = True

        if genome is not None:
            self.dapfit = genome.dapfit
        else:
            self.dapfit = None

        self.symbols_added = False

        self.renamed_features = []

        self.small_cds_removed = False
        if "_minus_small_CDSs" in annot_file_path:
            self.small_cds_removed = True

        self.promoter_size = 3000
        
        self.generated_all_sequences = False
        self.contains_all_sequences = False        
        self.generated_CDS_sequences = False
        self.contains_CDS_sequences = False
        self.generated_protein_sequences = False
        self.contains_protein_sequences = False

        self.promoter_types = "standard"
        if original_annotation != None:
            self.liftover = True

        keys = [
            "repeated_gene_IDs", "1bp_gene", "1bp_transcript", "subfeature_to_gene",
            "transcript_to_more_than_1_gene", "transcript_to_gene_other_chr", "subfeature_to_transcript_other_chr", "subfeature_to_gene_other_chr"
        ]
        # Create a dictionary where every value is an empty set
        self.errors = {key: set() for key in keys}
        
        keys = [
            "1bp_exon", "1bp_CDS", "1bp_UTR", "1bp_other_subfeature", "decreasing_coordinates",
            "missing_subfeature_parent", "transcript_to_inexistent_gene", "transcript_with_no_parent",
            "missing_subfeature_parent_liftover", "multiple_CDSs_per_transcript",
            "possible_policistronic_transcript", "transcript_with_no_exons",
            "gene_with_no_transcripts", "subfeature_with_no_parent", "subfeature_to_gene", "repeat_transcript_different_genes", "repeat_transcript_same_gene"
        ]

        self.warnings = {key: set() for key in keys}

        del keys
        
        encoding = read_file_with_fallback(self.file)
        file_format = detect_file_format(self.file, encoding)
        if not quiet:
            print(f"{file_format} file format detected for file='{self.file}'")

        if file_format == 'gtf':
            gff_file = f"{self.name}_{pid}.tmp"
            convert_gtf_to_gff3(self.file, gff_file, encoding, quiet=quiet)
        else:
            gff_file = self.file

        staging = self.load_data(gff_file, encoding=encoding, chosen_chromosomes=chosen_chromosomes, chosen_coordinates=chosen_coordinates, skip_atypical_features=skip_atypical_features, quiet=quiet)

        self.stats = AnnotationStats(self)
        self.export = AnnotationExport(self)
        self.motifs = AnnotationMotifs(self)
        self.overlaps = AnnotationOverlaps(self)
        self.redundancy = AnnotationRedundancy(self)

        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected or quite mode
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False

        self._gene_info = {}
        self._transcript_info = {}

        batch_size = 5000

        for n, stage in enumerate(staging):

            progress_bar = tqdm(total=len(staging[stage]), disable=disable, bar_format=(
                f'\033[1;{Annotation.bar_colors[n]}mAdding {stage}s:\033[0m '
                '{percentage:3.0f}%|'
                f'\033[1;{Annotation.bar_colors[n]}m{{bar}}\033[0m| '
                '{n}/{total} [{elapsed}<{remaining}]'))
            
            count = 0

            for entry in staging[stage]:
                if stage in ["UTR", "CDS", "exon", "other_subfeature"]:
                    self._add_subfeature(entry, ft_level=stage, infer_gene_and_transcript_from_subfeatures=infer_genes_from_subfeatures, skip_orphaned_features=skip_orphaned_features, quiet=quiet)
                elif stage == "transcript":
                    self._add_transcript(entry, rename_repeated_id=incorporate_and_rename_repeated_ids, infer_gene_from_transcript=infer_genes_from_transcripts, skip_orphaned_features=skip_orphaned_features, quiet=quiet)
                else:
                    self._add_gene(entry, rename_repeated_id=incorporate_and_rename_repeated_ids, quiet=quiet)

                count += 1
                if count >= batch_size:
                    progress_bar.update(count)
                    count = 0

            if count > 0:
                progress_bar.update(count)

            progress_bar.close()
            staging[stage] = []
            gc.collect()

        del self._gene_info
        del self._transcript_info
        del staging

        gc.collect()

        if self.file != gff_file:
            os.remove(gff_file)

        now = time.time()
        lapse = now - start_time

        if not quiet:
            print(f"\nCreating {self.id} annotation object took {round(lapse/60, 1)} minutes\n")

        self.update(original_annotation=original_annotation, genome=genome, sort_processes=sort_processes, define_synteny=define_synteny, rename_features=rename_features, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id, quiet=quiet, consider_polycistronic=consider_polycistronic, consider_read_utrs=consider_read_utrs, collapse_exons=collapse_exons, collapse_CDSs=collapse_CDSs, standardise_features=standardise_features)

        if (rework_all_CDSs or work_out_missing_CDSs) and genome:
            self.rework_CDSs(genome, override=rework_all_CDSs, quiet=quiet)
            self.update(original_annotation=original_annotation, genome=genome, sort_processes=sort_processes, define_synteny=define_synteny, rename_features=rename_features, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id, quiet=quiet, consider_polycistronic=consider_polycistronic, consider_read_utrs=consider_read_utrs, collapse_exons=collapse_exons, collapse_CDSs=collapse_CDSs, standardise_features=standardise_features)

    def load_data(self, gff_file, encoding, chosen_chromosomes:tuple[str, ...]|None=None, chosen_coordinates:tuple[int, int]|None=None, skip_atypical_features:bool=False, quiet:bool=False):
        
        staging = {"gene": [], "transcript": [], "exon": [], "CDS": [], "UTR": [], "other_subfeature": []}

        chromosomes_t = set()

        with open(gff_file, encoding=encoding) as f:

            for line in f:
                line = line.strip()

                if not line:
                    continue

                if line[0] == "#":
                    if line == "###" or line == "#":
                        continue
                    self.gff_header.append(line)
                    continue

                parts = line.split("\t")
                if len(parts) <= 2:
                    self.gff_header.append(line)
                    continue

                if len(parts) < 9:
                    continue

                entry = parse_gff_line(line)

                ID = entry["id"]
                ft = entry["feature"]
                ch = entry["ch"]
                start = entry["start"]
                end = entry["end"]

                if chosen_chromosomes is not None:
                    if ch not in chosen_chromosomes:
                        continue

                if chosen_coordinates is not None:

                    if start < chosen_coordinates[0]:
                        continue
                    if end > chosen_coordinates[1]:
                        continue

                cat = default_features_r.get(ft, "atypical")

                if not ID:
                    if not quiet:
                        print(f"{self.id} Warning: ID not found for {ft}: {line}")

                source = entry["source"]
                strand = entry["strand"]
                score = entry["score"]
                parents = entry["parents"]
                
                attributes = entry["attributes"]
                decreasing_coordinates = entry["decreasing_coordinates"]

                if decreasing_coordinates:
                    self.warnings["decreasing_coordinates"].add(ID)
                    if not quiet:
                        print(f"{self.id} Warning: Decreasing coordinates for {ft} {ID}")

                if cat == "atypical":
                    if not skip_atypical_features:
                        if ft not in self.features:
                            self.features[ft] = 1
                        else:
                            self.features[ft] += 1
                        chromosomes_t.add(ch)
                        self.atypical_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))
                else:
                    if ft not in self.features:
                        self.features[ft] = 1
                    else:
                        self.features[ft] += 1
                    chromosomes_t.add(ch)
                    staging[cat].append(entry)

        chromosomes_t = sorted(list(chromosomes_t))
        for ch in chromosomes_t:
            self.chrs[ch] = {}

        return staging

    def _add_gene(self, entry, rename_repeated_id:bool=False, quiet:bool=False):
        ID = entry["id"]
        transposable = entry["transposable"]
        pseudogene = entry["pseudogene"]
        ch = entry["ch"]
        source = entry["source"]
        ft = entry["feature"]
        strand = entry["strand"]
        start = entry["start"]
        end = entry["end"]
        score = entry["score"]
        attributes = entry["attributes"]

        if start == end:
            self.errors["1bp_gene"].add(ID)
            if not quiet:
                print(f"{self.id} Error: 1bp gene {ID} feature")

        if ID not in self.all_gene_ids:
            self.all_gene_ids[ID] = ch
            self.chrs[ch][ID] = Gene(pseudogene, transposable, ID, ch, source, ft, strand, start, end, score, ".", [], attributes)
        else:
            self.errors["repeated_gene_IDs"].add(ID)
            if not rename_repeated_id:
                if not quiet:
                    print(f"{self.id} Error: repeated gene ID {ID} will not be added based on the chosen parameter (rename_repeated_id={rename_repeated_id})")
            else:
                renamed_id = self._get_unique_gene_id(ID)
                if not quiet:
                    print(f"{self.id} Warning: repeated gene id {ID} was renamed to {renamed_id}")
                self.all_gene_ids[renamed_id] = ch
                self.chrs[ch][renamed_id] = Gene(pseudogene, transposable, renamed_id, ch, source, ft, strand, start, end, score, ".", [], attributes)
                if ID not in self._gene_info:
                    self._gene_info[ID] = {ID, renamed_id}
                else:
                    self._gene_info[ID].add(renamed_id)

    def _add_transcript(self, entry, rename_repeated_id:bool=False, infer_gene_from_transcript:bool=False, skip_orphaned_features:bool=False, quiet:bool=False):
        ID = entry["id"]
        ch = entry["ch"]
        source = entry["source"]
        ft = entry["feature"]
        strand = entry["strand"]
        start = entry["start"]
        end = entry["end"]
        score = entry["score"]
        attributes = entry["attributes"]
        parents = entry["parents"]

        if rename_repeated_id:
            unique_t_id = self._get_unique_transcript_id(ID)
        else:
            unique_t_id = ID

        if start == end:
            if not quiet:
                print(f"{self.id} Error: 1bp transcript level {ID} feature")
            self.errors["1bp_transcript"].add(ID)

        if len(parents) > 1:
            """
            So far, this does not happen, so no action taken yet.
            """
            if not quiet:
                print(f"{self.id} Error: {ID} transcript refers to more than 1 gene. Parents = {parents}. The transcript will not be added to any of the parental genes")
            self.errors["transcript_to_more_than_1_gene"].add(ID)
            if not skip_orphaned_features:
                self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))

        elif not parents:
            """
            The current assumption is that there will not be more than one parentless transcript actually referring to the same gene.
            """
            self.warnings["transcript_with_no_parent"].add(ID)
            if infer_gene_from_transcript:

                created_gene_id = f"{unique_t_id}_gene"

                if rename_repeated_id:
                    count = 0
                    while created_gene_id in self.all_gene_ids:
                        count += 1
                        unique_t_id = f"{unique_t_id}_{count}"
                        while unique_t_id in self.all_transcript_ids:
                            count += 1
                            unique_t_id = f"{unique_t_id}_{count}"
                        created_gene_id = f"{unique_t_id}_gene"

                if created_gene_id not in self.all_gene_ids and unique_t_id not in self.all_transcript_ids:
                    entry["id"] = created_gene_id
                    entry["feature"] = "gene"
                    self._add_gene(entry, quiet=quiet)

                    self.all_transcript_ids[unique_t_id] = (ch, created_gene_id)
                    self.chrs[ch][created_gene_id].transcripts[unique_t_id] = Transcript(unique_t_id, ch, source, ft, strand, start, end, score, ".", [created_gene_id], attributes)
                    if unique_t_id != ID:
                        if ID in self._transcript_info:
                            self._transcript_info[ID].add(unique_t_id)
                        else:
                            self._transcript_info[ID] = {ID, unique_t_id}
                        if not quiet:
                            print(f"{self.id} Warning: No parent provided so the following {created_gene_id} gene was created for {ID} transcript which was renamed to {unique_t_id} to avoid repeat ids")
                    elif not quiet:
                        print(f"{self.id} Warning: No parent provided so the following {created_gene_id} gene was created for {ID} transcript")
                            
                elif not quiet:
                    print(f"{self.id} Warning: {ID} transcript could not be assigned to a gene with the chosen parameter (rename_repeated_id={rename_repeated_id}) as no parent was provided and transcript id was not unique")
            elif not quiet:
                print(f"{self.id} Warning: {ID} transcript could not be assigned to a gene with the chosen parameter (infer_genes_from_transcripts={infer_gene_from_transcript}) as no parent was provided")

        else:
            parent = parents[0]
            # This is an output in some lifton/liftoff cases and they are repeated entries, so ignoring them is appropriate
            if parent == ID:
                if not quiet:
                    print(f"{self.id} Warning: Transcript {ID} has its own id as a parent so this transcript feature was ignored")
                if not skip_orphaned_features:
                    self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))

            elif parent not in self.all_gene_ids:

                self.warnings["transcript_to_inexistent_gene"].add(ID)

                if infer_gene_from_transcript:

                    if unique_t_id not in self.all_transcript_ids:
                        entry["id"] = parent
                        entry["feature"] = "gene"
                        self._add_gene(entry, quiet=quiet)

                        self.all_transcript_ids[unique_t_id] = (ch, parent)
                        self.chrs[ch][parent].transcripts[unique_t_id] = Transcript(unique_t_id, ch, source, ft, strand, start, end, score, ".", parents, attributes)
                        if unique_t_id != ID:
                            if ID in self._transcript_info:
                                self._transcript_info[ID].add(unique_t_id)
                            else:
                                self._transcript_info[ID] = {ID, unique_t_id}
                                print(f"{self.id} Warning: Transcript referred to an inexistent gene so the following {parent} gene was created for {ID} transcript")
                        elif not quiet:
                            print(f"{self.id} Warning: Transcript referred to an inexistent gene so the following {parent} gene was created for {ID} transcript which was renamed to {unique_t_id} to avoid repeat ids")
                    elif not quiet:
                        print(f"{self.id} Warning: {ID} transcript could not be assigned to a gene with the chosen parameter (rename_repeated_id={rename_repeated_id}) as transcript referred to an inexistent gene parent and the provided transcript id was not unique")
                elif not quiet:
                    print(f"{self.id} Warning: {ID} transcript could not be assigned to a gene with the chosen parameter (infer_genes_from_transcripts={infer_gene_from_transcript}) as transcript referred to an inexistent gene parent")

            elif self.all_gene_ids[parent] != ch:
                self.errors["transcript_to_gene_other_chr"].add(ID)
                if not quiet:
                    print(f"{self.id} Error: {ID} transcript refers to a gene in a different chromosome, it could not be assigned to a gene")
                if not skip_orphaned_features:
                    self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))

            else:
            
                if parent not in self._gene_info:

                    if unique_t_id not in self.all_transcript_ids:
                        self.all_transcript_ids[unique_t_id] = (ch, parent)
                        self.chrs[ch][parent].transcripts[unique_t_id] = Transcript(unique_t_id, ch, source, ft, strand, start, end, score, ".", parents, attributes)
                        if unique_t_id != ID:
                            if ID in self._transcript_info:
                                self._transcript_info[ID].add(unique_t_id)
                            else:
                                self._transcript_info[ID] = {ID, unique_t_id}

                            if parent == self.all_transcript_ids[ID][1]:
                                self.warnings["repeat_transcript_same_gene"].add(ID)
                                if not quiet:
                                    print(f"{self.id} Warning: {ID} transcript was renamed to {unique_t_id} to avoid repeat ids for the same {parent} gene")
                            else:
                                self.warnings["repeat_transcript_different_genes"].add(ID)
                                if not quiet:
                                    print(f"{self.id} Warning: {ID} transcript was renamed to {unique_t_id} to avoid repeat transcript ids across genes ({parent} and {self.all_transcript_ids[ID][1]})")
                    elif not quiet:
                        print(f"{self.id} Warning: {ID} transcript could not be assigned to a gene with the chosen parameter (rename_repeated_id={rename_repeated_id}) as the provided transcript id was not unique")

                else:

                    possible_parents = self._gene_info[parent]

                    found = False

                    for parent in possible_parents:

                        if parent not in self.chrs[ch]:
                            continue
                        if self.chrs[ch][parent].start <= start:
                            continue
                        if self.chrs[ch][parent].end >= end:
                            continue

                        found = True
                        if unique_t_id not in self.all_transcript_ids:
                            self.all_transcript_ids[unique_t_id] = (ch, parent)
                            self.chrs[ch][parent].transcripts[unique_t_id] = Transcript(unique_t_id, ch, source, ft, strand, start, end, score, ".", [parent], attributes)
                            if unique_t_id != ID:
                                if ID in self._transcript_info:
                                    self._transcript_info[ID].add(unique_t_id)
                                else:
                                    self._transcript_info[ID] = {ID, unique_t_id}

                                if parent == self.all_transcript_ids[ID][1]:
                                    self.warnings["repeat_transcript_same_gene"].add(ID)
                                    if not quiet:
                                        print(f"{self.id} Warning: {ID} transcript was renamed to {unique_t_id} to avoid repeat ids for the same {parent} gene")
                                else:
                                    self.warnings["repeat_transcript_different_genes"].add(ID)
                                    if not quiet:
                                        print(f"{self.id} Warning: {ID} transcript was renamed to {unique_t_id} to avoid repeat transcript ids across genes ({parent} and {self.all_transcript_ids[ID][1]})")

                        elif not quiet:
                            print(f"{self.id} Warning: {ID} transcript could not be assigned to a gene with the chosen parameter (rename_repeated_id={rename_repeated_id}) as the provided transcript id was not unique")
                        break

                    if not found and not quiet:
                        if unique_t_id != ID:
                            print(f"{self.id} Error: {ID} transcript, renamed to {unique_t_id}, could not be assigned to any gene. Possibly due to unforseen id clash issue")
                        else:
                            print(f"{self.id} Error: {ID} transcript could not be assigned to any gene. Possibly due to unforseen id clash issue")
                        if not skip_orphaned_features:
                            self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))

    def _add_subfeature(self, entry, ft_level, infer_gene_and_transcript_from_subfeatures:bool=False, only_infer_if_none_of_the_parents_exist:bool=True, skip_orphaned_features:bool=False, liftover_exception:bool=True, quiet:bool=False):
        ID = entry["id"]
        ch = entry["ch"]
        source = entry["source"]
        ft = entry["feature"]
        strand = entry["strand"]
        start = entry["start"]
        end = entry["end"]
        phase = entry["phase"]
        score = entry["score"]
        attributes = entry["attributes"]
        parents = entry["parents"]

        if start == end:
            self.warnings[f"1bp_{ft_level}"].add(ID)

        if not parents:
            self.warnings["subfeature_with_no_parent"].add(ID)

            if ft_level != "CDS":

                if not quiet:
                    print(f"{self.id} Warning: No parent provided so the following {ft} subfeature {ID} could not be assigned to any transcript")

                if not skip_orphaned_features:
                    self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))
            
            elif infer_gene_and_transcript_from_subfeatures:

                created_gene = f"{ID}_gene"
                created_transcript = f"{ID}_transcript"

                if created_gene not in self.all_gene_ids:
                    if not quiet and ft != "nucleotide_to_protein_match":
                        print(f"{self.id} Warning: {ft} {ID} subfeature was assigned to created transcript {created_transcript} and created gene {created_gene}")
                    entry["id"] = created_gene
                    entry["feature"] = "gene"
                    self._add_gene(entry, quiet=quiet)

                if created_transcript not in self.all_transcript_ids:
                    entry["id"] = created_transcript
                    entry["feature"] = "transcript"
                    entry["parents"] = [created_gene]
                    self._add_transcript(entry, quiet=quiet, skip_orphaned_features=skip_orphaned_features)

                self.chrs[ch][created_gene].transcripts[created_transcript].temp_CDSs.append(Feature(ID, ch, source, ft, strand, start, end, score, phase, parents, attributes))

            else:

                if not quiet:
                    print(f"{self.id} Warning: No parent provided and infer_gene_and_transcript_from_subfeatures={infer_gene_and_transcript_from_subfeatures} so the following {ft} subfeature {ID} could not be assigned to a created transcript")

                if not skip_orphaned_features:
                    self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))

        true_orphans = True

        for parent in parents:
            if parent in self.all_transcript_ids:
                true_orphans = False
            elif parent in self.all_gene_ids:
                true_orphans = False

        for parent in parents:
            if parent not in self._transcript_info:
                if parent in self.all_transcript_ids:
                    if self.all_transcript_ids[parent][0] != ch:
                        self.errors["subfeature_to_transcript_other_chr"].add(ID)
                        if not quiet:
                            print(f"{self.id} Error: {ID} subfeature refers to a transcript in a different chromosome, it could not be assigned to its parent")
                        if not skip_orphaned_features:
                            self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))
                        continue

                    gene_parent = self.all_transcript_ids[parent][1]

                    if ft_level == "CDS":
                        self.chrs[ch][gene_parent].transcripts[parent].temp_CDSs.append(Feature(ID, ch, source, ft, strand, start, end, score, phase, parents, attributes))
                    elif ft_level == "exon":
                        self.chrs[ch][gene_parent].transcripts[parent].exons.append(Exon(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))
                    elif ft_level == "UTR":
                        self.chrs[ch][gene_parent].transcripts[parent].temp_UTRs.append(UTR(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))
                    else:
                        self.chrs[ch][gene_parent].transcripts[parent].miRNAs.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))

                elif parent in self.all_gene_ids:
                    if self.all_gene_ids[parent] != ch:
                        self.errors["subfeature_to_gene_other_chr"].add(ID)
                        if not quiet:
                            print(f"{self.id} Error: {ID} subfeature refers to a gene in a different chromosome, it could not be assigned to its parent")
                        if not skip_orphaned_features:
                            self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))
                        continue

                    # cases where subfeature directly references a pseudogene creating a pseudotranscript for the pseudogene
                    if self.chrs[ch][parent].pseudogene:
                        if parent[0:4] == "gene":
                            pseudo_t = f"pseudo_t_{parent[5:]}"
                        else:
                            pseudo_t = f"pseudo_t_{parent}"

                        if self.chrs[ch][parent].transcripts == {}:
                            # PSEUDOGENES referred to by transcript subfeatures are given a single pseudotranscript
                            self.chrs[ch][parent].transcripts[pseudo_t] = Transcript(pseudo_t, self.chrs[ch][parent].ch, self.chrs[ch][parent].source, "pseudotranscript", self.chrs[ch][parent].strand, self.chrs[ch][parent].start, self.chrs[ch][parent].end, self.chrs[ch][parent].score, ".", [parent])
                            self.all_transcript_ids[pseudo_t] = (ch, parent)
                            
                        if pseudo_t in self.chrs[ch][parent].transcripts:
                            if ft_level == "CDS":
                                self.chrs[ch][parent].transcripts[pseudo_t].temp_CDSs.append(Feature(ID, ch, source, ft, strand, start, end, score, phase, [pseudo_t], attributes))
                            elif ft_level == "exon":
                                self.chrs[ch][parent].transcripts[pseudo_t].exons.append(Exon(ID, ch, source, ft, strand, start, end, score, ".", [pseudo_t], attributes))
                            elif ft_level == "UTR":
                                self.chrs[ch][parent].transcripts[pseudo_t].temp_UTRs.append(UTR(ID, ch, source, ft, strand, start, end, score, ".", [pseudo_t], attributes))
                            else:
                                self.chrs[ch][parent].transcripts[pseudo_t].miRNAs.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", [pseudo_t], attributes))
                        else:
                            if not skip_orphaned_features:
                                self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))
                            if not quiet:
                                print(f"{self.id} Error: {parent} pseudogene already had a transcript which {ft} subfeature {ID} ignores")

                    # cases where subfeature directly points to a gene that is not a pseudogene
                    else:
                        # gene without transcripts pointed to then auto-create a transcript
                        if self.chrs[ch][parent].transcripts == {}:
                            if not infer_gene_and_transcript_from_subfeatures:
                                self.warnings["subfeature_to_gene"].add(ID)
                                if not skip_orphaned_features:
                                    self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))
                                if not quiet:
                                    print(f"{self.id} Warning: {ft} subfeature {ID} references {parent} which is a gene, but since infer_gene_and_transcript_from_subfeatures is False, the gene and transcript were not auto-created.")
                            else:
                                count = 1
                                auto_t_id = f"{parent}_t{count}"
                                while auto_t_id in self.all_transcript_ids:
                                    count += 1
                                    auto_t_id = f"{parent}_t{count}"

                                gene_obj = self.chrs[ch][parent]
                                self.all_transcript_ids[auto_t_id] = (ch, parent)
                                self.chrs[ch][parent].transcripts[auto_t_id] = Transcript(
                                    auto_t_id, gene_obj.ch, gene_obj.source, "mRNA",
                                    gene_obj.strand, gene_obj.start, gene_obj.end,
                                    gene_obj.score, ".", [parent])
                                if not quiet:
                                    print(f"{self.id} Warning: {ft} subfeature {ID} references {parent} gene which has no transcripts, auto-created transcript {auto_t_id}")
                                self.warnings["subfeature_to_gene"].add(ID)
                        # correctly linking the subfeature to the single transcript that exists
                        if len(self.chrs[ch][parent].transcripts) == 1:
                            temp_id = list(self.chrs[ch][parent].transcripts.keys())[0]
                            if ft_level == "CDS":
                                self.chrs[ch][parent].transcripts[temp_id].temp_CDSs.append(Feature(ID, ch, source, ft, strand, start, end, score, phase, [temp_id], attributes))
                            elif ft_level == "exon":
                                self.chrs[ch][parent].transcripts[temp_id].exons.append(Exon(ID, ch, source, ft, strand, start, end, score, ".", [temp_id], attributes))
                            elif ft_level == "UTR":
                                self.chrs[ch][parent].transcripts[temp_id].temp_UTRs.append(UTR(ID, ch, source, ft, strand, start, end, score, ".", [temp_id], attributes))
                            else:
                                self.chrs[ch][parent].transcripts[temp_id].miRNAs.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", [temp_id], attributes))  
                        elif len(self.chrs[ch][parent].transcripts) > 1:
                            if not quiet:
                                print(f"{self.id} Error: {ft} subfeature {ID} references {parent} gene which is not a pseudogene and has multiple transcripts")
                            self.errors["subfeature_to_gene"].add(ID)

                # parent not found within created genes or transcripts
                else:
                    
                    if liftover_exception and self.liftover:
                        if not quiet:
                            print(f"{self.id} Warning: {ft} subfeature {ID} references {parent} which is not found in the annotation, possibly not transfered during liftover")
                        self.warnings["missing_subfeature_parent_liftover"].add(ID)
                        if not skip_orphaned_features:
                            self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))
                        continue

                    else:
                        if not infer_gene_and_transcript_from_subfeatures:
                            self.warnings["missing_subfeature_parent"].add(ID)
                            if not skip_orphaned_features:
                                self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))
                            if not quiet:
                                print(f"{self.id} Warning: {ft} subfeature {ID} references {parent} which is not found in the annotation")

                        elif only_infer_if_none_of_the_parents_exist and not true_orphans:
                            self.warnings["missing_subfeature_parent"].add(ID)
                            if not skip_orphaned_features:
                                self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))
                            print(f"{self.id} Warning: {ft} subfeature {ID} references {parent} which is not found in the annotation, although at least one of the listed parents was found")

                        else:
                            self.warnings["missing_subfeature_parent"].add(ID)
                            if not quiet:
                                print(f"{self.id} Warning: {ft} subfeature {ID} references {parent} which is not found in the gff and both transcript and gene will be generated")

                            inferred_g_id = f"{parent}_gene"
                            count = 0
                            while inferred_g_id in self.all_gene_ids:
                                count += 1
                                inferred_g_id = f"{parent}_{count}_gene"

                            entry["id"] = inferred_g_id
                            entry["feature"] = "gene"
                            self._add_gene(entry, quiet=quiet)

                            entry["id"] = parent
                            entry["feature"] = "transcript"
                            entry["parents"] = [inferred_g_id]
                            self._add_transcript(entry, quiet=quiet, skip_orphaned_features=skip_orphaned_features)

                            if true_orphans:
                                temp_parents = parents
                            else:
                                temp_parents = [parent]


                            if ft_level == "CDS":
                                self.chrs[ch][inferred_g_id].transcripts[parent].temp_CDSs.append(Feature(ID, ch, source, ft, strand, start, end, score, phase, temp_parents, attributes))
                            elif ft_level == "exon":
                                self.chrs[ch][inferred_g_id].transcripts[parent].exons.append(Exon(ID, ch, source, ft, strand, start, end, score, ".", temp_parents, attributes))
                            elif ft_level == "UTR":
                                self.chrs[ch][inferred_g_id].transcripts[parent].temp_UTRs.append(UTR(ID, ch, source, ft, strand, start, end, score, ".", temp_parents, attributes))
                            else:
                                self.chrs[ch][inferred_g_id].transcripts[parent].miRNAs.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", temp_parents, attributes)) 

            # Cases where the transcript parent was repeated by id in the gff
            else:
                
                possible_parents = self._transcript_info[parent]

                found = False

                for parent in possible_parents:

                    gene_parent = self.all_transcript_ids[parent][1]

                    if gene_parent not in self.chrs[ch]:
                        continue
                    if self.chrs[ch][parent].start <= start:
                        continue
                    if self.chrs[ch][parent].end >= end:
                        continue

                    found = True

                    if ft_level == "CDS":
                        self.chrs[ch][gene_parent].transcripts[parent].temp_CDSs.append(Feature(ID, ch, source, ft, strand, start, end, score, phase, [parent], attributes))
                    elif ft_level == "exon":
                        self.chrs[ch][gene_parent].transcripts[parent].exons.append(Exon(ID, ch, source, ft, strand, start, end, score, ".", [parent], attributes))
                    elif ft_level == "UTR":
                        self.chrs[ch][gene_parent].transcripts[parent].temp_UTRs.append(UTR(ID, ch, source, ft, strand, start, end, score, ".", [parent], attributes))
                    else:
                        self.chrs[ch][gene_parent].transcripts[parent].miRNAs.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", [parent], attributes)) 

                    break

                if not found and not quiet:
                    print(f"{self.id} Error: {ID} {ft} feature could not be assigned to any transcript. Possibly due to unforseen id clash issue")
                    if not skip_orphaned_features:
                        self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, ".", parents, attributes))

    def _get_unique_transcript_id(self, t_id):
        unique_id = t_id
        count = 0
        while unique_id in self.all_transcript_ids:
            count += 1
            unique_id = f"{t_id}_{count}"

        return unique_id

    def _get_unique_gene_id(self, id):
        unique_id = id
        count = 0
        while unique_id in self.all_gene_ids:
            count += 1
            unique_id = f"{id}_{count}"

        return unique_id

    def copy(self):
        return copy.deepcopy(self)

    @property
    def summary(self) -> dict:
        return self.stats.data
    
    def update(self, original_annotation:Annotation|None=None, rename_features:list=[], keep_existing_ids_if_derived_from_base_id:bool=False, genome:Genome|None=None, define_synteny:bool=False, sort_processes:int=1, quiet:bool=False, consider_polycistronic:bool=False, consider_read_utrs:bool=False, collapse_exons:bool=True, collapse_CDSs:bool=True, standardise_features:bool=False):
        start_time = time.time()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False

        batch_size = 1000
        count = 0

        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;62mUpdating {self.id} genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;62m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                count += 1
                if count >= batch_size:
                    progress_bar.update(count)
                    count = 0
                g.collapse_subfeatures(exons=collapse_exons, CDSs=collapse_CDSs)
                for t in g.transcripts.values():
                    t.update(quiet=quiet, consider_polycistronic=consider_polycistronic, consider_read_utrs=consider_read_utrs)
                    if t.polycistronic == "no":
                        continue
                    elif t.polycistronic == "maybe":
                        self.warnings["possible_policistronic_transcript"].add(t.id) 
                    elif t.polycistronic == "yes":
                        self.warnings["multiple_CDSs_per_transcript"].add(t.id)
                g.update()

        if count > 0:
            progress_bar.update(count)
        
        progress_bar.close()
        self.update_features(standardise=standardise_features, quiet=quiet)
        
        if rename_features != []:
            self.rename_ids(features=rename_features, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id, quiet=quiet)
        self.remove_missing_transcript_parent_references()
        self.homogenise_parents_for_shared_exons_utrs()
        self.correct_gene_transcript_and_subfeature_coordinates()
        if not self.sorted:
            self.sort_genes(processes=sort_processes)
        if define_synteny:
            self.define_synteny(original_annotation=original_annotation, sort_processes=sort_processes)
        if self.liftover and original_annotation != None:
            for g_id in original_annotation.all_gene_ids:
                #improve at some point, slow with massive lists
                if g_id not in self.all_gene_ids:
                    self.unmapped.append(g_id)
        self.stats.update(genome=genome)

        self.update_suffixes()

        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nWhole update process for {self.id} annotation object took {round(lapse/60, 1)} minutes\n")

    def update_suffixes(self, quiet:bool=True):
        if not quiet:
            print(f"\nUpdating suffixes for {self.id}")
        self.feature_suffix = ""
        self.suffix = ""

        if self.combined:
            self.feature_suffix += "_combined"
        if self.small_cds_removed:
            self.feature_suffix += "_minus_small_CDSs"
        if self.coding_removed:
            self.feature_suffix += "_minus_coding"
        if self.non_coding_removed:
            self.feature_suffix += "_minus_non_coding"
        if self.transposable_removed:
            self.feature_suffix += "_minus_TE"
        if self.non_transposable_removed:
            self.feature_suffix += "_minus_non_TE"

        if self.full_renamed_ids:
            self.feature_suffix += "_full_renamed_ids"

        self.suffix = self.feature_suffix

        if self.confrenamed:
            self.suffix += "_confrenamed"
        if self.featurecounts:
            self.suffix += "_fcounts"
        if self.dapmod:
            self.suffix += "_dapmod"
        elif self.dapfit:
            self.suffix += "_dapfit"
        if self.symbols_added:
            self.suffix += "_plus_symbols"
        if not quiet:
            print(f"\nUpdated suffixes for {self.id}")

    def update_features(self, standardise=False, quiet:bool=True):
        if not quiet:
            print(f"\nUpdating features for {self.id}")

        self.features = {}
        for genes in self.chrs.values():
            for g in genes.values():
                if standardise:
                    g.feature = default_features_r[g.feature]
                if g.feature not in self.features:
                    self.features[g.feature] = 1
                else:
                    self.features[g.feature] += 1
                for t in g.transcripts.values():
                    if standardise:
                        if t.feature not in default_noncoding_transcripts:
                            t.feature = "mRNA"
                    if t.feature not in self.features:
                        self.features[t.feature] = 1
                    else:
                        self.features[t.feature] += 1
                    for e in t.exons:
                        if standardise:
                            e.feature = "exon"
                        if e.feature not in self.features:
                            self.features[e.feature] = 1
                        else:
                            self.features[e.feature] += 1
                    for c in t.CDSs.values():
                        if standardise:
                            c.feature = "CDS"
                            for cs in c.CDS_segments:
                                cs.feature = "CDS"
                        if c.feature not in self.features:
                            self.features[c.feature] = 1
                        else:
                            self.features[c.feature] += 1

        for ft in self.atypical_features:
            if ft.feature not in self.features:
                self.features[ft.feature] = 1
            else:
                self.features[ft.feature] += 1
        for ft in self.orphaned_features:
            if ft.feature not in self.features:
                self.features[ft.feature] = 1
            else:
                self.features[ft.feature] += 1

        if not quiet:
            print(f"\nUpdated features for {self.id}")

    def mark_transposable_element_genes(self, TE_genes_file):
        TE_genes = set()
        f_in = open(TE_genes_file)
        for line in f_in:
            line = line.strip()
            if line != "":
                TE_genes.add(line)
        f_in.close()
        for genes in self.chrs.values():
            for g in genes.values():
                if g.id in TE_genes:
                    g.transposable = True
        self.update()

    def mark_rRNA_transcripts(self, rRNA_transcripts_file, clean:bool=True):
        rRNA_transcripts = set()
        f_in = open(rRNA_transcripts_file)
        for line in f_in:
            line = line.strip()
            if line != "":
                rRNA_transcripts.add(line)
        f_in.close()
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    if t.id in rRNA_transcripts:
                        t.feature = "rRNA"
                        t.CDSs = {}
                        t.coding = False
                        t.frame = "."
                        t.blast_hits = []
                        t.main = False
                        t.temp_CDSs = []
                        t.temp_UTRs = []
        if clean:
            self.remove_other_mRNA_transcripts_from_rRNA_genes()
            self.update(rename_features=["transcript", "CDS", "exon", "UTR"])
        else:
            self.update()
    
    def remove_other_mRNA_transcripts_from_rRNA_genes(self):
        for chrom, genes in self.chrs.items():
            for g in genes.values():
                rRNA_positive = False
                mRNA_positive = False
                for t in g.transcripts.values():
                    if t.feature == "rRNA":
                        rRNA_positive = True
                    elif t.feature == "mRNA":
                        mRNA_positive = True
                if rRNA_positive == True and mRNA_positive == True:
                    mRNA_transcripts_to_remove = []
                    for t in g.transcripts.values():
                        if t.feature == "mRNA":
                            mRNA_transcripts_to_remove.append(t.id)
                    for t_id in mRNA_transcripts_to_remove:
                        del self.chrs[chrom][g.id].transcripts[t_id]   

    def correct_gene_transcript_and_subfeature_coordinates(self, quiet:bool=True):
        if not quiet:
            print(f"Correcting feature coordinates for {self.id}")

        for genes in self.chrs.values():
            for g in genes.values():
                earliest_start = None
                latest_end = None
                for t in g.transcripts.values():
                    if t.exons:
                        for c in t.CDSs.values():
                            if c.start < t.exons[0].start:
                                if not quiet:
                                    print(f"Warning: {c.id} start should not be earlier than for first {t.id} exon, proceeding to fix {self.id}")
                                t.exons[0].start = c.start
                                t.exons[0].update_size()
                                self.sorted = False
                            if c.end > t.exons[-1].end:
                                if not quiet:
                                    print(f"Warning: {c.id} end should not extend beyond the last {t.id} exon, proceeding to fix {self.id}")
                                t.exons[-1].end = c.end
                                t.exons[-1].update_size()
                                self.sorted = False

                        if t.exons[0].start < t.start:
                            if not quiet:
                                print(f"First exon start should not be earlier than for {t.id}, proceeding to fix {self.id}")
                        elif t.exons[0].start > t.start:
                            if not quiet:
                                print(f"First exon should not start later than {t.id}, proceeding to fix {self.id}")
                        if t.exons[-1].end < t.end:
                            if not quiet:
                                print(f"Last exon should not finish earlier than {t.id}, proceeding to fix {self.id}")
                        elif t.exons[-1].end > t.end:
                            if not quiet:
                                print(f"Last exon should not finish later than {t.id}, proceeding to fix {self.id}")
                        if t.start != t.exons[0].start or t.end != t.exons[-1].end:
                            t.start = t.exons[0].start
                            t.end = t.exons[-1].end
                            t.update_size()
                            self.sorted = False

                    if t.start < g.start:
                        if not quiet:
                            print(f"{t.id} start should not be earlier than for {g.id}, proceeding to fix {self.id}")
                        g.start = t.start
                        self.sorted = False
                    if t.end > g.end:
                        if not quiet:
                            print(f"{t.id} end should not extend beyond {g.id}, proceeding to fix {self.id}")
                        g.end = t.end
                        self.sorted = False
                    if earliest_start is None or latest_end is None:
                        earliest_start = t.start
                        latest_end = t.end
                    else:
                        if t.start < earliest_start:
                            earliest_start = t.start
                        if t.end > latest_end:
                            latest_end = t.end

                if earliest_start is not None and latest_end is not None:
                    if g.start != earliest_start or g.end != latest_end:
                        if not quiet:
                            print(f"{g.id} was too long and had to be trimmed to longest transcript ({self.id})")
                        g.start = earliest_start
                        g.end = latest_end
                        g.update_size()
                        self.sorted = False
        if not quiet:
            print(f"Corrected feature coordinates for {self.id}")

    def generate_sequences(self, genome:Genome, just_CDSs:bool=False, quiet:bool=False):
        start_time = time.time()
        for o in self.atypical_features:
            o.generate_sequence(genome)
        for o in self.orphaned_features:
            o.generate_sequence(genome)
        if not just_CDSs:
            for genes in self.chrs.values():
                for g in genes.values():
                    g.generate_sequence(genome)
                    for t in g.transcripts.values():
                        t.generate_sequence(genome)
                        for c in t.CDSs.values():
                            c.generate_sequence(genome)
            self.generated_all_sequences = True
            self.contains_all_sequences = True 
        else:
            for genes in self.chrs.values():
                for g in genes.values():
                    for t in g.transcripts.values():
                        for c in t.CDSs.values():
                            c.generate_sequence(genome)
        self.generated_CDS_sequences = True
        self.contains_CDS_sequences = True
        self.generated_protein_sequences = True
        self.contains_protein_sequences = True

        now = time.time()
        lapse = now - start_time

        if not quiet:
            print(f"\nGenerating sequences for {self.id} annotation object took {round(lapse/60, 1)} minutes\n")

    def clear_sequences(self, just_hard=False, keep_proteins:bool=False, quiet:bool=True):
        start_time = time.time()
        for o in self.atypical_features:
            o.clear_sequence(just_hard=just_hard)
        for o in self.orphaned_features:
            o.clear_sequence(just_hard=just_hard)
        for genes in self.chrs.values():
            for g in genes.values():
                g.clear_sequence(just_hard=just_hard)
                for t in g.transcripts.values():
                    t.clear_sequence(just_hard=just_hard)
                    for c in t.CDSs.values():
                        c.clear_sequence(just_hard=just_hard, keep_proteins=keep_proteins)
        if not keep_proteins:
            self.contains_protein_sequences = False
        self.contains_CDS_sequences = False
        self.contains_all_sequences = False

        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"Clearing sequences of {self.id} annotation object took {round(lapse, 1)} seconds\n")

    def generate_promoters(self, genome:Genome, promoter_size:int=2000, promoter_type:str = "standard", generate_sequence:bool=False):
        """
        promoter_type (str): Defines the reference point for the promoters.
            - standard (default): Promoter based on 'promoter_size' is generated upstream of the transcript's start site (TSS)
            - upstream_ATG : Promoter based on 'promoter_size' is generated upstream of the main CDS's start codon (ATG). If no CDS, falls back to standard.
            - standard_plus_up_to_ATG : Promoter based on 'promoter_size' is generated upstream of the transcript's start site (TSS) and any gene sequence up to the start codon (ATG) is also added. If no CDS, falls back to standard.
        """
        self.promoter_types = promoter_type
        self.promoter_size = promoter_size

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    t.generate_promoter(promoter_size, genome.scaffolds[t.ch].size, promoter_type)
                    if generate_sequence:
                        t.promoter.generate_sequence(genome)

    def return_random_gene_ids(self, number:int=1, to_avoid:list=[], coding:bool=True):
        random_ids = []
        while len(random_ids) < number:
            r = random.choice(list(self.all_gene_ids.keys()))
            ch = self.all_gene_ids[r]
            if coding:
                if self.chrs[ch][r].coding and r not in random_ids and r not in to_avoid:
                    random_ids.append(r)
            elif r not in random_ids and r not in to_avoid:
                random_ids.append(r)

        return random_ids

    def combine_transcripts(self, genome:Genome, respect_non_coding:bool=False):
        for genes in self.chrs.values():
            for g in genes.values():
                g.combine_transcripts(genome, respect_non_coding=respect_non_coding)
        self.sorted = False
        self.update(rename_features=["transcript", "CDS", "exon", "UTR"])
        self.combined = True
        self.update_suffixes()

    def sort_genes(self, processes:int=2, quiet:bool=True, noisy:bool=False):
        if not quiet:
            print(f"\nSorting genes for {self.id}")
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                        bar_format=(
            f'\033[38;2;210;180;140m\033[1mSorting {self.id} genes:\033[0m '
            '{percentage:3.0f}%|'
            f'\033[38;2;210;180;140m{{bar}}\033[0m| '
            '{n}/{total} [{elapsed}<{remaining}]'))

        if processes > 1:
            if not quiet:
                print(f"Parallel mode chosen for {self.id}")
            with Pool(processes=processes) as pool:
                if not quiet:
                    print(f"Parallel sorting started pool for {self.id}")
                for chrom, sorted_genes in pool.starmap(sort_and_update_genes, [(chrom, genes) for chrom, genes in self.chrs.items()]):
                    if not quiet and noisy:
                        print(f"Parallel sorting {chrom} genes for {self.id}")
                    self.chrs[chrom] = sorted_genes
                    progress_bar.update(len(sorted_genes))
        else:
            if not quiet:
                print(f"Sequential mode chosen for {self.id}")
            for chrom, genes in self.chrs.items():
                sorted_genes = sorted(genes.values())
                if not quiet and noisy:
                    print(f"Sequential sorting {chrom} genes for {self.id}")
                self.chrs[chrom] = {g.id: g for g in sorted_genes}
                progress_bar.update(len(sorted_genes))

        progress_bar.close()
        self.sorted = True
        if not quiet:
            print(f"Sorted genes for {self.id}")

    def define_synteny(self, original_annotation:Annotation|None=None, sort_processes:int=1, quiet:bool=True):
        if not quiet:
            print(f"\nDefining synteny for {self.id} annotation genes")
        start_time = time.time()
        if not self.sorted:
            self.sort_genes(processes=sort_processes)
        # defining synteny
        for genes in self.chrs.values():
            gene_list = list(genes.values())
            for n, g in enumerate(gene_list):
                g.synteny_order = n
                # This works even if only a single gene has been annotated
                # in a chromosome or scaffold
                if len(genes) == 1:
                    g.previous_gene = False
                    g.next_gene = False
                elif n == 0:
                    g.previous_gene = False
                    g.next_gene = gene_list[n+1].id
                elif n != (len(genes) - 1):
                    g.previous_gene = gene_list[n-1].id
                    g.next_gene = gene_list[n+1].id
                else:
                    g.previous_gene = gene_list[n-1].id
                    g.next_gene = False

        if self.liftover and original_annotation:
            for genes in self.chrs.values():
                for g_id, g in genes.items():
                    # this extra bit is for extra liftover copies
                    if g_id not in original_annotation.all_gene_ids:
                        g_id = "_".join(g_id.split("_")[:-1])
                        if g_id not in original_annotation.all_gene_ids:
                            continue

                    g.old_previous_gene = original_annotation.chrs[original_annotation.all_gene_ids[g_id]][g_id].previous_gene
                    g.old_next_gene = original_annotation.chrs[original_annotation.all_gene_ids[g_id]][g_id].next_gene                
                    g.old_synteny_order = original_annotation.chrs[original_annotation.all_gene_ids[g_id]][g_id].synteny_order

                    if g.old_previous_gene == g.previous_gene and g.old_next_gene == g.next_gene:
                        g.conserved_synteny = True
                    else: 
                        g.conserved_synteny = False
        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nDefining synteny for {self.id} annotation genes took {round(lapse, 1)} seconds\n")

    def homogenise_parents_for_shared_exons_utrs(self):
        for genes in self.chrs.values():
            for g in genes.values():
                exon_groups: dict[tuple[int, int, str, str], list[Exon]] = {}
                for t in g.transcripts.values():
                    for e in t.exons:
                        key = (e.start, e.end, e.ch, e.strand)
                        if key not in exon_groups:
                            exon_groups[key] = [e]
                        else:
                            exon_groups[key].append(e)

                for group in exon_groups.values():
                    if len(group) > 1:
                        all_parents = set()
                        for e in group:
                            all_parents.update(e.parents)
                        merged_parents = sorted(all_parents)
                        for e in group:
                            e.parents = merged_parents[:]
                    else:
                        group[0].parents.sort()

                utr_groups: dict[tuple[int, int, str, str], list[UTR]] = {}
                for t in g.transcripts.values():
                    for c in t.CDSs.values():
                        for u in c.UTRs:
                            key = (u.start, u.end, u.ch, u.strand)
                            if key not in utr_groups:
                                utr_groups[key] = [u]
                            else:
                                utr_groups[key].append(u)

                for group in utr_groups.values():
                    if len(group) > 1:
                        all_parents = set()
                        for u in group:
                            all_parents.update(u.parents)
                        merged_parents = sorted(all_parents)
                        for u in group:
                            u.parents = merged_parents[:]
                    else:
                        group[0].parents.sort()

        self.shared_exons = True
        self.shared_UTRs = True

    def single_parent_for_exons_utrs(self):
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    for e in t.exons:
                        e.parents = [t.id]
                    for c in t.CDSs.values():
                        for u in c.UTRs:
                            u.parents = [t.id]

        self.shared_exons = False
        self.shared_UTRs = False

    def add_aliases(self, overlap_threshold:int=6):
        for genes in self.chrs.values():
            for g in genes.values():
                for name, hits in g.overlaps.items():
                    if name == "other":
                        for hit in hits:
                            if hit.score >= overlap_threshold:
                                if hit.id not in g.aliases:
                                    g.aliases.append(hit.id)
    
    def clear_aliases(self):
        for genes in self.chrs.values():
            for g in genes.values():
                g.aliases = []

    def CDS_to_CDS_segment_ids(self, override:bool=False):
        repeat_CDS_segment_id = False

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    for j, c in enumerate(t.CDSs.values()):
                        for x1, cs1 in enumerate(c.CDS_segments):
                            for x2, cs2 in enumerate(c.CDS_segments):
                                if x1 == x2:
                                    continue
                                if cs1.id == cs2.id:
                                    repeat_CDS_segment_id = True
                                    break

        if repeat_CDS_segment_id or override:
            for genes in self.chrs.values():
                for g in genes.values():
                    for t in g.transcripts.values():
                        count = 1
                        for j, c in enumerate(t.CDSs.values()):
                            if not c.main:
                                count += 1
                                c.id = f"{t.id}_CDS{count}"
                            else:
                                c.id = f"{t.id}_CDS1"
                            for x, cs in enumerate(c.CDS_segments):
                                cs.id = f"{c.id}_{x+1}"

    def CDS_segment_to_CDS_ids(self, override:bool=False):
        common_protein_CDS_ids = True

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    for c in t.CDSs.values():
                        for x1, cs1 in enumerate(c.CDS_segments):
                            for x2, cs2 in enumerate(c.CDS_segments):
                                if x1 == x2:
                                    continue
                                if cs1 != cs2:
                                    common_protein_CDS_ids = False
                                    break

        if not common_protein_CDS_ids or override:
            for genes in self.chrs.values():
                for g in genes.values():
                    for t in g.transcripts.values():
                        count = 1
                        for c in t.CDSs.values():
                            if not c.main:
                                count += 1
                                c.id = f"{t.id}_CDS{count}"
                            else:
                                c.id = f"{t.id}_CDS1"
                            for cs in c.CDS_segments:
                                cs.id = c.id

    def merge(self, other:Annotation, max_cds_overlap:int|float=100, max_exon_overlap:int|float=100, max_gene_overlap:int|float=100, features_to_rename:list=["gene", "transcript", "CDS", "exon", "UTR"], quiet:bool=False):
        """
        Priority is given to self annotation
        """
        start_time = time.time()
        self.update(quiet=quiet)
        other.update(quiet=quiet)

        if max_cds_overlap != 100 or max_exon_overlap != 100 or max_gene_overlap != 100:
            self.overlaps.detect(other, quiet=quiet)

        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False

        progress_bar = tqdm(total=len(other.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[38;2;46;204;113m\033[1mMerging {other.id} and {self.id} annotations:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[38;2;46;204;113m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        if self.merged:
            if "_..._" in self.name:
                first_name = self.name.split("_..._")[0]
            else:
                first_name = self.name.split("_merged_")[0]
            self.name = f"{first_name}_..._{other.name}"
        else:
            self.name = f"{self.name}_merged_{other.name}"
        if self.genome != None:
            self.id = f"{self.name}_on_{self.genome}"
        else:
            self.id = self.name
        count = 0
        if max_cds_overlap == 100 and max_exon_overlap == 100 and max_gene_overlap == 100:
            for chr, genes in other.chrs.items():
                if chr not in self.chrs:
                    self.chrs[chr] = {}
                for g in genes.values():
                    progress_bar.update(1)
                    count = 0
                    temp_id = g.id
                    while temp_id in self.all_gene_ids:
                        count += 1
                        temp_id = f"{g.id}_{count}"
                    self.chrs[chr][temp_id] = g.copy()
                    self.chrs[chr][temp_id].id = temp_id
                    self.all_gene_ids[temp_id] = chr

        else:
            for chr, genes in other.chrs.items():
                if chr not in self.chrs:
                    self.chrs[chr] = {}
                for g in genes.values():
                    progress_bar.update(1)
                    count = 0
                    if g.overlaps["other"] == []:
                        temp_id = g.id
                        while temp_id in self.all_gene_ids:
                            count += 1
                            temp_id = f"{g.id}_{count}"
                        self.chrs[chr][temp_id] = g.copy()
                        self.chrs[chr][temp_id].id = temp_id
                        self.all_gene_ids[temp_id] = chr
                    else:
                        cds_scores = [0]
                        exon_scores = [0]
                        gene_scores = [0]
                        for o in g.overlaps["other"]:
                            if o.min_CDS_percent != None:
                                cds_scores.append(o.min_CDS_percent)
                            if o.min_exon_percent != None:
                                exon_scores.append(o.min_exon_percent)
                            if o.min_gene_percent != None:
                                gene_scores.append(o.min_gene_percent)

                        check_cdss = any(ov.CDSs_in_both for ov in g.overlaps["other"])
                        check_exons = any(ov.exons_in_both for ov in g.overlaps["other"])

                        # these conditionals are placed here since some genes may lack CDSs, others may even lack exons and hence some of the thresholds may not be applicable
                        
                        if check_cdss:
                            if max(cds_scores) <= max_cds_overlap and max(exon_scores) <= max_exon_overlap and max(gene_scores) <= max_gene_overlap:
                                temp_id = g.id
                                while temp_id in self.all_gene_ids:
                                    count += 1
                                    temp_id = f"{g.id}_{count}"
                                self.chrs[chr][temp_id] = g.copy()
                                self.chrs[chr][temp_id].id = temp_id
                                self.all_gene_ids[temp_id] = chr
                        
                        elif check_exons:
                            if max(exon_scores) <= max_exon_overlap and max(gene_scores) <= max_gene_overlap:
                                temp_id = g.id
                                while temp_id in self.all_gene_ids:
                                    count += 1
                                    temp_id = f"{g.id}_{count}"
                                self.chrs[chr][temp_id] = g.copy()
                                self.chrs[chr][temp_id].id = temp_id
                                self.all_gene_ids[temp_id] = chr

                        elif max(gene_scores) <= max_gene_overlap:
                            temp_id = g.id
                            while temp_id in self.all_gene_ids:
                                count += 1
                                temp_id = f"{g.id}_{count}"
                            self.chrs[chr][temp_id] = g.copy()
                            self.chrs[chr][temp_id].id = temp_id
                            self.all_gene_ids[temp_id] = chr

        self.merged = True
        self.sorted = False
        self.generated_all_sequences = False
        self.contains_all_sequences = False        
        self.generated_CDS_sequences = False
        self.contains_CDS_sequences = False
        self.generated_protein_sequences = False
        self.contains_protein_sequences = False
        progress_bar.close()
        self.update_gene_and_transcript_list(quiet=quiet)
        self.update(rename_features=features_to_rename, quiet=quiet)
        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nMerging {self.id} and {other.id} annotations took {round(lapse/60, 1)} minutes")

    def remove_exons_with_unmatched_strand(self, quiet:bool=False):
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    new_e = []
                    for e in t.exons:
                        if e.strand == t.strand:
                            new_e.append(e)
                    t.exons = new_e

        self.remove_transcripts_with_no_exons(quiet=quiet)
        self.remove_genes_with_no_transcripts(quiet=quiet)

    def remove_transcripts_with_no_exons(self, quiet:bool=False):
        transcripts_to_remove = []
        for chrom, genes in self.chrs.items():
            for g in genes.values():
                for t in g.transcripts.values():
                    if t.exons == []:
                        transcripts_to_remove.append((chrom, g.id, t.id))
        for chrom, g_id, t_id in transcripts_to_remove:
            del self.chrs[chrom][g_id].transcripts[t_id]
            if not quiet:
                print(f"{t_id} Warning: transcript with no exons which could have been removed because of strand inconsistencies")
            self.warnings["transcript_with_no_exons"].add(t_id)

    def remove_transcripts(self, to_remove:set, quiet:bool=False):
        for t in to_remove:
            if t in self.all_transcript_ids:
                chrom = self.all_transcript_ids[t][0]
                g_id = self.all_transcript_ids[t][1]
                del self.chrs[chrom][g_id].transcripts[t]
            elif not quiet:
                warnings.warn(f"Transcript level id {t} is not present in annotation {self.id}", category=UserWarning)

        self.remove_genes_with_no_transcripts(quiet=quiet)

    def remove_genes_with_no_transcripts(self, quiet:bool=False):
        genes_to_remove = []
        for chrom, genes in self.chrs.items():
            for g in genes.values():
                if g.transcripts == {}:
                    genes_to_remove.append((chrom, g.id))
        for chrom, g_id in genes_to_remove:
            del self.chrs[chrom][g_id]
            if not quiet:
                print(f"{g_id} Warning: gene with no transcripts which could have been removed because of strand inconsistencies of its exons")
            self.warnings["gene_with_no_transcripts"].add(g_id)

        self.update_gene_and_transcript_list(quiet=quiet)

    def remove_missing_transcript_parent_references(self, quiet:bool=True):
        if not quiet:
            print(f"Removing missing transcript parent references for {self.id} annotation.")
        self.remove_transcripts_with_no_exons(quiet=quiet)
        self.remove_genes_with_no_transcripts(quiet=quiet)

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    for e in t.exons:
                        new_parents = []
                        for p in e.parents:
                            if p in self.all_transcript_ids:
                                if g.transcripts[p].strand == e.strand:
                                    new_parents.append(p)
                        e.parents = new_parents.copy()
                        e.parents.sort()
                    for c in t.CDSs.values():
                        new_parents = []
                        for p in c.parents:
                            if p in self.all_transcript_ids:
                                new_parents.append(p)
                        c.parents = new_parents.copy()
                        c.parents.sort()
                        for cs in c.CDS_segments:
                            new_parents = []
                            for p in cs.parents:
                                if p in self.all_transcript_ids:
                                    new_parents.append(p)
                            cs.parents = new_parents.copy()
                            cs.parents.sort()
        if not quiet:
            print(f"Removed missing transcript parent references for {self.id} annotation.")

    def rework_CDSs(self, genome:Genome, override:bool=True, low_memory:bool=True, coding_ratio_threshold:float=0.8, quiet:bool=False):
        start_time = time.time()
        if low_memory:
            self.clear_sequences()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mReworking {self.id} CDSs:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
    
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                for t in g.transcripts.values():
                    if t.coding and not override:
                        continue
                    t.generate_sequence(genome, low_memory)
                    t.generate_best_protein(genome)
                    t.generate_CDSs_based_on_ORF(low_memory)
                    for c in t.CDSs.values():
                        c.generate_sequence(genome, low_memory)
                    t.update()
                    if t.coding_ratio < coding_ratio_threshold:
                        t.generate_sequence(genome, low_memory)
                        t.generate_best_protein(genome, must_have_stop=False)
                        t.generate_CDSs_based_on_ORF(low_memory)
                        for c in t.CDSs.values():
                            c.generate_sequence(genome, low_memory)
                    t.update()

        progress_bar.close()
        self.update(rename_features=["CDS", "exon", "UTR"], quiet=quiet)
        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nReworking CDSs for {self.id} took {round(lapse/60, 1)} minutes")

    def update_gene_and_transcript_list(self, quiet:bool=True):
                # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False

        total = 0
        for genes in self.chrs.values():
            total += len(genes.keys())

        progress_bar = tqdm(total=total, disable=disable,
                                bar_format=(
                    f'\033[1;95mUpdating {self.id} annotation gene and transcript lists:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;95m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        self.all_gene_ids = {}
        self.all_transcript_ids = {}
        for chr, genes in self.chrs.items():
            for g in genes.values():
                progress_bar.update(1)
                self.all_gene_ids[g.id] = chr
                for t in g.transcripts.values():
                    self.all_transcript_ids[t.id] = (chr, g.id)
        progress_bar.close()

    def make_alternative_transcripts_into_genes(self, quiet:bool=False):
        new_chrs = {}
        for chrom, genes in self.chrs.items():
            new_genes = {}
            for g in genes.values():
                new_transcripts = {}
                for t in g.transcripts.values():
                    if t.main:
                        new_transcripts[t.id] = t.copy()
                    else:
                        g_new = g.copy()
                        count = 0
                        while g_new.id in self.all_gene_ids:
                            count += 1
                            g_new.id = f"{g.id}_{count}"
                        
                        g_new.transcripts = {t.id : t.copy()}
                        g_new.start = t.start
                        g_new.end = t.end
                        new_genes[g_new.id] = g_new.copy()
                        self.all_gene_ids[g_new.id] = chrom

                g.transcripts = new_transcripts.copy()
            new_chrs[chrom] = new_genes.copy()

        for nc in new_chrs:
            for g in new_chrs[nc].values():
                self.chrs[nc][g.id] = g.copy()

        del new_chrs

        self.update_gene_and_transcript_list(quiet=quiet)
        self.update(rename_features=["gene", "transcript", "CDS", "exon", "UTR"], quiet=quiet)

    def rename_ids(self, custom_path:str="", features:list[str]=["gene", "transcript", "CDS", "exon", "UTR"], keep_existing_ids_if_derived_from_base_id:bool=False, remove_point_suffix:bool=False, strip_gene_tag:bool=False, keep_subfeature_numbers:bool=False, cds_segment_ids:bool=False, prefix:str="", suffix:str="", spacer:int=100, sep:str="_", g_id_digits:int=5, t_id_digits:int=3, correspondences:bool=False, quiet:bool=False):

        acceptable_features = ["gene", "transcript", "CDS", "exon", "UTR"]

        for f in features:
            if f not in acceptable_features:
                raise ValueError(f"Incorrect feature '{f}' chosen. Select from: {acceptable_features}.")
        
        if features == []:
            raise ValueError(f"Rename ids was called but no feature levels were chosen. Select from: {acceptable_features}.")
        
        if prefix:
            if keep_existing_ids_if_derived_from_base_id or features != acceptable_features or remove_point_suffix:
                ignored_options = []
                if keep_existing_ids_if_derived_from_base_id:
                    ignored_options.append("keep_existing_ids_if_derived_from_base_id")
                if features != acceptable_features:
                    ignored_options.append("features")
                if remove_point_suffix:
                    ignored_options.append("remove_point_suffix")
                if strip_gene_tag:
                    ignored_options.append("strip_gene_tag")
                warnings.warn(f"Providing a prefix '{prefix}' means all features will be renamed based on the prefix, the following provided options are to be ignored: {ignored_options}.", category=UserWarning)

        elif suffix:
            warnings.warn(f"Provided suffix={suffix} will have no effect as no prefix was provided and custom renaming will therefore be skipped.", category=UserWarning)

        if cds_segment_ids and "CDS" not in features:
            warnings.warn("CDS features will be changed if need be since cds_segment_ids have been requested.", category=UserWarning)

        start_time = time.time()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        changed_features = set()
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[38;2;156;42;42m\033[1mRenaming {self.id} Gene Ids:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[38;2;156;42;42m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

        correspondences_d = {}

        for genes in self.chrs.values():
            g_count = 0
            for g in genes.values():
                progress_bar.update(1)
                g_count += spacer

                if "gene" in features or prefix:
                    g.rename(g_count, sep=sep, digits=g_id_digits, prefix=prefix, suffix=suffix, base_id_as_id=strip_gene_tag, remove_point_suffix=remove_point_suffix)
                    if g.renamed:
                        changed_features.add("gene")

                    correspondences_d[g.original_id] = g.id

                tmains = 0

                base_id_present = False
                base_id_missing = False

                for t in g.transcripts.values():
                    if t.main:
                        tmains += 1
                    if g.base_id in t.id:
                        base_id_present = True
                    else:
                        base_id_missing = True
                    
                if base_id_present and base_id_missing:
                    warnings.warn(f"{self.id} annotation gene {g.original_id} has a mix of transcript id formats and renaming errors could occur!", category=UserWarning)
                
                if tmains > 1:
                    raise ValueError(f"There should not be more than one main transcript for {g.original_id} in {self.id} annotation.")


                t_count = 1

                for t in g.transcripts.values():

                    t.parents = [g.id]

                    if not t.main:
                        t_count += 1

                    if "transcript" in features or prefix:

                        if prefix or (base_id_present and base_id_missing):
                            t.rename(base_id=g.base_id, sep=sep, count=t_count, digits=t_id_digits, keep_numbering=keep_subfeature_numbers)
                        else:
                            t.rename(base_id=g.base_id, sep=sep, count=t_count, digits=t_id_digits, keep_numbering=keep_subfeature_numbers, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id)

                        if t.renamed:
                            changed_features.add("transcript")

                    cmains = 0
                    base_id_present = False
                    base_id_missing = False

                    for c in t.CDSs.values():
                        if c.main:
                            cmains += 1
                        if g.base_id in c.id:
                            base_id_present = True
                        else:
                            base_id_missing = True

                        for cs in c.CDS_segments:
                            if g.base_id in cs.id:
                                base_id_present = True
                            else:
                                base_id_missing = True

                    if base_id_present and base_id_missing:
                        warnings.warn(f"{self.id} annotation transcript {t.original_id} has a mix of CDS id formats and renaming errors could occur!", category=UserWarning)

                    if cmains > 1:
                        raise ValueError(f"There shouldn't be more than one main CDS for transcript {t.original_id} in {self.id} annotation.")

                    c_count = 1
                    for c in t.CDSs.values():
                        c.parents = [t.id]

                        if not c.main:
                            c_count += 1

                        if "CDS" in features or prefix:

                            if prefix or (base_id_present and base_id_missing):
                                c.rename(base_id=t.id, base_gene_id=g.base_id, count=c_count, sep=sep, digits=t_id_digits, keep_numbering=keep_subfeature_numbers, keep_existing_ids_if_derived_from_base_id=False, cds_segment_ids=cds_segment_ids)
                            else:
                                c.rename(base_id=t.id, base_gene_id=g.base_id, count=c_count, sep=sep, digits=t_id_digits, keep_numbering=keep_subfeature_numbers, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id, cds_segment_ids=cds_segment_ids)

                            if c.renamed:
                                changed_features.add("CDS")

                        for cs in c.CDS_segments:
                            cs.parents = [t.id]


                    base_id_present = False
                    base_id_missing = False

                    for e in t.exons:
                        e.parents = [t.id]
                        if g.base_id in e.id:
                            base_id_present = True
                        else:
                            base_id_missing = True

                    if base_id_present and base_id_missing:
                        warnings.warn(f"{self.id} annotation transcript {t.original_id} has a mix of exon id formats and renaming errors could occur!", category=UserWarning)


                    base_id_present = False
                    base_id_missing = False

                    for c in t.CDSs.values():
                        for u in c.UTRs:
                            u.parents = [t.id]
                            if g.base_id in u.id:
                                base_id_present = True
                            else:
                                base_id_missing = True

                    if base_id_present and base_id_missing:
                        warnings.warn(f"{self.id} annotation transcript {t.original_id} has a mix of UTR id formats and renaming errors could occur!", category=UserWarning)


                if "exon" in features or prefix:
                    g.rename_exons(sep=sep, digits=t_id_digits, keep_numbering=keep_subfeature_numbers, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id)

                    if g.renamed_exons:
                        changed_features.add("exon")

                    if g.renamed_utrs:
                        changed_features.add("UTR")


        progress_bar.close()

        self.homogenise_parents_for_shared_exons_utrs()
        self.update_keys(quiet=quiet)
        self.update_gene_and_transcript_list(quiet=quiet)

        self.renamed_features = changed_features

        if features == ["gene", "transcript", "CDS", "exon", "UTR"]:
            self.full_renamed_ids = True

        if correspondences:

            export_folder = Path(custom_path or self.path) / "out_gffs"
            export_folder.mkdir(parents=True, exist_ok=True)
            export_folder = str(export_folder) + "/"

            if "gene" in changed_features:
                out_text = ["old_gene_id\tnew_gene_id"]
                for k, v in correspondences_d.items():
                    out_text.append(f"{v}\t{k}")
                f_out = open(f"{export_folder}{self.id}{self.feature_suffix}_rename_eqs.tsv", "w", encoding="utf-8")
                out_text = "\n".join(out_text)
                f_out.write(out_text)
                f_out.close()

            else:
                warnings.warn(f"Correspondences on gene ids were requested, however gene ids remained unchanged.", category=UserWarning)

        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nRenaming {self.id} ids with prefix='{prefix}', changing={self.renamed_features} features took {round(lapse/60, 1)} minutes")        

    def update_keys(self, quiet:bool=True):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys())*3, disable=disable,
                                bar_format=(
                    f'\033[1;95mUpdating {self.id} dictionary keys:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;95m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))


        # fixing CDS keys
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                for t in g.transcripts.values():
                    new_CDSs = {}
                    for c in t.CDSs.values():
                        if c.id not in new_CDSs:
                            new_CDSs[c.id] = c.copy()
                        else:
                            print(f"Error: Repeated {c.id} CDS in {t.id} transcript in {self.id} annotation when updating keys.")
                    t.CDSs = new_CDSs.copy()

        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                new_transcripts = {}
                for t in g.transcripts.values():
                    if t.id not in new_transcripts:
                        new_transcripts[t.id] = t.copy()
                    else:
                        print(f"Error: Repeated {t.id} transcript in {g.id} gene in {self.id} annotation when updating keys.")
                g.transcripts = new_transcripts.copy()


        for chrom, genes in self.chrs.items():
            new_genes = {}
            for g in genes.values():
                progress_bar.update(1)
                new_genes[g.id] = g.copy()
            self.chrs[chrom] = new_genes.copy()
            del new_genes

        progress_bar.close()

    def create_featurecounts_ids(self):
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    t.gene_id = g.id
                    for c in t.CDSs.values():
                        for cs in c.CDS_segments:
                            cs.gene_id = g.id
                        for u in c.UTRs:
                            u.gene_id = g.id

    def create_gtf_attributes(self):

        for genes in self.chrs.values():
            for g in genes.values():

                g.gtf_attributes = [f'gene_id "{g.id}"']
            
                if g.symbols != []:
                    name_string = ",".join(g.symbols)
                    g.gtf_attributes.append(f'gene_name "{name_string}"')
                elif g.names != []:
                    name_string = ",".join(g.names)
                    g.gtf_attributes.append(f'gene_name "{name_string}"')

                for t in g.transcripts.values():
                    t.gtf_attributes = [f'gene_id "{g.id}"', f'transcript_id "{t.id}"']
                    for x, e in enumerate(t.exons):
                        e.gtf_attributes = [f'gene_id "{g.id}"', f'transcript_id "{t.id}"', f'exon_number "{x+1}"']

                    for c in t.CDSs.values():
                        c.gtf_attributes = [f'gene_id "{g.id}"', f'transcript_id "{t.id}"']
                        for cs in c.CDS_segments:
                            cs.gtf_attributes = [f'gene_id "{g.id}"', f'transcript_id "{t.id}"']
                        for u in c.UTRs:
                            u.gtf_attributes = [f'gene_id "{g.id}"', f'transcript_id "{t.id}"']

    def add_blast_hits(self, source, blastfile, mode:str="protein"):
        """
        Blast results could be at a gene, transcript or protein level. In future,
        other features can also be easily included.

        This function for now assumes one hit per protein!
        """
        accepted_levels = ["protein", "CDS", "transcript", "gene"]

        if mode in accepted_levels:
            hits = {}
            f_in = open(blastfile, encoding="utf-8")
            for line in f_in:
                line = line.strip().split("\t")
                ID = line[0]
                evalue = float(line[-2])
                bitscore = float(line[-1])
                if ID not in hits:
                    hits[ID] = BlastHit(source, bitscore, evalue)
                else:
                    print(f"Warning: for now the add blast hits method does not accept more than one blast hit per ID ({ID})")
            f_in.close()
            if mode == "protein":
                for protein_id in self.protein_equivalences:
                    chrom, g, t, c = self.all_protein_ids[protein_id]
                    prot = self.chrs[chrom][g].transcripts[t].CDSs[c].protein
                    if protein_id in hits:
                        if prot is not None:
                            prot.blast_hits.append(hits[protein_id])
                    for protein_id_copy in self.protein_equivalences[protein_id]:
                        chrom, g, t, c = self.all_protein_ids[protein_id_copy]
                        prot = self.chrs[chrom][g].transcripts[t].CDSs[c].protein
                        if protein_id_copy in self.chrs[chrom] and protein_id in hits:
                            if prot is not None:
                                prot.blast_hits.append(hits[protein_id])
            else:
                "Adding blast hits is not available for gene, transcripts and CDSs yet."
        else:
            print(f"Warning: {mode} chosen is not in accepted list of choices=['protein', 'CDS', 'transcript', 'gene]")   

    def remove_chromosomes_from_header(self):
        new_header = []
        ft_to_keep = set(self.chrs)
        
        for line in self.gff_header:
            if line.startswith("##sequence-region"):
                keep_header = False
                for ft in ft_to_keep:
                    if ft in line:
                        keep_header = True
                if keep_header:
                    new_header.append(line)

        self.gff_header = new_header.copy()

    def subset(self, chosen_features, gene_cap:int=3000, common_chromosomes:set|None=None, min_genes:int=1500, quiet:bool=False):

        initial_chosen_features = chosen_features.copy()

        for chosen_feature in chosen_features:
            if chosen_feature not in self.chrs:
                raise ValueError(f"Chosen scaffold/chromosome {chosen_feature} is not in {self.name} genome. Choose from '{self.chrs.keys()}'")

        if common_chromosomes is None:
            total_chromosomes = set(self.chrs)

        else:
            total_chromosomes = common_chromosomes.copy()

        if min_genes > 0:

            num_genes_in_chosen_features = 0
            for ft in chosen_features:
                num_genes_in_chosen_features += len(self.chrs[ft])

            remaining_to_chose_from = total_chromosomes - chosen_features
            chr_cap_overriden = False

            while num_genes_in_chosen_features < min_genes and remaining_to_chose_from:
                chr_cap_overriden = True
                
                chosen_features.add(random.choice(list(remaining_to_chose_from)))

                remaining_to_chose_from = total_chromosomes - chosen_features

                num_genes_in_chosen_features = 0
                for ft in chosen_features:
                    num_genes_in_chosen_features += len(self.chrs[ft])

            if chr_cap_overriden:
                print(f"Chromosome/scaffold cap of {len(initial_chosen_features)} was overriden by min_genes = {min_genes} parameter as not enough genes were present in {initial_chosen_features}. The final selection includes {len(chosen_features)} features: {chosen_features}")

        features_to_remove = set(self.chrs) - chosen_features
        genes_to_keep_per_chromosome = math.ceil(gene_cap / len(chosen_features))

        self.remove_chromosomes(features_to_remove, update=False, quiet=quiet)

        genes_to_remove = set()

        total_deficit = 0

        for genes in self.chrs.values():
            deficit = genes_to_keep_per_chromosome - len(genes)
            if deficit < 0:
                deficit = 0
            total_deficit += deficit

        for genes in self.chrs.values():

            gene_list = list(genes)
            surplus = len(genes) - genes_to_keep_per_chromosome

            if surplus > 0 :

                contribution_to_cover_deficit = min(surplus, total_deficit)

                final_genes_to_keep = genes_to_keep_per_chromosome + contribution_to_cover_deficit
                
                if len(gene_list) > final_genes_to_keep:
                    genes_to_keep_sample = set(random.sample(gene_list, final_genes_to_keep))
                    
                    genes_to_remove_from_this_chr = set(gene_list) - genes_to_keep_sample
                    genes_to_remove.update(genes_to_remove_from_this_chr)
                    
                total_deficit -= contribution_to_cover_deficit

        if genes_to_remove:
            self.remove_genes(genes_to_remove, quiet=quiet)
        else:
            warnings.warn(f"The cap value {gene_cap} was not enforced as there are not enough genes in the subset chromosomes in annotation {self.id}.", category=UserWarning)

        self.update(quiet=quiet)

        return chosen_features

    def filter_by_rna_class(self, rna_classes=['mRNA'], quiet:bool=False):

        transcript_to_remove = set()

        for genes in self.chrs.values():

            for g in genes.values():

                for t in g.transcripts.values():

                    if t.feature not in rna_classes:

                        transcript_to_remove.add(t.id)

        self.remove_transcripts(transcript_to_remove, quiet=quiet)
    
        self.update(quiet=quiet)
        

    def remove_chromosomes(self, features_to_remove:set, update:bool=True, quiet:bool=False):
        if features_to_remove:
            for ft in features_to_remove:
                del self.chrs[ft]
            self.remove_chromosomes_from_header()

        if update:
            self.update(quiet=quiet)

    def remove_genes(self, to_remove:set|None=None, override_rescue:bool=False, quiet:bool=False):

        if to_remove is None:
            to_remove = set()

        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False

        total_count = 0

        for genes in self.chrs.values():
            total_count += len(genes)
        
        progress_bar = tqdm(total=total_count, disable=disable,
                                bar_format=(
                    f'\033[1;91mRemoving {self.id} genes:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        
        for gene in to_remove:
            if gene in self.all_gene_ids:
                chrom = self.all_gene_ids[gene]
                self.chrs[chrom][gene].remove = True
                if override_rescue:
                    self.chrs[chrom][gene].rescue = False
            else:
                warnings.warn(f"Gene {gene} is not present in annotation {self.id}.", category=UserWarning)

        genes_to_remove = set()
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if g.remove and not g.rescue:
                    genes_to_remove.add(g.id)

        for g_id in genes_to_remove:
            chrom = self.all_gene_ids[g_id]
            del self.chrs[chrom][g_id]
            del self.all_gene_ids[g_id]
        progress_bar.close()
        
        self.remove_missing_genes_in_overlaps(quiet=quiet)
        self.update_gene_and_transcript_list(quiet=quiet)

    def remove_missing_genes_in_overlaps(self, quiet:bool=True):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False

        total_count = 0
        for genes in self.chrs.values():
            total_count += len(genes)

        progress_bar = tqdm(total=total_count, disable=disable,
                                bar_format=(
                    f'\033[1;91mRemoving {self.id} missing genes in overlaps:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                new_overlaps = []
                for o in g.overlaps["self"]:
                    if o.id in self.all_gene_ids:
                        new_overlaps.append(o)
                g.overlaps["self"] = new_overlaps
        progress_bar.close()              

    def remove_duplicate_transcripts(self, quiet:bool=False):
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        progress_bar = tqdm(total=len(self.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[1;91mRemoving repeat transcripts per gene of {self.id}:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[1;91m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                to_eliminate = []
                for t1 in g.transcripts.values():
                    if t1.id not in to_eliminate:
                        for t2 in g.transcripts.values():
                            if t2.id not in to_eliminate:
                                if t1.id == t2.id:
                                    continue
                                if t1.almost_equal(t2):
                                    to_eliminate.append(t1.id)
                to_eliminate = set(to_eliminate)
                for t_eliminate in to_eliminate:
                    del g.transcripts[t_eliminate]
        progress_bar.close()
        self.update(rename_features=["transcript", "CDS", "exon", "UTR"], quiet=quiet)

    def make_alternative_genes_into_transcripts(self, quiet:bool=False):

        correspondence = {}

        to_remove = []
        for chrom, genes in self.chrs.items():
            for g in genes.values():
                if (not g.remove or g.rescue) and (g.id not in correspondence):
                    for o in g.overlaps["self"]:
                        if o.score >= 5:

                            if o.id in correspondence:

                                correspondence_gene = correspondence[o.id]
                                correspondence[g.id] = correspondence_gene
                                to_remove.append((g.id, chrom))
                                for t in self.chrs[chrom][g.id].transcripts.values():
                                    self.chrs[chrom][correspondence_gene].transcripts[t.id] = t.copy()

                            else:

                                correspondence[o.id] = g.id
                                to_remove.append((o.id, chrom))
                                for t in self.chrs[chrom][o.id].transcripts.values():
                                    g.transcripts[t.id] = t.copy()

        for g, chrom in to_remove:
            if g in self.chrs[chrom]:
                del self.chrs[chrom][g]

        self.update(rename_features=["gene", "transcript", "CDS", "exon", "UTR"], quiet=quiet)

    def remove_genes_with_small_CDSs(self, CDS_threshold:int=200, quiet:bool=False):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                remove = False
                for t in g.transcripts.values():
                    if t.main:
                        for c in t.CDSs.values():
                            if c.main:
                                if c.size < CDS_threshold:
                                    remove = True
                                    removed_any = True
                if remove:
                    g.rescue = False
                    g.remove = True
        self.remove_genes(quiet=quiet)

        if removed_any:
            self.small_cds_removed = True
        self.update(quiet=quiet)

    def remove_TE_genes(self, quiet:bool=False):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if g.transposable:
                    g.rescue = False
                    g.remove = True
                    removed_any = True

        self.remove_genes(quiet=quiet)

        if removed_any:
            self.transposable_removed = True
        self.update(quiet=quiet)

    def remove_non_TE_genes(self, quiet:bool=False):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if not g.transposable:
                    g.rescue = False
                    g.remove = True
                    removed_any = True

        self.remove_genes(quiet=quiet)
        
        if removed_any:
            self.non_transposable_removed = True

        self.update(quiet=quiet)

    def remove_non_coding_genes_and_transcripts(self, quiet:bool=False):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if g.noncoding == True and g.coding == False:
                    g.rescue = False
                    g.remove = True
                    removed_any = True

        self.remove_genes(quiet=quiet)

        removed_any = self.remove_non_coding_transcripts_from_coding_genes(removed_any, False, quiet=quiet)

        if removed_any:
            self.non_coding_removed = True

        self.update(quiet=quiet)

    def remove_coding_genes_and_transcripts(self, quiet:bool=False):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if g.noncoding == False and g.coding == True:
                    g.rescue = False
                    g.remove = True
                    removed_any = True

        self.remove_genes(quiet=quiet)

        removed_any = self.remove_coding_transcripts_from_non_coding_genes(removed_any, False, quiet=quiet)

        if removed_any:
            self.coding_removed = True

        self.update(quiet=quiet)

    def remove_coding_transcripts_from_non_coding_genes(self, removed_any:bool=False, update=True, quiet:bool=False):
        transcripts_to_remove = []

        for chrom, genes in self.chrs.items():
            for g in genes.values():
                if g.coding:
                    for t in g.transcripts.values():
                        if t.coding == True:
                            transcripts_to_remove.append((chrom, g.id, t.id))
        
        for chrom, g_id, t_id in transcripts_to_remove:
            del self.chrs[chrom][g_id].transcripts[t_id]

        if transcripts_to_remove != []:
            removed_any = True

        if removed_any:
            self.coding_removed = True
        
        self.overlaps.clear()

        if update:
            self.update(rename_features=["transcript", "CDS", "exon", "UTR"], quiet=quiet)

        return removed_any

    def remove_non_coding_transcripts_from_coding_genes(self, removed_any:bool=False, update=True, quiet:bool=False):
        transcripts_to_remove = []

        for chrom, genes in self.chrs.items():
            for g in genes.values():
                if g.noncoding:
                    for t in g.transcripts.values():
                        if t.coding == False:
                            transcripts_to_remove.append((chrom, g.id, t.id))
        
        for chrom, g_id, t_id in transcripts_to_remove:
            del self.chrs[chrom][g_id].transcripts[t_id]

        if transcripts_to_remove != []:
            removed_any = True

        if removed_any:
            self.non_coding_removed = True

        self.overlaps.clear()

        if update:
            self.update(rename_features=["transcript", "CDS", "exon", "UTR"], quiet=quiet)

        return removed_any

    def clear_gene_names_and_symbols(self, quiet:bool=False):
        for genes in self.chrs.values():
            for g in genes.values():
                g.names = []
                g.symbols = []
                g.synonyms = []
        self.update(quiet=quiet)

    def remove_genes_without_symbols(self, quiet:bool=False):
        for genes in self.chrs.values():
            for g in genes.values():
                if g.symbols == []:
                    g.remove = True
        self.remove_genes(quiet=quiet)
        self.update(quiet=quiet)

    def rename_chromosomes(self, equivalences, dap:bool=False, quiet:bool=False):

        renamed_scaffolds = False
        for old, new in equivalences.items():
            if new != old and old in self.chrs:
                renamed_scaffolds = True

        self.chrs = {equivalences.get(k, k): v for k, v in self.chrs.items()}

        for chrom, genes in self.chrs.items():
            for g in genes.values():
                g.ch = chrom
                for t in g.transcripts.values():
                    t.ch = chrom
                    if hasattr(t, "promoter"):
                        t.promoter.ch = chrom
                    for c in t.CDSs.values():
                        c.ch = chrom
                        for cs in c.CDS_segments:
                            cs.ch = chrom
                        if hasattr(c, "protein"):
                            if c.protein != None:
                                c.protein.ch = chrom
                        for u in c.UTRs:
                            u.ch = chrom
                    for e in t.exons:
                        e.ch = chrom
                    for i in t.introns:
                        i.ch = chrom

        for o in self.atypical_features:
            if o.ch in equivalences:
                o.ch = equivalences[o.ch]

        for o in self.orphaned_features:
            if o.ch in equivalences:
                o.ch = equivalences[o.ch]

        if dap:
            self.dapfit = True
            if renamed_scaffolds:
                self.dapmod = True
        elif renamed_scaffolds:
            self.confrenamed = True

        self.update_suffixes(quiet=quiet)

    def add_gene_symbols_pseudogenes(self, file_path:str, just_gene_names:bool=True, clear:bool=True, header:bool=False, sep:str="\t"):
        if clear:
            self.clear_gene_names_and_symbols()

        if header:
            if file_path.endswith(".xlsx"):
                df = pd.read_excel(file_path, skiprows=1, dtype=str)
            else:
                df = pd.read_csv(file_path, skiprows=1, sep=sep, dtype=str)
        else:
            if file_path.endswith(".xlsx"):
                df = pd.read_excel(file_path, dtype=str)
            else:
                df = pd.read_csv(file_path, sep=sep, dtype=str)

        df = df.fillna("")

        pseudogene_col_exists = False
        synonym_col_exists = False

        if "pseudogene" in df.columns:
            pseudogene_col_exists = True
        
        if "gene name synonym(s)" in df.columns:
            synonym_col_exists = True

        for _, row in df.iterrows():
            gene_name = str(row.iloc[0])
            gene_id = str(row.iloc[1])
            if pseudogene_col_exists:
                pseudogene = str(row["pseudogene"])
            if synonym_col_exists:
                gene_synonyms = str(row["gene name synonym(s)"]).split("; ")
            ch = self.all_gene_ids[gene_id]
            self.chrs[ch][gene_id].symbols.append(gene_name)
            if not just_gene_names:
                if synonym_col_exists:
                    self.chrs[ch][gene_id].synonyms += gene_synonyms
                if pseudogene_col_exists:
                    if pseudogene == "pseudogene" or pseudogene == "Pseudogene":
                        self.chrs[ch][gene_id].pseudogene = True
                    elif pseudogene == "gene" or pseudogene == "Gene":
                        self.chrs[ch][gene_id].pseudogene = False

        self.symbols_added = True

    def add_gene_symbols(self, file_path:str, clear:bool=True, header:bool=False, sep:str="\t"):

        if clear:
            self.clear_gene_names_and_symbols()

        if header:
            if file_path.endswith(".xlsx"):
                df = pd.read_excel(file_path, skiprows=1, dtype=str)
            else:
                df = pd.read_csv(file_path, skiprows=1, sep=sep, dtype=str)
        else:
            if file_path.endswith(".xlsx"):
                df = pd.read_excel(file_path, dtype=str)
            else:
                df = pd.read_csv(file_path, sep=sep, dtype=str)

        df = df.fillna("")

        for i in range(len(df)):
            gene_name = str(df.iat[i, 1])
            gene_id = str(df.iat[i, 0])
            ch = self.all_gene_ids[gene_id]
            self.chrs[ch][gene_id].symbols.append(gene_name)

        self.symbols_added = True

    def release(self, name, id, source_name, id_prefix, spacer:int=10, suffix:str="", custom_path:str="", tag:str=".gff3", skip_atypical_fts:bool=True, main_only:bool=False, UTRs:bool=True, clear_aliases=True, extra_attributes=False, quiet:bool=False):
        if clear_aliases:
            self.clear_aliases()
        self.name = name
        self.id = id
        for genes in self.chrs.values():
            for g in genes.values():
                g.source = source_name
                for t in g.transcripts.values():
                    t.source = source_name
                    for e in t.exons:
                        e.source = source_name
                    for c in t.CDSs.values():
                        c.source = source_name
                        for cs in c.CDS_segments:
                            cs.source = source_name
                        for u in c.UTRs:
                            u.source = source_name

        self.update(quiet=quiet)
        self.rename_ids(prefix=id_prefix, spacer=spacer, suffix=suffix, features=["gene", "transcript", "CDS", "exon", "UTR"], quiet=quiet)
        self.update(quiet=quiet)
        self.export.gff(custom_path=custom_path, extra_attributes=extra_attributes, tag=tag, skip_atypical_fts=skip_atypical_fts, main_only=main_only, UTRs=UTRs, quiet=quiet)

    def __str__(self):
        return str(self.id)
