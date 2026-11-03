"""
Created on Fri Oct 7 15:09:46 2022

Module defining several genomic classes.

@authors: David Navarro, Antonio Santiago
"""

from __future__ import annotations
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from .genome import Genome
    from .utils.gtf_gff import GffEntry

import copy
import random
import time
import pandas as pd
import os
import warnings
import math
import re

from multiprocessing import Pool
from pathlib import Path

from .feature import Feature
from .genome import classify_compartment
from .gene import Gene
from .transcript import Transcript
from .subfeatures import Exon, UTR, CDS
from .hits import BlastHit
from .utils.genefunctions import sort_and_update_genes, TAXONOMY_ORGANELLE_CODES, resolve_taxonomy_tables, NCBI_GENETIC_CODES
from .utils.misc import read_file_with_fallback, open_file, start_progress_bar, unambiguous_numbers
from .utils.gtf_gff import parse_gff_parts, convert_gtf_to_gff3, detect_file_format, gtf_quote, INT_NAME, INT_ALIAS, MAPPED_ATTRIBUTES, GffEntry
from .annotation_components.stats import AnnotationStats
from .annotation_components.export import AnnotationExport
from .annotation_components.motifs import AnnotationMotifs
from .annotation_components.overlaps import AnnotationOverlaps
from .annotation_components.redundancy import AnnotationRedundancy
from .conf import default_noncoding_transcripts, default_features_r, default_features, default_codons

class Annotation():

    WARNING_KEYS = (
        "1bp_exon", "1bp_CDS", "1bp_UTR", "1bp_miRNA", "decreasing_coordinates",
        "missing_subfeature_parent", "transcript_to_inexistent_gene", "transcript_with_no_parent",
        "missing_subfeature_parent_liftover", "overlapping_CDS_segments",
        "possible_polycistronic_transcript", "independent_CDS_segment_ids",
        "overlapping_exons", "abutting_exons",
        "transcript_with_no_exons",
        "gene_with_no_transcripts", "subfeature_with_no_parent", "subfeature_to_gene", "repeat_transcript_different_genes", "repeat_transcript_same_gene",
        "feature_exceeds_scaffold_length", "chromosome_not_in_genome",
        "phase_mismatch_across_intron",
        "cds_segment_coordinate_mismatch", "cds_exceeds_exon", "cds_segment_in_intron",
        "utr_exceeds_exon", "exon_transcript_boundary_mismatch", "transcript_exceeds_gene",
        "gene_exceeds_transcripts", "read_utrs_inconsistent", "read_utrs_differ_from_derived",
        "repeated_subfeature_ids", "self_parented_transcript",
        "codon_added_to_CDS", "codon_not_incorporated",
        "linked_by_shared_attribute", "shared_attribute_not_linked", "repeated_feature_line",
        "pseudogene_line_read_as_transcript", "linked_by_file_order", "linkable_by_file_order"
    )

    # Warnings that describe the state of the annotation object and can therefore be rechecked
    # after modifications; the rest are only detectable while parsing the input file
    RECHECKABLE_WARNING_KEYS = (
        "phase_mismatch_across_intron", "overlapping_exons", "abutting_exons",
        "cds_segment_coordinate_mismatch", "cds_exceeds_exon", "cds_segment_in_intron", "utr_exceeds_exon",
        "exon_transcript_boundary_mismatch", "transcript_exceeds_gene", "gene_exceeds_transcripts",
        "gene_with_no_transcripts", "transcript_with_no_exons",
    )

    @staticmethod
    def _warning_ids_hint(warning_key:str) -> str:
        return f"; to list their ids, run `aegis summary` on this annotation, which exports them in the '{warning_key}' column of its _warnings.csv"

    @classmethod
    def _new_warnings_dict(cls, keys:tuple[str, ...]|None=None) -> dict[str, set[str]]:
        return {key: set() for key in (cls.WARNING_KEYS if keys is None else keys)}

    _stats: AnnotationStats
    _export: AnnotationExport
    _motifs: AnnotationMotifs
    _overlaps: AnnotationOverlaps
    _redundancy: AnnotationRedundancy

    entry: GffEntry
    
    genome: Genome | None
    original_annotation: Annotation | None
    all_gene_ids:dict[str, str]
    # the tuple values consist of transcript id keys consist of (chromosome, gene_id)
    all_transcript_ids:dict[str, tuple[str, str]]
    _gene_info:dict[str, set[str]]
    _transcript_info:dict[str, set[str]]
    _miRNA_info:set[str]
    contains_promoters:bool

    tags:set[str]

    chrs:dict[str, dict[str, Gene]]
    
    bar_colors = ["31", "32", "34", "33", "33", "33"]
    overlapped_annotations: set[str]

    tags_to_detect:set[str] = { "clean", "dapfit", "dapmod", "confrenamed", "plus_symbols", "standardised_features"}
    feature_tags_to_detect:set[str] = {"minus_TE", "minus_non_TE", "minus_coding", "minus_non_coding", "minus_pseudogene", "minus_non_pseudogene", "minus_small_CDS", "combined", "full_renamed_ids"}
    # Tags of renamed chromosomes only change exports that list chromosome names
    chromosome_tags:set[str] = {"confrenamed", "dapmod"}

    def __init__(self, annot_file_path:str, name:str|None=None, genome:Genome|None=None, hard_masked_genome:Genome|None=None, original_annotation:Annotation|None=None, target:bool=False, to_overlap:bool=True, rework_all_CDSs:bool=False, work_out_missing_CDSs:bool=False, fallback_to_trim:bool=False, chosen_chromosomes:tuple[str, ...]|None=None, chosen_coordinates:tuple[int, int]|None=None, sort_processes:int=1, define_synteny=False, rename_features:tuple[str,...]=(), keep_existing_ids_if_derived_from_base_id:bool=False, quiet:bool=False, consider_read_utrs:bool=False, infer_genes_from_transcripts:bool=True, infer_genes_from_subfeatures:bool=True, skip_orphaned_features:bool=True, skip_atypical_features:bool=True, incorporate_and_rename_repeated_ids:bool=True, collapse_exons:bool=True, collapse_CDSs:bool=True, standardise_features:bool=False, remove_missing_transcript_parent_references:bool=False, remove_transcripts_with_no_exons:bool=False, remove_genes_with_no_transcripts:bool=False, remove_genes_with_no_transcripts_even_if_pseudogene:bool=False, rename_source:str="", adjust_internal_shifts:Literal["intra_exon", "all", "none"]|bool="intra_exon", taxonomy:Literal["plant", "vertebrate", "invertebrate", "yeast"]|str="plant", table:int|str=1, auto_organelle_detection:bool=True, mito_table:int|str|None=None, plastid_table:int|str|None=None, mitochondria_chroms:list[str]|tuple[str, ...]|str|None=None, chloroplast_chroms:list[str]|tuple[str, ...]|str|None=None, coordinate_polishing:Literal["none", "contain", "trim"]="contain", coding_ratio_threshold:float=0.7, allow_internal_stops:bool=True, allow_partial:bool=True, enforce_start_codon:bool=True, orf_choice_mode:Literal["longest", "earliest"]="longest", min_codon_len:int=2, recalculate_phases:bool=False, reset_phases_zero:bool=False, initiator_methionine:bool=True, alt_start_methionine:bool=True, list_load_warnings:bool=False, link_by:tuple[str, ...]|list[str]|str|None=("locus_tag",), link_by_file_order:bool=False):
        
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
        self.overlapped_annotations = set()
        self.merged = False
        self.contains_promoters = False

        if coordinate_polishing not in ("none", "contain", "trim"):
            raise ValueError(f"Invalid coordinate polishing mode: '{coordinate_polishing}'. Choose from: 'none', 'contain', 'trim'.")
        self.coordinate_polishing = coordinate_polishing
        self.coding_ratio_threshold = coding_ratio_threshold
        self.allow_internal_stops = allow_internal_stops
        self.allow_partial = allow_partial
        self.enforce_start_codon = enforce_start_codon
        self.orf_choice_mode = orf_choice_mode
        self.min_codon_len = min_codon_len
        self.recalculate_phases = recalculate_phases
        self.reset_phases_zero = reset_phases_zero
        # Followed by every later update unless it is given its own
        self.collapse_exons = collapse_exons
        self.collapse_CDSs = collapse_CDSs
        self._protein_qc_summary = None

        self.adjust_internal_shifts = adjust_internal_shifts
        self.initiator_methionine = initiator_methionine
        self.alt_start_methionine = alt_start_methionine
        self._store_genetic_codes(taxonomy=taxonomy, table=table, mito_table=mito_table, plastid_table=plastid_table, auto_organelle_detection=auto_organelle_detection, mitochondria_chroms=mitochondria_chroms, chloroplast_chroms=chloroplast_chroms)

        # list_load_warnings prints every feature behind the load-time warnings that are otherwise printed as counts
        self.list_load_warnings = list_load_warnings
        self._load_note_counts = {}
        # self-parented transcripts are only resolved once all transcript rows are loaded
        self._defer_self_parented_transcripts = True
        self._self_parented_transcripts = []

        self.genome = genome
        
        if self.genome is not None:
            self.id = f"{self.name}_on_{self.genome.name}"
            Feature._ACTIVE_GENOME = self.genome
        else:
            self.id = self.name

        if hard_masked_genome is not None:
            Feature._ACTIVE_HARD_GENOME = hard_masked_genome

        self.original_annotation = original_annotation

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
        # sequences marked Is_circular=true in a GFF region line
        self.circular_chromosomes = set()

        # optional default tag system for output gffs/gtfs, detected from the file name (not its folders)
        file_name = Path(annot_file_path).name
        self.tags = {tag for tag in self.tags_to_detect if self._has_tag(file_name, tag)}
        self.feature_tags = {tag for tag in self.feature_tags_to_detect if self._has_tag(file_name, tag)}

        self.renamed_features = []

        self.promoter_size = 3000

        self.contains_protein_sequences = False

        self.promoter_type = "standard"

        keys = [
            "repeated_gene_IDs", "1bp_gene", "1bp_transcript", "subfeature_to_gene",
            "transcript_to_more_than_1_gene", "transcript_to_gene_other_chr", "subfeature_to_transcript_other_chr", "subfeature_to_gene_other_chr"
        ]
        # Create a dictionary where every value is an empty set
        self.errors = {key: set() for key in keys}
        
        # Warnings describe the input annotation: they are recorded while loading (including
        # the load-time update) and are not modified by later operations. See recheck_warnings()
        # for the state of the annotation object after modifications.
        self.warnings = self._new_warnings_dict()
        self._recording_input_warnings = True

        del keys
        
        encoding = read_file_with_fallback(self.file)
        file_format = detect_file_format(self.file, encoding)
        self.input_format = file_format # 'gff3' or 'gtf'
        if not quiet:
            print(f"{file_format} file format detected for file='{self.file}'")

        if file_format == 'gtf':
            gff_file = f"{self.name}_{pid}.tmp"
            convert_gtf_to_gff3(self.file, gff_file, encoding, quiet=quiet)
        else:
            gff_file = self.file

        staging = self.load_data(gff_file, encoding=encoding, chosen_chromosomes=chosen_chromosomes, chosen_coordinates=chosen_coordinates, skip_atypical_features=skip_atypical_features, quiet=quiet)

        if isinstance(link_by, str):
            link_by = tuple(k.strip() for k in link_by.split(",") if k.strip())
        self._link_by_shared_attribute(staging, tuple(link_by or ()), create_features=infer_genes_from_subfeatures, quiet=quiet)

        # Gene and transcript lines with neither ID nor Parent (e.g. MITOS: one 'Name=' line per record) are
        # named after their location, so they do not all share the ID '' and overwrite each other
        for stage in ("gene", "transcript"):
            for entry in staging[stage]:
                if not entry.id:
                    entry.id = f"{entry.feature}_{entry.ch}_{entry.start}_{entry.end}"

        # Last resort, opt-in: lines still without Parent join the gene/transcript they follow in the file
        self._link_by_file_order(staging, apply=link_by_file_order, quiet=quiet)

        # Features crossing the origin of circular sequences are read across it (Feature.seq)
        if self.genome is not None:
            for ch in self.circular_chromosomes:
                scaffold = self.genome.get_scaffold(ch)
                if scaffold is not None:
                    scaffold.circular = True

        self._gene_info = {}
        self._transcript_info = {}
        self._miRNA_info = set()

        batch_size = 5000

        subfeatures = {"exon","CDS", "UTR", "miRNA"}

        for n, stage in enumerate(staging):

            if len(staging[stage]) < batch_size:
                corrected_batch_size = len(staging[stage])
            else:
                corrected_batch_size = batch_size

            progress_bar = start_progress_bar(total=len(staging[stage]), description=f"Adding {self.id} {stage}s", colour=Annotation.bar_colors[n], quiet=quiet)
            
            count = 0

            for entry in staging[stage]:
                if stage in subfeatures:
                    self._add_subfeature(entry, ft_level=stage, infer_gene_and_transcript_from_subfeatures=infer_genes_from_subfeatures, skip_orphaned_features=skip_orphaned_features, quiet=quiet)
                elif stage == "transcript":
                    self._add_transcript(entry, rename_repeated_id=incorporate_and_rename_repeated_ids, infer_gene_from_transcript=infer_genes_from_transcripts, skip_orphaned_features=skip_orphaned_features, quiet=quiet)
                else:
                    self._add_gene(entry, rename_repeated_id=incorporate_and_rename_repeated_ids, quiet=quiet)

                count += 1
                if count >= corrected_batch_size:
                    progress_bar.update(count)
                    count = 0

            if count > 0:
                progress_bar.update(count)

            progress_bar.close()
            staging[stage] = []

            if stage == "transcript":
                self._add_self_parented_transcripts(skip_orphaned_features=skip_orphaned_features, quiet=quiet)

        self._incorporate_codons(skip_orphaned_features=skip_orphaned_features, quiet=quiet)

        del self._gene_info
        del self._transcript_info
        del self._miRNA_info
        del staging

        if self.file != gff_file:
            os.remove(gff_file)

        now = time.time()
        lapse = now - start_time

        if genome is not None:
            if genome.dapfit:
                self.tags.add("dapfit")
            else:
                rogue_chromosome_format = False
                for chrom in self.chrs:
                    if chrom.startswith("chr"):
                        number_str = chrom[3:]
                        if not number_str.isdigit():
                            rogue_chromosome_format = True
                            break
                    else:
                        rogue_chromosome_format = True
                        break
                if not rogue_chromosome_format:
                    self.tags.add("dapfit")

        # Validated on the parsed coordinates, before any update or polishing, as a load-time check
        self._genome_validation = self._validate_against_genome(quiet=quiet)

        misc_attributes = False
        for genes in self.chrs.values():
            for g in genes.values():
                if g.misc_attributes:
                    misc_attributes = True
                    break
                for t in g.transcripts.values():
                    if t.misc_attributes:
                        misc_attributes = True
                        break
                    for e in t.exons:
                        if e.misc_attributes:
                            misc_attributes = True
                            break
                    for c in t.CDSs.values():
                        if c.misc_attributes:
                            misc_attributes = True
                            break
                        for u in c.UTRs:
                            if u.misc_attributes:
                                misc_attributes = True
                                break
                    for ft in t.miRNAs:
                        if ft.misc_attributes:
                            misc_attributes = True
                            break
                    if misc_attributes:
                        break
                if misc_attributes:
                    break
            if misc_attributes:
                break
        for ft in self.atypical_features:
            if ft.misc_attributes:
                misc_attributes = True
                break
        for ft in self.orphaned_features:
            if ft.misc_attributes:
                misc_attributes = True
                break

        if not misc_attributes:
            self.tags.add("clean")

        if not quiet:
            for (note, warning_key), count in self._load_note_counts.items():
                print(f"{self.id} Warning: {count} {note}" + ("" if self.list_load_warnings else Annotation._warning_ids_hint(warning_key)))
            print(f"\nCreating {self.id} annotation object took {round(lapse/60, 1)} minutes\n")
        del self._load_note_counts

        self.detect_polycistronic_transcripts(quiet=quiet)

        if consider_read_utrs:
            for t in self.iter_transcripts():
                t.keep_read_UTRs()

        self.update(sort_processes=sort_processes, define_synteny=define_synteny, rename_features=rename_features, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id, quiet=quiet, collapse_exons=collapse_exons, collapse_CDSs=collapse_CDSs, standardise_features=standardise_features, remove_missing_transcript_parent_references=remove_missing_transcript_parent_references, remove_transcripts_with_no_exons=remove_transcripts_with_no_exons, remove_genes_with_no_transcripts=remove_genes_with_no_transcripts, remove_genes_with_no_transcripts_even_if_pseudogene=remove_genes_with_no_transcripts_even_if_pseudogene, coordinate_polishing=coordinate_polishing, recalculate_phases=recalculate_phases, reset_phases_zero=reset_phases_zero)
        self._recording_input_warnings = False

        self._validate_organelle_chroms()

        if (rework_all_CDSs or work_out_missing_CDSs) and genome:
            # rework_CDSs only gives new ids to the CDSs and UTRs it changes; rename_features then applies to all features
            self.rework_CDSs(
                override=rework_all_CDSs,
                coding_ratio_threshold=coding_ratio_threshold,
                fallback_to_trim=fallback_to_trim,
                min_codon_len=min_codon_len,
                quiet=quiet,
                allow_internal_stops=allow_internal_stops,
                allow_partial=allow_partial,
                enforce_start_codon=enforce_start_codon,
                orf_choice_mode=orf_choice_mode,
            )
            self.update(sort_processes=sort_processes, define_synteny=define_synteny, rename_features=rename_features, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id, quiet=quiet, collapse_exons=collapse_exons, collapse_CDSs=collapse_CDSs, standardise_features=standardise_features, remove_missing_transcript_parent_references=remove_missing_transcript_parent_references, remove_transcripts_with_no_exons=remove_transcripts_with_no_exons, remove_genes_with_no_transcripts=remove_genes_with_no_transcripts, remove_genes_with_no_transcripts_even_if_pseudogene=remove_genes_with_no_transcripts_even_if_pseudogene, coordinate_polishing=coordinate_polishing, recalculate_phases=recalculate_phases, reset_phases_zero=reset_phases_zero)

        if rename_source:
            self.rename_source(rename_source)

        if self.genome is not None:
            self._check_genome_compatibility(quiet=quiet)

    def _check_genome_compatibility(self, quiet: bool = False):
        if not self.genome or not getattr(self.genome, "scaffolds", None) or not self.chrs:
            return

        annot_chrs = set(self.chrs.keys())
        genome_scaffolds = set(self.genome.scaffolds.keys())
        overlap = annot_chrs.intersection(genome_scaffolds)

        if len(overlap) == 0:
            # Check if any annotation chromosome matches OriSeqID in genome scaffold descriptions
            desc_oriseq_ids = set()
            for s in self.genome.scaffolds.values():
                d = getattr(s, "description", "")
                m = re.search(r"(?:^|[\s;,])OriSeqID=([^\s;,]+)", d)
                if m:
                    desc_oriseq_ids.add(m.group(1))

            oriseq_overlap = annot_chrs.intersection(desc_oriseq_ids)
            sample_annot = sorted(list(annot_chrs))[:3]
            sample_scaffolds = sorted(list(genome_scaffolds))[:3]

            if oriseq_overlap:
                msg = (
                    f"None of the {len(annot_chrs)} chromosome IDs in annotation '{self.name}' (e.g. {sample_annot}) "
                    f"match the {len(genome_scaffolds)} scaffold IDs in genome '{self.genome.name}' (e.g. {sample_scaffolds}). "
                    f"However, {len(oriseq_overlap)} chromosome IDs match 'OriSeqID' found in FASTA header descriptions (Genome Warehouse / GWH format). "
                    f"Please use '--gwh' or '--header-id-tag OriSeqID' (or run 'aegis tidy-genome --gwh') to re-header the genome so features align."
                )
                raise ValueError(msg)
            else:
                msg = (
                    f"None of the {len(annot_chrs)} chromosome IDs in annotation '{self.name}' (e.g. {sample_annot}) "
                    f"match the {len(genome_scaffolds)} scaffold IDs in genome '{self.genome.name}' (e.g. {sample_scaffolds})."
                )
                raise ValueError(msg)
        elif len(overlap) < len(annot_chrs) and not quiet:
            missing = annot_chrs - genome_scaffolds
            if len(missing) <= 5:
                missing_str = ", ".join(sorted(missing))
            else:
                missing_str = f"{', '.join(sorted(list(missing))[:5])}... ({len(missing)} total)"
            print(f"Notice: {len(missing)} chromosomes in annotation '{self.name}' are not present in genome '{self.genome.name}': {missing_str}")

    @property
    def motifs(self) -> AnnotationMotifs:
        if not hasattr(self, '_motifs'):
            self._motifs = AnnotationMotifs(self)
        return self._motifs

    @property
    def overlaps(self) -> AnnotationOverlaps:
        if not hasattr(self, '_overlaps'):
            self._overlaps = AnnotationOverlaps(self)
        return self._overlaps

    @property
    def redundancy(self) -> AnnotationRedundancy:
        if not hasattr(self, '_redundancy'):
            self._redundancy = AnnotationRedundancy(self)
        return self._redundancy

    @property
    def stats(self) -> AnnotationStats:
        if not hasattr(self, '_stats'):
            self._stats = AnnotationStats(self)
        return self._stats

    @property
    def export(self) -> AnnotationExport:
        if not hasattr(self, '_export'):
            self._export = AnnotationExport(self)
        return self._export

    @property
    def summary(self) -> dict:
        return self.stats.data

    @property
    def genome_validation(self) -> dict:
        """
        Assembly-annotation validation findings for the input annotation, computed while loading.
        """
        return self._genome_validation

    def _validate_against_genome(self, quiet: bool = False) -> dict:
        """
        Validates the input annotation features against the associated genome assembly. Run once
        while loading, as its findings are recorded as input warnings.
        Distinguishes between:
          - Fatal mismatch: Zero contigs from annotation exist in genome assembly (leads to halt).
          - Warning discrepancies: Even a single contig missing or single feature out of bounds (leads to warning).
        Returns a dictionary summarizing findings.
        """
        results = {
            "missing_chromosomes": [],
            "out_of_bounds_features": [],
            "unannotated_scaffolds": [],
            "matching_chromosomes": [],
            "fatal_mismatch": False,
            "has_warning": False,
            "summary_message": "",
        }
        if self.genome is None or not getattr(self.genome, "scaffolds", None):
            return results

        genome_scaffolds = self.genome.scaffolds
        missing_chroms = [ch for ch in self.chrs if ch not in genome_scaffolds]
        matching_chroms = [ch for ch in self.chrs if ch in genome_scaffolds]
        results["matching_chromosomes"] = matching_chroms
        results["missing_chromosomes"] = missing_chroms

        for ch in missing_chroms:
            self.warnings["chromosome_not_in_genome"].add(ch)

        # Fatal mismatch check: absolutely zero contigs match between annotation and genome
        if len(self.chrs) > 0 and len(genome_scaffolds) > 0 and len(matching_chroms) == 0:
            results["fatal_mismatch"] = True
            results["summary_message"] = (
                f"Fatal mismatch: None of the {len(self.chrs)} contig(s) in annotation '{self.name}' "
                f"were found in genome '{self.genome.name}'. Aborting assembly comparison."
            )
            if not quiet:
                print(f"\n[ERROR] {results['summary_message']}")
            return results

        out_of_bounds = []
        for ch in matching_chroms:
            genes = self.chrs[ch]
            scf_len = genome_scaffolds[ch].size
            def exceeds(ft) -> bool:
                return (ft.start is not None and ft.start < 1) or (ft.end is not None and ft.end > scf_len)

            for g in genes.values():
                # Every feature level is checked, as in the input annotation subfeatures
                # may extend beyond the scaffold even if their parents do not
                exceeding = []
                if exceeds(g):
                    exceeding.append(g.id)
                for t in g.transcripts.values():
                    if exceeds(t):
                        exceeding.append(t.id)
                    exceeding.extend(e.id for e in t.exons if exceeds(e))
                    for c in t.CDSs.values():
                        if any(exceeds(cs) for cs in c.CDS_segments):
                            exceeding.append(c.id)
                        exceeding.extend(u.id for u in c.UTRs if exceeds(u))
                    exceeding.extend(u.id for u in (t.temp_UTRs or []) if exceeds(u))
                    exceeding.extend(cs.id for cs in (t.temp_CDSs or []) if exceeds(cs))
                if exceeding:
                    out_of_bounds.append(g.id)
                    self.warnings["feature_exceeds_scaffold_length"].update(exceeding)

        results["out_of_bounds_features"] = out_of_bounds

        # Genome scaffolds without genes
        unannotated = [scf_name for scf_name in genome_scaffolds if scf_name not in self.chrs]
        results["unannotated_scaffolds"] = unannotated

        # Any missing contig or even a single out-of-bounds feature triggers warning
        if len(out_of_bounds) > 0 or len(missing_chroms) > 0:
            results["has_warning"] = True
            msg_parts = []
            if missing_chroms:
                if len(missing_chroms) == 1:
                    msg_parts.append(f"1 contig in annotation not found in genome ('{missing_chroms[0]}')")
                else:
                    msg_parts.append(f"{len(missing_chroms)} contigs in annotation not found in genome ({', '.join(missing_chroms[:3])}{'...' if len(missing_chroms) > 3 else ''})")
            if out_of_bounds:
                if len(out_of_bounds) == 1:
                    msg_parts.append(f"1 feature has coordinates exceeding contig boundaries in genome ('{out_of_bounds[0]}')")
                else:
                    msg_parts.append(f"{len(out_of_bounds)} features have coordinates exceeding contig boundaries in genome")
            
            results["summary_message"] = (
                f"Genome-annotation discrepancies detected for '{self.name}' on '{self.genome.name}': "
                + "; ".join(msg_parts)
                + ". Annotation coordinates and assembly contig sizes should be verified."
            )
            if not quiet:
                print(f"\n[WARNING] {results['summary_message']}")

        return results

    def iter_genes(self):
        for chrom in self.chrs.values():
            yield from chrom.values()

    def iter_transcripts(self):
        for chrom in self.chrs.values():
            for g in chrom.values():
                yield from g.transcripts.values()

    def iter_exons(self):
        for chrom in self.chrs.values():
            for g in chrom.values():
                for t in g.transcripts.values():
                    yield from t.exons

    def iter_CDSs(self, segments: bool = True):
        if segments:
            for chrom in self.chrs.values():
                for g in chrom.values():
                    for t in g.transcripts.values():
                        for c in t.CDSs.values():
                            yield c
                            yield from c.CDS_segments
        else:
            for chrom in self.chrs.values():
                for g in chrom.values():
                    for t in g.transcripts.values():
                        yield from t.CDSs.values()

    def iter_UTRs(self):
        for chrom in self.chrs.values():
            for g in chrom.values():
                for t in g.transcripts.values():
                    for c in t.CDSs.values():
                        yield from c.UTRs

    def iter_introns(self):
        for chrom in self.chrs.values():
            for g in chrom.values():
                for t in g.transcripts.values():
                    if t.introns:
                        yield from t.introns

    def iter_miRNAs(self):
        for chrom in self.chrs.values():
            for g in chrom.values():
                for t in g.transcripts.values():
                    yield from t.miRNAs

    def iter_proteins(self):
        for chrom in self.chrs.values():
            for g in chrom.values():
                for t in g.transcripts.values():
                    for c in t.CDSs.values():
                        yield c.protein

    def iter_promoters(self):
        for chrom in self.chrs.values():
            for g in chrom.values():
                for t in g.transcripts.values():
                    if t.promoter is not None:
                        yield t.promoter

    def iter_features(self, include_atypical: bool = True, include_orphaned: bool = True):
        for chrom in self.chrs.values():
            for g in chrom.values():
                yield g
                for t in g.transcripts.values():
                    yield t
                    yield from t.exons
                    
                    for c in t.CDSs.values():
                        yield c
                        yield from c.CDS_segments
                        yield from c.UTRs
                    
                    if t.introns:
                        yield from t.introns
                    yield from t.miRNAs

        if include_atypical:
            yield from self.atypical_features
            
        if include_orphaned:
            yield from self.orphaned_features

    @staticmethod
    def _has_tag(text:str, tag:str) -> bool:
        # A tag is a whole part of a name: "_combined" in "X_combined_on_Y" or "X_combined.gff3", but not in "X_combinedV2"
        return re.search(rf"_{re.escape(tag)}(?=[_.]|$)", text) is not None

    def _suffix(self, tags:set[str]) -> str:
        # Tags already in the id (e.g. kept from the name of a reloaded file) are not repeated
        return "".join(f"_{tag}" for tag in sorted(tags) if not self._has_tag(self.id, tag))

    @property
    def feature_suffix(self) -> str:
        """Tags that change the features themselves, naming every export."""
        return self._suffix(self.feature_tags)

    @property
    def all_suffixes(self) -> str:
        """All tags, naming gff/gtf exports."""
        return self._suffix(self.feature_tags | self.tags)

    @property
    def chromosome_suffix(self) -> str:
        """Tags of renamed chromosomes, naming exports that list chromosome names."""
        return self._suffix(self.tags & self.chromosome_tags)

    def _resolve_output_path(self, filepath: str | None=None, output_dir: str | None=None, filename: str | None=None, 
        suffix:str = "", extra_suffixes: list[str] | None = None, 
        extension: str = ".fasta", use_annot_dir: bool = False, 
        subfolder_name: str = "features", subfolder: bool = False, use_name_not_id: bool = False, prefix:str ="", create_dir:bool=True
    ) -> Path:
        
        # Exact same logic as your original code
        if filepath:
            p = Path(filepath)
            p.parent.mkdir(parents=True, exist_ok=True)
            if (output_dir is not None) or subfolder or filename or use_annot_dir:
                warnings.warn(f"Exact output file path ({filepath}) was specified. Ignoring output_dir, use_annot_dir, subfolder, and filename arguments.")
            return p
        else:
            if use_annot_dir:
                export_folder = Path(self.path)
                if output_dir is not None:
                    warnings.warn(f"Both 'use_annot_dir={use_annot_dir}' and 'output_dir={output_dir}' were provided. Defaulting to the annotation directory ({self.path}).")
            else:
                export_folder = Path(output_dir or ".")

            if subfolder:
                export_folder = export_folder / subfolder_name.strip("/")


            if create_dir:

                export_folder.mkdir(parents=True, exist_ok=True)

            if extension and not extension.startswith("."):
                extension = f".{extension}"

            if not filename:

                if prefix and not prefix.endswith("_"):
                    prefix = f"{prefix}_"

                tag_str = "".join([f"_{t}" for t in (extra_suffixes or []) if t])
                if use_name_not_id:
                    filename = f"{prefix}{self.name}{suffix}{tag_str}{extension}"
                else:
                    filename = f"{prefix}{self.id}{suffix}{tag_str}{extension}"
            
            if not Path(filename).suffix:
                filename = f"{filename}{extension}"

            return export_folder / filename

    def _resolve_output_dir(self, output_dir: str | None = None, use_annot_dir: bool = False, subfolder_name: str = "features", subfolder: bool = False, create_dir:bool=True) -> Path:

        if use_annot_dir:
            export_folder = Path(self.path)
            if output_dir is not None:
                warnings.warn(f"Both 'use_annot_dir={use_annot_dir}' and 'output_dir={output_dir}' were provided. Defaulting to the annotation directory ({self.path}).")
        else:
            export_folder = Path(output_dir or ".")

        if subfolder:
            export_folder = export_folder / subfolder_name.strip("/")

        if create_dir:

            export_folder.mkdir(parents=True, exist_ok=True)

        return export_folder

    def load_data(self, gff_file, encoding, chosen_chromosomes:tuple[str, ...]|None=None, chosen_coordinates:tuple[int, int]|None=None, skip_atypical_features:bool=False, quiet:bool=False):
        
        staging = {"gene": [], "transcript": [], "miRNA": [],"exon": [], "CDS": [], "UTR": []}
        # start_codon/stop_codon rows are folded into the CDS rows once all subfeatures are added (see _incorporate_codons)
        self._codon_entries = []

        chromosomes_t = set()
        # gene and transcript lines by ID, to skip exact repeats (e.g. a whole gene block written twice)
        seen_lines = {}

        with open_file(gff_file, encoding=encoding) as f:
            for line_number, raw_line in enumerate(f, start=1):
                line = raw_line.strip()

                if not line:
                    continue

                # sequences appended after the features end the annotation (GFF3 spec)
                if line.startswith("##FASTA") or line[0] == ">":
                    break

                if line[0] == "#":
                    if line == "###" or line == "#":
                        continue
                    self.gff_header.append(line)
                    continue

                parts = line.split("\t")
                if len(parts) <= 2:
                    self.gff_header.append(line)
                    continue

                # strip() also took the tab before an empty attribute column, or after an empty seqid
                if len(parts) == 8:
                    raw_parts = raw_line.rstrip("\r\n").split("\t")
                    if len(raw_parts) >= 9 and not raw_parts[0].strip():
                        if not quiet:
                            print(f"{self.id} Warning: line {line_number} has no seqid and was skipped: {line}")
                        continue
                    parts.append("")

                if len(parts) < 9:
                    continue

                entry = parse_gff_parts(parts)
                entry.line = line_number

                ID = entry.id
                ft = entry.feature

                if ft == "region" and entry.attributes.get("Is_circular", "").lower() == "true":
                    self.circular_chromosomes.add(entry.ch)

                cat = "codon" if ft in default_codons else default_features_r.get(ft, "atypical")

                # Ensembl (GRCh37) writes the transcripts of pseudogenes as 'pseudogene' lines with a Parent gene
                if ft == "pseudogene" and entry.parents:
                    self._note("'pseudogene' lines with a Parent were read as pseudogenic_transcript", f"pseudogene line {entry.id} has Parent {','.join(entry.parents)} and was read as a pseudogenic_transcript", "pseudogene_line_read_as_transcript", quiet)
                    self.warnings["pseudogene_line_read_as_transcript"].add(entry.id)
                    ft = entry.feature = "pseudogenic_transcript"
                    cat = "transcript"

                if skip_atypical_features and cat == "atypical":
                    continue

                ch = entry.ch
                start = entry.start
                end = entry.end

                if chosen_chromosomes is not None:
                    if ch not in chosen_chromosomes:
                        continue

                if chosen_coordinates is not None:

                    if start < chosen_coordinates[0]:
                        continue
                    if end > chosen_coordinates[1]:
                        continue

                if not ID:
                    if not quiet:
                        print(f"{self.id} Warning: ID not found for {ft}: {line}")

                source = entry.source
                strand = entry.strand
                score = entry.score
                parents = entry.parents
                
                attributes = entry.attributes
                decreasing_coordinates = entry.decreasing_coordinates

                if decreasing_coordinates:
                    self.warnings["decreasing_coordinates"].add(ID)
                    if not quiet:
                        print(f"{self.id} Warning: Decreasing coordinates for {ft} {ID}")

                if cat == "atypical":
                    if ft not in self.features:
                        self.features[ft] = 1
                    else:
                        self.features[ft] += 1
                    chromosomes_t.add(ch)
                    self.atypical_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes))
                else:
                    if ft not in self.features:
                        self.features[ft] = 1
                    else:
                        self.features[ft] += 1
                    chromosomes_t.add(ch)
                    if cat in ("gene", "transcript") and ID:
                        line_hash = hash(line)
                        if seen_lines.get((cat, ID)) == line_hash:
                            self.warnings["repeated_feature_line"].add(ID)
                            self.features[ft] -= 1
                            self._note("exact repeats of gene or transcript lines were skipped", f"{ft} line of {ID} repeats an earlier line exactly and was skipped", "repeated_feature_line", quiet)
                            continue
                        seen_lines.setdefault((cat, ID), line_hash)
                    if cat == "codon":
                        self._codon_entries.append(entry)
                    else:
                        staging[cat].append(entry)

        chromosomes_t = sorted(list(chromosomes_t))
        for ch in chromosomes_t:
            self.chrs[ch] = {}

        return staging

    @staticmethod
    def _shared_attribute(entry, link_by:tuple[str, ...]) -> tuple[str, str] | None:
        """(attribute, value) of the first link_by attribute the entry carries as a single value, else None."""
        attrs = entry.attributes
        for key in link_by:
            value = attrs.get(key)
            if value is None:
                mapped = MAPPED_ATTRIBUTES.get(key.lower())
                if mapped is not None:
                    value = attrs.get(mapped)
            if isinstance(value, list):
                value = value[0] if len(value) == 1 else None
            if value:
                return key, value
        return None

    def _link_by_shared_attribute(self, staging:dict, link_by:tuple[str, ...], create_features:bool=True, quiet:bool=False):
        """
        Links transcripts and subfeatures that have no Parent attribute through an attribute they share
        with their gene (e.g. locus_tag in Prokka, NCBI and GenBank-derived files), before the features are
        added. Only lines without Parent are changed, so files with Parent attributes load as before.

        Features are grouped by (seqid, attribute, value). A group is linked only if it holds at most one
        gene and all its features are on one strand; otherwise nothing in it is changed and the IDs are
        listed under 'shared_attribute_not_linked'. In a linked group:
          * transcripts without Parent get the gene as Parent (a gene is created if there is none),
          * subfeatures without Parent get the group's transcript if it has exactly one; with several, the
            transcript whose ID appears in the subfeature ID, else the only one spanning it; with none,
            transcripts are created, a new one starting whenever a subfeature overlaps one of the same type
            already taken (an isoform repeating an exon), in file order.
        Features are created only if create_features (infer_genes_from_subfeatures) is True. Every linked
        feature is listed under 'linked_by_shared_attribute'.
        """
        # GFF3: rows sharing an ID are one feature, so a row without Parent takes the Parent of its
        # siblings (e.g. one CDS row of a multi-row CDS whose Parent was mistyped)
        for stage in ("miRNA", "exon", "CDS", "UTR"):
            orphan_ids = {e.id for e in staging[stage] if not e.parents and e.id}
            if not orphan_ids:
                continue
            sibling_parents = {}
            for e in staging[stage]:
                if e.parents and e.id in orphan_ids:
                    sibling_parents.setdefault((e.ch, e.feature, e.id), e.parents)
            n = 0
            for e in staging[stage]:
                if not e.parents and (e.ch, e.feature, e.id) in sibling_parents:
                    e.parents = list(sibling_parents[(e.ch, e.feature, e.id)])
                    self.warnings["linked_by_shared_attribute"].add(e.id)
                    n += 1
            if n and not quiet:
                print(f"{self.id} Warning: {n} {stage} rows without Parent took the Parent of rows sharing their ID" + ("" if self.list_load_warnings else self._warning_ids_hint("linked_by_shared_attribute")))

        if not link_by:
            return

        subfeature_entries = [e for stage in ("miRNA", "exon", "CDS", "UTR") for e in staging[stage] if not e.parents]
        subfeature_entries += [e for e in self._codon_entries if not e.parents]
        unparented_transcripts = [e for e in staging["transcript"] if not e.parents]
        if not subfeature_entries and not unparented_transcripts:
            return

        groups = {}
        def group(entry):
            shared = self._shared_attribute(entry, link_by)
            if shared is None:
                return None
            key = (entry.ch,) + shared
            if key not in groups:
                groups[key] = {"genes": [], "transcripts": [], "subfeatures": []}
            return groups[key]

        for e in unparented_transcripts:
            g = group(e)
            if g is not None:
                g["transcripts"].append(e)
        for e in subfeature_entries:
            g = group(e)
            if g is not None:
                g["subfeatures"].append(e)
        if not groups:
            return

        # genes and transcripts with Parent only matter for the groups found above
        for e in staging["gene"]:
            shared = self._shared_attribute(e, link_by)
            if shared is not None and (e.ch,) + shared in groups:
                groups[(e.ch,) + shared]["genes"].append(e)

        gene_ids = {e.id for e in staging["gene"] if e.id}
        transcript_ids = {e.id for e in staging["transcript"] if e.id}
        transcripts_by_gene = {}
        for e in staging["transcript"]:
            for p in e.parents:
                transcripts_by_gene.setdefault(p, []).append(e)
            if e.parents:
                shared = self._shared_attribute(e, link_by)
                if shared is not None and (e.ch,) + shared in groups:
                    groups[(e.ch,) + shared]["transcripts"].append(e)

        def unique_id(base, taken):
            new_id, n = base, 1
            while new_id in taken:
                n += 1
                new_id = f"{base}_{n}"
            taken.add(new_id)
            return new_id

        def new_entry(template, ID, feature, members, parents):
            return GffEntry(ch=template.ch, source=template.source, feature=feature,
                            start=min(m.start for m in members), end=max(m.end for m in members),
                            score=".", strand=template.strand, attributes={}, id=ID, parents=parents,
                            pseudogene=False, transposable=False, decreasing_coordinates=False,
                            line=0) # created, not read: no place in the file order

        linked, not_linked = [], []

        def link(entry, parent_id):
            entry.parents = [parent_id]
            if not entry.id:
                # as parse_gff_parts names lines that have a Parent but no ID
                entry.id = f"{parent_id}_{entry.feature}_{entry.start}_{entry.end}"
            linked.append(entry.id)

        for (ch, key, value), g in groups.items():
            members = g["genes"] + g["transcripts"] + g["subfeatures"]
            strands = {m.strand for m in members if m.strand in ("+", "-")}
            if len(g["genes"]) > 1 or len(strands) > 1:
                reason = f"{len(g['genes'])} genes share it" if len(g["genes"]) > 1 else "its features are on both strands"
                not_linked.extend(m.id or f"{m.feature}_{m.ch}_{m.start}_{m.end}" for m in g["transcripts"] + g["subfeatures"] if not m.parents)
                if not quiet and self.list_load_warnings:
                    print(f"{self.id} Warning: features with {key}={value} on {ch} were not linked through it, as {reason}")
                continue

            gene = g["genes"][0] if g["genes"] else None
            if gene is not None and not gene.id:
                gene.id = unique_id(value, gene_ids)

            def ensure_gene(template):
                nonlocal gene
                if gene is None and create_features:
                    gene = new_entry(template, unique_id(value, gene_ids), "gene", members, [])
                    gene.attributes = {key: value}
                    staging["gene"].append(gene)
                return gene

            transcripts = list(g["transcripts"])
            if gene is not None:
                transcripts += [t for t in transcripts_by_gene.get(gene.id, []) if t not in transcripts]

            for t in g["transcripts"]:
                if t.parents:
                    continue
                if ensure_gene(t) is not None:
                    link(t, gene.id)

            subfeatures = sorted(g["subfeatures"], key=lambda e: e.line)
            if not subfeatures:
                continue

            if len(transcripts) == 1:
                for sf in subfeatures:
                    link(sf, transcripts[0].id)

            elif transcripts:
                for sf in subfeatures:
                    named = [t for t in transcripts if t.id and re.search(rf"(?<![A-Za-z0-9]){re.escape(t.id)}(?![A-Za-z0-9])", sf.id)]
                    spanning = [t for t in transcripts if t.start <= sf.start and sf.end <= t.end]
                    chosen = named if len(named) == 1 else spanning
                    if len(chosen) == 1:
                        link(sf, chosen[0].id)
                    else:
                        not_linked.append(sf.id or f"{sf.feature}_{sf.ch}_{sf.start}_{sf.end}")

            elif ensure_gene(subfeatures[0]) is not None:
                isoforms = [[]]
                for sf in subfeatures:
                    if any(o.feature == sf.feature and o.start <= sf.end and sf.start <= o.end for o in isoforms[-1]):
                        isoforms.append([])
                    isoforms[-1].append(sf)
                for isoform in isoforms:
                    coding = any(sf.feature in default_features["CDS"] for sf in isoform)
                    t = new_entry(isoform[0], unique_id(f"{gene.id}_t1", transcript_ids), "mRNA" if coding else "transcript", isoform, [gene.id])
                    t.attributes = {key: value}
                    staging["transcript"].append(t)
                    for sf in isoform:
                        link(sf, t.id)

            else:
                not_linked.extend(sf.id or f"{sf.feature}_{sf.ch}_{sf.start}_{sf.end}" for sf in subfeatures)

        for ID in linked:
            self.warnings["linked_by_shared_attribute"].add(ID)
        for ID in not_linked:
            self.warnings["shared_attribute_not_linked"].add(ID)

        if not quiet:
            if linked:
                print(f"{self.id} Warning: {len(linked)} features without Parent were linked to their gene or transcript by a shared attribute ({', '.join(link_by)})" + ("" if self.list_load_warnings else self._warning_ids_hint("linked_by_shared_attribute")))
            if not_linked:
                print(f"{self.id} Warning: {len(not_linked)} features without Parent share a {'/'.join(link_by)} value but could not be linked unambiguously (several genes share it, mixed strands or several candidate transcripts)" + ("" if self.list_load_warnings else self._warning_ids_hint("shared_attribute_not_linked")))


    def _link_by_file_order(self, staging:dict, apply:bool=False, quiet:bool=False):
        """
        Last resort for lines that are still without Parent after Parent attributes, shared IDs and
        shared attributes (link_by): each joins the transcript, or else the gene, that it follows in the
        file. Only with apply (link_by_file_order=True); otherwise the lines it would link are listed
        under 'linkable_by_file_order' and a hint to opt in is printed.

        A line joins only if it lies within that transcript/gene on the same seqid and strand. A
        subfeature that overlaps one of the same type already in the transcript (an isoform, or a
        repeat) ends the chain: it and the lines after it stay unlinked until the next gene or
        transcript line. Lines whose shared attribute was ambiguous ('shared_attribute_not_linked') and
        features created while linking are never used. Lines not linked fall back to the usual handling
        of features without Parent.
        """
        subfeatures = [e for stage in ("miRNA", "exon", "CDS", "UTR") for e in staging[stage]] + self._codon_entries
        if not any(not e.parents for e in subfeatures) and not any(not e.parents for e in staging["transcript"]):
            return

        def name(e):
            return e.id or f"{e.feature}_{e.ch}_{e.start}_{e.end}"

        ambiguous = self.warnings["shared_attribute_not_linked"]
        genes_by_id = {e.id: e for e in staging["gene"]}
        transcripts_by_id = {e.id: e for e in staging["transcript"]}

        lines = [(e.line, 0, e) for e in staging["gene"] if e.line]
        lines += [(e.line, 1, e) for e in staging["transcript"] if e.line]
        lines += [(e.line, 2, e) for e in subfeatures if e.line]
        lines.sort(key=lambda x: x[0])

        def fits(e, container):
            if container is None or e.ch != container.ch:
                return False
            if e.strand in ("+", "-") and container.strand in ("+", "-") and e.strand != container.strand:
                return False
            return container.start <= e.start and e.end <= container.end

        candidates = []
        gene = transcript = None
        taken = []  # subfeatures already in the current transcript (or gene), to spot isoforms
        for _, kind, e in lines:
            if kind == 0:
                gene, transcript, taken = e, None, []
            elif kind == 1:
                if e.parents:
                    gene = genes_by_id.get(e.parents[0])
                elif name(e) not in ambiguous and fits(e, gene):
                    candidates.append((e, gene))
                else:
                    gene = None
                transcript, taken = e, []
            else:
                if e.parents:
                    parent = transcripts_by_id.get(e.parents[0])
                    if parent is not None and parent is not transcript:
                        transcript, taken = parent, []
                        gene = genes_by_id.get(parent.parents[0]) if parent.parents else None
                    taken.append(e)
                    continue
                if name(e) in ambiguous:
                    continue
                container = transcript if transcript is not None else gene
                if not fits(e, container):
                    continue
                if any(o.feature == e.feature and o.start <= e.end and e.start <= o.end for o in taken):
                    # an isoform or a repeat: file order cannot tell which transcript it belongs to
                    gene = transcript = None
                    taken = []
                    continue
                candidates.append((e, container))
                taken.append(e)

        if not candidates:
            return

        if apply:
            for e, container in candidates:
                e.parents = [container.id]
                if not e.id:
                    e.id = f"{container.id}_{e.feature}_{e.start}_{e.end}"
                self.warnings["linked_by_file_order"].add(e.id)
            if not quiet:
                print(f"{self.id} Warning: {len(candidates)} features without Parent were linked to the gene or transcript they follow in the file (link_by_file_order)" + ("" if self.list_load_warnings else self._warning_ids_hint("linked_by_file_order")))
        else:
            for e, _ in candidates:
                self.warnings["linkable_by_file_order"].add(name(e))
            if not quiet:
                print(f"{self.id} Warning: {len(candidates)} features without Parent lie within the gene or transcript they follow in the file, but were not linked to it; if the file lists each gene's features after it, try loading it with --link-by-file-order (link_by_file_order=True). Without it they are kept as features without Parent" + ("" if self.list_load_warnings else self._warning_ids_hint("linkable_by_file_order")))


    def _add_gene(self, entry, rename_repeated_id:bool=False, skip_orphaned_features:bool=False, quiet:bool=False):
        ID = entry.id
        transposable = entry.transposable
        pseudogene = entry.pseudogene
        ch = entry.ch
        source = entry.source
        ft = entry.feature
        strand = entry.strand
        start = entry.start
        end = entry.end
        score = entry.score
        attributes = entry.attributes

        if start == end:
            self.errors["1bp_gene"].add(ID)
            if not quiet:
                print(f"{self.id} Error: 1bp gene {ID} feature")

        if ID not in self.all_gene_ids:
            self.all_gene_ids[ID] = ch
            self.chrs[ch][ID] = Gene(pseudogene, transposable, ID, ch, source, ft, strand, start, end, score, [], attributes)
        else:
            self.errors["repeated_gene_IDs"].add(ID)
            if not rename_repeated_id:
                if not quiet:
                    print(f"{self.id} Error: repeated gene ID {ID} will not be added as a gene based on the chosen parameter (rename_repeated_id={rename_repeated_id})")
                if not skip_orphaned_features:
                    self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, [], attributes))
            else:
                renamed_id = self._get_unique_gene_id(ID)
                if not quiet:
                    print(f"{self.id} Warning: repeated gene id {ID} was renamed to {renamed_id}")
                self.all_gene_ids[renamed_id] = ch
                self.chrs[ch][renamed_id] = Gene(pseudogene, transposable, renamed_id, ch, source, ft, strand, start, end, score, [], attributes)
                if ID not in self._gene_info:
                    self._gene_info[ID] = {ID, renamed_id}
                else:
                    self._gene_info[ID].add(renamed_id)

    def _add_transcript(self, entry, rename_repeated_id:bool=False, infer_gene_from_transcript:bool=False, skip_orphaned_features:bool=False, quiet:bool=False):
        ID = entry.id
        ch = entry.ch
        source = entry.source
        ft = entry.feature
        strand = entry.strand
        start = entry.start
        end = entry.end
        score = entry.score
        attributes = entry.attributes
        parents = entry.parents

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
                self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes))

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
                    entry.id = created_gene_id
                    entry.feature = "gene"
                    self._add_gene(entry, quiet=quiet)

                    self.all_transcript_ids[unique_t_id] = (ch, created_gene_id)
                    self.chrs[ch][created_gene_id].transcripts[unique_t_id] = Transcript(unique_t_id, ch, source, ft, strand, start, end, score, [created_gene_id], attributes)
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
            # This is an output in some lifton/liftoff cases and they are repeated entries, so ignoring them is
            # appropriate. Other GFFs give the mRNA the id of its gene: those are kept once all transcripts are loaded.
            if parent == ID:
                if self._defer_self_parented_transcripts and ID in self.all_gene_ids:
                    self._self_parented_transcripts.append(entry)
                else:
                    self.warnings["self_parented_transcript"].add(ID)
                    self._note("transcripts had their own id as a parent so they were ignored", f"Transcript {ID} has its own id as a parent so this transcript feature was ignored", "self_parented_transcript", quiet)
                    if not skip_orphaned_features:
                        self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes))

            elif parent not in self.all_gene_ids:

                self.warnings["transcript_to_inexistent_gene"].add(ID)

                if infer_gene_from_transcript:

                    if unique_t_id not in self.all_transcript_ids:
                        entry.id = parent
                        entry.feature = "gene"
                        self._add_gene(entry, quiet=quiet)

                        self.all_transcript_ids[unique_t_id] = (ch, parent)
                        self.chrs[ch][parent].transcripts[unique_t_id] = Transcript(unique_t_id, ch, source, ft, strand, start, end, score, parents, attributes)
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

            elif parent in self._gene_info:

                possible_parents = self._gene_info[parent]

                found = False

                for parent in possible_parents:

                    if parent not in self.chrs[ch]:
                        continue

                    if self.chrs[ch][parent].end < start:
                        continue
                    if self.chrs[ch][parent].start > end:
                        continue

                    found = True
                    if unique_t_id not in self.all_transcript_ids:
                        self.all_transcript_ids[unique_t_id] = (ch, parent)
                        self.chrs[ch][parent].transcripts[unique_t_id] = Transcript(unique_t_id, ch, source, ft, strand, start, end, score, [parent], attributes)
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
                        self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes))

            elif self.all_gene_ids[parent] != ch:
                self.errors["transcript_to_gene_other_chr"].add(ID)
                if not quiet:
                    print(f"{self.id} Error: {ID} transcript refers to a gene in a different chromosome, it could not be assigned to a gene")
                if not skip_orphaned_features:
                    self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes))

            else:

                if unique_t_id not in self.all_transcript_ids:
                    self.all_transcript_ids[unique_t_id] = (ch, parent)
                    self.chrs[ch][parent].transcripts[unique_t_id] = Transcript(unique_t_id, ch, source, ft, strand, start, end, score, parents, attributes)
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

    def _note(self, summary:str, message:str, warning_key:str, quiet:bool=False):
        """
        Load-time note: counted under summary, which is printed once per annotation after loading
        with a pointer to the warning_key column of the aegis summary warnings export, and printed
        individually only with list_load_warnings.
        """
        if hasattr(self, "_load_note_counts"):
            key = (summary, warning_key)
            self._load_note_counts[key] = self._load_note_counts.get(key, 0) + 1
            if self.list_load_warnings and not quiet:
                print(f"{self.id} Warning: {message}")
        elif not quiet:
            print(f"{self.id} Warning: {message}")

    def _incorporate_codons(self, skip_orphaned_features:bool=False, quiet:bool=False):
        """
        Folds the start_codon/stop_codon rows read by load_data into the CDS rows of their transcripts
        (see Transcript.incorporate_codons), so CDSs annotated without their stop (or start) codon, as in
        GTF 2.2, cover it. Codon rows within the CDS rows change nothing. Run once all subfeatures are added,
        while repeated transcript ids can still be resolved.
        """
        codons_by_transcript = {}
        for entry in self._codon_entries:
            transcripts = []
            for parent in entry.parents:
                if parent in self._transcript_info:
                    candidates = self._transcript_info[parent]
                elif parent in self.all_transcript_ids:
                    candidates = (parent,)
                elif parent in self.all_gene_ids and parent in self.chrs.get(entry.ch, {}) and len(self.chrs[entry.ch][parent].transcripts) == 1:
                    candidates = tuple(self.chrs[entry.ch][parent].transcripts)
                else:
                    candidates = ()
                for t_id in candidates:
                    ch, g_id = self.all_transcript_ids[t_id]
                    if ch != entry.ch:
                        continue
                    t = self.chrs[ch][g_id].transcripts[t_id]
                    if len(candidates) > 1 and (t.end < entry.start or t.start > entry.end):
                        continue
                    transcripts.append(t)

            if not transcripts:
                self.warnings["missing_subfeature_parent"].add(entry.id)
                self._note("start/stop codon rows could not be assigned to a transcript", f"{entry.feature} {entry.id} could not be assigned to a transcript (parents: {entry.parents})", "missing_subfeature_parent", quiet)
                if not skip_orphaned_features:
                    self.orphaned_features.append(Feature(entry.id, entry.ch, entry.source, entry.feature, entry.strand, entry.start, entry.end, entry.score, entry.parents, entry.attributes))
                continue

            for t in transcripts:
                if t.id not in codons_by_transcript:
                    codons_by_transcript[t.id] = (t, [])
                codons_by_transcript[t.id][1].append((entry.start, entry.end))

        for t, codons in codons_by_transcript.values():
            changed, not_placed = t.incorporate_codons(codons)
            if changed:
                self.warnings["codon_added_to_CDS"].add(t.id)
                self._note("transcripts had start/stop codon rows outside their CDS rows, which now cover them", f"Transcript {t.id} CDS rows were extended to cover its start/stop codon rows", "codon_added_to_CDS", quiet)
            if not_placed:
                self.warnings["codon_not_incorporated"].add(t.id)
                if not t.temp_CDSs:
                    message = f"Transcript {t.id} has start/stop codon rows but no CDS rows, so they were ignored"
                else:
                    message = f"Transcript {t.id} start/stop codon rows {not_placed} neither lie within nor continue its CDS rows, so they were ignored"
                self._note("transcripts had start/stop codon rows that could not be added to their CDS", message, "codon_not_incorporated", quiet)

        del self._codon_entries

    def _add_self_parented_transcripts(self, skip_orphaned_features:bool=False, quiet:bool=False):
        """
        Transcripts whose parent is their own id (e.g. ID=Pda01g0001;Parent=Pda01g0001) are kept, renamed
        {gene id}_t1, as the transcript of the gene with that id if it is not a pseudogene and has no other
        transcripts. Otherwise they are repeated entries (as in some lifton/liftoff outputs) and are ignored.
        """
        self._defer_self_parented_transcripts = False
        for entry in self._self_parented_transcripts:
            ID = entry.id
            ch = entry.ch
            self.warnings["self_parented_transcript"].add(ID)
            gene_ch = self.all_gene_ids.get(ID)
            # pseudogenes get their pseudotranscript from their subfeatures instead
            if gene_ch == ch and not self.chrs[ch][ID].transcripts and not self.chrs[ch][ID].pseudogene:
                count = 1
                t_id = f"{ID}_t{count}"
                while t_id in self.all_transcript_ids:
                    count += 1
                    t_id = f"{ID}_t{count}"
                self.all_transcript_ids[t_id] = (ch, ID)
                self.chrs[ch][ID].transcripts[t_id] = Transcript(t_id, ch, entry.source, entry.feature, entry.strand, entry.start, entry.end, entry.score, [ID], entry.attributes)
                self._note("transcripts had their gene's id as their own id and parent, so they were kept renamed as {gene id}_t1", f"Transcript {ID} has its own id as a parent, kept as transcript {t_id} of gene {ID}", "self_parented_transcript", quiet)
            else:
                self._note("transcripts had their own id as a parent so they were ignored", f"Transcript {ID} has its own id as a parent so this transcript feature was ignored", "self_parented_transcript", quiet)
                if not skip_orphaned_features:
                    self.orphaned_features.append(Feature(ID, ch, entry.source, entry.feature, entry.strand, entry.start, entry.end, entry.score, entry.parents, entry.attributes))
        self._self_parented_transcripts = []

    def _add_subfeature(self, entry, ft_level, infer_gene_and_transcript_from_subfeatures:bool=False, only_infer_if_none_of_the_parents_exist:bool=True, skip_orphaned_features:bool=False, liftover_exception:bool=True, quiet:bool=False):
        ID = entry.id
        ch = entry.ch
        source = entry.source
        ft = entry.feature
        strand = entry.strand
        start = entry.start
        end = entry.end
        score = entry.score
        attributes = entry.attributes
        parents = entry.parents
        phase = entry.phase if (ft_level == "CDS" or ft in default_features["CDS"]) else None

        if start == end:
            self.warnings[f"1bp_{ft_level}"].add(ID)

        if not parents:
            self.warnings["subfeature_with_no_parent"].add(ID)

            if ft_level not in ("CDS", "exon"):

                if not quiet:
                    print(f"{self.id} Warning: No parent provided so the following {ft} subfeature {ID} could not be assigned to any transcript")

                if not skip_orphaned_features:
                    self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes, phase=phase))
            
            elif infer_gene_and_transcript_from_subfeatures:

                # Lines with neither ID nor Parent are named after their location, so that they do not all
                # end up in one '_gene' (CDS lines sharing an ID still form one transcript)
                base_id = ID or f"{ft}_{ch}_{start}_{end}"
                if not ID:
                    ID = base_id
                created_gene = f"{base_id}_gene"
                created_transcript = f"{base_id}_transcript"

                if created_gene not in self.all_gene_ids:
                    if not quiet and ft != "nucleotide_to_protein_match":
                        print(f"{self.id} Warning: {ft} {ID} subfeature was assigned to created transcript {created_transcript} and created gene {created_gene}")
                    entry.id = created_gene
                    entry.feature = "gene"
                    self._add_gene(entry, quiet=quiet)

                if created_transcript not in self.all_transcript_ids:
                    entry.id = created_transcript
                    entry.feature = "transcript"
                    entry.parents = [created_gene]
                    self._add_transcript(entry, quiet=quiet, skip_orphaned_features=skip_orphaned_features)

                if ft_level == "CDS":
                    self.chrs[ch][created_gene].transcripts[created_transcript].temp_CDSs.append(Feature(ID, ch, source, ft, strand, start, end, score, [created_transcript], attributes, phase=phase))
                else:
                    self.chrs[ch][created_gene].transcripts[created_transcript].exons.append(Exon(ID, ch, source, ft, strand, start, end, score, [created_transcript], attributes))

            else:

                if not quiet:
                    print(f"{self.id} Warning: No parent provided and infer_gene_and_transcript_from_subfeatures={infer_gene_and_transcript_from_subfeatures} so the following {ft} subfeature {ID} could not be assigned to a created transcript")

                if not skip_orphaned_features:
                    self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes, phase=phase))

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
                            self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes, phase=phase))
                        continue

                    gene_parent = self.all_transcript_ids[parent][1]

                    if ft_level == "CDS":
                        self.chrs[ch][gene_parent].transcripts[parent].temp_CDSs.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes, phase=phase))
                    elif ft_level == "exon":
                        self.chrs[ch][gene_parent].transcripts[parent].exons.append(Exon(ID, ch, source, ft, strand, start, end, score, parents, attributes))
                    elif ft_level == "UTR":
                        self.chrs[ch][gene_parent].transcripts[parent].temp_UTRs.append(UTR(ID, ch, source, ft, strand, start, end, score, parents, attributes))
                    else:
                        self.chrs[ch][gene_parent].transcripts[parent].miRNAs.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes))
                        self._miRNA_info.add(ID)

                elif parent in self.all_gene_ids:
                    if self.all_gene_ids[parent] != ch:
                        self.errors["subfeature_to_gene_other_chr"].add(ID)
                        if not quiet:
                            print(f"{self.id} Error: {ID} subfeature refers to a gene in a different chromosome, it could not be assigned to its parent")
                        if not skip_orphaned_features:
                            self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes, phase=phase))
                        continue

                    # cases where subfeature directly references a pseudogene creating a pseudotranscript for the pseudogene
                    if self.chrs[ch][parent].pseudogene:
                        if parent[0:4] == "gene":
                            pseudo_t = f"pseudo_t_{parent[5:]}"
                        else:
                            pseudo_t = f"pseudo_t_{parent}"

                        if self.chrs[ch][parent].transcripts == {}:
                            # PSEUDOGENES referred to by transcript subfeatures are given a single pseudotranscript
                            self.chrs[ch][parent].transcripts[pseudo_t] = Transcript(pseudo_t, self.chrs[ch][parent].ch, self.chrs[ch][parent].source, "pseudotranscript", self.chrs[ch][parent].strand, self.chrs[ch][parent].start, self.chrs[ch][parent].end, self.chrs[ch][parent].score, [parent])
                            self.all_transcript_ids[pseudo_t] = (ch, parent)
                            
                        if pseudo_t in self.chrs[ch][parent].transcripts:
                            if ft_level == "CDS":
                                self.chrs[ch][parent].transcripts[pseudo_t].temp_CDSs.append(Feature(ID, ch, source, ft, strand, start, end, score, [pseudo_t], attributes, phase=phase))
                            elif ft_level == "exon":
                                self.chrs[ch][parent].transcripts[pseudo_t].exons.append(Exon(ID, ch, source, ft, strand, start, end, score, [pseudo_t], attributes))
                            elif ft_level == "UTR":
                                self.chrs[ch][parent].transcripts[pseudo_t].temp_UTRs.append(UTR(ID, ch, source, ft, strand, start, end, score, [pseudo_t], attributes))
                            else:
                                self.chrs[ch][parent].transcripts[pseudo_t].miRNAs.append(Feature(ID, ch, source, ft, strand, start, end, score, [pseudo_t], attributes))
                                self._miRNA_info.add(ID)
                        else:
                            if not skip_orphaned_features:
                                self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes, phase=phase))
                            if not quiet:
                                print(f"{self.id} Error: {parent} pseudogene already had a transcript which {ft} subfeature {ID} ignores")

                    # cases where subfeature directly points to a gene that is not a pseudogene
                    else:
                        # gene without transcripts pointed to then auto-create a transcript
                        if self.chrs[ch][parent].transcripts == {}:
                            if not infer_gene_and_transcript_from_subfeatures:
                                self.warnings["subfeature_to_gene"].add(ID)
                                if not skip_orphaned_features:
                                    self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes, phase=phase))
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
                                    gene_obj.score, [parent])
                                self._note("genes with subfeatures but no transcripts got an auto-created transcript", f"{ft} subfeature {ID} references {parent} gene which has no transcripts, auto-created transcript {auto_t_id}", "subfeature_to_gene", quiet)
                                self.warnings["subfeature_to_gene"].add(ID)
                        # correctly linking the subfeature to the single transcript that exists
                        if len(self.chrs[ch][parent].transcripts) == 1:
                            temp_id = list(self.chrs[ch][parent].transcripts.keys())[0]
                            if ft_level == "CDS":
                                self.chrs[ch][parent].transcripts[temp_id].temp_CDSs.append(Feature(ID, ch, source, ft, strand, start, end, score, [temp_id], attributes, phase=phase))
                            elif ft_level == "exon":
                                self.chrs[ch][parent].transcripts[temp_id].exons.append(Exon(ID, ch, source, ft, strand, start, end, score, [temp_id], attributes))
                            elif ft_level == "UTR":
                                self.chrs[ch][parent].transcripts[temp_id].temp_UTRs.append(UTR(ID, ch, source, ft, strand, start, end, score, [temp_id], attributes))
                            else:
                                self.chrs[ch][parent].transcripts[temp_id].miRNAs.append(Feature(ID, ch, source, ft, strand, start, end, score, [temp_id], attributes))
                                self._miRNA_info.add(ID)  
                        elif len(self.chrs[ch][parent].transcripts) > 1:
                            if not quiet:
                                print(f"{self.id} Error: {ft} subfeature {ID} references {parent} gene which is not a pseudogene and has multiple transcripts")
                            self.errors["subfeature_to_gene"].add(ID)

                # if parent an miRNA
                elif parent in self._miRNA_info:
                    if ft_level == "exon":
                        if not skip_orphaned_features:
                            self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes))
                    elif not quiet:
                        print(f"Warning: {ft} subfeature {ID} references {parent} miRNA and this feature is not an exon.")

                # parent not found within created genes or transcripts or miRNAs
                else:
                    
                    if liftover_exception and self.original_annotation is not None:
                        if not quiet:
                            print(f"{self.id} Warning: {ft} subfeature {ID} references {parent} which is not found in the annotation, possibly not transfered during liftover")
                        self.warnings["missing_subfeature_parent_liftover"].add(ID)
                        if not skip_orphaned_features:
                            self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes))
                        continue

                    else:
                        if not infer_gene_and_transcript_from_subfeatures:
                            self.warnings["missing_subfeature_parent"].add(ID)
                            if not skip_orphaned_features:
                                self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes))
                            if not quiet:
                                print(f"{self.id} Warning: {ft} subfeature {ID} references {parent} which is not found in the annotation")

                        elif only_infer_if_none_of_the_parents_exist and not true_orphans:
                            self.warnings["missing_subfeature_parent"].add(ID)
                            if not skip_orphaned_features:
                                self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes))
                            if not quiet:
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

                            entry.id = inferred_g_id
                            entry.feature = "gene"
                            self._add_gene(entry, quiet=quiet)

                            entry.id = parent
                            entry.feature = "transcript"
                            entry.parents = [inferred_g_id]
                            self._add_transcript(entry, quiet=quiet, skip_orphaned_features=skip_orphaned_features)

                            if true_orphans:
                                temp_parents = parents
                            else:
                                temp_parents = [parent]


                            if ft_level == "CDS":
                                self.chrs[ch][inferred_g_id].transcripts[parent].temp_CDSs.append(Feature(ID, ch, source, ft, strand, start, end, score, temp_parents, attributes, phase=phase))
                            elif ft_level == "exon":
                                self.chrs[ch][inferred_g_id].transcripts[parent].exons.append(Exon(ID, ch, source, ft, strand, start, end, score, temp_parents, attributes))
                            elif ft_level == "UTR":
                                self.chrs[ch][inferred_g_id].transcripts[parent].temp_UTRs.append(UTR(ID, ch, source, ft, strand, start, end, score, temp_parents, attributes))
                            else:
                                self.chrs[ch][inferred_g_id].transcripts[parent].miRNAs.append(Feature(ID, ch, source, ft, strand, start, end, score, temp_parents, attributes))
                                self._miRNA_info.add(ID)

            # Cases where the transcript parent was repeated by id in the gff
            else:
                
                possible_parents = self._transcript_info[parent]

                found = False

                for parent in possible_parents:

                    gene_parent = self.all_transcript_ids[parent][1]

                    if gene_parent not in self.chrs[ch]:
                        continue
                    if self.chrs[ch][gene_parent].transcripts[parent].end < start:
                        continue
                    if self.chrs[ch][gene_parent].transcripts[parent].start > end:
                        continue

                    found = True

                    if ft_level == "CDS":
                        self.chrs[ch][gene_parent].transcripts[parent].temp_CDSs.append(Feature(ID, ch, source, ft, strand, start, end, score, [parent], attributes, phase=phase))
                    elif ft_level == "exon":
                        self.chrs[ch][gene_parent].transcripts[parent].exons.append(Exon(ID, ch, source, ft, strand, start, end, score, [parent], attributes))
                    elif ft_level == "UTR":
                        self.chrs[ch][gene_parent].transcripts[parent].temp_UTRs.append(UTR(ID, ch, source, ft, strand, start, end, score, [parent], attributes))
                    else:
                        self.chrs[ch][gene_parent].transcripts[parent].miRNAs.append(Feature(ID, ch, source, ft, strand, start, end, score, [parent], attributes))

                if not found and not quiet:
                    print(f"{self.id} Error: {ID} {ft} feature could not be assigned to any transcript. Possibly due to unforseen id clash issue")
                    if not skip_orphaned_features:
                        self.orphaned_features.append(Feature(ID, ch, source, ft, strand, start, end, score, parents, attributes, phase=phase))

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
    
    def update(self, rename_features:tuple[str,...]=(), keep_existing_ids_if_derived_from_base_id:bool=False, define_synteny:bool=False, sort_processes:int=1, quiet:bool=False, collapse_exons:bool|None=None, collapse_CDSs:bool|None=None, standardise_features:bool=False, remove_missing_transcript_parent_references:bool=False, remove_transcripts_with_no_exons:bool=False, remove_genes_with_no_transcripts:bool=False, remove_genes_with_no_transcripts_even_if_pseudogene:bool=False, update_gene_and_transcript_list:bool=False, coordinate_polishing:Literal["none", "contain", "trim"]|None=None, recalculate_phases:bool|None=None, reset_phases_zero:bool|None=None):
        start_time = time.time()

        self._protein_qc_summary = None
        if collapse_exons is None:
            collapse_exons = getattr(self, "collapse_exons", True)
        if collapse_CDSs is None:
            collapse_CDSs = getattr(self, "collapse_CDSs", True)
        if recalculate_phases is None:
            recalculate_phases = getattr(self, "recalculate_phases", False)
        if reset_phases_zero is None:
            reset_phases_zero = getattr(self, "reset_phases_zero", False)
        # Checks only record into the input warnings during the load-time update
        warnings_sink = self.warnings if getattr(self, "_recording_input_warnings", False) else self._new_warnings_dict()

        self.detect_polycistronic_transcripts(quiet=quiet, record_to=warnings_sink)

        batch_size = 1000
        count = 0

        progress_bar = start_progress_bar(total=len(self.all_gene_ids), description=f"Updating {self.id} genes", colour="62", quiet=quiet)

        for genes in self.chrs.values():
            for g in genes.values():
                count += 1
                if count >= batch_size:
                    progress_bar.update(count)
                    count = 0
                # UTRs read from the GFF are only staged until the first update, which decides whether they are used
                read_UTRs = {}
                for t in g.transcripts.values():
                    # Checked before collapsing, so the input exons are assessed
                    self._check_exon_overlaps(t, warnings_sink)
                    if t.temp_UTRs:
                        read_UTRs[t.id] = (t.read_UTRs_CDS is not None, {(u.start, u.end) for u in t.temp_UTRs})
                g.collapse_subfeatures(exons=collapse_exons, CDSs=collapse_CDSs)
                for t in g.transcripts.values():
                    t.update(quiet=quiet)
                    if t.id in read_UTRs:
                        self._check_read_UTRs(t, *read_UTRs[t.id], warnings_sink)

                    for c in t.CDSs.values():
                        # Checked before phases are recalculated, so the input phases are assessed
                        self._check_phase_mismatch_across_intron(c, warnings_sink)
                        if reset_phases_zero:
                            c.update_phase(override=True, full_override=True)
                            c.update_frame()
                        elif recalculate_phases:
                            c.update_phase(override=True, full_override=False)
                            c.update_frame()
                g.update(quiet=quiet)

        if count > 0:
            progress_bar.update(count)
        
        progress_bar.close()
        if standardise_features:
            self.update_features(standardise=standardise_features, quiet=quiet)
        
        # Removals happen before renaming so that numbering has no gaps
        self.detect_genes_with_no_transcripts(remove=remove_genes_with_no_transcripts, remove_pseudogene=remove_genes_with_no_transcripts_even_if_pseudogene, quiet=quiet, record_to=warnings_sink)
        self.detect_transcripts_with_no_exons(remove_transcripts=remove_transcripts_with_no_exons, remove_genes_accordingly=remove_genes_with_no_transcripts, quiet=quiet, record_to=warnings_sink)

        if coordinate_polishing is None:
            coordinate_polishing = getattr(self, "coordinate_polishing", "contain")

        # Coordinate inconsistencies are always checked (and recorded during loading) on the
        # unpolished coordinates, so warnings do not depend on the polishing mode
        self.correct_gene_transcript_and_subfeature_coordinates(mode="none", quiet=quiet, record_to=warnings_sink)
        if coordinate_polishing != "none":
            self.correct_gene_transcript_and_subfeature_coordinates(mode=coordinate_polishing, quiet=quiet)

        # Coordinates must be final and genes sorted before renaming, as ids are numbered by position
        self.sort_genes(processes=sort_processes)

        # Named before shared UTRs are homogenised, so a derived UTR takes the id of a read one with its coordinates
        for genes in self.chrs.values():
            for g in genes.values():
                provisional = [u for t in g.transcripts.values() for c in t.CDSs.values() for u in c.UTRs if u.provisional]
                if provisional:
                    g.name_new_UTRs(provisional)
        self.homogenise_shared_exons_utrs()
        # Resolved before renaming, so repeated input ids are recorded whether or not they are renamed
        self.resolve_subfeature_id_clashes(record_to=warnings_sink)

        if rename_features:
            self.rename_ids(features=rename_features, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id, quiet=quiet)

        if remove_missing_transcript_parent_references:
            self.remove_missing_transcript_parent_references(quiet=quiet)

        if update_gene_and_transcript_list:
            self.update_gene_and_transcript_list(quiet=quiet)

        if define_synteny:
            self.define_synteny(sort_processes=sort_processes)
        if self.original_annotation is not None:
            for g_id in self.original_annotation.all_gene_ids:
                if g_id not in self.all_gene_ids:
                    self.unmapped.append(g_id)

        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nWhole update process for {self.id} annotation object took {round(lapse/60, 1)} minutes\n")

    def update_features(self, standardise=False, quiet:bool=True):

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
                            if t.feature != "mRNA":
                                self.tags.add("standardised_features")
                                t.feature = "mRNA"
                    if t.feature not in self.features:
                        self.features[t.feature] = 1
                    else:
                        self.features[t.feature] += 1
                    for e in t.exons:
                        if standardise:
                            if e.feature != "exon":
                                self.tags.add("standardised_features")
                                e.feature = "exon"
                        if e.feature not in self.features:
                            self.features[e.feature] = 1
                        else:
                            self.features[e.feature] += 1
                    for c in t.CDSs.values():
                        if standardise:
                            c.feature = "CDS"
                            for cs in c.CDS_segments:
                                if cs.feature != "CDS":
                                    self.tags.add("standardised_features")
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
            if standardise:
                print(f"\nUpdated standardised features for {self.id}")
            else:
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
                        t.quality._blast_hits = None
                        t.main = False
                        t.temp_CDSs = []
                        t.temp_UTRs = []
                        t.read_UTRs_CDS = None

        if clean:
            self.remove_other_mRNA_transcripts_from_rRNA_genes()
            self.update(rename_features=("transcript", "CDS", "exon", "UTR"))
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

    def _check_exon_overlaps(self, t:Transcript, record_to:dict[str, set[str]]):
        """
        Flags transcripts with overlapping or directly adjacent (abutting) exons, which collapse_exons merges.
        """
        cur_end = None
        for e in sorted(t.exons):
            if cur_end is not None:
                if e.start <= cur_end:
                    record_to["overlapping_exons"].add(t.id)
                elif e.start == cur_end + 1:
                    record_to["abutting_exons"].add(t.id)
            cur_end = e.end if cur_end is None else max(cur_end, e.end)

    def _check_read_UTRs(self, t:Transcript, kept:bool, read:set[tuple[int, int]], record_to:dict[str, set[str]]):
        """
        Flags transcripts whose UTRs read from the GFF were replaced by UTRs derived from the exons: with
        consider_read_utrs because they do not fit the exons and CDS (read_utrs_inconsistent), and otherwise
        when the derived UTRs have other coordinates (read_utrs_differ_from_derived).
        """
        if kept:
            if t.read_UTRs_CDS is None:
                record_to["read_utrs_inconsistent"].add(t.id)
        elif read != {(u.start, u.end) for c in t.CDSs.values() for u in c.UTRs}:
            record_to["read_utrs_differ_from_derived"].add(t.id)

    def _check_phase_mismatch_across_intron(self, c:CDS, record_to:dict[str, set[str]]):
        """
        Flags CDSs whose segment phases are inconsistent across an intron (see CDS.inconsistent_junctions).
        """
        if any(not within_exon for _, within_exon in c.inconsistent_junctions()):
            record_to["phase_mismatch_across_intron"].add(c.id)

    def recheck_warnings(self, quiet:bool=True) -> dict[str, set[str]]:
        """
        Runs the checks that describe the annotation object (coordinate hierarchy, phases across
        introns, overlapping/abutting exons, genes with no transcripts and transcripts with no exons) on its current state and
        returns the findings, keyed as in self.warnings. Unlike self.warnings, which describes the
        input annotation, this reflects any modifications made since loading. self.warnings is not
        modified. Checks only possible while parsing the input (e.g. parent references) are not included.
        """
        rechecked = self._new_warnings_dict(self.RECHECKABLE_WARNING_KEYS)
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    self._check_exon_overlaps(t, rechecked)
                    for c in t.CDSs.values():
                        self._check_phase_mismatch_across_intron(c, rechecked)
        self.correct_gene_transcript_and_subfeature_coordinates(mode="none", quiet=quiet, record_to=rechecked)
        self.detect_genes_with_no_transcripts(quiet=quiet, record_to=rechecked)
        self.detect_transcripts_with_no_exons(quiet=quiet, record_to=rechecked)
        return rechecked

    def correct_gene_transcript_and_subfeature_coordinates(self, mode:Literal["none", "contain", "trim"]="contain", quiet:bool=True, record_to:dict[str, set[str]]|None=None):
        """
        Enforces the coordinate hierarchy gene >= transcript >= exons >= CDSs/UTRs.

        mode:
            "none": coordinates are left untouched; inconsistencies are recorded into record_to.
            "contain": CDSs and UTRs are kept within their transcript's exons (by extending
                the first/last exon, or the one exon a CDS row reaching into an intron overlaps;
                rows within an intron or spanning one are only reported; CDS rows beyond the
                first/last exon and separated from it by a gap, which the CDS defines as an intron,
                are given their own exon instead), transcripts span exactly
                their exons and genes are
                extended to cover all their transcripts. Gene regions extending beyond
                their transcripts (e.g. polyA sites or TF binding sites) are kept.
            "trim": same as "contain", but genes are also trimmed to span exactly from
                the start of their earliest transcript to the end of their latest one.

        record_to: warnings dict where inconsistencies are recorded in "none" mode (not
            recorded if None). Correcting modes do not record, as update() always checks
            the unpolished coordinates first.
        """
        if mode not in ("none", "contain", "trim"):
            raise ValueError(f"Invalid coordinate polishing mode: '{mode}'. Choose from: 'none', 'contain', 'trim'.")

        skip_correction = mode == "none"
        if record_to is None:
            record_to = self._new_warnings_dict()

        if not quiet:
            if skip_correction:
                print(f"Checking feature coordinates (skipping correction) for {self.id}")
            else:
                print(f"Correcting feature coordinates ({mode} mode) for {self.id}")

        # Features are counted under a summary printed once at the end, and printed individually only with list_load_warnings
        list_load_warnings = getattr(self, "list_load_warnings", False)
        flagged = {}
        def note(summary:str, warning_key:str, feature_id:str, message:str):
            if quiet:
                return
            flagged.setdefault((summary, warning_key), set()).add(feature_id)
            if list_load_warnings:
                print(message)

        for genes in self.chrs.values():
            for g in genes.values():
                earliest_start = None
                latest_end = None
                transcripts_moved = False
                extended_exons = []
                new_exons = []
                for t in g.transcripts.values():

                    for c in t.CDSs.values():
                        if hasattr(c, 'CDS_segments') and c.CDS_segments:
                            if len(c.CDS_segments) > 1:
                                c.CDS_segments.sort()
                            seg_start = c.CDS_segments[0].start
                            seg_end = max(cs.end for cs in c.CDS_segments)

                            if c.start != seg_start:
                                if skip_correction:
                                    record_to["cds_segment_coordinate_mismatch"].add(c.id)
                                    note("CDSs have a start or end that differs from their CDS_segments", "cds_segment_coordinate_mismatch", c.id, f"Warning: {c.id} start differs from its first CDS_segment ({self.id})")
                                else:
                                    note("CDSs had their start or end set to those of their CDS_segments", "cds_segment_coordinate_mismatch", c.id, f"Warning: {c.id} start differs from its first CDS_segment, proceeding to fix {self.id}")
                                    c.start = seg_start

                            if c.end != seg_end:
                                if skip_correction:
                                    record_to["cds_segment_coordinate_mismatch"].add(c.id)
                                    note("CDSs have a start or end that differs from their CDS_segments", "cds_segment_coordinate_mismatch", c.id, f"Warning: {c.id} end differs from its last CDS_segment ({self.id})")
                                else:
                                    note("CDSs had their start or end set to those of their CDS_segments", "cds_segment_coordinate_mismatch", c.id, f"Warning: {c.id} end differs from its last CDS_segment, proceeding to fix {self.id}")
                                    c.end = seg_end

                    if t.exons:
                        # CDS rows must not reach into introns (parts beyond the first/last exon are checked below as
                        # cds_exceeds_exon). When correcting, a row reaching past the one exon it overlaps extends that exon;
                        # rows within an intron, spanning one, or whose fix would make two exons touch are only reported
                        span_start = min(e.start for e in t.exons)
                        span_end = max(e.end for e in t.exons)
                        for c in t.CDSs.values():
                            for cs in c.CDS_segments:
                                start, end = max(cs.start, span_start), min(cs.end, span_end)
                                if start > end or any(e.start <= start and end <= e.end for e in t.exons):
                                    continue
                                overlapping = [e for e in t.exons if e.start <= end and start <= e.end]
                                if skip_correction:
                                    record_to["cds_segment_in_intron"].add(c.id)
                                    note("CDSs have segments reaching into an intron", "cds_segment_in_intron", c.id, f"Warning: {c.id} segment {cs.start}..{cs.end} reaches into a {t.id} intron ({self.id})")
                                    continue
                                if len(overlapping) == 1:
                                    exon = overlapping[0]
                                    new_start, new_end = min(exon.start, start), max(exon.end, end)
                                    # An extension touching another exon would remove a splice junction, so it is only reported
                                    if not any(e is not exon and e.start <= new_end + 1 and new_start - 1 <= e.end for e in t.exons):
                                        note("CDSs had segments reaching into an intron, fixed by extending the exon they overlap", "cds_segment_in_intron", c.id, f"Warning: {c.id} segment {cs.start}..{cs.end} reaches into a {t.id} intron, proceeding to extend exon {exon.id} {self.id}")
                                        exon.start, exon.end = new_start, new_end
                                        extended_exons.append(exon)
                                        continue
                                note("CDSs have segments reaching into an intron that cannot be fixed by extending one exon", "cds_segment_in_intron", c.id, f"Warning: {c.id} segment {cs.start}..{cs.end} reaches into a {t.id} intron and cannot be fixed by extending one exon ({self.id})")

                        if len(t.exons) > 1:
                            t.exons.sort()
                        # Exons are sorted by (start, end), so the first one has the earliest start,
                        # but with overlapping/nested exons the last one may not have the latest end
                        first_exon = t.exons[0]
                        last_exon = max(t.exons, key=lambda e: e.end)
                        for c in t.CDSs.values():
                            if skip_correction:
                                if c.start < first_exon.start:
                                    record_to["cds_exceeds_exon"].add(c.id)
                                    note("CDSs extend beyond the first or last exon of their transcript", "cds_exceeds_exon", c.id, f"Warning: {c.id} start is earlier than for first {t.id} exon ({self.id})")
                                if c.end > last_exon.end:
                                    record_to["cds_exceeds_exon"].add(c.id)
                                    note("CDSs extend beyond the first or last exon of their transcript", "cds_exceeds_exon", c.id, f"Warning: {c.id} end extends beyond the last {t.id} exon ({self.id})")
                            elif c.start < first_exon.start or c.end > last_exon.end:
                                # Segments next to the outermost exon extend it; segments beyond an intron (a gap
                                # between CDS segments) get their own exon, as extending would absorb the intron
                                extended, added = t.cover_CDS_segments_beyond_exons(c)
                                if extended:
                                    note("CDSs extended beyond the first or last exon of their transcript, which was extended to cover them", "cds_exceeds_exon", c.id, f"Warning: {c.id} should not extend beyond the first or last {t.id} exon, proceeding to extend it {self.id}")
                                    extended_exons.extend(extended)
                                if added:
                                    segments = ", ".join(f"{e.start}..{e.end}" for e in added)
                                    note("CDSs had segments beyond an intron past the first or last exon of their transcript, which were given their own exon", "cds_exceeds_exon", c.id, f"Warning: {c.id} segments {segments} lie beyond an intron past the first or last {t.id} exon, proceeding to add an exon for each {self.id}")
                                    new_exons.extend(added)
                                first_exon = t.exons[0]
                                last_exon = max(t.exons, key=lambda e: e.end)

                            for u in c.UTRs:
                                if u.start < first_exon.start:
                                    if skip_correction:
                                        record_to["utr_exceeds_exon"].add(u.id)
                                        note("UTRs extend beyond the first or last exon of their transcript", "utr_exceeds_exon", u.id, f"Warning: {u.id} start is earlier than for first {t.id} exon ({self.id})")
                                    else:
                                        note("UTRs extended beyond the first or last exon of their transcript, which was extended to cover them", "utr_exceeds_exon", u.id, f"Warning: {u.id} start should not be earlier than for first {t.id} exon, proceeding to fix {self.id}")
                                        first_exon.start = u.start
                                        extended_exons.append(first_exon)
                                if u.end > last_exon.end:
                                    if skip_correction:
                                        record_to["utr_exceeds_exon"].add(u.id)
                                        note("UTRs extend beyond the first or last exon of their transcript", "utr_exceeds_exon", u.id, f"Warning: {u.id} end extends beyond the last {t.id} exon ({self.id})")
                                    else:
                                        note("UTRs extended beyond the first or last exon of their transcript, which was extended to cover them", "utr_exceeds_exon", u.id, f"Warning: {u.id} end should not extend beyond the last {t.id} exon, proceeding to fix {self.id}")
                                        last_exon.end = u.end
                                        extended_exons.append(last_exon)

                        if skip_correction:
                            boundary_summary = "transcripts do not start and end with their first and last exons"
                        else:
                            boundary_summary = "transcripts were set to span exactly from their first to their last exon"
                        if first_exon.start < t.start:
                            note(boundary_summary, "exon_transcript_boundary_mismatch", t.id, f"First exon start should not be earlier than for {t.id}, proceeding to fix {self.id}" if not skip_correction else f"Warning: First exon start is earlier than for {t.id} ({self.id})")
                        elif first_exon.start > t.start:
                            note(boundary_summary, "exon_transcript_boundary_mismatch", t.id, f"First exon should not start later than {t.id}, proceeding to fix {self.id}" if not skip_correction else f"Warning: First exon starts later than for {t.id} ({self.id})")
                        if last_exon.end < t.end:
                            note(boundary_summary, "exon_transcript_boundary_mismatch", t.id, f"Last exon should not finish earlier than {t.id}, proceeding to fix {self.id}" if not skip_correction else f"Warning: Last exon finishes earlier than for {t.id} ({self.id})")
                        elif last_exon.end > t.end:
                            note(boundary_summary, "exon_transcript_boundary_mismatch", t.id, f"Last exon should not finish later than {t.id}, proceeding to fix {self.id}" if not skip_correction else f"Warning: Last exon finishes later than for {t.id} ({self.id})")
                        if t.start != first_exon.start or t.end != last_exon.end:
                            if skip_correction:
                                record_to["exon_transcript_boundary_mismatch"].add(t.id)
                            else:
                                t.start = first_exon.start
                                t.end = last_exon.end
                                transcripts_moved = True

                    if t.start < g.start:
                        if skip_correction:
                            record_to["transcript_exceeds_gene"].add(t.id)
                            note("transcripts extend beyond their gene", "transcript_exceeds_gene", t.id, f"Warning: {t.id} start is earlier than for {g.id} ({self.id})")
                        else:
                            note("transcripts extended beyond their gene, which was extended to cover them", "transcript_exceeds_gene", t.id, f"{t.id} start should not be earlier than for {g.id}, proceeding to fix {self.id}")
                            g.start = t.start
                    if t.end > g.end:
                        if skip_correction:
                            record_to["transcript_exceeds_gene"].add(t.id)
                            note("transcripts extend beyond their gene", "transcript_exceeds_gene", t.id, f"Warning: {t.id} end extends beyond {g.id} ({self.id})")
                        else:
                            note("transcripts extended beyond their gene, which was extended to cover them", "transcript_exceeds_gene", t.id, f"{t.id} end should not extend beyond {g.id}, proceeding to fix {self.id}")
                            g.end = t.end
                    if earliest_start is None or latest_end is None:
                        earliest_start = t.start
                        latest_end = t.end
                    else:
                        if t.start < earliest_start:
                            earliest_start = t.start
                        if t.end > latest_end:
                            latest_end = t.end

                if transcripts_moved and len(g.transcripts) > 1:
                    g.sort_transcripts()

                # Named first, so that extended exons given a new id below avoid the new exons' ids
                if new_exons:
                    g.name_new_exons(new_exons)
                # An extended exon may now differ from a copy shared with another transcript under the same id
                if extended_exons:
                    g.resolve_exon_id_clashes(extended_exons)

                if earliest_start is not None and latest_end is not None:
                    if g.start < earliest_start or g.end > latest_end:
                        if skip_correction:
                            record_to["gene_exceeds_transcripts"].add(g.id)
                            note("genes extend beyond their transcripts", "gene_exceeds_transcripts", g.id, f"Warning: {g.id} extends beyond its transcripts ({self.id})")
                        elif mode == "trim":
                            note("genes extended beyond their transcripts and were trimmed to them", "gene_exceeds_transcripts", g.id, f"{g.id} extends beyond its transcripts and was trimmed to them ({self.id})")
                            g.start = earliest_start
                            g.end = latest_end
                        else:
                            note("genes extend beyond their transcripts, keeping gene boundaries", "gene_exceeds_transcripts", g.id, f"{g.id} extends beyond its transcripts, keeping gene boundaries ({self.id})")
        if not quiet:
            # Only the check (mode "none") points to the ids, as update() runs it before any correction of the same features
            for (summary, warning_key), ids in flagged.items():
                if skip_correction:
                    print(f"Warning: {len(ids)} {summary} ({self.id})" + ("" if list_load_warnings else Annotation._warning_ids_hint(warning_key)))
                else:
                    print(f"{len(ids)} {summary} ({self.id})")
            if skip_correction:
                print(f"Checked feature coordinates (skipped correction) for {self.id}")
            else:
                print(f"Corrected feature coordinates ({mode} mode) for {self.id}")

    def generate_promoters(self, promoter_size:int=2000, promoter_type:str = "standard"):
        """
        promoter_type (str): Defines the reference point for the promoters.
            - standard (default): Promoter based on 'promoter_size' is generated upstream of the transcript's start site (TSS)
            - upstream_ATG : Promoter based on 'promoter_size' is generated upstream of the main CDS's start codon (ATG). If no CDS, falls back to standard.
            - standard_plus_up_to_ATG : Promoter based on 'promoter_size' is generated upstream of the transcript's start site (TSS) and any gene sequence up to the start codon (ATG) is also added. If no CDS, falls back to standard.
        """
        self.promoter_type = promoter_type
        self.promoter_size = promoter_size

        self.contains_promoters = True

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    t.generate_promoter(promoter_size, self.genome.scaffolds[t.ch].size, promoter_type)
    
    def clear_promoters(self):
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    t.clear_promoter()
        self.contains_promoters = False

    def _store_genetic_codes(
        self,
        taxonomy: Literal["plant", "vertebrate", "invertebrate", "yeast"] | str = "plant",
        table: int | str | None = None,
        mito_table: int | str | None = None,
        plastid_table: int | str | None = None,
        auto_organelle_detection: bool = True,
        mitochondria_chroms: list[str] | tuple[str, ...] | str | None = None,
        chloroplast_chroms: list[str] | tuple[str, ...] | str | None = None,
    ):
        self.taxonomy = taxonomy
        self.table, self.mito_table, self.plastid_table = resolve_taxonomy_tables(
            taxonomy=taxonomy,
            table=table,
            mito_table=mito_table,
            plastid_table=plastid_table,
        )
        self.auto_organelle_detection = auto_organelle_detection
        parsed = []
        for chroms in (mitochondria_chroms, chloroplast_chroms):
            if not chroms:
                parsed.append(())
            elif isinstance(chroms, str):
                parsed.append(tuple(c.strip() for c in chroms.split(",") if c.strip()))
            else:
                parsed.append(tuple(chroms))
        self.mitochondria_chroms, self.chloroplast_chroms = parsed

    def set_genetic_codes(
        self,
        taxonomy: Literal["plant", "vertebrate", "invertebrate", "yeast"] | str = "plant",
        table: int | str | None = None,
        mito_table: int | str | None = None,
        plastid_table: int | str | None = None,
        auto_organelle_detection: bool = True,
        mitochondria_chroms: list[str] | tuple[str, ...] | str | None = None,
        chloroplast_chroms: list[str] | tuple[str, ...] | str | None = None,
    ):
        """
        Replaces the genetic code configuration used to translate CDSs. Tables
        left as None follow the taxonomy preset. Existing proteins are cleared
        if the configuration changes, since they were translated with the old codes.
        """
        previous = self.genetic_codes
        self._store_genetic_codes(taxonomy=taxonomy, table=table, mito_table=mito_table, plastid_table=plastid_table, auto_organelle_detection=auto_organelle_detection, mitochondria_chroms=mitochondria_chroms, chloroplast_chroms=chloroplast_chroms)
        self._validate_organelle_chroms()
        if self.genetic_codes != previous:
            self.clear_proteins()

    @property
    def genetic_codes(self) -> dict:
        """The genetic code configuration used to translate CDSs."""
        return {
            "taxonomy": self.taxonomy,
            "table": self.table,
            "mito_table": self.mito_table,
            "plastid_table": self.plastid_table,
            "auto_organelle_detection": self.auto_organelle_detection,
            "mitochondria_chroms": self.mitochondria_chroms,
            "chloroplast_chroms": self.chloroplast_chroms,
        }

    def _validate_organelle_chroms(self):
        all_known = set(self.chrs.keys())
        if self.genome:
            all_known.update(self.genome.scaffolds.keys())
        for kind, chroms in (("mitochondrial", self.mitochondria_chroms), ("chloroplast", self.chloroplast_chroms)):
            for chrom in chroms:
                if chrom not in all_known:
                    raise ValueError(f"Specified {kind} chromosome '{chrom}' was not found in the annotation or genome.")

    def chromosome_compartment(self, chrom: str) -> Literal["nuclear", "mitochondria", "chloroplast"]:
        """
        Compartment of a chromosome, deciding its genetic code, its gene ids (see rename_ids) and how it is
        reported: user-listed organelle chromosomes first, then (with auto_organelle_detection) the
        classification of its genome scaffold or, without one, of its name (see classify_compartment).
        """
        scaffold = self.genome.get_scaffold(chrom) if self.genome is not None else None
        return classify_compartment(chrom, scaffold=scaffold, mitochondria_chroms=self.mitochondria_chroms, chloroplast_chroms=self.chloroplast_chroms, auto_detection=self.auto_organelle_detection)

    def translation_table(self, chrom: str, cds: CDS | None = None) -> int | str:
        """
        Returns the genetic code table used to translate CDSs on a chromosome. A CDS given with its
        own GFF transl_table attribute is translated with that table instead.
        """
        if cds is not None and cds.transl_table is not None:
            return cds.transl_table
        compartment = self.chromosome_compartment(chrom)
        if compartment == "mitochondria":
            return self.mito_table
        if compartment == "chloroplast":
            return self.plastid_table
        return self.table

    def _print_genetic_code_info(self):
        names = {}
        for key in ("table", "mito_table", "plastid_table"):
            t = getattr(self, key)
            names[key] = NCBI_GENETIC_CODES.get(t, {}).get("name", "Custom") if isinstance(t, int) else "Custom"
        print(f"Info: Genetic code configuration [taxonomy='{self.taxonomy}']: nuclear={self.table} ({names['table']}), mitochondrial={self.mito_table} ({names['mito_table']}), plastid={self.plastid_table} ({names['plastid_table']})")

    def generate_proteins(
        self,
        mode: Literal["start", "end", "orf", "orf_or_end", "orf_or_start"] = "end",
        max_nucleotide_trim: int | None = None,
        tolerated_stops: int | None = 0,
        orf_choice_mode: Literal["longest", "earliest"] = "longest",
        must_have_stop: bool = False,
        enforce_start_codon: bool = True,
        min_codon_len: int = 2,
        start_codons: tuple[str, ...] | None = None,
        stop_codons: tuple[str, ...] | None = None,
        correct_CDS: bool = False,
        always_resolve_strand: bool = True,
        ignore_ambiguous_strands: bool = False,
        quiet: bool = True,
    ):
        if not quiet and self.auto_organelle_detection:
            self._print_genetic_code_info()

        for chrom, genes in self.chrs.items():
            if self.genome is not None and chrom not in self.genome.scaffolds:
                continue

            for g in genes.values():
                for t in g.transcripts.values():
                    for c in t.CDSs.values():
                        c.generate_protein(
                            mode=mode,
                            max_nucleotide_trim=max_nucleotide_trim,
                            tolerated_stops=tolerated_stops,
                            orf_choice_mode=orf_choice_mode,
                            must_have_stop=must_have_stop,
                            enforce_start_codon=enforce_start_codon,
                            min_codon_len=min_codon_len,
                            start_codons=start_codons,
                            stop_codons=stop_codons,
                            correct_CDS=correct_CDS,
                            always_resolve_strand=always_resolve_strand,
                            ignore_ambiguous_strands=ignore_ambiguous_strands,
                            quiet=quiet,
                            adjust_internal_shifts=self.adjust_internal_shifts,
                            table=self.translation_table(chrom, c),
                            initiator_methionine=self.initiator_methionine,
                            alt_start_methionine=self.alt_start_methionine,
                        )
                    if correct_CDS:
                        t.update(quiet=quiet)
                if correct_CDS:
                    g.update(quiet=quiet)
        self.contains_protein_sequences = True

    def get_protein_qc_summary(self) -> dict[str, int]:
        """
        Computes summary metrics for CDSs and translated proteins across the annotation. Protein
        metrics are counted from Protein.tags, as written in detailed FASTA headers:
        - total_cds: Total CDS features.
        - total_proteins: CDSs with translated protein sequences.
        - complete_proteins: Translated proteins with valid start codon, stop codon and no early stop.
        - alt_start_proteins: Complete proteins whose start codon is an alternative (non-curated) NCBI initiator.
        - partial_proteins: Proteins missing a start codon or a stop codon.
        - partial_5prime: Proteins missing valid start codon (or with skipped 5' CDS bases).
        - partial_3prime: Proteins missing stop codon.
        - ambiguous_proteins: Proteins with residues translated from ambiguous codons (X).
        - truncated_proteins: Proteins with premature internal stop codons.
        - frameshifted_proteins: Proteins translated across a frameshift (annotated overlap or gap, or adjusted phases).
        - frameshifted_cds: CDSs with internal phase shifts / frameshifts (intra-exon or intron); needs no genome.
        - phase_mismatches: CDSs with phase mismatch across introns.
        """
        if hasattr(self, "_protein_qc_summary") and self._protein_qc_summary is not None:
            return self._protein_qc_summary

        if self.genome is not None and not getattr(self, "contains_protein_sequences", False):
            try:
                self.generate_proteins(quiet=True)
            except Exception:
                pass

        total_cds = 0
        total_proteins = 0
        complete_proteins = 0
        alt_start_proteins = 0
        ambiguous_proteins = 0
        partial_proteins = 0
        partial_5p = 0
        partial_3p = 0
        truncated_proteins = 0
        frameshifted_proteins = 0
        frameshifted_cds = 0

        intron_phase_mismatches = self.warnings.get("phase_mismatch_across_intron", set())

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    for c in t.CDSs.values():
                        total_cds += 1
                        if getattr(c, "has_internal_shift", False) or c.id in intron_phase_mismatches:
                            frameshifted_cds += 1
                        p = c.protein
                        if p is not None and p.seq:
                            total_proteins += 1
                            tags = p.tags
                            partial = tags.get("partial", "")
                            if partial:
                                partial_proteins += 1
                                if "5'" in partial:
                                    partial_5p += 1
                                if "3'" in partial:
                                    partial_3p += 1
                            if "internal_stops" in tags:
                                truncated_proteins += 1
                            elif not partial:
                                complete_proteins += 1
                                if "alt_start" in tags:
                                    alt_start_proteins += 1
                            if "ambiguous_residues" in tags:
                                ambiguous_proteins += 1
                            if "frameshifts" in tags:
                                frameshifted_proteins += 1

        summary = {
            "total_cds": total_cds,
            "total_proteins": total_proteins,
            "complete_proteins": complete_proteins,
            "alt_start_proteins": alt_start_proteins,
            "partial_proteins": partial_proteins,
            "partial_5prime": partial_5p,
            "partial_3prime": partial_3p,
            "truncated_proteins": truncated_proteins,
            "ambiguous_proteins": ambiguous_proteins,
            "frameshifted_proteins": frameshifted_proteins,
            "frameshifted_cds": frameshifted_cds,
            "phase_mismatches": len(intron_phase_mismatches),
        }
        self._protein_qc_summary = summary
        return summary

    def generate_protein_equivalences(
        self,
        quiet: bool = True,
    ):
        if not self.contains_protein_sequences:
            self.generate_proteins(quiet=quiet)

        all_protein_seqs = {}
        self.all_protein_ids = {}
        for chrom, genes in self.chrs.items():
            for g in genes.values():
                if g.coding:
                    for t in g.transcripts.values():
                        for c in t.CDSs.values():
                            if c.protein is not None:
                                if c.protein.seq != "":
                                    all_protein_seqs[c.protein.id] = c.protein.seq
                                    self.all_protein_ids[c.protein.id] = (chrom, g.id, t.id, c.id)

        progress_bar = start_progress_bar(total=len(all_protein_seqs.keys()), description=f"Determining unique {self.id} proteins", colour="91", quiet=quiet)

        unique_sequences = {}
        self.protein_equivalences = {}

        for protein_id, sequence in all_protein_seqs.items():
            progress_bar.update(1)
            if sequence not in unique_sequences:
                unique_sequences[sequence] = protein_id
                self.protein_equivalences[protein_id] = []
            else:
                first_protein_id = unique_sequences[sequence]
                self.protein_equivalences[first_protein_id].append(protein_id)

        progress_bar.close()

    def clear_proteins(self):
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    for c in t.CDSs.values():
                        c.clear_protein()
        self.contains_protein_sequences = False
        self.protein_equivalences = {}
        self._protein_qc_summary = None

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

    def combine_transcripts(self, respect_non_coding:bool=False, respect_non_combined:bool=False, redetect_CDS:bool=True, quiet:bool=False):
        for genes in self.chrs.values():
            for g in genes.values():
                g.combine_transcripts(respect_non_coding=respect_non_coding, respect_non_combined=respect_non_combined, redetect_CDS=redetect_CDS, quiet=quiet)
        self.update(rename_features=("transcript", "CDS", "exon", "UTR"), quiet=quiet)
        self.feature_tags.add("combined")

    def sort_genes(self, processes:int=2, quiet:bool=True, noisy:bool=False):
        if not quiet:
            print(f"\nSorting genes for {self.id}")

        progress_bar = start_progress_bar(total=len(self.all_gene_ids), description=f"Sorting {self.id} genes", colour="38;2;210;180;140", quiet=quiet)

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
        if not quiet:
            print(f"Sorted genes for {self.id}")

    def define_synteny(self, sort_processes:int=1, quiet:bool=True):
        if not quiet:
            print(f"\nDefining synteny for {self.id} annotation genes")
        start_time = time.time()
        self.sort_genes(processes=sort_processes)
        # defining synteny
        for genes in self.chrs.values():
            gene_list = list(genes.values())
            for n, g in enumerate(gene_list):
                g.synteny.order = n
                # This works even if only a single gene has been annotated
                # in a chromosome or scaffold
                if len(genes) == 1:
                    g.synteny.previous = None
                    g.synteny.next = None
                elif n == 0:
                    g.synteny.previous = None
                    g.synteny.next = gene_list[n+1].id
                elif n != (len(genes) - 1):
                    g.synteny.previous = gene_list[n-1].id
                    g.synteny.next = gene_list[n+1].id
                else:
                    g.synteny.previous = gene_list[n-1].id
                    g.synteny.next = None

        if self.original_annotation is not None:
            for genes in self.chrs.values():
                for g_id, g in genes.items():
                    # this extra bit is for extra liftover copies
                    if g_id not in self.original_annotation.all_gene_ids:
                        g_id = "_".join(g_id.split("_")[:-1])
                        if g_id not in self.original_annotation.all_gene_ids:
                            continue

                    g.synteny.old_previous = self.original_annotation.chrs[self.original_annotation.all_gene_ids[g_id]][g_id].synteny.previous
                    g.synteny.old_next = self.original_annotation.chrs[self.original_annotation.all_gene_ids[g_id]][g_id].synteny.next                
                    g.synteny.old_order = self.original_annotation.chrs[self.original_annotation.all_gene_ids[g_id]][g_id].synteny.order

                    if g.synteny.old_previous == g.synteny.previous and g.synteny.old_next == g.synteny.next:
                        g.synteny.liftover_conserved = True
                    else: 
                        g.synteny.liftover_conserved = False
        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nDefining synteny for {self.id} annotation genes took {round(lapse, 1)} seconds\n")

    def homogenise_shared_exons_utrs(self):
        """
        Exons (and UTRs) of a gene with the same coordinates are one shared feature: they take the id of
        the first one, whether ids were renamed or kept, and their parents are the transcripts of the gene
        containing them (references to transcripts outside the gene are kept, see
        remove_missing_transcript_parent_references), so features no longer shared after being merged
        or extended lose the transcripts that no longer contain them.
        Different features left with the same id (e.g. ids numbered per transcript, such as exon_3 being
        a different exon in each transcript) keep it for the first one; the others get <id>_<n>.
        """
        for genes in self.chrs.values():
            for g in genes.values():
                transcript_ids = {t.id for t in g.transcripts.values()}
                exon_groups: dict[tuple[int, int, str, str], list[tuple[str, Exon]]] = {}
                utr_groups: dict[tuple[int, int, str, str], list[tuple[str, UTR]]] = {}
                for t in g.transcripts.values():
                    for e in t.exons:
                        exon_groups.setdefault((e.start, e.end, e.ch, e.strand), []).append((t.id, e))
                    for c in t.CDSs.values():
                        for u in c.UTRs:
                            utr_groups.setdefault((u.start, u.end, u.ch, u.strand), []).append((t.id, u))

                for groups in (exon_groups, utr_groups):
                    for group in groups.values():
                        first = group[0][1]
                        containing = {t_id for t_id, _ in group}
                        outside_gene = {p for _, f in group for p in (f.parents or []) if p not in transcript_ids}
                        parents = sorted(containing | outside_gene)
                        for _, f in group:
                            f.parents = parents[:]
                            f.id, f.id_number = first.id, first.id_number

                    owners = {}
                    used_ids = {group[0][1].id for group in groups.values()}
                    for key, group in groups.items():
                        f_id = group[0][1].id
                        if owners.setdefault(f_id, key) == key:
                            continue
                        n = 2
                        while f"{f_id}_{n}" in used_ids:
                            n += 1
                        used_ids.add(f"{f_id}_{n}")
                        for _, f in group:
                            f.id = f"{f_id}_{n}"
                            f.update_numbering()

        self.shared_exons = True
        self.shared_UTRs = True

    def resolve_subfeature_id_clashes(self, record_to:dict[str, set[str]]|None=None):
        """
        Exon, CDS and UTR ids are only reconciled within a gene (see homogenise_shared_exons_utrs), but
        an id may only be used by one feature of the annotation. An id is kept by the first gene (in
        chromosome order) and feature type using it; later uses, and uses of a gene or transcript id, get
        <id>_<n>, shared by the features of that gene and type with that id. Repeated ids read from the
        input are recorded as repeated_subfeature_ids.
        """
        all_ids = set()
        for genes in self.chrs.values():
            for g in genes.values():
                all_ids.add(g.id)
                for t in g.transcripts.values():
                    all_ids.add(t.id)
                    all_ids.update(e.id for e in t.exons)
                    for c in t.CDSs.values():
                        all_ids.add(c.id)
                        all_ids.update(cs.id for cs in c.CDS_segments)
                        all_ids.update(u.id for u in c.UTRs)
        taken = {g.id for genes in self.chrs.values() for g in genes.values()}
        taken.update(t.id for genes in self.chrs.values() for g in genes.values() for t in g.transcripts.values())

        for genes in self.chrs.values():
            for g in genes.values():
                groups: dict[tuple[str, str], list[Feature]] = {}
                for t in g.transcripts.values():
                    for e in t.exons:
                        groups.setdefault(("exon", e.id), []).append(e)
                    for c in t.CDSs.values():
                        groups.setdefault(("CDS", c.id), []).append(c)
                        for cs in c.CDS_segments:
                            groups.setdefault(("CDS", cs.id), []).append(cs)
                        for u in c.UTRs:
                            groups.setdefault(("UTR", u.id), []).append(u)

                renamed_CDSs = False
                for (kind, f_id), features in groups.items():
                    if f_id not in taken:
                        taken.add(f_id)
                        continue
                    if record_to is not None and any(f.original_id == f_id for f in features):
                        record_to["repeated_subfeature_ids"].add(f_id)
                    n = 2
                    while f"{f_id}_{n}" in all_ids:
                        n += 1
                    new_id = f"{f_id}_{n}"
                    all_ids.add(new_id)
                    taken.add(new_id)
                    for f in features:
                        f.id = new_id
                        f.update_numbering()
                        if isinstance(f, CDS):
                            renamed_CDSs = True
                            if f.protein is not None:
                                f.protein.id = f"{new_id}.prot"

                if renamed_CDSs:
                    for t in g.transcripts.values():
                        t.CDSs = {c.id: c for c in t.CDSs.values()}

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
                                if g.aliases is None:
                                    g.aliases = []
                                if hit.id not in g.aliases:
                                    g.aliases.append(hit.id)
    
    def clear_aliases(self):
        for genes in self.chrs.values():
            for g in genes.values():
                g.aliases = None

    def CDS_to_CDS_segment_ids(self, override:bool=False):
        repeat_CDS_segment_id = False

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    for c in t.CDSs.values():
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
                        for c in t.CDSs.values():
                            if not c.main:
                                count += 1
                                c.id = f"{t.id}_CDS{count}"
                            else:
                                c.id = f"{t.id}_CDS1"
                            for x, cs in enumerate(c.CDS_segments):
                                cs.id = f"{c.id}_{x+1}"

        self.update_keys(CDS_keys=True, gene_keys=False, transcript_keys=False)

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
        
        self.update_keys(CDS_keys=True, gene_keys=False, transcript_keys=False)

    def merge(self, other:Annotation, max_cds_overlap:int|float=100, max_exon_overlap:int|float=100, max_gene_overlap:int|float=100, features_to_rename:tuple[str,...]=("gene", "transcript", "CDS", "exon", "UTR"), rename_clashing_ids:bool=True, quiet:bool=False):

        """
        Priority is given to self annotation
        """
        start_time = time.time()
        self.update(quiet=quiet)
        other.update(quiet=quiet)

        if rename_clashing_ids:

            total_features_to_rename = list(features_to_rename)
            
            total_features_to_rename.extend(["transcript", "CDS", "exon", "UTR"])

            features_to_rename = tuple(set(total_features_to_rename))

        if max_cds_overlap != 100 or max_exon_overlap != 100 or max_gene_overlap != 100:
            self.overlaps.detect(other, quiet=quiet)

        progress_bar = start_progress_bar(total=len(other.all_gene_ids), description=f"Adding {other.id} annotation to {self.id}", colour="38;2;46;204;113", quiet=quiet)
        
        if self.merged:
            if "_..._" in self.name:
                first_name = self.name.split("_..._")[0]
            else:
                first_name = self.name.split("_merged_")[0]
            self.name = f"{first_name}_..._{other.name}"
        else:
            self.name = f"{self.name}_merged_{other.name}"
        if self.genome is not None:
            self.id = f"{self.name}_on_{self.genome.name}"
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
                    if rename_clashing_ids:
                        while temp_id in self.all_gene_ids:
                            count += 1
                            temp_id = f"{g.id}_{count}"

                    if temp_id not in self.chrs[chr]:
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
                        if rename_clashing_ids:
                            while temp_id in self.all_gene_ids:
                                count += 1
                                temp_id = f"{g.id}_{count}"
                        if temp_id not in self.chrs[chr]:
                            self.chrs[chr][temp_id] = g.copy()
                            self.chrs[chr][temp_id].id = temp_id
                            self.all_gene_ids[temp_id] = chr
                    else:
                        cds_scores: list[float] = [0.0]
                        exon_scores: list[float] = [0.0]
                        gene_scores: list[float] = [0.0]
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
                                if rename_clashing_ids:
                                    while temp_id in self.all_gene_ids:
                                        count += 1
                                        temp_id = f"{g.id}_{count}"
                                if temp_id not in self.chrs[chr]:
                                    self.chrs[chr][temp_id] = g.copy()
                                    self.chrs[chr][temp_id].id = temp_id
                                    self.all_gene_ids[temp_id] = chr
                        
                        elif check_exons:
                            if max(exon_scores) <= max_exon_overlap and max(gene_scores) <= max_gene_overlap:
                                temp_id = g.id
                                if rename_clashing_ids:
                                    while temp_id in self.all_gene_ids:
                                        count += 1
                                        temp_id = f"{g.id}_{count}"
                                if temp_id not in self.chrs[chr]:
                                    self.chrs[chr][temp_id] = g.copy()
                                    self.chrs[chr][temp_id].id = temp_id
                                self.all_gene_ids[temp_id] = chr

                        elif max(gene_scores) <= max_gene_overlap:
                            temp_id = g.id
                            if rename_clashing_ids:
                                while temp_id in self.all_gene_ids:
                                    count += 1
                                    temp_id = f"{g.id}_{count}"
                            if temp_id not in self.chrs[chr]:
                                self.chrs[chr][temp_id] = g.copy()
                                self.chrs[chr][temp_id].id = temp_id
                            self.all_gene_ids[temp_id] = chr

        self.merged = True
  
        self.clear_proteins()
        progress_bar.close()
        self.update_gene_and_transcript_list(quiet=quiet)
        self.update(rename_features=features_to_rename, quiet=quiet)
        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nMerging {self.id} and {other.id} annotations took {round(lapse/60, 1)} minutes")

    def remove_exons_with_unmatched_strand(self, remove_transcripts_accordingly:bool=False, remove_genes_accordingly:bool=False, quiet:bool=False):
        for genes in self.chrs.values():
            genes_to_remove = set()
            for g in genes.values():
                transcripts_to_remove = set()
                for t in g.transcripts.values():
                    new_e = []
                    for e in t.exons:
                        if e.strand == t.strand:
                            new_e.append(e)
                    if len(new_e) > 0:
                        t.exons = new_e
                    else:
                        if len(t.exons) > 0 and remove_transcripts_accordingly:
                            transcripts_to_remove.add(t.id)
                if len(g.transcripts) == 0:
                    continue
                if len(transcripts_to_remove) == len(g.transcripts):
                    if remove_genes_accordingly:
                        genes_to_remove.add(g.id)
                for t_id in transcripts_to_remove:
                    del g.transcripts[t_id]
                    if not quiet:
                        print(f"{t_id} Warning: transcript with no exons was removed because of exon strand inconsistencies")
            for g_id in genes_to_remove:
                del genes[g_id]
                if not quiet:
                    print(f"{g_id} Warning: gene with no transcripts because of exon strand inconsistencies was removed")

        self.update_gene_and_transcript_list(quiet=quiet)
        self.remove_missing_transcript_parent_references(quiet=quiet)

    def detect_transcripts_with_no_exons(self, remove_transcripts:bool=False, remove_genes_accordingly:bool=False, quiet:bool=False, record_to:dict[str, set[str]]|None=None):
        """
        record_to: warnings dict where findings are recorded (not recorded if None).
        """
        if record_to is None:
            record_to = self._new_warnings_dict()
        removed = False
        for genes in self.chrs.values():
            genes_to_remove = set()
            for g in genes.values():
                transcripts_to_remove = set()
                for t in g.transcripts.values():
                    if t.exons == []:
                        record_to["transcript_with_no_exons"].add(t.id)
                        if remove_transcripts:
                            transcripts_to_remove.add(t.id)
                            removed = True
                if len(g.transcripts) == 0:
                    continue
                if len(transcripts_to_remove) == len(g.transcripts):
                    if remove_genes_accordingly:
                        genes_to_remove.add(g.id)
                for t_id in transcripts_to_remove:
                    del g.transcripts[t_id]
                    if not quiet:
                        print(f"{t_id} Warning: transcript with no exons was removed")
            for g_id in genes_to_remove:
                del genes[g_id]
                if not quiet:
                    print(f"{g_id} Warning: gene with no transcripts, deleted because of lack of exons, was removed")
                record_to["gene_with_no_transcripts"].add(g_id)

        if removed:
            self.update_gene_and_transcript_list(quiet=quiet)
            self.remove_missing_transcript_parent_references(quiet=quiet)

    def remove_transcripts(self, to_remove:set, remove_genes_accordingly:bool=False,quiet:bool=False):
        removed = False

        for t in to_remove:
            if t in self.all_transcript_ids:
                chrom = self.all_transcript_ids[t][0]
                g_id = self.all_transcript_ids[t][1]
                del self.chrs[chrom][g_id].transcripts[t]
                removed = True
                if remove_genes_accordingly:
                    if self.chrs[chrom][g_id].transcripts == {}:
                        del self.chrs[chrom][g_id]
            elif not quiet:
                warnings.warn(f"Transcript level id {t} is not present in annotation {self.id}", category=UserWarning)

        if removed:
            self.update_gene_and_transcript_list(quiet=quiet)
            self.remove_missing_transcript_parent_references(quiet=quiet)

    def detect_genes_with_no_transcripts(self, remove:bool=False, remove_pseudogene:bool=False, quiet:bool=False, record_to:dict[str, set[str]]|None=None):
        """
        record_to: warnings dict where findings are recorded (not recorded if None).
        """
        if record_to is None:
            record_to = self._new_warnings_dict()
        removed = False
        for genes in self.chrs.values():
            genes_to_remove = set()
            for g in genes.values():
                if g.transcripts == {}:
                    record_to["gene_with_no_transcripts"].add(g.id)
                    if remove:
                        removed = True
                        if g.pseudogene:
                            if remove_pseudogene:
                                genes_to_remove.add(g.id)
                        else:
                            genes_to_remove.add(g.id)

            for g_id in genes_to_remove:
                del genes[g_id]
                if not quiet:
                    print(f"{g_id} Warning: gene with no transcripts was removed")

        if removed:
            self.update_gene_and_transcript_list(quiet=quiet)

    def detect_polycistronic_transcripts(self, quiet:bool=False, record_to:dict[str, set[str]]|None=None):
        """
        Detects polycistronic transcripts, overlapping CDS segments, or independent CDS segment ID
        formatting across the annotation based on raw temp_CDSs before CDS generation, and populates
        self.warnings['overlapping_CDS_segments'], self.warnings['possible_polycistronic_transcript'],
        and self.warnings['independent_CDS_segment_ids'].

        Exemplified by canonical_gene_example.gff3 (the Sequence Ontology / GFF3 specification EDEN gene example),
        where mRNA00003 contains multiple logical CDSs (cds00003 and cds00004) under the same transcript.

        record_to: warnings dict where findings are recorded (defaults to self.warnings, as this is a load-time check).
        """
        if record_to is None:
            record_to = self.warnings
        list_load_warnings = getattr(self, "list_load_warnings", False)
        flagged = {}
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    if t.temp_CDSs:
                        warnings = t._detect_polycistronic()
                        for warning in warnings:
                            if t.id not in record_to[warning]:
                                flagged[warning] = flagged.get(warning, 0) + 1
                                if list_load_warnings and not quiet:
                                    print(f"Warning: Transcript {t.id} flagged as {warning} ({self.id})")
                            record_to[warning].add(t.id)
        if not quiet:
            for warning, count in flagged.items():
                print(f"Warning: {count} transcripts flagged as {warning} ({self.id})" + ("" if list_load_warnings else Annotation._warning_ids_hint(warning)))


    def remove_missing_transcript_parent_references(self, quiet:bool=True):

        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    for e in t.exons:
                        new_parents = []
                        if e.parents:
                            for p in e.parents:
                                if p in self.all_transcript_ids:
                                    if p not in g.transcripts or g.transcripts[p].strand == e.strand:
                                        new_parents.append(p)
                                elif not quiet:
                                    print(f"Transcript level id {p} parent reference for {e.id} is not present in annotation {self.id}")
                        e.parents = new_parents
                        e.parents.sort()
                    for c in t.CDSs.values():
                        new_parents = []
                        if c.parents:
                            for p in c.parents:
                                if p in self.all_transcript_ids:
                                    new_parents.append(p)
                            c.parents = new_parents
                            c.parents.sort()
                            for cs in c.CDS_segments:
                                new_parents = []
                                for p in (cs.parents or []):
                                    if p in self.all_transcript_ids:
                                        new_parents.append(p)
                                    elif not quiet:
                                        print(f"Transcript level id {p} parent reference for {cs.id} is not present in annotation {self.id}")
                                cs.parents = new_parents
                                cs.parents.sort()

    def rework_CDSs(
        self,
        override: bool = True,
        coding_ratio_threshold: float = 0.7,
        fallback_to_trim: bool = False,
        start_codons: tuple[str, ...] | None = None,
        stop_codons: tuple[str, ...] | None = None,
        min_codon_len: int = 2,
        quiet: bool = False,
        enforce_start_codon: bool = True,
        must_have_stop: bool = True,
        tolerated_stops: int | None = 0,
        orf_choice_mode: Literal["longest", "earliest"] = "longest",
        mode: Literal["orf", "orf_or_end", "orf_or_start"] = "orf",
        allow_internal_stops: bool = True,
        allow_partial: bool = True,
    ):
        """
        Replaces the CDSs of transcripts (all of them with override, otherwise only those of
        non-coding transcripts) by the best ORF across their exons.

        Only what the rework changes gets a new id: a reworked CDS with the segments and phases of
        the CDS it replaces is that CDS (keeping its id and attributes), and so are UTRs with the
        coordinates of replaced ones (see Transcript.keep_unchanged_subfeatures). New CDSs are named
        <transcript>_CDS001 and new UTRs <gene>_u<n> (by update, see Gene.name_new_UTRs); other
        features keep their ids.
        """
        start_time = time.time()

        if not quiet and self.auto_organelle_detection:
            self._print_genetic_code_info()

        progress_bar = start_progress_bar(total=len(self.all_gene_ids), description=f"Reworking {self.id} CDSs", colour="91", quiet=quiet)

        for chrom, genes in self.chrs.items():
            chrom_table = self.translation_table(chrom)

            for g in genes.values():
                progress_bar.update(1)
                for t in g.transcripts.values():
                    if t.coding and not override:
                        continue
                    old_CDSs = list(t.CDSs.values())
                    old_UTRs = [u for c in old_CDSs for u in c.UTRs]

                    t.generate_best_protein(
                        mode=mode,
                        start_codons=start_codons,
                        stop_codons=stop_codons,
                        min_codon_len=min_codon_len,
                        enforce_start_codon=enforce_start_codon,
                        must_have_stop=must_have_stop,
                        tolerated_stops=tolerated_stops,
                        orf_choice_mode=orf_choice_mode,
                        quiet=quiet,
                        table=chrom_table,
                        initiator_methionine=self.initiator_methionine,
                        alt_start_methionine=self.alt_start_methionine,
                    )
                    t.update(quiet=quiet)

                    if t.coding_ratio < coding_ratio_threshold and allow_internal_stops:
                        t.generate_best_protein(
                            mode=mode,
                            start_codons=start_codons,
                            stop_codons=stop_codons,
                            tolerated_stops=1 if (tolerated_stops is None or tolerated_stops == 0) else tolerated_stops,
                            min_codon_len=min_codon_len,
                            enforce_start_codon=enforce_start_codon,
                            must_have_stop=must_have_stop,
                            orf_choice_mode=orf_choice_mode,
                            quiet=quiet,
                            table=chrom_table,
                            initiator_methionine=self.initiator_methionine,
                            alt_start_methionine=self.alt_start_methionine,
                        )
                        t.update(quiet=quiet)

                    if t.coding_ratio < coding_ratio_threshold and allow_partial:
                        t.generate_best_protein(
                            mode=mode,
                            start_codons=start_codons,
                            stop_codons=stop_codons,
                            must_have_stop=False,
                            tolerated_stops=tolerated_stops,
                            min_codon_len=min_codon_len,
                            enforce_start_codon=enforce_start_codon,
                            orf_choice_mode=orf_choice_mode,
                            quiet=quiet,
                            table=chrom_table,
                            initiator_methionine=self.initiator_methionine,
                            alt_start_methionine=self.alt_start_methionine,
                        )
                        t.update(quiet=quiet)

                    if t.coding_ratio < coding_ratio_threshold and fallback_to_trim:
                        trim_mode = "orf_or_end" if mode == "orf" else mode
                        t.generate_best_protein(
                            mode=trim_mode,
                            start_codons=start_codons,
                            stop_codons=stop_codons,
                            must_have_stop=False,
                            min_codon_len=min_codon_len,
                            enforce_start_codon=enforce_start_codon,
                            orf_choice_mode=orf_choice_mode,
                            quiet=quiet,
                            table=chrom_table,
                            initiator_methionine=self.initiator_methionine,
                            alt_start_methionine=self.alt_start_methionine,
                        )
                        t.update(quiet=quiet)

                    t.keep_unchanged_subfeatures(old_CDSs, old_UTRs, quiet=quiet)

                g.update(quiet=quiet)

        progress_bar.close()
        # Names the new UTRs (see Gene.name_new_UTRs)
        self.update(quiet=quiet)
        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nReworking CDSs for {self.id} took {round(lapse/60, 1)} minutes")

    def update_gene_and_transcript_list(self, quiet:bool=True):

        total = 0
        for genes in self.chrs.values():
            total += len(genes)

        progress_bar = start_progress_bar(total=total, description=f"Updating {self.id} annotation gene and transcript lists", colour="95", quiet=quiet)

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
                # Trim the gene to its remaining transcripts, as its span may have come from those now split off
                if g.transcripts:
                    g.start = min(t.start for t in g.transcripts.values())
                    g.end = max(t.end for t in g.transcripts.values())
            new_chrs[chrom] = new_genes.copy()

        for nc in new_chrs:
            for g in new_chrs[nc].values():
                self.chrs[nc][g.id] = g.copy()

        del new_chrs

        self.update_gene_and_transcript_list(quiet=quiet)
        self.update(rename_features=("gene", "transcript", "CDS", "exon", "UTR"), quiet=quiet)

    def rename_source(self, new_source:str="aegis", atypical:bool=True, orphaned:bool=True):
        for genes in self.chrs.values():
            for g in genes.values():
                g.source = new_source
                for t in g.transcripts.values():
                    t.source = new_source
                    for e in t.exons:
                        e.source = new_source
                    if t.introns:
                        for i in t.introns:
                            i.source = new_source
                    for c in t.CDSs.values():
                        c.source = new_source
                        for cs in c.CDS_segments:
                            cs.source = new_source
                        for u in c.UTRs:
                            u.source = new_source
                    for ft in t.miRNAs:
                        ft.source = new_source

        if atypical:
            for a in self.atypical_features:
                a.source = new_source
        
        if orphaned:
            for o in self.orphaned_features:
                o.source = new_source
        
    def rename_ids(self, features:tuple[str,...]=("gene", "transcript", "CDS", "exon", "UTR"), keep_existing_ids_if_derived_from_base_id:bool=False, remove_point_suffix:bool=False, strip_gene_tag:bool=False, keep_subfeature_numbers:bool=False, cds_segment_ids:bool=False, prefix:str="", suffix:str="", spacer:int=100, sep:str="_", g_id_digits:int=5, t_id_digits:int=3, correspondences:bool=False, correspondences_output_dir:str | None = None, correspondences_extension=".tsv", correspondences_filename:str | None = None, correspondences_use_annot_dir:bool=False, correspondences_subfolder:bool=False, correspondences_subfolder_name:str="out_gffs", correspondences_filepath:str | None = None, old_id_names:Literal["alias", "keep"]="alias", quiet:bool=False,
    #deprecated arguments
    custom_path:str=""):
        """
        Renames the chosen feature levels from their gene ids.

        prefix: rename all genes as {prefix}{chromosome}g{number}{sep}{suffix}, numbered by position along
            each chromosome. Genes on organelles (see chromosome_compartment) use chrMT or chrCP as
            chromosome, numbered across all chromosomes of that organelle.
        keep_subfeature_numbers: reuse the trailing number of current transcript, CDS, exon and
            UTR ids. Exons, UTRs and CDSs only take a number written after their type
            ('AT4G00005:exon:1', 'g1_u002', 't1_CDS1'), so e.g. the accession version in
            'cds-NP_065184.2' is not read as one. Numbers are only kept when every feature of the
            group (transcripts or exons/UTRs of a gene) has one and no number is shared by
            different features; otherwise the group is numbered sequentially.
        cds_segment_ids: give each CDS row its own id (<CDS id>_<n>) instead of the CDS id.
        old_id_names: 'alias' moves a Name that is the feature's old ID (e.g. Araport 'ID=AT4G00005;Name=AT4G00005')
            to Alias, so renamed features are not labelled with an ID they no longer have; 'keep' leaves Names as read.
        """

        if old_id_names not in ("alias", "keep"):
            raise ValueError(f"old_id_names must be 'alias' or 'keep', not '{old_id_names}'.")

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'correspondences_output_dir' instead.", DeprecationWarning, stacklevel=2)
            if correspondences_output_dir is None:
                correspondences_output_dir = custom_path

        if correspondences_subfolder_name != "out_gffs":
            correspondences_subfolder = True

        if correspondences_output_dir is not None or correspondences_subfolder or correspondences_filepath is not None or correspondences_filename is not None:
            correspondences = True

        acceptable_features = ["gene", "transcript", "CDS", "exon", "UTR"]

        for f in features:
            if f not in acceptable_features:
                raise ValueError(f"Incorrect feature '{f}' chosen. Select from: {acceptable_features}.")
        
        if not features:
            raise ValueError(f"Rename ids was called but no feature levels were chosen. Select from: {acceptable_features}.")
        
        if prefix:
            ignored_options = []
            if keep_existing_ids_if_derived_from_base_id:
                ignored_options.append("keep_existing_ids_if_derived_from_base_id")
            if set(features) != set(acceptable_features):
                ignored_options.append("features")
            if remove_point_suffix:
                ignored_options.append("remove_point_suffix")
            if strip_gene_tag:
                ignored_options.append("strip_gene_tag")
            if ignored_options:
                warnings.warn(f"Providing a prefix '{prefix}' means all features will be renamed based on the prefix, the following provided options are to be ignored: {ignored_options}.", category=UserWarning)

        elif suffix:
            warnings.warn(f"Provided suffix={suffix} will have no effect as no prefix was provided and custom renaming will therefore be skipped.", category=UserWarning)

        if cds_segment_ids and "CDS" not in features:
            warnings.warn("CDS features will be changed if need be since cds_segment_ids have been requested.", category=UserWarning)

        start_time = time.time()

        changed_features = set()

        # Gene ids are numbered by position along each chromosome
        self.sort_genes(processes=1)

        progress_bar = start_progress_bar(total=len(self.all_gene_ids), description=f"Renaming {self.id} gene ids", colour="38;2;156;42;42", quiet=quiet)

        correspondences_d = {}

        # Organelle genes are labelled chrMT/chrCP instead of their chromosome name; chromosomes sharing
        # a label (e.g. several mitochondrial contigs) continue its gene numbering
        organelle_labels = {"mitochondria": "chrMT", "chloroplast": "chrCP"}
        g_counts: dict[str, int] = {}

        for chrom, genes in self.chrs.items():
            chrom_label = organelle_labels.get(self.chromosome_compartment(chrom), chrom)
            for g in genes.values():
                progress_bar.update(1)
                g_counts[chrom_label] = g_counts.get(chrom_label, 0) + spacer
                g_count = g_counts[chrom_label]

                if "gene" in features or prefix:
                    g.rename(g_count, sep=sep, digits=g_id_digits, prefix=prefix, suffix=suffix, base_id_as_id=strip_gene_tag, remove_point_suffix=remove_point_suffix, chrom_label=chrom_label)
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
                    
                # Id formats only decide what is renamed when existing ids may be kept
                if keep_existing_ids_if_derived_from_base_id and base_id_present and base_id_missing:
                    warnings.warn(f"{self.id} annotation gene {g.original_id} has a mix of transcript id formats and renaming errors could occur!", category=UserWarning)
                
                if tmains > 1:
                    raise ValueError(f"There should not be more than one main transcript for {g.original_id} in {self.id} annotation.")


                t_count = 1
                keep_transcript_numbers = keep_subfeature_numbers and unambiguous_numbers((t.id, t.id_number) for t in g.transcripts.values())

                for t in g.transcripts.values():

                    t.parents = [g.id]

                    if not t.main:
                        t_count += 1

                    if "transcript" in features or prefix:

                        if prefix or base_id_missing:
                            t.rename(base_id=g.base_id, sep=sep, count=t_count, digits=t_id_digits, keep_numbering=keep_transcript_numbers, keep_existing_ids_if_derived_from_base_id=False)
                        else:
                            t.rename(base_id=g.base_id, sep=sep, count=t_count, digits=t_id_digits, keep_numbering=keep_transcript_numbers, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id)

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

                    if keep_existing_ids_if_derived_from_base_id and base_id_present and base_id_missing:
                        warnings.warn(f"{self.id} annotation transcript {t.original_id} has a mix of CDS id formats and renaming errors could occur!", category=UserWarning)

                    if cmains > 1:
                        raise ValueError(f"There shouldn't be more than one main CDS for transcript {t.original_id} in {self.id} annotation.")

                    c_count = 1
                    for c in t.CDSs.values():
                        c.parents = [t.id]

                        if not c.main:
                            c_count += 1

                        if "CDS" in features or prefix:

                            if prefix or base_id_missing:
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

                    if keep_existing_ids_if_derived_from_base_id and base_id_present and base_id_missing:
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

                    if keep_existing_ids_if_derived_from_base_id and base_id_present and base_id_missing:
                        warnings.warn(f"{self.id} annotation transcript {t.original_id} has a mix of UTR id formats and renaming errors could occur!", category=UserWarning)


                if "exon" in features or prefix:
                    g.rename_exons(sep=sep, digits=t_id_digits, keep_numbering=keep_subfeature_numbers, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id)

                    if g.renamed_exons:
                        changed_features.add("exon")

                if "UTR" in features or prefix:
                    g.rename_utrs(sep=sep, digits=t_id_digits, keep_numbering=keep_subfeature_numbers, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id)

                    if g.renamed_utrs:
                        changed_features.add("UTR")


        progress_bar.close()

        self.homogenise_shared_exons_utrs()
        self.update_keys()
        self.update_gene_and_transcript_list(quiet=quiet)
        if old_id_names == "alias":
            self._move_stale_names_to_aliases()

        self.renamed_features = changed_features

        if set(features) == {"gene", "transcript", "CDS", "exon", "UTR"}:
            self.feature_tags.add("full_renamed_ids")

        if correspondences:
        
            if "gene" in changed_features:

                extra_suffixes = ["rename_eqs"]

                final_correspondences_path = self._resolve_output_path(output_dir=correspondences_output_dir, extension=correspondences_extension, extra_suffixes=extra_suffixes, suffix=self.feature_suffix, filename=correspondences_filename, use_annot_dir=correspondences_use_annot_dir, subfolder=correspondences_subfolder, filepath=correspondences_filepath)

                with open(str(final_correspondences_path), "w", encoding="utf-8") as f_out:
                    f_out.write("old_gene_id\tnew_gene_id\n")
                    for k, v in correspondences_d.items():
                        f_out.write(f"{k}\t{v}\n")

            else:
                warnings.warn(f"Correspondences on gene ids were requested, however gene ids remained unchanged.", category=UserWarning)

        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nRenaming {self.id} ids with prefix='{prefix}', changing={self.renamed_features} features took {round(lapse/60, 1)} minutes")        

    def _move_stale_names_to_aliases(self):
        """
        A renamed feature whose Name is its old ID (e.g. Araport 'ID=AT4G00005;Name=AT4G00005') would
        be labelled with an ID it no longer has: that Name is moved to Alias, the GFF3 attribute for
        secondary identifiers.
        """
        def move(f):
            if f.names and f.id != f.original_id and f.original_id in f.names:
                f.names = [n for n in f.names if n != f.original_id] or None
                if f.aliases is None:
                    f.aliases = []
                    f.place_attribute_like(INT_ALIAS, INT_NAME)
                if f.original_id not in f.aliases:
                    f.aliases.append(f.original_id)

        for f in self._iter_all_features():
            move(f)

    def standardise_attribute_keys(self):
        """Writes Name, Alias and symbol under their default keys instead of the spelling read (e.g. 'gene_synonym', 'Symbol')."""
        for f in self._iter_all_features():
            if f.attribute_keys:
                f.attribute_keys = None

    def reorder_attributes(self):
        """Writes Name, Alias and symbol first and in that order (AEGIS's default), instead of where they were in column 9."""
        for f in self._iter_all_features():
            if f.attribute_order:
                f.attribute_order = None

    def _iter_all_features(self):
        for genes in self.chrs.values():
            for g in genes.values():
                yield g
                for t in g.transcripts.values():
                    yield t
                    yield from t.exons
                    yield from t.miRNAs
                    for c in t.CDSs.values():
                        yield c
                        yield from c.CDS_segments
                        yield from c.UTRs
        yield from self.atypical_features
        yield from self.orphaned_features

    def update_keys(self, gene_keys:bool=True, transcript_keys:bool=True, CDS_keys:bool=True):

        for chrom, genes_dict in self.chrs.items():
            new_genes = {}
            
            for gene in genes_dict.values():
                new_transcripts = {}
                for transcript in gene.transcripts.values():
                    new_cdss = {cds.id: cds for cds in transcript.CDSs.values()}
                    if CDS_keys:
                        transcript.CDSs = new_cdss
                    new_transcripts[transcript.id] = transcript
                
                if transcript_keys:
                    gene.transcripts = new_transcripts

                new_genes[gene.id] = gene

            if gene_keys:
                self.chrs[chrom] = new_genes


    def create_featurecounts_ids(self):
        for genes in self.chrs.values():
            for g in genes.values():
                g.gene_id = g.id
                for t in g.transcripts.values():
                    t.gene_id = g.id
                    for c in t.CDSs.values():
                        for cs in c.CDS_segments:
                            cs.gene_id = g.id
                        for u in c.UTRs:
                            u.gene_id = g.id
                    for e in t.exons:
                        e.gene_id = g.id

    def create_gtf_attributes(self):
        """
        GTF attributes for every gene and subfeature; values are quoted with gtf_quote. As in Ensembl
        GTFs, gene_name is on every line of a gene and transcript_name on every line of a transcript,
        as tools that read names from exon lines (e.g. Cell Ranger) expect. Names are symbol_values(),
        else the Names as read.
        """
        for genes in self.chrs.values():
            for g in genes.values():

                gene_id = f'gene_id {gtf_quote(g.id)}'
                gene_names = g.symbol_values() or g.names
                gene_name = [f'gene_name {gtf_quote(",".join(gene_names))}'] if gene_names else []
                g.gtf_attributes = [gene_id] + gene_name

                for t in g.transcripts.values():
                    transcript_names = t.symbol_values() or t.names
                    transcript_name = [f'transcript_name {gtf_quote(",".join(transcript_names))}'] if transcript_names else []
                    ids = [gene_id, f'transcript_id {gtf_quote(t.id)}']
                    names = gene_name + transcript_name
                    t.gtf_attributes = ids + names
                    num_exons = len(t.exons)
                    for x, e in enumerate(t.exons):
                        exon_num = (num_exons - x) if t.strand == "-" else (x + 1)
                        e.gtf_attributes = ids + [f'exon_number "{exon_num}"'] + names

                    for c in t.CDSs.values():
                        c.gtf_attributes = ids + names
                        for cs in c.CDS_segments:
                            cs.gtf_attributes = ids + names
                        for u in c.UTRs:
                            u.gtf_attributes = ids + names

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
                if not self.protein_equivalences:
                    self.generate_protein_equivalences(quiet=True)
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

        self.gff_header = new_header

    def subset(
        self,
        chosen_features: set[str] | list[str] | tuple[str] | None = None,
        gene_cap: int | None = 3000,
        common_chromosomes: set | list | tuple | None = None,
        min_genes: int | None = None,
        quiet: bool = False,
        chr_cap: int | None = None,
        no_gene_cap: bool = False,
        no_min_genes: bool = False,
        seed: int | None = None,
    ):

        if seed is not None:
            random.seed(seed)

        if chosen_features is None:
            chosen_features = set()
        elif isinstance(chosen_features, (list, tuple)):
            chosen_features = set(chosen_features)
        else:
            chosen_features = chosen_features.copy()

        initial_chosen_features = chosen_features.copy()

        for chosen_feature in chosen_features:
            if chosen_feature not in self.chrs:
                raise ValueError(f"Chosen scaffold/chromosome {chosen_feature} is not in {self.name} genome. Choose from '{self.chrs.keys()}'")

        if common_chromosomes is None:
            total_chromosomes = set(self.chrs)
        elif isinstance(common_chromosomes, (list, tuple)):
            total_chromosomes = set(common_chromosomes)
        else:
            total_chromosomes = common_chromosomes.copy()

        if chr_cap is not None and chr_cap > 0 and not chosen_features:
            if chr_cap <= len(total_chromosomes):
                chosen_features = set(random.sample(list(total_chromosomes), chr_cap))
            else:
                chosen_features = total_chromosomes.copy()

        if not initial_chosen_features:
            if no_min_genes or (min_genes is not None and min_genes <= 0):
                effective_min_genes = 0
            elif min_genes is None:
                effective_min_genes = 1500
            else:
                effective_min_genes = min_genes
        else:
            if no_min_genes or min_genes is None or min_genes <= 0:
                effective_min_genes = 0
            else:
                effective_min_genes = min_genes

        if effective_min_genes > 0:

            num_genes_in_chosen_features = 0
            for ft in chosen_features:
                num_genes_in_chosen_features += len(self.chrs[ft])

            remaining_to_chose_from = total_chromosomes - chosen_features
            chr_cap_overriden = False

            while num_genes_in_chosen_features < effective_min_genes and remaining_to_chose_from:
                chr_cap_overriden = True
                
                chosen_features.add(random.choice(list(remaining_to_chose_from)))

                remaining_to_chose_from = total_chromosomes - chosen_features

                num_genes_in_chosen_features = 0
                for ft in chosen_features:
                    num_genes_in_chosen_features += len(self.chrs[ft])

            if chr_cap_overriden and not quiet:
                print(f"Chromosome/scaffold cap of {len(initial_chosen_features)} was overriden by min_genes = {effective_min_genes} parameter as not enough genes were present in {initial_chosen_features}. The final selection includes {len(chosen_features)} features: {chosen_features}")

        if not chosen_features:
            chosen_features = total_chromosomes.copy()

        features_to_remove = set(self.chrs) - chosen_features
        self.remove_chromosomes(features_to_remove, update=False, quiet=quiet)

        is_gene_cap_enabled = not no_gene_cap and gene_cap is not None and gene_cap > 0

        if is_gene_cap_enabled:
            genes_to_keep_per_chromosome = math.ceil(gene_cap / len(chosen_features)) if len(chosen_features) > 0 else 0 #type: ignore

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

        self.update(quiet=quiet, update_gene_and_transcript_list=True, remove_missing_transcript_parent_references=True)

        return chosen_features

    def filter_by_rna_class(self, rna_classes=['mRNA'], remove_genes_accordingly:bool=True, quiet:bool=False):

        transcript_to_remove = set()

        for genes in self.chrs.values():

            for g in genes.values():

                for t in g.transcripts.values():

                    if t.feature not in rna_classes:

                        transcript_to_remove.add(t.id)

        self.remove_transcripts(transcript_to_remove, remove_genes_accordingly=remove_genes_accordingly, quiet=quiet)
    
        self.update(quiet=quiet)

    def remove_chromosomes(self, features_to_remove:set, update:bool=True, quiet:bool=False):
        if features_to_remove:
            for ft in features_to_remove:
                del self.chrs[ft]
            self.remove_chromosomes_from_header()

        if update:
            self.update(quiet=quiet, update_gene_and_transcript_list=True, remove_missing_transcript_parent_references=True)

    def remove_genes(self, to_remove:set|None=None, override_rescue:bool=False, quiet:bool=False):

        if to_remove is None:
            to_remove = set()

        total_count = 0

        for genes in self.chrs.values():
            total_count += len(genes)

        progress_bar = start_progress_bar(total=total_count, description=f"Removing {self.id} genes", colour="91", quiet=quiet)
        
        for gene in to_remove:
            if gene in self.all_gene_ids:
                chrom = self.all_gene_ids[gene]
                self.chrs[chrom][gene].quality.remove = True
                if override_rescue:
                    self.chrs[chrom][gene].quality.rescue = False
            else:
                warnings.warn(f"Gene {gene} is not present in annotation {self.id}.", category=UserWarning)

        genes_to_remove = set()
        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                if g.quality.remove and not g.quality.rescue:
                    genes_to_remove.add(g.id)

        for g_id in genes_to_remove:
            chrom = self.all_gene_ids[g_id]
            del self.chrs[chrom][g_id]
            del self.all_gene_ids[g_id]
        progress_bar.close()
        
        self.remove_missing_genes_in_overlaps(quiet=quiet)
        self.update_gene_and_transcript_list(quiet=quiet)

    def remove_missing_genes_in_overlaps(self, quiet:bool=True):

        total_count = 0
        for genes in self.chrs.values():
            total_count += len(genes)

        progress_bar = start_progress_bar(total=total_count, description=f"Removing {self.id} missing genes in overlaps", colour="91", quiet=quiet)

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

        progress_bar = start_progress_bar(total=len(self.all_gene_ids), description=f"Removing {self.id} duplicate transcripts", colour="91", quiet=quiet)

        for genes in self.chrs.values():
            for g in genes.values():
                progress_bar.update(1)
                to_eliminate = set()
                transcripts = list(g.transcripts.values())
                for i in range(len(transcripts)):
                    t1 = transcripts[i]
                    if t1.id in to_eliminate:
                        continue
                    for j in range(i + 1, len(transcripts)):
                        t2 = transcripts[j]
                        if t2.id in to_eliminate:
                            continue
                        if t1.almost_equal(t2):
                            to_eliminate.add(t2.id)

                for t_eliminate in to_eliminate:
                    del g.transcripts[t_eliminate]
        progress_bar.close()
        self.update(rename_features=("transcript", "CDS", "exon", "UTR"), quiet=quiet)

    def make_alternative_genes_into_transcripts(self, quiet:bool=False):

        correspondence = {}

        to_remove = []
        for chrom, genes in self.chrs.items():
            for g in genes.values():
                if (not g.quality.remove or g.quality.rescue) and (g.id not in correspondence):
                    if g.overlaps:
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

        self.update(rename_features=("gene", "transcript", "CDS", "exon", "UTR"), quiet=quiet)

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
                    g.quality.rescue = False
                    g.quality.remove = True
        self.remove_genes(quiet=quiet)

        if removed_any:
            self.feature_tags.add("minus_small_CDS")
        self.update(quiet=quiet)

    def remove_TE_genes(self, quiet:bool=False):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if g.transposable:
                    g.quality.rescue = False
                    g.quality.remove = True
                    removed_any = True

        self.remove_genes(quiet=quiet)

        if removed_any:
            self.feature_tags.add("minus_TE")
        self.update(quiet=quiet)

    def remove_non_TE_genes(self, quiet:bool=False):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if not g.transposable:
                    g.quality.rescue = False
                    g.quality.remove = True
                    removed_any = True

        self.remove_genes(quiet=quiet)
        
        if removed_any:
            self.feature_tags.add("minus_non_TE")

        self.update(quiet=quiet)

    def remove_pseudogene_genes(self, quiet:bool=False):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if g.pseudogene:
                    g.quality.rescue = False
                    g.quality.remove = True
                    removed_any = True

        self.remove_genes(quiet=quiet)

        if removed_any:
            self.feature_tags.add("minus_pseudogene")
        self.update(quiet=quiet)

    def remove_non_pseudogene_genes(self, quiet:bool=False):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if not g.pseudogene:
                    g.quality.rescue = False
                    g.quality.remove = True
                    removed_any = True

        self.remove_genes(quiet=quiet)
        
        if removed_any:
            self.feature_tags.add("minus_non_pseudogene")

        self.update(quiet=quiet)

    def remove_non_coding_genes_and_transcripts(self, quiet:bool=False):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if g.noncoding == True and g.coding == False:
                    g.quality.rescue = False
                    g.quality.remove = True
                    removed_any = True

        self.remove_genes(quiet=quiet)

        removed_any = self.remove_non_coding_transcripts_from_coding_genes(removed_any, False, quiet=quiet)

        if removed_any:
            self.feature_tags.add("minus_non_coding")

        self.update(quiet=quiet)

    def remove_coding_genes_and_transcripts(self, quiet:bool=False):

        removed_any = False

        for genes in self.chrs.values():
            for g in genes.values():
                if g.noncoding == False and g.coding == True:
                    g.quality.rescue = False
                    g.quality.remove = True
                    removed_any = True

        self.remove_genes(quiet=quiet)

        removed_any = self.remove_coding_transcripts_from_non_coding_genes(removed_any, False, quiet=quiet)

        if removed_any:
            self.feature_tags.add("minus_coding")

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
            self.feature_tags.add("minus_coding")
        
        self.overlaps.clear()

        if update:
            self.update(rename_features=("transcript", "CDS", "exon", "UTR"), quiet=quiet)

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
            self.feature_tags.add("minus_non_coding")

        self.overlaps.clear()

        if update:
            self.update(rename_features=("transcript", "CDS", "exon", "UTR"), quiet=quiet)

        return removed_any

    def clear_gene_names_and_symbols(self):
        for genes in self.chrs.values():
            for g in genes.values():
                g.names = None
                g.symbols = None
        self.tags.discard("plus_symbols")

    def remove_genes_without_symbols(self, quiet:bool=False):
        """Removes genes with no symbol: no Name other than an ID and no Symbol (see Feature.symbol_values)."""
        for genes in self.chrs.values():
            for g in genes.values():
                if not g.symbol_values():
                    g.quality.remove = True
        self.remove_genes(quiet=quiet)
        self.update(quiet=quiet)

    def rename_chromosomes(self, equivalences, dap:bool=False):

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
                    if t.promoter is not None:
                        t.promoter.ch = chrom
                    for c in t.CDSs.values():
                        c.ch = chrom
                        for cs in c.CDS_segments:
                            cs.ch = chrom
                        if c.protein is not None:
                            c.protein.ch = chrom
                        for u in c.UTRs:
                            u.ch = chrom
                    for e in t.exons:
                        e.ch = chrom
                    if t.introns is not None:
                        for i in t.introns:
                            i.ch = chrom

        for o in self.atypical_features:
            if o.ch in equivalences:
                o.ch = equivalences[o.ch]

        for o in self.orphaned_features:
            if o.ch in equivalences:
                o.ch = equivalences[o.ch]
        if dap:
            self.tags.add("dapfit")
            if renamed_scaffolds:
                self.tags.add("dapmod")
                self.tags.discard("dapfit")
        elif renamed_scaffolds:
            self.tags.add("confrenamed")

    _GENE_ATTRIBUTE_TARGETS = {"name": "Name", "alias": "Alias", "symbol": "Symbol"}
    _GENE_ATTRIBUTE_PROPERTIES = {"Name": "names", "Alias": "aliases", "Symbol": "symbols"}

    @staticmethod
    def _attribute_values(raw) -> list[str]:
        """A string is one value (never split on ','); empty and NaN values are dropped; order is kept."""
        if raw is None:
            return []
        items = [raw] if isinstance(raw, str) else list(raw)
        values = []
        for v in items:
            if v is None or (isinstance(v, float) and math.isnan(v)):
                continue
            v = str(v).strip()
            if v and v not in values:
                values.append(v)
        return values

    _MATCH_FIELDS = ("id", "name", "alias")

    @classmethod
    def _target_attributes(cls, targets) -> list[str]:
        """'Name', 'Alias' and/or 'Symbol' (case-insensitive) as attribute names, without repeats."""
        if isinstance(targets, str):
            targets = (targets,)
        attributes = []
        for t in targets:
            attribute = cls._GENE_ATTRIBUTE_TARGETS.get(str(t).strip().lower())
            if attribute is None:
                raise ValueError(f"Unknown target '{t}'. Choose from: Name, Alias, Symbol.")
            if attribute not in attributes:
                attributes.append(attribute)
        if not attributes:
            raise ValueError("No targets chosen. Choose from: Name, Alias, Symbol.")
        return attributes

    @classmethod
    def _match_fields(cls, match_on) -> tuple[str, ...]:
        if isinstance(match_on, str):
            match_on = (match_on,)
        fields = tuple(str(m).strip().lower() for m in match_on)
        invalid = [m for m in fields if m not in cls._MATCH_FIELDS]
        if invalid or not fields:
            raise ValueError(f"Invalid match_on {invalid or list(fields)}. Choose from: {', '.join(cls._MATCH_FIELDS)}.")
        return fields

    @staticmethod
    def _check_choice(name, value, choices):
        if value not in choices:
            raise ValueError(f"{name} must be one of {', '.join(repr(c) for c in choices)}, not '{value}'.")

    def _gene_finder(self, match_on:tuple[str, ...], id_type:str):
        """Function from a key to the genes it matches: by the first of match_on fields that matches any."""
        indexes = []
        for field in match_on:
            index: dict[str, list[Gene]] = {}
            for genes in self.chrs.values():
                for g in genes.values():
                    if field == "id":
                        keys = (g.original_id if id_type == "original" else g.id,)
                    else:
                        keys = (g.names if field == "name" else g.aliases) or ()
                    for k in keys:
                        bucket = index.setdefault(k, [])
                        if not bucket or bucket[-1] is not g:
                            bucket.append(g)
            indexes.append(index)

        def find(key) -> list[Gene]:
            for index in indexes:
                genes = index.get(key)
                if genes:
                    return genes
            return []
        return find

    def _transcript_finder(self, id_type:str):
        if id_type == "original":
            by_original: dict = {}
            for genes in self.chrs.values():
                for g in genes.values():
                    for t in g.transcripts.values():
                        by_original.setdefault(t.original_id, t)
            return by_original.get

        def find(transcript_id):
            location = self.all_transcript_ids.get(transcript_id)
            if location is None:
                return None
            g = self.chrs.get(location[0], {}).get(location[1])
            return g.transcripts.get(transcript_id) if g is not None else None
        return find

    def _resolve_ids(self, mapping:dict, id_type:str, transcript_names:bool, match_on:tuple[str, ...]):
        """
        Mapping keys to features: ({id(gene): (gene, values)}, {id(transcript): (transcript, values)},
        ids not found, ids that matched several genes).
        """
        find_genes = self._gene_finder(match_on, id_type)
        find_transcript = self._transcript_finder(id_type) if transcript_names else None

        gene_values: dict[int, tuple[Gene, list[str]]] = {}
        transcript_values: dict[int, tuple[Transcript, list[str]]] = {}
        missing, shared = [], []
        for key, raw in mapping.items():
            key = str(key).strip()
            values = self._attribute_values(raw)
            genes = find_genes(key)
            if genes:
                if len(genes) > 1:
                    shared.append(key)
                found = [(gene_values, g) for g in genes]
            else:
                t = find_transcript(key) if find_transcript is not None else None
                if t is None:
                    missing.append(key)
                    continue
                found = [(transcript_values, t)]
            for store, feature in found:
                current = store.setdefault(id(feature), (feature, []))[1]
                for v in values:
                    if v not in current:
                        current.append(v)
        return gene_values, transcript_values, missing, shared

    def _report_ids(self, missing:list[str], shared:list[str], transcript_names:bool, on_missing:str, match_on:tuple[str, ...]):
        def examples(ids):
            return f"{', '.join(ids[:5])}{', ...' if len(ids) > 5 else ''}"
        if shared:
            warnings.warn(f"{len(shared)} ids matched several {self.id} genes by {'/'.join(f for f in match_on if f != 'id')} and were set on all of them: {examples(shared)}", category=UserWarning)
        if not missing or on_missing == "ignore":
            return
        message = f"{len(missing)} {'gene or transcript' if transcript_names else 'gene'} ids not found in {self.id} annotation: {examples(missing)}"
        if not transcript_names:
            transcript_ids = sum(1 for m in missing if m in self.all_transcript_ids)
            if transcript_ids:
                message += f" ({transcript_ids} of them are transcript ids, which are only used with transcript_names / --transcript-names)"
        if on_missing == "raise":
            raise KeyError(message)
        warnings.warn(message, category=UserWarning)

    def _apply_attribute_values(self, attributes:list[str], gene_values:dict, transcript_values:dict, mode:str, clear_unmapped:bool, transcript_names:bool):
        all_genes = [g for genes in self.chrs.values() for g in genes.values()]
        all_transcripts = [t for g in all_genes for t in g.transcripts.values()]
        gene_keys = {a: self._most_used_gff_key(all_genes, a) for a in attributes}
        transcript_keys = {a: self._most_used_gff_key(all_transcripts, a) or gene_keys[a] for a in attributes}

        def assign(feature, values, keys):
            for attribute in attributes:
                prop = self._GENE_ATTRIBUTE_PROPERTIES[attribute]
                current = getattr(feature, prop) or []
                updated = values[:] if mode == "replace" else current + [v for v in values if v not in current]
                if not current and updated:
                    feature.set_gff_key(attribute, keys[attribute])
                setattr(feature, prop, updated or None)

        def clear(feature):
            for attribute in attributes:
                setattr(feature, self._GENE_ATTRIBUTE_PROPERTIES[attribute], None)

        for g in all_genes:
            gene_entry = gene_values.get(id(g))
            if gene_entry is not None:
                assign(g, gene_entry[1], gene_keys)
            elif clear_unmapped:
                clear(g)
            for t in g.transcripts.values():
                transcript_entry = transcript_values.get(id(t))
                if transcript_entry is not None:
                    assign(t, transcript_entry[1], transcript_keys)
                elif gene_entry is not None:
                    if transcript_names:
                        assign(t, gene_entry[1], transcript_keys)
                    elif mode == "replace":
                        clear(t)
                elif clear_unmapped:
                    clear(t)

    def set_gene_attribute(self, mapping:dict, targets:str|tuple[str, ...]|list[str]=("Name",), mode:Literal["replace", "append"]="replace", clear_unmapped:bool=False, id_type:Literal["current", "original"]="current", on_missing:Literal["warn", "raise", "ignore"]="warn", transcript_names:bool=False, match_on:str|tuple[str, ...]|list[str]=("id",)) -> list[str]:
        """
        Sets gene names, aliases and/or symbols from a {gene id: value(s)} mapping.

        mapping: values are a string (one value, never split on ',') or a list of strings.
        targets: "Name" (labelled and indexed by genome browsers such as JBrowse), "Alias" and/or "Symbol",
            case-insensitive; every target gets the same values.
        mode: 'replace' the target values of mapped genes, or 'append' the new ones to them.
        clear_unmapped: also empty the targets of genes not in the mapping. By default they keep their values.
        id_type: mapping keys are 'current' ids, or 'original' ids as read (before rename_ids).
        on_missing: what to do with ids not in the annotation: 'warn', 'raise' a KeyError (before
            anything is changed) or 'ignore'.
        transcript_names: transcripts of mapped genes get their gene's values, and mapping keys may
            also be transcript ids, which set that transcript's values instead (not its gene's).
        match_on: what mapping keys are compared with to find genes, tried in order: 'id' (see
            id_type), 'name' and/or 'alias' (e.g. old ids moved to Alias by rename_ids). A key that
            matches several genes by name or alias sets all of them, with a warning.

        Transcripts follow their gene: without transcript_names, 'replace' empties the targets of the
        transcripts of mapped genes (they would still describe the gene's old values), and
        clear_unmapped those of the transcripts of unmapped genes too; 'append' leaves them as read.

        Values a feature did not have before are written with the key spelling most used in the
        annotation for that attribute (e.g. 'gene_synonym' for aliases), else Name, Alias or symbol.
        Returns the mapping ids not found.
        """
        attributes = self._target_attributes(targets)
        self._check_choice("mode", mode, ("replace", "append"))
        self._check_choice("id_type", id_type, ("current", "original"))
        self._check_choice("on_missing", on_missing, ("warn", "raise", "ignore"))
        match_on = self._match_fields(match_on)

        gene_values, transcript_values, missing, shared = self._resolve_ids(mapping, id_type, transcript_names, match_on)
        self._report_ids(missing, shared, transcript_names, on_missing, match_on)
        self._apply_attribute_values(attributes, gene_values, transcript_values, mode, clear_unmapped, transcript_names)
        return missing

    def _most_used_gff_key(self, features, attribute) -> str | None:
        """Key spelling most used for attribute among features that have values for it (None if none do)."""
        prop = self._GENE_ATTRIBUTE_PROPERTIES[attribute]
        counts: dict[str, int] = {}
        for f in features:
            if getattr(f, prop):
                key = f.gff_key(attribute)
                counts[key] = counts.get(key, 0) + 1
        return max(counts, key=lambda k: counts[k]) if counts else None

    @classmethod
    def read_symbols_table(cls, file_path:str, targets:str|tuple[str, ...]|list[str]=("Name",), id_column:int|str=1, header:bool=False, sep:str|None=None, split:str|None=None) -> list[tuple[list[str], dict[str, list[str]]]]:
        """
        Reads a TSV/CSV/XLSX table of ids and symbols for set_symbol_columns: one
        (attributes, {id: values}) pair per value column.

        targets: "Name", "Alias" and/or "Symbol", each optionally followed by ':column' naming the
            column with its values (default 2), e.g. ("Name:2", "Alias:3"). Targets sharing a column
            get the same values.
        id_column: column with the ids (default 1).
            Columns are 1-based numbers or, with header, column names.
        sep: column separator; None uses ',' for .csv files and tab otherwise (ignored for .xlsx/.xls).
        split: also split each cell on this separator (e.g. ';'); by default a cell is one value
            (never split on ',').
        Cells are read as text (e.g. 'NA' is a symbol, not a missing value). Rows with an empty id or
        value map nothing, so an empty cell does not clear a gene (see clear_unmapped). Repeated ids
        give a feature several values.
        """
        if isinstance(targets, str):
            targets = (targets,)
        columns: dict[str, list[str]] = {}
        for t in targets:
            target, _, column = str(t).partition(":")
            attribute = cls._target_attributes(target)[0]
            column = column.strip() or "2"
            if attribute not in columns.setdefault(column, []):
                columns[column].append(attribute)
        if not columns:
            raise ValueError("No targets chosen. Choose from: Name, Alias, Symbol.")

        lower = file_path.lower()
        if lower.endswith((".xlsx", ".xls")):
            try:
                df = pd.read_excel(file_path, header=0 if header else None, dtype=str, keep_default_na=False)
            except ImportError as e:
                package = "xlrd" if "xlrd" in str(e) else "openpyxl"
                raise ValueError(f"Reading {os.path.basename(file_path)} needs the '{package}' package (pip install {package}), or save the table as TSV/CSV.") from None
        else:
            if sep is None:
                sep = "," if lower.removesuffix(".gz").endswith(".csv") else "\t"
            df = pd.read_csv(file_path, header=0 if header else None, sep=sep, dtype=str, keep_default_na=False)
        df = df.fillna("")

        def position(column) -> int:
            column = str(column).strip()
            if column.isdigit() and 1 <= int(column) <= df.shape[1]:
                return int(column) - 1
            if header and column in [str(c) for c in df.columns]:
                return [str(c) for c in df.columns].index(column)
            available = ", ".join(str(c) for c in df.columns) if header else f"1-{df.shape[1]}"
            raise ValueError(f"Column '{column}' not found in {file_path} (columns: {available}; separator {sep!r}).")

        ids = df.iloc[:, position(id_column)]
        result = []
        for column, attributes in columns.items():
            mapping: dict[str, list[str]] = {}
            for feature_id, cell in zip(ids, df.iloc[:, position(column)]):
                feature_id = str(feature_id).strip()
                values = [v.strip() for v in (str(cell).split(split) if split else [str(cell)])]
                values = [v for v in values if v]
                if feature_id and values:
                    mapping.setdefault(feature_id, []).extend(values)
            result.append((attributes, mapping))
        return result

    def set_symbol_columns(self, columns:list[tuple[list[str], dict[str, list[str]]]], mode:Literal["replace", "append"]="replace", clear_unmapped:bool=False, id_type:Literal["current", "original"]="current", on_missing:Literal["warn", "raise", "ignore"]="warn", transcript_names:bool=False, match_on:str|tuple[str, ...]|list[str]=("id",)) -> list[str]:
        """
        Sets each (attributes, {id: values}) pair from read_symbols_table with set_gene_attribute (see
        it for the other arguments). Ids are checked for all columns first, so on_missing='raise'
        changes nothing. Returns the ids not found.
        """
        self._check_choice("mode", mode, ("replace", "append"))
        self._check_choice("id_type", id_type, ("current", "original"))
        self._check_choice("on_missing", on_missing, ("warn", "raise", "ignore"))
        match_on = self._match_fields(match_on)

        resolved = []
        missing, shared = [], []
        for attributes, mapping in columns:
            gene_values, transcript_values, column_missing, column_shared = self._resolve_ids(mapping, id_type, transcript_names, match_on)
            resolved.append((self._target_attributes(attributes), gene_values, transcript_values))
            missing.extend(m for m in column_missing if m not in missing)
            shared.extend(m for m in column_shared if m not in shared)
        self._report_ids(missing, shared, transcript_names, on_missing, match_on)

        for attributes, gene_values, transcript_values in resolved:
            self._apply_attribute_values(attributes, gene_values, transcript_values, mode, clear_unmapped, transcript_names)
        self.tags.add("plus_symbols")
        return missing

    def add_gene_symbols(self, file_path:str, targets:str|tuple[str, ...]|list[str]=("Name",), mode:Literal["replace", "append"]="replace", clear_unmapped:bool=False, header:bool=False, sep:str|None=None, id_type:Literal["current", "original"]="current", on_missing:Literal["warn", "raise", "ignore"]="warn", transcript_names:bool=False, id_column:int|str=1, split:str|None=None, match_on:str|tuple[str, ...]|list[str]=("id",)) -> list[str]:
        """
        Reads a table of ids and symbols (read_symbols_table: targets, id_column, header, sep, split)
        and sets them (set_symbol_columns / set_gene_attribute: the other arguments). Returns the ids
        not found.
        """
        columns = self.read_symbols_table(file_path, targets=targets, id_column=id_column, header=header, sep=sep, split=split)
        return self.set_symbol_columns(columns, mode=mode, clear_unmapped=clear_unmapped, id_type=id_type, on_missing=on_missing, transcript_names=transcript_names, match_on=match_on)

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
        self.rename_ids(prefix=id_prefix, spacer=spacer, suffix=suffix, features=("gene", "transcript", "CDS", "exon", "UTR"), quiet=quiet)
        self.update(quiet=quiet)
        self.export.gff(output_dir=custom_path, extra_attributes=extra_attributes, filename=tag, skip_atypical_fts=skip_atypical_fts, main_only=main_only, UTRs=UTRs, quiet=quiet, clean_attributes=True)

    def generate_introns(self):
        for genes in self.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    t.generate_introns()

    def __str__(self):
        return str(self.id)
