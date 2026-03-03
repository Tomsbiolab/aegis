from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from ..annotation import Annotation
    from ..genome import Genome

import time
import os
import warnings
import pandas as pd
import sys

from pathlib import Path
from tqdm import tqdm
from typing import Literal

class AnnotationExport:
    """
    Component for handling export methods for the Annotation class.
    Accessed via 'annotation_object.export'.
    """
    def __init__(self, annotation: Annotation):
        self._annot = annotation

    def all_features(self, feature_output: Literal["main", "all", "both"] = "main", promoters: bool = True, verbose: bool = True, path: str = "", most_specific_id_level = "promoter", quiet: bool = False):
        """
        The "output" parameter can be both, main or all. This parameter only 
        affects promoter, transcript, CDS and protein sequences. If both is selected
        a "main features" file and and "all features" file will be produced.
        """

        feature_output_choices = ["main", "all", "both"]
        if feature_output not in feature_output_choices:
            raise ValueError(f"feature_output={feature_output} is not amongst the feature_output_choices={feature_output_choices} to export all features.")

        most_specific_id_level_choices = ["gene", "transcript", "CDS", "protein", "promoter"]

        if most_specific_id_level not in most_specific_id_level_choices:
            raise ValueError(f"most_specific_id_level={most_specific_id_level} is not amongst the most_specific_id_level_choices={most_specific_id_level_choices} to export proteins.")

        start_time = time.time()

        self.genes(verbose, path)

        if feature_output == "both":
            modes = [True, False]
        elif feature_output == "main":
            modes = [True]
        else:
            modes = [False]

        for b in modes:

            if most_specific_id_level == "promoter" or most_specific_id_level == "protein":
                self.proteins(b, verbose, path)
                self.CDSs(b, verbose, path)
                self.transcripts(b, verbose, path)
                if promoters:
                    self.promoters(b, verbose, path)

            elif most_specific_id_level == "CDS":
                self.proteins(b, verbose, path, used_id="CDS")
                self.CDSs(b, verbose, path)
                self.transcripts(b, verbose, path)
                if promoters:
                    self.promoters(b, verbose, path)

            elif most_specific_id_level == "transcript":
                self.proteins(b, verbose, path, used_id="transcript")
                self.CDSs(b, verbose, path, used_id="transcript")
                self.transcripts(b, verbose, path)
                if promoters:
                    self.promoters(b, verbose, path, used_id="transcript")

            elif most_specific_id_level == "gene":
                self.proteins(b, verbose, path, used_id="gene")
                self.CDSs(b, verbose, path, used_id="gene")
                self.transcripts(b, verbose, path, used_id="gene")
                if promoters:
                    self.promoters(b, verbose, path, used_id="gene")

        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"Extracting {self._annot.id} annotation features took {round(lapse, 1)} seconds\n")

    def proteins(self, only_main: bool = True, verbose: bool = True, custom_path: str = "", used_id: str = "protein", unique_proteins_per_gene: bool = False, only_cds_main: bool = True):
        """
        Main proteins means only proteins obtained from the main CDSs of the
        main transcripts. This equates to one protein per gene.

        Verbose will include protein tag details and readthrough options, strand, chromosome, and coordinates.

        Proteins will be exported into fastas with their protein ids unless 
        CDS, transcript or gene has been selected as used_id. To avoid duplicate
        equal entries, when choosing gene or transcript only main proteins will
        be able to be output.

        valid_id_choices = ["gene", "transcript", "CDS", "protein"]
        """

        if custom_path:
            output_file = Path(custom_path)
        else:
            output_file = Path(self._annot.path) / "features"
        output_file.mkdir(parents=True, exist_ok=True)
        output_file = str(output_file) + "/"

        output_file += self._annot.id
        output_file += self._annot.feature_suffix
        output_file += "_proteins"

        out = ""

        if unique_proteins_per_gene:
            only_main = False
            only_cds_main = False
            if used_id == "gene":
                used_id = "protein"
                warnings.warn(f"Used id 'gene' has been changed to 'protein' as unique proteins per gene was selected.", category=UserWarning)

        if used_id == "gene":
            only_main = True
            only_cds_main = True
            output_file += "_g_id_main"
        else:
            if used_id == "transcript":
                only_cds_main = True
                output_file += "_t_id"
                if unique_proteins_per_gene:
                    warnings.warn(f"If more than one CDS exists per transcript (this is rarely the case), CDSs beyond the main CDS will not be considered, since 'transcript' was the used_id. Select 'CDS' or 'protein' if all CDSs are to be considered.", category=UserWarning)
            elif used_id == "CDS":
                output_file += "_c_id"
            elif used_id == "protein":
                output_file += "_p_id"

            if unique_proteins_per_gene:
                output_file += "_unique_per_gene"
            elif only_main:
                output_file += "_main"
            else:
                output_file += "_all"

        if verbose:
            output_file += "_coordinates"
            
        output_file += ".fasta"
        
        valid_id_choices = ["gene", "transcript", "CDS", "protein"]
        if used_id not in valid_id_choices:
            raise ValueError(f"used_id={used_id} is not amongst the valid_id_choices={valid_id_choices} to export proteins.")
        
        for genes in self._annot.chrs.values():
            for g in genes.values():
                temp_cs = []

                for t in g.transcripts.values():
                    if only_main:
                        if t.main:
                            for c in t.CDSs.values():
                                if c.seq != "":
                                    if only_cds_main:
                                        if c.main:
                                            temp_cs.append(c)
                                    else:
                                        temp_cs.append(c)
                    else:
                        for c in t.CDSs.values():
                            if c.seq != "":
                                if only_cds_main:
                                    if c.main:
                                        temp_cs.append(c)
                                else:
                                    temp_cs.append(c)

                if not unique_proteins_per_gene:
                    final_cs = temp_cs.copy()

                else:
                    final_cs = []
                    if len(temp_cs) > 0:
                        final_cs.append(temp_cs[0])

                    for i, c1 in enumerate(temp_cs):
                        if i > 0:
                            add = True
                            for c2 in final_cs:
                                if c1.equal_segments(c2):
                                    add = False
                            if add:
                                final_cs.append(c1)

                for c in final_cs:

                    if used_id == "protein":
                        out += f">{c.protein.id}"
                    elif used_id == "CDS":
                        out += f">{c.id}"
                    elif used_id == "transcript":
                        out += f">{c.parents[0]}"
                    elif used_id == "gene":
                        out += f">{g.id}"

                    if c.protein.summary_tag and verbose:
                        out += f"|{c.protein.summary_tag}"
                    if verbose:
                        out += f"|readthrough:{c.protein.readthrough}|{c.strand}|{c.ch}|{c.start}:{c.end}"

                    out += f"\n{c.protein.seq}\n"        
                     
        if out != "":
            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()
        else:
            print(f"Warning: Run self.generate_sequences(genome) on {self._annot.id}")

    def unique_proteins(self, genome: Genome | None = None, custom_path: str = "", quiet: bool = False):
        start_time = time.time()
        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False
        
        if not self._annot.contains_protein_sequences:
            if genome == None:
                print(f"You forgot to provide the genome for {self._annot.id} and no sequences exist for CDSs")
            else:
                self._annot.generate_sequences(genome, quiet=quiet)
                self.unique_proteins(custom_path=custom_path)
        else:
            if custom_path:
                output_file = Path(custom_path)
            else:
                output_file = Path(self._annot.path) / "features"
            output_file.mkdir(parents=True, exist_ok=True)
            output_file = str(output_file) + "/"
            output_file += f"{self._annot.id}{self._annot.feature_suffix}_unique_proteins.fasta"
            out = ""
            all_protein_seqs = {}
            self._annot.all_protein_ids = {}
            for chrom, genes in self._annot.chrs.items():
                for g in genes.values():
                    if g.coding:
                        for t in g.transcripts.values():
                            for c in t.CDSs.values():
                                if c.protein is not None:
                                    if c.protein.seq != "":
                                        all_protein_seqs[c.protein.id] = c.protein.seq
                                        self._annot.all_protein_ids[c.protein.id] = (chrom, g.id, t.id, c.id)

            progress_bar = tqdm(total=len(all_protein_seqs.keys()), disable=disable,
                            bar_format=(
                f'\033[1;91mDetermining and exporting unique {self._annot.id} proteins:\033[0m '
                '{percentage:3.0f}%|'
                f'\033[1;91m{{bar}}\033[0m| '
                '{n}/{total} [{elapsed}<{remaining}]'))
            
            
            unique_sequences = {}
            self._annot.protein_equivalences = {}

            for protein_id, sequence in all_protein_seqs.items():
                progress_bar.update(1)
                if sequence not in unique_sequences:
                    unique_sequences[sequence] = protein_id
                    self._annot.protein_equivalences[protein_id] = []
                else:
                    first_protein_id = unique_sequences[sequence]
                    self._annot.protein_equivalences[first_protein_id].append(protein_id)

            for sequence, protein_id in unique_sequences.items():
                out += f">{protein_id}\n{sequence}\n"

            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()

            progress_bar.close()

            now = time.time()
            lapse = now - start_time
            if not quiet:
                print(f"\nExporting unique {self._annot.id} proteins took {round(lapse/60, 1)} minutes")

    def CDSs(self, only_main: bool = True, verbose: bool = True, custom_path: str = "", used_id: str = "CDS", unique_CDSs_per_gene: bool = False, only_cds_main: bool = True):
        """
        Main CDSs means only CDS sequence obtained from the main CDS of the
        main transcripts.

        Verbose will include strand, chromosome, and coordinates.

        CDSs will be exported into fastas with their CDS ids unless 
        transcript or gene has been selected as used_id. To avoid duplicate
        equal entries, when choosing gene or transcript only main CDSs will
        be able to be output.

        valid_id_choices = ["gene", "transcript", "CDS"]
        """

        if custom_path:
            output_file = Path(custom_path)
        else:
            output_file = Path(self._annot.path) / "features"
        output_file.mkdir(parents=True, exist_ok=True)
        output_file = str(output_file) + "/"

        output_file += self._annot.id
        output_file += self._annot.feature_suffix
        output_file += "_CDSs"

        out = ""

        if unique_CDSs_per_gene:
            only_main = False
            only_cds_main = False
            if used_id == "gene":
                used_id = "CDS"
                warnings.warn(f"Used id 'gene' has been changed to 'CDS' as unique CDSs per gene was selected.", category=UserWarning)

        if used_id == "gene":
            only_main = True
            only_cds_main = True
            output_file += "_g_id_main"
        else:
            if used_id == "transcript":
                only_cds_main = True
                output_file += "_t_id"
                if unique_CDSs_per_gene:
                    warnings.warn(f"If more than one CDS exists per transcript (this is rarely the case), CDSs beyond the main CDS will not be considered, since 'transcript' was the used_id. Select 'CDS' if all CDSs are to be considered.", category=UserWarning)
            elif used_id == "CDS":
                output_file += "_c_id"

            if unique_CDSs_per_gene:
                output_file += "_unique_per_gene"
            elif only_main:
                output_file += "_main"
            else:
                output_file += "_all"

        if verbose:
            output_file += "_coordinates"
            
        output_file += ".fasta"
        
        valid_id_choices = ["gene", "transcript", "CDS"]
        if used_id not in valid_id_choices:
            raise ValueError(f"used_id={used_id} is not amongst the valid_id_choices={valid_id_choices} to export CDSs.")

        for genes in self._annot.chrs.values():
            for g in genes.values():
                temp_cs = []
                for t in g.transcripts.values():
                    if only_main:
                        if t.main:
                            for c in t.CDSs.values():
                                if c.seq != "":
                                    if only_cds_main:
                                        if c.main:
                                            temp_cs.append(c)
                                    else:
                                        temp_cs.append(c)
                    else:
                        for c in t.CDSs.values():
                            if c.seq != "":
                                if only_cds_main:
                                    if c.main:
                                        temp_cs.append(c)
                                else:
                                    temp_cs.append(c)

                if not unique_CDSs_per_gene:
                    final_cs = temp_cs.copy()
                else:
                    final_cs = []

                    if len(temp_cs) > 0:
                        final_cs.append(temp_cs[0])

                    for i, c1 in enumerate(temp_cs):
                        if i > 0:
                            add = True
                            for c2 in final_cs:
                                if c1.equal_segments(c2):
                                    add = False
                            if add:
                                final_cs.append(c1)

                for c in final_cs:

                    if used_id == "CDS":
                        out += f">{c.id}"
                    elif used_id == "transcript":
                        out += f">{c.parents[0]}"
                    elif used_id == "gene":
                        out += f">{g.id}"

                    if verbose:
                        out += f"|{c.strand}|{c.ch}|{c.start}:{c.end}"

                    out += f"\n{c.seq}\n"
                     
        if out != "":
            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()
        else:
            print(f"Warning: Run self.generate_sequences(genome) on {self._annot.id}")

    def transcripts(self, only_main: bool = True, verbose: bool = True, custom_path: str = "", used_id: str = "transcript", rna_classes: list = []):
        """
        Main means only main transcript sequences are exported.

        Verbose will include strand, chromosome, and coordinates.

        Transcripts will be exported into fastas with their transcript ids unless 
        gene has been selected as used_id. To avoid duplicate
        equal entries, when choosing gene only main transcripts will
        be able to be output.

        valid_id_choices = ["gene", "transcript"]
        """
        if custom_path:
            output_file = Path(custom_path)
        else:
            output_file = Path(self._annot.path) / "features"
        output_file.mkdir(parents=True, exist_ok=True)
        output_file = str(output_file) + "/"
        output_file += self._annot.id
        output_file += self._annot.feature_suffix
        output_file += "_transcripts"

        out = ""

        if used_id == "gene":
            only_main = True
            output_file += "_g_id_main"
        elif used_id == "transcript":
            output_file += "_t_id"
            if only_main:
                output_file += "_main"
            else:
                output_file += "_all"

        if verbose:
            output_file += "_coordinates"
            
        output_file += ".fasta"

        valid_id_choices = ["gene", "transcript"]
        if used_id not in valid_id_choices:
            raise ValueError(f"used_id={used_id} is not amongst the valid_id_choices={valid_id_choices} to export transcripts.")
        
        for genes in self._annot.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    if (t.feature in rna_classes) or (not rna_classes):
                        if t.seq != "":
                            if only_main and not t.main:
                                continue

                            if used_id == "transcript":    
                                out += f">{t.id}"
                            elif used_id == "gene":
                                out += f">{g.id}"

                            if verbose:
                                out += f"|{t.strand}|{t.ch}|{t.start}:{t.end}"
                            out += f"\n{t.seq}\n"

        if out != "":
            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()
        else:
            print(f"Warning: Run self.generate_sequences(genome) on {self._annot.id}")

    def genes(self, verbose: bool = True, custom_path: str = ""):
        if custom_path:
            output_file = Path(custom_path)
        else:
            output_file = Path(self._annot.path) / "features"
        output_file.mkdir(parents=True, exist_ok=True)
        output_file = str(output_file) + "/"
        output_file += self._annot.id
        output_file += self._annot.feature_suffix
        output_file += "_genes"


        if verbose:
            output_file += "_coordinates"
        output_file += ".fasta"
        out = ""

        for genes in self._annot.chrs.values():
            for g in genes.values():
                if g.seq != "":
                    out += f">{g.id}"
                    if verbose:
                        out += f"|{g.strand}|{g.ch}|{g.start}:{g.end}"
                    out += f"\n{g.seq}\n"
        if out != "":
            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()
        else:
            print(f"Warning: Run self.generate_sequences(genome) on {self._annot.id}")

    def promoters(self, only_main: bool = True, verbose: bool = True, custom_path: str = "", used_id: str = "promoter"):
        """

        Verbose will include promoter type, strand, chromosome, and coordinates.

        """
        if custom_path:
            output_file = Path(custom_path)
        else:
            output_file = Path(self._annot.path) / "features"
        output_file.mkdir(parents=True, exist_ok=True)
        output_file = str(output_file) + "/"
        output_file += self._annot.id
        output_file += self._annot.feature_suffix
        output_file += "_promoters"

        out = ""

        if used_id == "gene":
            only_main = True
            output_file += "_g_id_main"
        else:
            if used_id == "transcript":
                output_file += "_t_id"
            if used_id == "promoter":
                output_file += "_prom_id"
            if only_main:
                output_file += "_main"
            else:
                output_file += "_all"

        if verbose:
            output_file += "_coordinates"
            
        output_file += f"_{self._annot.promoter_size}_{self._annot.promoter_types}.fasta"

        valid_id_choices = ["gene", "transcript", "promoter"]

        if used_id not in valid_id_choices:
            raise ValueError(f"used_id={used_id} is not amongst the valid_id_choices={valid_id_choices} to export promoters.")

        for genes in self._annot.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    if only_main and not t.main:
                        continue

                    if hasattr(t, "promoter"):
                        if t.promoter.seq != "":

                            if used_id == "promoter":
                                out += f">{t.promoter.id}"
                            elif used_id == "transcript":
                                out += f">{t.id}"
                            elif used_id == "gene":
                                out += f">{g.id}"

                            if verbose:
                                out += f"|{t.promoter.type}|{t.strand}|{t.ch}|{t.promoter.start}:{t.promoter.end}"

                            out += f"\n{t.promoter.seq}\n"
                        else:
                            print(f"Warning: transcript {t.id} from annotation {self._annot.id} has no promoter sequence.\n")
                    else:
                        print(f"Warning: transcript {t.id} from annotation {self._annot.id} has no promoter.\n")

        if out != "":
            f_out = open(output_file, "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()
        else:
            print(f"Warning: Run self.generate_promoters(genome) and self.generate_sequences(genome) on {self._annot.id}")

    def for_dapseq(self, genome: Genome, genome_out_folder: str = "", gff_out_folder: str = "", tag: str = "_for_dap.gff3", skip_atypical_fts: bool = True, main_only: bool = False, UTRs: bool = False, exclude_non_coding: bool = False):
        equivalences = genome.rename_features_dap(output_folder=genome_out_folder, return_equivalences=True, export=True)
        self._annot.rename_chromosomes(equivalences)
        self.gff(custom_path=gff_out_folder, tag=tag, skip_atypical_fts=skip_atypical_fts, main_only=main_only, UTRs=UTRs, just_genes=exclude_non_coding)

    def gff(self, custom_path: str = "", tag: str = ".gff3", skip_atypical_fts: bool = False, main_only: bool = False, UTRs: bool = False, just_genes: bool = False, no_1bp_features: bool = False, repeat_exons_utrs: bool = False, subfolder: bool = True, quiet: bool = False, skip_orphaned_fts: bool = False):

        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False

        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[38;2;46;204;113m\033[1mExporting gff for {self._annot.id} annotation:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[38;2;46;204;113m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

        if subfolder:
            export_folder = Path(custom_path or self._annot.path) / "out_gffs"
        else:
            export_folder = Path(custom_path or self._annot.path)
        export_folder.mkdir(parents=True, exist_ok=True)
        export_folder = str(export_folder) + "/"

        output_suffix = ""

        if just_genes:
            output_suffix += "_just_genes"
        if no_1bp_features:
            output_suffix += "_for_lifton"

        if tag == ".gff3":
            tag = f"{self._annot.id}{self._annot.suffix}{output_suffix}{tag}"
            if not quiet:
                print(f"Exporting {self._annot.id} gff with tag='{tag}' which is dapfit={self._annot.dapfit} and dapmod={self._annot.dapmod} and combined={self._annot.combined}.")
        elif not quiet:
            print(f"Exporting {self._annot.id} gff to {export_folder}{tag}.")

        with open(f"{export_folder}{tag}", "w", encoding="utf-8") as f_out:
            if self._annot.clean or self._annot.dapmod:
                f_out.write("##gff-version 3\n")
            else:
                f_out.write("\n".join(self._annot.gff_header) + "\n")

            if repeat_exons_utrs:
                self._annot.single_parent_for_exons_utrs(quiet=quiet)

            for x1, genes in enumerate(self._annot.chrs.values()):
                progress_bar.update(len(genes))
                for x2, g in enumerate(genes.values()):
                    
                    if no_1bp_features:
                        gene_1bp_feature = False
                        for t in g.transcripts.values():
                            for e in t.exons:
                                if e.size == 1:
                                    gene_1bp_feature = True
                            for c in t.CDSs.values():
                                for cs in c.CDS_segments:
                                    if cs.size == 1:
                                        gene_1bp_feature = True
                                for u in c.UTRs:
                                    if u.size == 1:
                                        gene_1bp_feature = True

                        if gene_1bp_feature:
                            continue

                    f_out.write(g.print_gff())

                    if just_genes:
                        continue

                    if repeat_exons_utrs or main_only or len(g.transcripts) == 1:
                        for t in g.transcripts.values():
                            if main_only:
                                if not t.main:
                                    continue
                            f_out.write(t.print_gff())
                            for e in t.exons:
                                f_out.write(e.print_gff())

                            for m in t.miRNAs:
                                f_out.write(m.print_gff())

                            for c in t.CDSs.values():
                                if main_only:
                                    if not c.main:
                                        continue
                                for c_seg in c.CDS_segments:
                                    f_out.write(c_seg.print_gff())

                                if UTRs:
                                    if hasattr(c, "UTRs"):
                                        for u in c.UTRs:
                                            f_out.write(u.print_gff())

                    else:
                                
                        for t in g.transcripts.values():                                
                            f_out.write(t.print_gff())

                        exons = []
                        for t in g.transcripts.values():
                            exons.extend(t.exons)

                        exons.sort()
                        
                        unique_exons = []
                        if exons:
                            unique_exons.append(exons[0])
                            for i in range(1, len(exons)):
                                if not exons[i].equal_sequence(exons[i-1]):
                                    unique_exons.append(exons[i])

                        for e in unique_exons:
                            f_out.write(e.print_gff())

                        for m in t.miRNAs:
                            f_out.write(m.print_gff())

                        for t in g.transcripts.values():
                            for c in t.CDSs.values():
                                for c_seg in c.CDS_segments:
                                    f_out.write(c_seg.print_gff())


                        if UTRs:
                            utrs = []
                            for t in g.transcripts.values():
                                for c in t.CDSs.values():
                                    if hasattr(c, "UTRs"):
                                        utrs.extend(c.UTRs)
                            utrs.sort()

                            unique_utrs = []
                            if utrs:
                                unique_utrs.append(utrs[0])
                                for i in range(1, len(utrs)):
                                    if not utrs[i].equal_sequence(utrs[i-1]):
                                        unique_utrs.append(utrs[i])

                            for u in unique_utrs:
                                f_out.write(u.print_gff())

                    if x1 == (len(self._annot.chrs) - 1) and x2 == (len(genes) - 1):
                        continue
                    f_out.write("###\n")

            progress_bar.close()

            if not skip_atypical_fts:
                if not just_genes:
                    for x, ft in enumerate(self._annot.atypical_features):
                        if x == 0:
                            f_out.write("###\n")
                        f_out.write(ft.print_gff())
                        if x == (len(self._annot.atypical_features) - 1):
                            continue
                        f_out.write("###\n")
            if not skip_orphaned_fts:
                if not just_genes:
                    for x, ft in enumerate(self._annot.orphaned_features):
                        if x == 0:
                            f_out.write("###\n")
                        f_out.write(ft.print_gff())
                        if x == (len(self._annot.orphaned_features) - 1):
                            continue
                        f_out.write("###\n")

    def gtf(self, custom_path: str = "", tag: str = ".gtf", main_only: bool = False, UTRs: bool = False, just_genes: bool = False, no_1bp_features: bool = False, quiet: bool = False):

        self._annot.create_gtf_attributes()

        # Check if stdout or stderr are redirected to files
        stdout_redirected = not sys.stdout.isatty()
        stderr_redirected = not sys.stderr.isatty()

        # Disable tqdm if stdout or stderr are redirected
        if stdout_redirected or stderr_redirected or quiet:
            disable = True
        else:
            disable = False

        progress_bar = tqdm(total=len(self._annot.all_gene_ids.keys()), disable=disable,
                                bar_format=(
                    f'\033[38;2;46;204;113m\033[1mExporting gtf for {self._annot.id} annotation:\033[0m '
                    '{percentage:3.0f}%|'
                    f'\033[38;2;46;204;113m{{bar}}\033[0m| '
                    '{n}/{total} [{elapsed}<{remaining}]'))

        export_folder = Path(custom_path or self._annot.path) / "out_gtfs"
        export_folder.mkdir(parents=True, exist_ok=True)
        export_folder = str(export_folder) + "/"

        output_suffix = ""

        if just_genes:
            output_suffix += "_just_genes"
        if no_1bp_features:
            output_suffix += "_for_lifton"

        if tag == ".gtf":
            tag = f"{self._annot.id}{self._annot.suffix}{output_suffix}{tag}"
            if not quiet:
                print(f"Exporting {self._annot.id} gtf with tag='{tag}' which is dapfit={self._annot.dapfit} and dapmod={self._annot.dapmod} and combined={self._annot.combined}.")
        elif not quiet:
            print(f"Exporting {self._annot.id} gtf to {export_folder}{tag}.")

        with open(f"{export_folder}{tag}", "w", encoding="utf-8") as f_out:
            f_out.write("#gtf-version 2.2\n")

            for x1, genes in enumerate(self._annot.chrs.values()):
                progress_bar.update(len(genes))
                for x2, g in enumerate(genes.values()):

                    if no_1bp_features:
                        gene_1bp_feature = False
                        for t in g.transcripts.values():
                            for e in t.exons:
                                if e.size == 1:
                                    gene_1bp_feature = True
                            for c in t.CDSs.values():
                                for cs in c.CDS_segments:
                                    if cs.size == 1:
                                        gene_1bp_feature = True
                                for u in c.UTRs:
                                    if u.size == 1:
                                        gene_1bp_feature = True

                        if gene_1bp_feature:
                            continue

                    f_out.write(g.print_gtf())

                    if just_genes:
                        continue

                    for t in g.transcripts.values():
                        if main_only:
                            if not t.main:
                                continue
                        original_feature = t.feature
                        t.feature = "transcript"
                        f_out.write(t.print_gtf())
                        t.feature = original_feature
                        for e in t.exons:
                            f_out.write(e.print_gtf())
                        for c in t.CDSs.values():
                            if main_only:
                                if not c.main:
                                    continue
                            for c_seg in c.CDS_segments:
                                f_out.write(c_seg.print_gtf())
                            if UTRs:
                                if hasattr(c, "UTRs"):
                                    for u in c.UTRs:
                                        f_out.write(u.print_gtf())

                    if x1 == (len(self._annot.chrs) - 1) and x2 == (len(genes) - 1):
                        continue
                    f_out.write("###\n")

            progress_bar.close()

    def gene_list(self, custom_path: str = "", output_file: str = "", lengths: bool = False, coordinates: bool = False, chromosomes: bool = False, coding_info: bool = False, skip_coding: bool = False, skip_non_coding: bool = False, sep: str = "\t", skip_pseudogenes: bool = False, skip_transposables: bool = False, gene_symbols: bool = False):

        if not custom_path:
            export_folder = Path(self._annot.path) / "lists"
        else:
            export_folder = Path(custom_path)
        export_folder.mkdir(parents=True, exist_ok=True)
        export_folder = str(export_folder) + "/"

        if not output_file:
            output_file = f"{export_folder}{self._annot.name}_genes.txt"
        else:
            output_file = f"{export_folder}{output_file}"

        if skip_coding and skip_non_coding:
            print("Warning: Both skip_coding and skip_non_coding are set to True. No genes will be listed.")
            return

        with open(output_file, "w", encoding="utf-8") as f_out:

            header = ["gene_id"]
            if chromosomes or coordinates:
                header.append("scaffold")
            if coordinates:
                header.append("gene_start")
                header.append("gene_end")
            if lengths:
                header.append("gene_length")
            if coding_info:
                header.append("coding")
            if gene_symbols:
                header.append("gene_symbol")
            f_out.write(sep.join(header) + "\n")

            for chrom, genes in self._annot.chrs.items():
                for g in genes.values():
                    if skip_pseudogenes and g.pseudogene:
                        continue
                    if skip_transposables and g.transposable:
                        continue
                    if skip_coding and g.coding:
                        continue
                    if skip_non_coding and not g.coding:
                        continue

                    out = [g.id]
                    if chromosomes or coordinates:
                        out.append(chrom)
                    if coordinates:
                        out.append(str(g.start))
                        out.append(str(g.end))
                    if lengths:
                        out.append(str(g.size))
                    if coding_info:
                        out.append(str(g.coding))
                    if gene_symbols:
                        out.append("|".join(g.symbols))
                    f_out.write(sep.join(out) + "\n")

    def transcript_list(self, custom_path: str = "", output_file: str = "", lengths: bool = False, coordinates: bool = False, chromosomes: bool = False, coding_info: bool = False, skip_coding: bool = False, skip_non_coding: bool = False, sep: str = "\t", skip_pseudogenes: bool = False, skip_transposables: bool = False, gene_symbols: bool = False):
        if not custom_path:
            export_folder = Path(self._annot.path) / "lists"
        else:
            export_folder = Path(custom_path)
        export_folder.mkdir(parents=True, exist_ok=True)
        export_folder = str(export_folder) + "/"

        if not output_file:
            output_file = f"{export_folder}{self._annot.name}_transcripts.txt"
        else:
            output_file = f"{export_folder}{output_file}"

        if skip_coding and skip_non_coding:
            print("Warning: Both skip_coding and skip_non_coding are set to True. No transcripts will be listed.")
            return

        with open(output_file, "w", encoding="utf-8") as f_out:

            header = ["transcript_id"]
            if chromosomes or coordinates:
                header.append("scaffold")
            if coordinates:
                header.append("transcript_start")
                header.append("transcript_end")
            if lengths:
                header.append("transcript_length")
            if coding_info:
                header.append("coding")
            if gene_symbols:
                header.append("gene_symbol")
            f_out.write(sep.join(header) + "\n")

            for chrom, genes in self._annot.chrs.items():
                for g in genes.values():
                    if skip_pseudogenes and g.pseudogene:
                        continue
                    if skip_transposables and g.transposable:
                        continue

                    for t in g.transcripts.values():
                        if skip_coding and t.coding:
                            continue
                        if skip_non_coding and not t.coding:
                            continue

                        out = [t.id]
                        if chromosomes or coordinates:
                            out.append(chrom)
                        if coordinates:
                            out.append(str(t.start))
                            out.append(str(t.end))
                        if lengths:
                            out.append(str(t.size))
                        if coding_info:
                            out.append(str(t.coding))
                        if gene_symbols:
                            out.append("|".join(g.symbols))
                        f_out.write(sep.join(out) + "\n")


