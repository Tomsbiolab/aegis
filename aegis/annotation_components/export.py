from __future__ import annotations
from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from ..genome import Genome
    from ..subfeatures import CDS

import time
import warnings

from typing import Literal

from ..utils.misc import start_progress_bar
from .base import AnnotationComponent

class AnnotationExport(AnnotationComponent):
    """
    Component for handling export methods for the Annotation class.
    Accessed via 'annotation_object.export'.
    """

    def all_features(self, feature_output: Literal["main", "all", "both"] = "main", promoters: bool = True, verbose: bool = True, most_specific_id_level = "promoter", output_dir: str | None = None, use_annot_dir:bool=False, subfolder:bool=False, subfolder_name:str="features", extension:str=".fasta", quiet: bool = False,
        #deprecated arguments
        path:str=""):
        """
        The "output" parameter can be both, main or all. This parameter only 
        affects promoter, transcript, CDS and protein sequences. If both is selected
        a "main features" file and and "all features" file will be produced.
        """

        if path != "":
            warnings.warn("'path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = path

        if subfolder_name != "features":
            subfolder = True

        feature_output_choices = ["main", "all", "both"]
        if feature_output not in feature_output_choices:
            raise ValueError(f"feature_output={feature_output} is not amongst the feature_output_choices={feature_output_choices} to export all features.")

        most_specific_id_level_choices = ["gene", "transcript", "CDS", "protein", "promoter"]

        if most_specific_id_level not in most_specific_id_level_choices:
            raise ValueError(f"most_specific_id_level={most_specific_id_level} is not amongst the most_specific_id_level_choices={most_specific_id_level_choices} to export proteins.")

        start_time = time.time()

        self.genes(verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir)

        if feature_output == "both":
            modes = [True, False]
        elif feature_output == "main":
            modes = [True]
        else:
            modes = [False]

        for b in modes:

            if most_specific_id_level == "promoter" or most_specific_id_level == "protein":
                self.proteins(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir)
                self.CDSs(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir)
                self.transcripts(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir)
                if promoters:
                    self.promoters(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir)

            elif most_specific_id_level == "CDS":
                self.proteins(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir, used_id="CDS")
                self.CDSs(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir)
                self.transcripts(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir)
                if promoters:
                    self.promoters(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir)

            elif most_specific_id_level == "transcript":
                self.proteins(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir, used_id="transcript")
                self.CDSs(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir, used_id="transcript")
                self.transcripts(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir)
                if promoters:
                    self.promoters(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir, used_id="transcript")

            elif most_specific_id_level == "gene":
                self.proteins(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir, used_id="gene")
                self.CDSs(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir, used_id="gene")
                self.transcripts(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir, used_id="gene")
                if promoters:
                    self.promoters(only_main=b, verbose=verbose, output_dir=output_dir, subfolder=subfolder, subfolder_name=subfolder_name, extension=extension, use_annot_dir=use_annot_dir, used_id="gene")

        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"Extracting {self._annot.id} annotation features took {round(lapse, 1)} seconds\n")

    def proteins(self, only_main: bool = True, verbose: bool = True, used_id: str = "protein", unique_proteins_per_gene: bool = False, only_cds_main: bool = True, mode: Literal["start", "end", "orf", "orf_or_end"] = "end", use_name_not_id: bool = False, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "features", extension=".fasta",
        #deprecated arguments
        custom_filename: str="", custom_path:str=""):

        if custom_filename != "":
            warnings.warn("'custom_filename' is deprecated. Please use 'filename' instead.", DeprecationWarning, stacklevel=2)
            if filename is None:
                filename = custom_filename

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        if subfolder_name != "features":
            subfolder = True

        final_cs: list[CDS]

        valid_id_choices = ["gene", "transcript", "CDS", "protein"]
        if used_id not in valid_id_choices:
            raise ValueError(f"used_id={used_id} is not amongst the valid_id_choices={valid_id_choices} to export proteins.")

        if not self._annot.contains_protein_sequences:
            self._annot.generate_proteins(mode=mode)

        extra_suffixes = ["proteins"]

        if unique_proteins_per_gene:
            only_main = False
            only_cds_main = False
            if used_id == "gene":
                used_id = "protein"
                warnings.warn(f"Used id 'gene' has been changed to 'protein' as unique proteins per gene was selected.", category=UserWarning)

        if used_id == "gene":
            only_main = True
            only_cds_main = True
            extra_suffixes.append("g_id_main")
        else:
            if used_id == "transcript":
                only_cds_main = True
                extra_suffixes.append("t_id")
                if unique_proteins_per_gene:
                    warnings.warn(f"If more than one CDS exists per transcript (this is rarely the case), CDSs beyond the main CDS will not be considered, since 'transcript' was the used_id. Select 'CDS' or 'protein' if all CDSs are to be considered.", category=UserWarning)
            elif used_id == "CDS":
                extra_suffixes.append("c_id")
            elif used_id == "protein":
                extra_suffixes.append("p_id")

            if unique_proteins_per_gene:
                extra_suffixes.append("unique_per_gene")
            elif only_main:
                extra_suffixes.append("main")
            else:
                extra_suffixes.append("all")

        if verbose:
            extra_suffixes.append("coordinates")

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self._annot.feature_suffix, extension=extension, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder, extra_suffixes=extra_suffixes, use_name_not_id=use_name_not_id)
        
        with open(str(final_output_path), "w", encoding="utf-8") as f_out:
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
                            f_out.write(f">{c.protein.id}")
                        elif used_id == "CDS":
                            f_out.write(f">{c.id}")
                        elif used_id == "transcript":
                            f_out.write(f">{c.parents[0]}")
                        elif used_id == "gene":
                            f_out.write(f">{g.id}")

                        if c.protein.summary_tag and verbose:
                            f_out.write(f"|{c.protein.summary_tag}")
                        if verbose:
                            f_out.write(f"|readthrough:{c.protein.readthrough}|{c.strand}|{c.protein.ch}|{c.protein.start}:{c.protein.end}")

                        f_out.write(f"\n{c.protein.seq}\n")

    def unique_proteins(self, use_name_not_id: bool = False, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "features", extension=".fasta", quiet: bool = False, mode: Literal["start", "end", "orf", "orf_or_end"] = "end",
        #deprecated arguments
        custom_path:str=""):

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        start_time = time.time()

        extra_suffixes = ["unique_proteins"]

        if subfolder_name != "features":
            subfolder = True

        self._annot.generate_protein_equivalences(mode=mode, quiet=quiet)

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self._annot.feature_suffix, extension=extension, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder, extra_suffixes=extra_suffixes, use_name_not_id=use_name_not_id)

        with open(str(final_output_path), "w", encoding="utf-8") as f_out:
            for protein_id in self._annot.protein_equivalences:
                chrom, g, t, c = self._annot.all_protein_ids[protein_id]
                sequence = self._annot.chrs[chrom][g].transcripts[t].CDSs[c].protein.seq
                f_out.write(f">{protein_id}\n{sequence}\n")

        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nExporting unique {self._annot.id} proteins took {round(lapse/60, 1)} minutes")

    def unique_transcripts(self, use_name_not_id: bool = False, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "features", extension=".fasta", quiet: bool = False, rna_classes: list = [],
        #deprecated arguments
        custom_path:str=""):

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        if subfolder_name != "features":
            subfolder = True

        start_time = time.time()

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self._annot.feature_suffix, extension=extension, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder, extra_suffixes=["unique_transcripts"], use_name_not_id=use_name_not_id)
        
        all_transcript_seqs = {}
        for genes in self._annot.chrs.values():
            for g in genes.values():
                for t in g.transcripts.values():
                    if (t.feature in rna_classes) or (not rna_classes):
                        if t.seq != "":
                            all_transcript_seqs[t.id] = t.seq

        progress_bar = start_progress_bar(total=len(all_transcript_seqs.keys()), description=f"Exporting unique {self._annot.id} transcripts", quiet=quiet, colour="91")
        
        unique_sequences = {}

        for transcript_id, sequence in all_transcript_seqs.items():
            progress_bar.update(1)
            if sequence not in unique_sequences:
                unique_sequences[sequence] = transcript_id

        with open(str(final_output_path), "w", encoding="utf-8") as f_out:
            for sequence, transcript_id in unique_sequences.items():
                f_out.write(f">{transcript_id}\n{sequence}\n")

        progress_bar.close()

        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nExporting unique {self._annot.id} transcripts took {round(lapse/60, 1)} minutes")

    def unique_CDSs(self, use_name_not_id: bool = False, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "features", extension=".fasta", quiet: bool = False,
        #deprecated arguments
        custom_path:str=""):

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        if subfolder_name != "features":
            subfolder = True
        
        start_time = time.time()

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self._annot.feature_suffix, extension=extension, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder, extra_suffixes=["unique_CDSs"], use_name_not_id=use_name_not_id)

        all_CDS_seqs = {}
        for genes in self._annot.chrs.values():
            for g in genes.values():
                if g.coding:
                    for t in g.transcripts.values():
                        for c in t.CDSs.values():
                            if c.seq != "":
                                all_CDS_seqs[c.id] = c.seq

        progress_bar = start_progress_bar(total=len(all_CDS_seqs.keys()), description=f"Exporting unique {self._annot.id} CDSs", quiet=quiet, colour="91")
        
        unique_sequences = {}

        for CDS_id, sequence in all_CDS_seqs.items():
            progress_bar.update(1)
            if sequence not in unique_sequences:
                unique_sequences[sequence] = CDS_id

        with open(str(final_output_path), "w", encoding="utf-8") as f_out:
            for sequence, CDS_id in unique_sequences.items():
                f_out.write(f">{CDS_id}\n{sequence}\n")

        progress_bar.close()

        now = time.time()
        lapse = now - start_time
        if not quiet:
            print(f"\nExporting unique {self._annot.id} CDSs took {round(lapse/60, 1)} minutes")

    def CDSs(self, only_main: bool = True, verbose: bool = True, used_id: str = "CDS", unique_CDSs_per_gene: bool = False, only_cds_main: bool = True, use_name_not_id: bool = False, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "features", extension=".fasta",
        #deprecated arguments
        custom_filename: str="", custom_path:str=""):
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

        if custom_filename != "":
            warnings.warn("'custom_filename' is deprecated. Please use 'filename' instead.", DeprecationWarning, stacklevel=2)
            if filename is None:
                filename = custom_filename

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        valid_id_choices = ["gene", "transcript", "CDS"]
        if used_id not in valid_id_choices:
            raise ValueError(f"used_id={used_id} is not amongst the valid_id_choices={valid_id_choices} to export CDSs.")

        if subfolder_name != "features":
            subfolder = True

        extra_suffixes = ["CDSs"]

        if unique_CDSs_per_gene:
            only_main = False
            only_cds_main = False
            if used_id == "gene":
                used_id = "CDS"
                warnings.warn(f"Used id 'gene' has been changed to 'CDS' as unique CDSs per gene was selected.", category=UserWarning)

        if used_id == "gene":
            only_main = True
            only_cds_main = True
            extra_suffixes.append("g_id_main")
        else:
            if used_id == "transcript":
                only_cds_main = True
                extra_suffixes.append("t_id")
                if unique_CDSs_per_gene:
                    warnings.warn(f"If more than one CDS exists per transcript (this is rarely the case), CDSs beyond the main CDS will not be considered, since 'transcript' was the used_id. Select 'CDS' if all CDSs are to be considered.", category=UserWarning)
            elif used_id == "CDS":
                extra_suffixes.append("c_id")

            if unique_CDSs_per_gene:
                extra_suffixes.append("unique_per_gene")
            elif only_main:
                extra_suffixes.append("main")
            else:
                extra_suffixes.append("all")

        if verbose:
            extra_suffixes.append("coordinates")

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self._annot.feature_suffix, extension=extension, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder, extra_suffixes=extra_suffixes, use_name_not_id=use_name_not_id)

        with open(str(final_output_path), "w", encoding="utf-8") as f_out:
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
                            f_out.write(f">{c.id}")
                        elif used_id == "transcript":
                            f_out.write(f">{c.parents[0]}")
                        elif used_id == "gene":
                            f_out.write(f">{g.id}")

                        if verbose:
                            f_out.write(f"|{c.strand}|{c.ch}|{c.start}:{c.end}")

                        f_out.write(f"\n{c.seq}\n")


    def transcripts(self, only_main: bool = True, verbose: bool = True, used_id: str = "transcript", rna_classes: list = [], unique_transcripts_per_gene: bool = False, use_name_not_id: bool = False, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "features", extension=".fasta",
        #deprecated arguments
        custom_filename: str="", custom_path:str=""):
        """
        Main means only main transcript sequences are exported.

        Verbose will include strand, chromosome, and coordinates.

        Transcripts will be exported into fastas with their transcript ids unless 
        gene has been selected as used_id. To avoid duplicate
        equal entries, when choosing gene only main transcripts will
        be able to be output.

        valid_id_choices = ["gene", "transcript"]
        """

        if custom_filename != "":
            warnings.warn("'custom_filename' is deprecated. Please use 'filename' instead.", DeprecationWarning, stacklevel=2)
            if filename is None:
                filename = custom_filename

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        if subfolder_name != "features":
            subfolder = True

        valid_id_choices = ["gene", "transcript"]
        if used_id not in valid_id_choices:
            raise ValueError(f"used_id={used_id} is not amongst the valid_id_choices={valid_id_choices} to export transcripts.")

        extra_suffixes = ["transcripts"]

        if unique_transcripts_per_gene:
            only_main = False
            if used_id == "gene":
                used_id = "transcript"
                warnings.warn(f"Used id 'gene' has been changed to 'transcript' as unique transcripts per gene was selected.", category=UserWarning)

        if used_id == "gene":
            only_main = True
            extra_suffixes.append("g_id_main")
        elif used_id == "transcript":
            extra_suffixes.append("t_id")
            if unique_transcripts_per_gene:
                extra_suffixes.append("unique_per_gene")
            elif only_main:
                extra_suffixes.append("main")
            else:
                extra_suffixes.append("all")

        if verbose:
            extra_suffixes.append("coordinates")

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self._annot.feature_suffix, extension=extension, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder, extra_suffixes=extra_suffixes, use_name_not_id=use_name_not_id)

        with open(str(final_output_path), "w", encoding="utf-8") as f_out:

            for genes in self._annot.chrs.values():
                for g in genes.values():
                    temp_ts = []
                    for t in g.transcripts.values():
                        if (t.feature in rna_classes) or (not rna_classes):
                            if t.seq != "":
                                if only_main:
                                    if t.main:
                                        temp_ts.append(t)
                                else:
                                    temp_ts.append(t)

                    if not unique_transcripts_per_gene:
                        final_ts = temp_ts.copy()
                    else:
                        final_ts = []
                        if len(temp_ts) > 0:
                            final_ts.append(temp_ts[0])

                        for i, t1 in enumerate(temp_ts):
                            if i > 0:
                                add = True
                                for t2 in final_ts:
                                    if t1.seq == t2.seq:
                                        add = False
                                if add:
                                    final_ts.append(t1)

                    for t in final_ts:
                        if used_id == "transcript":    
                            f_out.write(f">{t.id}")
                        elif used_id == "gene":
                            f_out.write(f">{g.id}")

                        if verbose:
                            f_out.write(f"|{t.strand}|{t.ch}|{t.start}:{t.end}")
                        f_out.write(f"\n{t.seq}\n")

    def genes(self, verbose: bool = True, use_name_not_id: bool = False, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "features", extension=".fasta",
        #deprecated arguments
        custom_filename: str="", custom_path:str=""):

        if custom_filename != "":
            warnings.warn("'custom_filename' is deprecated. Please use 'filename' instead.", DeprecationWarning, stacklevel=2)
            if filename is None:
                filename = custom_filename

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        if subfolder_name != "features":
            subfolder = True

        extra_suffixes = ["genes"]
        
        if verbose:
            extra_suffixes.append("coordinates")

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self._annot.feature_suffix, extension=extension, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder, extra_suffixes=extra_suffixes, use_name_not_id=use_name_not_id)

        with open(str(final_output_path), "w", encoding="utf-8") as f_out:
            for genes in self._annot.chrs.values():
                for g in genes.values():
                    if g.seq != "":
                        f_out.write(f">{g.id}")
                        if verbose:
                            f_out.write(f"|{g.strand}|{g.ch}|{g.start}:{g.end}")
                        f_out.write(f"\n{g.seq}\n")

    def promoters(self, only_main: bool = True, verbose: bool = True, used_id: str = "promoter", promoter_size: int = 2000, promoter_type: str = "standard", use_name_not_id: bool = False, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "features", extension=".fasta", quiet:bool=False,
        #deprecated arguments
        custom_filename: str="", custom_path:str=""):
        """
        Verbose will include promoter type, strand, chromosome, and coordinates.
        """

        if custom_filename != "":
            warnings.warn("'custom_filename' is deprecated. Please use 'filename' instead.", DeprecationWarning, stacklevel=2)
            if filename is None:
                filename = custom_filename

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        if subfolder_name != "features":
            subfolder = True

        valid_id_choices = ["gene", "transcript", "promoter"]
        if used_id not in valid_id_choices:
            raise ValueError(f"used_id={used_id} is not amongst the valid_id_choices={valid_id_choices} to export promoters.")


        extra_suffixes = ["promoters"]

        if used_id == "gene":
            only_main = True
            extra_suffixes.append("g_id_main")
        else:
            if used_id == "transcript":
                extra_suffixes.append("t_id")
            if used_id == "promoter":
                extra_suffixes.append("prom_id")
            if only_main:
                extra_suffixes.append("main")
            else:
                extra_suffixes.append("all")

        if verbose:
            extra_suffixes.append("coordinates")

        extra_suffixes.append(str(self._annot.promoter_size))
        extra_suffixes.append(self._annot.promoter_type)

        if not self._annot.contains_promoters:
            self._annot.generate_promoters(promoter_size=promoter_size, promoter_type=promoter_type)

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self._annot.feature_suffix, extension=extension, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder, extra_suffixes=extra_suffixes, use_name_not_id=use_name_not_id)

        with open(str(final_output_path), "w", encoding="utf-8") as f_out:
            for genes in self._annot.chrs.values():
                for g in genes.values():
                    for t in g.transcripts.values():
                        if only_main and not t.main:
                            continue

                        if hasattr(t, "promoter"):
                            if t.promoter.seq != "":

                                if used_id == "promoter":
                                    f_out.write(f">{t.promoter.id}")
                                elif used_id == "transcript":
                                    f_out.write(f">{t.id}")
                                elif used_id == "gene":
                                    f_out.write(f">{g.id}")

                                if verbose:
                                    f_out.write(f"|{t.promoter.type}|{t.strand}|{t.ch}|{t.promoter.start}:{t.promoter.end}")

                                f_out.write(f"\n{t.promoter.seq}\n")
                            elif not quiet:
                                print(f"Warning: transcript {t.id} from annotation {self._annot.id} has no promoter sequence.\n")
                        elif not quiet:
                            print(f"Warning: transcript {t.id} from annotation {self._annot.id} has no promoter.\n")

    def for_dapseq(self, genome: Genome, genome_output_dir: str | None = None, gff_output_dir: str | None = None, genome_filepath: str | None = None, genome_filename: str | None = None, use_genome_dir: bool = False, genome_subfolder: bool = False, genome_subfolder_name: str = "out_genomes", genome_extension:str=".fasta", gff_filepath: str | None = None, gff_filename: str | None = None, use_annot_dir: bool = False, gff_subfolder: bool = False, gff_subfolder_name: str = "out_gffs", gff_extension:str=".gff3", gff_use_name_not_id:bool = False, skip_atypical_fts: bool = True, main_only: bool = False, UTRs: bool = False, exclude_non_coding: bool = False,
        #deprecated arguments
        genome_out_folder:str="", gff_out_folder:str="", tag:str="_for_dap.gff3"):


        if genome_out_folder != "":
            warnings.warn("'genome_out_folder' is deprecated. Please use 'genome_output_dir' instead.", DeprecationWarning, stacklevel=2)
            if genome_output_dir is None:
                genome_output_dir = genome_out_folder

        if gff_out_folder != "":
            warnings.warn("'gff_out_folder' is deprecated. Please use 'gff_output_dir' instead.", DeprecationWarning, stacklevel=2)
            if gff_output_dir is None:
                gff_output_dir = gff_out_folder

        if genome_subfolder_name != "out_genomes":
            genome_subfolder = True

        if gff_subfolder_name != "out_gffs":
            gff_subfolder = True

        if tag != "_for_dap.gff3":
            warnings.warn("'tag' is deprecated for naming files. Please use 'gff_filename' instead or 'gff_filepath' for the whole path.", DeprecationWarning, stacklevel=2)
            if gff_filename is None:
                gff_filename = tag

        equivalences = genome.rename_features_dap(output_dir=genome_output_dir, filename=genome_filename, filepath=genome_filepath, extension=genome_extension, use_genome_dir=use_genome_dir, subfolder_name=genome_subfolder_name, subfolder=genome_subfolder, return_equivalences=True, export=True)

        self._annot.rename_chromosomes(equivalences)

        extra_suffixes = ["for_dap"]

        final_output_path = self._resolve_output_path(filepath=gff_filepath, output_dir=gff_output_dir, filename=gff_filename, suffix=self._annot.all_suffixes, extension=gff_extension, use_annot_dir=use_annot_dir, subfolder_name=gff_subfolder_name, subfolder=gff_subfolder, extra_suffixes=extra_suffixes, use_name_not_id=gff_use_name_not_id)

        self.gff(filepath=str(final_output_path), skip_atypical_fts=skip_atypical_fts, main_only=main_only, UTRs=UTRs, just_genes=exclude_non_coding)

    def gff(self, 
            filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "out_gffs", extension:str=".gff3", skip_atypical_fts: bool = False, main_only: bool = False, UTRs: bool = False, just_genes: bool = False, no_1bp_features: bool = False, repeat_exons_utrs: bool = False, quiet: bool = False, skip_orphaned_fts: bool = False, featurecountsID: bool = False, extra_attributes: bool = False, clean_attributes: bool = True, aliases: bool = False, symbols: bool = False, symbols_as_description: bool = False, print_empty_attributes: bool = False, miRNAs: bool = True, clean_header: bool = False,
            # Deprecated arguments
            custom_path: str = "", tag: str = ".gff3"):

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        if tag != ".gff3":
            warnings.warn("'tag' is deprecated for naming files. Please use 'filename' instead or 'filepath' for the whole path.", DeprecationWarning, stacklevel=2)
            if filename is None:
                filename = tag

        progress_bar = start_progress_bar(total=len(self._annot.all_gene_ids.keys()), description=f"Exporting gff for {self._annot.id} annotation", quiet=quiet, colour="38;2;46;204;113")

        if subfolder_name != "out_gffs":
            subfolder = True

        if repeat_exons_utrs:
            self._annot.single_parent_for_exons_utrs()

        if featurecountsID:
            self._annot.create_featurecounts_ids()
            self._annot.tags.add("fcounts")
        else:
            self._annot.tags.discard("fcounts")
        
        if clean_attributes:
            self._annot.tags.add("clean")

        extra_suffixes = []
        if just_genes:
            extra_suffixes.append("just_genes")

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self._annot.all_suffixes, extra_suffixes=extra_suffixes, extension=extension, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder)

        if not quiet:
            print(f"Exporting {self._annot.id} gff to {final_output_path}.")

        with open(f"{final_output_path}", "w", encoding="utf-8") as f_out:
            if clean_header or "dapmod" in self._annot.tags or self._annot.gff_header == []:
                f_out.write("##gff-version 3\n")
            else:
                f_out.write("\n".join(self._annot.gff_header) + "\n")


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

                    f_out.write(g.print_gff(extra_attributes=extra_attributes, clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))

                    if just_genes:
                        continue

                    if repeat_exons_utrs or main_only or len(g.transcripts) == 1:
                        for t in g.transcripts.values():
                            if main_only:
                                if not t.main:
                                    continue
                            f_out.write(t.print_gff(featurecountsID=featurecountsID, clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))
                            for e in t.exons:
                                f_out.write(e.print_gff(featurecountsID=featurecountsID, clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))

                            if miRNAs:
                                for m in t.miRNAs:
                                    f_out.write(m.print_gff(featurecountsID=featurecountsID, clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))

                            for c in t.CDSs.values():
                                if main_only:
                                    if not c.main:
                                        continue
                                for c_seg in c.CDS_segments:
                                    f_out.write(c_seg.print_gff(featurecountsID=featurecountsID, clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))

                                if UTRs:
                                    if hasattr(c, "UTRs"):
                                        for u in c.UTRs:
                                            f_out.write(u.print_gff(featurecountsID=featurecountsID, clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))

                    else:
                                
                        for t in g.transcripts.values():                                
                            f_out.write(t.print_gff(featurecountsID=featurecountsID, clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))

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
                            f_out.write(e.print_gff(featurecountsID=featurecountsID, clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))


                        for t in g.transcripts.values():

                            for c in t.CDSs.values():
                                for c_seg in c.CDS_segments:
                                    f_out.write(c_seg.print_gff(featurecountsID=featurecountsID, clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))

                            if miRNAs:
                                for m in t.miRNAs:
                                    f_out.write(m.print_gff(featurecountsID=featurecountsID, clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))
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
                                f_out.write(u.print_gff(featurecountsID=featurecountsID, clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))

                    if x1 == (len(self._annot.chrs) - 1) and x2 == (len(genes) - 1):
                        continue
                    f_out.write("###\n")

            progress_bar.close()

            if not skip_atypical_fts:
                if not just_genes:
                    for x, ft in enumerate(self._annot.atypical_features):
                        if x == 0:
                            f_out.write("###\n")
                        f_out.write(ft.print_gff(clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))
                        if x == (len(self._annot.atypical_features) - 1):
                            continue
                        f_out.write("###\n")
            if not skip_orphaned_fts:
                if not just_genes:
                    for x, ft in enumerate(self._annot.orphaned_features):
                        if x == 0:
                            f_out.write("###\n")
                        f_out.write(ft.print_gff(clean=clean_attributes, aliases=aliases, symbols=symbols, symbols_as_description=symbols_as_description, print_empty_attributes=print_empty_attributes))
                        if x == (len(self._annot.orphaned_features) - 1):
                            continue
                        f_out.write("###\n")

    def gtf(self, 
            filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "out_gtfs", extension=".gtf", main_only: bool = False, UTRs: bool = False, just_genes: bool = False, no_1bp_features: bool = False, quiet: bool = False,
            # Deprecated arguments
            custom_path: str = "", tag: str = ".gtf"):

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        if tag != ".gtf":
            warnings.warn("'tag' is deprecated for naming files. Please use 'filename' instead or 'filepath' for the whole path.", DeprecationWarning, stacklevel=2)
            if filename is None:
                filename = tag

        if subfolder_name != "out_gtfs":
            subfolder = True

        self._annot.create_gtf_attributes()

        progress_bar = start_progress_bar(total=len(self._annot.all_gene_ids.keys()), description=f"Exporting gtf for {self._annot.id} annotation", quiet=quiet, colour="38;2;46;204;113")

        extra_suffixes = []

        if just_genes:
            extra_suffixes.append("just_genes")

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self._annot.all_suffixes, extra_suffixes=extra_suffixes, extension=extension, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder)

        if not quiet:
            print(f"Exporting {self._annot.id} gtf to {final_output_path}.")

        with open(f"{final_output_path}", "w", encoding="utf-8") as f_out:
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

    def gene_list(self, use_name_not_id: bool = False, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "lists", extension=".txt", lengths: bool = False, coordinates: bool = False, chromosomes: bool = False, coding_info: bool = False, skip_coding: bool = False, skip_non_coding: bool = False, sep: str = "\t", skip_pseudogenes: bool = False, skip_transposables: bool = False, gene_symbols: bool = False, include_header: bool = True, main_transcript_length_instead_of_gene_length: bool = False, main_gene_length_when_transcript_missing: bool = False, quiet:bool=False,
        #deprecated arguments
        custom_path: str = "", output_file: str = ""):

        if output_file != "":
            warnings.warn("'output_file' is deprecated. Please use 'filename' instead.", DeprecationWarning, stacklevel=2)
            if filename is None:
                filename = output_file

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        if subfolder_name != "lists":
            subfolder = True

        extra_suffixes = ["genes"]

        if skip_coding and skip_non_coding:
            print("Warning: Both skip_coding and skip_non_coding are set to True. No genes will be listed.")
            return

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, extra_suffixes=extra_suffixes, extension=extension, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder, use_name_not_id=use_name_not_id)

        if not quiet:
            print(f"Exporting {self._annot.id} gene list to {final_output_path}.")

        with open(final_output_path, "w", encoding="utf-8") as f_out:

            header = ["gene_id"]
            if chromosomes or coordinates:
                header.append("chromosome")
            if coordinates:
                header.append("gene_start")
                header.append("gene_end")
            if lengths:
                header.append("gene_length")
            if coding_info:
                header.append("coding")
            if gene_symbols:
                header.append("gene_symbol")
            if include_header:
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
                        if main_transcript_length_instead_of_gene_length:
                            t_size = 0
                            for t in g.transcripts.values():
                                if t.main:
                                    t_size = t.size
                                    break
                            if t_size == 0 and main_gene_length_when_transcript_missing:
                                out.append(str(g.size))
                            else:
                                out.append(str(t_size))
                        else:
                            out.append(str(g.size))
                    if coding_info:
                        out.append(str(g.coding))
                    if gene_symbols and g.symbols:
                        out.append("|".join(g.symbols))
                    f_out.write(sep.join(out) + "\n")

    def transcript_list(self, use_name_not_id: bool = False, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_annot_dir: bool = False, subfolder: bool = False, subfolder_name: str = "lists", extension=".txt", lengths: bool = False, coordinates: bool = False, chromosomes: bool = False, coding_info: bool = False, skip_coding: bool = False, skip_non_coding: bool = False, sep: str = "\t", skip_pseudogenes: bool = False, skip_transposables: bool = False, gene_symbols: bool = False, include_header: bool = True, quiet:bool=False,
        #deprecated arguments
        custom_path: str = "", output_file: str = ""):

        if output_file != "":
            warnings.warn("'output_file' is deprecated. Please use 'filename' instead.", DeprecationWarning, stacklevel=2)
            if filename is None:
                filename = output_file

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        if subfolder_name != "lists":
            subfolder = True

        extra_suffixes = ["transcripts"]

        if skip_coding and skip_non_coding:
            print("Warning: Both skip_coding and skip_non_coding are set to True. No transcripts will be listed.")
            return

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, extra_suffixes=extra_suffixes, extension=extension, use_annot_dir=use_annot_dir, subfolder_name=subfolder_name, subfolder=subfolder, use_name_not_id=use_name_not_id)

        if not quiet:
            print(f"Exporting {self._annot.id} transcript list to {final_output_path}.")

        with open(str(final_output_path), "w", encoding="utf-8") as f_out:

            header = ["transcript_id"]
            if chromosomes or coordinates:
                header.append("chromosome")
            if coordinates:
                header.append("transcript_start")
                header.append("transcript_end")
            if lengths:
                header.append("transcript_length")
            if coding_info:
                header.append("coding")
            if gene_symbols:
                header.append("gene_symbol")
            if include_header:
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
                        if gene_symbols and g.symbols:
                            out.append("|".join(g.symbols))
                        f_out.write(sep.join(out) + "\n")


