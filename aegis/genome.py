import time
import textwrap
import pandas as pd
import copy
import re
import random
import warnings
import hashlib

from typing import Literal

from os import system
from Bio import SeqIO
from pathlib import Path

from aegis.utils.misc import open_file
from aegis.utils.genefunctions import reverse_complement, sequence_hash

class Scaffold():
    # Compared in lower case, also as the prefix before the first "_" (e.g. chrUn_KI270302v1, chr00_scaffold_1)
    unknown_chromosome_names = {"chrun", "chr00", "chr0"}
    # Unlocalised, alternate and fix-patch scaffolds of a chromosome, e.g. chr1_random, chr4_KI270790v1_alt, chr1_NW_123.1
    chromosome_scaffold_name_pattern = re.compile(r"_(?:random|alt|fix)(?:_|$)|_[A-Z]{2}\d{6}v\d+|_N[TW]_\d", re.IGNORECASE)
    # Bare "m"/"c" are deliberately excluded: single-letter names are also used for nuclear chromosomes or linkage groups
    mitochondria_names = {"mt", "mito", "chrm", "chrmt", "mitochondria", "mitochondrion"}
    chloroplast_names = {"pt", "chrc", "chrpt", "pltd", "plastid", "chloroplast"}
    organelle_name_pattern = re.compile(r"^chr?[_.-]?([mc])$", re.IGNORECASE)
    mitochondria_description_pattern = re.compile(r"\bmitochondri(?:on|al|a)\b", re.IGNORECASE)
    chloroplast_description_pattern = re.compile(r"\b(?:chloroplast|plastid|apicoplast)\b", re.IGNORECASE)
    # DAP renaming (Genome.rename_features_dap) names organelles chr0001... and scaffolds chr00000001...
    dap_organelle_name_pattern = re.compile(r"^chr\d{4}$")
    dap_scaffold_name_pattern = re.compile(r"^chr\d{8}$")
    # e.g. "NC_000001.11 Homo sapiens chromosome 1, GRCh38.p14 Primary Assembly"
    chromosome_description_pattern = re.compile(r"\bchromosome\s+(\w+)", re.IGNORECASE)
    non_chromosome_description_pattern = re.compile(r"\b(?:unlocali[sz]ed|unplaced|scaffold|patch|contig|alternate)\b", re.IGNORECASE)

    @staticmethod
    def is_unknown_chromosome_name(name: str) -> bool:
        """True for unknown chromosome names, e.g. chrUn, chr00, chrUn_KI270302v1 or chr00_scaffold_1."""
        lower_name = name.lower()
        return lower_name in Scaffold.unknown_chromosome_names or lower_name.split("_")[0] in Scaffold.unknown_chromosome_names

    @staticmethod
    def description_chromosome(description: str) -> str | None:
        """
        Chromosome name (e.g. "1" or "X") given in a FASTA description, or None if there
        is none or the sequence is described as an unlocalized/unplaced scaffold, patch, etc.
        """
        if not description:
            return None
        match = Scaffold.chromosome_description_pattern.search(description)
        if not match or Scaffold.non_chromosome_description_pattern.search(description):
            return None
        if match.group(1).lower() in ("un", "unknown"):
            return None
        return match.group(1)

    @staticmethod
    def organelle_type(name: str, description: str = "") -> str | None:
        """
        Classifies a sequence as "mitochondria", "chloroplast" or None from its name
        (exact organelle names, or ch/chr followed by M or C) or, failing that, from
        organelle keywords in its FASTA description.
        """
        lower_name = name.lower()
        if lower_name in Scaffold.mitochondria_names:
            return "mitochondria"
        if lower_name in Scaffold.chloroplast_names:
            return "chloroplast"
        match = Scaffold.organelle_name_pattern.match(name)
        if match:
            return "mitochondria" if match.group(1).lower() == "m" else "chloroplast"
        if description:
            if Scaffold.mitochondria_description_pattern.search(description):
                return "mitochondria"
            if Scaffold.chloroplast_description_pattern.search(description):
                return "chloroplast"
        return None

    def __init__(self, name, sequence, original_name:str="", description:str="", dapmod:bool=False):

        self.name = name
        self.renamed = False
        self.dapfit = False
        self.dapmod = dapmod
        self.chromosome = False
        self.mitochondria = False
        self.chloroplast = False
        self.organelle = False
        # Set from Is_circular=true in a GFF region line; organelles are always treated as circular
        self.circular = False
        self.seq = sequence
        self.unknown_chromosome = False
        self.size = len(self.seq)
        self.description = description if description else name
        self._seq_hash: str | None = None
        self._rc_seq_hash: str | None = None

        if original_name:
            self.original_name = original_name
        else:
            self.original_name = self.name

        if self.original_name != self.name:
            self._rename_description(self.original_name, self.name)

        self.update()

    def _rename_description(self, old_name:str, new_name:str):
        """
        Gives the description (the FASTA header) new_name as its id, followed by the id of the input
        file (original_name) whenever the two differ, e.g. "NC_012007.3 Malus domestica chromosome 1"
        renamed chr01 becomes "chr01 NC_012007.3 Malus domestica chromosome 1". Descriptions not
        starting with old_name or original_name are left as they are.
        """
        parts = self.description.split(maxsplit=1)
        if not parts or parts[0] not in (old_name, self.original_name):
            return
        rest = parts[1] if len(parts) > 1 else ""
        # The original id is only given once, however many times the sequence is renamed
        rest_parts = rest.split(maxsplit=1)
        if rest_parts and rest_parts[0] == self.original_name:
            rest = rest_parts[1] if len(rest_parts) > 1 else ""
        # The input id of a reloaded _dapmod genome is not known (its placeholder is not given), but its
        # description already holds it if it was exported by aegis
        if new_name == self.original_name or self.original_name.startswith("unknown_dapmod_"):
            ids = [new_name]
        else:
            ids = [new_name, self.original_name]
        self.description = " ".join(ids + ([rest] if rest else []))

    def update(self, new_name:str=""):
        if new_name:
            if hasattr(self, "description") and self.description:
                self._rename_description(self.name, new_name)
            self.name = new_name

        if self.name != self.original_name:
            self.renamed = True

        match = re.search(r'(\d+(?:\.\d+)?)$', self.name)

        if match:
            self.number = match.group()
        else:
            self.number = ""

        # DAP names say what each sequence is, but scaffold names (chr00000001) would
        # otherwise be classified as chromosomes for starting with "ch". The organelle
        # type and unknown chromosomes can't be told from a DAP name, so they come from
        # the input id (original_name) when it is known.
        dap_scaffold = self.dapmod and Scaffold.dap_scaffold_name_pattern.match(self.name)
        dap_organelle = self.dapmod and Scaffold.dap_organelle_name_pattern.match(self.name)
        lower_name = self.name.lower()
        # Checked before the "ch" prefix, which would otherwise make all of them chromosomes
        chromosome_scaffold = Scaffold.chromosome_scaffold_name_pattern.search(self.name)
        if dap_scaffold:
            if Scaffold.is_unknown_chromosome_name(self.original_name):
                self.unknown_chromosome = True
        elif dap_organelle:
            self.chromosome = True
            self.organelle = True
        elif Scaffold.is_unknown_chromosome_name(self.name):
            self.unknown_chromosome = True
        elif chromosome_scaffold:
            pass
        elif self.renamed and not self.name.startswith("chr"):
            if len(self.number) == 2 or len(self.number) == 3:
                self.chromosome = True
        elif lower_name.startswith("ch"):
            self.chromosome = True

        # Accession-named chromosomes (e.g. NC_000001.11) are recognised from the description,
        # which also gives their number
        if not (self.chromosome and self.number.isdigit()) and not self.unknown_chromosome and not dap_scaffold and not chromosome_scaffold:
            description_chromosome = Scaffold.description_chromosome(self.description)
            if description_chromosome is not None:
                self.chromosome = True
                self.number = description_chromosome if description_chromosome.isdigit() else ""

        # Organelle flags are only ever set, never cleared: renamed organelles
        # (e.g. DAP renaming to chr0001) keep the classification of their original name.
        organelle = Scaffold.organelle_type(self.name, self.description)
        if organelle is None and dap_organelle:
            organelle = Scaffold.organelle_type(self.original_name)
        if organelle is not None:
            self.chromosome = True
            self.organelle = True
            if organelle == "mitochondria":
                self.mitochondria = True
            else:
                self.chloroplast = True

        if self.name.startswith("chr"):
            number_str = self.name[3:]
            if number_str.isdigit():
                self.dapfit = True

    @property
    def seq_hash(self) -> str:
        """Checksum of the sequence (cached)."""
        if self._seq_hash is None:
            self._seq_hash = sequence_hash(self.seq)
        return self._seq_hash

    @property
    def rc_seq_hash(self) -> str:
        """Checksum of the reverse-complemented sequence (cached)."""
        if self._rc_seq_hash is None:
            self._rc_seq_hash = sequence_hash(reverse_complement(self.seq))
        return self._rc_seq_hash

    @property
    def upper_seq(self) -> str:
        """Uppercase version of sequence."""
        return self.seq.upper()

    @property
    def soft_masked_bp(self) -> int:
        """Count of lowercase (soft-masked) nucleotides."""
        return sum(1 for ch in self.seq if ch.islower())

    @property
    def soft_masked_fraction(self) -> float:
        """Fraction of sequence that is soft-masked (0.0 to 1.0)."""
        return round(self.soft_masked_bp / self.size, 4) if self.size > 0 else 0.0

    def copy(self):
        return copy.deepcopy(self)

def classify_compartment(chrom: str, scaffold: Scaffold | None = None, mitochondria_chroms: tuple[str, ...] = (), chloroplast_chroms: tuple[str, ...] = (), auto_detection: bool = True) -> Literal["nuclear", "mitochondria", "chloroplast"]:
    """
    The one rule deciding whether a chromosome is "mitochondria", "chloroplast" or "nuclear": user-listed
    organelle chromosomes first; then, with auto_detection, the classification of its genome scaffold
    (name or FASTA description) or, without one, of its name (e.g. GFF column 1), both by
    Scaffold.organelle_type.
    """
    if chrom in mitochondria_chroms:
        return "mitochondria"
    if chrom in chloroplast_chroms:
        return "chloroplast"
    if not auto_detection:
        return "nuclear"
    if scaffold is not None:
        if scaffold.mitochondria:
            return "mitochondria"
        if scaffold.chloroplast:
            return "chloroplast"
        return "nuclear"
    organelle = Scaffold.organelle_type(chrom)
    if organelle == "mitochondria":
        return "mitochondria"
    if organelle == "chloroplast":
        return "chloroplast"
    return "nuclear"

class Genome():
    # FASTA files written by Genome.export in this run, to warn when a later export overwrites one
    _exported_files: set[str] = set()

    def __init__(self, name:str, genome_file_path:str, chromosome_dict:dict={}, rename_chromosomes:bool=False, quiet:bool=False,
                 header_id_tag:str|None=None, header_id_regex:str|None=None, gwh:bool=False, preserve_case:bool=True):
        start = time.time()
        self.name = name
        self.preserve_case = preserve_case

        self.file = str(Path(genome_file_path).resolve())
        self.path = str(Path(genome_file_path).resolve().parent) + "/"

        self.suffix = ""

        self.confrenamed = False

        self.dapmod = False

        self.dapfit = False

        # Set once DAP renaming has been applied, so that the genome stays "dapmod" even
        # after the renamed sequences (e.g. organelles) are removed
        self.dap_renamed = False

        self.unknown_chromosome = False

        # Output tags are read from the file name only, so that folder names (e.g. results_dapmod/) can't set them
        file_name = Path(self.file).name

        if "_dapmod" in file_name:
            self.dapmod = True
            self.dapfit = True
            self.dap_renamed = True

        if "_confrenamed" in file_name:
            self.confrenamed = True

        if "_scffree" in file_name:
            self.non_chromosomal_scaffolds = False
        
        if "_dapfit" in file_name:
            self.dapfit = True

        if "_organellefree" in file_name:
            self.organelles = False
            self.mitochondria = False
            self.chloroplast = False
        
        if "_chr00" in file_name:
            self.unknown_chromosome = True

        self.scaffolds = {}
        self.aliases = {}
        self.header_id_tag = header_id_tag
        self.header_id_regex = header_id_regex
        self.gwh = gwh

        if self.gwh and not self.header_id_tag and not self.header_id_regex:
            self.header_id_tag = "OriSeqID"

        # Creating a dictionary with the genome sequence of each chromosome or scaffold, still referred as chromosomes in the code
        self.chromosome_dict = chromosome_dict

        self.equivalences = {}

        count = 0
        with open_file(self.file, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                raw_id = record.id
                desc = record.description if record.description else raw_id
                scaffold_id = raw_id

                if self.header_id_regex:
                    match = re.search(self.header_id_regex, desc)
                    if match and match.groups():
                        scaffold_id = match.group(1).strip()
                elif self.header_id_tag:
                    match = re.search(rf"(?:^|[\s;,]){re.escape(self.header_id_tag)}=([^\s;,]+)", desc)
                    if match:
                        scaffold_id = match.group(1).strip()

                count += 1
                if scaffold_id in self.scaffolds:
                    print((f"Error: scaffold feature {scaffold_id} is repeated in {self.name}, genome (file: {self.file})"))
                raw_seq = str(record.seq) if self.preserve_case else str(record.seq).upper()
                if self.dapmod:
                    original_name = Genome._dap_original_name(scaffold_id, desc) or f"unknown_dapmod_{count}"
                    self.scaffolds[scaffold_id] = Scaffold(scaffold_id, raw_seq, original_name=original_name, description=desc, dapmod=True)
                    self.equivalences[scaffold_id] = scaffold_id
                elif rename_chromosomes and scaffold_id in self.chromosome_dict:
                    self.confrenamed = True
                    self.scaffolds[self.chromosome_dict[scaffold_id]] = Scaffold(self.chromosome_dict[scaffold_id], raw_seq, scaffold_id, description=desc)
                    self.equivalences[scaffold_id] = self.chromosome_dict[scaffold_id]
                else:
                    self.scaffolds[scaffold_id] = Scaffold(scaffold_id, raw_seq, original_name=raw_id, description=desc)
                    self.equivalences[scaffold_id] = scaffold_id

                if raw_id != scaffold_id:
                    self.aliases[raw_id] = scaffold_id
                    self.equivalences[raw_id] = scaffold_id

        self.update()

        now = time.time()
        lapse = now - start
        if not quiet:
            print(f"\n{self.name} genome chromosomes/scaffolds: {self.features[:30]} ...")
            print(f"\nCreating {self.name} genome object took {round(lapse, 1)} seconds\n")

    @staticmethod
    def _dap_original_name(name: str, description: str) -> str:
        """
        Input id of a sequence renamed by DAP renaming, read back from a _dapmod FASTA file: the second
        header word, which Genome.export writes after the new name (e.g. "chr00000001 scf1 some scaffold").
        Only DAP scaffold and organelle names (chr00000001, chr0001) always come from renaming; for any
        other name, or a header without a second word, "" is returned.
        """
        if not (Scaffold.dap_scaffold_name_pattern.match(name) or Scaffold.dap_organelle_name_pattern.match(name)):
            return ""
        parts = description.split()
        if len(parts) < 2 or parts[0] != name:
            return ""
        return parts[1]

    def get_scaffold(self, scaffold_id: str) -> Scaffold | None:
        if scaffold_id in self.scaffolds:
            return self.scaffolds[scaffold_id]
        if scaffold_id in self.aliases and self.aliases[scaffold_id] in self.scaffolds:
            return self.scaffolds[self.aliases[scaffold_id]]
        return None

    def __getitem__(self, scaffold_id: str) -> Scaffold:
        scf = self.get_scaffold(scaffold_id)
        if scf is not None:
            return scf
        raise KeyError(scaffold_id)

    def __contains__(self, scaffold_id: str) -> bool:
        return scaffold_id in self.scaffolds or (scaffold_id in self.aliases and self.aliases[scaffold_id] in self.scaffolds)

    def update(self, update_scaffolds:bool=False):

        self.suffix = ""

        self.dapfit = True
        self.mitochondria = False
        self.chloroplast = False
        self.organelles = False
        self.unknown_chromosome = False
        self.non_chromosomal_scaffolds = False
        # Non-chromosomal scaffolds besides unknown chromosomes (chr00), which are marked by _chr00 instead of preventing _scffree
        self.other_scaffolds = False
        self.chromosomes = False
        self.dapmod = False
        self.confrenamed = False

        self.features = list(self.scaffolds.keys())
        self.size = sum(scaffold.size for scaffold in self.scaffolds.values())
        self.chromosome_size = sum(scaffold.size for scaffold in self.scaffolds.values() if scaffold.chromosome)
        self.nuclear_chromosome_size = sum(scaffold.size for scaffold in self.scaffolds.values() if scaffold.chromosome and not scaffold.organelle)

        self.chromosome_names = set()
        self.scaffold_names = set()
        self.accessory_chromosome_names = set()

        for scaffold in self.scaffolds.values():

            if update_scaffolds:
                scaffold.update()

            if not scaffold.name.startswith("chr"):
                self.dapfit = False
            else:
                try:
                    _ = int(scaffold.name.split("r")[1])
                except:
                    self.dapfit = False

            if scaffold.dapmod:
                self.dapmod = True
            elif scaffold.renamed:
                self.confrenamed = True

            if scaffold.unknown_chromosome:
                self.unknown_chromosome = True

            if scaffold.chromosome:
                self.chromosomes = True
                if scaffold.mitochondria:
                    self.mitochondria = True
                    self.accessory_chromosome_names.add(scaffold.name)
                elif scaffold.chloroplast:
                    self.chloroplast = True
                    self.accessory_chromosome_names.add(scaffold.name)
                elif scaffold.organelle:
                    # organelle of unknown type, e.g. chr0001 in a reloaded DAP genome
                    self.organelles = True
                    self.accessory_chromosome_names.add(scaffold.name)
                else:
                    self.chromosome_names.add(scaffold.name)
            else:
                self.non_chromosomal_scaffolds = True
                self.scaffold_names.add(scaffold.name)
                if not scaffold.unknown_chromosome:
                    self.other_scaffolds = True

            if update_scaffolds:
                scaffold.update()

        if self.mitochondria or self.chloroplast:
            self.organelles = True

        if self.dap_renamed:
            self.dapmod = True

        if self.dapmod:
            self.suffix += "_dapmod"
        
        if self.confrenamed:
            self.suffix += "_confrenamed"

        if not self.other_scaffolds:
            self.suffix += "_scffree"

        if self.dapfit and not self.dapmod:
            self.suffix += "_dapfit"

        if not self.organelles:
            self.suffix += "_organellefree"

        if self.unknown_chromosome:
            self.suffix += "_chr00"

    def _resolve_output_path(self, filepath: str | None = None, output_dir: str | None = None, 
        filename: str | None = None, suffix: str = "", extra_suffixes: list[str] | None = None, extension: str = ".fasta", subfolder_name: str = "out_genomes", subfolder: bool = False, use_genome_dir: bool = False):

        if filepath:
            out_path = Path(filepath)
            out_path.parent.mkdir(parents=True, exist_ok=True)
            if (output_dir is not None) or subfolder or filename or use_genome_dir:
                warnings.warn(f"Exact output file path ({filepath}) was specified. Ignoring output_dir, use_genome_dir, subfolder, and filename arguments.")
            return out_path
        else:
            if use_genome_dir:
                export_folder = Path(self.path)
                if output_dir is not None:
                    warnings.warn(f"Both 'use_genome_dir={use_genome_dir}' and 'output_dir={output_dir}' were provided. Defaulting to the genome directory ({self.path}).")

            else:
                export_folder = Path(output_dir or ".")

            if subfolder:
                export_folder = export_folder / subfolder_name.strip("/")
            
            export_folder.mkdir(parents=True, exist_ok=True)

            if extension and not extension.startswith("."):
                extension = f".{extension}"

            if not filename:
                tag_str = "".join([f"_{t}" for t in (extra_suffixes or []) if t])
                filename = f"{self.name}{tag_str}{suffix}{extension}"

            if not Path(filename).suffix:
                filename = f"{filename}{extension}"

            out_path = export_folder / filename

            return out_path

    def export_feature_sizes(self, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_genome_dir: bool = False, subfolder: bool = False, subfolder_name: str = "coordinates", extension=".tsv",quiet:bool=False,
        #deprecated arguments
        custom_path:str=""):

        if custom_path != "":
            warnings.warn("'custom_path' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = custom_path

        self.update()

        extra_suffixes = ["genome_feature_sizes"]

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self.suffix, extra_suffixes=extra_suffixes, extension=extension, use_genome_dir=use_genome_dir, subfolder_name=subfolder_name, subfolder=subfolder)

        out = []
        for scaffold_id, scaffold in self.scaffolds.items():
            out.append(f"{scaffold_id}\t{scaffold.size}")
        if out:
            out[-1] += "\n"
            with open(str(final_output_path), "w", encoding="utf-8") as f_out:
                f_out.write("\n".join(out))
        elif not quiet:
            print(f"No scaffolds/chromosomes found in {self.name} genome object.")

    def rename_features_dap(self, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_genome_dir: bool = False, subfolder: bool = False, subfolder_name: str = "out_genomes", extension=".fasta", return_equivalences:bool=False, export:bool=False, keep_description:bool=True,
        #deprecated argument
        output_folder:str=""):
        """
        Renames scaffolds and chromosomes to become dapfit.
        """

        if output_folder != "":
            warnings.warn("'output_folder' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = output_folder

        if not self.dapfit:

            # checking to see if all chromosomes have unique numbers, else they
            # will be all renumbered

            use_chromosome_count = False
            chromosome_numbers = []
            nuclear_chromosome_found = False

            for scaffold in self.scaffolds.values():
                if scaffold.chromosome:
                    if not scaffold.mitochondria and not scaffold.chloroplast:
                        if not scaffold.organelle:
                            nuclear_chromosome_found = True
                        if not scaffold.number.isdigit():
                            use_chromosome_count = True
                        else:
                            if scaffold.number in chromosome_numbers:
                                use_chromosome_count = True
                            chromosome_numbers.append(scaffold.number)

            # Otherwise every nuclear sequence silently becomes a chr00000001-style scaffold
            if not nuclear_chromosome_found:
                warnings.warn(f"No recognised nuclear chromosomes in {self.name} genome: all {len(self.scaffolds)} sequences were classed as scaffolds or organelles, so DAP renaming gives every non-organelle sequence a scaffold name (chr00000001...). Name its chromosomes with a chromosome dictionary (rename_chromosomes=True).")

            # A single unknown chromosome (chrUn, chr00...) is named chr00, which keeps it recognisable;
            # several of them get scaffold names, as does one whose name a chromosome numbered 0 takes
            unknown_chromosome_count = sum(1 for scaffold in self.scaffolds.values() if scaffold.unknown_chromosome and not scaffold.chromosome)
            chromosome_zero = not use_chromosome_count and any(int(number) == 0 for number in chromosome_numbers)
            keep_unknown_chromosome = unknown_chromosome_count == 1 and not chromosome_zero

            organelle_count = 0
            scaffold_count = 0
            chromosome_count = 0

            for scaffold in self.scaffolds.values():

                if scaffold.mitochondria or scaffold.chloroplast:
                    organelle_count += 1
                    number = "{:04d}".format(organelle_count)
                    scaffold.number = number
                    self.equivalences[scaffold.original_name] = f"chr{scaffold.number}"

                elif scaffold.chromosome:
                    if use_chromosome_count:
                        chromosome_count += 1
                        scaffold.number = "{:02d}".format(chromosome_count)
                        self.equivalences[scaffold.original_name] = f"chr{scaffold.number}"
                    else:
                        scaffold.number = "{:02d}".format(int(scaffold.number))
                        self.equivalences[scaffold.original_name] = f"chr{scaffold.number}"
                elif keep_unknown_chromosome and scaffold.unknown_chromosome:
                    scaffold.number = "00"
                    self.equivalences[scaffold.original_name] = "chr00"
                else:
                    scaffold_count += 1
                    scaffold.number = "{:08d}".format(scaffold_count)
                    self.equivalences[scaffold.original_name] = f"chr{scaffold.number}"

            new_scaffolds = {}

            for scaffold in self.scaffolds.values():
                previous_name = scaffold.name
                was_chromosome = scaffold.chromosome
                new_name = self.equivalences[scaffold.original_name]
                scaffold.update(new_name=new_name)
                # update() classifies any name starting with "ch" as a chromosome
                scaffold.chromosome = was_chromosome
                if previous_name != scaffold.name:
                    scaffold.dapmod = True
                new_scaffolds[scaffold.name] = scaffold.copy()

            self.scaffolds = new_scaffolds.copy()
            del new_scaffolds

            self.dap_renamed = True
            self.update()

            if export:
                self.export(filepath=filepath, output_dir=output_dir, use_genome_dir=use_genome_dir, subfolder=subfolder, subfolder_name=subfolder_name, filename=filename, extension=extension, keep_description=keep_description)

            if return_equivalences:
                return self.equivalences
            
        else:

            print(f"Warning: {self.name} genome is already fit for DAP-Seq analysis, so it will not be modified.")

    def rename_features_from_dic(self, rename_map: dict) -> dict:

        new_scaffolds = {}
        final_equivalences = {}

        for scaffold in self.scaffolds.values():
            
            original_name = scaffold.original_name
            new_name = rename_map.get(scaffold.name, rename_map.get(original_name, scaffold.name))
            
            copied_scaffold = scaffold.copy()
            copied_scaffold.update(new_name=new_name)

            final_equivalences[original_name] = new_name
            final_equivalences[scaffold.name] = new_name

            new_scaffolds[new_name] = copied_scaffold

        self.scaffolds = new_scaffolds
        self.update()

        return final_equivalences

    def remove_scaffolds(self, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_genome_dir: bool = False, subfolder: bool = False, subfolder_name: str = "out_genomes", extension=".fasta", export:bool=False, remove_00:bool=True, remove_organelles:bool=False, keep_description:bool=True):

        sequence_count = len(self.scaffolds)

        if self.non_chromosomal_scaffolds:
            new_scaffolds = {}

            for scaffold_id, scaffold in self.scaffolds.items():
                if scaffold.chromosome:
                    new_scaffolds[scaffold_id] = scaffold.copy()
                elif not remove_00 and scaffold.unknown_chromosome:
                    new_scaffolds[scaffold_id] = scaffold.copy()

            self.scaffolds = new_scaffolds.copy()
            del new_scaffolds

            self.update()

            if not self.chromosome_names:
                kept = f"only organelles ({', '.join(sorted(self.accessory_chromosome_names))}) are left" if self.scaffolds else "nothing is left"
                warnings.warn(f"No recognised nuclear chromosomes in {self.name} genome: {sequence_count} sequences were classed as scaffolds or organelles, so after removing scaffolds {kept}. Name its chromosomes with a chromosome dictionary (rename_chromosomes=True).")
                # An organelle-only genome is not a scaffold-free version of the input
                if export:
                    warnings.warn(f"There was nothing to export for {self.name} genome.")
                    export = False

        if remove_organelles:
            self.remove_organelles()

        if export:
            # e.g. chr00 is the only scaffold and remove_00=False: the result is the full genome under the same name
            if len(self.scaffolds) == sequence_count:
                removed = "scaffolds or organelles were" if remove_organelles else "scaffolds were"
                version = "scaffold- and organelle-free" if remove_organelles else "scaffold-free"
                warnings.warn(f"No {removed} removed from {self.name} genome, so its {version} version was not exported.")
            else:
                self.export(filepath=filepath, output_dir=output_dir, use_genome_dir=use_genome_dir, subfolder=subfolder, subfolder_name=subfolder_name, filename=filename, extension=extension, keep_description=keep_description)

    def remove_organelles(self, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_genome_dir: bool = False, subfolder: bool = False, subfolder_name: str = "out_genomes", extension=".fasta", export:bool=False, remove_mitochondria:bool=True, remove_chloroplast:bool=True, keep_description:bool=True):

        sequence_count = len(self.scaffolds)
        new_scaffolds = {}

        for scaffold_id, scaffold in self.scaffolds.items():
            if remove_mitochondria:
                if scaffold.mitochondria:
                    continue
            if remove_chloroplast:
                if scaffold.chloroplast:
                    continue
            # organelles of unknown type are removed only when removing both kinds
            if remove_mitochondria and remove_chloroplast and scaffold.organelle:
                continue
            new_scaffolds[scaffold_id] = scaffold.copy()

        self.scaffolds = new_scaffolds.copy()
        del new_scaffolds

        self.update()

        if export:
            if len(self.scaffolds) == sequence_count:
                warnings.warn(f"No organelles were removed from {self.name} genome, so its organelle-free version was not exported.")
            else:
                self.export(filepath=filepath, output_dir=output_dir, use_genome_dir=use_genome_dir, subfolder=subfolder, subfolder_name=subfolder_name, filename=filename, extension=extension, keep_description=keep_description)

    def export(self, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_genome_dir: bool = False, subfolder: bool = False, subfolder_name: str = "out_genomes", extension=".fasta", quiet:bool=False, keep_description: bool = True,
        #deprecated arguments
        output_folder:str="", file:str=".fasta"):

        if output_folder != "":
            warnings.warn("'output_folder' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = output_folder

        if file != ".fasta":
            warnings.warn("'file' is deprecated. Please use 'filename' instead.", DeprecationWarning, stacklevel=2)
            if filename is None:
                filename = file

        self.update()

        # Tabs (e.g. GWH "OriSeqID=chr1\tLen=120") are handled badly by some tools
        headers = [scaffold.description.replace("\t", " ") if (keep_description and getattr(scaffold, "description", None)) else scaffold.name for scaffold in self.scaffolds.values()]
        # Headers that are only ids, whether descriptions were dropped or never given, are marked by _no_description
        suffix = self.suffix
        if headers and all(len(header.split()) == 1 for header in headers):
            suffix += "_no_description"

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=suffix, use_genome_dir=use_genome_dir, subfolder_name=subfolder_name, subfolder=subfolder, extension=extension)

        if not self.scaffolds:
            warnings.warn(f"There was nothing to export for {self.name} genome.")
            return

        resolved_output_path = str(Path(final_output_path).resolve())
        if resolved_output_path in Genome._exported_files:
            warnings.warn(f"{final_output_path} was already written earlier in this run and is being overwritten by {self.name} genome (suffix '{suffix}').")
        Genome._exported_files.add(resolved_output_path)

        if not quiet:
            print(f"Exporting {self.name} genome to {final_output_path}.")

        with open(str(final_output_path), "w", encoding="utf-8") as f_out:
            for scaffold, header in zip(self.scaffolds.values(), headers):
                f_out.write(f">{header}\n{scaffold.seq}\n")

    def copy(self):
        return copy.deepcopy(self)
        
    def extract_peak_sequences(self, DAPseq_output_file:str, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_genome_dir: bool = False, subfolder: bool = False, subfolder_name: str = "out_peak_seqs", extension=".fasta", top=600,
        #deprecated arguments
        output_file_name:str="", output_folder:str=""):

        if output_folder != "":
            warnings.warn("'output_folder' is deprecated. Please use 'output_dir' instead.", DeprecationWarning, stacklevel=2)
            if output_dir is None:
                output_dir = output_folder

        if output_file_name != "":
            warnings.warn("'output_file_name' is deprecated. Please use 'filename' instead.", DeprecationWarning, stacklevel=2)
            if filename is None:
                filename = output_file_name

        try:
            df = pd.read_csv(DAPseq_output_file, delimiter='\t', dtype=str)
            df.dropna(how='all', inplace=True)
        except FileNotFoundError:
            print(f"Error: DAPseq file not found: {DAPseq_output_file}")
            return

        if top == 'all':
            top_df = df
        else:
            df['score'] = pd.to_numeric(df['score'], errors='coerce')
            df.dropna(subset=['score'], inplace=True)
            top_df = df.nlargest(top, 'score')

        peaks = {}

        for index, row in top_df.iterrows():
            try:

                seq_name = row['seqnames']
                feature = row['feature']
                peak_str = row['peak']

                peak_position = int(peak_str.split(':')[1]) # type: ignore

                scaffold_seq = self.scaffolds[seq_name].seq
                scaffold_len = len(scaffold_seq)

                start_pos = max(0, peak_position - 100)
                end_pos = min(scaffold_len, peak_position + 100)
                
                if start_pos >= end_pos:
                    print(f"Warning: Invalid range for row {index}. Skipping.")
                    continue

                extracted_seq = scaffold_seq[start_pos:end_pos]
                header = f"{seq_name}_{feature}_{start_pos}:{end_pos}"

                peaks[header] = str(extracted_seq)

            except (IndexError, ValueError):

                print(f"Warning: Incorrect peak format in row {index} (value: '{row.get('peak', 'N/A')}'). Skipping.")
            except KeyError:

                print(f"Warning: Scaffold '{row.get('seqnames', 'N/A')}' not found in reference. Skipping row {index}.")
            except Exception as e:

                print(f"Error processing row {index}: {e}. Skipping.")


        if not peaks:
            print("Warning: No peaks extracted for any sequence.")
            return

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self.suffix, extension=extension, use_genome_dir=use_genome_dir, subfolder_name=subfolder_name, subfolder=subfolder)

        with open(str(final_output_path), "w", encoding="utf-8") as f_out:
            for header, seq in peaks.items():
                f_out.write(f'>{header}\n')
                f_out.write(f'{textwrap.fill(seq, width=60)}\n')

    def subset(self, chosen_features:set|list|tuple|None=None, cap:int=2, quiet:bool=False):

        if chosen_features is None:
            chosen_features = set()
        elif isinstance(chosen_features, (list, tuple)):
            chosen_features = set(chosen_features)

        if chosen_features:
            for chosen_feature in chosen_features:
                if chosen_feature not in self.scaffolds:
                    raise ValueError(f"Chosen scaffold/chromosome {chosen_feature} is not in {self.name} genome. Choose from '{self.scaffolds.keys()}'")
            scaffolds_to_remove = set(self.scaffolds) - chosen_features
        else:
            if cap > len(self.scaffolds) and not quiet:
                warnings.warn(f"Cap value {cap} exceeds the number of available scaffolds/chrosomomes ({len(self.scaffolds)}). No features removed in subset genome {self.name}.", category=UserWarning)
                return
            scaffolds_to_remove = set(self.scaffolds) - set(random.sample(list(self.scaffolds), cap))

        if scaffolds_to_remove:
            self.remove_features(scaffolds_to_remove)
        elif not quiet:
            print(f"No scaffolds/chromosomes removed from {self.name} genome.")

    def remove_features(self, features_to_remove:set):
        for ft in features_to_remove:
            del self.scaffolds[ft]
        self.update()

    def get_stats(self, estimated_genome_size: int | None = None) -> dict:
        """
        Calculates and returns assembly statistics for the genome.

        Args:
            estimated_genome_size (int, optional): Expected/estimated genome size (in bp)
                for calculating NG50, LG50, and auNG.
        """
        self.update()

        lengths = [s.size for s in self.scaffolds.values()]
        lengths_sorted = sorted(lengths, reverse=True)
        total_size = sum(lengths)

        # N50 and L50
        half_size = total_size / 2.0
        n50, l50 = 0, 0
        running_sum = 0
        for idx, length in enumerate(lengths_sorted, 1):
            running_sum += length
            if running_sum >= half_size:
                n50 = length
                l50 = idx
                break

        # N90 and L90
        target_90 = total_size * 0.9
        n90, l90 = 0, 0
        running_sum = 0
        for idx, length in enumerate(lengths_sorted, 1):
            running_sum += length
            if running_sum >= target_90:
                n90 = length
                l90 = idx
                break

        # auN (length-weighted mean contig length / area under Nx curve)
        sum_sq = sum(length ** 2 for length in lengths)
        au_n = round(sum_sq / total_size) if total_size > 0 else 0

        # NG50, LG50, auNG (if estimated_genome_size is provided)
        ng50, lg50, au_ng = None, None, None
        if estimated_genome_size is not None and estimated_genome_size > 0:
            ng_half = estimated_genome_size / 2.0
            running_sum = 0
            ng50, lg50 = 0, 0
            for idx, length in enumerate(lengths_sorted, 1):
                running_sum += length
                if running_sum >= ng_half:
                    ng50 = length
                    lg50 = idx
                    break
            au_ng = round(sum_sq / estimated_genome_size)

        # GC and Gap content
        total_gc = 0
        total_gap = 0
        total_soft_masked = 0
        for scf in self.scaffolds.values():
            seq = scf.seq
            total_gc += seq.count('G') + seq.count('C') + seq.count('g') + seq.count('c')
            total_gap += seq.count('N') + seq.count('n')
            total_soft_masked += scf.soft_masked_bp

        gc_content = round((total_gc / total_size * 100), 2) if total_size > 0 else 0.0
        gap_content = round((total_gap / total_size * 100), 2) if total_size > 0 else 0.0
        soft_masked_pct = round((total_soft_masked / total_size * 100), 2) if total_size > 0 else 0.0

        if self.scaffolds:
            longest = max(self.scaffolds.values(), key=lambda s: s.size)
            shortest = min(self.scaffolds.values(), key=lambda s: s.size)
            longest_scaffold = (longest.name, longest.size)
            shortest_scaffold = (shortest.name, shortest.size)
        else:
            longest_scaffold = ("", 0)
            shortest_scaffold = ("", 0)

        num_chromosomes = len(self.chromosome_names) + len(self.accessory_chromosome_names)
        if self.unknown_chromosome:
            num_chromosomes += sum(1 for s in self.scaffolds.values() if s.unknown_chromosome)

        scaffold_size = total_size - self.chromosome_size

        res = {
            "total_size": total_size,
            "chromosome_size": self.chromosome_size,
            "nuclear_chromosome_size": self.nuclear_chromosome_size,
            "scaffold_size": scaffold_size,
            "num_sequences": len(self.scaffolds),
            "num_chromosomes": num_chromosomes,
            "num_scaffolds": len(self.scaffolds) - num_chromosomes,
            "n50": n50,
            "l50": l50,
            "n90": n90,
            "l90": l90,
            "auN": au_n,
            "gc_content": gc_content,
            "gap_content": gap_content,
            "soft_masked_bp": total_soft_masked,
            "soft_masked_pct": soft_masked_pct,
            "longest_scaffold": longest_scaffold,
            "shortest_scaffold": shortest_scaffold,
        }
        if estimated_genome_size is not None and estimated_genome_size > 0:
            res["estimated_genome_size"] = estimated_genome_size
            res["ng50"] = ng50
            res["lg50"] = lg50
            res["auNG"] = au_ng

        return res

    @property
    def stats(self) -> dict:
        return self.get_stats()

    def get_sorted_features(self, chromosomes_only: bool = False, sort_by: str = "name") -> list[str]:
        """
        Returns a sorted list of scaffold/chromosome names.
        
        Args:
            chromosomes_only (bool): If True, only include chromosomes (nuclear, unknown, organellar).
            sort_by (str): 'name' (natural alphanumeric sort), 'size' (descending size), or 'order' (FASTA order).
        """
        self.update()

        candidates = []
        for name, scf in self.scaffolds.items():
            is_chr = scf.chromosome or scf.unknown_chromosome or scf.organelle
            if chromosomes_only and not is_chr:
                continue
            candidates.append(scf)

        if sort_by == "order":
            return [s.name for s in candidates]

        if sort_by == "size":
            return [s.name for s in sorted(candidates, key=lambda s: (-s.size, s.name))]

        # sort_by == "name": category hierarchy + natural sort
        def _sort_key(scf):
            if scf.chromosome and not scf.organelle and not scf.unknown_chromosome:
                category = 1
            elif scf.unknown_chromosome:
                category = 2
            elif scf.organelle:
                category = 3
            else:
                category = 4
            parts = [int(text) if text.isdigit() else text.lower() for text in re.split(r'(\d+)', scf.name)]
            return (category, parts)

        return [s.name for s in sorted(candidates, key=_sort_key)]