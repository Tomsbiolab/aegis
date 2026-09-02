import time
import textwrap
import pandas as pd
import copy
import re
import random
import warnings

from os import system
from Bio import SeqIO
from pathlib import Path

from aegis.utils.misc import open_file

class Scaffold():
    mitochondria_suffixes = ["m", "M"]
    chloroplast_suffixes = ["c", "C"]
    unknown_chromosome_names = ["chrUn", "chrun", "ChrUn", "Chrun", "chr00", "Chr00", "chr0", "Chr0"]
    organelle_suffixes = mitochondria_suffixes + chloroplast_suffixes
    def __init__(self, name, sequence, original_name:str="", description:str=""):

        self.name = name
        self.renamed = False
        self.dapfit = False
        self.dapmod = False
        self.chromosome = False
        self.mitochondria = False
        self.chloroplast = False
        self.organelle = False
        self.seq = sequence
        self.unknown_chromosome = False
        self.size = len(self.seq)
        self.description = description if description else name

        if original_name:
            self.original_name = original_name
        else:
            self.original_name = self.name

        if self.description and self.original_name != self.name:
            if self.description.startswith(self.original_name):
                self.description = self.name + self.description[len(self.original_name):]

        self.update()

    def update(self, new_name:str=""):
        if new_name:
            if hasattr(self, "description") and self.description:
                if self.description.startswith(self.name):
                    self.description = new_name + self.description[len(self.name):]
                elif hasattr(self, "original_name") and self.description.startswith(self.original_name):
                    self.description = new_name + self.description[len(self.original_name):]
            self.name = new_name

        if self.name != self.original_name:
            self.renamed = True

        match = re.search(r'(\d+(?:\.\d+)?)$', self.name)

        if match:
            self.number = match.group()
        else:
            self.number = ""

        if self.renamed and not self.name.startswith("chr"):
            if len(self.number) == 2 or len(self.number) == 3 or self.name[-1].lower() in Scaffold.organelle_suffixes:
                self.chromosome = True
                if self.name[-1].lower() in Scaffold.organelle_suffixes:
                    self.organelle = True
                    if self.name[-1].lower() in Scaffold.mitochondria_suffixes:
                        self.mitochondria = True
                    else:
                        self.chloroplast = True
        else:
            if self.name.lower() in Scaffold.unknown_chromosome_names:
                self.unknown_chromosome = True
            elif self.name.lower().startswith("ch"):
                self.chromosome = True
                if self.name[-1].lower() in Scaffold.organelle_suffixes:
                    self.organelle = True
                    if self.name[-1].lower() in Scaffold.mitochondria_suffixes:
                        self.mitochondria = True
                    else:
                        self.chloroplast = True

        if self.name.startswith("chr"):
            number_str = self.name[3:]
            if number_str.isdigit():
                self.dapfit = True

    def copy(self):
        return copy.deepcopy(self)

class Genome():
    def __init__(self, name:str, genome_file_path:str, chromosome_dict:dict={}, rename_chromosomes:bool=False, quiet:bool=False,
                 header_id_tag:str|None=None, header_id_regex:str|None=None, gwh:bool=False):
        start = time.time()
        self.name = name

        self.file = str(Path(genome_file_path).resolve())
        self.path = str(Path(genome_file_path).resolve().parent) + "/"

        self.suffix = ""

        self.confrenamed = False

        self.dapmod = False

        self.dapfit = False

        self.unknown_chromosome = False

        if "_dapmod" in self.file:
            self.dapmod = True
            self.dapfit = True

        if "_confrenamed" in self.file:
            self.confrenamed = True

        if "_scffree" in self.file:
            self.non_chromosomal_scaffolds = False
        
        if "_dapfit" in self.file:
            self.dapfit = True

        if "_organellefree" in self.file:
            self.organelles = False
            self.mitochondria = False
            self.chloroplast = False
        
        if "_chr00" in self.file:
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
                if self.dapmod:
                    self.scaffolds[scaffold_id] = Scaffold(scaffold_id, str(record.seq).upper(), original_name=f"unknown_dapmod_{count}", description=desc)
                    self.equivalences[scaffold_id] = scaffold_id
                elif rename_chromosomes and scaffold_id in self.chromosome_dict:
                    self.confrenamed = True
                    self.scaffolds[self.chromosome_dict[scaffold_id]] = Scaffold(self.chromosome_dict[scaffold_id], str(record.seq).upper(), scaffold_id, description=desc)
                    self.equivalences[scaffold_id] = self.chromosome_dict[scaffold_id]
                else:
                    self.scaffolds[scaffold_id] = Scaffold(scaffold_id, str(record.seq).upper(), original_name=raw_id, description=desc)
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
                else:
                    self.chromosome_names.add(scaffold.name)
            else:
                self.non_chromosomal_scaffolds = True
                self.scaffold_names.add(scaffold.name)

            if update_scaffolds:
                scaffold.update()

        if self.mitochondria or self.chloroplast:
            self.organelles = True

        if self.dapmod:
            self.suffix += "_dapmod"
        
        if self.confrenamed:
            self.suffix += "_confrenamed"

        if not self.non_chromosomal_scaffolds:
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

    def rename_features_dap(self, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_genome_dir: bool = False, subfolder: bool = False, subfolder_name: str = "out_genomes", extension=".fasta", return_equivalences:bool=False, export:bool=False, 
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

            for scaffold in self.scaffolds.values():
                if scaffold.chromosome:
                    if not scaffold.mitochondria and not scaffold.chloroplast:
                        if not scaffold.number:
                            use_chromosome_count = True
                        else:
                            if scaffold.number in chromosome_numbers:
                                use_chromosome_count = True
                            chromosome_numbers.append(scaffold.number)

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
                else:
                    scaffold_count += 1
                    scaffold.number = "{:08d}".format(scaffold_count)
                    self.equivalences[scaffold.original_name] = f"chr{scaffold.number}"

            new_scaffolds = {}

            for scaffold in self.scaffolds.values():
                previous_name = scaffold.name
                new_name = self.equivalences[scaffold.original_name]
                scaffold.update(new_name=new_name)
                if previous_name != scaffold.name:
                    scaffold.dapmod = True
                new_scaffolds[scaffold.name] = scaffold.copy()

            self.scaffolds = new_scaffolds.copy()
            del new_scaffolds
            
            self.update()

            if export:
                self.export(filepath=filepath, output_dir=output_dir, use_genome_dir=use_genome_dir, subfolder=subfolder, subfolder_name=subfolder_name, filename=filename, extension=extension)

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

    def remove_scaffolds(self, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_genome_dir: bool = False, subfolder: bool = False, subfolder_name: str = "out_genomes", extension=".fasta", export:bool=False, remove_00:bool=True, remove_organelles:bool=False):
        
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

            if remove_organelles:
                self.remove_organelles(export=export, output_dir=output_dir)

            elif export:
                self.export(filepath=filepath, output_dir=output_dir, use_genome_dir=use_genome_dir, subfolder=subfolder, subfolder_name=subfolder_name, filename=filename, extension=extension)

        elif remove_organelles:
            self.remove_organelles(export=export, output_dir=output_dir)

    def remove_organelles(self, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_genome_dir: bool = False, subfolder: bool = False, subfolder_name: str = "out_genomes", extension=".fasta", export:bool=False, remove_mitochondria:bool=True, remove_chloroplast:bool=True):

        new_scaffolds = {}

        for scaffold_id, scaffold in self.scaffolds.items():
            if remove_mitochondria:
                if scaffold.mitochondria:
                    continue
            if remove_chloroplast:
                if scaffold.chloroplast:
                    continue
            new_scaffolds[scaffold_id] = scaffold.copy()

        self.scaffolds = new_scaffolds.copy()
        del new_scaffolds

        self.update()

        if export:
            self.export(filepath=filepath, output_dir=output_dir, use_genome_dir=use_genome_dir, subfolder=subfolder, subfolder_name=subfolder_name, filename=filename, extension=extension)

    def export(self, filepath: str | None = None, output_dir: str | None = None, filename: str | None = None, use_genome_dir: bool = False, subfolder: bool = False, subfolder_name: str = "out_genomes", extension=".fasta", quiet:bool=False, keep_description: bool = False,
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

        final_output_path = self._resolve_output_path(filepath=filepath, output_dir=output_dir, filename=filename, suffix=self.suffix, use_genome_dir=use_genome_dir, subfolder_name=subfolder_name, subfolder=subfolder, extension=extension)

        if not self.scaffolds:
            warnings.warn(f"There was nothing to export for {self.name} genome.")
            return

        if not quiet:
            print(f"Exporting {self.name} genome to {final_output_path}.")
        
        with open(str(final_output_path), "w", encoding="utf-8") as f_out:
            for scaffold in self.scaffolds.values():
                header = scaffold.description if (keep_description and getattr(scaffold, "description", None)) else scaffold.name
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