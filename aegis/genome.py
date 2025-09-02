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


class Scaffold():
    mitochondria_suffixes = ["m", "M"]
    chloroplast_suffixes = ["c", "C"]
    unknown_chromosome_names = ["chrUn", "chrun", "ChrUn", "Chrun", "chr00", "Chr00", "chr0", "Chr0"]
    organelle_suffixes = mitochondria_suffixes + chloroplast_suffixes
    def __init__(self, name, sequence, original_name:str=""):

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

        if original_name:
            self.original_name = original_name
        else:
            self.original_name = self.name

        self.update()

    def update(self, new_name:str=""):
        if new_name:
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
    def __init__(self, name:str, genome_file_path:str, chromosome_dict:dict={}, rename_chromosomes:bool=False, quiet:bool=False):
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

        # Creating a dictionary with the genome sequence of each chromosome or scaffold, still referred as chromosomes in the code
        self.chromosome_dict = chromosome_dict

        self.equivalences = {}

        count = 0
        with open (self.file, "r", encoding="utf-8") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                scaffold_id = record.id
                count += 1
                if scaffold_id in self.scaffolds:
                    print((f"Error: scaffold feature {scaffold_id} is repeated in {self.name}, genome (file: {self.file})"))
                if self.dapmod:
                    self.scaffolds[scaffold_id] = Scaffold(scaffold_id, str(record.seq).upper(), original_name=f"unknown_dapmod_{count}")
                    self.equivalences[scaffold_id] = scaffold_id
                elif rename_chromosomes and scaffold_id in self.chromosome_dict:
                    self.confrenamed = True
                    self.scaffolds[self.chromosome_dict[scaffold_id]] = Scaffold(self.chromosome_dict[scaffold_id], str(record.seq).upper(), scaffold_id)
                    self.equivalences[scaffold_id] = self.chromosome_dict[scaffold_id]
                else:
                    self.scaffolds[scaffold_id] = Scaffold(scaffold_id, str(record.seq).upper())
                    self.equivalences[scaffold_id] = scaffold_id

        self.update()

        now = time.time()
        lapse = now - start
        if not quiet:
            print(f"\n{self.name} genome chromosomes/scaffolds: {self.features[:30]} ...")
            print(f"\nCreating {self.name} genome object took {round(lapse, 1)} seconds\n")

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

        for scaffold_id, scaffold in self.scaffolds.items():

            if update_scaffolds:
                scaffold.update()

            if not scaffold.name.startswith("chr"):
                self.dapfit = False
            else:
                try:
                    number = int(scaffold.name.split("r")[1])
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
                elif scaffold.chloroplast:
                    self.chloroplast = True
            else:
                self.non_chromosomal_scaffolds = True

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
        
    def export_feature_sizes(self, custom_path:str=""):

        self.update()

        if custom_path == "":
            export_folder = self.path + "genome_feature_sizes/"
        else:
            export_folder = custom_path
        if export_folder[-1] != "/":
            export_folder += "/"
        system(f"mkdir -p {export_folder}")

        tag = f"{self.name}_genome_feature_sizes{self.suffix}.tsv"

        out = ""
        for scaffold_id, scaffold in self.scaffolds.items():
            out += f"{scaffold_id}\t{scaffold.size}\n"
        f_out = open(f"{export_folder}{tag}", "w", encoding="utf-8")
        f_out.write(out)
        f_out.close()
    
    def rename_features_dap(self, output_folder:str="", return_equivalences:bool=False, export:bool=False, chromosome_dict:dict={}):
        """
        Renames scaffolds and chromosomes to become dapfit.
        """

        if not self.dapfit:

            # checking to see if all chromosomes have unique numbers, else they
            # will be all renumbered

            use_chromosome_count = False
            chromosome_numbers = []

            for scaffold_id, scaffold in self.scaffolds.items():
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

            for scaffold_id, scaffold in self.scaffolds.items():

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

            for scaffold_id, scaffold in self.scaffolds.items():
                previous_name = scaffold.name
                scaffold.name = self.equivalences[scaffold.original_name]
                if previous_name != scaffold.name:
                    scaffold.dapmod = True
                new_scaffolds[scaffold.name] = scaffold.copy()

            self.scaffolds = new_scaffolds.copy()
            del new_scaffolds
            
            self.update()

            if export:
                self.export(output_folder=output_folder, file=".fasta")

            if return_equivalences:
                return self.equivalences
            
        else:

            print(f"Warning: {self.name} genome is already fit for DAP-Seq analysis, so it will not be modified.")

    def rename_features_from_dic(self, rename_map: dict) -> dict:

        new_scaffolds = {}
        final_equivalences = {}

        for scaffold_id, scaffold in self.scaffolds.items():
            
            original_name = scaffold.original_name
            new_name = rename_map.get(original_name, scaffold.name)
            scaffold.name = new_name

            final_equivalences[original_name] = new_name

            new_scaffolds[new_name] = copy.copy(scaffold)

        self.scaffolds = new_scaffolds
        self.update()

        return final_equivalences

    def remove_scaffolds(self, output_folder:str="", export:bool=False, chromosome_dict:dict={}, remove_00:bool=True, remove_organelles:bool=False):
        
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
                self.remove_organelles(export=export, output_folder=output_folder)

            elif export:
                self.export(output_folder=output_folder, file=".fasta")

        elif remove_organelles:
            self.remove_organelles(export=export, output_folder=output_folder)

    def remove_organelles(self, output_folder:str="", export:bool=False, remove_mitochondria:bool=True, remove_chloroplast:bool=True):

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
            self.export(output_folder=output_folder, file=".fasta")

    def export(self, output_folder:str="", file:str=".fasta", quiet:bool=False):

        self.update()

        if not output_folder:
            export_folder = self.path + "out_genomes/"
        else:
            export_folder = output_folder
        if export_folder[-1] != "/":
            export_folder += "/"

        system(f"mkdir -p {export_folder}")

        if file == ".fasta":
            file = f"{self.name}{self.suffix}{file}"

        out = ""

        for scaffold in self.scaffolds.values():
            out += f">{scaffold.name}\n{scaffold.seq}\n"

        if out:
            if not quiet:
                print(f"Exporting {self.name} genome to {file}.")
            f_out = open(f"{export_folder}{file}", "w", encoding="utf-8")
            f_out.write(out)
            f_out.close()
        else:
            print(f"Warning: there was nothing to export for {file} genome")

    def copy(self):
        return copy.deepcopy(self)
        
    def extract_peak_sequences(self, output_file_name, DAPseq_output_file, output_folder: str = "", top=600):

        export_folder = Path(output_folder or Path(self.path) / "out_peak_seqs/")
        export_folder.mkdir(parents=True, exist_ok=True)

        try:
            df = pd.read_csv(DAPseq_output_file, delimiter='\t', dtype=str)
            df.dropna(how='all', inplace=True)
        except FileNotFoundError:
            print(f"Error: DAPseq file not found: {DAPseq_output_file}")
            return

        # --- Selection of 'top' peaks ---
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

                peak_position = int(peak_str.split(':')[1])

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

        output_path = export_folder / output_file_name
        with open(output_path, "w", encoding="utf-8") as f_out:
            for header, seq in peaks.items():
                f_out.write(f'>{header}\n')
                f_out.write(f'{textwrap.fill(seq, width=60)}\n')

    def subset(self, chosen_features:set=None, cap:int=2):

        if chosen_features is None:
            chosen_features = set()

        if chosen_features:
            for chosen_feature in chosen_features:
                if chosen_feature not in self.scaffolds:
                    raise ValueError(f"Chosen scaffold/chromosome {chosen_feature} is not in {self.name} genome. Choose from '{self.scaffolds.keys()}'")
            scaffolds_to_remove = set(self.scaffolds) - chosen_features
        else:
            if cap > len(self.scaffolds):
                warnings.warn(f"Cap value {cap} exceeds the number of available scaffolds/chrosomomes ({len(self.scaffolds)}). No features removed in subset genome {self.id}.", category=UserWarning)
                return
            scaffolds_to_remove = set(self.scaffolds) - set(random.sample(list(self.scaffolds), cap))

        if scaffolds_to_remove:
            self.remove_features(scaffolds_to_remove)
        else:
            print(f"No scaffolds/chromosomes removed from {self.id} genome.")

    def remove_features(self, features_to_remove:set):
        for ft in features_to_remove:
            del self.scaffolds[ft]
        self.update()