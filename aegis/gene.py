from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .genome import Genome
    from .transcript import Transcript
    from .hits import OverlapHit
    from .subfeatures import UTR

from .feature import Feature
from .subfeatures import Exon

class Gene(Feature):

    __slots__ = (
        'pseudogene', 'transposable', 'transcripts', 'previous_gene',
        'next_gene', 'synteny_order', 'old_previous_gene', 'old_next_gene',
        'old_synteny_order', 'conserved_synteny', 'coding', 'noncoding',
        'overlaps', 'remove', 'rescue', 'reliable', 'reliable_score',
        'overlap_reliable', 'unrescuable', 'overlap_with_selected_CDS',
        'overlap_with_selected_exon', 'alternative_transcript_rescue',
        'intron_nested', 'intron_nested_fully_contained', 'intron_nested_single',
        'UTR_intron_nested', 'transcriptomic_evidence', 'abinitio_evidence',
        'base_id', 'original_base_id', 'renamed_exons', 'renamed_utrs'
    )

    transcripts:dict[str, Transcript]
    synteny_order:int|None
    old_synteny_order:int|None
    previous_gene:str|None|bool
    next_gene:str|None|bool
    old_previous_gene:str|None|bool
    old_next_gene:str|None|bool
    conserved_synteny:bool|None
    alternative_transcript_rescue:list
    overlaps:dict[str, list[OverlapHit]]|None
    
    def __init__(self, pseudogene:bool, transposable:bool, feature_id:str, 
                 ch:str, source:str, feature:str, strand:str,
                 start:int, end:int, score:str, phase:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end,
                         score, phase, parents, attributes)
        self.pseudogene = pseudogene
        self.transposable = transposable
        # transcripts will be added as {"transcript_id" : transcript_object}
        self.transcripts = {}
        # order within the chromosome
        self.previous_gene = None
        self.next_gene = None
        self.synteny_order = None
        self.old_previous_gene = None
        self.old_next_gene = None
        self.old_synteny_order = None
        self.conserved_synteny = None
        # A gene can have both coding and noncoding transcripts
        self.coding = False
        self.noncoding = False

        self.overlaps = None

        self.remove = False
        self.rescue = False
        self.reliable = False
        self.reliable_score = 0
        self.overlap_reliable = False
        self.unrescuable = False

        self.overlap_with_selected_CDS = False
        self.overlap_with_selected_exon = False
        self.alternative_transcript_rescue = []

        self.intron_nested = False
        self.intron_nested_fully_contained = False
        self.intron_nested_single = False
        self.UTR_intron_nested = False

        self.transcriptomic_evidence = False
        self.abinitio_evidence = False

        self.renamed_exons = False
        self.renamed_utrs = False

        self.obtain_base_id(original=True)


    def update(self, quiet:bool=False):
        self.update_size()
        self.sort_transcripts()
        self.coding = False
        self.noncoding = False

        generated_exons = False
        for t in self.transcripts.values():
            if t.generated_exons:
                generated_exons = True
            if t.coding:
                self.coding = True
            else:
                self.noncoding = True
            if self.strand == ".":
                if t.strand != ".":
                    self.strand = t.strand

        if self.coding:
            best_id = None
            best_cds_size = -1
            best_exon_size = -1
            for t in self.transcripts.values():
                t.main = False
                cds_size = 0
                if t.CDSs != {}:
                    t.determine_main_CDS()
                    for c in t.CDSs.values():
                        if c.main:
                            cds_size = c.size
                exon_size = t.size
                # Select by largest CDS, break ties by largest exon span
                if (cds_size > best_cds_size or (cds_size == best_cds_size and exon_size > best_exon_size)):
                    best_cds_size = cds_size
                    best_exon_size = exon_size
                    best_id = t.id
            if best_id:
                self.transcripts[best_id].main = True

        elif self.noncoding:
            best_id = None
            best_exon_size = -1
            for t in self.transcripts.values():
                t.main = False
                if t.size > best_exon_size:
                    best_exon_size = t.size
                    best_id = t.id
            if best_id:
                self.transcripts[best_id].main = True
        elif not quiet and not self.pseudogene:
            print(f"Warning: gene {self.id} has no transcripts annotated")

        self.homogenise_exon_scores()

        if generated_exons:
            self.rename_exons()

    def obtain_base_id(self, original:bool=False):

        if self.id.endswith("gene"):
            self.base_id = self.id[:-4]
        elif self.id.startswith("gene"):
            self.base_id = self.id[4:]
        else:
            self.base_id = self.id

        self.base_id = self.base_id.strip("_,.-/:;")

        if original:
            self.original_base_id = self.base_id

    def rename(self, count:int, sep:str="_", digits:int=5, prefix:str="", suffix:str="", base_id_as_id:bool=False, remove_point_suffix:bool=False):

        if remove_point_suffix:
            if "." in self.id:
                self.id = self.id.split(".")[0]
                if self.original_id != self.id:
                    self.renamed = True
                    self.obtain_base_id()

        if base_id_as_id:
            if self.base_id != self.id:
                self.id = self.base_id
                if self.original_id != self.id:
                    self.renamed = True

        if prefix:

            if suffix:
                self.id = f"{prefix}{self.ch}g{count:0{digits}d}{sep}{suffix}"
            else:
                self.id = f"{prefix}{self.ch}g{count:0{digits}d}"
    
            if self.original_id != self.id:
                self.renamed = True
                self.obtain_base_id()
            
        if self.renamed:
            self.update_numbering()

    def sort_transcripts(self):
        sorted_transcripts = sorted(self.transcripts.values())
        self.transcripts = {t.id: t for t in sorted_transcripts}

    def homogenise_exon_scores(self):
        all_exons = {}
        for t in self.transcripts.values():
            for e in t.exons:
                tag = (e.start, e.end, e.strand)
                if tag not in all_exons:
                    all_exons[tag] = e.score
                elif e.score != ".":
                    if all_exons[tag] == ".":
                        all_exons[tag] = e.score
                    elif float(e.score) > float(all_exons[tag]):
                        all_exons[tag] = e.score

        for t in self.transcripts.values():
            for e in t.exons:
                e.score = all_exons[(e.start, e.end, e.strand)]

    def clear_UTRs(self):
        for t in self.transcripts.values():
            t.clear_UTRs()
        self.update()

    def combine_transcripts(self, genome:Genome, low_memory:bool=True, respect_non_coding:bool=False, quiet:bool=False):
        """
        Useful for RNA-Seq read counting for transcript variants as "one" gene.
        """
        temp_fts = []
        for t in self.transcripts.values():
            for e in t.exons:
                temp_fts.append(e.copy())
        self.transcripts = {}
        if temp_fts != []:
            still_overlapping_fts = True
            while still_overlapping_fts:
                features_to_remove = []
                features_to_add = []
                for i, tempft1 in enumerate(temp_fts):
                    for j, tempft2 in enumerate(temp_fts):
                        if i != j:
                            small = min(tempft1.start, tempft2.start)
                            large = max(tempft1.end, tempft2.end)

                            overlapping, _ = tempft1.overlap(tempft2)

                            if overlapping:
                                temp = Exon("combined", self.ch, self.source, "exon", self.strand, small, large, self.score, ".", [self.id])
                                add = True
                                # this is to avoid adding a same overlap twice
                                for f in features_to_add:
                                    if temp.equal_coordinates(f):
                                        add = False
                                        break
                                if add:
                                    features_to_add.append(temp)
                                if tempft1 not in features_to_remove:
                                    features_to_remove.append(tempft1)
                                if tempft2 not in features_to_remove:
                                    features_to_remove.append(tempft2)
                for sub_to_add in features_to_add:
                    temp_fts.append(sub_to_add)
                for sub_to_rem in features_to_remove:
                    temp_fts.remove(sub_to_rem)
                if features_to_remove == [] and features_to_add == []:
                    still_overlapping_fts = False
            temp_fts.sort()

            temp_coding_feature = "mRNA"
            if not self.coding:
                temp_coding_feature = "lncRNA"
            self.transcripts[f"{self.id}_t001"] = Transcript(f"{self.id}_t001", self.ch, self.source, temp_coding_feature, self.strand, temp_fts[0].start, temp_fts[-1].end, self.score, ".", [self.id])
            self.transcripts[f"{self.id}_t001"].exons = temp_fts.copy()
            counter = 0
            for e in self.transcripts[f"{self.id}_t001"].exons:
                counter += 1
                e.feature = "exon"
                e.id = f"{self.id}_generated_exon_{counter}"
                e.parents = [f"{self.id}_t001"]

        for t in self.transcripts.values():
            t.update(consider_polycistronic=False, consider_read_utrs=False, quiet=quiet)
            if respect_non_coding:
                if not self.coding:
                    continue
            t.generate_sequence(genome, low_memory)
            t.generate_best_protein(genome)
            t.generate_CDSs_based_on_ORF(low_memory)
            for c in t.CDSs.values():
                c.generate_sequence(genome, low_memory)
            if t.coding_ratio < 0.80:
                t.generate_sequence(genome, low_memory)
                t.generate_best_protein(genome, must_have_stop=False)
                t.generate_CDSs_based_on_ORF(low_memory)
                for c in t.CDSs.values():
                    c.generate_sequence(genome, low_memory)
            t.update(consider_polycistronic=False, consider_read_utrs=False, quiet=quiet)

    def longer_CDS(self, other:Gene):
        for t1 in self.transcripts.values():
            if t1.main:
                for c1 in t1.CDSs.values():
                    if c1.main:
                        for t2 in other.transcripts.values():
                            if t2.main:
                                for c2 in t2.CDSs.values():
                                    if c2.main:
                                        if c1.longer(c2):
                                            return True
                                        else:
                                            return False
                                        
    def compare_protein_blast_hits(self, other:Gene, source_priority:list):
        """
        Method required to deal with the fact that a gene may have several blast hits due to several proteins...
        """
        self_proteins = []
        for t in self.transcripts.values():
            for c in t.CDSs.values():
                self_proteins.append(c.protein)

        other_proteins = []
        for t in other.transcripts.values():
            for c in t.CDSs.values():
                other_proteins.append(c.protein)

        best_self_protein = None
        if self_proteins == []:
            print(f"Warning: {self.id} gene has no associated proteins.")
        elif len(self_proteins) == 1:
            best_self_protein = self_proteins[0]
        else:
            best_self_protein = self_proteins[0]
            for n, p in enumerate(self_proteins):
                if n > 0:
                    current_best = best_self_protein.compare_blast_hits(p, source_priority)
                    if not current_best:
                        best_self_protein = p

        best_other_protein = None
        if other_proteins == []:
            print(f"Warning: {other.id} gene has no associated proteins.")
        elif len(other_proteins) == 1:
            best_other_protein = other_proteins[0]
        else:
            best_other_protein = other_proteins[0]
            for n, p in enumerate(other_proteins):
                if n > 0:
                    current_best = best_other_protein.compare_blast_hits(p, source_priority)
                    if not current_best:
                        best_other_protein = p
        
        if best_self_protein != None and best_other_protein != None:
            query_best = best_self_protein.compare_blast_hits(best_other_protein, source_priority)

            return query_best
        
        elif best_other_protein == None and best_other_protein == None:
            print(f"Warning: {self.id} or {other.id} genes have no associated proteins.")
            
            return None
        
        elif best_other_protein == None:

            return True

        elif best_self_protein == None:
            
            return False

    def get_main_CDS_range(self):
        for t in self.transcripts.values():
            if t.main:
                for c in t.CDSs.values():
                    if c.main:
                        return c.start, c.end
                break

        return None

    def rename_exons(self, base_id:str="", sep:str="_", digits:int=3, keep_numbering:bool=False, keep_existing_ids_if_derived_from_base_id:bool=False, name_exons_independently_for_each_transcript:bool=False):
        """
        Deals with gene level naming of exons, making sure shared exons receive the same numbering, also considering strand orientation.
        """

        if base_id == "":
            base_id = self.base_id

        rename = False
        # as long as any of the exons does not contain the base id, rename is triggered
        if keep_existing_ids_if_derived_from_base_id:
            for t in self.transcripts.values():
                for e in t.exons:
                    if base_id not in e.id:
                        rename = True
        else:
            rename = True

        if rename:

            if name_exons_independently_for_each_transcript:
                for t in self.transcripts.values():
                    t.rename_exons(base_id=self.id, sep=sep, digits=digits, keep_numbering=keep_numbering, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id)
                    if t.renamed_exons:
                        self.renamed_exons = True

            else:
                
                exon_names: dict[tuple[int, int, str, str], str] = {}

                exon_d:dict[str, list[Exon]] = {"+":[], "-":[], ".": []}

                for t in self.transcripts.values():
                    for e in t.exons:
                        exon_d[t.strand].append(e)
                for exons in exon_d.values():
                    exons.sort()

                x = 0

                for e in exon_d["+"]:
                    key = (e.start, e.end, e.ch, e.strand)
                    if key not in exon_names:
                        if keep_numbering and e.id_number != None:
                            exon_names[key] = f"{base_id}{sep}e{e.id_number:0{digits}d}"
                        else:
                            x += 1
                            exon_names[key] = f"{base_id}{sep}e{x:0{digits}d}"


                rev_exon_names = {}
                for e in exon_d["-"]:

                    key = (e.start, e.end, e.ch, e.strand)
                    if key not in rev_exon_names:
                        if keep_numbering and e.id_number != None:
                            rev_exon_names[key] = str(e.id_number)
                        else:
                            rev_exon_names[key] = ""

                x += len(rev_exon_names) + 1
                total = x

                for key in rev_exon_names:
                    if rev_exon_names[key] == "":
                        x -= 1
                        exon_names[key] = f"{base_id}{sep}e{x:0{digits}d}"
                    else:
                        exon_names[key] = f"{base_id}{sep}e{rev_exon_names[key]:0{digits}d}"

                x = total

                for e in exon_d["."]:
                    key = (e.start, e.end, e.ch, e.strand)
                    if key not in exon_names:
                        if keep_numbering and e.id_number != None:
                            exon_names[key] = f"{base_id}{sep}e{e.id_number:0{digits}d}"
                        else:
                            x += 1
                            exon_names[key] = f"{base_id}{sep}e{x:0{digits}d}"

                for t in self.transcripts.values():
                    for e in t.exons:
                        key = (e.start, e.end, e.ch, e.strand)
                        e.id = exon_names[key]

                        if e.original_id != e.id:
                            self.renamed_exons = True
                            t.renamed_exons = True
                            e.update_numbering()



    def rename_utrs(self, base_id:str="", sep:str="_", digits:int=3, keep_numbering:bool=False, keep_existing_ids_if_derived_from_base_id:bool=False, name_exons_independently_for_each_transcript:bool=False):
        """
        Deals with gene level naming of utrs, making sure shared utrs receive the same numbering, also considering strand orientation.
        """

        if base_id == "":
            base_id = self.base_id

        rename = False
        # as long as any of the exons does not contain the base id, rename is triggered
        if keep_existing_ids_if_derived_from_base_id:
            for t in self.transcripts.values():
                for c in t.CDSs.values():
                    for u in c.UTRs:
                        if base_id not in u.id:
                            rename = True
        else:
            rename = True

        if rename:

            if name_exons_independently_for_each_transcript:
                for t in self.transcripts.values():
                    t.rename_utrs(base_id=self.id, sep=sep, digits=digits, keep_numbering=keep_numbering, keep_existing_ids_if_derived_from_base_id=keep_existing_ids_if_derived_from_base_id)
                    if t.renamed_utrs:
                        self.renamed_utrs = True

            else:
                
                utr_names: dict[tuple[int, int, str, str], str] = {}

                utr_d:dict[str, list[UTR]] = {"+":[], "-":[], ".": []}

                for t in self.transcripts.values():
                    for c in t.CDSs.values():
                        for u in c.UTRs:
                            utr_d[t.strand].append(u)
                for utrs in utr_d.values():
                    utrs.sort()

                x = 0

                for u in utr_d["+"]:
                    key = (u.start, u.end, u.ch, u.strand)
                    if key not in utr_names:
                        if keep_numbering and u.id_number != None:
                            utr_names[key] = f"{base_id}{sep}u{u.id_number:0{digits}d}"
                        else:
                            x += 1
                            utr_names[key] = f"{base_id}{sep}u{x:0{digits}d}"

                rev_utr_names = {}
                for u in utr_d["-"]:

                    
                    key = (u.start, u.end, u.ch, u.strand)
                    if key not in rev_utr_names:
                        if keep_numbering and u.id_number != None:
                            rev_utr_names[key] = str(u.id_number)
                        else:
                            rev_utr_names[key] = ""

                x += len(rev_utr_names) + 1
                total = x

                for key in rev_utr_names:
                    if rev_utr_names[key] == "":
                        x -= 1
                        utr_names[key] = f"{base_id}{sep}u{x:0{digits}d}"
                    else:
                        utr_names[key] = f"{base_id}{sep}u{rev_utr_names[key]:0{digits}d}"

                x = total

                for u in utr_d["."]:
                    key = (u.start, u.end, u.ch, u.strand)
                    if key not in utr_names:
                        if keep_numbering and u.id_number != None:
                            utr_names[key] = f"{base_id}{sep}u{u.id_number:0{digits}d}"
                        else:
                            x += 1
                            utr_names[key] = f"{base_id}{sep}u{x:0{digits}d}"

                for t in self.transcripts.values():
                    for c in t.CDSs.values():
                        for u in c.UTRs:
                            key = (u.start, u.end, u.ch, u.strand)
                            u.id = utr_names[key]

                            if u.original_id != u.id:
                                self.renamed_utrs = True
                                t.renamed_utrs = True
                                u.update_numbering()

    def collapse_subfeatures(self, exons:bool=True, CDSs:bool=True):

        collapsed_exons = False
        collapsed_CDS_segments = False
        for t in self.transcripts.values():
            if exons:
                t.collapse_exons()
                if t.collapsed_exons == True:
                    collapsed_exons = True
            if CDSs:
                t.collapse_CDS_segments()
                if t.collapsed_CDS_segments:
                    collapsed_CDS_segments = True

        if collapsed_exons or collapsed_CDS_segments:
            if collapsed_exons:
                self.rename_exons()
            if collapsed_CDS_segments:
                self.rename_utrs()
            self.update()

    def print_gff(self, clean:bool=False, names:bool=False, symbols:bool=False, aliases:bool=False, symbols_as_description:bool=False, extra_attributes:bool=False, print_empty_attributes:bool=False):

        temp_attributes = [f"ID={self.id}"]
        if print_empty_attributes:
            if self.symbols is None:
                self.symbols = []
            if self.names is None:
                self.names = []
            if self.aliases is None:
                self.aliases = []
            if self.misc_attributes is None:
                self.misc_attributes = []

        if symbols and self.symbols is not None:
            symbol_string = ",".join(self.symbols)
            if symbol_string or print_empty_attributes:
                temp_attributes.append(f"Symbol={symbol_string}")

        if symbols_as_description and self.symbols is not None:
            symbol_string = ",".join(self.symbols)
            if symbol_string or print_empty_attributes:
                temp_attributes.append(f"Description={symbol_string}")

        if names and self.names is not None:
            name_string = ",".join(self.names)
            if name_string or print_empty_attributes:
                temp_attributes.append(f"Name={name_string}")

        if aliases and self.aliases is not None:
            alias_string = ",".join(self.aliases)
            if alias_string or print_empty_attributes:
                temp_attributes.append(f"Alias={alias_string}")

        if not clean and self.misc_attributes is not None:
            temp_attributes.extend(self.misc_attributes)

        if self.pseudogene:
            temp_attributes.append(f"Pseudogene={self.pseudogene}")
        if self.transposable:
            temp_attributes.append(f"Transposable={self.transposable}")

        if extra_attributes:

            temp_attributes.append(f"reliable_score={self.reliable_score}")
            temp_attributes.append(f"remove={self.remove}")
            temp_attributes.append(f"rescue={self.rescue}")
            blasts = []
            if self.blast_hits is not None:
                for b in self.blast_hits:
                    blasts.append(f"{b.source}_{b.score}")
            blasts = ",".join(blasts)
            if blasts:
                temp_attributes.append(f"blasts={blasts}")
            alternative_transcript_rescue = ",".join(self.alternative_transcript_rescue)
            if alternative_transcript_rescue:
                temp_attributes.append(f"alternative_transcript_rescue={alternative_transcript_rescue}")
            overlaps = []
            if self.overlaps is not None:
                for o in self.overlaps["self"]:
                    if o.score >= 5:
                        overlaps.append(o.id)
            overlaps = ",".join(overlaps)
            if overlaps:
                temp_attributes.append(f"CDS_orientated_overlaps={overlaps}")

            gene_masked_fraction = self.masked_fraction
            transcript_masked_fraction = 0
            CDS_masked_fraction = 0
            gene_GC_content = self.gc_content
            transcript_GC_content = 0
            CDS_GC_content = 0

            for t in self.transcripts.values():
                if t.main:
                    transcript_masked_fraction = t.masked_fraction
                    transcript_GC_content = t.gc_content
                    for c in t.CDSs.values():
                        if c.main:
                            CDS_masked_fraction = c.masked_fraction
                            CDS_GC_content = c.gc_content

            temp_attributes.append(f"gene_masked_fraction={gene_masked_fraction}")
            temp_attributes.append(f"transcript_masked_fraction={transcript_masked_fraction}")
            temp_attributes.append(f"CDS_masked_fraction={CDS_masked_fraction}")
            temp_attributes.append(f"gene_GC_content={gene_GC_content}")
            temp_attributes.append(f"transcript_GC_content={transcript_GC_content}")
            temp_attributes.append(f"CDS_GC_content={CDS_GC_content}")
            temp_attributes.append(f"intron_nested={self.intron_nested}")
            temp_attributes.append(f"intron_nested_fully_contained={self.intron_nested_fully_contained}")
            temp_attributes.append(f"intron_nested_single={self.intron_nested_single}")
            temp_attributes.append(f"intron_UTR_nested={self.UTR_intron_nested}")

        attribute_string = ";".join(temp_attributes)

        return(f"{self.ch}\t{self.source}\t{self.feature}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t{attribute_string}\n")


    def __str__(self):
        if self.symbols:
            return f"{self.id}: {self.symbols}"
        elif self.names:
            return f"{self.id}: {self.names}"
        else:
            return f"{self.id}"
