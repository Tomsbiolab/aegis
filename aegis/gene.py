from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .hits import OverlapHit
    from .subfeatures import UTR

from .feature import Feature
from .subfeatures import Exon
from .other_components import GeneSynteny
from .transcript import Transcript

class Gene(Feature):

    __slots__ = ('pseudogene', 'transposable', 'transcripts', 'noncoding', '_overlaps', '_synteny', '_quality', 'original_base_id', 'renamed_exons', 'renamed_utrs')

    transcripts:dict[str, Transcript]
    alternative_transcript_rescue:list
    _overlaps:dict[str, list[OverlapHit]]|None
    _synteny: GeneSynteny | None
    
    def __init__(self, pseudogene:bool, transposable:bool, feature_id:str, ch:str, source:str, feature:str, strand:str, start:int, end:int, score:str, parents:list[str]=[], attributes:dict={}):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, parents, attributes)
        self.pseudogene = pseudogene
        self.transposable = transposable
        # transcripts will be added as {"transcript_id" : transcript_object}
        self.transcripts = {}
        # A gene can have both coding and noncoding transcripts -> coding attribute is present in Feature already
        self.noncoding = False

        self._overlaps = None

        self.renamed_exons = False
        self.renamed_utrs = False

        self.original_base_id = self.base_id

        self._synteny = None


    def update(self, quiet:bool=False):
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

    def rename(self, count:int, sep:str="_", digits:int=5, prefix:str="", suffix:str="", base_id_as_id:bool=False, remove_point_suffix:bool=False):

        if remove_point_suffix:
            if "." in self.id:
                self.id = self.id.split(".")[0]
                if self.original_id != self.id:
                    self.renamed = True

        if base_id_as_id:
            if self.base_id != self.id:
                self.id = self.base_id
                if self.original_id != self.id:
                    self.renamed = True

        if prefix:

            if suffix:
                suffix = f"{sep}{suffix}"

            self.id = f"{prefix}{self.ch}g{count:0{digits}d}{suffix}"
    
            if self.original_id != self.id:
                self.renamed = True
            
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

    def combine_transcripts(self, respect_non_coding:bool=False, respect_non_combined:bool=False, redetect_CDS:bool=False, quiet:bool=False):
        """
        Useful for RNA-Seq read counting for transcript variants as "one" gene.
        """
        if len(self.transcripts) > 1:

            temp_fts = []
            coordinate_tuples = set()
            for t in self.transcripts.values():
                for e in t.exons:
                    if (e.start, e.end) not in coordinate_tuples:
                        temp_fts.append(e)
                        coordinate_tuples.add((e.start, e.end))

            collapsed_exons = False
            if temp_fts != []:

                temp_fts.sort()

                for t in self.transcripts.values():
                    if t.main:
                        new_transcript = t

                new_transcript.id = f"{self.id}_t001"
                new_transcript.exons = temp_fts

                self.transcripts = {new_transcript.id: new_transcript}

                new_transcript.collapse_exons()
                if new_transcript.collapsed_exons == True:
                    collapsed_exons = True

                if collapsed_exons:
                    self.rename_exons()

                self.update(quiet=quiet)

            for t in self.transcripts.values():

                if not redetect_CDS:
                    continue

                if respect_non_combined:
                    if not collapsed_exons:
                        continue

                if respect_non_coding:
                    if not self.coding:
                        continue

                t.generate_best_protein(quiet=quiet)
                t.update(quiet=quiet)

                if t.coding_ratio < 0.80:
                    t.generate_best_protein(tolerated_stops=1, quiet=quiet)
                t.update(quiet=quiet)

                if t.coding_ratio < 0.80:
                    t.generate_best_protein(must_have_stop=False, quiet=quiet)
                t.update(quiet=quiet)
            
            self.update(quiet=quiet)

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
        self_proteins = [c.protein for t in self.transcripts.values() for c in t.CDSs.values() if c.protein]
        other_proteins = [c.protein for t in other.transcripts.values() for c in t.CDSs.values() if c.protein]

        if not self_proteins and not other_proteins:
            print(f"Warning: {self.id} and {other.id} genes have no associated proteins.")
            return None
        elif not self_proteins:
            print(f"Warning: {self.id} gene has no associated proteins.")
            return False
        elif not other_proteins:
            print(f"Warning: {other.id} gene has no associated proteins.")
            return True

        if not source_priority:
            return True

        best_self_protein = min(self_proteins, key=lambda p: p.get_blast_hits_key(source_priority))
        best_other_protein = min(other_proteins, key=lambda p: p.get_blast_hits_key(source_priority))

        return best_self_protein.compare_blast_hits(best_other_protein, source_priority)

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

    def collapse_subfeatures(self, exons:bool=True, CDSs:bool=True, quiet:bool=False):

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
            self.update(quiet=quiet)

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
            extra = self.quality.get_attributes()
            temp_attributes.extend(extra)

        attribute_string = ";".join(temp_attributes)
        phase = self.phase if self.phase is not None else "."
        return(f"{self.ch}\t{self.source}\t{self.feature}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{phase}\t{attribute_string}\n")

    def __str__(self):
        if self.symbols:
            return f"{self.id}: {self.symbols}"
        elif self.names:
            return f"{self.id}: {self.names}"
        else:
            return f"{self.id}"

    # all of these gene properties are cached-like -> compatible with slots
    @property
    def overlaps(self) -> dict[str, list[OverlapHit]]:
        if self._overlaps is None:
            self._overlaps = {"self": [], "other":[]}
        return self._overlaps

    @property
    def synteny(self) -> GeneSynteny:
        if self._synteny is None:
            self._synteny = GeneSynteny()
        return self._synteny

    @property
    def base_id(self) -> str:

        base = self.id
        
        if base.endswith("gene"):
            base = base[:-4]
        elif base.startswith("gene"):
            base = base[4:]

        return base.strip("_,.-/:;")

    @property
    def alternative_transcript_rescue(self) -> list:
        return self.quality.alternative_transcript_rescue

