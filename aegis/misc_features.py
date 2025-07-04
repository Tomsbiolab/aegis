import copy
from .genefunctions import translate
from .feature import Feature

class Protein():
    def __init__(self, prot_id:str, nucleotides:str, chrom:str, readthrough:str="both"):
        self.id = prot_id
        self.ch = chrom
        self.summary_tag = ""
        self.readthrough = readthrough
        self.blast_hits = []
        # readthrough can be start, end, both or none 
        self.start, self.end_stop, self.early_stop, self.nucleotide_surplus, self.gaps, self.seq, self.coding_start, self.coding_end = translate(nucleotides, readthrough)
        if self.start == "late" or self.start == "absent" or self.end_stop == False or self.nucleotide_surplus or self.gaps:
            self.partial = True
            self.summary_tag += "partial"
        else:
            self.partial = False
        if self.early_stop and self.end_stop:
            self.truncated = True
            if self.summary_tag == "":
                self.summary_tag += "truncated"
            else:
                self.summary_tag += "_truncated"
        else:
            self.truncated = False
        self.size = len(self.seq)

    def copy(self):
        return copy.deepcopy(self)

    def compare_blast_hits(self, other:object, source_priority:list):
        compared = False
        query_best = True
        while not compared:
            for s in source_priority:
                query_evalue = float(2)
                query_bitscore = float(-1)
                target_evalue = float(2)
                target_bitscore = float(-1)
                
                for b in self.blast_hits:
                    if b.source == s:
                        if b.evalue < query_evalue:
                            query_evalue = b.evalue
                        if b.score > query_bitscore:
                            query_bitscore = b.score

                for b in other.blast_hits:
                    if b.source == s:
                        if b.evalue < target_evalue:
                            target_evalue = b.evalue
                        if b.score > target_bitscore:
                            target_bitscore = b.score

                if query_evalue < target_evalue:
                    compared = True
                    break
                elif query_evalue > target_evalue:
                    query_best = False
                    compared = True
                    break
                elif query_bitscore > target_bitscore:
                    compared = True
                    break
                elif query_bitscore < target_bitscore:
                    query_best = False
                    compared = True
                    break
                elif s == source_priority[-1]:
                    compared = True

        return query_best

class Promoter(Feature):
    def __init__(self, promoter_type, feature_id:str, ch:str, source:str, 
                 feature:str, strand:str, start:int, end:int, score:str,
                 phase:str, attributes:str):
        super().__init__(feature_id, ch, source, feature, strand, start, end, score, phase, attributes)
        self.type = promoter_type