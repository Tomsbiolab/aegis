class OverlapHit():
    def __init__(self, ID, origin, orientation, gene_query_percent,
                 gene_target_percent, exons_in_both, exon_query_percent,
                 exon_target_percent, CDSs_in_both, CDS_query_percent, 
                 CDS_target_percent, protein_query_percent, protein_target_percent, target_synteny_conserved):
        self.full_exon_overlaps = 0
        self.full_protein_overlaps = 0
        self.full_CDS_overlaps = 0
        self.full_UTR_overlaps = 0
        self.score = 1
        self.id = ID
        self.origin = origin
        if gene_query_percent is not None:
            self.gene_query_percent = round(gene_query_percent, 1)
        else:
            self.gene_query_percent = gene_query_percent
        if gene_target_percent is not None:
            self.gene_target_percent = round(gene_target_percent, 1)
        else:
            self.gene_target_percent = gene_target_percent
        if exon_query_percent is not None:
            self.exon_query_percent = round(exon_query_percent, 1)
        else:
            self.exon_query_percent = exon_query_percent
        if exon_target_percent is not None:
            self.exon_target_percent = round(exon_target_percent, 1)
        else:
            self.exon_target_percent = exon_target_percent
        if CDS_query_percent is not None:
            self.CDS_query_percent = round(CDS_query_percent, 1)
        else:
            self.CDS_query_percent = CDS_query_percent
        if CDS_target_percent is not None:
            self.CDS_target_percent = round(CDS_target_percent, 1)
        else:
            self.CDS_target_percent = CDS_target_percent
        if protein_query_percent is not None:
            self.protein_query_percent = round(protein_query_percent, 1)
        else:
            self.protein_query_percent = protein_query_percent
        if protein_target_percent is not None:
            self.protein_target_percent = round(protein_target_percent, 1)
        else:
            self.protein_target_percent = protein_target_percent
        self.exons_in_both = exons_in_both
        self.CDSs_in_both = CDSs_in_both
        self.orientation = orientation
        self.min_gene_percent = min(self.gene_query_percent, self.gene_target_percent)
        if self.exons_in_both:
            self.min_exon_percent = min(self.exon_query_percent, self.exon_target_percent)
        else:
            self.min_exon_percent = None
        if CDSs_in_both:
            self.min_CDS_percent = min(self.CDS_query_percent, self.CDS_target_percent)
        else:
            self.min_CDS_percent = None

        if self.protein_query_percent != None:
            self.min_protein_percent = min(self.protein_query_percent, self.protein_target_percent)
        self.target_synteny_conserved = target_synteny_conserved
        
        """
        Classifying the overlaps into categories based on a reliability score.
        Some things must be considered, if either gene has no exons annotated
        or if they are coding genes (i.e. CDSs annotated).
        """
        if self.orientation:
            # where either of the query and target don't have CDSs
            if not self.CDSs_in_both:
                # where either of the query and target don't have exons
                if not self.exons_in_both:
                    if self.min_gene_percent >= 100:
                        self.score = 9                        
                    elif self.min_gene_percent >= 90:
                        self.score = 8
                    elif self.min_gene_percent >= 70:
                        self.score = 7
                    elif self.min_gene_percent >= 50:
                        self.score = 6                
                    elif self.min_gene_percent >= 30:
                        self.score = 5
                    elif self.min_gene_percent >= 10:
                        self.score = 4
                    else:
                        self.score = 3
                # Exons in both
                else:
                    if self.min_exon_percent == 0:
                        self.score = 2
                    elif self.min_exon_percent >= 100:
                        self.score = 10
                    elif self.min_exon_percent >= 90:
                        self.score = 9
                    elif self.min_exon_percent >= 70:
                        self.score = 8
                    elif self.min_exon_percent >= 50:
                        self.score = 7                
                    elif self.min_exon_percent >= 30:
                        self.score = 6    
                    elif self.min_exon_percent >= 10:
                        self.score = 5                             
                    else:
                        self.score = 4

            # coding query vs coding target
            else:
                if self.min_CDS_percent == 0:
                    self.score = 3
                elif self.min_CDS_percent >= 100:
                    self.score = 11                        
                elif self.min_CDS_percent >= 90:
                    self.score = 10
                elif self.min_CDS_percent >= 70:
                    self.score = 9
                elif self.min_CDS_percent >= 50:
                    self.score = 8                
                elif self.min_CDS_percent >= 30:
                    self.score = 7          
                elif self.min_CDS_percent >= 10:
                    self.score = 6
                else:
                    self.score = 5
        else:
            # where either of the query and target don't have CDSs
            if not self.CDSs_in_both:
                # where either of the query and target don't have exons
                if not self.exons_in_both:
                    if self.min_gene_percent >= 100:
                        self.antiscore = 9                        
                    elif self.min_gene_percent >= 90:
                        self.antiscore = 8
                    elif self.min_gene_percent >= 70:
                        self.antiscore = 7
                    elif self.min_gene_percent >= 50:
                        self.antiscore = 6                
                    elif self.min_gene_percent >= 30:
                        self.antiscore = 5
                    elif self.min_gene_percent >= 10:
                        self.antiscore = 4
                    else:
                        self.antiscore = 3
                # Exons in both
                else:
                    if self.min_exon_percent == 0:
                        self.antiscore = 2
                    elif self.min_exon_percent >= 100:
                        self.antiscore = 10
                    elif self.min_exon_percent >= 90:
                        self.antiscore = 9
                    elif self.min_exon_percent >= 70:
                        self.antiscore = 8
                    elif self.min_exon_percent >= 50:
                        self.antiscore = 7                
                    elif self.min_exon_percent >= 30:
                        self.antiscore = 6    
                    elif self.min_exon_percent >= 10:
                        self.antiscore = 5                             
                    else:
                        self.antiscore = 4

            # coding query vs coding target
            else:
                if self.min_CDS_percent == 0:
                    self.antiscore = 3
                elif self.min_CDS_percent >= 100:
                    self.antiscore = 11                        
                elif self.min_CDS_percent >= 90:
                    self.antiscore = 10
                elif self.min_CDS_percent >= 70:
                    self.antiscore = 9
                elif self.min_CDS_percent >= 50:
                    self.antiscore = 8                
                elif self.min_CDS_percent >= 30:
                    self.antiscore = 7          
                elif self.min_CDS_percent >= 10:
                    self.antiscore = 6
                else:
                    self.antiscore = 5

class BlastHit():
    def __init__(self, source, score, evalue):
        self.source = source
        self.score = score
        self.evalue = evalue