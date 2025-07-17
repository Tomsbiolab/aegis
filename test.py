from aegis.genome import Genome
from aegis.annotation import Annotation

# genome = Genome('rice_genome', '/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/rice_genome.fasta')
# anno = Annotation('rice_annot', '/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/Oryza_sativa.IRGSP-1.0.61.gtf', genome)

# genome = Genome('arabidopsis_genome', '/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/thale_cress/TAIR10_genome/TAIR10_genome.fasta')
# anno = Annotation('arabidopsis_annot', '/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/Araport11_chr1_chr2_subset.gff', genome)
genome = Genome(name="genome", genome_file_path=None)
annotation = Annotation(name="annotation", annot_file_path='/media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/thale_cress/TAIR10_genome/TAIR10_genome.fasta', genome=None)

annotation.self_overlap_genes()

annotation.export_equivalences(output_file='/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/equivalences.csv', export_self=True, export_csv=True, return_df=False)