from aegis.genome import Genome
from aegis.annotation import Annotation
import os

genome = Genome('rice_genome', '/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/rice_genome.fasta')
anno1 = Annotation('/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/Oryza_sativa.IRGSP-1.0.61.gtf', 'rice_annot', genome)

os.system(f'docker run --rm -v /media:/media -w /media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis quay.io/biocontainers/liftoff:1.6.3--pyhdfd78af_0 liftoff /media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/rice_genome.fasta /media/tomslab2/Storage/Antonio/TOMSBioLab/genomes_and_annotation/thale_cress/TAIR10_genome/TAIR10_genome.fasta -g /media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/Araport11_chr1_chr2_subset.gff -o /media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/Araport11_chr1_chr2_subset_liftoff.gff')

anno2 = Annotation('/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/Araport11_chr1_chr2_subset_liftoff.gff', 'arabidopsis_annot', genome)

anno1.detect_gene_overlaps(anno2)

anno1.export_equivalences(output_file='/media/tomslab2/Storage/Antonio/TOMSBioLab/Paper_aegis/equivalences.csv', export_self=True, export_csv=True, return_df=False)
