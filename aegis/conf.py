
default_features = {}
default_features["gene"] = {"gene", "pseudogene", "transposable_element_gene"}
default_coding_transcripts = {"mRNA"}
default_generic_transcripts = {"transcript", "transcript_region", "primary_transcript", "pseudotranscript", "pseudogenic_transcript", "mRNA_TE_gene"}
default_noncoding_transcripts = {"antisense_lncRNA", "antisense_RNA", "miRNA_primary_transcript", "ncRNA", "lncRNA", "lnc_RNA", "pseudogenic_tRNA", "rRNA", "snoRNA", "snRNA", "tRNA", "pre_miRNA", "tRNA_pseudogene", "SRP_RNA", "RNase_MRP_RNA", "Y_RNA", "YRNA", "scaRNA", "vault_RNA", "telomerase_RNA", "scRNA", "RNase_P_RNA", "RNase_MRP_RNA"}
default_other_transcript_level_features = {"V_gene_segment", "D_gene_segment", "J_gene_segment", "C_gene_segment"}
default_codons = {"start_codon", "stop_codon"}
default_introns = {"intron"}
# Some features are clearly transcript level features but they cannot be
# classed as coding/noncoding just by looking at the name   
default_features["transcript"] = default_generic_transcripts.union(default_coding_transcripts).union(default_noncoding_transcripts).union(default_other_transcript_level_features)
default_features["UTR"] = {"UTR", "three_prime_UTR", "five_prime_UTR", "five_prime_utr", "three_prime_utr"}
default_features["exon"] = {"exon", "pseudogenic_exon"}
default_features["CDS"] = {"CDS", "nucleotide_to_protein_match"}
default_features["miRNA"] = {"miRNA"}

default_subfeatures = default_features["UTR"].union(default_features["exon"]).union(default_features["CDS"]).union(default_codons).union(default_features["miRNA"]).union(default_introns)

default_features_r = {}
for key, values in default_features.items():
    for value in values:
        default_features_r[value] = key
