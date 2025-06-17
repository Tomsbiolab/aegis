
# Standard library packages
import os
import time

# Import Numpy, Pandas and Seaborn
import numpy as np
import pandas as pd

# Import Biopython tools for running local BLASTX
from Bio.Blast.Applications import NcbiblastpCommandline

from config.cannabis import *

cores = 8

start = time.time()

# by default 1-e05
evalue = 0.00001
coverage = 30
identity = 30
max_hsps = 1

aegis_output = f"{root}aegis_output/"
blast_path = f"{aegis_output}blasts/"

os.system(f"mkdir -p {aegis_output}")
os.system(f"mkdir -p {blast_path}")
os.system(f"mkdir -p {blast_path}rbh/")

features_path_1 = f"{path}aegis_output/features/"

from config.PN40024 import *

features_path_2 = f"{path}aegis_output/features/"

protein_file_tag = "_proteins_g_id_main.fasta"

annotations_1 = ["CBDRx_on_CBDRx", "Finola_on_Finola", "JamaicanLion_on_JamaicanLion", "PurpleKush_on_PurpleKush"]

annotations_2 = ["5.1_on_T2T_ref"]

fastas_1 = []
for annotation in annotations_1:
    fastas_1.append(f"{features_path_1}{annotation}{protein_file_tag}")

fastas_2 = []
for annotation in annotations_2:
    fastas_2.append(f"{features_path_2}{annotation}{protein_file_tag}")

checked_pairs = []

for x1, f1 in enumerate(fastas_1):
    for x2, f2 in enumerate(fastas_2):
        if f1 == f2:
            continue
        pair = [f1, f2]
        pair.sort()
        if pair in checked_pairs:
            continue
        checked_pairs.append(pair)

        pair_start = time.time()

        s1 = os.path.abspath(f1)
        s2 = os.path.abspath(f2)

        t1 = annotations_1[x1]
        t2 = annotations_2[x2]

        print(f"Starting Blasts, RBH, and RBBHs with {t1} and {t2} proteins.")

        fwd_out = os.path.abspath(f'{blast_path}fwd_{t1}_against_{t2}.csv')
        rev_out = os.path.abspath(f'{blast_path}rev_{t2}_against_{t1}.csv')

        fwd_single_out = os.path.abspath(f'{blast_path}fwd_{t1}_against_{t2}_single.csv')
        rev_single_out = os.path.abspath(f'{blast_path}rev_{t2}_against_{t1}_single.csv')

        rbh_out = os.path.abspath(f'{blast_path}rbh/rbh_{t1}_against_{t2}.csv')
        rbbh_out = os.path.abspath(f'{blast_path}rbh/rbbh_{t1}_against_{t2}.csv')

        if not os.path.isfile(fwd_out):
            fwd_blastp = NcbiblastpCommandline(query=s1, subject=s2, out=fwd_out,
                                            outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                            evalue=evalue, num_threads=cores, max_hsps=max_hsps, qcov_hsp_perc=coverage)
            print("Forward: %s" % fwd_blastp)
            fwd_stdout, fwd_stderr = fwd_blastp()

        if not os.path.isfile(rev_out):
            rev_blastp = NcbiblastpCommandline(query=s2, subject=s1, out=rev_out,
                                            outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                            evalue=evalue, num_threads=cores, max_hsps=max_hsps, qcov_hsp_perc=coverage)
            print("Reverse: %s" % rev_blastp)
            rev_stdout, rev_stderr = rev_blastp()

        if not os.path.isfile(fwd_single_out):
            fwd_blastp_single = NcbiblastpCommandline(query=s1, subject=s2, out=fwd_single_out,
                                            outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                            max_target_seqs=1, num_threads=cores, max_hsps=max_hsps)
            print("Forward: %s" % fwd_blastp_single)
            fwd_stdout_single, fwd_stderr_single = fwd_blastp_single()

        if not os.path.isfile(rev_single_out):
            rev_blastp_single = NcbiblastpCommandline(query=s2, subject=s1, out=rev_single_out,
                                            outfmt="6 qseqid sseqid pident qcovs qlen slen length bitscore evalue",
                                            max_target_seqs=1, num_threads=cores, max_hsps=max_hsps)
            print("Reverse: %s" % rev_blastp_single)
            rev_stdout_single, rev_stderr_single = rev_blastp_single()

        if not os.path.isfile(rbh_out):

            fwd_results = pd.read_csv(fwd_out, sep="\t", header=None)
            rev_results = pd.read_csv(rev_out, sep="\t", header=None)

            headers = ["query", "subject", "identity", "coverage", "qlength", "slength", "alength", "bitscore", "E-value"]
            fwd_results.columns = headers
            rev_results.columns = headers

            fwd_results = fwd_results[(fwd_results["identity"] >= identity)]
            rev_results = rev_results[(rev_results["identity"] >= identity)]


            # Create a new column in both dataframes: normalised bitscore
            fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength
            rev_results['norm_bitscore'] = rev_results.bitscore/rev_results.qlength

            # Create query and subject coverage columns in both dataframes
            fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
            rev_results['qcov'] = rev_results.alength/rev_results.qlength
            fwd_results['scov'] = fwd_results.alength/fwd_results.slength
            rev_results['scov'] = rev_results.alength/rev_results.slength

            # Clip maximum coverage values at 1.0
            fwd_results['qcov'] = fwd_results['qcov'].clip(upper=1)
            rev_results['qcov'] = rev_results['qcov'].clip(upper=1)
            fwd_results['scov'] = fwd_results['scov'].clip(upper=1)
            rev_results['scov'] = rev_results['scov'].clip(upper=1)

            # Merge forward and reverse results
            rbh = pd.merge(fwd_results, rev_results, left_on=['subject', 'query'], right_on=['query', 'subject'], how='inner')

            rbh.to_csv(rbh_out, sep = '\t')

            del rbh

        if not os.path.isfile(rbbh_out):

            fwd_results = pd.read_csv(fwd_single_out, sep="\t", header=None)
            rev_results = pd.read_csv(rev_single_out, sep="\t", header=None)

            headers = ["query", "subject", "identity", "coverage", "qlength", "slength", "alength", "bitscore", "E-value"]
            fwd_results.columns = headers
            rev_results.columns = headers

            # Create a new column in both dataframes: normalised bitscore
            fwd_results['norm_bitscore'] = fwd_results.bitscore/fwd_results.qlength
            rev_results['norm_bitscore'] = rev_results.bitscore/rev_results.qlength

            # Create query and subject coverage columns in both dataframes
            fwd_results['qcov'] = fwd_results.alength/fwd_results.qlength
            rev_results['qcov'] = rev_results.alength/rev_results.qlength
            fwd_results['scov'] = fwd_results.alength/fwd_results.slength
            rev_results['scov'] = rev_results.alength/rev_results.slength

            # Clip maximum coverage values at 1.0
            fwd_results['qcov'] = fwd_results['qcov'].clip(upper=1)
            rev_results['qcov'] = rev_results['qcov'].clip(upper=1)
            fwd_results['scov'] = fwd_results['scov'].clip(upper=1)
            rev_results['scov'] = rev_results['scov'].clip(upper=1)

            # Merge forward and reverse results
            rbbh = pd.merge(fwd_results, rev_results, left_on=['subject', 'query'], right_on=['query', 'subject'], how='inner')

            rbbh.to_csv(rbbh_out, sep = '\t')

            duplicates = rbbh[rbbh.duplicated(subset=['query_x', 'subject_x'], keep=False)]
            if not duplicates.empty:
                print(f"\nWarning: Duplicate rows found based on ['query_x', 'subject_x']: for {f1} and {f2} RBBHs")
                print(duplicates)

            del rbbh

        pair_end = time.time()

        lapse = pair_end - pair_start

        print(f"\nRunning Blasts, RBHs, and RBBHs on {t1} vs {t2} proteins took {round((lapse/60)/60, 1)} hours\n")

end = time.time()
lapse = end - start
print(f"\nRunning RBHs between all these proteins={annotations_1} and {annotations_2} took {round((lapse/60)/60, 1)} hours\n")