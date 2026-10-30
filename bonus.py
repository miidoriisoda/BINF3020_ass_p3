from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

from helper_functions import kimura_two_parameter

aligned_seq = Path("data/p1-sars-cov-2-variants.fasta")
ref_fasta = Path("data/p1-NC_045512-ref.fasta")
ref_seq = str(SeqIO.read(ref_fasta, "fasta").seq)

PURINES = {"A", "G"}

def find_date(description):
    # date is stored in the 3rd column of fasta description
    return pd.to_datetime(description.split("|")[2])


records = list(SeqIO.parse(aligned_seq, "fasta"))
dates = [find_date(r.description) for r in records]
# starting point = earliest collection date
start_date = min(dates)

# include reference itself at the origin (0, 0)
days = [0]
genome_dist = [0.0]

# for record, date in zip(records, dates):
#     days.append((date - start_date).days)
#     genome_dist.append(kimura_two_parameter(ref_seq, str(record.seq)))

# days = np.array(days)
# genome_dist = np.array(genome_dist)

# # linear regression
# slope, intercept = np.polyfit(days, genome_dist, 1)

# sars_cov_2_per_genome_per_day = slope
# sars_cov_2_per_nt_per_day = slope / len(ref_seq)

# print(sars_cov_2_per_genome_per_day)
# print(sars_cov_2_per_nt_per_day)

# --- calculate transversion vs transitions ---
def count_ts_tv(ref, seq):
    ts = tv = 0
    for a, b in zip(ref.upper(), seq.upper()):
        if a in "ACGT" and b in "ACGT" and a != b:
            if (a in PURINES) == (b in PURINES): ts += 1
            else: tv += 1
    return ts, tv

other_seq = []
for r in records:
    other_seq.append(r.seq)
ts_tv = [count_ts_tv(ref_seq, s) for s in other_seq]
mean_ts = sum(t for t, _ in ts_tv) / len(ts_tv)
mean_tv = sum(v for _, v in ts_tv) / len(ts_tv)

print(f"mean transitions: {mean_ts}")
print(f"mean traversions: {mean_tv}")