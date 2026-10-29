from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

from helper_functions import jukes_cantor

aligned_seq = Path("data/p1-sars-cov-2-variants.fasta")
ref_fasta = Path("data/p1-NC_045512-ref.fasta")

ref_seq = str(SeqIO.read(ref_fasta, "fasta").seq)


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

for record, date in zip(records, dates):
    days.append((date - start_date).days)
    genome_dist.append(jukes_cantor(ref_seq, str(record.seq)))

days = np.array(days)
genome_dist = np.array(genome_dist)

# linear regression
slope, intercept = np.polyfit(days, genome_dist, 1)

sars_cov_2_per_genome_per_day = slope
sars_cov_2_per_nt_per_day = slope / len(ref_seq)

print(sars_cov_2_per_genome_per_day)
print(sars_cov_2_per_nt_per_day)

# plot
x_fit = np.linspace(days.min(), days.max(), 100)
plt.figure(figsize=(7, 5))
plt.scatter(days, genome_dist, alpha=0.6, label="Genomes")
plt.plot(x_fit, slope * x_fit + intercept, color="red",
         label=f"Fit: {slope:.3f} mutations/genome/day")
plt.xlabel("Elasped Time (Days)")
plt.ylabel("JC-Corrected Genetic Distance")
plt.legend()
plt.savefig("problem1.svg")
plt.show()
