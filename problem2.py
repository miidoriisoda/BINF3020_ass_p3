from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import SeqIO

from helper_functions import jukes_cantor


def estimate_mutation_rate(fasta_path, include_reference=False):
    records = list(SeqIO.parse(fasta_path, "fasta"))

    # first genome in the file is the reference (earliest known virus)
    ref = records[0]
    ref_seq = str(ref.seq)
    # format="mixed" copes with dates like 2014-03 or 2014 alongside full dates
    ref_date = pd.to_datetime(ref.description.split("|")[2], format="mixed")

    days, dist = [], []
    for record in records[1:]:
        date = pd.to_datetime(record.description.split("|")[2], format="mixed")
        days.append((date - ref_date).days)
        dist.append(jukes_cantor(ref_seq, str(record.seq)))

    days = np.array(days)
    dist = np.array(dist)

    if include_reference:
        days = np.insert(days, 0, 0)
        dist = np.insert(dist, 0, 0.0)

    slope, intercept = np.polyfit(days, dist, 1)
    per_genome = slope
    per_nt = slope / len(ref_seq)  # normalise by alignment length
    return days, dist, slope, intercept, per_genome, per_nt


def plot_fit(days, dist, slope, intercept, title, outfile):
    x = np.linspace(0, days.max(), 100)
    plt.figure(figsize=(7, 5))
    plt.scatter(days, dist, alpha=0.6, label="Genomes")
    plt.plot(x, slope * x + intercept, color="red", label=f"Fit: slope = {slope:.4f}")
    plt.scatter([0], [0], color="black", zorder=3, label="Reference")
    plt.xlabel("Days Since Reference")
    plt.ylabel("JC-Corrected Distance")
    plt.title(title)
    plt.legend()
    plt.savefig(outfile)
    plt.show()


# Ebola
d, g, s, i, ebola_per_genome_per_day, ebola_per_nt_per_day = estimate_mutation_rate(
    Path("data/p2-ebola-viruses.fasta")
)
plot_fit(d, g, s, i, "Zaire ebolavirus", "problem2_ebola.svg")

# Zika
d, g, s, i, zika_per_genome_per_day, zika_per_nt_per_day = estimate_mutation_rate(
    Path("data/p2-zika-viruses.fasta")
)
plot_fit(d, g, s, i, "Zika virus", "problem2_zika.svg")

print(ebola_per_genome_per_day, ebola_per_nt_per_day)
print(zika_per_genome_per_day, zika_per_nt_per_day)