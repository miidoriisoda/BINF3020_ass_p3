import re
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO

gene_locations = {
    'S': (21462, 25284),
    'E': (26144, 26372),
    'M': (26422, 27091),
    'N': (28173, 29433)
}

alpha_variants  = ["EPI_ISL_3072041","EPI_ISL_1010703","EPI_ISL_1121977","EPI_ISL_1061036","EPI_ISL_1184096",
                   "EPI_ISL_1184504","EPI_ISL_3568432","EPI_ISL_3568427","EPI_ISL_1184097","EPI_ISL_1447289",
                   "EPI_ISL_1615593","EPI_ISL_1660402","EPI_ISL_1660401","EPI_ISL_2404987","EPI_ISL_2462429", 
                   "EPI_ISL_593696"]

delta_variants =  ["EPI_ISL_1904460", "EPI_ISL_1904466", "EPI_ISL_1904469", "EPI_ISL_1904454", "EPI_ISL_1756027", 
                   "EPI_ISL_2404986", "EPI_ISL_2107448", "EPI_ISL_2107529", "EPI_ISL_2650457", "EPI_ISL_2716414", 
                   "EPI_ISL_2633525", "EPI_ISL_2786933", "EPI_ISL_2811809", "EPI_ISL_2695796", "EPI_ISL_2650470", 
                   "EPI_ISL_2811992", "EPI_ISL_2695794", "EPI_ISL_2695799", "EPI_ISL_3281261", "EPI_ISL_3186692",
                   "EPI_ISL_3186471", "EPI_ISL_3186636", "EPI_ISL_3186262", "EPI_ISL_3186487", "EPI_ISL_3280996", 
                   "EPI_ISL_3305513", "EPI_ISL_3072100", "EPI_ISL_6757593", "EPI_ISL_3426420", "EPI_ISL_3539677", 
                   "EPI_ISL_4055434", "EPI_ISL_3834282", "EPI_ISL_3643716", "EPI_ISL_3501929", "EPI_ISL_5535311", 
                   "EPI_ISL_3502476", "EPI_ISL_3502311", "EPI_ISL_3834174", "EPI_ISL_9964573", "EPI_ISL_4747964",
                   "EPI_ISL_4205585", "EPI_ISL_4552276", "EPI_ISL_4205719", "EPI_ISL_4205539", "EPI_ISL_5689721", 
                   "EPI_ISL_4205197", "EPI_ISL_4552130", "EPI_ISL_4747905", "EPI_ISL_5347813", "EPI_ISL_5534640", 
                   "EPI_ISL_5936432", "EPI_ISL_5778590", "EPI_ISL_5534654", "EPI_ISL_5778467", "EPI_ISL_5778492", 
                   "EPI_ISL_6252077", "EPI_ISL_5535134", "EPI_ISL_5305653", "EPI_ISL_7264600", "EPI_ISL_6585779", 
                   "EPI_ISL_6586228", "EPI_ISL_6251927", "EPI_ISL_6585956", "EPI_ISL_7264849", "EPI_ISL_6586173", 
                   "EPI_ISL_6586136", "EPI_ISL_7264591", "EPI_ISL_6757382", "EPI_ISL_7621131", "EPI_ISL_7621166", 
                   "EPI_ISL_8901286", "EPI_ISL_7748984", "EPI_ISL_7621531", "EPI_ISL_7621692", "EPI_ISL_7621124", 
                   "EPI_ISL_7748623", "EPI_ISL_7748687", "EPI_ISL_7621974", "EPI_ISL_8901160", "EPI_ISL_10116719", 
                   "EPI_ISL_8901318", "EPI_ISL_8901337", "EPI_ISL_10116652", "EPI_ISL_8900986", "EPI_ISL_8901124",
                   "EPI_ISL_8901168", "EPI_ISL_8901097", "EPI_ISL_10116931"]

# --- load alignment ---
records = list(SeqIO.parse("data/p1-sars-cov-2-variants.fasta", "fasta"))
seqs = {}
for r in records:
    m = re.search(r"EPI_ISL_\d+", r.id)
    seqs[m.group(0) if m else r.id] = str(r.seq).upper()

ref_seq = str(SeqIO.read("data/p1-NC_045512-ref.fasta", "fasta").seq)
ref = np.array(list(ref_seq.upper()))

def to_array(ids):
    return np.array([list(seqs[i]) for i in ids])

alpha = to_array(alpha_variants)
delta = to_array(delta_variants)

# --- mutation frequency per position ---
# note: gaps = mismatches, N != mismatches
def mutation_freq(arr):
    mism = (arr != ref) & (arr != "N")
    return mism.mean(axis=0)

alpha_freq = mutation_freq(alpha)
delta_freq = mutation_freq(delta)

# --- plot ---
colors = {"S": "#e07a5f", "E": "#3d405b", "M": "#81b29a", "N": "#f2cc8f"}
fig, axes = plt.subplots(2, 1, figsize=(14, 7), sharex=True, sharey=True)
start = 20000
x = np.arange(len(ref))

for ax, freq, name in zip(axes, [alpha_freq, delta_freq],
                          ["Alpha (20I, V1)", "Delta (21A)"]):
    ax.plot(x[start:], freq[start:], color="black", lw=0.9)
    for g, (s, e) in gene_locations.items():
        ax.axvspan(s, e, color=colors[g], alpha=0.35)
        ax.text((s + e) / 2, 1.04, g, ha="center", va="bottom",
                fontweight="bold", color=colors[g])
    ax.set_title(name, loc="left", fontsize=12)
    ax.set_ylabel("Mutation occurrence")
    ax.set_ylim(0, 1.15)
    ax.set_xlim(start, 29600)
    ax.spines[["top", "right"]].set_visible(False)

axes[-1].set_xlabel("Sequence Index (Reference NC_045512.2)")
fig.suptitle("Average Mutation Occurrence of SARS-CoV-2 Genome", fontsize=14)
fig.tight_layout()
fig.savefig("problem3.svg", format="svg")

# --- vital mutations (>0.5) in S ---
s0, s1 = gene_locations["S"]
def vital_sites(freq):
    return {i for i in range(s0, s1 + 1) if freq[i] > 0.5}

def consensus_base(arr, i):
    vals, counts = np.unique(arr[:, i], return_counts=True)
    return vals[counts.argmax()]

vital_alpha, vital_delta = vital_sites(alpha_freq), vital_sites(delta_freq)
mutations_in_both_variants = sorted(int(i) for i in vital_alpha & vital_delta)

same = [f"{ref[i]}{i}{consensus_base(alpha, i)}"
        for i in mutations_in_both_variants
        if consensus_base(alpha, i) == consensus_base(delta, i)]
same_mutation = same[0]

print(mutations_in_both_variants, same_mutation)
print(sorted(vital_alpha)); print(sorted(vital_delta))