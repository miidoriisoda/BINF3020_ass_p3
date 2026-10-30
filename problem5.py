import json, collections
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import SeqIO

# ---------- load ----------
def read_fasta(path):
    recs, h, s = [], None, []
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            if h is not None: recs.append((h, "".join(s)))
            h, s = line[1:], []
        else:
            s.append(line)
    recs.append((h, "".join(s)))
    return recs

recs = read_fasta("data/p5-aus-spike-proteins.fasta")
variants = json.load(open("data/variants.json"))
L = len(recs[0][1])

ids   = [h.split("|")[2] for h, _ in recs]
dates = [h.split("|")[1] for h, _ in recs]
seqs  = [s for _, s in recs]

# ---------- reference spike ----------
ref_seq = str(SeqIO.read("data/p1-NC_045512-ref.fasta", "fasta").seq)
ref = np.array(list(ref_seq.upper()))

# ---------- classify ----------
def classify(seq):
    # all differences vs reference (ignore X)
    if all(a == b or a == "X" for a, b in zip(seq, ref)):
        return "NO_MUTATIONS"
    best, best_frac = None, 0.0
    # only consider non-synonymous spike proteins
    for name, v in variants.items():
        muts = [m for m in v["nonsynonymous"] if m["gene"] == "S"]
        hits = sum(seq[m["pos"] - 1] == m["right"] for m in muts)
        frac = hits / len(muts)
        if frac > best_frac:
            best, best_frac = name, frac
    return best if best_frac >= 0.51 else "UNKNOWN"

labels = [classify(s) for s in seqs]
pd.DataFrame({"id": ids, "variant": labels}).to_csv(
    "problem5-classification.csv", header=False, index=False)

# ---------- monthly distribution ----------
df = pd.DataFrame({"date": pd.to_datetime(dates), "variant": labels})
counts = (df.groupby([pd.Grouper(key="date", freq="MS"), "variant"])
            .size().unstack(fill_value=0))
# months with no sequences skipped
counts = counts[counts.sum(axis=1) > 0]
pct = counts.div(counts.sum(axis=1), axis=0) * 100

order = [c for c in pct.columns if c not in ("UNKNOWN", "NO_MUTATIONS")]
order = ["NO_MUTATIONS"] * ("NO_MUTATIONS" in pct.columns) + order + ["UNKNOWN"] * ("UNKNOWN" in pct.columns)
pct = pct[order]

cmap = plt.get_cmap("tab20")
colors = []
for i, c in enumerate(order):
    colors.append("#bdbdbd" if c == "UNKNOWN" else "#f0f0f0" if c == "NO_MUTATIONS" else cmap(i % 20))

fig, ax = plt.subplots(figsize=(13, 5.5))
ax.stackplot(pct.index, pct.T.values, labels=order, colors=colors, alpha=0.95)
ax.set_ylim(0, 100); ax.set_xlim(pct.index.min(), pct.index.max())
ax.set_ylabel("Percentage Distribution of Variant")
ax.set_xlabel("Time")
ax.set_title("SARS-CoV-2 Spike Variants in Australia Over Time")
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
fig.savefig("problem5.svg", format="svg")

print(pd.Series(labels).value_counts())