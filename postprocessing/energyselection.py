#!/usr/bin/env python

# %%
basepath = "/mnt/c/Users/guido/data/spectrumscan-40d3089e55ff4be5a902cac1844bc22e"
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

# %%
def get_energies(matching):
    folder = "NON_MATCHING"
    if matching:
        folder = "MATCHING"

    with open(f"{basepath}/{folder}/fingerprints") as fh:
        fingerprints = fh.readlines()

    rows = []
    for fn in glob.glob(f"{basepath}/{folder}/energies*dat"):
        with open(fn) as fh:
            for line in fh:
                try:
                    fpidx, energy = line.strip().split()
                except:
                    continue
                fp = fingerprints[int(fpidx)].strip()
                energy = float(energy)
                rows.append({"fp": fp, "energy": energy})
    return pd.DataFrame(rows)


matchingE = get_energies(True)
nonmatchingE = get_energies(False)

# %%
def get_groups():
    rows = []
    for folder in "MATCHING NON_MATCHING".split():
        with open(f"{basepath}/{folder}/groups") as fh:
            for line in fh:
                count, fp = line.strip().split(maxsplit=1)
                rows.append({"fp": fp, "size": int(count)})
    return pd.DataFrame(rows)


groups = get_groups()

# %%
plt.rc("font", size=12)
s = pd.merge(nonmatchingE.groupby("fp").min().reset_index(), groups)
plt.hist(s.energy, bins=50, histtype="step", label="non-matching spectra")
s = pd.merge(matchingE.groupby("fp").min().reset_index(), groups)
plt.hist(s.energy, bins=50, histtype="step", label="matching spectra")
plt.tick_params(axis="y", which="both", left=False, labelleft=False)
plt.ylabel("Frequency [a.u]")
plt.xlabel("Total energy [Ha]")
plt.legend()
plt.savefig("energies.pdf", bbox_inches="tight")


# %%
for fp in (
    pd.merge(nonmatchingE.groupby("fp").min().reset_index(), groups)
    .sort_values("size")
    .tail(5)
    .fp
):
    print(nonmatchingE.query("fp == @fp").sort_values("energy").head(1))

# %%

# %%
