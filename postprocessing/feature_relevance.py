#!/usr/bin/env python
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %%
FEATURENAMES = "OO OC OH CC CH O7 O6 O5 O4 O3 O2 =O Aromatic C=C (C=C-C)1 (C=C-C)2 (C=C-C)3 (C=C-C)4 (C=C-C)5 C=C=C C-O-C'".split()
df = pd.read_csv(
    "data/fullcounts",
    sep="\s+",
    names=["count"] + FEATURENAMES + "cat1 cat2 cat3 cat4".split(),
)

#%%
def do_plot(spectrumid, ylim):
    plt.rc("font", size=14)
    f, axs = plt.subplots(1, 4, figsize=(15, 4.5), sharey=True)
    axBonds, axChains, axDouble, axConjugated = axs

    def _plot_into(axs, markers, features, title, loc):
        for fidx, feature in enumerate(features.split()):
            s = df.groupby(f"{feature} cat{spectrumid}".split()).sum().reset_index()
            s = s.pivot(
                index=feature, columns=f"cat{spectrumid}", values="count"
            ).reset_index()
            ys = np.array(list(s.M.values / (s.M.values + s.S.values + s.B.values)))
            if "$" in title:
                label = feature[-1]
            else:
                label = feature
            axs.plot(s[feature].values, ys * 100, label=label, marker=markers[fidx])

        axs.set_xlabel("Frequency / Molecule")
        axs.set_title(title)
        axs.legend(frameon=False, loc="upper left", ncol=2, columnspacing=0.3)

    _plot_into(
        axBonds, "ovs>v<^", "OO OC OH CC CH C-O-C'", "Bond frequencies", "upper left"
    )
    _plot_into(
        axChains,
        "ovs>v<^",
        "O2 O3 O4 O5 O6 O7",
        "Oxygen-Oxygen chains\nof length $n$",
        "lower left",
    )
    _plot_into(
        axDouble, "ovs>", "=O C=C C=C=C Aromatic", "Double bond features", "upper left"
    )
    _plot_into(
        axConjugated,
        "o+vsp",
        "(C=C-C)1 (C=C-C)2 (C=C-C)3 (C=C-C)4 (C=C-C)5",
        "Conjugated double\nbond chains of length $n$",
        "upper left",
    )
    axBonds.set_ylabel("Share matching spectrum [%]")
    axBonds.set_ylim(0, ylim)
    plt.subplots_adjust(hspace=0, wspace=0)
    plt.savefig(f"features-{spectrumid}.pdf", bbox_inches="tight")


do_plot(1, 50)
do_plot(2, 80)
do_plot(3, 60)
do_plot(4, 60)

#%%
s = df
total = s["count"].sum()
s = s.query("cat1 != 'B'")
stable = s["count"].sum()
s = s.query("OH != 0")
hasOH = s["count"].sum()
s = s.query("O7 == 0 & O6 == 0 &O5 == 0 &O4 == 0 &O3 == 0 ")
nolongchains = s["count"].sum()
s = s[s["C-O-C'"] > 0]
hasbridge = s["count"].sum()
s = s.query("Aromatic>0")
isaromatic = s["count"].sum()
total, stable, hasOH, nolongchains, hasbridge, isaromatic

# region
