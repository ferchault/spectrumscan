#!/usr/bin/env python
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# %%
FEATURENAMES = "OO OC OH CC CH O7 O6 O5 O4 O3 O2 =O Aromatic C=C=C (C=C-C)5 (C=C-C)4 (C=C-C)3 (C=C-C)2 (C=C-C)1 C=C".split()

# read and merge folders A...ZE
df = pd.read_csv("data/FTCOUNT", sep=" ", names="feature value state count".split())
df = df.groupby("feature value state".split()).sum().reset_index()
# give columns for states
df = df.pivot_table(
    index="feature value".split(), columns="state", values="count"
).reset_index()
df = df.fillna(0)

# %%
df["total"] = df.B + df.M + df.S
df.B /= df.total
df.M /= df.total
df.S /= df.total
df.M *= 100

# %%
df

# %%

q = df.sort_values(by="M", ascending=False)
xs = range(len(q.M.values))
# plt.fill_between(xs, q.M.values*0, q.M.values*100, alpha=0.3, label="Matching")
plt.plot(q.M.values * 100, label="Matching")
# plt.fill_between(xs, q.M.values*100, (q.M.values+q.S.values)*100, alpha=0.3, label="Stable")
# plt.fill_between(xs, (q.M.values+q.S.values)*100, q.M.values*0+100, alpha=0.1, label="Possible", color="grey")
avg = (df.M.values * df.total.values).sum() / df.total.sum() * 100
deviation = 0.3
plt.axhline(
    avg * (1 + deviation), label="average + %d%%" % (deviation * 100), color="C1"
)
plt.axhline(
    avg * (1 - deviation), label="average - %d%%" % (deviation * 100), color="C2"
)
plt.ylabel("Matching spectrum [%]")
plt.xlabel("Sorted feature index")
plt.legend(frameon=False)

# %%
plt.figure(figsize=(8, 4))
above = q[(q.M.values > avg * (1 + deviation) / 100)]
below = q[(q.M.values < avg * (1 - deviation) / 100)]
labels = []
for feature, value in zip(above.feature.values, above.value.values):
    labels.append(f"{FEATURENAMES[feature]}: {int(value)}")
for feature, value in zip(below.feature.values, below.value.values):
    labels.append(f"{FEATURENAMES[feature]}: {int(value)}")

values = np.zeros(len(labels))
values[: len(above.M.values)] = above.M.values
plt.bar(
    labels, values * 100, label="More likely", alpha=1, facecolor="none", edgecolor="C0"
)
plt.bar(labels, values * 100, alpha=0.2, facecolor="C0", edgecolor="none")

values = np.zeros(len(labels))
values[len(above.M.values) :] = below.M.values
plt.bar(
    labels, values * 100, label="Less likely", alpha=1, facecolor="none", edgecolor="C1"
)
plt.bar(labels, values * 100, alpha=0.2, facecolor="C1", edgecolor="none")
plt.xticks(rotation=90)
plt.legend(frameon=False)
plt.ylabel("Matching spectrum [%]")
plt.xlabel("Feature")


# %%
plt.rc("font", size=14)
f, axs = plt.subplots(1, 4, figsize=(15, 4), sharey=True)
axBonds, axChains, axDouble, axConjugated = axs

# Bonds
markers = "v>v<o"
THISPANEL = "OO OC OH CC CH".split()
xs, ys, labels = [], [], []
for fidx, feature in enumerate(THISPANEL):
    featureidx = FEATURENAMES.index(feature)
    s = df.query("feature == @featureidx").sort_values("value")
    if feature in "OC CH".split():
        s = s.query("value > 0")
    if feature in "CC".split():
        s = s.query("value > 1")
    ys += list(s.M.values)
    labels += [f"{feature}: {int(_)}" for _ in s.value.values]
    axBonds.plot(s.value.values, s.M.values, label=feature, marker=markers[fidx])
xs = list(range(len(ys)))
axBonds.legend(frameon=False, loc="upper left", ncol=3, columnspacing=0.3)
axBonds.set_xlabel("Frequency / Molecule")
axBonds.set_title("Bond frequencies")
axBonds.set_ylim(0, 50)
axBonds.set_ylabel("Share matching spectrum [%]")

# O chains
markers = "ovs>v<^"
THISPANEL = "O2 O3 O4 O5 O6 O7".split()
xs, ys, labels = [], [], []
for fidx, feature in enumerate(THISPANEL):
    featureidx = FEATURENAMES.index(feature)
    s = df.query("feature == @featureidx").sort_values("value")
    ys += list(s.M.values)
    labels += [f"{feature}: {int(_)}" for _ in s.value.values]
    axChains.plot(
        s.value.values, s.M.values, label=f"$n$={feature[-1:]}", marker=markers[fidx]
    )
xs = list(range(len(ys)))
axChains.legend(frameon=False, loc="lower left", ncol=2, columnspacing=0.3)
axChains.set_title("Oxygen-Oxygen chains\nof length $n$")
axChains.set_xlabel("Frequency / Molecule")

# double bond features
markers = "ovs>"
THISPANEL = "=O C=C C=C=C Aromatic".split()
for fidx, feature in enumerate(THISPANEL):
    featureidx = FEATURENAMES.index(feature)
    s = df.query("feature == @featureidx").sort_values("value")
    if len(s) == 0:
        featureidx = -FEATURENAMES[::-1].index(feature) - 1
        s = df.query("feature == @featureidx").sort_values("value")
    ys += list(s.M.values)
    labels += [f"{feature}: {int(_)}" for _ in s.value.values]
    axDouble.plot(s.value.values, s.M.values, label=feature, marker=markers[fidx])
xs = list(range(len(ys)))
axDouble.legend(frameon=False, ncol=2, columnspacing=0.3, loc="upper left")
axDouble.set_title("Double bond features")
axDouble.set_xlabel("Frequency / Molecule")

markers = "o+vsp"
THISPANEL = "(C=C-C)5 (C=C-C)4 (C=C-C)3 (C=C-C)2 (C=C-C)1".split()[::-1]
for fidx, feature in enumerate(THISPANEL):
    featureidx = FEATURENAMES.index(feature)
    s = df.query("feature == @featureidx").sort_values("value")
    if len(s) == 0:
        featureidx = -FEATURENAMES[::-1].index(feature) - 1
        s = df.query("feature == @featureidx").sort_values("value")
    ys += list(s.M.values)
    labels += [f"{int(_)}" for _ in s.value.values]
    axConjugated.plot(s.value.values, s.M.values, alpha=0.5, color=f"C{fidx}")
    axConjugated.plot(
        s.value.values,
        s.M.values,
        label=f"$n$={feature[-1:]}",
        marker=markers[fidx],
        markersize=8,
        alpha=1,
        lw=0,
        color=f"C{fidx}",
    )
xs = list(range(len(ys)))
axConjugated.legend(
    frameon=False, ncol=3, columnspacing=0.0, handletextpad=0.1, loc="upper left"
)
axConjugated.set_title("Conjugated double\nbond chains of length $n$")
axConjugated.set_xlabel("Frequency / Molecule")
plt.subplots_adjust(hspace=0, wspace=0)
# %%
