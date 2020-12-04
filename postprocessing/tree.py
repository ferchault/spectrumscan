#!/usr/bin/env python
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import string

# %%
FEATURENAMES = "OO OC OH CC CH O7 O6 O5 O4 O3 O2 =O Aromatic C=C (C=C-C)1 (C=C-C)2 (C=C-C)3 (C=C-C)4 (C=C-C)5 C=C=C C-O-C'".split()
df = pd.read_csv(
    "data/fullcounts",
    sep="\s+",
    names=["count"] + FEATURENAMES + "cat1 cat2 cat3 cat4".split(),
)

# region
def imbalance(df, colname):
    a = df[df[colname] > 0]
    b = df[df[colname] == 0]
    return np.abs(a.sum()["count"] - b.sum()["count"]), a, b


def choose_next(df):
    imbalances = {_: imbalance(df, _) for _ in FEATURENAMES}
    best_feature = sorted(imbalances.keys(), key=lambda _: imbalances[_][0])[0]
    left = imbalances[best_feature][1]
    right = imbalances[best_feature][2]
    return best_feature, imbalances[best_feature], left, right


# region
def explain(df, label=""):
    suf = "".join(
        random.choice(string.ascii_uppercase + string.digits) for _ in range(10)
    )
    total_molecules = df.sum()["count"]
    if total_molecules == 0:
        return
    broken = df[df["cat1"] == "B"].sum()["count"]
    percentage = 1 - broken / total_molecules
    matching = [
        df[df[f"cat{cat}"] == "M"].sum()["count"] / total_molecules
        for cat in (1, 2, 3, 4)
    ]
    grouppattern = f"""
    <g
       id="g2359{suf}">
       <text
         xml:space="preserve"
         id="text941-4-7{suf}2"
         style="font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;font-size:4.23333px;line-height:1.25;font-family:'Linux Libertine';-inkscape-font-specification:'Linux Libertine';letter-spacing:0px;word-spacing:0px;white-space:pre;shape-inside:url(#rect943-4-9);"
         transform="translate(-33.020172,-87.690285)"><tspan
           x="38.419922"
           y="73.338984"><tspan
             style="shape-inside:url(#rect943-4-9)">{label}</tspan></tspan></text>
      <text
         xml:space="preserve"
         id="text941-4-7{suf}"
         style="font-style:normal;font-variant:normal;font-weight:normal;font-stretch:normal;font-size:4.23333px;line-height:1.25;font-family:'Linux Libertine';-inkscape-font-specification:'Linux Libertine';letter-spacing:0px;word-spacing:0px;white-space:pre;shape-inside:url(#rect943-4-9);"
         transform="translate(-33.020172,-87.690285)"><tspan
           x="38.419922"
           y="73.338984"><tspan
             style="shape-inside:url(#rect943-4-9)">{total_molecules:,}</tspan></tspan></text>
      <path
         style="fill:none;stroke:#000000;stroke-width:0.264583px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
         d="m 0.486041,-13.813886 v 8.084905"
         id="path947-45-6{suf}" />
      <path
         id="path947-4-8-6{suf}"
         d="m 29.752054,-13.813886 v 8.084905"
         style="fill:none;stroke:#000000;stroke-width:0.264583px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1" />
      <rect
         style="fill:#0260bb;fill-opacity:1;stroke:none;stroke-width:0.352777;stroke-linecap:round;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1"
         id="rect966-49-7{suf}"
         width="{29*matching[0]}"
         height="1.3421344"
         x="0.62443095"
         y="-13.270615" />
      <rect
         y="-11.385208"
         x="0.62443101"
         height="1.3421344"
         width="{29*matching[1]}"
         id="rect966-3-6-7{suf}"
         style="fill:#0260bb;fill-opacity:1;stroke:none;stroke-width:0.352777;stroke-linecap:round;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1" />
      <rect
         y="-9.4998016"
         x="0.62443101"
         height="1.3421344"
         width="{29*matching[2]}"
         id="rect966-4-3-9{suf}"
         style="fill:#0260bb;fill-opacity:1;stroke:none;stroke-width:0.352777;stroke-linecap:round;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1" />
      <rect
         y="-7.6143951"
         x="0.62443101"
         height="1.3421344"
         width="{29*matching[3]}"
         id="rect966-7-2-6{suf}"
         style="fill:#0260bb;fill-opacity:1;stroke:none;stroke-width:0.352777;stroke-linecap:round;stroke-miterlimit:4;stroke-dasharray:none;stroke-opacity:1" />
      <path
         style="fill:none;stroke:#ff0000;stroke-width:0.264583px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1"
         d="m {percentage*29},-13.813886 v 8.084905"
         id="path947-4-9-4-1{suf}" />
    </g>"""
    return grouppattern


# region
f, i, l, r = choose_next(df)
print("#META", f)
# print(explain(l))
# print(explain(r))

#%%
def nested(max_level, l_in, r_in):
    if max_level == 0:
        return
    f, i, l_A, l_B = choose_next(l_in)
    if l_A.sum()["count"] < l_B.sum()["count"]:
        l_out = l_A
        label = f"has {f}"
    else:
        l_out = l_B
        label = f"no {f}"
    print("#META", max_level, "left column", f, l_out.sum()["count"])
    print(explain(l_out, label))
    f, i, r_A, r_B = choose_next(r_in)
    if r_A.sum()["count"] < r_B.sum()["count"]:
        r_out = r_A
        label = f"has {f}"
    else:
        r_out = r_B
        label = f"no {f}"
    print("#META", max_level, "right column", f, r_out.sum()["count"])
    print(explain(r_out, label))
    nested(max_level - 1, l_out, r_out)


nested(9, l, r)
