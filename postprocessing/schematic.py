#!/usr/bin/env python
#%%
import matplotlib.pyplot as plt
import numpy as np
import scipy.interpolate as sci

plt.rc("font", size=14)

# meta parameters
LINEWIDTH = 0.3
EXPERIMENT_UNCERTAIN_SCALE = 10

# reference data
grid_short = np.array(
    [
        0.138912671079212,
        0.136826895237182,
        0.134802828739590,
        0.132837772927060,
        0.130929184235579,
        0.129074663212412,
        0.127271944452462,
        0.125518887366340,
        0.123813467701037,
        0.122153769742578,
        0.120537979137517,
        0.118964376276714,
        0.117431330190674,
        0.115937292910894,
        0.114480794256235,
        0.113060437007398,
        0.111674892436229,
        0.110322896159762,
        0.109003244291822,
        0.107714789867569,
        0.106456439518648,
        0.105227150378710,
        0.104025927200871,
        0.102851819670387,
        0.101703919897280,
        0.100581360075014,
        0.099483310292536,
        0.098408976488081,
        0.097357598534149,
    ]
)

experiment_short = np.array(
    [
        [16062.12376, 9009.305174, 9103.97388, 8220.399292],
        [16022.67847, 8788.411527, 8835.74588, 7534.051174],
        [14294.97458, 7431.49341, 7210.599763, 5837.903528],
        [13537.62494, 7605.052704, 6910.815527, 5664.344234],
        [11707.36329, 6816.146822, 5159.444469, 5159.444469],
        [8788.411527, 5254.113175, 3187.179764, 3297.626587],
        [6153.465881, 5222.55694, 1956.486588, 2224.714588],
        [3234.514117, 3644.745175, 0, 757.3496468],
        [2193.158352, 3692.079528, 0, 0],
        [0, 2445.608235, 0, 0],
        [0, 2840.061176, 0, 0],
        [0, 2571.833176, 0, 0],
        [0, 2335.161411, 0, 0],
        [0, 2185.269293, 0, 0],
        [0, 2966.286117, 0, 0],
        [0, 3376.517175, 0, 0],
        [0, 3360.739058, 0, 0],
        [0, 4717.657175, 0, 0],
        [0, 4812.325881, 0, 0],
        [0, 5064.775763, 0, 0],
        [0, 4638.766587, 0, 0],
        [0, 7344.713763, 0, 0],
        [0, 8393.958586, 0, 0],
        [0, 7431.49341, 0, 0],
        [0, 8835.74588, 0, 0],
        [0, 9837.65635, 0, 0],
        [0, 10019.1047, 0, 0],
        [0, 10618.67317, 0, 0],
        [0, 10886.90117, 0, 0],
    ]
)
mae = np.array(
    (
        0.009367729738238392,
        0.007644026369760024,
        0.006634207502480805,
        0.03514044031325359,
        0.037518177682891826,
        0.035671701408550395,
    )
)
# Definition of spectrum
def fspectrum(values):
    sigma_inv = 1 / (LINEWIDTH * 0.036749304951208)
    const = 16836.70285
    dec = 0.4342944819032518
    return (
        const
        * dec
        * (0.05746523219)
        * sigma_inv
        * (
            values[3] * np.exp(-0.5 * ((grid_short - values[0]) * sigma_inv) ** 2)
            + values[4] * np.exp(-0.5 * ((grid_short - values[1]) * sigma_inv) ** 2)
            + values[5] * np.exp(-0.5 * ((grid_short - values[2]) * sigma_inv) ** 2)
        )
    )


def get_limits(EIDX):
    sigma = mae / np.sqrt(2 / np.pi)

    # limit specifications
    upper_limit = fspectrum((0.145, 1, 1, 0.5, 0, 0))
    factor = experiment_short[0, EIDX] / max(upper_limit)
    upper_limit = upper_limit * factor
    upper_limit = np.maximum(upper_limit, experiment_short[:, EIDX] + 300)
    # upper_limit = np.maximum(upper_limit, experiment_short * 2)
    upper_limit[:10] += experiment_short[0, EIDX] / 4

    lower_limit = experiment_short[:, EIDX]
    lower_limit = np.minimum(
        lower_limit, experiment_short[:, EIDX] - experiment_short[0, EIDX] / 4
    )
    lower_limit = np.maximum(lower_limit, 0)

    if EIDX == 1:
        upper_limit[12:] = 200000
        lower_limit[12:] = 0

    return lower_limit, upper_limit


# figure
# f, axs = plt.subplots(1, 5, figsize=(12, 5))
f = plt.figure(constrained_layout=True)
gs = f.add_gridspec(4, 2, hspace=0.0, wspace=0)
meta_axs = f.add_subplot(gs[:, 0])
meta_axs.tick_params(
    labelcolor="none", top=False, bottom=False, left=False, right=False
)
meta_axs.spines["right"].set_visible(False)
meta_axs.spines["top"].set_visible(False)
meta_axs.spines["left"].set_visible(False)
meta_axs.spines["bottom"].set_visible(False)
axs = [f.add_subplot(gs[_, 0]) for _ in range(4)]
schematic = f.add_subplot(gs[:, 1])

for eidx in range(4):
    spectrum = axs[eidx]
    lower_limit, upper_limit = get_limits(eidx)

    spectrum.fill_between(
        grid_short, lower_limit, upper_limit, alpha=0.2, label="Tolerance"
    )
    spectrum.plot(grid_short, lower_limit, color="C0")
    spectrum.plot(grid_short, upper_limit, color="C0")
    spectrum.plot(
        grid_short,
        experiment_short[:, eidx],
        "o-",
        color="C1",
        label="Experiment",
        markersize=3,
    )

    if eidx == 3:
        spectrum.set_xlabel("Energy [Ha]")
    if eidx == 2:
        spectrum.legend(frameon=False)
    meta_axs.set_ylabel("Absorbance [a.u.]", labelpad=-20)
    if eidx != 3:
        spectrum.set_xticklabels([])
        spectrum.set_xticks([])
    spectrum.set_yticklabels([])
    spectrum.set_yticks([])
    spectrum.set_ylim(0, 20000)

schematic.set_xlabel("S1 [Ha]")
schematic.set_ylabel("S2 [Ha]")
xs = np.linspace(0.05, 0.2, 40)
xv, yv = np.meshgrid(xs, xs)
center = (max(xs) + min(xs)) / 2
zs = np.exp(-0.5 * ((xv - center) / sigma[3]) ** 2) + np.exp(
    -0.5 * ((yv - center) / sigma[3]) ** 2
)
schematic.contourf(xv, yv, zs)
schematic.scatter(center, center, color="red")
schematic.set_yticklabels([])
schematic.set_yticks([])
schematic.set_xticklabels([])
schematic.set_xticks([])

# region
angle = np.linspace(0, 2 * np.pi)
radius = 0.03 * np.ones(len(angle))

np.random.seed(0)
xs = np.linspace(0, 2 * np.pi, 10)
ys = np.array((0, 1, 1, 0)) * 0.01
ys = np.random.normal(size=10) * 0.005
ys[-1] = ys[0]
spl = sci.CubicSpline(xs, ys, bc_type="periodic")

radius += spl(angle)

center = 0.1
xs = center + np.sin(angle) * radius
ys = center + np.cos(angle) * radius - 0.01
schematic.plot(xs, ys, color="white")
schematic.fill(xs, ys, color="white", alpha=0.4)
schematic.annotate(xy=(center - 0.02, center - 0.03), s="$\Omega$", fontsize=20)

f.align_xlabels(axs)

# export
plt.savefig("schematic.pdf", bbox_inches="tight")

# %%
