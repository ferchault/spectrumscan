#!/usr/bin/env python
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
        16062.12376,
        16022.67847,
        14294.97459,
        13537.62494,
        11707.36329,
        8788.411529,
        6153.465882,
        3234.514118,
        2193.158353,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
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
sigma = mae / np.sqrt(2 / np.pi)

# Definition of spectrum
def spectrum(values):
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


# limit specifications
upper_limit = spectrum((0.145, 1, 1, 0.5, 0, 0))
upper_limit = np.maximum(upper_limit, experiment_short + 300)
# upper_limit = np.maximum(upper_limit, experiment_short * 2)
upper_limit[:10] += 4000

lower_limit = experiment_short
lower_limit = np.minimum(lower_limit, experiment_short - 4000)
lower_limit = np.maximum(lower_limit, 0)

# figure
f, axs = plt.subplots(1, 2)
spectrum, schematic = axs

spectrum.fill_between(
    grid_short, lower_limit, upper_limit, alpha=0.2, label="Tolerance"
)
spectrum.plot(grid_short, lower_limit, color="C0")
spectrum.plot(grid_short, upper_limit, color="C0")
spectrum.plot(
    grid_short, experiment_short, "o-", color="C1", label="Experiment", markersize=3
)
spectrum.legend(frameon=False)
spectrum.set_xlabel("Energy [Ha]")
spectrum.set_ylabel("Absorbance [a.u.]")
spectrum.set_yticklabels([])
spectrum.set_yticks([])

schematic.set_xlabel("S1 [Ha]")
schematic.set_ylabel("S2 [Ha]")
xs = np.linspace(0.05, 0.2, 40)
xv, yv = np.meshgrid(xs, xs)
center = (max(xs) + min(xs)) / 2
zs = np.exp(-0.5 * ((xv - center) / sigma[3]) ** 2) + np.exp(
    -0.5 * ((yv - center) / sigma[3]) ** 2
)
schematic.contourf(xv, yv, zs)
plt.scatter(center, center, color="red")
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
