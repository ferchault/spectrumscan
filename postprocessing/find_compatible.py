#!/usr/bin/env python
# %%
import numpy as np
import numba
import scipy.optimize as sco

# meta parameters
LINEWIDTH = 0.3

# reference data
grid_short = np.array([ 0.138912671079212, 0.136826895237182, 0.134802828739590, 0.132837772927060, 0.130929184235579, 0.129074663212412, 0.127271944452462, 0.125518887366340, 0.123813467701037, 0.122153769742578, 0.120537979137517, 0.118964376276714, 0.117431330190674, 0.115937292910894, 0.114480794256235, 0.113060437007398, 0.111674892436229, 0.110322896159762, 0.109003244291822, 0.107714789867569, 0.106456439518648, 0.105227150378710, 0.104025927200871, 0.102851819670387, 0.101703919897280, 0.100581360075014, 0.099483310292536, 0.098408976488081, 0.097357598534149])
experiment_short = np.array([16062.12376,16022.67847,14294.97459,13537.62494,11707.36329,8788.411529,6153.465882,3234.514118,2193.158353,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0])

# %%
# Definition of spectrum
@numba.jit(nopython=True)
def spectrum(values):
    sigma_inv=1/(LINEWIDTH*0.036749304951208)
    const=16836.70285
    dec=0.4342944819032518
    return const*dec*(0.05746523219) * sigma_inv*(values[3]*np.exp(-0.5* ((grid_short-values[0]) * sigma_inv) ** 2) +values[4]*np.exp(-0.5* ((grid_short-values[1]) * sigma_inv) ** 2) +values[5]*np.exp(-0.5* ((grid_short-values[2]) * sigma_inv) ** 2))

# %%
# Define the "acceptable" range
import matplotlib.pyplot as plt
plt.plot(grid_short, experiment_short, label='Experimental')
upper_limit = spectrum((0.147, 1, 1, 0.5, 0, 0))
upper_limit = np.maximum(upper_limit, experiment_short + 300)
upper_limit = np.maximum(upper_limit, experiment_short * 2)

lower_limit = experiment_short * 0.5
lower_limit = np.minimum(lower_limit, experiment_short - 2000)
lower_limit = np.maximum(lower_limit, 0)
plt.fill_between(grid_short, lower_limit, upper_limit, alpha=0.2, label='Accepted')
plt.legend(loc="upper left")

# %%
@numba.jit(nopython=True)
def gaussian_pdf(mu, sigma, value):
    prefactor = 1/(sigma * np.sqrt(2*np.pi))
    exponent = ((value - mu) / sigma)**2*(-0.5)
    return prefactor * np.exp(exponent)

def likelihood(centers, sigmas, candidate_parameters):
    """ Given a probability distribution for each observable and concrete manifestations, calculate the likelihood of observing the experimental spectrum."""
    prefactor = 1/(sigmas * np.sqrt(2*np.pi))
    exponent = ((candidate_parameters - centers) / sigmas)**2*(-0.5)
    return prefactor * np.exp(exponent)

def band_target(candidate):
    """ Validity prefactor: 1 if within band, 0 otherwise."""

def objective(candidate):
    return band_target(candidate) * likelihood(centers, variances, candidate)

#def optimize_one(s1, s2)

def optimize(value):
    
    sigma = mae/np.sqrt(2/np.pi)
    lower_bounds = np.maximum(value - sigma, 0)
    upper_bounds = value + sigma
    return sco.differential_evolution(residuals, bounds=list(zip(lower_bounds, upper_bounds)), tol=1).x

def monte_carlo_likelihood():
    pass

# %%
dbg = np.array((0.10377130098640919,0.14916887134313583,0.15750830247998238, 0.013235272252773295,0.010341944675863707,0.012630209337737078))
mae = np.array((0.009367729738238392, 0.007644026369760024, 0.006634207502480805, 0.03514044031325359, 0.037518177682891826, 0.035671701408550395))
sigma = mae/np.sqrt(2/np.pi)
success = 0
Nrounds = 100000

osz1 = np.random.normal(dbg[3], sigma[3], Nrounds)
osz1 = osz1[osz1 > 0]
osz2 = np.random.normal(dbg[4], sigma[4], len(osz1))
osz2 = osz2[osz2 > 0]
osz3 = np.random.normal(dbg[5], sigma[5], len(osz2))
osz3 = osz3[osz3 > 0]
print (len(osz1), len(osz2), len(osz3))

# %%
for i in range(100000):
    candidate = np.random.normal(dbg, sigma)
    if min(candidate[3:]) < 0:
        continue
    s = spectrum(candidate)
    rescale = max(s) / max(experiment_short)
    if rescale > 5 or rescale <0.2:
        continue
    s /= rescale
    if np.all(s >= lower_limit) and np.all(s <= upper_limit):
        success += 1
success

# %%
xs = np.linspace(-5, 5, 100)
ys = likelihood(0, 0.1, xs)
plt.plot(xs, ys)

# %%
plt.plot(grid_short, spectrum(dbg))

# %%
sigma

# %%
1/(sigma * np.sqrt(2*np.pi))

# %%
from scipy.special import erf
from scipy.stats import truncnorm
def random_osz(mu, sigma, N):
    larger_than_0 = 1-(1 + erf(-mu/(sigma * np.sqrt(2))))*0.5
    larger_than_0 = int(larger_than_0 * N)
    a, b = (0 - mu) / sigma, (100 - mu) / sigma
    return truncnorm.rvs(a, np.inf, loc=mu, scale=sigma, size=larger_than_0)


plt.hist(random_osz(dbg[4], sigma[4], Nrounds), histtype='step', range=(0, 0.2), bins=100)
osz1 = np.random.normal(dbg[4], sigma[4], Nrounds)
osz1 = osz1[osz1 > 0]
plt.hist(osz1, histtype="step", range=(0, 0.2), bins=100)
import timeit
timeit.timeit('random_osz(dbg[4], sigma[4], Nrounds)', globals=globals(), number=1)


# %%
def numpysample(mu, sigma, N):
    osz1 = np.random.normal(mu, sigma, N)
    osz1 = osz1[osz1 > 0]
    return osz1
timeit.timeit('numpysample(dbg[4], sigma[4], Nrounds)', globals=globals(), number=1)

# %%
len(random_osz(dbg[4], sigma[4], Nrounds))

# %%
len(numpysample(dbg[4], sigma[4], Nrounds))

# %%
