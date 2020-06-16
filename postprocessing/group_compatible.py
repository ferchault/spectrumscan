#%%
import matplotlib.pyplot as plt
import numpy as np

#%%
def read_hist():
    lines = open("../SO.count").readlines()
    values, weights = [], []
    for line in lines:
        count, what = line.strip().split()
        if what == "--":
            continue
        values.append(float(what))
        weights.append(int(count))
    return values, weights


cs = np.array(read_hist())

# %%
sigma1 = np.exp(-0.5) ** 6 * max(cs[0])
plt.hist(cs[0], weights=cs[1], bins=50, log=True, histtype="step")
plt.axvline(sigma1, color="red")
cs[1][cs[0] > sigma1].sum()


# %%
sigma1 * 10000

# %%
sigma1

# %%
