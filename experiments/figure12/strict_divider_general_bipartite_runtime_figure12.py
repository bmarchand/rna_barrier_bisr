import sys
import matplotlib.pylab as plt
import numpy as np
import seaborn as sns
import json
from rna_barrier_bisr.general_bipartite import random_bipartite
from networkx.algorithms.bipartite.matching import maximum_matching
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats

DEGREE = 5
N_LWB = 4
N_UPB = 41
N_STEP = 6 
N_SEEDS = 200

fname = 'general_bipartite_strict_divider_results_'+str(DEGREE)+'.json'

with open(fname, 'r') as f:
    res_dict = json.load(f)

avg = []
err = []

ks_per_n = {}

for n in range(N_LWB, N_UPB, N_STEP):

    rtimes = []

    ks_per_n[n] = []

    for seed in range(1, N_SEEDS+1, 1):

        key = str(tuple((n,DEGREE,seed)))
        res = res_dict[key]

        rtimes.append(float(res.split('..')[1])) 
        ks_per_n[n].append(int(res.split('..')[0]))

    bsr = bs.bootstrap(np.array(rtimes), stat_func=bs_stats.mean)
    m = bsr.value
    lb = m - bsr.lower_bound
    ub = bsr.upper_bound - m
    avg.append(m)
    err.append((lb,ub))
    print("n, m, lb, ub, bsr.lb, bsr.ub", n, m, lb, ub, bsr.lower_bound, bsr.upper_bound)

fig, (ax1, ax2) = plt.subplots(2, figsize=(8, 4.8))

maxsize = -1

for n in ks_per_n.keys():
    if len(ks_per_n[n]) > maxsize:
        maxsize = len(ks_per_n[n])

for n, li in ks_per_n.items():

    new_li = []
    new_li += li
    for _ in range(len(li), maxsize, 1):
        new_li.append(np.nan)
        
    ks_per_n[n] = new_li

print(ks_per_n)
import pandas as pd
df = pd.DataFrame(data = ks_per_n)#{n:ks_per_n[n] for n in sorted(ks_per_n.keys())})
print(df)

sns.boxplot(data=df, ax=ax2, color='mediumorchid',native_scale=True)

ax1.set_yscale("log")
ax1.set_xticks([])
ax1.set_ylabel("run-time (s, log-scale)")
ax2.set_xlabel("number of vertices")
ax2.set_ylabel(r'$\rho$ value')

ax1.errorbar(list(range(N_LWB, N_UPB, N_STEP)), 
             avg, 
             yerr=np.array(err).transpose(), 
             fmt='o-', 
             ecolor='black',
             barsabove=True,
             label='Divider-Schedule', 
             capsize = 2.)

fig.legend()
fig.savefig("general_bipartite_strict_divider_figure12.pdf")

plt.show()
