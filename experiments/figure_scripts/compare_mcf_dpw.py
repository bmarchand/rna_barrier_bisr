import sys
import matplotlib.pylab as plt
import numpy as np
import seaborn as sns
import json
from general_bipartite import random_bipartite, vertex_cover_from_matching
from networkx.algorithms.bipartite.matching import maximum_matching
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats


fname1 = sys.argv[-2]
fname2 = sys.argv[-1]

if fname2.split('_')[3]=='mcf':
    fname1, fname2 = fname2, fname1

DEGREE = 5
THETA = 5
N_LWB = 10
N_UPB = 160
N_STEP = 10 
N_SEEDS = 500
N_MAX_MCF = 90
RNA = True

with open(fname1, 'r') as f:
    res_dict1 = json.load(f)

with open(fname2, 'r') as f:
    res_dict2 = json.load(f)

avg1 = []
avg2 = []

err1 = []
err2 = []

ks_per_n = {}

for n in range(N_LWB, N_UPB, N_STEP):

    rtimes1 = []
    rtimes2 = []

    ks_per_n[n] = []

    for seed in range(1, N_SEEDS+1, 1):

        key = str(tuple((n,DEGREE,seed)))


        try:
            res2 = res_dict2[key]
            print('res', res2)
            rtimes2.append(float(res2.split('..')[1])) 
            try:
                dpw = int(res2.split('..')[2])
            except IndexError:
                continue
        except KeyError:
            print(key, "not in ", fname2)
            continue

        ks_per_n[n].append(dpw)
       
        if n < N_MAX_MCF:
            try: 
                res1 = res_dict1[key]
                rtimes1.append(float(res1.split('..')[1])) 
            except KeyError:    
                continue

    bsr = bs.bootstrap(np.array(rtimes2), stat_func=bs_stats.mean)
    m = bsr.value
    lb = m - bsr.lower_bound
    ub = bsr.upper_bound - m
    avg2.append(m)
    err2.append((lb,ub))
    print("n, m, lb, ub, bsr.lb, bsr.ub", n, m, lb, ub, bsr.lower_bound, bsr.upper_bound)

    if n < N_MAX_MCF:
        bsr = bs.bootstrap(np.array(rtimes1), stat_func=bs_stats.mean)
        m = bsr.value
        lb = m - bsr.lower_bound 
        ub = bsr.upper_bound - m
        avg1.append(m)
        err1.append((lb,ub))
        print("lb,ub,m", lb, m, ub)

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

import pandas as pd
df = pd.DataFrame(data = {n:ks_per_n[n] for n in sorted(ks_per_n.keys())})

sns.boxplot(data=df, ax=ax2, color='mediumorchid')

ax1.set_yscale("log")
ax1.set_xticks([])
ax1.set_ylabel("run-time (s, log-scale)")
ax2.set_xlabel("number of nucleotides")
ax2.set_ylabel(r'$\rho$ value')

print(len(list(range(N_LWB, N_MAX_MCF, N_STEP))))
print(len(avg1))
print(np.array(err1).shape)
print("len1", len(list(range(N_LWB, N_MAX_MCF, N_STEP))))
print("len2",len(avg1))
ax1.errorbar(list(range(N_LWB, N_MAX_MCF, N_STEP)), 
             avg1, 
             yerr=np.array(err1).transpose(), 
             fmt='o-', 
             ecolor='black',
             barsabove=True,
             label='m-MIS', 
             capsize = 2.)
ax1.errorbar(list(range(N_LWB, N_UPB, N_STEP)), 
             avg2, 
             yerr=np.array(err2).transpose(), 
             fmt='o-', 
             ecolor='black',
             barsabove=True,
             label=fname2.split('_')[3],
             capsize=2.)

fig.legend()

plt.show()
