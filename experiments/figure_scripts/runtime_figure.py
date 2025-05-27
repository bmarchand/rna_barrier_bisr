import sys
import matplotlib.pylab as plt
import numpy as np
import seaborn as sns
import json
from general_bipartite import random_bipartite, vertex_cover_from_matching
from networkx.algorithms.bipartite.matching import maximum_matching

fname = sys.argv[-1]

with open(fname, 'r') as f:
    res_dict = json.load(f)

results_per_k = {}

ks = set([])

for key, value in res_dict.items():

    n = int(key[1:].split(',')[0])
    degree = int(key.split(',')[1])
    seed = int(key[:-1].split(',')[2])

    B, R, G = random_bipartite(n, degree, seed=seed)

    M = maximum_matching(G, top_nodes=B)
    K = vertex_cover_from_matching(B,R,G,M)
    I = set(G) - K

    k = int(value.split('..')[0]) + len(I) - len(B)

    ks.add(k)

    rt = float(value.split('..')[1])

    if k not in results_per_k.keys():

        results_per_k[k] = {}
        results_per_k[k][n] = [rt]

    else:
        if n not in results_per_k[k].keys():
            results_per_k[k][n] = [rt]
        else:  
            results_per_k[k][n].append(rt)

fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8, 4.8))

ns = set([])
ks_per_n = {}

for k in sorted(results_per_k.keys(), key=lambda x: int(x)):
    for n in sorted(results_per_k[k].keys()):
        ns.add(n)
        if n not in ks_per_n.keys():
            ks_per_n[n] = [k]
        else:
            ks_per_n[n].append(k)

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
df = pd.DataFrame(data = {2*n:ks_per_n[n] for n in sorted(ks_per_n.keys())})

sns.boxplot(data=df, ax=ax3, color='mediumorchid')

import matplotlib.cm as cm


N_c = len(results_per_k.keys())

for c, k in enumerate(sorted(results_per_k.keys(), key=lambda x: int(x))):

    xs = []
    ys = []
    err_ys = []

    for n in sorted(results_per_k[k].keys()):
        xs.append(2*n)
        ys.append(np.mean(results_per_k[k][n]))
        err_ys.append(np.std(results_per_k[k][n])/float(np.sqrt(len(results_per_k[k][n]))))    

    ax1.errorbar(xs,ys, yerr=err_ys, fmt='o-', color=cm.viridis(float(c)/float(N_c-1)))
    ax2.errorbar(xs,ys, label=r'$\rho$='+str(k), fmt='o-', color=cm.viridis(float(c)/float(N_c-1)))



fig.legend()
ax1.set_ylabel("runtime (s)")
ax1.set_xticks([])
ax2.set_xticks([])
ax2.set_yscale("log")
ax2.set_ylabel("runtime (s)\n (log scale)")
ax3.set_ylabel(r'average $\rho$ value')

plt.show()
