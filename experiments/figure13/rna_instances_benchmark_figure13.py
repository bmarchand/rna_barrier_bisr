import matplotlib.pylab as plt
import numpy as np
import seaborn as sns
import json

with open('rho_data.json') as f:
    rho_dict = json.load(f)
with open('run_time_data.json') as f:
    run_time_dict = json.load(f)
with open('run_time_data_dpw.json') as f:
    run_time_dict_dpw = json.load(f)

fig, (ax1, ax2, ax3) = plt.subplots(3, figsize=(8,6))

ax1.set_yscale("log")
ax1.set_xticks(list(range(10,60,5)))
ax1.set_ylabel('run-time(secs, log-scale)')
ax2.set_ylabel(r'$\Phi$ value')
ax3.set_ylabel(r'$\rho$ value')
ax3.set_xlabel("number of nucleotides")

ns = set([])

fname = 'random_structures_data_set.csv'
for line in open(fname).readlines()[1:]:
    ns.add(int(line.split(',')[1]))

ns = sorted(list(ns))

rhos_per_n = {}
phis_per_n = {}
rt_per_n = {}
rt_per_n_dpw = {}
for n in ns:
    rhos_per_n[n] = []
    phis_per_n[n] = []
    rt_per_n[n] = []
    rt_per_n_dpw[n] = []

for line in open(fname).readlines()[1:]:
    phi = int(line.split(',')[-1].rstrip('\n'))
    index = line.split(',')[0]
    n = int(line.split(',')[1])
    rt = run_time_dict[index]
    phis_per_n[n].append(phi)
    rhos_per_n[n].append(rho_dict[index])
    rt_per_n[n].append(rt)

    rt_dpw = run_time_dict_dpw[index]
    rt_per_n_dpw[n].append(rt_dpw)

ax1.errorbar(ns,
             [np.mean(rt_per_n[n]) for n in ns],
             yerr=[np.std(rt_per_n[n])/np.sqrt(len(rt_per_n[n])) for n in ns],
             fmt='o-',
             ecolor='black',
             barsabove=True,
             label='arboricity',
             capsize=2.)

print([np.std(rt_per_n_dpw[n]) for n in ns])
ax1.errorbar(ns,
             [np.mean(rt_per_n_dpw[n]) for n in ns],
             yerr=[np.std(rt_per_n_dpw[n])/np.sqrt(len(rt_per_n_dpw[n])) for n in ns],
             fmt='o-',
             ecolor='black',
             barsabove=True,
             label='dpw',
             capsize=2.)

import pandas as pd
df = pd.DataFrame(data = {n:phis_per_n[n] for n in ns})
sns.boxplot(data=df, ax=ax2, color='mediumorchid')

df = pd.DataFrame(data = {n:rhos_per_n[n] for n in ns})
sns.boxplot(data=df, ax=ax3, color='mediumorchid')

fig.legend()
fig.savefig('left_panel_bench.pdf')
plt.show()
