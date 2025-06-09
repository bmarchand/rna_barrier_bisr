import matplotlib.pylab as plt
import numpy as np
import seaborn as sns
import json

with open('rho_data.json') as f:
    rho_dict = json.load(f)
with open('run_time_data.json') as f:
    run_time_dict = json.load(f)
with open('run_time_data_divide.json') as f:
    run_time_dict_divide = json.load(f)

fig, (ax1, ax2) = plt.subplots(2, figsize=(8,6))

MAX_LENGTH = 40

ax1.set_yscale("log")
ax1.set_xticks(list(range(10,60,5)))
ax1.set_ylabel('run-time(secs, log-scale)')
ax2.set_ylabel(r'$\Phi$ value')
#ax3.set_ylabel(r'$\rho$ value')
ax2.set_xlabel("number of vertices in conflict graph")

ns = set([])
rna_lengths = set([])
ns_mcf = set([])
rna_length_mcf = set([])

fname = 'random_structures_data_set_figure13.csv'
for line in open(fname).readlines()[1:]:
    rna_length = int(line.split(',')[1])
    num_vertices_bipartite_graph = int(line.split(',')[-1]) 
    if num_vertices_bipartite_graph==0:
        continue
    ns.add(num_vertices_bipartite_graph)
    rna_lengths.add(rna_length)
    if rna_length < MAX_LENGTH:
        ns_mcf.add(num_vertices_bipartite_graph)
        rna_length_mcf.add(rna_length)

ns = sorted(list(ns))
ns_mcf = sorted(list(ns_mcf))


rhos_per_n = {}
phis_per_n = {}
rt_per_n = {}
rt_per_n_divide = {}
for n in ns:
    rhos_per_n[n] = []
    phis_per_n[n] = []
    rt_per_n[n] = []
    rt_per_n_divide[n] = []

for line in open(fname).readlines()[1:]:
    phi = int(line.split(',')[-2])
    index = line.split(',')[0]
    n = int(line.split(',')[-1].rstrip('\n'))
    if n==0 or phi==0:
        continue
    rna_length = int(line.split(',')[1])
    rt = run_time_dict[index]
    phis_per_n[n].append(phi)
    rt_per_n[n].append(rt)

    if rna_length < MAX_LENGTH:
        
        rt_divide = run_time_dict_divide[index]
        rt_per_n_divide[n].append(rt_divide)
        rhos_per_n[n].append(rho_dict[index])

print(phis_per_n)
input()

ax1.errorbar(ns,
             [np.mean(rt_per_n[n]) for n in ns],
             yerr=[np.std(rt_per_n[n])/np.sqrt(len(rt_per_n[n])) for n in ns],
             fmt='o-',
             ecolor='black',
             barsabove=True,
             label='Subtree-Schedule',
             capsize=2.)

print([np.std(rt_per_n_divide[n]) for n in ns])
ax1.errorbar(ns,
             [np.mean(rt_per_n_divide[n]) for n in ns],
             yerr=[np.std(rt_per_n_divide[n])/np.sqrt(len(rt_per_n_divide[n])) for n in ns],
             fmt='o-',
             ecolor='black',
             barsabove=True,
             label='Divide-Schedule',
             capsize=2.)

ax2.errorbar(ns,
             [np.mean(phis_per_n[n]) for n in ns],
             yerr=[np.std(phis_per_n[n])/np.sqrt(len(phis_per_n[n])) for n in ns],
             fmt='o-',
             ecolor='black',
             barsabove=True,
             capsize=2.)

import pandas as pd
#df = pd.DataFrame(data = {n:phis_per_n[n] for n in ns})
#sns.boxplot(data=df, ax=ax2, color='mediumorchid')

#df = pd.DataFrame(data = {n:rhos_per_n[n] for n in ns if n < MAX_LENGTH})
#sns.boxplot(data=df, ax=ax3, color='mediumorchid')

fig.legend()
fig.savefig('compare_divide_subtree_figure13.pdf')
plt.show()
