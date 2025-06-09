from rna_barrier_bisr.rna_interface import solve_mcf
from rna_barrier_bisr.dpw_interface import secondary_structures_to_G
import time
from networkx.algorithms.bipartite.matching import maximum_matching
from rna_barrier_bisr.general_bipartite import random_bipartite, vertex_cover_from_matching
from rna_barrier_bisr.random_rna_structures import sstopairs, ssparse
from networkx.algorithms.bipartite.matching import maximum_matching

fname = 'random_structures_data_set_figure13.csv'

run_time_dict = {}
rho_dict = {}

N = len(open(fname).readlines())

MAX_LENGTH = 40

for line in open(fname).readlines()[1:]:
    index = line.split(',')[0]
    print(index, 'out of',N)
    length = int(line.split(',')[1])
    if length >= MAX_LENGTH:
        continue
    s1 = line.split(',')[2] 
    s2 = line.split(',')[3] 
    print(s1,'->',s2)

    sizeL = s1.count("(")
    B,R,G = secondary_structures_to_G(ssparse(s1),ssparse(s2))
    M = maximum_matching(G, top_nodes=B)
    K = vertex_cover_from_matching(B,R,G,M)
    I = set(G) - K
    alpha=len(I)


    t0 = time.time()
    k, _ = solve_mcf(s1, s2)
    run_time = time.time() - t0
    run_time_dict[index] = run_time


    rho_dict[index] = k+alpha-sizeL

import json
with open('run_time_data_divide.json','w') as f:
    json.dump(run_time_dict,f)

with open('rho_data.json','w') as f:
    json.dump(rho_dict,f)
