from rna_barrier_subtree_solver.subtree_solver import solve_subtree
import time

fname = 'random_structures_data_set.csv'

run_time_dict = {}

N = len(open(fname).readlines())

for line in open(fname).readlines()[1:]:
    index = line.split(',')[0]
    print(index, 'out of',N)
    s1 = line.split(',')[2] 
    s2 = line.split(',')[3] 
    print(s1,'->',s2)
    t0 = time.time()
    _, _, _ = solve_subtree(s1, s2)
    run_time = time.time() - t0
    run_time_dict[index] = run_time

import json
with open('run_time_data.json','w') as f:
    json.dump(run_time_dict,f)
