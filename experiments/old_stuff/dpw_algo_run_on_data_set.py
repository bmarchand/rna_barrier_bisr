from bisr_dpw.dpw_interface import solve_dpw
import time

fname = 'random_structures_data_set.csv'

run_time_dict = {}
rho_dict = {}

N = len(open(fname).readlines())

for line in open(fname).readlines()[1:]:
    index = line.split(',')[0]
    print(index, 'out of',N)
    s1 = line.split(',')[2] 
    s2 = line.split(',')[3] 
    print(s1,'->',s2)
    t0 = time.time()
    rho, _, _ = solve_dpw(s1, s2)
    run_time = time.time() - t0
    run_time_dict[index] = run_time
    rho_dict[index] = rho

import json
with open('run_time_data_dpw.json','w') as f:
    json.dump(run_time_dict,f)

with open('rho_data.json','w') as f:
    json.dump(rho_dict,f)
