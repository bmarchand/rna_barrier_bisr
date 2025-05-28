import random
import networkx as nx
import json
import time
from rna_barrier_bisr.general_bipartite import realize, random_bipartite

DEGREE = 5
N_LWB = 4
N_UPB = 41
N_STEP = 6
N_SEEDS = 200
RECOMPUTE = False

fname = 'general_bipartite_strict_divider_results_'+str(DEGREE)+'.json'

try:
    with open(fname,'r') as f:
        res_dict = json.load(f)
except:
    res_dict = {}

for n in range(N_LWB, N_UPB, N_STEP):
    for seed in range(1, N_SEEDS+1, 1):

        print("n, DEGREE, seed : ", n, DEGREE, seed)       
 
        if str(tuple((n, DEGREE, seed))) in res_dict.keys():
            if not RECOMPUTE:
                print("SKIPPED.")
                continue

        B, R, G = random_bipartite(int(n/2), DEGREE, seed=seed)

        t0 = time.time()

        k = 0

        while True:
            P = realize(B,R,G,k)
            if P:
                break
            k += 1

        rt = time.time() - t0

        res_dict[str(tuple((n, DEGREE, seed)))] = str(k)+'..'+str(rt)

        with open(fname, 'w') as f:
            json.dump(res_dict, f)
