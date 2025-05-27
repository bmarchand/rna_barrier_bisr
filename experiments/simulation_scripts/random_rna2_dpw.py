from barriers_py3 import randomDistantStructPairs, ssstring, realize, ssparse, sstopairs, ssrandom
import json
from networkx.algorithms.bipartite import maximum_matching
from dpw_interface import secondary_structures_to_G, extend_to_pm, construct_digraph
from dpw_algorithm import DPW 
import time

THETA = 5
N_LWB = 150
N_UPB = 160
N_STEP = 10
N_SEEDS = 30
RECOMPUTE = False


with open('experiments/data/random_rna2_results_dpw_theta'+str(THETA)+'.json','r') as f:
    res_dict = json.load(f)


for n in range(N_LWB, N_UPB, N_STEP):
    for seed in range(1, N_SEEDS, 1):

        print("n, THETA, seed : ", n, THETA, seed)       
        with open('experiments/data/random_rna2_results_dpw_theta'+str(THETA)+'.json','r') as f:
            res_dict = json.load(f)
 
        if str(tuple((n, THETA, seed))) in res_dict.keys():
            if not RECOMPUTE:
                print("SKIPPED.")
                continue
        
        if (n, seed) == (120,53):
            continue

        
        count = {}

        s1 = ssrandom(n, count, seed=2*seed, THETA=THETA)
        s2 = ssrandom(n, count, seed=2*seed+1, THETA=THETA)

        if s1 is None or s2 is None:    
            print("None, continue")
            continue

        if s2.count("(") > s1.count("("):
            s1, s2 = s2, s1

        ss1 = ssparse(s1)
        ss2 = ssparse(s2)
        bps1 = sstopairs(ss1)
        bps2 = sstopairs(ss2)

        # Remove common base pairs
        commonbps = bps1 & bps2
        if len(commonbps) > 0:
            for (i, j) in commonbps:
                (ss1[i], ss1[j]) = (-1, -1)
                (ss2[i], ss2[j]) = (-1, -1)
                s1 = ssstring(ss1)
                s2 = ssstring(ss2)
                bps1 = sstopairs(ss1)
                bps2 = sstopairs(ss2)
        
        t0 = time.time()

        k = 0
    
        if len(bps1)==0:
            rt = time.time() - t0
            res_dict[str(tuple((n, THETA, seed)))] = str(k)+'..'+str(rt)
            continue
        if len(bps2)==0:
            rt = time.time() - t0
            k = len(bps1) 
            res_dict[str(tuple((n, THETA, seed)))] = str(k)+'..'+str(rt)
            continue

        B, R, G = secondary_structures_to_G(ss1, ss2) 
        M = maximum_matching(G, top_nodes = B)

        n_ur = 0
        for u in R:
            if u not in M.keys():
                n_ur += 1 

        B, R, G, M = extend_to_pm(B,R,G, M)

        H, int_to_e, e_to_int = construct_digraph(B,R,G,M)

        dpw = 0

        while True:
            P = DPW(H, dpw, break_into_scc=True) 
            if P:
                break
            dpw += 1

        rt = time.time() - t0

        res_dict[str(tuple((n, THETA, seed)))] = str(dpw-n_ur+1)+'..'+str(rt)+'..'+str(dpw)

        with open('experiments/data/random_rna2_results_dpw_theta'+str(THETA)+'.json', 'w') as f:
            json.dump(res_dict, f)
