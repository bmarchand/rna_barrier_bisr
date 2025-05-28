from barriers_py3 import ssrandom, ssstring, realize, ssparse, sstopairs
import json
import time

THETA = 5
N_LWB = 10
N_UPB = 90
N_STEP = 10
N_SEEDS = 500
RECOMPUTE = True


for n in range(N_LWB, N_UPB, N_STEP):
    for seed in range(50, N_SEEDS+1, 1):

        with open('experiments/data/random_rna2_results_mcf_theta'+str(THETA)+'.json','r') as f:
            res_dict = json.load(f)

        print("n, THETA, seed : ", n, THETA, seed)       
 
        if str(tuple((n, THETA, seed))) in res_dict.keys():
            if not RECOMPUTE:
                print("SKIPPED.")
                continue
        
        count = {}

        s1 = ssrandom(n, count, seed=2*seed)
        s2 = ssrandom(n, count, seed=2*seed+1)

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

        while True:
            P = realize(ss1,ss2,k)
            print("k:",k)
            if P:
                break
            k += 1


        rt = time.time() - t0

        res_dict[str(tuple((n, THETA, seed)))] = str(k)+'..'+str(rt)

        with open('experiments/data/random_rna2_results_mcf_theta'+str(THETA)+'.json', 'w') as f:
            json.dump(res_dict, f)
