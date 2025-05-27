from barriers_py3 import randomDistantStructPairs, ssstring, realize, ssparse, sstopairs
import json
import time

THETA = 3
N_LWB = 4
N_UPB = 60
N_STEP = 10
N_SEEDS = 10
RECOMPUTE = True


with open('experiments/data/random_rna1_results_mcf_theta'+str(THETA)+'.json','r') as f:
    res_dict = json.load(f)


for n in range(N_LWB, N_UPB, N_STEP):
    for seed in range(1, N_SEEDS+1, 1):
        with open('experiments/data/random_rna1_results_mcf_theta'+str(THETA)+'.json','r') as f:
            res_dict = json.load(f)

        print("n, THETA, seed : ", n, THETA, seed)       
 
        if str(tuple((n, THETA, seed))) in res_dict.keys():
            if not RECOMPUTE:
                print("SKIPPED.")
                continue
        
        _, s1, s2 = randomDistantStructPairs(n, theta=THETA, seed=seed)

        bps1 = sstopairs(s1)
        bps2 = sstopairs(s2)

        # Remove common base pairs
        commonbps = bps1 & bps2
        if len(commonbps) > 0:
            for (i, j) in commonbps:
                (s1[i], s1[j]) = (-1, -1)
                (s2[i], s2[j]) = (-1, -1)
                bps1 = sstopairs(s1)
                bps2 = sstopairs(s2)
        
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
            P = realize(s1,s2,k)
            print("k:",k)
            if P:
                break
            k += 1


        rt = time.time() - t0

        res_dict[str(tuple((n, THETA, seed)))] = str(k)+'..'+str(rt)

        with open('experiments/data/random_rna1_results_mcf_theta'+str(THETA)+'.json', 'w') as f:
            json.dump(res_dict, f)
