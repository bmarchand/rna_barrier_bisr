import random
from dpw_algorithm import DPW
from graph_classes import Digraph
import networkx as nx
from general_bipartite import random_bipartite
from dpw_interface import extend_to_pm, construct_digraph
from networkx.algorithms.bipartite import sets, maximum_matching, to_vertex_cover
import json
import time

DEGREE = 5
N_LWB = 4
N_UPB = 60
N_STEP = 6
N_SEEDS = 200
RECOMPUTE = False


for n in range(N_LWB, N_UPB, N_STEP):
    for seed in range(1, N_SEEDS+1, 1):
        with open('experiments/data/general_bipartite_results_dpw_degree'+str(DEGREE)+'.json','r') as f:
            res_dict = json.load(f)

        print("n, DEGREE, seed : ", n, DEGREE, seed)       
 
        if str(tuple((n, DEGREE, seed))) in res_dict.keys():
            if not RECOMPUTE:
                print("SKIPPED.")
                continue


        B, R, G = random_bipartite(int(n/2), DEGREE, seed=seed)

        t0 = time.time()

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

        res_dict[str(tuple((n, DEGREE, seed)))] = str(dpw-n_ur+1)+'..'+str(rt)

        with open('experiments/data/general_bipartite_results_dpw_degree'+str(DEGREE)+'.json', 'w') as f:
            json.dump(res_dict, f)
