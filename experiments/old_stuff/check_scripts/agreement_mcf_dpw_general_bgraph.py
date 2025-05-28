import random
import networkx as nx
from general_bipartite import realize, random_bipartite
from dpw_interface import construct_digraph
from dpw_interface import extend_to_pm
from dpw_algorithm import DPW
from graph_classes import Digraph
from maximum_matching import maximum_matching
import json
import time

DEGREE = 3
N_LWB = 1
N_UPB = 20
N_STEP = 1
N_SEEDS = 300


for n in range(N_LWB, N_UPB, N_STEP):
    for seed in range(1, N_SEEDS+1, 1):

        print("n, DEGREE, seed : ", n, DEGREE, seed)       
        
        B, R, G = random_bipartite(n, DEGREE, seed=seed)
        
        mcf_k = 0

        while True:
            P = realize(B,R,G,mcf_k)
            if P:
                break
            mcf_k += 1

        t0 = time.time()

        M = maximum_matching(G)

        n_ur = 0
        for u in R:
            if u not in M.keys():
                n_ur += 1 

        G, M = extend_to_pm(G, M)

        H, int_to_e, e_to_int = construct_digraph(B,R,G,M)
        
        dpw = 0

        print("n_ur, ", n_ur)

        while True:
            dpw_P = DPW(H, dpw, break_into_scc=True) 
            if dpw_P:
                break
            dpw += 1

        print("dpw: ", dpw)

        dpw_k = dpw - n_ur + 1

        try:
            assert(dpw_k==mcf_k)
        except AssertionError:
            print("MCF P: ", P)
            print("dpw_k: ", dpw_k, "mcf_k: ", mcf_k)
            B, R, G = random_bipartite(n, DEGREE, seed=seed)

            for u in G.nodes:
                print(u, ": ", [v for v in G.neighbors(u)])

            M = maximum_matching(G, top_nodes = B)
            n_ur = 0
            for u in R:
                if u not in M.keys():
                    n_ur += 1 
            print("maximum matching ", M)
            print("n_ur", n_ur)
            B, R, G, M = extend_to_pm(B,R,G, M)
            print("extnded matching ", M)
            print("new B: ", B, "new R: ", R)
            H = Digraph()

            e_to_int = {}
            int_to_e = {}

            for l, u in enumerate(B):
                e_to_int[(u, M[u])] = l
                int_to_e[l] = (u, M[u])
                H.add_node(l)

            
            for u in B:
                for v in B:
                    if u != v:
                        if G.has_edge(u, M[v]):
                            H.add_edge(e_to_int[(u,M[u])], e_to_int[(v,M[v])])

            for u in H.out_ngbh.keys():
                print(u, ": ", [v for v in H.out_ngbh[u]])

            print("mapping: ", int_to_e)

            raise AssertionError
