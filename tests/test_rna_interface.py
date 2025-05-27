from bisr_dpw.barriers_py3 import ssrandom
from bisr_dpw.rna_interface import random_conflict_graph, layout_to_pathway
from bisr_dpw.dpw_interface import extend_to_pm
from bisr_dpw.graph_classes import Digraph
from bisr_dpw.maximum_matching import maximum_matching
from bisr_dpw.dpw_algorithm import DPW

def test_random_conflict_graph():

    n = 100

    random_conflict_graph(n, 3, 0.1, seed=124)
    
def test_small_instance():

    n = 10

    G = random_conflict_graph(n, seed=4, theta=1, unpaired_weight=0.2)

    print(G.ngbh)

    assert(len(G.ngbh.keys())==8)
    assert(sorted(G.ngbh[(1,3)]) == [(1,8), (2,7)])

def test_layout_to_pathway():

    ns = [15, 20, 25]

    nruns = 20

    seed = 2021

    for n in ns:
        for inc in range(nruns):
            k, n_ur, int_to_e, H, seq = dpw_solution(n, seed+inc)
    
            G = random_conflict_graph(n, 
                                      seed=seed+inc, 
                                      theta=1, 
                                      unpaired_weight=0.2)
            
            M = maximum_matching(G)

            G, M = extend_to_pm(G, M)

            check_solution(k, n_ur, int_to_e, H, seq, G)

def check_solution(k, n_ur, int_to_e, H, seq, G):

    new_H = Digraph()

    for u in H.in_ngbh.keys():
        new_H.add_node(int_to_e[u])

    for u in H.out_ngbh.keys():
        for v in H.out_ngbh[u]:
            new_H.add_edge(int_to_e[u],int_to_e[v])

    if len(seq) > 0:
        if seq[0] == -1:
            seq = seq[1:]

    new_seq = []

    for u in seq:
        new_seq.append(int_to_e[u])

    pathway = layout_to_pathway(new_H, new_seq)

    print("new_seq ", new_seq)
    print("G.ngbh ", G.ngbh)
    print("G.left", G.left)
    print("G.right", G.right)
    print("pathway ", pathway)

    delta = 0 # must remain \leq k+1-n_ur

    active_set = set(G.left) 

    for u in pathway:
        print("processing ", u)
        if G.side[u] == 0:
            active_set.remove(u)
            delta += 1
            assert(delta <= k+1)
        elif G.side[u] == 1:
            active_set.add(u)
            delta -= 1
            for v in G.ngbh[u]:
                assert(v not in active_set)

        print("delta", delta)

    assert(active_set==set(G.right))        

    

def dpw_solution(n, seed):

    G = random_conflict_graph(n, seed=seed, theta=1, unpaired_weight=0.2)

    M = maximum_matching(G)
    
    matched = {}

    for u, v in M:
        matched[u] = True
        matched[v] = True

    for u in G.side.keys():
        try:
            matched[u]
        except KeyError:
            matched[u] = False

    n_ur = 0
    for u in matched.keys():
        if G.side[u] == 1 and matched[u] == False:
            n_ur += 1 

    G, M = extend_to_pm(G, M)

    H = Digraph()

    e_to_int = {}
    int_to_e = {}

    for k, (u, v) in enumerate(M):
        e_to_int[(u,v)] = k
        int_to_e[k] = (u,v)
        H.add_node(k)

    
    for u,v in M:
        for w, x in M:
            if x in G.ngbh[u]:
                H.add_edge(e_to_int[(u,v)], e_to_int[(w,x)])

    k = 0

    while k <= n:
    
        seq = DPW(H, k, break_into_scc=True)      
 
        if seq is not None:
            break

        k+= 1
    
    return k, n_ur, int_to_e, H, seq 

if __name__=='__main__':
    test_small_instance()
