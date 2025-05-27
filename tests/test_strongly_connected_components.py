from networkx.algorithms.components import number_strongly_connected_components
import networkx as nx
import random
from bisr_dpw.dpw_algorithm import Digraph
from bisr_dpw_cpp_routines import strongly_connected_components

def test_2components():

    out_ngbh = {}

    out_ngbh[0] = [1]
    out_ngbh[1] = [2]
    out_ngbh[2] = [0]
    
    out_ngbh[3] = [4]
    out_ngbh[4] = [5]
    out_ngbh[5] = [3]

    out_ngbh[1].append(3)

    in_ngbh = {}

    for key in out_ngbh.keys():
        in_ngbh[key] = []

    for key in out_ngbh.keys():
        for neigh in out_ngbh[key]:
            in_ngbh[neigh].append(key)    

    sccs = strongly_connected_components(out_ngbh,
                                         in_ngbh)

    assert(len(sccs)==2)

def test_dag():

    H = Digraph()

    for u in range(1, 5, 1):
        H.add_edge(u,u+1)

    sccs = strongly_connected_components(H.out_ngbh, H.in_ngbh)

    assert(len(sccs)==5)   
    assert(sccs==[[1],[2],[3],[4],[5]])

def test_random_graphs():

    nruns = 20

    N = 100

    random.seed(2021)
 
    for _ in range(nruns):

        H = Digraph()
        nx_H = nx.DiGraph()

        for u in range(1, N+1, 1):
            for v in range(1, N+1, 1):
                if u != v:
                    if random.random() < 0.5:
                        H.add_edge(u, v)
                        nx_H.add_edge(u, v)

        sccs = strongly_connected_components(H.out_ngbh, H.in_ngbh) 

        nx_nsccs = number_strongly_connected_components(nx_H)

        assert(len(sccs)==nx_nsccs) 
