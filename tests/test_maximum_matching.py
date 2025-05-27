from bisr_dpw.maximum_matching import maximum_matching
from bisr_dpw.graph_classes import BipartiteGraph
from bisr_dpw.dpw_interface import extend_to_pm

def test_complete_bipartite():

    N = 5

    G = BipartiteGraph()

    for i in range(N):
        G.add_node(i, 0)

    for i in range(N, N+N, 1):
        G.add_node(i, 1)

    for i in range(N):
        for j in range(N, N+N, 1):
            G.add_edge(i,j)

    M = maximum_matching(G)

    assert(len(M)==5)

def test_asymmetric_complete_bipartite():

    N = 5
    M = 4
    
    G = BipartiteGraph()

    for i in range(N):
        G.add_node(i, 0)

    for i in range(N, N+M, 1):
        G.add_node(i, 1)

    for i in range(N):
        for j in range(N, N+M, 1):
            G.add_edge(i,j)

    Matching = maximum_matching(G)

    assert(len(Matching)==4)

def test_slightly_more_complex_instance():

    G = BipartiteGraph()

    G.add_node(0, 0)
    G.add_node(1, 0)
    G.add_node(2, 0)

    G.add_node(3, 1)
    G.add_node(4, 1)
    G.add_node(5, 1)

    G.add_edge(0, 3)
    G.add_edge(1, 3)
    G.add_edge(0, 5)
    G.add_edge(2, 4)
    G.add_edge(2, 5)

    M = maximum_matching(G)

    assert(len(M)==3)
    assert(M == [(0, 5), (1, 3), (2, 4)])
    
    G, M = extend_to_pm(G, M)

    assert(len(M)==3)

def test_extend_matching():

    n = 5
    m = 4
    
    G = BipartiteGraph()

    for i in range(n):
        G.add_node(i, 0)

    for i in range(n, n+m, 1):
        G.add_node(i, 1)

    for i in range(n):
        for j in range(n, n+m, 1):
            G.add_edge(i,j)

    M = maximum_matching(G)

    assert(len(M)==4)
    
    G, M = extend_to_pm(G, M)

    assert(len(M) == 5)
