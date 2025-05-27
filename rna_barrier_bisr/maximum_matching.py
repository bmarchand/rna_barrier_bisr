from bisr_dpw_cpp_routines import hopcroft_karp
from .graph_classes import BipartiteGraph


def maximum_matching(G):

    new_G = BipartiteGraph()

    int_map = {}
    inv_map = {}

    for k, u in enumerate(G.left):

        new_G.add_node(k, 0)

        int_map[k] = u
        inv_map[u] = k

    for l in range(len(G.left), len(G.left)+len(G.right), 1):

        new_G.add_node(l, 1)

        int_map[l] = G.right[l-len(G.left)]
        inv_map[G.right[l-len(G.left)]] = l

    for k, u in enumerate(G.left):
        for v in G.ngbh[u]:
            new_G.add_edge(k, inv_map[v])

    matching = hopcroft_karp(new_G.ngbh, new_G.left, new_G.right)

    actual_matching = []

    for k, l in matching:
        actual_matching.append((int_map[k], int_map[l]))

    return actual_matching
