from .graph_classes import Digraph
from .dpw_algorithm import DPW
from .maximum_matching import maximum_matching
from .random_rna_structures import ssparse, sstopairs
from .rna_interface import structures_to_conflict_graph
import networkx as nx

def solve_dpw(s1, s2):
    G = structures_to_conflict_graph(s1, s2)

    k, n_ur, schedule = solve_graph_dpw(G)

    return k-n_ur+1, schedule, k

def solve_graph_dpw(G):

    if len(G.left) == 0:
        return -1, 0, list(G.right)

    M = maximum_matching(G)

    M_dict = {}

    for u, v in M:
        M_dict[u] = v
        M_dict[v] = u

    n_ur = 0
    for u in G.right:
        if u not in M_dict.keys():
            n_ur += 1

    G_ext, M_full = extend_to_pm(G, M)

    M_dict = {}

    for u, v in M_full:
        M_dict[u] = v
        M_dict[v] = u

    H, int_to_e, e_to_int = construct_digraph(G_ext.left, G_ext.right, G_ext, M_dict)

    k = 0

    while k <= G_ext.n_nodes:

        seq = DPW(H, k, break_into_scc=True)

        if seq is not None:
            break

        k += 1

    return k, n_ur, [int_to_e[i] for i in seq]


def extend_to_pm(G, M):

    matched = {}

    for u, v in M:
        matched[u] = True
        matched[v] = True

    G_ext = G.copy()
    M_full = [e for e in M] 

    for u in G_ext.side.keys():
        try:
            matched[u]
        except KeyError:
            matched[u] = False

    for u in matched.keys():
        if not matched[u]:
            G_ext.add_node(G_ext.n_nodes, 1 - G_ext.side[u])
            for v in G_ext.ngbh.keys():
                if G_ext.side[v] == G_ext.side[u]:
                    G_ext.add_edge(v, G_ext.n_nodes - 1)

            if G_ext.side[u] == 0:
                M_full.append((u, G_ext.n_nodes - 1))
            if G_ext.side[u] == 1:
                M_full.append((G_ext.n_nodes - 1, u))

    return G_ext, M_full


def secondary_structures_to_G(s1, s2):

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

    G = nx.Graph()
    B = set([])
    R = set([])

    for u in bps1:
        G.add_node(u)
        B.add(u)
        for v in bps2:
            G.add_node(v)
            R.add(v)

            conflict = False

            if len(set([u[0], v[0], u[1], v[1]])) != 4:
                conflict = True

            if (u[0] < v[0] < u[1] < v[1]) or (v[0] < u[0] < v[1] < u[1]):
                conflict = True

            if conflict:
                G.add_edge(u, v)

    return B, R, G


def construct_digraph(B, R, G, M):

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
                    H.add_edge(e_to_int[(u, M[u])], e_to_int[(v, M[v])])

    return H, int_to_e, e_to_int


def extend_to_pm_nx(B, R, G, M):

    N = G.order()

    for b in B:
        if b not in M.keys():
            G.add_node(N)
            R.add(N)
            for b2 in B:
                G.add_edge(b2, N)
            M[N] = b
            M[b] = N

            N += 1

    for r in R:
        if r not in M.keys():
            G.add_node(N)
            B.add(N)
            for r2 in R:
                G.add_edge(r2, N)
            M[N] = r
            M[r] = N

            N += 1

    return B, R, G, M
