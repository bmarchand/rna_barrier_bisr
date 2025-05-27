from .barriers_py3 import ssrandom, ssparse, sstopairs, ssstring, realize
from .graph_classes import Digraph, BipartiteGraph

import random

def solve_mcf(s1, s2, theta=1):

    bps1, bps2, ss1, ss2, inversion = pre_processing(s1, s2)

    k = 0

    print("ss1,ss2", ss1, ss2)

    while k < len(s1):

        p = realize(ss1, ss2, k, 1, theta=theta)
        if p is not None:
            break
        k += 1

    if inversion:
        return k-(len(bps1)-len(bps2)), p[::-1]
    else:
        return k, p


def pre_processing(s1, s2):

    inversion = False
    if s2.count("(") > s1.count("("):
        s1, s2 = s2, s1
        inversion = True

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

    return bps1, bps2, ss1, ss2, inversion


def structures_to_conflict_graph(s1, s2):

    bps1, bps2, ss1, ss2, inversion = pre_processing(s1, s2)

    G = BipartiteGraph()

    if inversion:
        side1 = 1
        side2 = 0
    else:
        side1 = 0
        side2 = 1

    for u in bps1:
        G.add_node(u, side1)
        for v in bps2:
            G.add_node(v, side2)

            conflict = False

            if len(set([u[0], v[0], u[1], v[1]])) != 4:
                conflict = True

            if (u[0] < v[0] < u[1] < v[1]) or (v[0] < u[0] < v[1] < u[1]):
                conflict = True

            if conflict:
                G.add_edge(u, v)

    return G


def random_conflict_graph(n, theta, unpaired_weight, seed=None):

    if seed:
        random.seed(seed)

    count = {}

    s1 = ssrandom(n, count, theta=theta, unpaired_weight=unpaired_weight)
    s2 = ssrandom(n, count, theta=theta, unpaired_weight=unpaired_weight)

    return structures_to_conflict_graph(s1, s2)


def layout_to_pathway(H, seq):
    """
    Conversion from a "Tamaki-style" layout (seq) to a reconfiguration
    pathway of equivalent quality for the corresponding underlying
    conflict graph G.

    An "interval representation" of H will be used as an intermediary
    structure.

    Args:
        - H: directed graph, whose vertices are the edges of a perfect matching
    of the conflict graph G. If G did not have a perfect matching, then
    artificial nodes have been added to it.

        - seq: an ordering of the vertices of H.

    Warning:
        - the code relies on the vertices of H being 2-uples of objects.
    (as they are supposed to be edges of a perfect matching)
    """

    opened = {}  # dictionary for keeping track of which intervals are open.
    closed = {}

    when_opened = {}  # dictionary for storing when interval starts
    when_closed = {}  # "       "     "     "      "     ends

    prefix_in_ngbh = set([])  # in neighbors of current prefix in loop

    for u in seq:
        opened[u] = False
        closed[u] = False

    for k, u in enumerate(seq):

        print("adding u", u)

        # if not opened, open it just before closing it
        if not opened[u]:
            when_opened[u] = k

        # when we find a vertex, we close its interval
        when_closed[u] = k
        opened[u] = False
        closed[u] = True

        # update in_ngbh
        prefix_in_ngbh.discard(u)
        for v in H.in_ngbh[u]:
            if not closed[v]:
                prefix_in_ngbh.add(v)

        # opening newly added vertices
        for v in prefix_in_ngbh:
            if not opened[v]:
                opened[v] = True
                when_opened[v] = k

        print("n(sigma)- ", prefix_in_ngbh)

    print("when_opened", when_opened)
    print("when_closed", when_closed)

    # converting to reconfiguration pathway

    pathway = []  # At first, will contain at position k a list of
    #  vertices processed at step k. Then will be unrolled.

    for _ in range(len(seq)):
        pathway.append([])

    # opening events
    for u in seq:
        pathway[when_opened[u]].append(u[0])  # relying on u being 2-uple

    # closing events
    for u in seq:
        pathway[when_closed[u]].append(u[1])

    print("pathway brut", pathway)

    result = []

    # unrolling
    for sub_list in pathway:
        result += sub_list

    return result
