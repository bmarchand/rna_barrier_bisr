import itertools
import random
from bisr_dpw.dpw_algorithm import DPW,  valid_layout
from bisr_dpw.graph_classes import Digraph


def test_dag():

    N = 5

    H = Digraph()

    for u in range(1, 5, 1):
        H.add_edge(u, u+1)

    seq = None
    k = 0

    while not seq:
        seq = DPW(H, k)
        if k > 90:
            break
        if not seq:
            k += 1

    assert(k == 0)


def test_dpw_equal_1():

    print("hello")
    H = Digraph()

    for u in range(1, 5, 1):
        H.add_edge(u, u+1)

    H.add_edge(5, 1)

    print(H.in_ngbh)
    print(H.out_ngbh)

    seq = None
    k = 0

    seq = DPW(H, k)

    print("I am printing result ", seq)
    assert(seq is None)

    k = 1

    seq = DPW(H, k)

    print("working seq, ", seq)
    assert(seq is not None)


def test_dpw_equal_2():

    H = Digraph()

    for u in range(1, 5, 1):
        H.add_edge(u, u + 1)

    H.add_edge(1, 3)
    H.add_edge(1, 4)
    H.add_edge(1, 5)
    H.add_edge(2, 4)
    H.add_edge(2, 5)
    H.add_edge(3, 5)

    H.add_edge(4, 1)
    H.add_edge(5, 2)

    k = 0

    seq = DPW(H, k)
    assert(seq is None)

    k = 1
    seq = DPW(H, k)
    assert(seq is None)

    k = 2
    seq = DPW(H, k)
    print(seq)
    assert(seq is not None)


def test_small_graph0():

    H = Digraph()

    H.add_edge(1, 2)
    H.add_edge(1, 3)
    H.add_edge(2, 1)
    H.add_edge(2, 5)
    H.add_edge(3, 2)
    H.add_edge(3, 5)
    H.add_edge(4, 2)
    H.add_edge(4, 3)
    H.add_edge(5, 1)

    k = 0

    seq = DPW(H, k)
    assert(seq is None)

    k = 1

    seq = DPW(H, k)
    assert(seq is not None)

def test_small_graph1():

    H = Digraph()

    for u in range(5):
        H.add_node(u)

    H.add_edge(0, 1)
    H.add_edge(0, 2)
    H.add_edge(0, 3)
    H.add_edge(0, 5)
    
    H.add_edge(1, 3)
    H.add_edge(1, 5)

    H.add_edge(2, 0)
    H.add_edge(2, 1)
    H.add_edge(2, 3)
    H.add_edge(2, 5)

    H.add_edge(3, 2)
    H.add_edge(3, 5)

    H.add_edge(4, 5)

    k = 0
    seq = DPW(H, k, break_into_scc=True)
    assert(seq is None)

    k = 1

    seq = DPW(H, k, break_into_scc=True)
    assert(seq is not None)

    print("seq", seq)

def test_random_small_graphs():

    random.seed(2020)

    N = 5

    nruns = 20

    for _ in range(nruns):

        H = Digraph()

        for u in range(1, N+1, 1):
            for v in range(1, N+1, 1):
                if u != v:
                    if random.random() < 0.5:
                        H.add_edge(u, v)

        seq = None
        k = 0

        print(H.out_ngbh)

        while not seq:
            print("k", k)
            seq = DPW(H, k)
            if k > 90:
                break
            if not seq:
                k += 1

        assert(valid_layout(H, seq, k))

        best_k = 10**9
        best_sigma = seq
        for sigma in itertools.permutations(list(range(1, N+1))):

            kk = 0

            while not valid_layout(H, sigma, kk):
                kk += 1
                if kk > 90:
                    break

            if kk < best_k:
                best_k = kk
                best_sigma = sigma

        print("best found ", best_sigma, "best_k =", best_k)
        if k == best_k:
            print("not better than what DPW gets")
        else:
            print("better than DPW solution: ", seq[1:], "k =", k)
        assert(k == best_k)


def test_suppression():

    H = Digraph()

    H.add_edge(1, 2)
    H.add_edge(1, 3)
    H.add_edge(1, 4)
    H.add_edge(1, 5)
    H.add_edge(2, 3)
    H.add_edge(2, 4)
    H.add_edge(2, 5)

    H.add_edge(3, 4)
    H.add_edge(4, 5)
    H.add_edge(5, 3)

    k = 0

    seq = DPW(H, k)

    assert(seq is None)

    k = 1

    seq, trie = DPW(H, k, full_results=True, break_into_scc=False)

    assert(seq is not None)

    assert(len(trie.generate_sequences(depth_limit=5)) == 2)


if __name__ == '__main__':
    pass
#    test_dag()
#    test_dpw_equal_1()
    test_dpw_equal_2()
#    test_random_small_graphs()
#    test_small_graph0()
#    test_small_graph1()
#    test_suppresion()
