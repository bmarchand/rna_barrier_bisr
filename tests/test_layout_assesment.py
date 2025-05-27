from bisr_dpw.dpw_algorithm import valid_layout, Digraph


def test_valid_layout():

    def test_graph1():

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

        seq = [1, 2, 3, 4, 5]

        assert(valid_layout(H, seq, 2))

    def test_graph2():

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

        print("test 2 - call 1")
        assert(valid_layout(H, [4, 3, 2, 1, 5], 1))
        print("test 2 - call 2")
        assert(valid_layout(H, [4, 3, 2, 1, 5], 2))

    test_graph1()
    test_graph2()
