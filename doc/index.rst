.. bisr-dpw documentation master file, created by
   sphinx-quickstart on Wed Apr 28 17:24:54 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the documentation of bisr-dpw!
=========================================

Source code documentation for the implementation of the algorithms used
in `this paper <https://hal.inria.fr/hal-03272963>`_.

The following code example shows how to import and use the **directed pathwidth**
algorithm:

.. code-block:: python

    
    from graph_classes import Digraph
    from dpw_algorithm import DPW

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

    seq = None
    k = 0

    while not seq:
        seq = DPW(H, k)
        if not seq:
            k += 1

    print("directed pathwidth value: ", k)
    print("layout: ", seq)

The output should be:

.. code-block:: python

    >> directed pathwidth value: 2
    >> layout: [1, 3, 4, 2, 5]    

The following code example shows how to import and use our algorithm
for **bipartite independent set reconfiguration**:

.. code-block:: python

    import networkx as nx
    from general_bipartite import realize

    G = nx.Graph()
    
    for u in range(7):
        G.add_node(u)

    G.add_edge(0,4)
    G.add_edge(1,3)
    G.add_edge(1,4)
    G.add_edge(1,5)
    G.add_edge(2,4)
    G.add_edge(2,5)
    G.add_edge(2,6)
    G.add_edge(1,6)
    
    B, R = sets(G)
    
    k=0
    while True:
        #try with k
        P = realize(B,R,G,k)
        if P:
            break
        k += 1

    print("found k=", k)
    print("found P=", P)

.. autofunction:: dpw_algorithm.DPW

.. autofunction:: general_bipartite.realize

.. autoclass:: graph_classes.Digraph
    :special-members: __init__
    :members:
