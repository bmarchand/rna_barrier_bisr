import copy

class BipartiteGraph:

    def __init__(self):

        self.n_nodes = 0
        self.ngbh = {}

        self.left = []
        self.right = []
        self.side = {}

    def add_node(self, u, side):

        if side not in [0, 1]:
            raise "specify a side with 0 or 1. 0 = left, 1 = right"

        try:
            self.ngbh[u]
        except KeyError:
            self.n_nodes += 1
            self.ngbh[u] = []
            self.side[u] = side
            if side == 0:
                self.left.append(u)
            if side == 1:
                self.right.append(u)

    def add_edge(self, u, v):

        self.ngbh[v].append(u)
        self.ngbh[u].append(v)

    def has_edge(self, u, v):

        return u in self.ngbh[v]

    def copy(self):

        G_copy = BipartiteGraph()

        G_copy.n_nodes = copy.copy(self.n_nodes)
        G_copy.ngbh = {copy.copy(key):[e for e in val] for key, val in self.ngbh.items()} 
        G_copy.left = [copy.copy(v) for v in self.left]
        G_copy.right = [copy.copy(v) for v in self.right]
        G_copy.side = {copy.copy(key):val for key, val in self.side.items()}
    
        return G_copy 

class Digraph:

    def __init__(self):
        """
        Initializes attributes (n_nodes, in_nghbh, out_ngbh) to empty/zero default values
        """

        #: (**integer**) - initialized at 0 by __init__. Number of nodes in the graph
        self.n_nodes = 0

        #: (**dict**) - initialized at {} by __init__. In-adjacency dictionary. Mapping vertices (which are integers) to the list of their in-neighbors.
        self.in_ngbh = {}

        #: (**dict**) - initialized at {} by __init__. In-adjacency dictionary. Mapping vertices (which are integers) to the list of their in-neighbors.
        self.out_ngbh = {}

    def add_node(self, u):
        """
        Method to add a node to the graph.

        :param u: integer describing the new vertex. If already in the graph nothing happens
        :type u: **int**
        """

        try:
            self.in_ngbh[u]
        except KeyError:
            self.n_nodes += 1
            self.in_ngbh[u] = []
            self.out_ngbh[u] = []

    def add_edge(self, u, v):
        """
        Method to add an edge to the graph, from u to v. If one of the vertices does not exist, it is added to the graph.

        :param u: integer describing a vertex, source of the edge.
        :type u: **int**

        :param v: integer describing a vertex, tip of the edge.
        :type v: **int**
        """

        if u == v:
            return
        try:
            self.in_ngbh[v].append(u)
            self.out_ngbh[u].append(v)
        except KeyError:
            self.add_node(u)
            self.add_node(v)
            self.in_ngbh[v].append(u)
            self.out_ngbh[u].append(v)

    def subgraph(self, seq):

        sub_G = Digraph()

        for u in seq:
            sub_G.add_node(u)

        for u in seq:
            for v in self.in_ngbh[u]:
                try:
                    sub_G.in_ngbh[v]
                    sub_G.add_edge(v, u)
                except KeyError:
                    pass

        return sub_G
