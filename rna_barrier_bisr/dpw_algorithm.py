from bisr_dpw_cpp_routines import strongly_connected_components
import networkx as nw
from .graph_classes import BipartiteGraph, Digraph


class Node:
    def __init__(self, value, int_attribute, in_ngbh=set([])):
        """
        Args:
            - value: will store the label of a node
            - int_attribute: meant to store the in-degree of the sequence
            represented by the node, in implementation of DPW below.
        """

        self.children = []  # lexicographically (< only) sorted
        self.parent = None
        self.int_attribute = int_attribute
        self.in_ngbh = in_ngbh  # The set of in-ngbh not in seq=prefix+node.val
        self.value = value


class Trie:

    def __init__(self):

        self.root = Node(-1, 0)

    def insert(self, seq):

        node = self.root

        for val in seq:
            found = False
            for child in node.children:
                if val == child.value:
                    found = True
                    node = child
                    break

            if found:
                continue

            new_node = Node(val, -1)

            wh_to_insert = 0
            for ind, child in enumerate(node.children):
                if child.value < val:
                    wh_to_insert = ind+1
                    break

            node.children.insert(wh_to_insert, new_node)
            node = new_node

    def breadth(self):

        def dfs(node):
            if len(node.children)==0:
                return 1

            n_leaves = 0
            for child in node.children:
                n_leaves += dfs(child)

            return n_leaves

        return dfs(self.root)

    def print(self):

        def dfs(node, prefix):

            if len(node.children) == 0:
                print(prefix + [node.value])

            for child in node.children:
                dfs(child, prefix + [node.value])

        dfs(self.root, [])

    def generate_sequences(self, depth_limit=None):

        sequences = []

        def dfs(node, prefix, depth):

            if len(node.children) == 0:
                if depth_limit:
                    if depth >= depth_limit:
                        sequences.append(prefix + [node.value])
                else:
                    sequences.append(prefix + [node.value])

            for child in node.children:
                dfs(child, prefix + [node.value], depth + 1)

        dfs(self.root, [], 0)

        return sequences


def seq_to_int(seq):

    if seq[0] == -1:
        seq = seq[1:]

    n = ""
    for s in seq:
        n += str(s)

    return int(n)


def DPW(G, k, full_results=False, break_into_scc=True):
    """
    Implementation of Tamaki's algorithm for directed pathwidth.

    :param G: Digraph instance, directed graph to solve.
    :type G: Digraph

    :param k: directed pathwidth decision value: the function must decide whether the input Digraph has directed pathwidth <= k or not.
    :type k: int

    :param full_results: The first step of the algorithm is (unless break_into_scc is set to False) to decompose the input Digraph into its strongly connected components, which are solved separately and independently. when set to **True**, it returns a list of individual solutions for each SCC. defaults to **False**.
    :type full_results: bool

    :param break_into_scc: Boolean for specifying whether the input Digraph should be broken into strongly connected components, which would then be solved independently. Such a pre-processing step can always be carried out. It was made specifiable for debug purposes. defaults to **True**
    :type break_into_scc: bool

    :return: **seq**: either equal to **None** if the input Digraph has directed pathwidth > k or to a valid layout, in the form of a **list** of vertices.
    :rtype: **list** or **None**
    """

    if not break_into_scc:
        return DPW_scc(G, k, full_results=full_results)

    sccs = strongly_connected_components(G.out_ngbh, G.in_ngbh)

    total_seq = []
    for scc in sccs:
        sub_G = G.subgraph(scc)
        seq = DPW_scc(sub_G, k, full_results=full_results)

        if seq is None:
            return None

        else:
            if full_results:
                total_seq += [seq]
            else:
                total_seq += seq[1:]

    return total_seq


def DPW_scc(G, k, full_results=False):


    # DFS AUXILIARY FUNCTIONS
    def dfs1(node, prefix, depth):
        """
        DFS function whose role is to generate all immediate extensions.
        (T_i set in paper)
        They are added to the trie. Some will potentially be removed
        in later steps.
        """

        # will we find at least one extension ?
        found_immediate_extension = False

        # If reached bottom of tree: generate immediate extensions
        if len(node.children) == 0 and depth == i-1:

            seq = prefix + [node.value]  # the sequence will will extend.

            banned = {}  # to quickly check, afterwards, what is in seq.
            for u in seq:
                banned[u] = True

            # loop over all vertices. tagging whether \in seq or not.
            for v in G.in_ngbh.keys():
                try:
                    banned[v]
                except KeyError:
                    banned[v] = False

            # loop for generating immediate extensions
            for v in G.in_ngbh.keys():

                # not banned = not in seq = we can extend seq with v
                if not banned[v]:

                    # collecting in-ngbhs of v from V \ seq
                    in_ngbh = set([])
                    for w in G.in_ngbh[v]:
                        if not banned[w]:
                            in_ngbh.add(w)

                    # form in-ngbh of extended sequence
                    in_ngbh = in_ngbh.union(node.in_ngbh)

                    # if v was an in-ngbh of seq remove it
                    try:
                        in_ngbh.remove(v)
                    except KeyError:
                        pass

                    # num of in-ngbh of new seq
                    new_seq_deg = len(in_ngbh)

                    # condidition for immediate extension acceptance
                    if new_seq_deg <= k:
                        new_node = Node(v, new_seq_deg, in_ngbh=in_ngbh)

                        # need to figure out where to insert (lexico-order)
                        wh_to_insert = 0
                        for ind, child in enumerate(node.children):
                            wh_to_insert = ind+1
                            if child.value > new_node.value:
                                wh_to_insert = ind
                                break

                        # actual inserting
                        node.children.insert(wh_to_insert, new_node)
                        new_node.parent = node

                        # book keeping
                        found_immediate_extension = True

            return found_immediate_extension

        else:
            # haven't reach leaves. recursive calls.
            for child in node.children:
                found_here = dfs1(child, prefix + [node.value], depth + 1)
                found_immediate_extension = found_immediate_extension or \
                    found_here

            return found_immediate_extension


    def dfs2(node, prefix, depth):
        """
        DFS function, which stops depth exploration at level (i-1), i.e
        just before the leaves.

        At this level, for each eta (paper notation), we look among the
        children for most preferable non-expanding immediate extension
        (definition of T_i').

        If we find one, we set tagged[eta] = True,
        and tagged[extension_in_question] = True as well.

        Otherwise, we set tagged[eta] = False. Same for all its immediate
        extensions.
        """

        # stoping when depth counter reaches level just before leaves.
        if depth == i-1 and i > 1:
            seq = prefix + [node.value]

            tagged[node] = False

            min_degree = 10**9
            # among non-expanding children of \eta, choose the best.
            for child in node.children:
                # if non-expanding.
                if child.int_attribute <= node.int_attribute:
                    tagged[node] = True
                    if child.int_attribute < min_degree:
                        # strictness of "<" + ordering
                        # ensures selection of most preferable
                        min_degree = child.int_attribute
                        lucky_child = child

            if tagged[node]:
                # then lucky child is defined, and tagged as well.
                tagged[lucky_child] = True

            # registering ancestor/inheritor:
            if tagged[node]:
                sig_degree = lucky_child.int_attribute
                eta_degree = node.int_attribute

                # node is a valid ancestor, with lucky_child as inheritor
                # but we want to find the shortest prefix such that this
                # is True. Two cases are distinguished: sig_degree=eta_deg
                # and sig_degree < eta_degree

                # in this case the ancestor can only be eta.
                if sig_degree == eta_degree:
                    furthest_ancestor = node

                # else: eta is valid, but iterate further back for other
                # possibilities.
                else:

                    ancestor = node
                    furthest_ancestor = node
                    min_since_eta = eta_degree
                    while ancestor.parent.value != -1:
                        ancestor = ancestor.parent

                        # if lower than sigma, break
                        if ancestor.int_attribute < sig_degree:
                            break

                        # higher than sigma, but lower than anything seen
                        if ancestor.int_attribute < min_since_eta:
                            min_since_eta = ancestor.int_attribute
                            furthest_ancestor = ancestor

                # register inheritor. There is always one here.
                try:
                    # if ancestor already had inheritor, compare
                    other_child = inheritor[furthest_ancestor]
                    if lucky_child.int_attribute < \
                            other_child.int_attribute:
                        # if new is better then:
                        inheritor[furthest_ancestor] = lucky_child
                except KeyError:
                    # else, just register.
                    inheritor[furthest_ancestor] = lucky_child

        else:
            # haven't reached bottom, recursive call
            for child in node.children:
                dfs2(child, prefix + [node.value], depth + 1)

    def dfs3(node, prefix, depth):
        """
        DFS to identify and remove suppressed elements.

        When we reach level i-1, just before the leaves, we stop
        and look for suppressed elements.
        """

        # stopping at level just before leaves.
        if depth == i-1 and i > 1:
            seq = prefix + [node.value]

            if not tagged[node]:
                # all or nothing case: if a prefix has an inheritor,
                # drop everything
                # else keep everything (all immediate extensions)

                # iterating up-tree over prefixes to see if ancestor.
                ancestor = node

                found_ancestor = False
                while ancestor.value != -1:
                    try:
                        inheritor[ancestor]
                        found_ancestor = True
                        break
                    except KeyError:
                        ancestor = ancestor.parent

                # if found ancestor, prune all leaves.
                if found_ancestor:
                    node.children = []

            # else \eta has a child in T_i'.
            else:
                # identifying tagged child:
                for child in node.children:
                    try:
                        if tagged[child]:
                            lucky_child = child
                    except KeyError:
                        continue
                        # there should be a tagged child at some point.

                # children are pruned, will see if later if lucky_child
                # added
                node.children = []

                ancestor = node
                while ancestor.value != -1:
                    try:
                        if inheritor[ancestor] == lucky_child:
                            if len(node.children) == 0:
                                node.children.append(lucky_child)
                        else:
                            node.children = []
                            break
                        ancestor = ancestor.parent
                    except KeyError:
                        ancestor = ancestor.parent

        else:
            # have not reached level i-1, going deeper.
            for child in node.children:
                dfs3(child, prefix + [node.value], depth + 1)

    def dfs4(node, prefix, depth):
        """
        Final DFS to return solutions if there are any.
        """

        if depth == G.n_nodes:
            return prefix + [node.value]

        for child in node.children:
            seq = dfs4(child, prefix + [node.value], depth+1)
            if seq:
                return seq


    ## ACTUAL ALGORITHM ##

    trie = Trie()  # Will store all the S_i \forall i.

    i = 0

    while i < G.n_nodes:  # purpose is to try and generate S_n

        i += 1

        # STEP 1: immediate extensions
        found_immediate_extension = dfs1(trie.root, [], 0)

        if not found_immediate_extension:
            return None

        # STEP 2: shortest non-expanding extensions, ancestors, inheritors.
        # (tagging elements of T_iprime, recording ancestors)
        # storage of Boolean for elements of T_i' and their parent.

        tagged = {}
        inheritor = {}  # mapping ancestors to inheritors.

        dfs2(trie.root, [], 0)

        # STEP 3: FILTERING OUT SUPPRESSED ELEMENT.
        dfs3(trie.root, [], 0)

    # final dfs to build solution if there is one
    seq = dfs4(trie.root, [], 0)

    if not full_results:
        return seq
    else:
        return seq, trie


def valid_layout(G, seq, k):

    valid = True

    if seq[0] == -1:
        seq = seq[1:]

    for i in range(1, G.n_nodes+1, 1):

        sigma = seq[:i]

        in_neighbors = set([])

        for u in sigma:
            for v in G.in_ngbh[u]:
                if v not in sigma:
                    in_neighbors.add(v)

        if len(in_neighbors) > k:
            valid = False

    return valid
