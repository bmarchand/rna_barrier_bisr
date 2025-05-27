from rna_barrier_subtree_solver.utilities import num_leaves

def test_num_leaves1():
    s = '.....(.).'
    assert(num_leaves(s)==2)

def test_num_leaves2():
    s = '(..)..(..).'
    assert(num_leaves(s)==2)

def test_num_leaves3():
    s = '...........'
    assert(num_leaves(s)==0)
