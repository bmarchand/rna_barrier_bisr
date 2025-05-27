from bisr_dpw.rna_interface import random_conflict_graph, structures_to_conflict_graph, solve_mcf
from bisr_dpw.dpw_interface import construct_digraph, solve_dpw
from bisr_dpw.barriers_py3 import realize, ssrandom
import random

import matplotlib.pylab as plt

seed = 2021
theta = 1
unpaired_weight = .5

def test_agreement_dpw_mcf_small_graphs():

    ns = [10, 15, 20]

    nruns = 10

    for n in ns:
        for inc in range(nruns):
            s1 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=inc)
            s2 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=inc)
            k1 ,schedule = solve_dpw(s1, s2)
            k2 ,schedule = solve_mcf(s1, s2)
            assert(k1==k2)

def test_small_instance1():

    n = 20
    s1 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=seed)
    s2 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=seed)

    k1 ,schedule = solve_dpw(s1, s2)
    k2 ,schedule = solve_mcf(s1, s2)

    assert(k1==k2)

def test_small_instance2():

    n = 20
    s1 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=seed+1)
    s2 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=seed+1)

    k1 ,schedule = solve_dpw(s1, s2)
    k2 ,schedule = solve_mcf(s1, s2)

    assert(k1==k2)

def test_small_instance3():

    n = 20

    s1 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=seed+62)
    s2 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=seed+62)
    k1 ,schedule = solve_dpw(s1, s2)
    k2 ,schedule = solve_mcf(s1, s2)

    assert(k1==k2)

def test_small_instance4():

    n = 20

    s1 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=seed+62)
    s2 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=seed+62)
    k1 ,schedule = solve_dpw(s1, s2)
    k2 ,schedule = solve_mcf(s1, s2)

    assert(k1==k2)


def test_small_instance5():

    s1 = '.(.(...)).'
    s2 = '(..(...).)'

    # solve dpw
    k1 ,schedule = solve_dpw(s1, s2)

    # solve mcf
    k2 ,schedule = solve_mcf(s1, s2)
    
    assert(k1==k2)

def test_small_instance6():
    s1 =".((...(((.)))))"
    s2 ="((((.)((.).))))"

    # solve dpw
    k1 ,schedule = solve_dpw(s1, s2)

    # solve mcf
    k2 ,schedule = solve_mcf(s1, s2)
    
    assert(k1==k2)

def test_small_instance7():
    s1 = "...(.)...."
    s2 = ".(..(.).)."

    # solve dpw
    k1 ,schedule = solve_dpw(s1, s2)

    # solve mcf
    k2 ,schedule = solve_mcf(s1, s2)
    
    assert(k1==k2)

def test_small_instance8():
    s1 = "((((..))))"
    s2 = "(((....)))"

    # solve dpw
    k1 ,schedule = solve_dpw(s1, s2)

    # solve mcf
    k2 ,schedule = solve_mcf(s1, s2)
    
    assert(k1==k2)
    

def test_small_instance9():
    s1 = ".((...(((.)))))"
    s2 = "((((.)((.).))))"

    # solve dpw
    k1 ,schedule = solve_dpw(s1, s2)

    # solve mcf
    k2 ,schedule = solve_mcf(s1, s2)
    
    assert(k1==k2)

def test_small_instance10():
    s1 = ".(.(...))."
    s2 = "((..)(..))"

    # solve dpw
    k1 ,schedule = solve_dpw(s1, s2)

    # solve mcf
    k2 ,schedule = solve_mcf(s1, s2)
    
    assert(k1==k2)


if __name__=='__main__':
    test_small_instance7()
