from rna_barrier_subtree_solver.subtree_solver import solve_subtree
from bisr_dpw.barriers_py3 import ssrandom
from rna_barrier_subtree_solver.utilities import num_leaves, filter_common_bps

theta = 4
unpaired_weight=0.3
n = 100
seed = 2023

s1 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=seed)
s2 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=seed+n)
s1, s2 = filter_common_bps(s1, s2)

k, sched = solve_subtree(s1, s2)

print(k,sched)
