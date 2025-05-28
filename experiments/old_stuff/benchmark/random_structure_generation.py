from bisr_dpw.barriers_py3 import ssrandom
from rna_barrier_subtree_solver.utilities import num_leaves, filter_common_bps, list_bps

theta = 1
unpaired_weight=.1

cnt = 0

f = open('random_structures_data_set.csv','w')
print('index,rna_length,s1,s2,arboricity',file=f)

for n in range(10,60,5):
    for inc in range(100):
        s1 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=2*cnt)
        s2 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=2*cnt+1)
        s1, s2 = filter_common_bps(s1, s2)
        param = min(num_leaves(s1), num_leaves(s2))
        print(','.join(map(str,[cnt,n,s1,s2,param])), file=f)
        cnt += 1

f.close()
