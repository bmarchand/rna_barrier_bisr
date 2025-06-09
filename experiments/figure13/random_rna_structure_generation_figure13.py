from rna_barrier_bisr.random_rna_structures import ssrandom
from rna_barrier_bisr.utilities import num_leaves, filter_common_bps, list_bps

theta = 1
unpaired_weight=.1

cnt = 0

f = open('random_structures_data_set_figure13.csv','w')
print('index,rna_length,s1,s2,arboricity,num_vertices_bipartite_graph',file=f)

for n in range(10,100,5):
    for inc in range(300):
        s1 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=2*cnt)
        s2 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=2*cnt+1)
        s1, s2 = filter_common_bps(s1, s2)
        num_vertices_bip_graph = s1.count('(')+s2.count('(')
        param = min(num_leaves(s1), num_leaves(s2))
        print(','.join(map(str,[cnt,n,s1,s2,param,num_vertices_bip_graph])), file=f)
        cnt += 1

f.close()
