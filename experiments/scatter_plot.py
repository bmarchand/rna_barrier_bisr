from rna_barrier_subtree_solver.subtree_solver import solve_subtree
from rna_barrier_subtree_solver.utilities import num_leaves, filter_common_bps
from bisr_dpw.barriers_py3 import ssrandom
import time
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


color_list = sns.color_palette('hls',8)

xs = []
ys = []
cs = []

theta = 1
unpaired_weight=0.3

params = set([])

for n in range(10,50,2):
    for inc in range(50):
        s1 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=inc)
        s2 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=inc+n)
        s1, s2 = filter_common_bps(s1, s2)
        print(s1,'->',s2)
        t0 = time.time()
        k3, schedule3 = solve_subtree(s1, s2)
        run_time = time.time() - t0
        param = min(num_leaves(s1), num_leaves(s2))
        if param==0:
            continue

        xs.append(n)
        ys.append(run_time)
        cs.append(color_list[param-2])
        params.add(param)

plt.scatter(xs, ys, c=cs)

custom_lines = [Line2D([0], [0], color=c, lw=4) for c in color_list[:max(params)-2]]

plt.legend(custom_lines, [r'$\phi$='+str(param) for param in range(min(params), max(params)+1,1)])
plt.show()
