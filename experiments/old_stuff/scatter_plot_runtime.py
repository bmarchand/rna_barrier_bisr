#!/usr/bin/env python
# coding: utf-8

# In[1]:


from rna_barrier_subtree_solver.subtree_solver import solve_subtree
from rna_barrier_subtree_solver.utilities import num_leaves, filter_common_bps, list_bps
from bisr_dpw.barriers_py3 import ssrandom
import time
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


color_list = sns.color_palette('hls',8)


# In[2]:


xs = []
ys = []
cs = []

theta = 1
unpaired_weight=0.1

params = set([])

for n in range(10,60,4):
    for inc in range(60):
        s1 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=inc)
        s2 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=inc+n)
        s1, s2 = filter_common_bps(s1, s2)
        print(s1,'->',s2)
        t0 = time.time()
        _, _ = solve_subtree(s1, s2)
        run_time = time.time() - t0
        param = min(num_leaves(s1), num_leaves(s2))
        assert param!=1, "phi=1 for "+s1+'->'+s2
        if param==0:
            continue
        if param==num_leaves(s1):
            xs.append(len(list_bps(s1)))
        else:
            xs.append(len(list_bps(s2)))
        
        ys.append(run_time)
        cs.append(color_list[param-2])
        params.add(param)


# In[3]:


print(params)
plt.figure(dpi=200)
plt.scatter(xs, ys, c=cs, s=[8 for _ in xs])
plt.yscale('log')
custom_lines = [Line2D([0], [0], color=c, lw=4) for c in color_list[:max(params)-2]]

plt.xlabel('number of base pairs in initial structure')
plt.ylabel('run-time (seconds)')
plt.legend(custom_lines, [r'$\phi$='+str(param) for param in range(min(params), max(params)+1,1)])
plt.show()


# In[22]:


import rna_barrier_subtree_solver.subtree_solver as subtree_solver
subtree_solver.re_init_dp_table()
from bisr_dpw.dpw_interface import solve_dpw

ts_dpw = []
ts_subtree = []

theta = 3
unpaired_weight=0.5

for n in range(10,50,5):
    for inc in range(60):
        s1 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=inc)
        s2 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=inc+n)
        s1, s2 = filter_common_bps(s1, s2)
        print(s1,'->',s2)
        t0 = time.time()
        k1, _ = solve_subtree(s1, s2)
        ts_subtree.append(time.time() - t0)
        
        t0 = time.time()
        k2, _ = solve_dpw(s1, s2)
        ts_dpw.append(time.time()-t0)
        assert(k1==k2)
        


# In[23]:


plt.scatter(ts_dpw, ts_subtree)
plt.xlabel('dpw solve time')
plt.ylabel('subtree method solve time')
M = max(max(ts_dpw), max(ts_subtree))
plt.plot([0,M],[0,M], c='r')


# In[17]:


def biclique(n):
    s1 = n*'('+n*'.'+n*'('+n*'.'
    s2 = n*'.'+n*'('+n*'.'+n*')'
    print(s1)
    print(s2)
    return s1,s2


# In[20]:


ns = range(1,200,10)
ts_subtree = []
ts_dpw = []
for n in ns:
    s1, s2 = biclique(n)
    
    t0 = time.time()
    _, _ = solve_subtree(s1, s2)
    ts_subtree.append(time.time() - t0)
    
    t0 = time.time()
    _, _ = solve_dpw(s1, s2)
    ts_dpw.append(time.time()-t0)


# In[21]:


plt.plot(ns, ts_dpw,'o-',label='dpw')
plt.plot(ns, ts_subtree,'o-',label='subtree')
plt.legend()
plt.show()


# In[ ]:




