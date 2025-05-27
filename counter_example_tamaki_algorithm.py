from rna_barrier_bisr.dpw_algorithm import DPW
from rna_barrier_bisr.graph_classes import Digraph
from itertools import groupby
  

"""
Counter example for DPW/Tamaki

Target solution (S)
A5   B4  C7     D5   X3 Y3 Z4  E3
 _               _          _           
/ \       ___   / \     /\ / \         
|  \  ___/   \  |  \  __| \|  \        
|   \/        \ |   \/         \       
|              \|               \  
|                                \ 
          x Kill                       
S is pruned at branch C7 by an other branch T below (SNEKFE Y shorter than C)

T::
A5   B4  Y3 Z4  C7     D5   X3 E3
 _           _          _                      
/ \      /\ / \  ___   / \          
|  \  ___| \|  \/   \  |  \  __       
|   \/               \ |   \/  \ 
|                     \|        \ 
|                                \ 
                 x Kill            
T is pruned at branch B4 by an other branch U (extension XBC is shorter than BYZC, but XBC is not a SNEKFE: BC is the true "lost" SNEKFE )

U::
A5   X3 B4  C7    
 _           ___      
/ \      ___/   \ 
|  \  __/        \ 
|   \/            \. 
|               
U does not lead to a solution: there is no way to place D5 after X3
NB: here U is not killed by a branch with YZ, since there is not enough budget left after AXB (and YZ can only appear after B)

Constraint on the ordering of the bumps:
        >Y->Z 
       /      \ 
   A->B->C->D->E
    \         / 
      -> X ->


When running the script the following happens: a graph whose solution space consists of bumps, with ordering
constraints outlined as above, is created. It has dpw=5. Yet, launching the algorithm (DPW function)
on it returns false (for the reasons described above). Adding a single edge, we force C7 to be before
Z4, which (here, for this graph) does not change the dpw, but makes Tamaki's algorithm find it.
"""


    
def bump(H,F, up, flat, down, names, thisname):
    pref=[]
    out=[]
    if up>0:
        flat+=1
    for _ in range(up):
        x=H.n_nodes
        F.append(x)
        H.add_node(x)
        for y in pref:
            H.add_edge(y,x)
            H.add_edge(x,y)
        pref.append(x)
    for _ in range(flat):
        x=H.n_nodes
        H.add_node(x)
        for y in pref:
            H.add_edge(y,x)
        pref.append(x)
        out.append(x)
        names[x]=thisname 
    for _ in range(down):
        assert(F)
        x=F.pop(-1)
        for y in pref:
            H.add_edge(y,x)
        pref.append(x)
        out.append(x)
        names[x]=thisname
    return out 

def linklist(H, sets):
    pref=[]
    for s in sets:
        for p in pref:
            for u in s:
                H.add_edge(p,u)
        pref.extend(s)

def counterex(addExtraEdge, verbose=False):
    H = Digraph()
    k=5
    Free=[]
    Names={}
    A=bump(H,Free,5,1,3,Names,'A')
    B=bump(H,Free,1,3,0,Names,'B')
    C=bump(H,Free,1,3,3,Names,'C')
    D=bump(H,Free,4,1,3,Names,'D')
    X=bump(H,Free,1,2,0,Names,'X')
    E=bump(H,Free,0,0,3,Names,'E')
    Y=bump(H,Free,2,0,2,Names,'Y')
    Z=bump(H,Free,2,1,2,Names,'Z')
    assert(not Free)
    linklist(H, [A,B,C,D,E])
    linklist(H, [A,X,E])
    linklist(H, [B,Y,Z,E])
    if addExtraEdge:
        H.add_edge(12,32)
        print('( add an edge ',Names[12],'->', Names[32],')')
    if verbose:
        print(Names)

    
    result = DPW(H, k, full_results=True, break_into_scc=True)      
    
    if result is not None:
        print(result)
        seq,trie=result[0] 
        if verbose:
            for i, trie in enumerate(rest):
                print("level i =", i+1)
                current=set()   
                for smallseq, val in trie:         
                    current.add(str(val)+' '+' '.join( [x+':'+str(sum(1 for _ in group)) for x, group in groupby([Names[i] for i in smallseq[1:]])]))
                print('\n'.join(list(current)))



        print(f"Found sequence (k={k}): ", seq[1:])
        count_dups = [x+str(sum(1 for _ in group)) for x, group in groupby([Names[i] for i in seq[1:]])]
        print (' '.join(count_dups))
    else:
        print(f'No solution (k={k})')

print('** graph H: ')
counterex(False, verbose=False)
print()
print('** graph H + one edge: ')
counterex(True, verbose=False)
