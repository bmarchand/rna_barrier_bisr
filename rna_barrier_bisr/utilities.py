
def list_bps(s):
    stack = []
    bps = []
    for k, c in enumerate(s):
        if c=='(':
            stack.append(k)
        if c==')':
            bps.append((stack.pop(), k))
    return bps

def num_leaves(s):
    """
        s: structure given as dot bracket notation

        the number of leaves is the number of times
        the pattern (...) (with an arbitrary number
        of points) appears in the structure,
        plus an additional +1 if the root
        of the tree is of degree 1.
    """
    # number of leaves
    n_leaves = s.replace('.','').count('()')
    if n_leaves==0:
    # case of an empty struct
        return 0

    # +1 for the root if deg(root)==1
    for k in range(1,len(s)-1,1):
        if s[k]=='.':
            if s[:k-1].count('(') > 0 and s[k+1:].count('(')==s[k+1:].count(')') > 0:
                if s[:k-1].count('(')==s[:k-1].count(')') and s[k+1:].count('(')==s[k+1:].count(')'):
                    return n_leaves

    return n_leaves+1

def filter_common_bps(s1, s2):
    bps1 = list_bps(s1)
    bps2 = list_bps(s2)

    s1 = list(s1)
    s2 = list(s2)

    for bp in bps1:
        if bp in bps2:
            i,j = bp
            s1[i] = '.'
            s1[j] = '.'
            s2[i] = '.'
            s2[j] = '.'

    return "".join(s1), "".join(s2)

if __name__=='__main__':
    # quick checks
    assert(num_leaves('(...((..))..((..))...)')==3)
    assert(num_leaves('...((..))..((..))...')==2)
