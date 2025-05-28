import random
import sys
import os
import getopt
import time
from .random_rna_structures import sstopairs, ssparse, compatiblepairs, pairstoss

DEBUG=False
SMALLDEBUG = True
LEAVEINCLUDES = True
LEAVEBLACKLIST = True
RUNBRUTEFORCE = False
TIME = False

def boolTuples(dim):
    if dim <= 0:
        return [()]
    else:
        res = []
        for t in boolTuples(dim-1):
            for b in [True, False]:
                res.append(tuple([b]+list(t)))
        return res



def printpath(p, B, R):
    bpsB = sstopairs(B)
    bpsR = sstopairs(R)
    tmp = [c for c in ssstring(B)]
    h = 0
    print("".join(tmp), "h=", h)
    for (i, j) in p:
        if (i, j) in bpsB:
            (tmp[i], tmp[j]) = ('.', '.')
            h += 1
        if (i, j) in bpsR:
            (tmp[i], tmp[j]) = ('[', ']')
            h -= 1
        print("".join(tmp), "h=", h)


def IInverse(Sep, B, R):
    # Function I is its own inverse
    return I(Sep, B, R)


def I(Sep, B, R):
    bpsB = sstopairs(B)
    bpsR = sstopairs(R)
    return (bpsB-Sep, bpsR & Sep)

def MCF(B, R, theta=1):
    # Returns max #pairs, and matching structure in well-parenthesized
    # notation, composed of compatible pairs in B + R (=E)
    # restricted to interval [i,j], and further constrained to
    # featuring >0 pair for B (iff alpha=True) or R (iff beta=True)
    # mat is just a caching dictionary for memoization
    # Note: It takes strength and discipline to be *this* lazy...
    mat = {}
    n = len(B)
    # Init for intervals smaller than THETA nucleotides
    for m in range(1, theta+1):
        for alpha, beta in boolTuples(2):
            for i in range(0, n-m+1):
                j = i+m-1
                if (alpha, beta) == (False, False):
                    # If interval too small, and no need to include
                    # pairs from either B or R (alpha,beta == False, False)
                    # Then the best structure has 0 pairs
                    mat[(i, j, alpha, beta)] = (0, "."*m)
                else:
                    # If interval too small, but we still need to include
                    # some pair from B or R (alpha,beta != False, False)
                    # Then there is no allowed structure (-infty pairs)
                    mat[(i, j, alpha, beta)] = (-sys.maxsize, "X"*m)
    for m in range(theta+1, n+1):
        for alpha, beta in boolTuples(2):
            for i in range(0, n-m+1):
                j = i+m-1

                # Case 1: First position unpaired, max #pairs is found by
                # recurrence on [i+1,j], forwarding the need to feature pair
                # from B (iff alpha=True) or R (iff beta=True)
                (k, ss) = mat[(i+1, j, alpha, beta)]
                mat[(i, j, alpha, beta)] = (k, "."+ss)

                # Case 2: First position paired, max #pairs is found by
                # recurrence on [i+1,j], forwarding the need to feature pair
                # from B (iff alpha=True) or R (iff beta=True)

                # Check if pairs from B and R at pos i can be considered
                # for the paired case (ie must be fully contained in [i,j])
                partners = []
                if B[i] > i and B[i] <= j:
                    partners.append((i, B[i]))
                if R[i] > i and R[i] <= j:
                    partners.append((i, R[i]))
                # Go over suitable partners for pos. i (within R or B)
                for (a, k) in partners:
                    assert(a == i)
                    # Consider all possible ways of "forwarding" constraints
                    # (having >0 pair from B, ie alpha=True, or R, ie beta=True),
                    # to subintervals [i+1,k-1] (alpha1,beta1) and [k+1,j] (alpha2,beta2).
                    for (alpha1, beta1, alpha2, beta2) in boolTuples(4):
                        # If alpha=True, check that "someone" satisfies constraint
                        # of having >0 pair from B (ie either through pair (a,k),
                        # recursion over [i+1,k-1], or rec. over [k+1,k])
                        if not alpha or (alpha1 or alpha2 or B[i] == k):
                            # Same for R w.r.t. beta1 and beta2
                            if not beta or (beta1 or beta2 or R[i] == k):
                                # Recurse over subintervals with suitable constraints
                                if i+1 <= k-1:
                                    (termA, ssA) = mat[(
                                        i+1, k-1, alpha1, beta1)]
                                elif not alpha1 and not beta1:
                                    (termA, ssA) = (0, "")
                                else:
                                    (termA, ssA) = (-sys.maxsize, "X")
                                if k < j:
                                    (termB, ssB) = mat[(
                                        k+1, j, alpha2, beta2)]
                                elif not alpha2 and not beta2:
                                    (termB, ssB) = (0, "")
                                else:
                                    (termB, ssB) = (-sys.maxsize, "X")
                                # Accumulate max #pairs, plus one for current pair (a,k)
                                mat[(i, j, alpha, beta)] = max(
                                    mat[(i, j, alpha, beta)], (1+termA+termB, "("+ssA+")"+ssB))
    return mat[(0, n-1, True, True)]


def separate(B, R, theta=1):
    assert(len(B) == len(R))
    # Call to rec. function, considering full interval [i,j] = [0,n-1],
    # and forcing >0 pair from both B and R with (alpha,beta) = (True,True)
    (k, sep) = MCF(B, R, theta=theta)
    if DEBUG:
        print("            MCF", k, sep)

    # I constraints cannot be satisfied, return nothing
    if k < 0:
        return None
    # Convert well-parenthesized separator struct into arc-annotated seq...
    ss = ssparse(sep)
    # ... and then to set of base pairs
    bps = sstopairs(ss)
    # Check that the prediction max #pairs and effective #pairs match
    assert(k == len(bps))
    if k < max(len(sstopairs(B)), len(sstopairs(R))):
        return None
    else:
        return IInverse(bps, B, R)


def delta(B, R):
    return len(B)-len(R)


def builddependencies(B, R):
    res = []
    for p1 in R:
        res.append((p1, [p2 for p2 in B if not compatiblepairs(p1, p2)]))
    tmp = sorted([(len(deps), p, deps) for p, deps in res])
    # print tmp
    return [(p, set(deps)) for (m, p, deps) in tmp]


def skipInclusions(allDeps):
    # assume allDeps is sorted by deps len,
    # filter out (p1,d1) if it contains (p2,d2)
    # (can be assumed to not appear first)
    if LEAVEINCLUDES:
        return allDeps
    res = []
    for (p1, deps1) in allDeps:
        for (p2, deps2) in res:
            if deps2.issubset(deps1):
                break
        else:
            res.append((p1, deps1))
    return res


def skipBlacklist(bl, allDeps):
    if LEAVEBLACKLIST:
        return allDeps
    res = []
    for (p, deps) in allDeps:
        if p not in bl:
            res.append((p, deps))
    return res


def no_conflicts(bpsB, bpsR):
    for i, j in bpsB:
        for k, l in bpsR:
            if (i <= k and k <= j and j <= l) or (k <= i and i <= l and l <= j):
                return False
    return True

def list_deps(bpR, bpsB):
    i,j = bpR
    dependencies = []
    for k,l in bpsB:
        if (i <= k and k <= j and j <= l) or (k <= i and i <= l and l <= j):
            dependencies.append((k,l))
    return dependencies

def realize(B, R, k, depth=0, firstRblacklist=set([]), lastBblacklist=set([]), theta=1):
    def indent(depth):
        return "  "*depth

    if DEBUG:
        print(indent(depth)+"realize (k=%s):" % k)
        print(indent(depth+1)+"B=%s" % (ssstring(B)))
        print(indent(depth+1)+"R=%s" % (ssstring(R)))

    n = len(B)
    assert(n == len(R))
    if k < 0:
        return None

    bpsB = sstopairs(B)
    bpsR = sstopairs(R)

    if len(bpsB) == 0:
        if DEBUG:
            print(indent(depth+1)+"|B|=0   -> Return R")
        return list(bpsR)

    for bpR in bpsR:
        if len(list_deps(bpR, bpsB))==0:
            stock = realize(B, 
                            pairstoss(bpsR-set([bpR]), n), 
                            k+1, 
                            depth=depth+1, 
                            firstRblacklist=firstRblacklist, 
                            lastBblacklist=lastBblacklist, 
                            theta=theta)
            if stock:
                return [bpR]+stock
            else:
                return None #Case where only red nodes remain already treated before

#    if no_conflicts(bpsB, bpsR):
#        return list(bpsR)+list(bpsB)

    if k == 0:
        return None

    if len(bpsR) == 0:
        if DEBUG:
            print(indent(depth+1)+"|R|=0   -> Return B")
        if len(bpsB) > k:
            return None
        else:
            return list(bpsB)

    if len(bpsB) <= k:
        if DEBUG:
            print(indent(depth+1)+"|B|<=k   -> Return B.R")
        return list(bpsB) + list(bpsR)

    if len(bpsB) < len(bpsR):
        if DEBUG:
            print(indent(depth+1)+"|B|<|R| -> Swap B<->R")
        k += (len(bpsR)-len(bpsB))
        res = realize(R, B, k, depth, lastBblacklist,
                      firstRblacklist, theta=theta)
        if res is None:
            return None
        else:
            return list(reversed(res))

    # Check if separator can be found to divide and conquer
    sep = separate(B, R, theta=theta)
    if sep is not None:
        (bpsBp, bpsRp) = sep
        (bpsBi, bpsRi) = I(bpsBp | bpsRp, B, R)
        if DEBUG:
            print(indent(depth+1)+"Sep!=0  -> I=",
                  ssstring(pairstoss(bpsBi | bpsRi, n)))
        ssBLeft = pairstoss(bpsBp, n)
        ssRLeft = pairstoss(bpsRp, n)
        ssBRight = pairstoss(bpsB-bpsBp, n)
        ssRRight = pairstoss(bpsR-bpsRp, n)
        resL = realize(ssBLeft, ssRLeft, k, depth+1,
                       firstRblacklist, set([]), theta=theta)
        if resL is None:
            return None
        resR = realize(ssBRight, ssRRight, k-delta(bpsBp, bpsRp),
                       depth+1, set([]), lastBblacklist, theta=theta)
        if resL is not None and resR is not None:
            return resL + resR
        else:
            return None
    else:  # sep is None
        nextBlackList = firstRblacklist.copy()
        if DEBUG:
            print(indent(depth+1)+"Blacklists: first R= %s  last B = %s" %
                  (firstRblacklist, lastBblacklist))
        dependencies = builddependencies(bpsB, bpsR)
        if DEBUG:
            print(indent(depth+1)+" tries before filter %s" %
                  (len(dependencies)))
        dependencies = skipBlacklist(
            firstRblacklist, skipInclusions(dependencies))
        if DEBUG:
            print(indent(depth+1)+" tries after filter %s" %
                  (len(dependencies)))

        if DEBUG:
            print(indent(depth+1)+"Sep=0   -> Choosing first elem in R")
        for (p, deps) in dependencies:
            assert(len(deps) > 1)  # otherwise |B|<=k or sep is not None
            (i, j) = p
            (bpsnB, bpsnR) = (bpsB-deps, bpsR - set([(i, j)]))
            if DEBUG:
                print(indent(depth+1) +
                      "Try (%s,%s) in R -> Remove %s from B" % (i, j, deps))
            if len(deps) <= k:
                #print(indent(depth+1)+"nextBL = %s " % nextBlackList)
                recres = realize(pairstoss(bpsnB, n), pairstoss(
                    bpsnR, n), k+1-len(deps), depth+2, nextBlackList, lastBblacklist, theta=theta)
                if recres is not None:
                    return list(deps) + [(i, j)] + recres
                nextBlackList.add(p)
    return None


def printRevPath(bpsC, p, n):
    ss2 = pairstoss(bpsC, n)
    print(bpsC, ssstring(ss2))
    for (i, j) in reversed(p):
        if (i, j) in bpsC:
            bpsC = bpsC - set([(i, j)])
        else:
            bpsC = bpsC | set([(i, j)])
    ss1 = pairstoss(bpsC, n)
    printpath(p, ss1, ss2)


MIN_OF_MAX = 1


def bruteforce(C, B, R, n, P=[], curr_h=0, curr_max_h=0, h={}):
    # print len(C),len(B),len(R),len(P),curr_h,curr_max_h

    if MIN_OF_MAX not in h:
        h[MIN_OF_MAX] = +sys.maxsize

# if DEBUG:
##        ssC = pairstoss(C,n)
##        ssB = pairstoss(B,n)
##        ssR = pairstoss(R,n)
# print ssstring(ssC),ssstring(ssB),ssstring(ssR),P,curr_h,curr_max_h,h

    if curr_max_h > h[MIN_OF_MAX]:
        return (sys.maxsize, [])

    if len(R) == len(B) == 0:
        if curr_max_h < h[MIN_OF_MAX]:
            ssR = pairstoss(R, n)
            ssC = pairstoss(C, n)
            x = h[MIN_OF_MAX]
            if h[MIN_OF_MAX] == sys.maxsize:
                x = "+infty"
            else:
                x = "%s" % (x)
            print("[Update min h: %s -> %s]" % (x, curr_max_h))
            h[MIN_OF_MAX] = min(h[MIN_OF_MAX], curr_max_h)
        return (curr_max_h, P)

    # Free elements from R should be added as early as possible (no combinatorial exploration)
    for p1 in R:
        canBeAdded = True
        for p2 in C:
            if not compatiblepairs(p1, p2):
                canBeAdded = False
        if canBeAdded:
            addset = set([p1])
            return bruteforce(C | addset, B, R-addset, n, P+[p1], curr_h-1, curr_max_h, h)

    # If no elements left in R, remove pairs from B in any order
    if len(R) == 0:
        ncurr_h = curr_h+len(B)
        return bruteforce(C-B, set([]), set([]), n, P+list(B), ncurr_h, max(curr_max_h, ncurr_h), h)

    # An element is left in R, but not free, investigate all possibilities
    bestPath = (sys.maxsize, [])
    for (p1, deps) in builddependencies(C, R):
        # Destroy |deps| pairs in B (+|deps|), create one pair in B
        ncurr_h = curr_h+len(deps)-1
        # Max height reached after destroyung |deps| pairs in B
        local_max_h = curr_h+len(deps)
        res = bruteforce(C-deps | set([p1]), B-deps, R-set([p1]), n,
                         P+list(deps)+[p1], ncurr_h, max(curr_max_h, local_max_h), h)
        (height, p) = res
        bestPath = min(bestPath, (height, p))
    return bestPath

