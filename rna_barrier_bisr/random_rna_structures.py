import random
import sys
import os
import getopt
import time

## Parameters for random generation ##
MULTIPLE_HELIX_LENGTH = 1

DEBUG = False
SMALLDEBUG = True
LEAVEINCLUDES = True
LEAVEBLACKLIST = True
RUNBRUTEFORCE = False
TIME = False

maxk = 0
############# Basic utilities #################


def boolTuples(dim):
    if dim <= 0:
        return [()]
    else:
        res = []
        for t in boolTuples(dim-1):
            for b in [True, False]:
                res.append(tuple([b]+list(t)))
        return res


############# RNA structure stuffs #################
# Mostly "house-keeping" functions (parsing, conversions,
# uniform random generation...)

def ssparse(seq):
    """
    from well-parenthesized to table with -1 for unpaired and ( or )
    otherwise
    """
    res = [-1 for c in seq]
    p = []
    for i, c in enumerate(seq):
        if c == '(':
            p.append(i)
        elif c == ')':
            j = p.pop()
            (res[i], res[j]) = (j, i)
    return res


def ssstring(ss):
    """
    from table with -1 on unpaired position and ( ) elsewhere
    to well-parenthesized string
    """
    res = ["." for i in ss]
    for i in range(len(ss)):
        if ss[i] == -1:
            pass
        elif ss[i] > i:
            res[i] = '('
        elif ss[i] < i:
            res[i] = ')'
    return "".join(res)


def sstopairs(ss):
    """
    well-parenthesized string to list of base pairs
    """
    if ss is None:
        return None
    res = []
    for i in range(len(ss)):
        if ss[i] > i:
            j = ss[i]
            res.append((i, j))
    return set(res)


def pairstoss(bps, n):
    res = [-1 for i in range(n)]
    for (i, j) in bps:
        res[i], res[j] = j, i
    return res


def sscount(n, count, theta=3, unpaired_weight=0.1):
    if n not in count:
        if n <= 0:
            count[n] = 1.
        else:
            count[n] = unpaired_weight*sscount(n-1, count, theta=theta)
            for i in range(theta+2, n+1):
                count[n] += sscount(i-2, count, theta=theta) * \
                    sscount(n-i, count, theta=theta)
    return count[n]


# Taille 50

BASE_PAIRS = set([("A", "U"), ("U", "A"), ("G", "U"),
                  ("U", "G"), ("C", "G"), ("G", "C")])


def fillmatrix(seq, theta, taboo):
    n = len(seq)
    tab = [[-10000 for i in range(n)] for j in range(n)]
    for m in range(1, theta+2):
        for i in range(0, n-m+1):
            j = i + m-1
            tab[i][j] = 0
    for m in range(theta+2, n+1):
        for i in range(0, n-m+1):
            j = i + m-1
            tab[i][j] = tab[i+1][j]
            for k in range(i+theta+1, j+1):
                if ((i, k) not in taboo) and (seq[i], seq[k]) in BASE_PAIRS:
                    t2 = 0
                    if (k < j):
                        t2 = tab[k+1][j]
                    tab[i][j] = max(tab[i][j], 1+tab[i+1][k-1]+t2)
    return tab


def backtrack(seq, i, j, tab, theta, taboo):
    n = len(seq)
    m = j-i+1
    if m < theta+2:
        return []
    else:
        if tab[i][j] == tab[i+1][j]:
            return backtrack(seq, i+1, j, tab, theta, taboo)
        for k in range(i+theta+1, j+1):
            if ((i, k) not in taboo) and (seq[i], seq[k]) in BASE_PAIRS:
                t2 = 0
                if (k < j):
                    t2 = tab[k+1][j]
                if tab[i][j] == 1+tab[i+1][k-1]+t2:
                    return [(i, k)] + backtrack(seq, i+1, k-1, tab, theta, taboo) + backtrack(seq, k+1, j, tab, theta, taboo)
    raise Exception


def randomUniformSeq(n, seed=None):
    if seed:
        random.seed(seed)
    return "".join([random.choice(['A', 'C', 'G', 'U']) for i in range(n)])


def randomDistantStructPairs(n, theta=3, seed=None):
    rseq = randomUniformSeq(n, seed=seed)
    taboo = set()
    tab1 = fillmatrix(rseq, theta, taboo)
    bps1 = backtrack(rseq, 0, n-1, tab1, theta, taboo)
    taboo = set(bps1)
    tab2 = fillmatrix(rseq, theta, taboo)
    bps2 = backtrack(rseq, 0, n-1, tab2, theta, taboo)
    return (rseq, pairstoss(bps1, n), pairstoss(bps2, n))


def ssrandom(n, count, seed=None, theta=1, unpaired_weight=0.1):
    if seed:
        random.seed(seed)
    if n <= 0:
        return ''
    else:
        r = random.random()*sscount(n, count, theta=theta, unpaired_weight=unpaired_weight)
        r -= unpaired_weight * \
            sscount(n-1, count, theta=theta, unpaired_weight=unpaired_weight)
        if r < 0:
            return '.'+ssrandom(n-1, count, theta=theta, unpaired_weight=unpaired_weight)
        for i in range(theta+2, n+1):
            r -= sscount(i-2, count, theta=theta, unpaired_weight=unpaired_weight) * \
                sscount(n-i, count, theta=theta,
                        unpaired_weight=unpaired_weight)
            if r < 0:
                return '('+str(ssrandom(i-2, count, theta=theta, unpaired_weight=unpaired_weight))+')'+str(ssrandom(n-i, count, theta=theta, unpaired_weight=unpaired_weight))


def compatiblepairs(p1, p2):
    (i, j) = p1
    (k, l) = p2
    return (i < j < k < l) or (i < k < l < j) or (k < l < i < j) or (k < i < j < l)


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


def delta(B, R):
    return len(B)-len(R)


#if __name__ == "__main__":
#    # Set random seed for sake of reproducibility
#    seed = 0
#    n = 50
#    maxAllowedK = -1
#    try:
#        opts, args = getopt.getopt(sys.argv[1:], "n:s:idtfbk:")
#    except getopt.GetoptError:
#        print('-n <size> -s <seed> -i (leave includes) -d (debug) -t (run timer, silent mode) -f (brute force) -b (leave blacklists) -k <max allowed k>')
#        sys.exit(2)
#    for opt, arg in opts:
#        if opt == '-n':
#            n = int(arg)
#        elif opt == '-s':
#            seed = int(arg)
#        elif opt == '-k':
#            maxAllowedK = int(arg)
#        elif opt == '-i':
#            LEAVEINCLUDES = True
#        elif opt == '-b':
#            LEAVEBLACKLIST = True
#        elif opt == '-f':
#            RUNBRUTEFORCE = True
#        elif opt == '-d':
#            DEBUG = True
#        elif opt == '-t':
#            TIME = True
#    # print csv format: n, seed, options, k,maxk,<empty>,<empty>, time
#    sys.stdout.write("%s,%s,v3b_%s," % (
#        n, seed, ("skipIncludes" if not LEAVEINCLUDES else "base")+("_skipBL" if not LEAVEBLACKLIST else "")))
#    if TIME:
#        saveStdout = sys.stdout
#        f = open(os.devnull, 'w')
#        sys.stdout = f
#        start_time = time.time()
#
#    random.seed(seed)
#    count = {}
#
#    s1 = ssrandom(n, count)
#    s2 = ssrandom(n, count)
#    if s2.count("(") > s1.count("("):
#        s1, s2 = s2, s1
#
#    ss1 = ssparse(s1)
#    ss2 = ssparse(s2)
#    bps1 = sstopairs(ss1)
#    bps2 = sstopairs(ss2)
#
#    # Remove common base pairs
#    commonbps = bps1 & bps2
#    if len(commonbps) > 0:
#        for (i, j) in commonbps:
#            (ss1[i], ss1[j]) = (-1, -1)
#            (ss2[i], ss2[j]) = (-1, -1)
#        s1 = ssstring(ss1)
#        s2 = ssstring(ss2)
#        bps1 = sstopairs(ss1)
#        bps2 = sstopairs(ss2)
#
#    print(s1, len(sstopairs(ss1)))
#    print(s2, len(sstopairs(ss2)))
#
#    if maxAllowedK < 0:
#        maxAllowedK = n
#    print("[Running XP algo]")
#    for k in range(0, maxAllowedK+1):
#        print("  [k=%s]" % k)
#        p = realize(ss1, ss2, k, 1)
#        if p is not None:
#            print(s1, "Origin")
#            printpath(p, ss1, ss2)
#            print(s2.replace("(", "[").replace(")", "]"), "Destination")
#            break
#        else:
#            print("    No solution found")
#            if k == maxAllowedK:
#                k = k+0.5
#
#    if TIME:
#        sys.stdout = saveStdout
#        print("%s,%s,,,%s" % (k, maxk, time.time() - start_time))
#
#    if RUNBRUTEFORCE:
#        print("[Running branch and bound algo]")
#        (minh, p) = bruteforce(bps1, bps1, bps2, n)
#        print(s1, "Origin")
#        printpath(p, ss1, ss2)
#        print(s2.replace("(", "[").replace(")", "]"), "Destination")
#        # sys.exit()


# for i in {0..15}; do echo $i; python ./barriers.py -n 50 -t -s $i>>times;   python ./barriers.py -n 50 -ti -s $i >>times; done
