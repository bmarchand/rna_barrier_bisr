from .utilities import num_leaves, list_bps

debug = True

def solve_subtree(s1, s2):

    # TODO ADD REMOVING OF COMMON BP + PRINT WARNING

    # reinit DP table:
    global M
    M = {}
        
    # switching to list of bps representation
    l1 = list_bps(s1)
    l2 = list_bps(s2)

    # length of rna string: always same thing
    global N
    N = len(s1)
    assert(len(s1)==len(s2))

    if num_leaves(s2) < num_leaves(s1):
        k, schedule = DP_solver(l2, l1)
        return k-len(l2)+len(l1), schedule[::-1], num_leaves(s2)

    k, sched = DP_solver(l1, l2)

    return k, sched, num_leaves(s1) 

M = {}

def re_init_dp_table():
    global M
    M = {}

def key(l1,l2):
    s1 = ['.' for _ in range(N)]
    s2 = ['.' for _ in range(N)]
    for i,j in l1:
        s1[i] = '('
        s1[j] = ')'
    for i,j in l2:
        s2[i] = '('
        s2[j] = ')'
    return '.'.join(s1)+'->'+''.join(s2)

def conflict(bp1, bp2):
    i1, j1 = bp1
    i2, j2 = bp2
    return (i1 <= i2 and i2 <= j1 and j2 >= j1) or (i2 <= i1 and i1 <= j2 and j1 >= j2)

def neighbor_dict(bps1, bps2):
    
    N = {}

    # init keys
    for bp1 in bps1:
        N[bp1] = []
    for bp2 in bps2:
        N[bp2] = []

    for bp1 in bps1:
        for bp2 in bps2:
            if conflict(bp1, bp2):
                N[bp1].append(bp2)
                N[bp2].append(bp1)

    return N

def delete_bp(bps, bp):
    """
    should return a copy of bps in which bp has been removed, not
    modify bps in place
    """
    return [x for x in bps if x!=bp]

def barrier(schedule, bps1, bps2):
    """
    computes the barrier of schedule, which
    transforms s1 into s2
    """
    # a silly check
    if debug:
        assert len(schedule)==len(bps1)+len(bps2), "len schedule is "+str(len(schedule))+" and len(1) and len(2):"+str(len(bps1))+" and "+str(len(bps2))

    worst = 0
    cur = 0 # ``current'' energy value

    for bp in schedule:
        if bp in bps1:
            cur += 1
        if bp in bps2:
            cur -= 1
        if cur > worst:
            worst = cur

    return worst

def split(bps, bp):
    """
    splits s into s_in, s_out
    w.r.t bp.
    base pairs from s in conflict with
    bp are ignored
    """
    i,j = bp

    bps_in = []
    bps_out = []

    for k,l in bps:
        if (k < i and l > j) or (k < i and l < i) or (k > j and l > j):
            bps_out.append((k,l))
        if k > i and k < j and l > i and l < j:
            bps_in.append((k,l))

    return bps_in, bps_out 

def detect_sus(schedule, bps1, bps2):
    """
    It is assumed that schedule is in canonical 
    form. If an SUS exsists, schedule starts with
    an SUS, which can simply be detected by looking
    at the first point of the schedule where cost < 0

    Important: this function may return the entire schedule,
    if it does not contain any SUS
    """

    cost = 0
    worst = 0
    for k, bp in enumerate(schedule):
        if debug:
            assert(bp in bps1 or bp in bps2)
        # if immediate gain: greedy and recurse
        if bp in bps2:
            cost -= 1
            continue
        # else do I stop now or look at base pair ?
        if cost < 0 and k > 0:
            return schedule[:k], worst
        if bp in bps1:
            cost += 1
        if cost > worst:
            worst = cost

    # complete instance case
    if cost < 0:
        return schedule, worst

    return None, float('inf')

def filter_bps(bps, sub_sched):
    """
    s is a structure and sub_sched a subschedule, usually
    a SUS
    """
    for bp in sub_sched:
        if bp in bps:
            bps = delete_bp(bps, bp)

    return bps

def merge(schedule_in, schedule_out, bps1_in, bps2_in, bps1_out, bps2_out):
    """
    The crux of the algorithm.
    schedule1 and schedule2 are assumed to be in canonical form.
    Their SUS and reverse-SUSs are detected and put in front/rear
    """
    if debug:
        assert(len(schedule_in)==len(bps1_in)+len(bps2_in))
        assert(len(schedule_out)==len(bps1_out)+len(bps2_out))
    
    if len(schedule_in)==0:
        return schedule_out
    if len(schedule_out)==0:
        return schedule_in

    # detecting SUS in schedule_in
    sus_in, worst_in = detect_sus(schedule_in, bps1_in, bps2_in) 
    
    # detecting SUS in schedule_out
    sus_out, worst_out = detect_sus(schedule_out, bps1_out, bps2_out) 

    # putting it first
    if sus_in or sus_out:
        if worst_out < worst_in:
            schedule_in, schedule_out = schedule_out, schedule_in
            bps1_in, bps1_out = bps1_out, bps1_in
            bps2_in, bps2_out = bps2_out, bps2_in
            sus_in, sus_out = sus_out, sus_in
            
        bps1_in = filter_bps(bps1_in, sus_in)
        bps2_in = filter_bps(bps2_in, sus_in)
        
        schedule = sus_in+merge(schedule_in[len(sus_in):], schedule_out, bps1_in, bps2_in, bps1_out, bps2_out)
        if debug:
            assert(len(schedule)==len(schedule_in)+len(schedule_out))
        return schedule
        #return sus_in+merge(schedule_in[len(sus_in):], schedule_out, s1_in, s2_in, s1_out, s2_out) 

    # detecting reverse SUS in schedule_in
    rsus_in, worst_in = detect_sus(schedule_in[::-1], bps2_in, bps1_in) 

    # detecting SUS in schedule_out
    rsus_out, worst_out = detect_sus(schedule_out[::-1], bps2_out, bps1_out) 

    # putting it first
    if rsus_in or rsus_out:
        if worst_out < worst_in:
            schedule_in, schedule_out = schedule_out, schedule_in
            bps1_in, bps1_out = bps1_out, bps1_in
            bps2_in, bps2_out = bps2_out, bps2_in
            rsus_in, sus_out = rsus_out, sus_in
            
        rsus_in.reverse()
        bps1_in = filter_bps(bps1_in, rsus_in)
        bps2_in = filter_bps(bps2_in, rsus_in)
        schedule = merge(schedule_in[:-len(rsus_in)], schedule_out, bps1_in, bps2_in, bps1_out, bps2_out)+rsus_in
        if debug:
            assert(len(schedule)==len(schedule_in)+len(schedule_out))
        return schedule
        #return merge(schedule_in[:-len(rsus_in)], schedule_out, s1_in, s2_in, s1_out, s2_out)+rsus_in
    
    # then, no SUS in any direction: delta=0 and can concateneate:
    return schedule_in + schedule_out

def budget_vector(schedule, bps1, bps2):
    """
    Returns a vector (list) containing bg(schedule)
    pref[i] and suff[i] for all i
    """

    pref = {}  # storing pref_i, the budget to go to level i <= 0
    suff = {}  # storing suff_i, the budget to go to level i <= 0 from end

    # initializing these values at +infty
    for i in range(0,-len(bps1)-len(bps2),-1):
        pref[i] = float('inf')
        suff[i] = float('inf')

    # looping over schedule1 and tracking energy
    cost = 0
    worst = 0

    # pref dictionary
    for bp in schedule:
        if bp in bps1:
            cost += 1
        if bp in bps2:
            cost -= 1
        if cost <= 0 and not pref[cost] < float('inf'):
            pref[cost] = worst
        if cost > worst:
            worst = cost

    cost = 0
    worst_rev = 0
    
    # suff dictionary
    for bp in schedule[::-1]:
        if bp in bps2:
            cost += 1
        if bp in bps1:
            cost -= 1
        if cost <= 0 and not suff[cost] < float('inf'):
            suff[cost] = worst_rev
        if cost > worst_rev:
            worst_rev = cost

    return [worst] + [val for _, val in pref.items()] + [val for _, val in suff.items()]

def preferable(schedule1, schedule2, bps1, bps2):
    """
    is schedule1 preferable to schedule2 for going
    from s1 to s2, return True. False if less preferable
    or not comparable
    """
    budget_vector1 = budget_vector(schedule1, bps1, bps2)
    budget_vector2 = budget_vector(schedule2, bps1, bps2)
    if debug:
        assert(len(budget_vector1)==len(budget_vector2))
    
    for k in range(len(budget_vector1)):
        if budget_vector1[k] > budget_vector2[k]:
            return False

    return True

def DP_solver(bps1, bps2, depth=0):

    # checking if already computed
    if key(bps1, bps2) in M:
        return M[key(bps1,bps2)]
        
    # adjacency dictionary computation: (possible optim: do once outside and subsample)
    N = neighbor_dict(bps1, bps2)
    
    ### TERMINAL CASES ###

    # no-dependency final element
    for bp2 in bps2:
        if len(N[bp2])==0:
            # solving and re-adding red
            k, schedule = DP_solver(bps1, delete_bp(bps2, bp2), depth=depth+1)
            k = max(0, k-1)
            schedule = [bp2] + schedule

            # storing and returning
            if debug:
                assert(len(bps1)+len(bps2)==len(schedule))
            M[key(bps1, bps2)] = k, schedule
            return M[key(bps1,bps2)]

    # no-reward initial element
    for bp1 in bps1:
        if len(N[bp1])==0:
            # solving and re-adding blue
            N = len(bps1)
            assert(len(delete_bp(bps1, bp1))==len(bps1)-1)
            assert(N==len(bps1))
            k, schedule = DP_solver(delete_bp(bps1, bp1), bps2, depth=depth+1)
            schedule = schedule + [bp1]
            k = barrier(schedule, bps1, bps2)

            # storing and returning
            if debug:
                assert(len(bps1)+len(bps2)==len(schedule))
            M[key(bps1, bps2)] = k, schedule
            return M[key(bps1,bps2)]

    ### GENERAL CASE ###
    # dummy schedule to start:
    M[key(bps1,bps2)] = len(bps1), bps1+bps2

    for bp1 in bps1[::-1]:
        # is bp1 the last ?
        bps1_in, bps1_out = split(bps1, bp1)
        bps2_in, bps2_out = split(bps2, bp1)

        _, schedule_in = DP_solver(bps1_in, bps2_in, depth=depth+1)
        _, schedule_out = DP_solver(bps1_out, bps2_out, depth=depth+1)
        assert(len(schedule_in)==len(bps1_in)+len(bps2_in))
        assert(len(schedule_out)==len(bps1_out)+len(bps2_out))

        schedule = merge(schedule_in, schedule_out, bps1_in, bps2_in, bps1_out, bps2_out)
        schedule = schedule + [bp1] + N[bp1]
        if debug:
            assert(len(bps1)+len(bps2)==len(schedule))
        if preferable(schedule, M[key(bps1, bps2)][1], bps1, bps2):
            M[key(bps1, bps2)] = barrier(schedule, bps1, bps2), schedule

    return M[key(bps1, bps2)]
