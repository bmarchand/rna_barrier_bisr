from rna_barrier_subtree_solver.subtree_solver import solve_subtree, list_bps, detect_sus, merge, barrier
import rna_barrier_subtree_solver.subtree_solver as subtree_solver
from rna_barrier_subtree_solver.utilities import filter_common_bps
from bisr_dpw.rna_interface import solve_mcf
from bisr_dpw.dpw_interface import solve_dpw
from bisr_dpw.barriers_py3 import ssrandom

seed = 2021
theta = 1
unpaired_weight = .5

def test_reinit_table():
    s1 = '(((....).)).....(((...)))...'
    s2 = '...((((.(..)))))...(((...)))'
    
    _, schedule = solve_subtree(s1, s2)

    assert(len(subtree_solver.M.keys()) > 0)
    subtree_solver.re_init_dp_table()
    assert(len(subtree_solver.M.keys()) == 0)

def test_detect_sus():

    s1 = '(((...).))....(((...)))...'
    s2 = '...(((.(..))))...(((...)))'
#         012345678910  14  18  22
#                    11  15  19  23
#                     12  16  20  24
#                      13  17  21  25

    schedule = [(0,9),(1,8),(7,10),(2,6),(3,13),(4,12),(5,11),(14,22),(15,21),(16,20),(17,25),(18,24),(19,23)]

    bps1 = list_bps(s1)
    bps2 = list_bps(s2)

    print(detect_sus(schedule, bps1, bps2)[0])
    assert(detect_sus(schedule, bps1, bps2)[0]==[(0,9),(1,8),(7,10),(2,6),(3,13),(4,12),(5,11)])
    assert(detect_sus(schedule, bps1, bps2)[1]==2)

def test_detect_sus2():
    
    s1 = '(((....).)).....(((...)))...'
    s2 = '...((((.(..)))))...(((...)))'

    bps1 = list_bps(s1)
    bps2 = list_bps(s2)
    
    _, schedule = solve_subtree(s1, s2)
    assert(len(detect_sus(schedule, bps1, bps2)[0])==8)
    assert(detect_sus(schedule, bps1, bps2)[1]==2)


def test_detect_sus3():

    s2 = '......................................................................................(((((......)))))..'
    s1 = '...............................................................................((((((......)))))).......'

    bps1 = list_bps(s1)
    bps2 = list_bps(s2)

    schedule = [(86, 101), (87, 100), (88, 99), (89, 98), (90, 97), (84, 91), (83, 92), (82, 93), (81, 94), (80, 95), (79, 96)]
    assert(len(detect_sus(schedule, bps2, bps1)[0])==11)
    assert(detect_sus(schedule, bps2, bps1)[1]==5)

def test_small_instances():

    for n in [10,20]:
        for inc in range(5):
            s1 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=inc)
            s2 = ssrandom(n, {}, theta=theta, unpaired_weight=unpaired_weight, seed=inc)
            s1, s2 = filter_common_bps(s1, s2)

            k1, schedule1 = solve_dpw(s1, s2)
            k2, schedule2 = solve_mcf(s1, s2)
            k3, schedule3 = solve_subtree(s1, s2)

            assert(k1==k2)
            assert(k1==k3)

def test_counter_example_merge():

    s1='.((.......((((((........(((((((((..............((()))))))))))))))(((((....(((((((......))))))))))))............)))((((((((((((..........)))))))).......))))))...'
    s2='(..(((((((......((((((((.........)))))))))))))).......................((((.......((((((............))))))))))))...............((((((((((........)))))))......)))'

    #s1 =".((......(((((.......((((((((............((()))))))))))))((((...((((((.....))))))))))..........))).((((((((((........))))))).....))))).."
    #s2 ="(..((((((.....(((((((........))))))))))))....................(((......(((((..........))))))))))...(..........)(((((((.......))))).....))"
    s1, s2 = filter_common_bps(s1, s2)
    assert(len(s1)==len(s2))
    k1, schedule1 = solve_dpw(s1, s2)
    k3, schedule3 = solve_subtree(s1, s2)
    assert(k3==9)
    assert(k3==k1)

def test_counter_example_merge2():

    s1 =".((......(((((.......((((((((............((()))))))))))))((((...((((((.....))))))))))..........))).((((((((((........))))))).....))))).."
    s2 ="(..((((((.....(((((((........))))))))))))....................(((......(((((..........))))))))))...(..........)(((((((.......))))).....))"
    s1, s2 = filter_common_bps(s1, s2)
    assert(len(s1)==len(s2))
    k1, schedule1 = solve_dpw(s1, s2)
    k3, schedule3 = solve_subtree(s1, s2)
    assert(k3==k1)
    assert(k3==9)

def test_counter_example_merge3():
    s1 = "(((......((((((((((((..........)))))))........))))))))..((...((((......))))))((((((((......)))))).....))"
    s2 = "...((((((............))))((((((.......))))))((........))..(((....)))))(..............)(((((......))))).."
    an = "EEEffccccCCCCCAAAAAAAccccaaaaaaAAAAAAAaaaaaaeeCCCCCEEEeeFFdddDDDDdddffgDDDDFFGGBBBBBBgbbbbbBBBBBBbbbbbGG"

    s1, s2 = filter_common_bps(s1, s2)
    assert(len(s1)==len(s2))
    k1, schedule1 = solve_dpw(s1, s2)
    k3, schedule3 = solve_subtree(s1, s2)
    print(schedule3)
    print(schedule1)
    assert(k3==k1)

def test_counter_example_merge4():
    s1 = "(((......((((((((((((..........)))))))........))))))))..((...((((.....))))))"
    s2 = "...((((((............))))((((((.......))))))((........))..(((....)))))......"
    an = "EEEffccccCCCCCAAAAAAAccccaaaaaaAAAAAAAaaaaaaeeCCCCCEEEeeFFdddDDDDdddffDDDDFF"

    s1, s2 = filter_common_bps(s1, s2)
    assert(len(s1)==len(s2))
    k1, schedule1 = solve_dpw(s1, s2)
    k3, schedule3 = solve_subtree(s1, s2)
    print(schedule3)
    for bp in schedule3:
        if bp in list_bps(s1):
            print('+1', end=' ')
        if bp in list_bps(s2):
            print('-1', end=' ')
    assert(k3==k1)

def test_counter_example_merge5():
    s1 = "(((....((((((((((((..........)))))))........)))))))).."
    s2 = "...((((............))))((((((.......))))))((........))"
    an = "EEEccccCCCCCAAAAAAAccccaaaaaaAAAAAAAaaaaaaeeCCCCCEEEee"

    s1, s2 = filter_common_bps(s1, s2)
    assert(len(s1)==len(s2))
    k1, schedule1 = solve_dpw(s1, s2)
    k3, schedule3 = solve_subtree(s1, s2)
    print(schedule3)
    assert(k3==k1)

def test_counter_example_merge6():
    s1 = "((...((((...))))))"
    s2 = "..(((....)))......"
    an = "FFdddDDDDdddDDDDFF"

    s1, s2 = filter_common_bps(s1, s2)
    assert(len(s1)==len(s2))
    k1, schedule1 = solve_dpw(s1, s2)
    k3, schedule3 = solve_subtree(s1, s2)
    print(schedule3)
    for bp in schedule3:
        if bp in list_bps(s1):
            print('+1', end=' ')
        if bp in list_bps(s2):
            print('-1', end=' ')
    assert(k3==k1)

def test_merge1():
    s1_out = "(((....((((((((((((..........)))))))........))))))))......................"
    s2_out = "...((((............))))((((((.......))))))((........))...................."
    an     = "EEEccccCCCCCAAAAAAAccccaaaaaaAAAAAAAaaaaaaeeCCCCCEEEee...................."
    _, schedule_out = solve_subtree(s1_out, s2_out)
    
    s1_in =  "........................................................((...((((...))))))"
    s2_in =  "..........................................................(((....)))......"
    an    =  "........................................................FFdddDDDDdddDDDDFF"
    _, schedule_in = solve_subtree(s1_in, s2_in)

    merge_sched = merge(schedule_in, schedule_out, list_bps(s1_in), list_bps(s2_in), list_bps(s1_out), list_bps(s2_out))
    
    s1 =     "(((....((((((((((((..........)))))))........))))))))....((...((((...))))))"
    s2 =     "...((((............))))((((((.......))))))((........))....(((....)))......" 
    print(s1)
    print(s2)
    assert(barrier(merge_sched, list_bps(s1), list_bps(s2))==7)

def test_merge2():

    s1_out = "(((......((((((((((((..........)))))))........))))))))..((...((((......))))))..........................."
    s2_out = "...((((((............))))((((((.......))))))((........))..(((....))))).................................."
    an     = "EEEffccccCCCCCAAAAAAAccccaaaaaaAAAAAAAaaaaaaeeCCCCCEEEeeFFdddDDDDdddffgDDDDFFGGBBBBBBgbbbbbBBBBBBbbbbbGG"
    _, schedule_out = solve_subtree(s1_out, s2_out)
    
    s1_in  = ".............................................................................((((((((......)))))).....))"
    s2_in  = "......................................................................................(((((......))))).."
    an     = "EEEffccccCCCCCAAAAAAAccccaaaaaaAAAAAAAaaaaaaeeCCCCCEEEeeFFdddDDDDdddffgDDDDFFGGBBBBBBgbbbbbBBBBBBbbbbbGG"
    _, schedule_in = solve_subtree(s1_in, s2_in)

    merge_sched = merge(schedule_in, schedule_out, list_bps(s1_in), list_bps(s2_in), list_bps(s1_out), list_bps(s2_out))
    
    s1 =     "(((......((((((((((((..........)))))))........))))))))..((...((((......))))))((((((((......)))))).....))"
    s2 =     "...((((((............))))((((((.......))))))((........))..(((....)))))................(((((......))))).."
    print(s1)
    print(s2)
    assert(barrier(merge_sched, list_bps(s1), list_bps(s2))==7)

def test_reverse_sus1():
    s1_in =  "........................................................((...((((...))))))"
    s2_in =  "..........................................................(((....)))......"
    an    =  "........................................................FFdddDDDDdddDDDDFF"
    _, schedule_in = solve_subtree(s1_in, s2_in)

    assert(len(detect_sus(schedule_in[::-1], list_bps(s2_in), list_bps(s1_in))[0])==2)

def test_small_instance5():
    
    s1 ="...(..)..."
    s2 =".........."
    s1, s2 = filter_common_bps(s1, s2)

    k1, schedule1 = solve_dpw(s1, s2)
    k2, schedule2 = solve_mcf(s1, s2)
    k3, schedule3 = solve_subtree(s1, s2)

    assert(k1==k2)
    assert k1==k3, "dpw said "+str(k1)+" while subtree said "+str(k3)

def test_small_instance6():
    
    s1 =".((...(((.)))))"
    s2 ="((((.)((.).))))"
    s1, s2 = filter_common_bps(s1, s2)

    k1, schedule1 = solve_dpw(s1, s2)
    k2, schedule2 = solve_mcf(s1, s2)
    k3, schedule3 = solve_subtree(s1, s2)

    assert(k1==k2)
    assert(k1==k3)

def test_small_instance7():
    
    s1 = '.(.(...)).'
    s2 = '(..(...).)'
    s1, s2 = filter_common_bps(s1, s2)

    k1, schedule1 = solve_dpw(s1, s2)
    k2, schedule2 = solve_mcf(s1, s2)
    k3, schedule3 = solve_subtree(s1, s2)

    assert(k1==k2)
    assert(k1==k3)


def test_small_instance8():
    
    s1= "........()"
    s2= "(...())..."
    s1, s2 = filter_common_bps(s1, s2)

    k1, schedule1 = solve_dpw(s1, s2)
    k2, schedule2 = solve_mcf(s1, s2)
    k3, schedule3 = solve_subtree(s1, s2)

    assert(k1==k2)
    assert(k1==k3)

def test_small_instance9():
    
    s1="(.((.......................(.((..(.(...))....)))..))......................)..."
    s2=".(..(((((((((((((((((((((((.(..((.(.(((..))))...))..)))))))))))))))))))))).)))"
    s1, s2 = filter_common_bps(s1, s2)

    k1, schedule1 = solve_dpw(s1, s2)
    k2, schedule2 = solve_mcf(s1, s2)
    k3, schedule3 = solve_subtree(s1, s2)

    assert(k1==k2)
    assert(k1==k3)

def test_small_instance11():
    s1 = '..((.))(.)'
    s2 = '(.(((.))))'

    s1, s2 = filter_common_bps(s1, s2)

    k1, schedule1 = solve_dpw(s1, s2)
    k2, schedule2 = solve_mcf(s1, s2)
    k3, schedule3 = solve_subtree(s1, s2)

    assert(k1==k2)
    assert(k1==k3)

if __name__=='__main__':
    test_merge2()
