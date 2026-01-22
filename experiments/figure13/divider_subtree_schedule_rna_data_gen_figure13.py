import time
import csv
import argparse
import os
from rna_barrier_bisr.rna_interface import solve_mcf, pre_processing
from rna_barrier_bisr.dpw_interface import secondary_structures_to_G
from rna_barrier_bisr.utilities import filter_common_bps
from networkx.algorithms.bipartite.matching import maximum_matching
from rna_barrier_bisr.general_bipartite import vertex_cover_from_matching
from rna_barrier_bisr.random_rna_structures import ssparse
from rna_barrier_bisr.subtree_solver import solve_subtree, DP_solver
from rna_barrier_bisr.utilities import num_leaves, list_bps
tab = []

def compute_barrier_and_infos(name, s1full, s2full, norangealgo=0, noarboricityalgo=0):
    if len(s1full)!= len(s2full) or len(["(" in s1full]) != len([")" in s1full]) or len(["(" in s2full]) != len([")" in s2full]):
        print("Input Issue\n")
        return
    s1, s2 = filter_common_bps(s1full, s2full)
    sizeL = s1.count("(")
    sizeR = s2.count(")")
    B,R,G = secondary_structures_to_G(ssparse(s1),ssparse(s2))
    M = maximum_matching(G, top_nodes=B)
    K = vertex_cover_from_matching(B,R,G,M)
    I = set(G) - K


    alpha=len(I)
    if not norangealgo:
        t0 = time.time()
        k, _ = solve_mcf(s1, s2)
        time_range_algo = time.time() - t0
        rho = k+alpha-sizeL
        barrier = k
    else:
        time_range_algo = -1
        rho = -1
        barrier = -1
    if not noarboricityalgo:
        t0 = time.time()
        phi1, phi2 = (s1.replace('.','').count('()'), s2.replace('.','').count('()'))
        k, _, _ = solve_subtree(s1, s2)
        time_arboricity_algo = time.time() - t0
    else:
        time_arboricity_algo = -1
        phi1, phi2 = -1, -1
        k = -1
    if k != barrier and barrier !=-1 and k !=-1:
        print("Algorithm Issue\n",[s1full, s2full, len(s1full), len([l for l in s1full if l == "("]), len([l for l in s2full if l == "("]), s1, s2, barrier, rho, phi1, phi2, min(phi1, phi2), time_range_algo, time_arboricity_algo])
        return
    elif barrier == -1:
        barrier = k
    return [name, s1full, s2full, len(s1full), len([l for l in s1full if l == "("]), len([l for l in s2full if l == "("]), s1, s2, barrier, rho, phi1, phi2, min(phi1, phi2), time_range_algo, time_arboricity_algo]


fname = os.getcwd() + '/experiments/NewTest2/random_structures_data_set_figure13.csv'
N = len(open(fname).readlines())
MAX_LENGTH = 101

def build_csv_random_rna_instances(restart=1, norangealgo=0, last_index=-1):
    c = 'a'
    if restart:
        c = 'w'
    with open('ResultsforRNARandomBenchmark.csv', c, newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        if restart:
            elem = ["name", "s1", "s2", "size_si", "nbBPs1", "nbBPs2", "s1nosameBPs", "s2nosameBPs", "barrier", "range", "arboricityL","arboricityR","arboricity", "time_range_algo(secs)", "time_arboricity_algo(secs)"]
            writer.writerow(elem)
        #length_seen = -1
        for line in open(fname).readlines()[1:]:
            index = line.split(',')[0]
            length = int(line.split(',')[1])
            if int(index) <= last_index:
                continue
            s1full = line.split(',')[2] 
            s2full = line.split(',')[3] 
            print("Now working on ", index, 'out of',N, "of size", length, "\n")
            print(s1full,'->',s2full)
            elem = compute_barrier_and_infos(index, s1full, s2full, norangealgo=norangealgo)
            writer.writerow(elem)

def build_csv_random_rna_instances_full(restart=1,last_index=-1):
    c = 'a'
    if restart:
        c = 'w'
    with open('ResultsforRNARandomBenchmark.csv', c, newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        if restart:
            elem = ["name", "s1", "s2", "size_si", "nbBPs1", "nbBPs2", "s1nosameBPs", "s2nosameBPs", "barrier", "range", "arboricityL","arboricityR","arboricity", "time_range_algo(secs)", "time_arboricity_algo(secs)"]
            writer.writerow(elem)
        for line in open(fname).readlines()[1:]:
            index = line.split(',')[0]
            length = int(line.split(',')[1])
            number = int(index)%200

            if int(index) <= last_index or number >=100:
                continue
            noarboricityalgo=0
            norangealgo=0
            s1full = line.split(',')[2] 
            s2full = line.split(',')[3] 
            print("Now working on ", index, 'out of',N, "of size", length, "\n")
            print(s1full,'->',s2full)
            
            elem = compute_barrier_and_infos(index, s1full, s2full, norangealgo=norangealgo, noarboricityalgo=noarboricityalgo)
            writer.writerow(elem)



parser = argparse.ArgumentParser()

parser.add_argument('--restart', type=str, required=False)
parser.add_argument('--lastindex', type=str, required=False)
args = parser.parse_args()


restart = 0
if args.restart:
    restart = int(args.restart)
last_index = -1
if args.lastindex:
    last_index = int(args.lastindex)

build_csv_random_rna_instances_full(restart=restart, last_index=last_index)