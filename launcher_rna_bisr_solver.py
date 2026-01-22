
import argparse
from rna_barrier_bisr.rna_interface import solve_mcf
from rna_barrier_bisr.utilities import filter_common_bps
from rna_barrier_bisr.subtree_solver import solve_subtree


def compute_barrier_and_infos(s1full, s2full, rangealgo=0):
    if len(s1full)!= len(s2full) or len(["(" in s1full]) != len([")" in s1full]) or len(["(" in s2full]) != len([")" in s2full]):
        print("\ns1 or s2 is not a secondary structure or have a different size\n")
        return
    s1, s2 = filter_common_bps(s1full, s2full)
    if rangealgo:
        k, _ = solve_mcf(s1, s2)
    else:
        k, _, _ = solve_subtree(s1, s2)
    return k





parser = argparse.ArgumentParser()

parser.add_argument('--method', type=str, required=True, help='Either divide or subtree')
parser.add_argument('--s1', type=str, required=True, help='The starting secondary structure for the reconfiguration')
parser.add_argument('--s2', type=str, required=True, help='The ending secondary structure for the reconfiguration')
args = parser.parse_args()

if args.method == "divide":
    print("\nBarrier:", compute_barrier_and_infos(args.s1, args.s2, rangealgo=1), "\n")
elif args.method == "subtree":
    print("\nBarrier:", compute_barrier_and_infos(args.s1, args.s2, rangealgo=0), "\n")
else:
    print("\nProvided method is not supported\n")
