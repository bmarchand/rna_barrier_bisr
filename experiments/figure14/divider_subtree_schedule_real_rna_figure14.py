import time
import csv
from rna_barrier_bisr.rna_interface import solve_mcf, pre_processing
from rna_barrier_bisr.dpw_interface import secondary_structures_to_G
from rna_barrier_bisr.utilities import filter_common_bps
from networkx.algorithms.bipartite.matching import maximum_matching
from rna_barrier_bisr.general_bipartite import vertex_cover_from_matching
from rna_barrier_bisr.random_rna_structures import ssparse
from rna_barrier_bisr.subtree_solver import solve_subtree

tab = []
tab.append(
    ("adenine_riboswitch_of_add_gene_V.vulnificus",
     "............((((((.........))))))........((((((.......)))))).((((((.((((.((.((((((..........)))))).))))))))))))..",
     "(((((((((...((((((.........))))))........((((((.......))))))..))))))))).........................................."
    )    
)

tab.append(
    ("lysine_riboswitch_of_lysC_gene",
     "..((((((..(.((((.((((........((((((...(((.......)))......)))))).......))))..)))))..(((((((((.(((.....))).)))))))))((((.((((((.....)))))).)))).((((((((....))))))))....))))))......((((((((((.....)))))))..)))...((((((((((.....))))))))))",
     ".....(((..(.((((.((((........((((((...(((.......)))......)))))).......))))..)))))..(((((((((.(((.....))).)))))))))((((.((((((.....)))))).)))).((((((((....))))))))....)))(((((((((((((((((((.....)))))))..))...))))))))))................"
    )    
)

#tab.append(
#    ("5prime_leader_of_MS2",
#     "((((((((((((....))))))))))))(((((((((((..(((((...((((...))))...)))))..)))(((((((((((....))))))))))).(((((...))))).))))))))............",
#     "XXXXXXXXXXXXXXXXXXXX...((((((((((((..((((.(((..((((....))))..))).))))..))))))((((((.(((....))).).)))))(((((((....)))))))))))))........"
#    )    
#)

tab.append(
    ("TPP_riboswitch_of_thiamine_gene",
     "...(((((((.(.....).))).))))..(((((.(((((((((((...)))))).....((((((.....)))))).))))).....((((..(((((....)))))..))))..)))))...(((((..(((......))).)))))......(((((((((((((....)))))))))))))",
     "...........................(((((((.(((((((((((...)))))).....((((((.....)))))).))))).......((((....))))..((((((..(..((((((...((.(((.((((.........)))).))).))...))))))..)..))))))..)))))))."
    )    
)

tab.append(
    ("adenine_riboswitch_of_ydhL_gene",
     "................((((((.........))))))..................((((((((((((((((((..((((...))))..))))))))))))))))))....",
     "........(((((...((((((.........))))))........(((((.........)))))..)))))....((((...))))........................"
    )    
)

tab.append(
    ("cyclic-di-GMP_riboswitch_of_tfoX_gene",
     "((((((.....((...((((.((....))))))...))...(((.((((((......((((....))))...)))..))))))...))))))................((((((....))))))",
     "...........((...((((.((....))))))...))..............................(((((.(((((...((((.(((((....))))).)))))))))).))))......."
    )    
)

tab.append(
    ("SAM_riboswitch_of_metE_gene",
     "((((((((....(.(((...(((.....)))......))))(((..((((((....))).))).)))........(((((......)))))...)))))))).........((((((((.......))))))))",
     "............(.(((...(((.....)))......))))(((..((((((....))).))).)))........(((((......)))))(((((.((((...........)))).)))))............"
    )    
)

#tab.append(
#    ("SV-11_ground_state",
#     "(((((((((((...)))))))..((((((((((....))))))))))........)))).....((((((((.(....).)))))))).((((((((.....)))))))).....",
#     "(((...((((((.....))))))((.((((((((((....))))))))))))......)))......((((((((((....))))))))))..(((((((...)))))))....."
#    )    
#)

tab.append(
    ("SV-11_metastable_plus_to_stable_plus",
     "(((((((((((...)))))))..((((((((((....))))))))))........)))).....((((((((.(....).)))))))).((((((((.....)))))))).....",
     "(((.(((((((.(((((((((((((((((((...((((((((((((((..(((.....)))..))))))))))))))...)))))))))))))))))))...)))))))..)))."
    )    
)

tab.append(
    ("guanine_riboswitch_of_xptpubX_gene",
     "((((((((...(.(((((.......))))).)........((((((.......))))))..))))))))........(((((........)))))............((((((((((((((.......))))))))))))))......",
     ".....(((...(.(((((.......))))).)........((((((.......))))))..)))((((((((((((.(((((........)))))..............))))))))))))..........................."
    )    
)

Active =     "......................................................(((((((.........(((........))))))))))......((((((((((((.....)))))))).)))).............."
Inhibitory = "..............................((((((((((........(((...(((((((......(((((......))))).))))))))))...((((((((((((.....)))))))).))))..)).))))))))." 
Permissive = "((((((.((.(((((.........))))))))))))).................(((((((.........(((........))))))))))......((((((((((((.....)))))))).)))).............."
#tab.append(
#    ("HDV_ribozyme_active_to_inhibitory", #Active C Permissive
#     Active,
#     Inhibitory
#    )    
#)

#tab.append(
#    ("HDV_ribozyme_active_to_permissive", #Active C Permissive
#     Active,
#     Permissive
#    )    
#)


tab.append(
    ("HDV_ribozyme_inhibitory_to_permissive", #Pseudoknot ignored, non specified pairs considered as unpaired.
     Inhibitory,
     Permissive
    )    
)

tab.append(
    ("5prime_UTR_of_the_Leviviridae_Levivirus",
     "....((((((((....))))))))((((.((.....)).)))).........................................................................................",
     "((((((((((((....))))))))))))(((((((((((..(((((...((((...))))...)))))..)))(((((((((((....))))).)))))).(((((...))))).))))))))........."
    )    
)


tab.append(
    ("tryphtophan_operon_leader",
     "....((((((((.(((((.(((.....))).)))))..)))))).))...................((((((((....).))))))).......", 
     "..........................(((((.(.((((......................)))).).)))))......................"
    )    
)

tab.append(
    ("ZTP_riboswitch_IH1_to_P1",
     ".....(((((..)))))....((((.....))))....................", #non specified pairs considered as unpaired.
     ".(((((((((...........((((.....))))......)))))))))....."
    )    
)

def compute_barrier_and_infos(name, s1full, s2full):
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
    t0 = time.time()
    k, _ = solve_mcf(s1, s2)
    time_range_algo = time.time() - t0
    rho = k+alpha-sizeL
    barrier = k#+alpha-sizeL
    t0 = time.time()
    phi1, phi2 = (s1.replace('.','').count('()'), s2.replace('.','').count('()'))
    k, _, _ = solve_subtree(s1, s2)
    time_arboricity_algo = time.time() - t0
    if k != barrier:
        #if k - sizeL + sizeR != barrier :
        print(k, barrier, sizeL, sizeR)
        print("Algorithm Issue\n",[s1full, s2full, len(s1full), len([l for l in s1full if l == "("]), len([l for l in s2full if l == "("]), s1, s2, barrier, rho, phi1, phi2, min(phi1, phi2), time_range_algo, time_arboricity_algo])
        return
    return [name, s1full, s2full, len(s1full), len([l for l in s1full if l == "("]), len([l for l in s2full if l == "("]), s1, s2, barrier, rho, phi1, phi2, min(phi1, phi2), time_range_algo, time_arboricity_algo]



def build_csv_nine_real_instances():
    resu = [["name", "s1", "s2", "size_si", "nbBPs1", "nbBPs2", "s1nosameBPs", "s2nosameBPs", "barrier", "range", "arboricityL","arboricityR","arboricity", "time_range_algo(secs)", "time_arboricity_algo(secs)"]]

    for (name, s1full, s2full) in tab:
        print("Now working on ", name, "\n")
        resu.append(compute_barrier_and_infos(name, s1full, s2full))

    with open('BarriersForRNAExplorerData.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile, delimiter=' ',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
        for elem in resu:
            writer.writerow(elem)

build_csv_nine_real_instances()