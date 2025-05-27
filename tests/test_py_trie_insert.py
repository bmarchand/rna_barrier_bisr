from bisr_dpw.dpw_algorithm import Trie, seq_to_int


def test_insert_print():

    seq1 = [1, 3, 4, 2]
    seq2 = [1, 3, 5, 2]

    trie = Trie()

    trie.insert(seq1)
    trie.insert(seq2)

    trie.print()


def test_dfs_order():

    trie = Trie()

    trie.insert([1, 2, 3])

    trie.insert([1, 2, 4])

    sequences = trie.generate_sequences()

    n1 = seq_to_int(sequences[0])
    n2 = seq_to_int(sequences[1])

    assert(n1 < n2)
