

"""
Need to know: (Investigate, then ask William)

* Is A:C pair equivalent to C:A ? Or "slightly equivalent"?

* Do we need to design pair-sequences explicitly to not dimerize,
    or will the anti-nucleation hairpins prevent this?

* Do we need to check for complementarity on non-perpendicular angles,
    e.g. a bundle that hybridizes out of lattice?

* Which criteria do we select against? And where to we set the criteria cutoff? (If using binary criteria)
    - Sliding
    - Reversability
    - Self dimerization? Palindromic sequences
    - Hetero dimerization? Cross check
    - Uniqueness to other parallel bundles? (Cross check, forward and reverse)
    - Uniqueness to complementary bundles? (Cross check, complement, forward and reverse)

* Do we need to have a scoring algorithm, or just binary criteria? (accept/reject)

* Is random generation of bundle-sequences followed by selection "good enough",
    or do we need to implement evolution-type optimization of the sequences?
    The sequences all depend on each other, so evolution could be tricky or produce weird results...

* How much should we pay attention to "ensemble" weaknesses,
    e.g. the fact that one weak bundle itself doesn't matter too much (e.g. sliding),
    but if the whole ensemble/unit can slide (or reverse), then that is much worse.
    Gut feeling: We don't have to worry too much at first, since we can have enough sequence space
    to produce good sequences for individual bundles.
    (And if they are all good, then ensemble should also be OK.)

"""


from random import choice


WC = dict(zip("ATGC", "TACG"))


def complement(seq):
    return "".join(WC[n] for n in seq)

def make_pair():
    first = choice("ATGC")
    opts = "".join(n for n in "ATGC" if WC[n] != first)
    second = choice(opts)
    return (first, second)

def make_sequence(length=8):
    return tuple(make_pair() for i in range(length))


#def hairpin_combinations(seq):
#    """
#    Return a generator of all possible palindromic combinations of seq.
#    There are both even and odd combinations.
#    """
#    # Even:
#    for i in range(2, len(seq)-2):
#        # even:
#        yield zip(seq[:i], seq[i:])
#        # Odd:
#        yield zip(seq[:i], seq[i+1:])
#
#def calc_hairpin_frac(seq_sets):
#    pali_combs = list(palindrome_combinations(seq_sets))
#    pali_parities = [[pair1 == pair2 for pair1, pair2 in comb] for comb in pali_combs]
#    pali_fracs = [sum(palindrome_parity)/len(seq_sets) for palindrome_parity in pali_parities]
#    #palindrome_parity = [pair1 == pair2 for pair1, pair2 in zip(seq_sets, reversed(seq_sets))]
#    #palindrome_frac = sum(palindrome_parity)/len(sequence)
#    max_frac = max(pali_fracs)
#    max_pali = pali_combs[pali_fracs.index(max_frac)]
#    return max_frac, max_pali


def calc_sliding_frac(seq, seq2=None):
    """ Calculate sliding scores. Bulges etc are not allowed. """
    if seq2 is None:
        seq2 = seq
    i_max = min(len(seq), len(seq2)) - 1
    combs = [zip(seq, seq2[i:]) for i in range(1, i_max)]
    parities = [[pair1 == pair2 for pair1, pair2 in comb] for comb in combs]
    fracs = [sum(parity)/len(seq) for parity in parities]
    max_frac = max(fracs)
    max_comb = combs[fracs.index(max_frac)]
    return max_frac, max_comb


def sequence_is_ok(sequence):
    """
    Sequence is a list of pairs: [(A, G), (A, A), ...]
    # Alternatively, have a scoring method
    """
    seq_sets = tuple(set(pair) for pair in sequence)

    # Check for sliding:
    max_frac, max_comb = calc_sliding_frac(seq_sets)
    if max_frac > 0.2:
        print("Disregarding high-sliding ({:.2f}) seq:".format(max_frac), seq_repr(sequence))
        return False

    # Check for palindromes/partial palindromes
    max_pali_frac, max_comb = calc_sliding_frac(seq_sets, seq_sets[::-1])
    if max_pali_frac > 0.2:
        print("Disregarding palindrome-like ({:.2f}) seq:".format(max_pali_frac), seq_repr(sequence))
        return False

    return True


def cross_check_seqs(seq1, seq2):
    """
    Crosscheck whether two sequences will interact.
    We need to check:
    1) Compare the generated sequences to each other.
    2) Compare the complementary bundles to the generated sequences, and each other.

    """
    # seq_sets = tuple(set(pair) for pair in sequence)
    pass


def gen_N_sequences(N=10, length=8):
    """ Generate exactly N number of sequences of length <length>. """
    seqs = set()
    while len(seqs) < N:
        seq = make_sequence(length)
        if sequence_is_ok(seq):
            # TODO: Cross-check seq against existing sequences.
            seqs.add(seq)
    return seqs


def seq_repr(seq):
    return ", ".join("{}:{}".format(*pair) for pair in seq)


def main(args):
    seqs = gen_N_sequences(args.get("length", 8))
    print("Sequences:")
    print("\n".join(seq_repr(seq) for seq in seqs))
    with open(args.get("outputfn", "cc_sequences.txt"), "w") as fd:
        fd.write("\n".join(str(seq) for seq in seqs))



if __name__ == '__main__':
    main({})
