
# pylint: disable=C0111,C0103,W0142
"""
Need to know: (Investigate, then ask William)

* Is A:C pair equivalent to C:A ? Or "slightly equivalent"?
    - No, it is not considered identical. The 'A' is at the 5' end of one strand
    and C at the 3' end of another strand. The complement, G:T will bind to A:C,
    but not C:A.
    The 5' end of one DNA strand does not interact with the 5' end of another strand,
    i.e. a dsDNA with 5' sticky end / overhang will not interact with another
    dsDNA with a C overhang at the 5' end.

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



import os
import yaml
from random import choice
from itertools import product
from copy import deepcopy
from datetime import datetime
from math import sqrt


WC = dict(zip("ATGC", "TACG"))

VERBOSE = 1

bases = set("ATGC")


def complement(seq):
    return "".join(WC[n] for n in seq)

def make_pair():
    first = choice("ATGC")  # random.choice returns seq[_randbelow(len(seq))] so must be indexable.
    # If the first is a G (A), the second must NOT be a C (T):
    opts = [n for n in "ATGC" if WC[n] != first]
    second = choice(opts)
    return [first, second]  # use tuple if you want to have sequences as a set, set(seqs)

def make_sequence1(length=8):
    # use tuple if you want to have sequences as a set, set(seqs)
    return [make_pair() for i in range(length)]

def make_sequence2(length=8):
    """
    Combined to eliminate a single function call; increases performance by 10-20 %

    python -m timeit -s "from seqgen import make_sequence, make_sequence2;" "make_sequence()"
    10000 loops, best of 3: 81.8 usec per loop

    python -m timeit -s "from seqgen import make_sequence, make_sequence2;" "make_sequence2()"
    10000 loops, best of 3: 70.6 usec per loop

    Also: making tuples or lists are about the same.
    Note: For pypy, the difference between make_sequence1 and make_sequence2
        is completely lost. - 16-17 usec/loop in both cases.

    """
    # use tuple if you want to have sequences as a set, set(seqs)
    seq_first_bases = [choice("ATGC") for i in range(length)]
    #pairs = [[first, choice(list(bases-{WC[first]}))] # 80 usec/loop (20 usec/loop in pypy)
    pairs = [[first, choice([n for n in "ATGC" if WC[n] != first])]  # 65 usec/loop (15 with pypy)
             for first in seq_first_bases]
    return pairs

make_random_sequence = make_sequence2


def gen_nonrandom_sequence(length=8):
    """
    all possible N-length combinations of
    all_possible_pairs =
    [('A', 'A'),
     ('A', 'T'),
     ('A', 'G'),
     ('A', 'C'),
     ...
     ('C', 'G'),
     ('C', 'C')]
    """
    all_possible_pairs = [list(pair) for pair in product("ATGC", repeat=2)]    # 16
    # Exclude A:T, T:A, G:C, C:G -- 12 pairs left
    all_allowed_pairs = [pair for pair in all_possible_pairs if pair[1] != WC[pair[0]]]
    all_possible_seqs = product(all_allowed_pairs, repeat=length) # 16**length
    return all_possible_seqs


def rcompl_pair(pair):
    return [WC[n] for n in reversed(pair)]


def rcompl_bundle_seq(seq):
    return [rcompl_pair(pair) for pair in reversed(seq)]


def gc_content(seq):
    n_gc = sum(n in "GC" for pair in seq for n in pair)
    return n_gc/(len(seq)*2)


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


def calc_sliding_frac(seq, seq2=None, offset=1):
    """
    Calculate sliding scores. Bulges etc are not allowed.
    Permutations:
    A1 B1 C1 D2 E1 F1    A1 B1 C1 D2 E1 F1	A1 B1 C1 D2 E1 F1	 A1 B1 C1 D2 E1 F1  A1 B1 C1 D2 E1 F1  A1 B1 C1 D2 E1 F1
    A2 B2 C2 D2 E2 F2    B2 C2 D2 E2 F2		C2 D2 E2 F2		 D2 E2 F2           E2 F2              F2
and
    A1 B1 C1 D1 E1 F1    B1 C1 D1 E1 F1		C1 D1 E1 F1		 D1 E1 F1           E1 F1              F1
    A2 B2 C2 D2 E2 F2    A2 B2 C2 D2 E2 F2	A2 B2 C2 D2 E2 F2	 A2 B2 C2 D2 E2 F2  A2 B2 C2 D2 E2 F2  A2 B2 C2 D2 E2 F2

    """
    if seq2 is None:
        seq2 = seq
    i_max = min(len(seq), len(seq2)) - 1
    combs = [zip(seq, seq2[i:]) for i in range(offset, i_max)] \
          + [zip(seq[i:], seq2) for i in range(offset, i_max)]
    parities = [[pair1 == pair2 for pair1, pair2 in comb] for comb in combs]
    fracs = [sum(parity)/len(seq) for parity in parities]
    max_frac = max(fracs)
    max_comb = combs[fracs.index(max_frac)]
    return max_frac, max_comb


def sequence_is_ok(sequence, max_self_sim=0.2, max_pali_sim=0.2, gc_range=(0.2, 0.8)):
    """
    Sequence is a list of pairs: [(A, G), (A, A), ...]
    # Alternatively, have a scoring method
    This just tests whether one sequence will slide or has palindromes.
    It does not compare for sequence complementarity with other sequences.
    """
    #seq_sets = tuple(set(pair) for pair in sequence)
    # the order *does* matter -- i.e. we can use a tuple instead of sets:

    # GC content between a defined fraction:
    gc_frac = gc_content(sequence)
    if not gc_range[0] < gc_frac < gc_range[1]:
        if VERBOSE > 1:
            print("Disregarding seq with high or low GC content ({:.2f}), seq:".format(gc_frac), seq_repr(sequence))
        return False

    # Check for sliding:
    # max_frac, max_comb = calc_sliding_frac(seq_sets)
    # use offset=1, because we dont want to compare seq directly with it self without any sliding.
    max_frac, _ = calc_sliding_frac(sequence, offset=1)
    if max_frac > max_self_sim:
        if VERBOSE > 1:
            print("Disregarding high-sliding ({:.2f}) seq:".format(max_frac), seq_repr(sequence))
        return False

    # Check for palindromes/partial palindromes
    max_pali_frac, _ = calc_sliding_frac(sequence, sequence[::-1], offset=0)
    # calc_sliding_frac doesn't check for complete identity (without sliding)
    # Note that the may be usecases for actually having palindromes.
    if sequence == sequence[::-1]:
        if VERBOSE > 1:
            print("Detected perfect palindrome seq:", seq_repr(sequence))
        return False
        #max_pali_frac = 1
    if max_pali_frac > max_pali_sim:
        if VERBOSE > 1:
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
    similarity = sum(seq1pair == seq2pair for seq1pair, seq2pair in zip(seq1, seq2))
    sim_frac = similarity/len(seq1)
    # Does the reversed complement matter?
    # We have perpendicular binding, not axially, so it really shouldn't.

    # However, we do want to check for "sliding" complementarity:
    max_slide_frac, _ = calc_sliding_frac(seq1, seq2, offset=0)
    # Also check the palindrome?
    max_pali_frac, _ = calc_sliding_frac(seq1, seq2[::-1], offset=0)

    # score:
    sim_score = sim_frac + max((max_slide_frac, max_pali_frac))
    return sim_score


def check_sequence_vs_existing_set(sequence, existing, method="max"):
    """
    Check a single sequence against a set of sequences.
    Will check the sequence for:
    1) Similarity
    -- 2) Complementarity -- Edit: No.
    """
    if not existing:
        return 0
    if method == "max":
        max_sim = max(cross_check_seqs(sequence, seq2) for seq2 in existing)
        return max_sim
    # root-sum-of-squares:
    elif method == "rss":
        rss_sim = sqrt(sum(cross_check_seqs(sequence, seq2)**2 for seq2 in existing))
        return rss_sim


def gen_N_sequences(N=10, length=8, existing=None, randomly=True, max_other_sim=0.2,
                    max_self_sim=0.2, max_pali_sim=0.2, gc_range=(0.2, 0.8),
                    dont_update_existing=False):
    """
    Generate exactly N number of sequences of length <length>.
    """
    if existing is None:
        seqs = []
    else:
        if isinstance(existing, str):
            print("Loading existing sequences from file", existing)
            with open(existing) as fp:
                seqs = yaml.load(fp)
        else:
            seqs = existing.copy()
        print("Basing on %s existing sequences" % len(seqs))

    if dont_update_existing:
        # Newly accepted sequences will not be added to the list of existing sequences
        # This means that future sequences will not be compared against other sequences produced in this run.
        # Useful if you want to generate several candidates to fit into an existing set,
        # but you want to do the actual selection of these by manual review.
        existing = deepcopy(seqs)
    else:
        # The list of existing sequences is the same the seqs list to which newly accepted sequences are added to.
        existing = seqs

    print("seqs:")
    print("\n".join(seq_repr(seq) for seq in seqs))
    N_identical = 0
    N_bad_seq = 0
    N_similar = 0
    n_rounds = 0

    if randomly:
        def sequence_generator():
            while True:
                yield make_random_sequence(length)
        seqgen = sequence_generator()
    else:
        #def sequence_generator():
        #    all_possible_seqs_gen = gen_nonrandom_sequence(length=length)
        #    for seq in all_possible_seqs_gen:
        #        yield seq
        seqgen = gen_nonrandom_sequence(length=length)
    try:
        while len(seqs) < N:
            n_rounds += 1
            if n_rounds % 100000 == 0:
                if n_rounds % 1000000 == 0:
                    print("- %s million rounds, %s seqs..." % (n_rounds/1000000, len(seqs)))
                elif n_rounds < 1000000:
                    print("- %s rounds, %s seqs..." % (n_rounds, len(seqs)))
            #seq = make_sequence(length)
            seq = next(seqgen)
            # Check if sequence is identical to an existing sequence:
            if seq in seqs or seq[::-1] in seqs:
                N_identical += 1
                continue

            # We have at least two further checks:
            # - One is whether the sequence in itself is OK.
            # - The other is whether the sequence is compatible with the existing sequences (i.e. not too similar)
            if not sequence_is_ok(seq, max_self_sim=0.2, max_pali_sim=0.2, gc_range=(0.2, 0.8)):
                N_bad_seq += 1
                continue

            # TODO: Cross-check seq against existing sequences.
            sim_score = check_sequence_vs_existing_set(seq, existing)
            if sim_score > max_other_sim:
                N_similar += 1
                if VERBOSE > 1:
                    print("Disregarding similar (%0.02f) sequence %s (%s seqs in existing set)" %
                          (sim_score, seq_repr(seq), len(seqs)))
                continue

            seqs.append(seq)
            print("- %s seqs..." % len(seqs))
    except StopIteration:
        print("Sequential sequence generator exhausted - %s sequences generated\n" % len(seqs))
    except KeyboardInterrupt:
        print("\nSequence generation aborted! - %s sequences generated\n" % len(seqs))
    print("N identical:", N_identical)
    print("N bad seq:", N_bad_seq)
    print("N similar:", N_similar)
    print("N rounds:", n_rounds)
    return seqs




def seq_repr(seq):
    return ", ".join("{}:{}".format(*pair) for pair in seq)


def main(args):

    fn_N = 10

    existing = None
    #existing = args.get("existing", "cc_sequences5.yaml")
    #existing = "cc_sequences_17.yaml"
    randomly = True
    #randomly = args.get("randomly", False)
    length = args.get("length", 8)
    target_count = args.get("target_count", 40)

    outputfn = args.get("outputfn", "cc_sequences.yaml")

    output_format = args.get("output_format", "yaml")
    uniquefn = args.get("uniquefn", True)

    if output_format == 'yaml' and 'yaml' not in os.path.splitext(outputfn)[1]:
        outputfn += ".yaml"

    if uniquefn:
        fnroot, ext = os.path.splitext(outputfn)
        fnfmt = fnroot + "_%s" + ext
        while True:
            if os.path.exists(outputfn):
                outputfn = fnfmt % fn_N
                fn_N += 1
            else:
                print("Found unique filename:", outputfn)
                break

    print("Outputfn:", outputfn)
    with open(outputfn, 'w') as fp:
        os.utime(outputfn, None) # times must be given specifically as None on pytho 3.2 (=pypy)
    print("Generating sequences... (%s)" % ("randomly" if randomly else "sequential combinations"))

    ### Acceptance criteria:
    criteria = dict(max_other_sim=0.25,
                    max_self_sim=0.25,
                    max_pali_sim=0.25,
                    gc_range=(0.26, 0.74))
    target_number_of_sequences = 50

    paramsfn = os.path.splitext(outputfn)[0] + ".params.yaml"
    with open(paramsfn, 'w') as fp:
        params = {'randomly': randomly,
                  'existing': existing,
                  'outputfn': outputfn,
                  'length': length,
                  'target_count': target_count,
                  'starttime': datetime.now(),
                  'criteria': criteria,
                  'target_number_of_sequences': target_number_of_sequences,
                  'Note': "Updated slide calc (sliding both ways now) and excluding perfect palindromes."
                 }
        yaml.dump(params, fp)

    #seqs = gen_N_sequences(N=40, length=length,
    #                       existing=existing, randomly=randomly)
    seqs = gen_N_sequences(N=target_number_of_sequences,
                           length=length,
                           existing=existing, randomly=randomly,
                           dont_update_existing=False,
                           **criteria)


    print("Sequences:")
    print("\n".join(seq_repr(seq) for seq in seqs))

    print("Saving data to file:", outputfn)

    #outputdir = os.path.dirname(outputfn)
    if not output_format == 'yaml':
        with open(outputfn, "w") as fd:
            fd.write("\n".join(str(seq) for seq in seqs))
    else:
        with open(outputfn, "w") as fd:
            yaml.dump(seqs, fd)


if __name__ == '__main__':
    main({})
