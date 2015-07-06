


from random import choice


WC = dict(zip("ATGC", "TACG"))


def complement(seq):
    return "".join(WC[n] for n in seq)

def make_pair():
    first = choice("ATGC")
    opts = "".join(n for n in "ATGC" if WC[n] != first)
    second = choice(opts)
    return (first, second)

def make_sequence(lenght=8):
    return tuple(make_pair() for i in range(lenght))

def sequence_is_ok(sequence):
    """ Sequence is a list of pairs: [(A, G), (A, A), ...] """
    seq_sets = tuple(set(pair) for pair in sequence)
    # Check for palindromes
    # or partial palindromes

    # check for sliding

    # Alternatively, have a scoring method

    return True


def main(args):
    seqs = set()
    while len(seqs) < 10:
        seq = make_sequence()
        if sequence_is_ok(seq):
            seqs.add(seq)
    print("Sequences:")
    print("\n".join(", ".join("{}:{}".format(*pair) for pair in seq) for seq in seqs))
    with open(args.get("outputfn", "cc_sequences.txt"), "w") as fd:
        fd.write("\n".join(str(seq) for seq in seqs))



if __name__ == '__main__':
    main({})
