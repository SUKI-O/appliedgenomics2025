#!/usr/bin/env python
"""
Compute Jaccard coefficient between two sequences (FASTA or plain) and report:
- Jaccard coefficient
- ANI by kth-root
- Approximate ANI using natural log

Usage: python compute_jaccard.py -a seqA.fa -b seqB.fa -k 21
"""
import argparse
import math
import sys

VALID = set("ACGT")


def read_sequence(path):
    s = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            s.append(line)
    seq = "".join(s).upper()  # convert to uppercase
    # convert non-ACGT to N
    seq = "".join(ch if ch in VALID else "N" for ch in seq)
    return seq


def kmers(seq, k):
    if k <= 0:
        return set()
    n = len(seq)
    out = set()
    for i in range(0, n - k + 1):
        kmer = seq[i : i + k]
        if "N" in kmer:
            continue
        out.add(kmer)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-a", required=True, help="Sequence A (FASTA or plain)")
    p.add_argument("-b", required=True, help="Sequence B (FASTA or plain)")
    p.add_argument("-k", required=True, type=int, help="k-mer length")
    args = p.parse_args()

    A = read_sequence(args.a)
    B = read_sequence(args.b)
    k = args.k
    if k <= 0:
        print("k must be positive", file=sys.stderr)
        sys.exit(1)

    setA = kmers(A, k)
    setB = kmers(B, k)

    inter = setA & setB
    union = setA | setB
    len_inter = len(inter)
    len_union = len(union)

    if len_union == 0:
        jaccard = 0.0
    else:
        jaccard = len_inter / len_union

    # ANI by kth-root
    if jaccard <= 0.0:
        ani_root = 0.0
        ani_ln = 0.0
    elif jaccard >= 1.0:
        ani_root = 1.0
        ani_ln = 1.0
    else:
        ani_root = jaccard ** (1.0 / k)
        # approximate via natural log (numerically equivalent, shown separately)
        ani_ln = math.exp(math.log(jaccard) / k)

    print(f"K={k}")
    print(f"|A k-mers| = {len(setA)}")
    print(f"|B k-mers| = {len(setB)}")
    print(f"|A ∩ B|   = {len_inter}")
    print(f"|A ∪ B|   = {len_union}")
    print(f"Jaccard  = {jaccard:.6f}")
    print(f"ANI (kth-root) = {ani_root:.6f}")
    print(f"ANI (ln approx) = {ani_ln:.6f}")


if __name__ == "__main__":
    main()
