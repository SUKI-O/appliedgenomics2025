import argparse
import math
import sys
import zlib

"""
Compute Jaccard coefficient between two sequences (FASTA or plain) using modimizers and report:
- Jaccard coefficient (modimizers only)
- ANI by kth-root
- Approximate ANI using natural log
- Number of modimizers found in each sequence

Usage: python compute_jaccard_modimizer.py -a seqA.fa -b seqB.fa -k 21 -m 100
"""

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
    seq = "".join(s).upper()
    seq = "".join(ch if ch in VALID else "N" for ch in seq)
    return seq


def modimizers(seq, k, m):
    """
    Return set of k-mers whose CRC32 hash mod m == 0 (modimizers).
    """
    if k <= 0 or m <= 0:
        return set()
    n = len(seq)
    out = set()
    for i in range(0, n - k + 1):
        kmer = seq[i : i + k]
        if "N" in kmer:
            continue
        h = zlib.crc32(kmer.encode('utf-8')) & 0xffffffff
        if h % m == 0:
            out.add(kmer)
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("-a", required=True, help="Sequence A (FASTA or plain)")
    p.add_argument("-b", required=True, help="Sequence B (FASTA or plain)")
    p.add_argument("-k", required=True, type=int, help="k-mer length")
    p.add_argument("-m", required=True, type=int, help="modimizer modulus (sampling level)")
    args = p.parse_args()

    A = read_sequence(args.a)
    B = read_sequence(args.b)
    k = args.k
    m = args.m
    if k <= 0:
        print("k must be positive", file=sys.stderr)
        sys.exit(1)
    if m <= 0:
        print("m must be positive", file=sys.stderr)
        sys.exit(1)

    setA = modimizers(A, k, m)
    setB = modimizers(B, k, m)

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
        ani_ln = math.exp(math.log(jaccard) / k)

    print(f"K={k}, M={m}")
    print(f"|A modimizers| = {len(setA)}")
    print(f"|B modimizers| = {len(setB)}")
    print(f"|A ∩ B|   = {len_inter}")
    print(f"|A ∪ B|   = {len_union}")
    print(f"Jaccard  = {jaccard:.6f}")
    print(f"ANI (kth-root) = {ani_root:.6f}")
    print(f"ANI (ln approx) = {ani_ln:.6f}")


if __name__ == "__main__":
    main()
