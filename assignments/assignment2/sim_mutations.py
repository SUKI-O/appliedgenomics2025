import argparse
import random

def read_fasta(filename):
    with open(filename) as f:
        header = ''
        seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    yield header, seq
                header = line
                seq = ''
            else:
                seq += line
        if header:
            yield header, seq

def mutate_sequence(seq, mutation_rate):
    bases = ['A', 'C', 'G', 'T']
    seq_list = list(seq)
    n_mutations = int(len(seq) * mutation_rate)
    positions = random.sample(range(len(seq)), n_mutations)
    for pos in positions:
        original = seq_list[pos].upper()
        if original not in bases:
            continue
        choices = [b for b in bases if b != original]
        seq_list[pos] = random.choice(choices)
    return ''.join(seq_list)

def write_fasta(filename, records):
    with open(filename, 'w') as f:
        for header, seq in records:
            f.write(f"{header}\n")
            for i in range(0, len(seq), 60):
                f.write(seq[i:i+60] + '\n')

def main():
    parser = argparse.ArgumentParser(description="Simulate random substitutions in a FASTA file.")
    parser.add_argument('-i', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output FASTA file')
    parser.add_argument('-m', '--mutation', type=float, required=True, help='Mutation rate (e.g., 0.015 for 1.5%)')
    parser.add_argument('-s', '--seed', type=int, required=True, help='Random seed')
    args = parser.parse_args()

    random.seed(args.seed)
    mutated_records = []
    for header, seq in read_fasta(args.input):
        mutated_seq = mutate_sequence(seq, args.mutation)
        mutated_records.append((header, mutated_seq))
    write_fasta(args.output, mutated_records)

if __name__ == '__main__':
    main()