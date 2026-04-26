import random
import sys

def generate_dna(num_seq, length):
    bases = ['A', 'C', 'G', 'T']
    sequences = []
    counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}

    for i in range(num_seq):
        base_list=random.choices(bases, k=length)
        seq = "".join(base_list)
        sequences.append(seq)
        for base in seq:
            counts[base] += 1
        print(seq)


    total = num_seq * length
    print("Frequencies")
    for base in bases:
        freq = counts[base] / total
        print(f"{base}: {freq:.4f}")

    return sequences

if __name__ == "__main__":
    generate_dna(int(sys.argv[1]), int(sys.argv[2]))

