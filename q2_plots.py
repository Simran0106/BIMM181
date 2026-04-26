import matplotlib.pyplot as plt
from q1_locAL import align
from q2_randomDNA import generate_dna

def plot():
    num_pairs=500
    seq_length=1000

    lengths_p1=[]
    lengths_p2=[]

    print("Generating 1000 sequences with frequencies of bases")
    all_seqs = generate_dna(num_pairs * 2, seq_length)

    print("Starting alignment")

    for i in range(num_pairs):
        seq_q = all_seqs[i * 2]
        seq_d=all_seqs[(i * 2)+1]

        res_p1 = align(seq_q, seq_d, 1, -30, 0)
        lengths_p1.append(res_p1['length'])
        
        
        res_p2 = align(seq_q, seq_d, 1, -30, -20)
        lengths_p2.append(res_p2['length'])
        
        if (i + 1) % 50 == 0:
            print(f"Completed")

    mean_p1 = sum(lengths_p1) / len(lengths_p1)
    mean_p2 = sum(lengths_p2) / len(lengths_p2)
    
    print(f"\nMean Length P1: {mean_p1:.2f}")
    print(f"Mean Length P2: {mean_p2:.2f}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    ax1.hist(lengths_p1, bins=20, color='blue', edgecolor='black')
    ax1.axvline(mean_p1, color='red', linestyle='dashed', linewidth=2, label=f'Mean: {mean_p1:.1f}')
    ax1.set_title('P1: Match 1, Mismatch -30, Indel 0')
    ax1.set_xlabel('Alignment Length')
    ax1.set_ylabel('Frequency')
    ax1.legend()

    ax2.hist(lengths_p2, bins=20, color='green', edgecolor='black')
    ax2.axvline(mean_p2, color='red', linestyle='dashed', linewidth=2, label=f'Mean: {mean_p2:.1f}')
    ax2.set_title('P2: Match 1, Mismatch -30, Indel -20')
    ax2.set_xlabel('Alignment Length')
    ax2.set_ylabel('Frequency')
    ax2.legend()

    plt.tight_layout()
    plt.show()


plot()


