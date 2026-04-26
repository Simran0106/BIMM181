import tracemalloc
import matplotlib.pyplot as plt
from q1_locAL import align
import sys
import numpy as np

def read_fasta(filename):
    seq = []
    with open(filename, 'r') as f:
        for line in f:
            if not line.startswith(">"):
                seq.append(line.strip())
    return "".join(seq)

def run_linear_space(seq_q, seq_d, match, mismatch, indel, threshold):
    len_q, len_d = len(seq_q), len(seq_d)
    
    score_curr = np.zeros(len_q + 1, dtype=np.int32)
    score_prev = np.zeros(len_q + 1, dtype=np.int32)
    
    start_i_curr = np.zeros(len_q + 1, dtype=np.int32) 
    start_j_curr = np.zeros(len_q + 1, dtype=np.int32)
    start_i_prev = np.zeros(len_q + 1, dtype=np.int32)
    start_j_prev = np.zeros(len_q + 1, dtype=np.int32)

    
    q_ord = np.frombuffer(seq_q.encode(), dtype=np.int8)
    d_ord = np.frombuffer(seq_d.encode(), dtype=np.int8)
    j_range = np.arange(len_q, dtype=np.int32)   
    hits_dict = {}

    print(f"Starting Alignment: {len_q} x {len_d}")

    for i in range(1, len_d + 1):
        
        match_scores = np.where(q_ord == d_ord[i-1], match, mismatch).astype(np.int32)
        diag_pot = score_prev[:-1] + match_scores   
        up_pot   = score_prev[1:]  + indel           

        score_curr[0] = 0
        score_curr[1:] = np.maximum(0, np.maximum(diag_pot, up_pot))

        
        use_diag = (score_curr[1:] > 0) & (diag_pot >= up_pot)
        fresh    = score_prev[:-1] == 0             # diagonal came from a zero cell

        start_i_curr[1:] = np.where(use_diag,
                               np.where(fresh, i - 1, start_i_prev[:-1]),
                               np.where(score_curr[1:] > 0, start_i_prev[1:], 0))
        start_j_curr[1:] = np.where(use_diag,
                               np.where(fresh, j_range, start_j_prev[:-1]),
                               np.where(score_curr[1:] > 0, start_j_prev[1:], 0))
        start_i_curr[0] = start_j_curr[0] = 0

        
        for j in range(1, len_q + 1):
            left_val = score_curr[j-1] + indel
            if left_val > score_curr[j]:
                score_curr[j]    = max(0, left_val)
                start_i_curr[j] = start_i_curr[j-1]
                start_j_curr[j] = start_j_curr[j-1]

        
        if int(score_curr.max()) >= threshold:
            idxs = np.where(score_curr >= threshold)[0]
            for j in idxs:
                j   = int(j)
                res = int(score_curr[j])
                coord = (int(start_i_curr[j]), int(start_j_curr[j]))
                if coord not in hits_dict or res > hits_dict[coord]['score']:
                    hits_dict[coord] = {
                        'score':   res,
                        'start_d': int(start_i_curr[j]),
                        'end_d':   i,
                        'start_q': int(start_j_curr[j]),
                        'end_q':   j,
                    }

        
        score_prev,   score_curr   = score_curr,   score_prev
        start_i_prev, start_i_curr = start_i_curr, start_i_prev
        start_j_prev, start_j_curr = start_j_curr, start_j_prev

        if i % 10000 == 0:
            sys.stdout.write(f"\rProgress: {(i/len_d)*100:.1f}% | Unique Hits: {len(hits_dict)}")
            sys.stdout.flush()

    print("\nAlignment Complete.")
    return list(hits_dict.values())

def prune_hits(hits):
    sorted_hits = sorted(hits, key=lambda x: x['score'], reverse=True)
    final_hits = []
    for hit in sorted_hits:
        keep = True
        for s in final_hits:
            overlap_start = max(s['start_d'], hit['start_d'])
            overlap_end = min(s['end_d'], hit['end_d'])
            overlap_len = max(0, overlap_end - overlap_start)
            s_len = s['end_d'] - s['start_d']
            if s_len > 0 and (overlap_len / s_len) >= 0.5:
                keep = False
                break
        if keep:
            final_hits.append(hit)
    return final_hits

def main():
    tracemalloc.start()
    seq_q = read_fasta("query.fa") 
    seq_d = read_fasta("database.fa") 
    
    tracemalloc.reset_peak()
    raw_hits = run_linear_space(seq_q, seq_d, 1, -2, -2, threshold=15)
    
    current, peak = tracemalloc.get_traced_memory()
    print(f"\nMemory: {current / 10**6:.2f} MB | Peak: {peak / 10**6:.2f} MB")
    tracemalloc.stop()

    final_hits = prune_hits(raw_hits)
    
    if not final_hits:
        return

    scores = [h['score'] for h in final_hits]
    thresholds = list(range(15, max(scores) + 2))
    hit_counts = [sum(1 for s in scores if s >= x) for x in thresholds]

    plt.figure(figsize=(8, 5))
    plt.plot(thresholds, hit_counts, marker='o', color='purple')
    plt.title('Number of Distinct Hits vs Score Threshold')
    plt.xlabel('Score Threshold')
    plt.ylabel('Number of Hits')
    plt.grid(True)
    plt.savefig('q4_hits_plot.png')

    for idx, hit in enumerate(final_hits[:15]):
        sub_q = seq_q[hit['start_q'] : hit['end_q']]
        sub_d = seq_d[hit['start_d'] : hit['end_d']]
        res = align(sub_q, sub_d, 1, -2, -2)
        print(f"\nRank {idx+1} | Score: {hit['score']} | Query: {hit['start_q']}-{hit['end_q']} | DB: {hit['start_d']}-{hit['end_d']}")
        print(f"Query : {res['q_res']}")
        print(f"        {res['mid']}")
        print(f"DB    : {res['db_res']}")

if __name__ == "__main__":
    main()