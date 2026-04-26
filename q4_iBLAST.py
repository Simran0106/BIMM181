import tracemalloc
import matplotlib.pyplot as plt
from q1_locAL import align
import numpy as np

def read_fasta(filename):
    
    with open(filename, 'r') as f:
        lines = f.readlines()
    return "".join([line.strip() for line in lines if not line.startswith(">")])



def run_linear_space(seq_q, seq_d, match, mismatch, indel, threshold):
    len_q = len(seq_q)
    len_d = len(seq_d)
    
    
    score_matrix = np.zeros((2, len_q + 1), dtype=np.int32)
    
   
    start_i = np.zeros((2, len_q + 1), dtype=np.int32)
    start_j = np.zeros((2, len_q + 1), dtype=np.int32)
    
    hits_dict = {}

    print(f"Aligning {len_q} x {len_d} matrix...")

    for i in range(1, len_d + 1):
        i2 = i % 2
        i1 = (i - 1) % 2
        char_d = seq_d[i-1]
        
        for j in range(1, len_q + 1):
            match_score = match if char_d == seq_q[j-1] else mismatch
            
            s_diag = score_matrix[i1, j-1] + match_score
            s_up   = score_matrix[i1, j] + indel
            s_left = score_matrix[i2, j-1] + indel
            
            max_s = max(0, s_diag, s_up, s_left)
            score_matrix[i2, j] = max_s
            
            
            if max_s == 0:
                start_i[i2, j], start_j[i2, j] = i, j
            elif max_s == s_diag:
                if score_matrix[i1, j-1] > 0:
                    start_i[i2, j] = start_i[i1, j-1]
                    start_j[i2, j] = start_j[i1, j-1]
                else:
                    start_i[i2, j], start_j[i2, j] = i-1, j-1
            elif max_s == s_up:
                start_i[i2, j] = start_i[i1, j]
                start_j[i2, j] = start_j[i1, j]
            else: 
                start_i[i2, j] = start_i[i2, j-1]
                start_j[i2, j] = start_j[i2, j-1]
                
            if max_s >= threshold:
                si, sj = start_i[i2, j], start_j[i2, j]
                coord = (si, sj)
                if coord not in hits_dict or max_s > hits_dict[coord]['score']:
                    hits_dict[coord] = {
                        'score': max_s,
                        'start_d': si, 'end_d': i,
                        'start_q': sj, 'end_q': j
                    }
        
       
        if i % 20000 == 0:
            print(f"Progress: {(i/len_d)*100:.1f}% done...")

    return list(hits_dict.values())
            
          

def prune_hits(hits):
    
    sorted_hits = sorted(hits, key=lambda x: x['score'], reverse=True)
    final_hits = []
    
    
    for hit in sorted_hits:
        keep = True
        for s in final_hits:
            s_start, s_end = s['start_d'], s['end_d']
            hit_start, hit_end = hit['start_d'], hit['end_d']
            
            
            overlap_start = max(s_start, hit_start)
            overlap_end = min(s_end, hit_end)
            overlap_len = max(0, overlap_end - overlap_start)
            
            s_len = s_end - s_start
            
            
            if s_len > 0 and (overlap_len / s_len) >= 0.5:
                keep = False
                break
                
        if keep:
            final_hits.append(hit)
            
    return final_hits

def main():
    
    tracemalloc.start()
    print("Reading files")
    
    
    seq_q = read_fasta("query.fa") 
    seq_d = read_fasta("database.fa") 
    
    current, peak = tracemalloc.get_traced_memory()
    print(f"Memory after reading: Current = {current / 10**6:.2f} MB, Peak = {peak / 10**6:.2f} MB")
    
    
    print("\nStarting Linear Space Alignment ")
    tracemalloc.reset_peak()
    
    
    raw_hits = run_linear_space(seq_q, seq_d, 1, -2, -2, threshold=15)
    
    current, peak = tracemalloc.get_traced_memory()
    print(f"Memory during alignment: Current = {current / 10**6:.2f} MB, Peak = {peak / 10**6:.2f} MB")
    tracemalloc.stop()

    
    final_hits = prune_hits(raw_hits)
    print(f"\nFiltered down to {len(final_hits)} distinct final hits.")

    
    thresholds = list(range(15, max(h['score'] for h in final_hits) + 2))
    hit_counts = [sum(1 for h in final_hits if h['score'] >= x) for x in thresholds]

    plt.figure(figsize=(8, 5))
    plt.plot(thresholds, hit_counts, marker='o', color='purple')
    plt.title('Number of Distinct Hits $\ge$ Score Threshold')
    plt.xlabel('Score Threshold ($x$)')
    plt.ylabel('Number of Hits')
    plt.grid(True)
    plt.show() 

    
    print("\n TOP 15 ALIGNMENTS ")
    for idx, hit in enumerate(final_hits[:15]):
        
        sub_q = seq_q[hit['start_q'] : hit['end_q']]
        sub_d = seq_d[hit['start_d'] : hit['end_d']]
        
        
        res = align(sub_q, sub_d, 1, -2, -2)
        
        print(f"\nRank {idx+1} | Score: {hit['score']} | Query: {hit['start_q']} to {hit['end_q']} | DB: {hit['start_d']} to {hit['end_d']}")
        
        
        print(f"Query : {res['q_res']}")
        print(f"        {res['mid']}")
        print(f"DB    : {res['db_res']}")

if __name__ == "__main__":
    main()