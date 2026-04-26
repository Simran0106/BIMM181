import argparse
import numpy as np


def read_fasta(file_path):
    with open(file_path, 'r') as f:
        return "".join(line.strip() for line in f if not line.startswith('>'))

def align(seq_q, seq_d, match, mismatch, indel):

    n, m = len(seq_q), len(seq_d)

    score_matrix=np.zeros((n+1, m+1) , dtype=int)

    max_score=0
    max_pos=(0,0)

    for i in range(1, n+1):
        for j in range(1, m+1):

            match_score= match if seq_q[i-1]==seq_d[j-1] else mismatch

            option=[
                0,
                score_matrix[i-1][j-1]+match_score,
                score_matrix[i][j-1]+ indel,
                score_matrix[i-1][j]+indel

            ]
            curr_val= max(option)
            score_matrix[i,j]=curr_val

            if curr_val>= max_score:
                max_score=curr_val
                max_pos=(i,j)
    
    align_q, align_db, mid_line = "", "", ""
    curr_i, curr_j = max_pos
    end_q, end_db = curr_i - 1, curr_j - 1

    while curr_i>0 and curr_j>0 and score_matrix[curr_i, curr_j]>0:
        curr_val=score_matrix[curr_i, curr_j]

        diag_score=match if seq_q[curr_i-1]==seq_d[curr_j-1] else mismatch
        if curr_val == score_matrix[curr_i-1, curr_j-1] + diag_score:
            align_q = seq_q[curr_i-1] + align_q
            align_db = seq_d[curr_j-1] + align_db
            mid_line = ("|" if seq_q[curr_i-1] == seq_d[curr_j-1] else " ") + mid_line
            curr_i -= 1
            curr_j -= 1
        elif curr_val == score_matrix[curr_i-1, curr_j] + indel:
            align_q = seq_q[curr_i-1] + align_q
            align_db = "-" + align_db
            mid_line = " " + mid_line
            curr_i -= 1
        else:
            align_q = "-" + align_q
            align_db = seq_d[curr_j-1] + align_db
            mid_line = " " + mid_line
            curr_j-= 1
        
    return {
            "begin_q": curr_i, "end_q": end_q,
            "begin_db": curr_j, "end_db": end_db,
            "score": max_score, "length": len(align_q),
            "q_res": align_q, "db_res": align_db, "mid": mid_line
        }



def print_alignment(res):
    
    for i in range(0, len(res['q_res']), 60):
        
        q_start = res['begin_q'] + i
        db_start = res['begin_db'] + i
        
        print(f"\nQuery  {q_start:<5} {res['q_res'][i:i+60]}")
        print(f"             {res['mid'][i:i+60]}")
        print(f"Sbjct  {db_start:<5} {res['db_res'][i:i+60]}")

def main():
    parser = argparse.ArgumentParser(description="locAL: Local Sequence Alignment")
    parser.add_argument("-q", "--query", required=True, help="Path to query FASTA")
    parser.add_argument("-d", "--db", required=True, help="Path to database FASTA")
    parser.add_argument("-m", "--match", type=int, required=True, help="Match score")
    parser.add_argument("-s", "--mismatch", type=int, required=True, help="Mismatch score")
    parser.add_argument("-i", "--indel", type=int, required=True, help="Indel penalty")
    parser.add_argument("-a", action="store_true", help="Option to show alignment")

    args = parser.parse_args()
    seq_query = read_fasta(args.query)
    seq_db = read_fasta(args.db)
    result = align(seq_query, seq_db, args.match, args.mismatch, args.indel)

    print(f"begin-query:{result['begin_q']}, end-query:{result['end_q']}, begin-db:{result['begin_db']}, end-db:{result['end_db']}, score:{result['score']}, length:{result['length']}")


    if args.a:
        print_alignment(result)

if __name__ == "__main__":
    main()
