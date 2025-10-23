#!/usr/bin/env python3
"""
Sequence alignment algorithms implementation in Python
Supports: Global, Local, Semi-global (Fitting), and Affine Gap alignments
Author: Week 4 Deliverable
"""

def read_fasta(filename):
    """Read sequences from FASTA file"""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:].split()[0]  # Take only first word as ID
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences


def global_alignment(seq1, seq2, match=3, mismatch=-3, gap=-2):
    """
    Needleman-Wunsch global alignment algorithm
    Returns the alignment score
    """
    m = len(seq1)
    n = len(seq2)
    
    # Initialize DP matrix
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    
    # Initialize first row and column
    for i in range(m + 1):
        dp[i][0] = i * gap
    for j in range(n + 1):
        dp[0][j] = j * gap
    
    # Fill DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                score_match = match
            else:
                score_match = mismatch
            
            dp[i][j] = max(
                dp[i-1][j-1] + score_match,  # diagonal (match/mismatch)
                dp[i-1][j] + gap,             # up (gap in seq2)
                dp[i][j-1] + gap              # left (gap in seq1)
            )
    
    return dp[m][n]


def local_alignment(seq1, seq2, match=3, mismatch=-3, gap=-2):
    """
    Smith-Waterman local alignment algorithm
    Returns the best local alignment score
    """
    m = len(seq1)
    n = len(seq2)
    
    # Initialize DP matrix
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    
    max_score = 0
    
    # Fill DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                score_match = match
            else:
                score_match = mismatch
            
            dp[i][j] = max(
                0,                            # no negative scores
                dp[i-1][j-1] + score_match,  # diagonal
                dp[i-1][j] + gap,             # up
                dp[i][j-1] + gap              # left
            )
            
            if dp[i][j] > max_score:
                max_score = dp[i][j]
    
    return max_score


def semiglobal_alignment(seq1, seq2, match=3, mismatch=-3, gap=-2):
    """
    Semi-global (fitting) alignment algorithm
    Fits seq2 into seq1 without penalizing end gaps in seq1
    Returns the best fitting alignment score
    """
    m = len(seq1)
    n = len(seq2)
    
    # Initialize DP matrix
    dp = [[0 for _ in range(n + 1)] for _ in range(m + 1)]
    
    # Initialize first column (no penalty for gaps at start of seq1)
    for i in range(m + 1):
        dp[i][0] = 0
    
    # Initialize first row (penalty for gaps in seq2)
    for j in range(1, n + 1):
        dp[0][j] = j * gap
    
    # Fill DP matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                score_match = match
            else:
                score_match = mismatch
            
            dp[i][j] = max(
                dp[i-1][j-1] + score_match,
                dp[i-1][j] + gap,
                dp[i][j-1] + gap
            )
    
    # Find max score in last column (no penalty for gaps at end of seq1)
    max_score = dp[0][n]
    for i in range(1, m + 1):
        if dp[i][n] > max_score:
            max_score = dp[i][n]
    
    return max_score


def affine_gap_alignment(seq1, seq2, match=3, mismatch=-3, gap_open=-5, gap_extend=-1):
    """
    Global alignment with affine gap penalty
    Gap cost = gap_open + k * gap_extend (where k is gap length)
    Returns the alignment score
    """
    m = len(seq1)
    n = len(seq2)
    
    NEG_INF = -999999999
    
    # Three matrices: M (match), I (gap in seq1), J (gap in seq2)
    M = [[NEG_INF for _ in range(n + 1)] for _ in range(m + 1)]
    I = [[NEG_INF for _ in range(n + 1)] for _ in range(m + 1)]
    J = [[NEG_INF for _ in range(n + 1)] for _ in range(m + 1)]
    
    # Initialize
    M[0][0] = 0
    
    # First row - gaps in seq1
    for j in range(1, n + 1):
        J[0][j] = gap_open + j * gap_extend
    
    # First column - gaps in seq2
    for i in range(1, m + 1):
        I[i][0] = gap_open + i * gap_extend
    
    # Fill matrices
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i-1] == seq2[j-1]:
                score_match = match
            else:
                score_match = mismatch
            
            # M: match/mismatch (can come from any state)
            M[i][j] = max(
                M[i-1][j-1] + score_match,
                I[i-1][j-1] + score_match,
                J[i-1][j-1] + score_match
            )
            
            # I: gap in seq1 (vertical - can extend or open)
            I[i][j] = max(
                M[i-1][j] + gap_open + gap_extend,
                I[i-1][j] + gap_extend
            )
            
            # J: gap in seq2 (horizontal - can extend or open)
            J[i][j] = max(
                M[i][j-1] + gap_open + gap_extend,
                J[i][j-1] + gap_extend
            )
    
    # Return best score from any of the three states
    return max(M[m][n], I[m][n], J[m][n])


def run_alignment(method, seq1, seq2, **kwargs):
    """Run specified alignment method"""
    if method == 'global':
        return global_alignment(seq1, seq2)
    elif method == 'local':
        return local_alignment(seq1, seq2)
    elif method == 'semiglobal':
        return semiglobal_alignment(seq1, seq2)
    elif method == 'affine':
        return affine_gap_alignment(seq1, seq2, **kwargs)
    else:
        raise ValueError(f"Unknown method: {method}")


import sys
import time

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python3 alignment.py <method> <file1> <file2>")
        print("Methods: global, local, semiglobal, affine")
        sys.exit(1)
    
    method = sys.argv[1]
    file1 = sys.argv[2]
    file2 = sys.argv[3]
    
    # Read sequences
    try:
        seqs1 = read_fasta(file1)
        seqs2 = read_fasta(file2)
    except (FileNotFoundError, IOError) as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    
    if not seqs1 or not seqs2:
        print("Error: No sequences found in input files", file=sys.stderr)
        sys.exit(1)
    
    # Get first sequence from each file
    seq1 = list(seqs1.values())[0]
    seq2 = list(seqs2.values())[0]
    
    # Run alignment
    start = time.time()
    
    if method == 'affine':
        score = run_alignment(method, seq1, seq2, gap_open=-5, gap_extend=-1)
    else:
        score = run_alignment(method, seq1, seq2)
    
    elapsed = (time.time() - start) * 1000  # Convert to milliseconds
    
    print(f"Score: {score}")
    print(f"Time: {elapsed:.2f}ms")
    sys.stdout.flush()