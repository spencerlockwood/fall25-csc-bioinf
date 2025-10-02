from typing import List, Dict, Tuple, Set, Optional
import numpy as np

# Import Python-only modules if needed
# from python import Counter as PyCounter
# from python import defaultdict as PyDefaultdict

DP_MODULE: str = "DP"

def is_valid_sequence(sequence: str) -> bool:
    """Check if sequence contains only ACGT"""
    DNA_CHARS = {'A', 'C', 'G', 'T'}
    for c in sequence:
        if c not in DNA_CHARS:
            return False
    return True

def get_motifs_from_visited_states_and_region(visited_states: List[str], region: str) -> List[str]:
    """Placeholder - needs full implementation from utils"""
    # This would need the full implementation from utils.py
    return []

class Decomposer:
    mode: str
    
    def __init__(self, mode: str = DP_MODULE):
        if mode == "DP_CY":
            self.mode = mode
        elif mode == "DP":
            self.mode = mode
        elif mode == "HMM":
            self.mode = mode
        else:
            raise ValueError(f"{mode} is invalid mode for tandem repeat decomposer.")

    @staticmethod
    def refine(decomposed_trs: List[List[str]], verbose: bool = False) -> List[List[str]]:
        """
        Refine the decomposed TRs to remove redundant motifs.
        """
        motif_pair_counter: Dict[Tuple[str, str], int] = {}
        motif_pair_str_counter: Dict[str, int] = {}
        motif_pair_str_to_motif_pair: Dict[str, Set[Tuple[str, str]]] = {}

        # Count motif pairs
        for tr in decomposed_trs:
            for i in range(len(tr) - 1):
                first_motif: str = tr[i]
                second_motif: str = tr[i + 1]
                motif_pair: Tuple[str, str] = (first_motif, second_motif)
                motif_pair_str: str = first_motif + second_motif

                if motif_pair not in motif_pair_counter:
                    motif_pair_counter[motif_pair] = 0
                motif_pair_counter[motif_pair] += 1
                
                if motif_pair_str not in motif_pair_str_counter:
                    motif_pair_str_counter[motif_pair_str] = 0
                motif_pair_str_counter[motif_pair_str] += 1
                
                if motif_pair_str not in motif_pair_str_to_motif_pair:
                    motif_pair_str_to_motif_pair[motif_pair_str] = set()
                motif_pair_str_to_motif_pair[motif_pair_str].add(motif_pair)

        refined_trs: List[List[str]] = []
        for tr in decomposed_trs:
            for i in range(len(tr) - 1):
                first_motif: str = tr[i]
                second_motif: str = tr[i + 1]
                motif_pair: Tuple[str, str] = (first_motif, second_motif)

                if motif_pair not in motif_pair_counter or motif_pair_counter[motif_pair] == 0:
                    continue

                pair_str: str = first_motif + second_motif
                if motif_pair_counter[motif_pair] < motif_pair_str_counter[pair_str]:
                    max_frequency: int = 0
                    max_frequency_motif_pair: Tuple[str, str] = ("", "")
                    
                    for mp in motif_pair_str_to_motif_pair[pair_str]:
                        if motif_pair_counter[mp] > max_frequency:
                            max_frequency = motif_pair_counter[mp]
                            max_frequency_motif_pair = mp
                    
                    if verbose:
                        print("Multiple pairs found", motif_pair_str_to_motif_pair[pair_str])
                        print("Max frequency motif pair:", max_frequency_motif_pair)

                    tr[i] = max_frequency_motif_pair[0]
                    tr[i + 1] = max_frequency_motif_pair[1]
            refined_trs.append(tr)

        return refined_trs

    def decompose(self, sequence: str, motifs, 
                 match_score: float = 1.0,
                 mismatch_score: float = -1.0,
                 insertion_score: float = -1.0,
                 deletion_score: float = -1.0,
                 min_score_threshold: float = float("-inf"),
                 verbose: bool = False) -> List[str]:
        """
        Decompose sequence into motifs
        """
        if not isinstance(sequence, str):
            raise TypeError("Sequence must be a string")
        
        motif_list: List[str]
        if isinstance(motifs, str):
            motif_list = [motifs]
        elif isinstance(motifs, list):
            motif_list = motifs
        else:
            raise TypeError("Motifs must be a list of strings")

        sequence = sequence.upper()
        motif_list = [m.upper() for m in motif_list]

        if not is_valid_sequence(sequence):
            raise ValueError(f"Sequence has invalid characters: {sequence}")

        for motif in motif_list:
            if not is_valid_sequence(motif):
                raise ValueError(f"The motif has invalid characters: {motif}")

        if self.mode == "DP_CY":
            raise NotImplementedError("DP_CY mode not available in Codon")
        elif self.mode == "DP":
            return self._decompose_dp(sequence, motif_list, 
                                     match_score, mismatch_score, 
                                     insertion_score, deletion_score,
                                     min_score_threshold, verbose)
        else:
            return self._decompose_hmm(sequence, motif_list, verbose)

    @staticmethod
    def _decompose_dp(
            sequence: str,
            motifs: List[str],
            match_score: float = 1.0,
            mismatch_score: float = -1.0,
            insertion_score: float = -1.0,
            deletion_score: float = -1.0,
            min_score_threshold: float = float("-inf"),
            verbose: bool = False
    ) -> List[str]:
        """
        Decompose sequence into motifs using dynamic programming
        """
        # Get max motif length
        max_motif_length: int = 0
        for motif in motifs:
            if len(motif) > max_motif_length:
                max_motif_length = len(motif)
        
        if verbose:
            motif_str = ','.join(motifs)
            print("Motifs used for decomposition: " + motif_str)
            print("Max motif length " + str(max_motif_length))

        # Initialize DP arrays
        seq_len: int = len(sequence)
        num_motifs: int = len(motifs)
        
        s = np.zeros((seq_len + 1, num_motifs, max_motif_length + 1), dtype=np.float64)
        backtrack_i = np.zeros((seq_len + 1, num_motifs, max_motif_length + 1), dtype=np.int32)
        backtrack_m = np.zeros((seq_len + 1, num_motifs, max_motif_length + 1), dtype=np.int32)
        backtrack_j = np.zeros((seq_len + 1, num_motifs, max_motif_length + 1), dtype=np.int32)

        # Boundary cases
        for m in range(num_motifs):
            motif: str = motifs[m]
            motif_len: int = len(motif)
            
            for i in range(seq_len + 1):
                for j in range(motif_len + 1):
                    if i == 0 and j == 0:
                        s[0, m, 0] = 0.0
                        backtrack_i[0, m, j] = 0
                        backtrack_m[0, m, j] = m
                        backtrack_j[0, m, j] = 0
                    elif i == 0 and j != 0:
                        s[0, m, j] = s[0, m, j - 1] + insertion_score
                        backtrack_i[0, m, j] = 0
                        backtrack_m[0, m, j] = m
                        backtrack_j[0, m, j] = j - 1
                    elif i != 0 and j == 0:
                        s[i, m, 0] = s[i - 1, m, 0] + insertion_score
                        backtrack_i[i, m, 0] = i - 1
                        backtrack_m[i, m, 0] = m
                        backtrack_j[i, m, 0] = 0

        # Normal cases
        for i in range(1, seq_len + 1):
            for m in range(num_motifs):
                motif: str = motifs[m]
                motif_len: int = len(motif)
                
                for j in range(1, motif_len + 1):
                    if j == 1:
                        if i == 1:
                            match_val: float = match_score if sequence[i - 1] == motif[j - 1] else mismatch_score
                            from_diagonal: float = s[i - 1, m, j - 1] + match_val
                            from_m_left: float = s[i - 1, m, j] + insertion_score
                            from_m_up: float = s[i, m, j - 1] + deletion_score

                            s[i, m, j] = max(from_diagonal, max(from_m_left, from_m_up))
                            
                            if from_diagonal >= from_m_left and from_diagonal >= from_m_up:
                                backtrack_i[i, m, j] = 0
                                backtrack_m[i, m, j] = m
                                backtrack_j[i, m, j] = 0
                            elif from_m_left >= from_m_up:
                                backtrack_i[i, m, j] = i - 1
                                backtrack_m[i, m, j] = m
                                backtrack_j[i, m, j] = j
                            else:
                                backtrack_i[i, m, j] = i
                                backtrack_m[i, m, j] = m
                                backtrack_j[i, m, j] = 0
                        else:
                            max_motif_val: float = float('-inf')
                            max_m_index: int = -1
                            max_j_of_max_m: int = -1
                            
                            for mi in range(num_motifs):
                                ms: str = motifs[mi]
                                m_end: float = s[i - 1, mi, len(ms)]
                                if m_end > max_motif_val:
                                    max_motif_val = m_end
                                    max_m_index = mi
                                    max_j_of_max_m = len(ms)

                            match_val2: float = match_score if sequence[i - 1] == motif[0] else mismatch_score
                            max_from_motif_end: float = max_motif_val + match_val2
                            from_m_left2: float = s[i - 1, m, 1] + insertion_score
                            from_m_up2: float = s[i, m, 0] + deletion_score

                            s[i, m, j] = max(max_from_motif_end, max(from_m_left2, from_m_up2))
                            
                            if max_from_motif_end >= from_m_left2 and max_from_motif_end >= from_m_up2:
                                backtrack_i[i, m, j] = i - 1
                                backtrack_m[i, m, j] = max_m_index
                                backtrack_j[i, m, j] = max_j_of_max_m
                            elif from_m_left2 >= from_m_up2:
                                backtrack_i[i, m, j] = i - 1
                                backtrack_m[i, m, j] = m
                                backtrack_j[i, m, j] = 1
                            else:
                                backtrack_i[i, m, j] = i
                                backtrack_m[i, m, j] = m
                                backtrack_j[i, m, j] = 0
                    else:
                        match_val3: float = match_score if sequence[i - 1] == motif[j - 1] else mismatch_score
                        diagonal: float = s[i - 1, m, j - 1] + match_val3
                        from_left: float = s[i - 1, m, j] + insertion_score
                        from_up: float = s[i, m, j - 1] + deletion_score

                        s[i, m, j] = max(diagonal, max(from_left, from_up))
                        
                        if diagonal >= from_left and diagonal >= from_up:
                            backtrack_i[i, m, j] = i - 1
                            backtrack_m[i, m, j] = m
                            backtrack_j[i, m, j] = j - 1
                        elif from_left >= from_up:
                            backtrack_i[i, m, j] = i - 1
                            backtrack_m[i, m, j] = m
                            backtrack_j[i, m, j] = j
                        else:
                            backtrack_i[i, m, j] = i
                            backtrack_m[i, m, j] = m
                            backtrack_j[i, m, j] = j - 1

        # Backtracking
        backtrack_max: float = min_score_threshold
        backtrack_start_i: int = -1
        backtrack_start_m: int = -1
        backtrack_start_j: int = -1
        
        for m in range(num_motifs):
            motif: str = motifs[m]
            if backtrack_max < s[seq_len, m, len(motif)]:
                backtrack_max = s[seq_len, m, len(motif)]
                backtrack_start_i = seq_len
                backtrack_start_m = m
                backtrack_start_j = len(motif)

        if backtrack_start_i == -1:
            err_msg = "No good match greater than score threshold of " + str(min_score_threshold)
            raise ValueError(err_msg)

        if verbose:
            print("Best score: " + str(backtrack_max))

        # Reconstruct path
        i: int = backtrack_start_i
        m: int = backtrack_start_m
        j: int = backtrack_start_j
        prev_i: int = -1
        prev_j: int = -1
        decomposed_motif: str = ""
        decomposed_motifs: List[str] = []

        while True:
            if verbose:
                print("Backtrack pointer " + str(i) + " " + str(m) + " " + str(j))

            if prev_j == 1 and j != 1:
                if verbose:
                    print("Decomposed motif: " + decomposed_motif[::-1])
                decomposed_motifs.append(decomposed_motif[::-1])
                decomposed_motif = ""
            
            if prev_i != i and i != 0:
                decomposed_motif += sequence[i - 1]

            if i == 0 and j == 0:
                break
            
            next_i: int = int(backtrack_i[i, m, j])
            next_m: int = int(backtrack_m[i, m, j])
            next_j: int = int(backtrack_j[i, m, j])
            
            prev_i = i
            prev_j = j
            i = next_i
            m = next_m
            j = next_j

        if verbose:
            input_str = ''.join(decomposed_motifs[::-1])
            decomp_str = ' '.join(decomposed_motifs[::-1])
            print("Input     : " + input_str)
            print("Decomposed: " + decomp_str)

        return decomposed_motifs[::-1]

    def _decompose_hmm(self, sequence: str, consensus_motif: List[str], verbose: bool = False) -> List[str]:
        """
        Decompose sequence into motifs using a HMM
        Note: This requires pomegranate from Python
        """
        from python import pomegranate
        
        consensus_motif_str: str = consensus_motif[0] if len(consensus_motif) > 0 else ""
        repeat_count: int = int(round(float(len(sequence)) / float(len(consensus_motif_str))))

        if verbose:
            print("Estimated repeat count " + str(repeat_count))
            print("Building HMM...")
            print("motif " + consensus_motif_str)
            print("sequence " + sequence)

        # HMM building would go here - requires pomegranate
        raise NotImplementedError("HMM decomposition requires full pomegranate integration")