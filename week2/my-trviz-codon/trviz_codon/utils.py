# utils.py

from typing import List, Dict, Tuple, Set

# Determine if we're running in Codon
try:
    __codon__
    CODON = True
except NameError:
    CODON = False

# Conditionally import Python-only libraries
if CODON:
    from python import SeqIO as PySeqIO
    from python import itertools as PyItertools
else:
    from Bio import SeqIO as PySeqIO
    import itertools as PyItertools

# Constants
LOWERCASE_LETTERS: str = "abcdefghijklmnopqrstuvwxyz"
UPPERCASE_LETTERS: str = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
DIGITS: str = "0123456789"
DNA_CHARACTERS: Set[str] = {'A', 'C', 'G', 'T'}



def get_sample_and_sequence_from_fasta(fasta_file: str) -> Tuple[List[str], List[str]]:
    """Read fasta file and output headers and sequences"""
    headers: List[str] = []
    sequences: List[str] = []
    
    # Use Python's SeqIO
    with open(fasta_file) as handle:
        for record in PySeqIO.parse(handle, "fasta"):
            headers.append(str(record.id))
            sequences.append(str(record.seq).upper())

    return headers, sequences


def get_motif_counter(decomposed_vntrs: List[List[str]]) -> Dict[str, int]:
    """Return a counter for each motif"""
    motif_counter: Dict[str, int] = {}
    
    for decomposed_vntr in decomposed_vntrs:
        for motif in decomposed_vntr:
            if motif not in motif_counter:
                motif_counter[motif] = 0
            motif_counter[motif] += 1

    return motif_counter


def is_emitting_state(state_name: str) -> bool:
    """Check if the given state is emitting state, that is insertion or matching state"""
    if state_name.startswith('M') or state_name.startswith('I') or \
       state_name.startswith('start_random_matches') or \
       state_name.startswith('end_random_matches'):
        return True
    return False


def get_repeating_pattern_lengths(visited_states: List[str]) -> List[int]:
    """Get lengths of repeating patterns from visited states"""
    lengths: List[int] = []
    prev_start: int = -1
    
    for i in range(len(visited_states)):
        if visited_states[i].startswith('unit_end') and prev_start != -1:
            current_len: int = 0
            for j in range(prev_start, i):
                if is_emitting_state(visited_states[j]):
                    current_len += 1
            lengths.append(current_len)
        if visited_states[i].startswith('unit_start'):
            prev_start = i
    
    return lengths


def get_motifs_from_visited_states_and_region(visited_states: List[str], region: str) -> List[str]:
    """Extract motifs from visited states and region"""
    lengths: List[int] = get_repeating_pattern_lengths(visited_states)
    repeat_segments: List[str] = []
    added: int = 0
    
    for l in lengths:
        repeat_segments.append(region[added:added + l])
        added += l
    
    return repeat_segments


def is_valid_sequence(sequence: str) -> bool:
    """Check if the given sequence is DNA sequence"""
    for s in sequence:
        if s not in DNA_CHARACTERS:
            return False
    return True


def sort_by_manually(aligned_vntrs: List[str], sample_ids: List[str], 
                    sample_order_file: str) -> Tuple[List[str], List[str]]:
    """Sort the aligned and encoded tandem repeats based on the given order"""
    if sample_order_file is None:
        print("Sample order file is not provided. Follow the given sample order.")
        return sample_ids, aligned_vntrs

    sample_order: List[str] = []
    with open(sample_order_file) as f:
        for line in f:
            sample_order.append(line.strip())
    
    sorted_sample_ids: List[str] = []
    sorted_aligned_vntrs: List[str] = []
    
    for sample_id in sample_order:
        if sample_id in sample_ids:
            idx: int = sample_ids.index(sample_id)
            sorted_sample_ids.append(sample_id)
            sorted_aligned_vntrs.append(aligned_vntrs[idx])

    return sorted_sample_ids, sorted_aligned_vntrs


def sort(aligned_vntrs: List[str], sample_ids: List[str], 
        symbol_to_motif: Dict[str, str], sample_order_file: str, 
        method: str = 'motif_count') -> Tuple[List[str], List[str]]:
    """Sort the aligned and encoded tandem repeats"""
    if method == 'name':
        # Create pairs and sort
        pairs: List[Tuple[str, str]] = []
        for i in range(len(sample_ids)):
            pairs.append((sample_ids[i], aligned_vntrs[i]))
        pairs.sort(key=lambda x: x[0])
        
        new_ids: List[str] = []
        new_vntrs: List[str] = []
        for p in pairs:
            new_ids.append(p[0])
            new_vntrs.append(p[1])
        return new_ids, new_vntrs
        
    elif method == 'motif_count':
        pairs2: List[Tuple[str, str, int]] = []
        for i in range(len(sample_ids)):
            count: int = len(aligned_vntrs[i].replace('-', ''))
            pairs2.append((sample_ids[i], aligned_vntrs[i], count))
        pairs2.sort(key=lambda x: x[2], reverse=True)
        
        new_ids2: List[str] = []
        new_vntrs2: List[str] = []
        for p in pairs2:
            new_ids2.append(p[0])
            new_vntrs2.append(p[1])
        return new_ids2, new_vntrs2
        
    elif method == 'simulated_annealing':
        return sort_by_simulated_annealing_optimized(aligned_vntrs, sample_ids, symbol_to_motif)
    elif method == 'manually':
        return sort_by_manually(aligned_vntrs, sample_ids, sample_order_file)
    else:
        raise ValueError("Please check the rearrangement method. {}".format(method))


def get_levenshtein_distance(s1: str, s2: str) -> int:
    """
    Calculate the Levenshtein distance between two strings.
    The minimum number of single-character edits (insertions, deletions or substitutions)
    required to change one string into the other.
    """
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances: List[int] = list(range(len(s1) + 1))
    
    for i2 in range(len(s2)):
        c2: str = s2[i2]
        distances_: List[int] = [i2 + 1]
        
        for i1 in range(len(s1)):
            c1: str = s1[i1]
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                min_val: int = min(distances[i1], min(distances[i1 + 1], distances_[-1]))
                distances_.append(1 + min_val)
        distances = distances_
    
    return distances[-1]


def _calculate_cost(seq1: str, seq2: str, alphabet_to_motif: Dict[str, str]) -> int:
    """Calculate cost between two sequences"""
    if len(seq1) != len(seq2):
        raise ValueError("The length of two sequences should be identical.")

    cost: int = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            if seq1[i] != '-' and seq2[i] != '-':
                s1: str = alphabet_to_motif[seq1[i].lower()]
                s2: str = alphabet_to_motif[seq2[i].lower()]
                cost += get_levenshtein_distance(s1, s2)
            else:
                if seq1[i] == '-':
                    cost += len(alphabet_to_motif[seq2[i].lower()])
                else:
                    cost += len(alphabet_to_motif[seq1[i].lower()])
    return cost


def calculate_cost_with_dist_matrix(aligned_encoded_vntr1: str, 
                                   aligned_encoded_vntr2: str, 
                                   dist_matrix: Dict[str, Dict[str, int]], 
                                   allow_copy_change: bool = False) -> int:
    """Calculate cost using precomputed distance matrix"""
    if len(aligned_encoded_vntr1) != len(aligned_encoded_vntr2):
        raise ValueError("The length of two sequences should be identical.")
    
    cost: int = 0
    for i in range(len(aligned_encoded_vntr1)):
        symbol1: str = aligned_encoded_vntr1[i]
        symbol2: str = aligned_encoded_vntr2[i]
        
        if symbol1 != symbol2:
            if symbol1 != '-' and symbol2 != '-':
                cost += dist_matrix[symbol1][symbol2]
            else:
                if symbol1 == '-':
                    if allow_copy_change:
                        cost += 1
                    else:
                        cost += dist_matrix[symbol2][symbol2]
                else:
                    if allow_copy_change:
                        cost += 1
                    else:
                        cost += dist_matrix[symbol1][symbol1]

    return cost


def calculate_cost(aligned_vntrs: List[str], alphabet_to_motif: Dict[str, str]) -> int:
    """Calculate total cost for a list of aligned VNTRs"""
    total_cost: int = 0
    for i in range(len(aligned_vntrs) - 1):
        total_cost += _calculate_cost(aligned_vntrs[i], aligned_vntrs[i + 1], alphabet_to_motif)

    return total_cost


def get_distance_matrix(symbol_to_motif: Dict[str, str], score: bool = False) -> Dict[str, Dict[str, int]]:
    """
    Stores the edit distance between a motif and another motif.
    if two motifs are the same (e.g. dist_matrix[motif_x][motif_x]) it stores the length of the motif.
    """
    dist_matrix: Dict[str, Dict[str, int]] = {}
    max_score: int = 5

    for symbol1 in symbol_to_motif:
        dist_matrix[symbol1] = {}
        for symbol2 in symbol_to_motif:
            motif_seq1: str = symbol_to_motif[symbol1]
            motif_seq2: str = symbol_to_motif[symbol2]
            
            if symbol1 == symbol2:
                dist_matrix[symbol1][symbol2] = len(motif_seq1)
            else:
                edit_dist: int = get_levenshtein_distance(motif_seq1, motif_seq2)
                if score:
                    max_len: int = max(len(motif_seq1), len(motif_seq2))
                    dist_matrix[symbol1][symbol2] = int(max_score * (1.0 - float(edit_dist) / float(max_len)))
                else:
                    dist_matrix[symbol1][symbol2] = edit_dist

    return dist_matrix


def get_score_matrix(symbol_to_motif: Dict[str, str],
                     match_score: int = 2,
                     mismatch_score_for_edit_dist_of_1: int = -1,
                     mismatch_score_for_edit_dist_greater_than_1: int = -2,
                     gap_open_penalty: float = 1.5,
                     gap_extension_penalty: float = 0.6) -> Dict[str, float]:
    """Get scoring matrix for alignment"""
    score_matrix: Dict[str, float] = {}
    score_matrix['gap_open'] = gap_open_penalty
    score_matrix['gap_extension'] = gap_extension_penalty

    for symbol1 in symbol_to_motif:
        score_matrix[symbol1] = {}
        for symbol2 in symbol_to_motif:
            motif_seq1: str = symbol_to_motif[symbol1]
            motif_seq2: str = symbol_to_motif[symbol2]
            
            if symbol1 == symbol2:
                score_matrix[symbol1][symbol2] = float(match_score)
            else:
                edit_dist: int = get_levenshtein_distance(motif_seq1, motif_seq2)
                edit_dist_cutoff: int = 1
                
                if abs(len(motif_seq1) - len(motif_seq2)) <= 1:
                    max_len: int = max(len(motif_seq2), len(motif_seq1))
                    edit_dist_cutoff += max_len // 30
                
                if edit_dist <= edit_dist_cutoff:
                    score_matrix[symbol1][symbol2] = float(mismatch_score_for_edit_dist_of_1)
                else:
                    score_matrix[symbol1][symbol2] = float(mismatch_score_for_edit_dist_greater_than_1)

    return score_matrix


def calculate_total_cost(aligned_vntrs: List[str], dist_matrix: Dict[str, Dict[str, int]]) -> int:
    """Calculate total cost for aligned VNTRs using distance matrix"""
    total_cost: int = 0
    for i in range(len(aligned_vntrs) - 1):
        total_cost += calculate_cost_with_dist_matrix(aligned_vntrs[i], aligned_vntrs[i + 1], dist_matrix)

    return total_cost


def sort_by_simulated_annealing_optimized(seq_list: List[str], sample_ids: List[str], 
                                         symbol_to_motif: Dict[str, str]) -> Tuple[List[str], List[str]]:
    """Sort sequences using simulated annealing optimization"""
    dist_matrix: Dict[str, Dict[str, int]] = get_distance_matrix(symbol_to_motif)

    initial_cost: int = calculate_total_cost(seq_list, dist_matrix)
    initial_seq_list: List[str] = seq_list.copy()
    initial_sample_ids: List[str] = sample_ids.copy()

    T: float = 1000.0
    DECAY: float = 0.9
    iteration: int = 0

    # Generate all index pairs
    all_index_pairs: List[Tuple[int, int]] = []
    for i in range(len(seq_list)):
        for j in range(i + 1, len(seq_list)):
            all_index_pairs.append((i, j))

    while True:
        iteration += 1
        print("T:", T)
        if T <= 1e-2:
            break
        print("iteration", iteration)

        for pair in all_index_pairs:
            index_1: int = pair[0]
            index_2: int = pair[1]
            
            current_cost: int = 0
            after_cost: int = 0

            # Flanking cost for the index_1 sequence - Right side
            current_cost += calculate_cost_with_dist_matrix(seq_list[index_1], seq_list[index_1 + 1], dist_matrix)
            if index_1 + 1 == index_2:
                after_cost += calculate_cost_with_dist_matrix(seq_list[index_2], seq_list[index_1], dist_matrix)
            else:
                after_cost += calculate_cost_with_dist_matrix(seq_list[index_2], seq_list[index_1 + 1], dist_matrix)
            
            if index_1 != 0:  # has left side
                current_cost += calculate_cost_with_dist_matrix(seq_list[index_1], seq_list[index_1 - 1], dist_matrix)
                after_cost += calculate_cost_with_dist_matrix(seq_list[index_2], seq_list[index_1 - 1], dist_matrix)

            # Flanking cost for the index_2 sequence - Left side
            current_cost += calculate_cost_with_dist_matrix(seq_list[index_2], seq_list[index_2 - 1], dist_matrix)
            if index_2 - 1 == index_1:
                after_cost += calculate_cost_with_dist_matrix(seq_list[index_1], seq_list[index_2], dist_matrix)
            else:
                after_cost += calculate_cost_with_dist_matrix(seq_list[index_1], seq_list[index_2 - 1], dist_matrix)
            
            if index_2 != len(seq_list) - 1:  # Right side
                current_cost += calculate_cost_with_dist_matrix(seq_list[index_2], seq_list[index_2 + 1], dist_matrix)
                after_cost += calculate_cost_with_dist_matrix(seq_list[index_1], seq_list[index_2 + 1], dist_matrix)

            if after_cost == current_cost:
                continue
            elif after_cost < current_cost:
                print("Swap occurred at {} and {}".format(index_1, index_2), "after cost", after_cost, "cur cost", current_cost)
                seq_list[index_1], seq_list[index_2] = seq_list[index_2], seq_list[index_1]
                sample_ids[index_1], sample_ids[index_2] = sample_ids[index_2], sample_ids[index_1]
            else:
                prob: float = np.exp(-(float(after_cost - current_cost)) * 10.0 / T)
                rand_val: float = np.random.uniform(0.0, 1.0)
                if prob > rand_val:
                    print("Swap occurred {}".format(prob), "after cost", after_cost, "cur cost", current_cost)
                    seq_list[index_1], seq_list[index_2] = seq_list[index_2], seq_list[index_1]
                    sample_ids[index_1], sample_ids[index_2] = sample_ids[index_2], sample_ids[index_1]

        T *= DECAY

    print("The initial cost", initial_cost)
    after_cost_final: int = calculate_total_cost(seq_list, dist_matrix)
    print("Cost after sorting", after_cost_final)

    if initial_cost < after_cost_final:
        return initial_sample_ids, initial_seq_list
    return sample_ids, seq_list


def add_padding(encoded_trs: List[str]) -> List[str]:
    """
    Add padding to encoded traces to make them all the same length.
    The padding is done by adding '-' to the end of each trace.
    """
    max_motif_count: int = 0
    for encoded_tr in encoded_trs:
        if len(encoded_tr) > max_motif_count:
            max_motif_count = len(encoded_tr)
    
    padded_trs: List[str] = []
    for encoded_tr in encoded_trs:
        padding_count: int = max_motif_count - len(encoded_tr)
        padded_trs.append(encoded_tr + '-' * padding_count)

    return padded_trs


def print_progress_bar(iteration: int, total: int, prefix: str = '', suffix: str = '', 
                      decimals: int = 1, length: int = 100, fill: str = '█', 
                      print_end: str = "\r"):
    """
    Call in a loop to create terminal progress bar
    """
    percent_val: float = 100.0 * (float(iteration) / float(total))
    percent_str: str = ("{0:." + str(decimals) + "f}").format(percent_val)
    filled_length: int = int(length * iteration // total)
    bar: str = fill * filled_length + '-' * (length - filled_length)
    print(f'\r{prefix} |{bar}| {percent_str}% {suffix}', end=print_end)
    
    if iteration == total:
        print()


def get_motif_marks(sample_ids: List[str], decomposed_trs: List[List[str]], 
                   region_prediction_file: str) -> Dict[str, str]:
    """
    Parse the region prediction file and store the result in a dictionary
    The format of the file:
    >sample_id
    region1,start,end \t region2,start,end
    """
    region_prediction: Dict[str, List[Tuple[str, int, int]]] = {}
    
    with open(region_prediction_file, 'r') as f:
        current_sample: str = ""
        for line in f:
            line_stripped: str = line.strip()
            if line_stripped.startswith('>'):
                current_sample = line_stripped[1:]  # remove the '>' sign
                region_prediction[current_sample] = []
            else:
                regions: List[str] = line_stripped.split('\t')
                for region in regions:
                    parts: List[str] = region.split(',')
                    region_name: str = parts[0]
                    start: int = int(parts[1])
                    end: int = int(parts[2])
                    region_prediction[current_sample].append((region_name, start, end))

    # get the motif marks
    motif_marks: Dict[str, str] = {}
    
    for idx in range(len(sample_ids)):
        sample_id: str = sample_ids[idx]
        decomposed_tr: List[str] = decomposed_trs[idx]
        motif_mark: str = ""
        cumulative_length: int = 0

        if sample_id not in region_prediction:
            print("Sample {} is not in the region prediction file".format(sample_id))
            continue

        for motif_seq in decomposed_tr:
            motif_start: int = cumulative_length
            motif_end: int = cumulative_length + len(motif_seq)
            found_intron: bool = False
            
            for region_tuple in region_prediction[sample_id]:
                region_name: str = region_tuple[0]
                start: int = region_tuple[1]
                end: int = region_tuple[2]
                
                if region_name == "intron":
                    # Check if the motif sequence has overlap with any intron regions
                    if (start <= motif_start < end) or (start <= motif_end < end) or \
                       (motif_start <= start and end < motif_end):
                        motif_mark += "I"
                        found_intron = True
                        break
            
            if not found_intron:
                motif_mark += "X"

            cumulative_length += len(motif_seq)

        motif_marks[sample_id] = motif_mark
        if len(motif_mark) != len(decomposed_tr):
            raise ValueError("The length of the motif mark is not equal to the number of motifs")
    
    return motif_marks


def get_sample_to_population(population_data: str, sep: str = '\t', 
                            sample_index: int = 0, 
                            population_index: int = 5) -> Dict[str, str]:
    """
    Parse the population data file and store the result in a dictionary
    """
    sample_to_population: Dict[str, str] = {}
    
    with open(population_data, "r") as f:
        for line in f:
            split: List[str] = line.strip().split(sep)
            sample: str = split[sample_index]
            population: str = split[population_index]
            sample_to_population[sample] = population

    return sample_to_population