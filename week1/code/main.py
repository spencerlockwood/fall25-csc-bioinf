from dbg import DBG
from utils import read_data
import sys
import os

import faulthandler
faulthandler.enable()
sys.setrecursionlimit(1000000)


def calculate_n50(lengths: list[int]) -> int:
    """Compute the N50 statistic from a list of contig lengths."""
    if not lengths:
        return 0
    total = sum(lengths)
    half = total / 2
    running = 0
    for length in sorted(lengths, reverse=True):
        running += length
        if running >= half:
            return length
    return 0


if __name__ == "__main__":
    argv = sys.argv
    short1, short2, long1 = read_data(argv[1])

    k = 25
    dbg = DBG(k=k, data_list=[short1, short2, long1])

    contig_lengths = []
    contig_file = os.path.join('./', argv[1], 'contig.fasta')

    with open(contig_file, 'w') as f:
        for i in range(20):
            c = dbg.get_longest_contig()
            if c is None:
                break
            contig_lengths.append(len(c))
            print(i, len(c))
            f.write(f'>contig_{i}\n')
            f.write(c + '\n')

    n50 = calculate_n50(contig_lengths)
    print(f"N50 = {n50}")
