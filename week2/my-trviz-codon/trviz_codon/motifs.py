# motifs.py
# from decomposer import Decomposer, is_valid_sequence
# from trviz_codon.decomposer import Decomposer, is_valid_sequence
from .decomposer import Decomposer, is_valid_sequence
from typing import List

class Motif:
    def __init__(self, sequence: str):
        self.sequence = sequence.upper()
        if not is_valid_sequence(self.sequence):
            raise ValueError(f"Invalid motif sequence: {self.sequence}")

def parse_motif(sequence: str) -> Motif:
    """Create a Motif object from a sequence string."""
    return Motif(sequence)

def validate_sequence(sequence: str):
    """Validate DNA sequence contains only A, C, G, T."""
    if not is_valid_sequence(sequence.upper()):
        raise ValueError(f"Sequence has invalid characters: {sequence}")

def decompose(sequence: str, motifs: List[str], mode: str = "DP", **kwargs) -> List[str]:
    """Decompose the sequence into motifs using the Decomposer."""
    decomposer = Decomposer(mode=mode)
    return decomposer.decompose(sequence, motifs, **kwargs)

# For convenience, expose a minimal set of the original Biopython-like API
class Motifs:
    def __init__(self, motifs: List[str]):
        self.motifs = [Motif(m) for m in motifs]

    def decompose(self, sequence: str, **kwargs) -> List[str]:
        motif_seqs = [m.sequence for m in self.motifs]
        return decompose(sequence, motif_seqs, **kwargs)
