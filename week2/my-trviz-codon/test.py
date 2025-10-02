# test.py

import trviz_codon.motifs as motifs
import trviz_codon.decomposer as decomposer


def test_motif_creation():
    print("Testing motif creation...", end=" ")
    try:
        m = motifs.Motif("ACGT")
        assert m.sequence == "ACGT", "Motif sequence should be uppercase and preserved"
        
        # Test invalid motif
        try:
            motifs.Motif("ACGX")
            assert False, "Invalid motif sequence should raise ValueError"
        except ValueError:
            pass
            
        print("PASSED")
        return True
    except Exception as e:
        print(f"FAILED: {e}")
        return False

def test_validate_sequence():
    print("Testing sequence validation...", end=" ")
    try:
        valid_seq = "ACGTACGT"
        motifs.validate_sequence(valid_seq)  # Should pass silently
        
        # Test invalid sequence
        try:
            motifs.validate_sequence("ACGTX")
            assert False, "Invalid sequence should raise ValueError"
        except ValueError:
            pass
            
        print("PASSED")
        return True
    except Exception as e:
        print(f"FAILED: {e}")
        return False

def test_decompose_simple():
    print("Testing simple decomposition...", end=" ")
    try:
        seq = "ACGTACGTAC"
        motif_list = ["AC", "GT"]
        # Decompose with default mode DP
        result = motifs.decompose(seq, motif_list)
        assert isinstance(result, list), "Decompose result should be a list"
        for r in result:
            assert isinstance(r, str), "Each decomposed motif should be a string"
            
        print("PASSED")
        return True
    except Exception as e:
        print(f"FAILED: {e}")
        return False

def test_decompose_invalid_sequence():
    print("Testing decomposition with invalid sequence...", end=" ")
    try:
        try:
            motifs.decompose("ACGTX", ["AC", "GT"])
            assert False, "Invalid input sequence should raise ValueError"
        except ValueError:
            pass
            
        print("PASSED")
        return True
    except Exception as e:
        print(f"FAILED: {e}")
        return False

def test_motifs_class():
    print("Testing Motifs class...", end=" ")
    try:
        m = motifs.Motifs(["AC", "GT"])
        seq = "ACGTAC"
        decomposed = m.decompose(seq)
        assert isinstance(decomposed, list), "Motifs.decompose should return a list"
        for item in decomposed:
            assert isinstance(item, str)
            
        print("PASSED")
        return True
    except Exception as e:
        print(f"FAILED: {e}")
        return False

def main():
    print("Starting all tests...")
    print("=" * 50)
    
    test_results = []
    
    test_results.append(test_motif_creation())
    test_results.append(test_validate_sequence())
    test_results.append(test_decompose_simple())
    test_results.append(test_decompose_invalid_sequence())
    test_results.append(test_motifs_class())
    
    print("=" * 50)
    
    passed = sum(test_results)
    total = len(test_results)
    
    if passed == total:
        print(f"All tests passed! ({passed}/{total})")
    else:
        print(f"Some tests failed! ({passed}/{total} passed)")
        
    return passed == total

if __name__ == "__main__":
    main()