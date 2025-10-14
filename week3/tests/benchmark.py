# week3/tests/benchmark.py
import time

def run_benchmark():
    """Run the tests and return runtime in milliseconds"""
    start_time = time.time()
    
    # Direct imports for Codon compatibility
    from sequence.test_phylo import run_all_tests
    run_all_tests()
    
    end_time = time.time()
    return int((end_time - start_time) * 1000)

if __name__ == "__main__":
    runtime = run_benchmark()
    print(runtime)