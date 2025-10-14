# week3/phylo_test.py
import time
import numpy as np

def upgma_simple(distances):
    """A simplified UPGMA that performs the clustering algorithm"""
    MAX_FLOAT = np.finfo(np.float64).max
    n = distances.shape[0]
    
    # Track cluster assignments
    cluster_assignments = [[i] for i in range(n)]
    cluster_sizes = np.ones(n, dtype=int)
    distances_v = distances.astype(float, copy=True)
    
    # Perform clustering iterations
    for step in range(n - 1):
        # Find minimum distance
        dist_min = MAX_FLOAT
        i_min = -1
        j_min = -1
        
        for i in range(n):
            for j in range(i):
                if distances_v[i, j] < dist_min:
                    dist_min = distances_v[i, j]
                    i_min = i
                    j_min = j
        
        if i_min == -1:
            break
        
        # Merge clusters
        cluster_assignments[i_min].extend(cluster_assignments[j_min])
        cluster_sizes[i_min] += cluster_sizes[j_min]
        cluster_assignments[j_min] = []
        
        # Update distance matrix using UPGMA formula
        for k in range(n):
            if k != i_min and k != j_min and cluster_sizes[k] > 0:
                mean_val = (
                    distances_v[i_min, k] * cluster_sizes[i_min] +
                    distances_v[j_min, k] * cluster_sizes[j_min]
                ) / (cluster_sizes[i_min] + cluster_sizes[j_min])
                distances_v[i_min, k] = mean_val
                distances_v[k, i_min] = mean_val
        
        # Mark merged cluster as invalid
        distances_v[j_min, :] = MAX_FLOAT
        distances_v[:, j_min] = MAX_FLOAT
    
    return len(cluster_assignments[0])  # Return size of root cluster

def run_benchmark():
    """Run the UPGMA benchmark"""
    # Test with substantial data
    n = 200
    distances = np.random.rand(n, n) * 50
    distances = (distances + distances.T) / 2
    np.fill_diagonal(distances, 0.0)
    
    start_time = time.time()
    result = upgma_simple(distances)
    end_time = time.time()
    
    runtime_ms = int((end_time - start_time) * 1000)
    return runtime_ms

if __name__ == "__main__":
    runtime = run_benchmark()
    print(runtime)