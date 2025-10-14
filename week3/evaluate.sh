#!/bin/bash

# Check if we're already in week3 directory, if not cd into it
if [ -d "week3" ]; then
    cd week3
elif [ -f "evaluate.sh" ]; then
    # We're already in week3 directory
    echo "Already in week3 directory"
else
    echo "Error: Cannot find week3 directory"
    exit 1
fi

echo "Language    Runtime"
echo "-------------------"

# Run Python version
echo "Running Python tests..."
python_output=$(python -c "
import time
import numpy as np

def upgma_simple(distances):
    MAX_FLOAT = float('inf')
    n = distances.shape[0]
    
    cluster_assignments = [[i] for i in range(n)]
    cluster_sizes = np.ones(n, dtype=int)
    distances_v = distances.astype(float, copy=True)
    
    for step in range(n - 1):
        dist_min = MAX_FLOAT
        i_min, j_min = -1, -1
        
        for i in range(n):
            for j in range(i):
                if distances_v[i, j] < dist_min:
                    dist_min = distances_v[i, j]
                    i_min, j_min = i, j
        
        if i_min == -1:
            break
        
        cluster_assignments[i_min].extend(cluster_assignments[j_min])
        cluster_sizes[i_min] += cluster_sizes[j_min]
        cluster_assignments[j_min] = []
        
        for k in range(n):
            if k != i_min and k != j_min and cluster_sizes[k] > 0:
                mean_val = (
                    distances_v[i_min, k] * cluster_sizes[i_min] +
                    distances_v[j_min, k] * cluster_sizes[j_min]
                ) / (cluster_sizes[i_min] + cluster_sizes[j_min])
                distances_v[i_min, k] = mean_val
                distances_v[k, i_min] = mean_val
        
        distances_v[j_min, :] = MAX_FLOAT
        distances_v[:, j_min] = MAX_FLOAT
    
    return len(cluster_assignments[0])

n = 200
distances = np.random.rand(n, n) * 50
distances = (distances + distances.T) / 2
np.fill_diagonal(distances, 0.0)

start_time = time.time()
result = upgma_simple(distances)
end_time = time.time()

print(int((end_time - start_time) * 1000))
")

python_time=$python_output

# Run Codon version
echo "Running Codon tests..."
if [ -f "phylo_test.py" ]; then
    codon_output=$(codon run --release phylo_test.py 2>&1)
    codon_time=$codon_output
else
    echo "Error: phylo_test.py not found"
    codon_time="ERROR"
fi

# Validate results are numbers
if ! [[ "$python_time" =~ ^[0-9]+$ ]]; then
    python_time="ERROR"
fi

if ! [[ "$codon_time" =~ ^[0-9]+$ ]]; then
    codon_time="ERROR"
fi

echo "python      ${python_time}ms"
echo "codon       ${codon_time}ms"