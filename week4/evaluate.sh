#!/bin/bash
# Main evaluation script for alignment algorithms

# Check if data directory exists
if [ ! -d "data" ]; then
    echo "Error: data directory not found!"
    exit 1
fi

# Check if all required data files exist
required_files=("MT-human.fa" "MT-orang.fa" "q1.fa" "q2.fa" "q3.fa" "q4.fa" "q5.fa" "t1.fa" "t2.fa" "t3.fa" "t4.fa" "t5.fa")
for file in "${required_files[@]}"; do
    if [ ! -f "data/$file" ]; then
        echo "Error: data/$file not found!"
        exit 1
    fi
done

# Create output header
echo "Method            Language    Runtime"
echo "--------------------------------------"

# Function to run a single test and return time in ms
run_single_test() {
    local method=$1
    local lang=$2
    local file1=$3
    local file2=$4
    
    if [ "$lang" = "python" ]; then
        output=$(python3 alignment.py "$method" "$file1" "$file2" 2>&1)
    else
        output=$(codon run alignment.codon "$method" "$file1" "$file2" 2>&1)
    fi
    
    # Extract time from output
    time_ms=$(echo "$output" | grep "Time:" | awk '{print $2}' | sed 's/ms//')
    
    if [ -z "$time_ms" ]; then
        echo "0"
    else
        echo "$time_ms"
    fi
}

# Test global alignment on MT-human
echo -n "Testing global-mt_human with python... " >&2
python_time=$(run_single_test "global" "python" "data/MT-human.fa" "data/MT-orang.fa")
printf "global-mt_human   python      %6.0fms\n" "$python_time"

echo -n "Testing global-mt_human with codon... " >&2
codon_time=$(run_single_test "global" "codon" "data/MT-human.fa" "data/MT-orang.fa")
printf "global-mt_human   codon       %6.0fms\n" "$codon_time"

# Test local alignment on MT
echo -n "Testing local-mt_human with python... " >&2
python_time=$(run_single_test "local" "python" "data/MT-human.fa" "data/MT-orang.fa")
printf "local-mt_human    python      %6.0fms\n" "$python_time"

echo -n "Testing local-mt_human with codon... " >&2
codon_time=$(run_single_test "local" "codon" "data/MT-human.fa" "data/MT-orang.fa")
printf "local-mt_human    codon       %6.0fms\n" "$codon_time"

# Test semiglobal alignment on MT
echo -n "Testing semiglobal-mt_human with python... " >&2
python_time=$(run_single_test "semiglobal" "python" "data/MT-human.fa" "data/MT-orang.fa")
printf "semiglobal-mt_human python    %6.0fms\n" "$python_time"

echo -n "Testing semiglobal-mt_human with codon... " >&2
codon_time=$(run_single_test "semiglobal" "codon" "data/MT-human.fa" "data/MT-orang.fa")
printf "semiglobal-mt_human codon     %6.0fms\n" "$codon_time"

# Test affine alignment on MT
echo -n "Testing affine-mt_human with python... " >&2
python_time=$(run_single_test "affine" "python" "data/MT-human.fa" "data/MT-orang.fa")
printf "affine-mt_human   python      %6.0fms\n" "$python_time"

echo -n "Testing affine-mt_human with codon... " >&2
codon_time=$(run_single_test "affine" "codon" "data/MT-human.fa" "data/MT-orang.fa")
printf "affine-mt_human   codon       %6.0fms\n" "$codon_time"

# Test all q/t pairs
for i in {1..5}; do
    # Global
    echo -n "Testing global-q${i} with python... " >&2
    python_time=$(run_single_test "global" "python" "data/q${i}.fa" "data/t${i}.fa")
    printf "global-q${i}        python      %6.0fms\n" "$python_time"
    
    echo -n "Testing global-q${i} with codon... " >&2
    codon_time=$(run_single_test "global" "codon" "data/q${i}.fa" "data/t${i}.fa")
    printf "global-q${i}        codon       %6.0fms\n" "$codon_time"
    
    # Local
    echo -n "Testing local-q${i} with python... " >&2
    python_time=$(run_single_test "local" "python" "data/q${i}.fa" "data/t${i}.fa")
    printf "local-q${i}         python      %6.0fms\n" "$python_time"
    
    echo -n "Testing local-q${i} with codon... " >&2
    codon_time=$(run_single_test "local" "codon" "data/q${i}.fa" "data/t${i}.fa")
    printf "local-q${i}         codon       %6.0fms\n" "$codon_time"
    
    # Semiglobal
    echo -n "Testing semiglobal-q${i} with python... " >&2
    python_time=$(run_single_test "semiglobal" "python" "data/q${i}.fa" "data/t${i}.fa")
    printf "semiglobal-q${i}    python      %6.0fms\n" "$python_time"
    
    echo -n "Testing semiglobal-q${i} with codon... " >&2
    codon_time=$(run_single_test "semiglobal" "codon" "data/q${i}.fa" "data/t${i}.fa")
    printf "semiglobal-q${i}    codon       %6.0fms\n" "$codon_time"
    
    # Affine
    echo -n "Testing affine-q${i} with python... " >&2
    python_time=$(run_single_test "affine" "python" "data/q${i}.fa" "data/t${i}.fa")
    printf "affine-q${i}        python      %6.0fms\n" "$python_time"
    
    echo -n "Testing affine-q${i} with codon... " >&2
    codon_time=$(run_single_test "affine" "codon" "data/q${i}.fa" "data/t${i}.fa")
    printf "affine-q${i}        codon       %6.0fms\n" "$codon_time"
done

echo "" >&2
echo "✓ Evaluation complete!" >&2