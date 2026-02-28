#!/bin/bash
# Usage: ./run_threads.sh <graph_file> <source> <delta>
# Runs the "delta" executable with different thread counts and reports execution time.

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 graph_file source delta"
    exit 1
fi

GRAPH=$1
SRC=$2
DELTA=$3

# make sure executable exists
if [ ! -x ./delta ]; then
    echo "delta executable not found or not executable. Compile it first."
    exit 2
fi

for threads in 2 4 6 8; do
    echo "Running with $threads threads..."
    # run program and capture its output
    output=$(./delta "$GRAPH" "$SRC" "$DELTA" $threads 2>&1)
    # extract the timing line assuming format "Parallel Delta-Stepping took: X seconds"
    elapsed=$(echo "$output" | grep -oP '(?<=took: )[^ ]+')
    echo "$elapsed"
done
