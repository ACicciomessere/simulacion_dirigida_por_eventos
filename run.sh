#!/bin/bash
set -euo pipefail

Ns=(300 400 500 600 700)
RUNS=3
ROOT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Compile once
cd "$ROOT_DIR/java"
javac *.java
cd "$ROOT_DIR"

ALL_FILES=()
ALL_NS=()

# Build file lists upfront
for N in "${Ns[@]}"; do
    for (( i=1; i<=RUNS; i++ )); do
        ALL_FILES+=("$ROOT_DIR/runs/N${N}/run_${i}/output.txt")
        ALL_NS+=("$N")
    done
done

# Run all simulations in parallel — each writes directly to its final destination
PIDS=()
for N in "${Ns[@]}"; do
    for (( i=1; i<=RUNS; i++ )); do
        DIR="$ROOT_DIR/runs/N${N}/run_${i}"
        mkdir -p "$DIR"
        (
            cd "$ROOT_DIR/java"
            java -server Main "$N" "$DIR/output.txt"
            echo "  done N=$N run=$i"
        ) &
        PIDS+=($!)
    done
done

# Wait for all and check exit codes
FAILED=0
for pid in "${PIDS[@]}"; do
    wait "$pid" || { echo "A run failed (pid $pid)"; FAILED=1; }
done
if [[ $FAILED -ne 0 ]]; then
    echo "One or more simulations failed." >&2
    exit 1
fi

NS_STRING=$(IFS=, ; echo "${ALL_NS[*]}")

python3 python/analyze.py \
    "${ALL_FILES[@]}" \
    --Ns "$NS_STRING" \
    --out analisis_total
