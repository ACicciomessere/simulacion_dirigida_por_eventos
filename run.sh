#!/bin/bash
Ns=(300 400 500 600 700)
RUNS=3

cd java
javac *.java
cd ..

ALL_FILES=()
ALL_NS=()

# ── 1.1: timing ──────────────────────────────────────────────
# TIMING_FILE="runs/timing.txt"
# mkdir -p runs
# echo "# N elapsed_ms" > "$TIMING_FILE"

# for N in "${Ns[@]}"; do
#     echo "Timing N=$N ..."
#     cd java
#     START_MS=$(date +%s%3N)
#     java Main $N timing        # modo timing: corre tf=5s sin escribir output
#     END_MS=$(date +%s%3N)
#     cd ..
#     ELAPSED=$(( END_MS - START_MS ))
#     echo "$N $ELAPSED" >> "$TIMING_FILE"
#     echo "  → ${ELAPSED} ms"
# done

# ── 1.2/1.3/1.4: múltiples realizaciones ────────────────────
for N in "${Ns[@]}"; do
    for (( i=1; i<=RUNS; i++ )); do
        DIR="runs/N${N}/run_${i}"
        mkdir -p "$DIR"
        echo "N=$N run=$i"
        cd java
        java Main $N
        cd ..
        cp java/output.txt "$DIR/output.txt"
        ALL_FILES+=("$DIR/output.txt")
        ALL_NS+=("$N")
    done
done

NS_STRING=$(IFS=, ; echo "${ALL_NS[*]}")

python3 python/analyze.py \
    "${ALL_FILES[@]}" \
    --Ns "$NS_STRING" \
    --out analisis_total