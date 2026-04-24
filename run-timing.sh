#!/bin/bash

# ── Configuración ─────────────────────────────────────────────
Ns=(200 300 400 500 600)
TIMING_RUNS=3

# ── Compilar Java ─────────────────────────────────────────────
cd java
javac *.java
if [ $? -ne 0 ]; then
    echo "Error: compilación Java fallida."
    exit 1
fi
cd ..

mkdir -p runs

# ── Timing: TIMING_RUNS mediciones por N ──────────────────────
TIMING_FILE="runs/timing.txt"
echo "# N elapsed_ms" > "$TIMING_FILE"

for N in "${Ns[@]}"; do
    echo "Timing N=$N  (${TIMING_RUNS} repeticiones) ..."
    for (( r=1; r<=TIMING_RUNS; r++ )); do
        cd java
        START_MS=$(date +%s%3N)
        java -Xmx2g Main $N timing
        END_MS=$(date +%s%3N)
        cd ..
        ELAPSED=$(( END_MS - START_MS ))
        echo "$N $ELAPSED" >> "$TIMING_FILE"
        echo "  run $r/${TIMING_RUNS} → ${ELAPSED} ms"
    done
done

echo ""
echo "Timing completo. Resultados en $TIMING_FILE"
echo "Generando gráfico..."

python3 python/analyze_timing.py \
    --timing_file runs/timing.txt \
    --out timing_analysis