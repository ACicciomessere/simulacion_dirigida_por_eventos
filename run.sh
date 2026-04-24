# !/bin/bash

# # --- Configuration Variables ---
# N=300                  # Number of particles
# M=10                   # Number of cell matrix divisions (L/M > rc)
# RC=1.0                 # Interaction radius
# PERIODIC="false"       # Whether boundaries should wrap around the box
# ITERATIONS=500         # Number of frames to simulate
# CIRCLE_LEADER="false"  # Whether to spawn a leader particle (requires orbit implementation)
# R0=40.0                # Inner circular swarm boundary radius limitation

# echo "--> Compiling Java code..."
# # Compile all source files inside java directory
# javac java/*.java

# if [ $? -ne 0 ]; then
#     echo "Error: Java compilation failed."
#     exit 1
# fi

# echo "--> Cleaning up old output..."
# rm -f particles_frames.txt

# echo "--> Running Simulation with R0=$R0"
# cd java
# # Run the application with configuration args
# # java App $N $M $RC $PERIODIC $ITERATIONS $CIRCLE_LEADER $R0
# java Main
# cd ..

# echo "--> Launching Python Visualizer..."
# # Make sure python depends are met (matplotlib, numpy)
# if [ "$1" == "T" ]; then
#     rm -rf output.gif
#     python3 python/vis_thom.py $R0 particles_frames.txt output.gif
# else
#     python3 python/visualize.py $R0 particles_frames.txt
# fi

#!/bin/bash
Ns=(200 300 400 500 600)
RUNS=3

cd java
javac *.java
cd ..

ALL_FILES=()
ALL_NS=()

# ── 1.1: timing ──────────────────────────────────────────────
TIMING_FILE="runs/timing.txt"
mkdir -p runs
echo "# N elapsed_ms" > "$TIMING_FILE"

for N in "${Ns[@]}"; do
    echo "Timing N=$N ..."
    cd java
    START_MS=$(python3 -c 'import time; print(int(time.time() * 1000))')
    java -Xmx2g Main $N timing        # modo timing: corre tf=5s sin escribir output
    END_MS=$(python3 -c 'import time; print(int(time.time() * 1000))')
    cd ..
    ELAPSED=$(( END_MS - START_MS ))
    echo "$N $ELAPSED" >> "$TIMING_FILE"
    echo "  → ${ELAPSED} ms"
done

# ── 1.2/1.3/1.4: múltiples realizaciones ────────────────────
for N in "${Ns[@]}"; do
    for (( i=1; i<=RUNS; i++ )); do
        DIR="runs/N${N}/run_${i}"
        mkdir -p "$DIR"
        echo "N=$N run=$i"
        cd java
        java -Xmx2g Main $N
        cd ..
        cp outputs/sim_circular/output.txt "$DIR/output.txt"
        ALL_FILES+=("$DIR/output.txt")
        ALL_NS+=("$N")
    done
done

NS_STRING=$(IFS=, ; echo "${ALL_NS[*]}")

python3 python/analyze.py \
    "${ALL_FILES[@]}" \
    --Ns "$NS_STRING" \
    --out analisis_total

# ── Análisis de timing (escala en función de N) ────────────────
echo ""
echo "Analizando timing (tiempo de ejecución vs N)..."
python3 python/analyze_timing.py --timing_file runs/timing.txt --out timing_analysis