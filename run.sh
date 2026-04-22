#!/bin/bash

# --- Configuration Variables ---
N=300                  # Number of particles
M=10                   # Number of cell matrix divisions (L/M > rc)
RC=1.0                 # Interaction radius
PERIODIC="false"       # Whether boundaries should wrap around the box
ITERATIONS=500         # Number of frames to simulate
CIRCLE_LEADER="false"  # Whether to spawn a leader particle (requires orbit implementation)
R0=40.0                # Inner circular swarm boundary radius limitation

echo "--> Compiling Java code..."
# Compile all source files inside java directory
javac java/*.java

if [ $? -ne 0 ]; then
    echo "Error: Java compilation failed."
    exit 1
fi

echo "--> Cleaning up old output..."
rm -f particles_frames.txt

echo "--> Running Simulation with R0=$R0"
cd java
# Run the application with configuration args
# java App $N $M $RC $PERIODIC $ITERATIONS $CIRCLE_LEADER $R0
java Main
cd ..

echo "--> Launching Python Visualizer..."
# Make sure python depends are met (matplotlib, numpy)
if [ "$1" == "T" ]; then
    rm -rf output.gif
    python3 python/vis_thom.py $R0 particles_frames.txt output.gif
else
    python3 python/visualize.py $R0 particles_frames.txt
fi