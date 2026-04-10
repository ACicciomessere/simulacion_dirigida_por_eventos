#!/bin/bash

# --- Configuration Variables ---
N=300                  # Number of particles
L=20                   # Bounding box size (L x L)
M=10                   # Number of cell matrix divisions (L/M > rc)
RC=1.0                 # Interaction radius
PERIODIC="false"       # Whether boundaries should wrap around the box
ITERATIONS=500         # Number of frames to simulate
ETA=0.5                # Noise factor
CIRCLE_LEADER="false"  # Whether to spawn a leader particle (requires orbit implementation)
R0=8.0                 # Inner circular swarm boundary radius limitation

echo "--> Compiling Java code..."
# Compile all source files inside java directory
javac java/*.java

if [ $? -ne 0 ]; then
    echo "Error: Java compilation failed."
    exit 1
fi

echo "--> Cleaning up old output..."
rm -f particles_frames.txt

echo "--> Running Simulation with R0=$R0 and L=$L..."
cd java
# Run the application with configuration args
java App $N $L $M $RC $PERIODIC $ITERATIONS $ETA $CIRCLE_LEADER $R0
cd ..

echo "--> Launching Python Visualizer..."
# Make sure python depends are met (matplotlib, numpy)
python3 python/visualizer.py $L $R0 particles_frames.txt
