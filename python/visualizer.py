import sys
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

def parse_frames(filename):
    frames = []
    if not os.path.exists(filename):
        print(f"Error: File {filename} not found.")
        sys.exit(1)
        
    with open(filename, 'r') as f:
        while True:
            line_N = f.readline()
            if not line_N:
                break
            
            try:
                N = int(line_N.strip())
            except ValueError:
                continue
                
            frame_label = f.readline().strip()
            
            x = []
            y = []
            vx = []
            vy = []
            radius = []
            is_leader = []
            
            for _ in range(N):
                parts = f.readline().strip().split()
                if len(parts) >= 7:
                    x.append(float(parts[1]))
                    y.append(float(parts[2]))
                    vx.append(float(parts[3]))
                    vy.append(float(parts[4]))
                    radius.append(float(parts[5]))
                    is_leader.append(int(parts[6]))
            
            frames.append({
                'x': np.array(x),
                'y': np.array(y),
                'vx': np.array(vx),
                'vy': np.array(vy),
                'radius': np.array(radius),
                'is_leader': np.array(is_leader)
            })
            
    return frames

def main():
    if len(sys.argv) < 3:
        print("Usage: python visualizer.py <L> <r0> [filename]")
        sys.exit(1)
        
    r0 = float(sys.argv[1])
    filename = sys.argv[2] if len(sys.argv) > 3 else "../particles_frames.txt"
    r_inner= 1.0
    
    frames = parse_frames(filename)
    if not frames:
        print(f"No frames found in {filename}.")
        sys.exit(1)
        
    fig, ax = plt.subplots(figsize=(8, 8))
    
    # Setup axes limits
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_aspect('equal')
    
    # Draw circular boundary
    circle_boundary = plt.Circle((50, 50), r0, color='g', fill=False, linestyle='-', linewidth=2)
    ax.add_patch(circle_boundary)

    #Draw inner circular boundary
    inner_circle_boundary = plt.Circle((50, 50), r_inner, color='r', fill=True, linestyle='--', linewidth=1)
    ax.add_patch(inner_circle_boundary)
    
    # Use quiver to draw arrows
    quiver = ax.quiver([], [], [], [], color=[], scale=1, scale_units='xy', angles='xy')
    
    def init():
        return quiver,
        
    def update(frame_idx):
        frame = frames[frame_idx]
        
        # Colors: leader is orange, normal are blue
        colors = ['orange' if leader == 1 else 'blue' for leader in frame['is_leader']]
        
        quiver.set_offsets(np.c_[frame['x'], frame['y']])
        
        # Normalize and scale vectors for visualization
        norm = np.sqrt(frame['vx']**2 + frame['vy']**2)
        # Avoid division by zero
        norm[norm == 0] = 1 
        
        ux = frame['vx'] / norm
        uy = frame['vy'] / norm

        quiver.set_color(colors)
        
        ax.set_title(f"Simulation Frame: {frame_idx}")
        return quiver,
        
    ani = animation.FuncAnimation(fig, update, frames=len(frames), init_func=init, blit=False, interval=50)
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
