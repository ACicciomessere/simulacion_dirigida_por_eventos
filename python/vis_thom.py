import matplotlib
matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import sys
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
            frame_data = []

            for _ in range(N):
                parts = f.readline().strip().split()
                if len(parts) >= 7:
                    x = float(parts[1])
                    y = float(parts[2])
                    vx = float(parts[3])
                    vy = float(parts[4])
                    radius = float(parts[5])
                    is_leader = int(parts[6])
                    frame_data.append([x, y, vx, vy, radius, is_leader])

            if frame_data:
                frames.append(np.array(frame_data))

    return frames


def main():
    # Usage: python vis_thom.py <L> <r0> [filename] [output_gif]
    if len(sys.argv) < 3:
        print("Usage: python vis_thom.py <L> <r0> [filename] [output_gif]")
        sys.exit(1)

    r0 = float(sys.argv[1])
    filename = sys.argv[2] if len(sys.argv) > 3 else "particles_frames.txt"
    output_gif = sys.argv[3] if len(sys.argv) > 4 else "output.gif"
    r_inner = 1.0  

    frames = parse_frames(filename)
    if not frames:
        print(f"No frames found in {filename}.")
        sys.exit(1)

    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_xlim(0, L)
    ax.set_ylim(0, L)
    ax.set_aspect('equal')

    circle_boundary = plt.Circle((50, 50), r0, color='g', fill=False, linestyle='-', linewidth=2)
    ax.add_patch(circle_boundary)

    inner_circle_boundary = plt.Circle((50, 50), r_inner, color='r', fill=True, linestyle='--', linewidth=1)
    ax.add_patch(inner_circle_boundary)

    init_data = frames[0]
    colors = ['orange' if c == 1 else 'blue' for c in init_data[:, 5]]
    scatter = ax.scatter(init_data[:, 0], init_data[:, 1], s=30, color=colors, zorder=2)
    quiver = ax.quiver(init_data[:, 0], init_data[:, 1], np.zeros_like(init_data[:, 0]), np.zeros_like(init_data[:, 1]),
                       color=colors, scale=1, scale_units='xy', angles='xy', headwidth=3, headlength=4, zorder=1)

    def init():
        return scatter

    def update(frame_idx):
        data = frames[frame_idx]
        x = data[:, 0]
        y = data[:, 1]
        vx = data[:, 2]
        vy = data[:, 3]
        is_leader = data[:, 5]

        offsets = np.c_[x, y]
        scatter.set_offsets(offsets)
        scatter.set_color(['orange' if c == 1 else 'blue' for c in is_leader])
        quiver.set_offsets(offsets)

        norm = np.sqrt(vx ** 2 + vy ** 2)
        norm[norm == 0] = 1
        ux = vx / norm
        uy = vy / norm
        quiver.set_color(['orange' if c == 1 else 'blue' for c in is_leader])

        ax.set_title(f"Simulation Frame: {frame_idx}")
        return scatter

    ani = animation.FuncAnimation(fig, update, frames=len(frames), init_func=init, blit=False, interval=100)
    plt.tight_layout()
    ani.save(output_gif, writer='pillow', fps=10)
    print(f"Saved to {output_gif}")


if __name__ == '__main__':
    main()
