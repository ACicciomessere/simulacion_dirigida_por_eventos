import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import os

# --- parámetros del sistema ---
R_OUTER = 80
R_INNER = 1
CX = 0
CY = 0

# --- archivo ---
base_dir = os.path.dirname(__file__)
path = os.path.join(base_dir, "../outputs/sim_circular/output.txt")

# --- leer archivo ---
def load_data(filename):
    frames = []
    with open(filename, 'r') as f:
        lines = f.readlines()

    i = 0
    while i < len(lines):
        if lines[i].startswith("#"):
            i += 1
            continue

        t = float(lines[i].strip())
        i += 1

        particles = []
        while i < len(lines):
            if lines[i].strip() == "":
                i += 1
                continue

            parts = lines[i].split()
            if len(parts) != 4:
                break

            x, y, vx, vy = map(float, parts)
            particles.append([x, y, vx, vy])
            i += 1

        frames.append(np.array(particles))

    return frames

# --- cargar datos ---
print("Cargando datos...")
frames = load_data(path)
print("Frames originales:", len(frames))

# 🔥 CLAVE: reducir frames (evita que se trabe)
frames = frames[::10]
print("Frames usados:", len(frames))

N = len(frames[0])

# --- setup gráfico ---
fig, ax = plt.subplots()
ax.set_aspect('equal')
ax.set_xlim(CX - R_OUTER, CX + R_OUTER)
ax.set_ylim(CY - R_OUTER, CY + R_OUTER)

# círculos
outer_circle = plt.Circle((CX, CY), R_OUTER, fill=False)
inner_circle = plt.Circle((CX, CY), R_INNER, fill=False)
ax.add_patch(outer_circle)
ax.add_patch(inner_circle)

# partículas (más chicas = más rápido)
scat = ax.scatter([], [], s=3)

# --- estado ---
prev_frame = None
collision_timer = np.zeros(N)

# --- update ---
def update(frame):
    global prev_frame, collision_timer

    xs = frame[:, 0]
    ys = frame[:, 1]

     # default: aparecen frescas
    colors = np.array(['black'] * len(frame))

    if prev_frame is not None:
        dv = np.sqrt((frame[:,2] - prev_frame[:,2])**2 +
                     (frame[:,3] - prev_frame[:,3])**2)

        collisions = dv > 1e-4
    #     collision_timer[collisions] = 100

    # collision_timer[:] = np.maximum(collision_timer - 1, 0)

    # colors = np.where(collision_timer > 0, 'red', 'blue')
        for i in np.where(collisions)[0]:
            dx = xs[i] - CX
            dy = ys[i] - CY
            dist = np.sqrt(dx*dx + dy*dy)

            # 🟠 choque con círculo interno
            if abs(dist - R_INNER) < 2:
                colors[i] = 'red'

            # 🟢 choque con pared externa
            elif abs(dist - R_OUTER) < 2:
                colors[i] = 'green'

            else:
                colors[i] = 'red'  # fallback (por si acaso)

    scat.set_offsets(np.column_stack((xs, ys)))
    scat.set_color(colors)

    prev_frame = frame.copy()
    return scat,


# --- animación ---
ani = animation.FuncAnimation(
    fig,
    update,
    frames=frames,
    interval=500,
    blit=True
)

print("Guardando video...")

ani.save(
    "simulacion.mp4",
    fps=20,
    dpi=150
)

plt.show()