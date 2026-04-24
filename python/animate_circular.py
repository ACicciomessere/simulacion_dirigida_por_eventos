"""
animate_circular.py
====================
Animación del Sistema 1: Scanning rate en recinto cerrado con obstáculo fijo.

Lee el archivo de salida generado por la simulación Java (Collisions.java)
y produce una animación con matplotlib.

Formato esperado del archivo de salida (outputs/sim_circular/output.txt):
    # t
    # x y vx vy
    <tiempo>
    x1 y1 vx1 vy1
    x2 y2 vx2 vy2
    ...

Uso:
    python animate_circular.py                         # usa ruta por defecto
    python animate_circular.py ruta/al/output.txt     # ruta personalizada
    python animate_circular.py --save video.mp4       # guarda video
    python animate_circular.py --save anim.gif        # guarda GIF

Lógica de colores (estados de partículas):
    - FRESCA  (verde)  : estado inicial; también tras rebotar en pared externa
    - USADA   (violeta): tras colisionar con el obstáculo central
"""

import sys
import os
import argparse
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib.collections import PathCollection

# ── Parámetros físicos del sistema (deben coincidir con la simulación) ──────
R_OUTER   = 40.0   # radio recinto externo  (L=80 → r=40)
R_INNER   = 0.5    # radio obstáculo central
PARTICLE_R = 1.0   # radio de cada partícula
CX, CY    = 0.0, 0.0

# ── Colores ──────────────────────────────────────────────────────────────────
COLOR_FRESH  = "#2ecc71"   # verde
COLOR_USED   = "#9b59b6"   # violeta
COLOR_EDGE   = "#2ecc71"    # borde recinto
COLOR_OBS    = "#c0392b"   # obstáculo central
BG_COLOR     = "#fcfcff"   # fondo oscuro
TEXT_COLOR   = "#1a1a2e"

# ─────────────────────────────────────────────────────────────────────────────

def parse_output(filepath: str):
    """
    Lee el archivo output.txt y devuelve una lista de frames.
    Cada frame es un dict:
        { 'time': float, 'x': np.array, 'y': np.array,
          'vx': np.array, 'vy': np.array }
    """
    frames = []
    current_time = None
    xs, ys, vxs, vys = [], [], [], []

    with open(filepath, "r") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split()

            # ¿Es una línea de tiempo? (sólo 1 valor)
            if len(parts) == 1:
                # guardar frame anterior si existe
                if current_time is not None and xs:
                    frames.append({
                        "time": current_time,
                        "x":  np.array(xs,  dtype=float),
                        "y":  np.array(ys,  dtype=float),
                        "vx": np.array(vxs, dtype=float),
                        "vy": np.array(vys, dtype=float),
                    })
                current_time = float(parts[0])
                xs, ys, vxs, vys = [], [], [], []

            # ¿Es una línea de partícula? (4 valores: x y vx vy)
            elif len(parts) == 4:
                xs.append(float(parts[0]))
                ys.append(float(parts[1]))
                vxs.append(float(parts[2]))
                vys.append(float(parts[3]))

    # último frame
    if current_time is not None and xs:
        frames.append({
            "time": current_time,
            "x":  np.array(xs,  dtype=float),
            "y":  np.array(ys,  dtype=float),
            "vx": np.array(vxs, dtype=float),
            "vy": np.array(vys, dtype=float),
        })

    return frames


def classify_particles(frames):
    """
    Determina el estado (fresca/usada) de cada partícula en cada frame.

    Reglas simplificadas:
      - Inicialmente TODAS son frescas (verde).
      - Si colisionó con el obstáculo central → cambia a usada (violeta).
      - Una vez usada, permanece usada (nunca vuelve a cambiar).

    Heurística: detectamos cambio de velocidad entre frames consecutivos.
    Si cambia velocidad Y está muy cerca del centro, es colisión con obstáculo.
    """
    if not frames:
        return []

    n_particles = len(frames[0]["x"])
    # estado[i] = True → usada (violeta), False → fresca (verde)
    state = np.zeros(n_particles, dtype=bool)
    all_states = []

    prev = frames[0]
    all_states.append(state.copy())

    for frame in frames[1:]:
        x, y   = frame["x"],  frame["y"]
        vx, vy = frame["vx"], frame["vy"]
        px, py   = prev["x"],  prev["y"]
        pvx, pvy = prev["vx"], prev["vy"]

        # velocidad cambió → hubo colisión
        speed_changed = (np.abs(vx - pvx) > 1e-9) | (np.abs(vy - pvy) > 1e-9)

        dist_center = np.sqrt(x**2 + y**2)

        # colisión con obstáculo central: velocidad cambió Y está muy cerca del centro
        hit_inner = speed_changed & (dist_center < (R_INNER + PARTICLE_R + 2.0))

        # aplicar cambio de estado: fresca → usada si choca con obstáculo
        # Una vez usada, se queda usada (nunca vuelve a cambiar)
        state[hit_inner & ~state] = True

        all_states.append(state.copy())
        prev = frame

    return all_states


def build_animation(frames, all_states, interval_ms=50):
    """Construye y devuelve la FuncAnimation."""

    fig, ax = plt.subplots(figsize=(8, 8), facecolor=BG_COLOR)
    ax.set_facecolor(BG_COLOR)
    ax.set_aspect("equal")
    margin = R_OUTER * 1.08
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin, margin)
    ax.axis("off")

    # ── Decoración estática ──────────────────────────────────────────────────
    # Recinto externo
    outer_ring = plt.Circle((CX, CY), R_OUTER,
                             fill=False, edgecolor=COLOR_EDGE,
                             linewidth=2.5, zorder=2)
    ax.add_patch(outer_ring)

    # Obstáculo central (fijo)
    inner_obs = plt.Circle((CX, CY), R_INNER,
                            fill=True, facecolor=COLOR_USED,
                            edgecolor="#954bb5" , linewidth=1.5, zorder=5)
    ax.add_patch(inner_obs)

    # Título y etiquetas
    title_text = ax.text(0, margin * 0.92, "Sistema 1 – Scanning Rate",
                         ha="center", va="center", fontsize=14,
                         color=TEXT_COLOR, fontweight="bold")
    time_text  = ax.text(-margin * 0.97, -margin * 0.97, "",
                         ha="left", va="bottom", fontsize=11,
                         color=TEXT_COLOR, family="monospace")
    info_text  = ax.text( margin * 0.97, -margin * 0.97, "",
                         ha="right", va="bottom", fontsize=10,
                         color=TEXT_COLOR, family="monospace")

    # Leyenda manual
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=COLOR_FRESH, markersize=10, label="Fresca"),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=COLOR_USED,  markersize=10, label="Usada"),
        Line2D([0], [0], marker='o', color='w',
               markerfacecolor=COLOR_OBS,   markersize=10, label="Obstáculo"),
    ]
    ax.legend(handles=legend_elements, loc="upper right",
              facecolor="#2c2c4a", edgecolor="#555577",
              labelcolor=TEXT_COLOR, fontsize=9)

    # ── Partículas (scatter) ─────────────────────────────────────────────────
    # tamaño visual: convertir radio físico → puntos^2
    # El radio en unidades de datos es PARTICLE_R; necesitamos el equivalente en pts
    def radius_to_pts2(r_data):
        """Convierte radio en unidades de datos a tamaño scatter (pts^2)."""
        fig_width_in  = fig.get_figwidth()
        ax_width_frac = ax.get_position().width
        ax_width_in   = fig_width_in * ax_width_frac
        data_range    = ax.get_xlim()[1] - ax.get_xlim()[0]
        dpi           = fig.dpi
        ax_width_pts  = ax_width_in * dpi
        pts_per_data  = ax_width_pts / data_range
        return (r_data * pts_per_data) ** 2

    pt_size = radius_to_pts2(PARTICLE_R) * 0.8  # pequeño margen visual

    frame0  = frames[0]
    state0  = all_states[0]
    colors0 = [COLOR_USED if s else COLOR_FRESH for s in state0]

    scat = ax.scatter(frame0["x"], frame0["y"],
                      s=pt_size, c=colors0,
                      zorder=4, linewidths=0.4,
                      edgecolors="white", alpha=0.92)

    # ── Función de actualización ─────────────────────────────────────────────
    def update(i):
        frame  = frames[i]
        state  = all_states[i]
        colors = [COLOR_USED if s else COLOR_FRESH for s in state]

        scat.set_offsets(np.column_stack([frame["x"], frame["y"]]))
        scat.set_facecolor(colors)

        n_used  = int(np.sum(state))
        n_fresh = len(state) - n_used
        time_text.set_text(f"t = {frame['time']:.3f} s")
        info_text.set_text(f"Frescas: {n_fresh}  Usadas: {n_used}")

        return scat, time_text, info_text

    anim = FuncAnimation(fig, update,
                         frames=len(frames),
                         interval=interval_ms,
                         blit=True)
    return fig, anim


# ─────────────────────────────────────────────────────────────────────────────
#  Estadísticas rápidas (opcional, se imprimen en consola)
# ─────────────────────────────────────────────────────────────────────────────

def print_stats(frames, all_states):
    n = len(frames[0]["x"])
    print(f"\n{'─'*50}")
    print(f"  Partículas : {n}")
    print(f"  Frames     : {len(frames)}")
    if frames:
        t0 = frames[0]["time"]
        tf = frames[-1]["time"]
        print(f"  Tiempo     : {t0:.3f} s → {tf:.3f} s")
    # fracción usada promedio en el último 20% de frames
    cutoff = max(1, int(len(all_states) * 0.8))
    frac = np.mean([np.mean(s) for s in all_states[cutoff:]])
    print(f"  F_est (≈)  : {frac:.3f}  (fracción usada en últimos 20% frames)")
    print(f"{'─'*50}\n")


# ─────────────────────────────────────────────────────────────────────────────
#  Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Animación Sistema 1 – Scanning Rate (recinto circular)")
    parser.add_argument(
        "output_file", nargs="?",
        default="outputs/sim_circular/output.txt",
        help="Ruta al archivo output.txt generado por la simulación Java")
    parser.add_argument(
        "--save", metavar="FILE",
        help="Guardar animación en archivo (.mp4 o .gif)")
    parser.add_argument(
        "--interval", type=int, default=25,
        help="Intervalo entre frames en ms (default: 25)")
    parser.add_argument(
        "--fps", type=int, default=120,
        help="FPS para el video guardado (default: 60)")
    parser.add_argument(
        "--skip", type=int, default=2,
        help="Mostrar 1 de cada N frames (default: 2)")
    args = parser.parse_args()

    # ── Cargar datos ────────────────────────────────────────────────────────
    if not os.path.isfile(args.output_file):
        print(f"[ERROR] No se encontró el archivo: {args.output_file}")
        print("  Asegurate de correr primero la simulación Java.")
        print("  El archivo debe estar en: outputs/sim_circular/output.txt")
        sys.exit(1)

    print(f"[INFO] Leyendo: {args.output_file} …")
    frames = parse_output(args.output_file)

    if not frames:
        print("[ERROR] El archivo no contiene frames válidos.")
        sys.exit(1)

    # Sub-muestreo de frames
    if args.skip > 1:
        frames = frames[::args.skip]

    # ── Clasificar estados ──────────────────────────────────────────────────
    print("[INFO] Clasificando estados de partículas …")
    all_states = classify_particles(frames)

    print_stats(frames, all_states)

    # ── Animación ───────────────────────────────────────────────────────────
    print("[INFO] Construyendo animación …")
    fig, anim = build_animation(frames, all_states, interval_ms=args.interval)

    if args.save:
        ext = os.path.splitext(args.save)[1].lower()
        print(f"[INFO] Guardando en '{args.save}' (fps={args.fps}) …")
        if ext == ".gif":
            writer = PillowWriter(fps=args.fps)
        else:
            writer = FFMpegWriter(fps=args.fps, bitrate=1800,
                                  extra_args=["-vcodec", "libx264"])
        anim.save(args.save, writer=writer, dpi=120)
        print("[INFO] ¡Animación guardada!")
    else:
        plt.tight_layout()
        plt.show()


if __name__ == "__main__":
    main()