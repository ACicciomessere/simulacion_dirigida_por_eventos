"""
analyze_timing.py  –  Punto 1.1: tiempo de ejecución vs N

Lee runs/timing.txt donde cada línea es  "N  elapsed_ms"
(puede haber varias líneas por N, una por repetición del timing).

Para cada N calcula:
    - media del tiempo (ms y s)
    - desvío estándar
    - error estándar de la media (SEM)

Grafica t_mean vs N en escala LINEAL con barras de error (±std).

Uso:
    python analyze_timing.py --timing_file runs/timing.txt --out timing_analysis
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict


# ──────────────────────────────────────────────────────────────
# 1. Lectura del archivo
# ──────────────────────────────────────────────────────────────

def load_timing(filepath):
    """
    Devuelve un dict  {N: [t1_ms, t2_ms, ...]}
    Ignora líneas que empiezan con '#' o están vacías.
    """
    data = defaultdict(list)
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            N  = int(parts[0])
            ms = float(parts[1])
            data[N].append(ms)
    return data


# ──────────────────────────────────────────────────────────────
# 2. Plot
# ──────────────────────────────────────────────────────────────

def plot_timing(data, out_prefix="timing_analysis"):
    """
    Grafica t_mean (s) vs N en escala lineal.
    Barras de error = desvío estándar de las repeticiones.
    Anota el número de repeticiones sobre cada punto.
    """
    Ns      = sorted(data.keys())
    means   = []    # segundos
    stds    = []    # segundos
    sems    = []    # error estándar de la media
    n_reps  = []

    for N in Ns:
        arr = np.array(data[N]) / 1000.0   # ms → s
        means.append(arr.mean())
        stds.append(arr.std(ddof=1) if len(arr) > 1 else 0.0)
        sems.append(arr.std(ddof=1) / np.sqrt(len(arr)) if len(arr) > 1 else 0.0)
        n_reps.append(len(arr))

    means = np.array(means)
    stds  = np.array(stds)
    sems  = np.array(sems)

    fig, ax = plt.subplots(figsize=(8, 5))
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')

    # Puntos con barras de error = std (variabilidad real de las mediciones)
    ax.errorbar(
        Ns, means, yerr=stds,
        fmt='o-',
        color='#457b9d',
        ecolor='#e63946',
        elinewidth=1.8,
        capsize=6,
        capthick=1.8,
        linewidth=2,
        markersize=8,
        markerfacecolor='white',
        markeredgewidth=2,
        label='Tiempo medio ± std'
    )

    # Anotar cada punto con su valor medio y nro de repeticiones
    for N, t, s, nr in zip(Ns, means, stds, n_reps):
        ax.annotate(
            f"{t:.2f} s\n(n={nr})",
            xy=(N, t),
            xytext=(0, 14),
            textcoords='offset points',
            ha='center',
            fontsize=8,
            color='#1d3557'
        )

    ax.set_xlabel("N  (número de partículas)", fontsize=12)
    ax.set_ylabel("Tiempo de ejecución  [s]",   fontsize=12)
    ax.set_title("1.1  Tiempo de ejecución vs N  (escala lineal)", fontsize=13)

    # Escala lineal explícita (no log)
    ax.set_xscale('linear')
    ax.set_yscale('linear')

    # Ajustar límites para que las anotaciones no queden cortadas
    ax.set_ylim(bottom=0, top=max(means) * 1.35)
    ax.set_xlim(Ns[0] - 20, Ns[-1] + 20)

    ax.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)
    ax.legend(fontsize=10, facecolor='white', edgecolor='#cccccc')

    for spine in ax.spines.values():
        spine.set_edgecolor('black')
    ax.tick_params(colors='black')

    plt.tight_layout()
    out_path = out_prefix + ".png"
    plt.savefig(out_path, dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    print(f"Figura guardada: {out_path}")

    # Imprimir tabla resumen en consola
    print("\n── Resumen timing ──────────────────────────────")
    print(f"{'N':>6}  {'n_reps':>6}  {'mean (s)':>10}  {'std (s)':>9}  {'SEM (s)':>9}")
    for N, t, s, sem, nr in zip(Ns, means, stds, sems, n_reps):
        print(f"{N:>6}  {nr:>6}  {t:>10.3f}  {s:>9.3f}  {sem:>9.4f}")

    plt.show()


# ──────────────────────────────────────────────────────────────
# 3. Entrypoint
# ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Análisis de timing: t vs N en escala lineal")
    parser.add_argument("--timing_file", default="runs/timing.txt",
                        help="Archivo con columnas N elapsed_ms")
    parser.add_argument("--out", default="timing_analysis",
                        help="Prefijo del archivo de salida")
    args = parser.parse_args()

    data = load_timing(args.timing_file)
    if not data:
        print(f"No se encontraron datos en {args.timing_file}")
        return

    plot_timing(data, out_prefix=args.out)


if __name__ == "__main__":
    main()