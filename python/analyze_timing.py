"""
analyze_timing.py – Análisis del tiempo de ejecución vs N para 500ms

Lee el archivo runs/timing.txt generado por run.sh y genera un gráfico mostrando
la escala de tiempo de ejecución con la cantidad de partículas N.

Uso:
    python analyze_timing.py [--timing_file runs/timing.txt] [--out timing_plot]
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats


def parse_timing_file(filepath):
    """
    Parsea archivo timing.txt con formato:
        # N elapsed_ms
        N1 time1
        N2 time2
        ...
    
    Devuelve:
        Ns: array de valores N
        times_ms: array de tiempos en ms
    """
    Ns, times_ms = [], []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) >= 2:
                try:
                    N = int(parts[0])
                    t = float(parts[1])
                    Ns.append(N)
                    times_ms.append(t)
                except ValueError:
                    continue
    return np.array(Ns), np.array(times_ms)


def plot_timing(Ns, times_ms, out_prefix="timing_plot"):
    """
    Grafica tiempo de ejecución vs N con regresión.
    """
    # Ajuste de potencia: t = a * N^b
    log_N = np.log(Ns)
    log_t = np.log(times_ms)
    slope, intercept, r_value, p_value, std_err = stats.linregress(log_N, log_t)
    
    # Parámetros: t = a * N^b
    b = slope
    a = np.exp(intercept)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.patch.set_facecolor('white')
    
    # ── Gráfico 1: escala lineal ──
    ax1 = axes[0]
    ax1.set_facecolor('white')
    ax1.plot(Ns, times_ms, 'o-', color='#58a6ff', lw=2.5, markersize=10, 
             label='Mediciones')
    
    # Curva de ajuste
    N_smooth = np.linspace(Ns.min(), Ns.max(), 100)
    t_smooth = a * (N_smooth ** b)
    ax1.plot(N_smooth, t_smooth, '--', color='#f0883e', lw=2, 
             label=f'Ajuste: t = {a:.3f}·N^{b:.2f}')
    
    ax1.set_xlabel('N (número de partículas)', fontsize=11, color='black')
    ax1.set_ylabel('Tiempo [ms]', fontsize=11, color='black')
    ax1.set_title('Tiempo de ejecución vs N (escala lineal)', fontsize=12, color='black')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=10, loc='upper left')
    ax1.tick_params(colors='black')
    for spine in ax1.spines.values():
        spine.set_edgecolor('black')
    
    # Anotaciones de puntos
    for N, t in zip(Ns, times_ms):
        ax1.annotate(f'{t:.0f}ms', xy=(N, t), xytext=(5, 5), 
                    textcoords='offset points', fontsize=9, color='black')
    
    # ── Gráfico 2: escala log-log ──
    ax2 = axes[1]
    ax2.set_facecolor('white')
    ax2.loglog(Ns, times_ms, 'o-', color='#58a6ff', lw=2.5, markersize=10, 
               label='Mediciones')
    ax2.loglog(N_smooth, t_smooth, '--', color='#f0883e', lw=2, 
               label=f'Ajuste: log(t) = {intercept:.3f} + {b:.2f}·log(N)')
    
    ax2.set_xlabel('N (log scale)', fontsize=11, color='black')
    ax2.set_ylabel('Tiempo [ms] (log scale)', fontsize=11, color='black')
    ax2.set_title(f'Tiempo de ejecución vs N (escala log-log) — R²={r_value**2:.4f}', 
                  fontsize=12, color='black')
    ax2.grid(True, alpha=0.3, which='both')
    ax2.legend(fontsize=10, loc='upper left')
    ax2.tick_params(colors='black')
    for spine in ax2.spines.values():
        spine.set_edgecolor('black')
    
    fig.suptitle('Análisis de escalabilidad: tiempo de ejecución (primera parte, 500ms simulación)', 
                 fontsize=13, color='black', y=0.98)
    plt.tight_layout()
    
    plt.savefig(out_prefix + ".png", dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    print(f"Figura guardada: {out_prefix}.png")
    
    # Imprimir estadísticas
    print("\n" + "="*60)
    print("ANÁLISIS DE TIMING (500ms de simulación)")
    print("="*60)
    print(f"Puntos medidos: {len(Ns)}")
    print(f"Rango N: [{Ns.min()}, {Ns.max()}]")
    print(f"Rango tiempo: [{times_ms.min():.0f}, {times_ms.max():.0f}] ms")
    print()
    print(f"Ajuste potencial: t = {a:.6f} × N^{b:.4f}")
    print(f"R² = {r_value**2:.6f}")
    print(f"Error estándar: {std_err:.6f}")
    print(f"p-value: {p_value:.2e}")
    print()
    print(f"Complejidad: O(N^{b:.2f})")
    
    # Predicciones
    print("\nPredicciones para otros valores de N:")
    for N_pred in [100, 400, 500, 600, 700, 1000]:
        t_pred = a * (N_pred ** b)
        print(f"  N={N_pred:4d} → {t_pred:8.1f} ms = {t_pred/1000:6.2f} s")
    
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Análisis de timing — tiempo de ejecución vs N")
    parser.add_argument("--timing_file", type=str, default="runs/timing.txt",
                       help="Archivo de timing a procesar")
    parser.add_argument("--out", type=str, default="timing_plot",
                       help="Prefijo para guardar figuras")
    args = parser.parse_args()
    
    try:
        Ns, times_ms = parse_timing_file(args.timing_file)
        if len(Ns) == 0:
            print(f"Error: No se encontraron datos en {args.timing_file}")
            sys.exit(1)
        plot_timing(Ns, times_ms, out_prefix=args.out)
    except FileNotFoundError:
        print(f"Error: Archivo {args.timing_file} no encontrado")
        sys.exit(1)


if __name__ == "__main__":
    main()
