"""
analyze.py  –  Puntos 1.2, 1.3 y 1.4 del Sistema 1 (Scanning rate en recinto circular)

Lee el archivo output.txt generado por la simulación Java original (sin modificar).
Formato esperado:
    # t
    # x y vx vy
    <tiempo>
    x1 y1 vx1 vy1
    x2 y2 vx2 vy2
    ...

Uso:
    python analyze.py output.txt [--N 200] [--r_outer 40] [--r_inner 1] [--radius 1]

Si hay múltiples realizaciones, pasar varios archivos:
    python analyze.py run1/output.txt run2/output.txt run3/output.txt
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

# ──────────────────────────────────────────────────────────────
# 1. Parser del output.txt original
# ──────────────────────────────────────────────────────────────

def parse_output(filepath):
    """
    Devuelve:
        times  : array (T,)            tiempos de los snapshots
        states : array (T, N, 4)       x, y, vx, vy  para cada snapshot
    """
    times, frames = [], []
    with open(filepath) as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]

    i = 0
    while i < len(lines):
        try:
            t = float(lines[i])
        except ValueError:
            i += 1
            continue
        i += 1
        frame = []
        while i < len(lines):
            parts = lines[i].split()
            if len(parts) == 4:
                frame.append([float(v) for v in parts])
                i += 1
            else:
                break
        if frame:
            times.append(t)
            frames.append(frame)

    times  = np.array(times)
    states = np.array(frames)   # (T, N, 4)
    return times, states


# ──────────────────────────────────────────────────────────────
# 2. Reconstruir estados fresca/usada y detectar colisiones
# ──────────────────────────────────────────────────────────────

def reconstruct_states(times, states, r_outer, r_inner, particle_radius, tol=0.5):
    """
    Para cada snapshot, infiere si cada partícula es 'fresca' (0) o 'usada' (1).

    Reglas (igual que la simulación):
      - Arranca todo el mundo fresca.
      - Si en cualquier snapshot una partícula está muy cerca del obstáculo interno
        (dist_to_center ≈ r_inner + particle_radius), marcamos ese evento como
        colisión con el centro → la partícula pasa a 'usada'.
      - Si una partícula usada llega muy cerca de la pared externa
        (dist_to_center ≈ r_outer - particle_radius), vuelve a ser 'fresca'.

    Devuelve:
        particle_state  : array (T, N)   0 = fresca, 1 = usada  en cada instante
        cfc_series      : array (T,)     conteo acumulado de transiciones fresh→used
        hit_center_mask : array (T, N)   True si en ese snapshot hubo colisión con centro
        hit_outer_mask  : array (T, N)   True si en ese snapshot hubo colisión con pared
    """
    T, N, _ = states.shape
    x  = states[:, :, 0]   # (T, N)
    y  = states[:, :, 1]
    vx = states[:, :, 2]
    vy = states[:, :, 3]

    dist = np.sqrt(x**2 + y**2)   # distancia al centro en cada snapshot

    # Detectar snapshots donde la partícula está "en contacto" con cada pared
    contact_inner = np.abs(dist - (r_inner + particle_radius)) < tol
    contact_outer = np.abs(dist - (r_outer - particle_radius)) < tol

    # Reconstruir estado snapshot a snapshot
    particle_state = np.zeros((T, N), dtype=int)  # 0=fresca, 1=usada
    current_state  = np.zeros(N, dtype=int)        # estado actual de cada partícula
    cfc = 0
    cfc_series = np.zeros(T)

    for t in range(T):
        for j in range(N):
            if contact_inner[t, j] and current_state[j] == 0:
                current_state[j] = 1   # fresca → usada
                cfc += 1
            elif contact_outer[t, j] and current_state[j] == 1:
                current_state[j] = 0   # usada → fresca
        particle_state[t] = current_state.copy()
        cfc_series[t] = cfc

    return particle_state, cfc_series, contact_inner, contact_outer


# ──────────────────────────────────────────────────────────────
# 3. Observables
# ──────────────────────────────────────────────────────────────

def compute_J(times, cfc_series):
    """Regresión lineal de Cfc(t) → pendiente = J (scanning rate)."""
    slope, intercept, r, p, se = stats.linregress(times, cfc_series)
    return slope, se

def compute_Fu(particle_state):
    """Fu(t) = fracción de partículas usadas en cada snapshot."""
    return particle_state.mean(axis=1)

def compute_radial_profiles(times, states, particle_state, r_inner, r_outer,
                             dS=0.2):
    """
    Para cada snapshot: selecciona partículas frescas con velocidad radial
    apuntando al centro (Rj·vj < 0), las agrupa por capa S, calcula
    densidad media y velocidad radial media.

    Devuelve arrays indexados por shell:
        S_centers, avg_density, avg_velocity, flux
    """
    num_shells = int(np.ceil((r_outer - r_inner) / dS))
    shell_edges = r_inner + np.arange(num_shells + 1) * dS

    sum_density  = np.zeros(num_shells)
    sum_velocity = np.zeros(num_shells)
    count_snaps  = np.zeros(num_shells, dtype=int)  # snapshots con datos en esa capa

    x  = states[:, :, 0]
    y  = states[:, :, 1]
    vx = states[:, :, 2]
    vy = states[:, :, 3]
    dist = np.sqrt(x**2 + y**2)

    for t in range(len(times)):
        fresh_mask = (particle_state[t] == 0)
        radial_vel = (x[t] * vx[t] + y[t] * vy[t]) / (dist[t] + 1e-15)
        toward_center = (radial_vel < 0) & fresh_mask

        for k in range(num_shells):
            in_shell = toward_center & (dist[t] >= shell_edges[k]) & (dist[t] < shell_edges[k+1])
            n_in = np.sum(in_shell)
            if n_in == 0:
                continue
            area = np.pi * (shell_edges[k+1]**2 - shell_edges[k]**2)
            sum_density[k]  += n_in / area
            sum_velocity[k] += np.mean(np.abs(radial_vel[in_shell]))
            count_snaps[k]  += 1

    valid = count_snaps > 0
    avg_density  = np.where(valid, sum_density  / np.where(count_snaps>0, count_snaps, 1), 0.0)
    avg_velocity = np.where(valid, sum_velocity / np.where(count_snaps>0, count_snaps, 1), 0.0)
    flux = avg_density * avg_velocity

    S_centers = r_inner + (np.arange(num_shells) + 0.5) * dS
    return S_centers, avg_density, avg_velocity, flux


# ──────────────────────────────────────────────────────────────
# 4. Plots
# ──────────────────────────────────────────────────────────────

COLORS = ['#e63946', '#457b9d', '#2a9d8f', '#e9c46a', '#f4a261']

def plot_all(results, r_outer=40, r_inner=1, particle_radius=1, out_prefix="plot"):
    """
    results: list of dicts, one per (N, realization):
        { 'N': int, 'times': ..., 'cfc_series': ..., 'Fu': ...,
          'S': ..., 'density': ..., 'velocity': ..., 'flux': ... }
    """
    Ns = sorted(set(r['N'] for r in results))

    # ── Agrupar por N ──
    by_N = {N: [r for r in results if r['N'] == N] for N in Ns}

    fig = plt.figure(figsize=(18, 14))
    fig.patch.set_facecolor('#0d1117')
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

    ax_cfc  = fig.add_subplot(gs[0, 0])   # 1.2a: Cfc(t)
    ax_J    = fig.add_subplot(gs[0, 1])   # 1.2b: <J>(N)
    ax_Fu   = fig.add_subplot(gs[1, 0])   # 1.3a: Fu(t)
    ax_Fest = fig.add_subplot(gs[1, 1])   # 1.3b: Fest vs N
    ax_Tss  = fig.add_subplot(gs[1, 2])   # 1.3c: T_ss vs N
    ax_prof = fig.add_subplot(gs[2, :2])  # 1.4a: perfiles radiales (N fijo)
    ax_Jin  = fig.add_subplot(gs[2, 2])   # 1.4b: Jin @ S≈2 vs N

    for ax in [ax_cfc, ax_J, ax_Fu, ax_Fest, ax_Tss, ax_prof, ax_Jin]:
        ax.set_facecolor('#161b22')
        ax.tick_params(colors='#c9d1d9')
        ax.xaxis.label.set_color('#c9d1d9')
        ax.yaxis.label.set_color('#c9d1d9')
        ax.title.set_color('#e6edf3')
        for spine in ax.spines.values():
            spine.set_edgecolor('#30363d')

    # ── 1.2a: Cfc(t) para cada N (primera realización) ──
    ax_cfc.set_title("1.2  Cfc(t) — conteo acumulado fresca→usada")
    ax_cfc.set_xlabel("t  [s]"); ax_cfc.set_ylabel("Cfc(t)")
    for i, N in enumerate(Ns):
        r0 = by_N[N][0]
        ax_cfc.plot(r0['times'], r0['cfc_series'], color=COLORS[i % len(COLORS)],
                    lw=1.5, label=f"N={N}")
        # línea de regresión
        J, _ = compute_J(r0['times'], r0['cfc_series'])
        t_line = np.linspace(r0['times'][0], r0['times'][-1], 100)
        ax_cfc.plot(t_line, J * t_line, color=COLORS[i % len(COLORS)],
                    ls='--', lw=1, alpha=0.6)
    ax_cfc.legend(fontsize=8, labelcolor='#c9d1d9', facecolor='#161b22', edgecolor='#30363d')

    # ── 1.2b: <J>(N) con barra de error ──
    ax_J.set_title("1.2  Scanning rate  ⟨J⟩  vs  N")
    ax_J.set_xlabel("N"); ax_J.set_ylabel("J  [colisiones/s]")
    J_means, J_stds, N_vals = [], [], []
    for N in Ns:
        Js = [compute_J(r['times'], r['cfc_series'])[0] for r in by_N[N]]
        J_means.append(np.mean(Js)); J_stds.append(np.std(Js)); N_vals.append(N)
    ax_J.errorbar(N_vals, J_means, yerr=J_stds, fmt='o-',
                  color='#58a6ff', ecolor='#f0883e', capsize=5, lw=2)

    # ── 1.3a: Fu(t) ──
    ax_Fu.set_title("1.3  Fu(t) — fracción de partículas usadas")
    ax_Fu.set_xlabel("t  [s]"); ax_Fu.set_ylabel("Fu(t) = Nu/N")
    for i, N in enumerate(Ns):
        for j, r in enumerate(by_N[N]):
            ax_Fu.plot(r['times'], r['Fu'], color=COLORS[i % len(COLORS)],
                       lw=1.2, alpha=0.5 if j > 0 else 1.0,
                       label=f"N={N}" if j == 0 else None)
    ax_Fu.legend(fontsize=8, labelcolor='#c9d1d9', facecolor='#161b22', edgecolor='#30363d')
    ax_Fu.set_ylim(-0.05, 1.05)

    # ── 1.3b: Fest vs N ──
    ax_Fest.set_title("1.3  Fest (valor estacionario) vs N")
    ax_Fest.set_xlabel("N"); ax_Fest.set_ylabel("Fest")
    Fest_vals = []
    for N in Ns:
        fests = []
        for r in by_N[N]:
            fu = r['Fu']
            # estacionario = promedio del último 30% del tiempo
            n_tail = max(1, int(len(fu) * 0.3))
            fests.append(np.mean(fu[-n_tail:]))
        Fest_vals.append(np.mean(fests))
    ax_Fest.plot(N_vals, Fest_vals, 'o-', color='#3fb950', lw=2)

    # ── 1.3c: T_ss vs N (tiempo al estacionario) ──
    ax_Tss.set_title("1.3  T_estacionario vs N")
    ax_Tss.set_xlabel("N"); ax_Tss.set_ylabel("T_ss  [s]")
    Tss_vals = []
    for N in Ns:
        tsss = []
        for r in by_N[N]:
            fu   = r['Fu']
            t    = r['times']
            fest = np.mean(fu[int(len(fu)*0.7):])
            # primer momento donde Fu supera 90% de Fest
            threshold = 0.9 * fest
            idx = np.where(fu >= threshold)[0]
            tsss.append(t[idx[0]] if len(idx) > 0 else t[-1])
        Tss_vals.append(np.mean(tsss))
    ax_Tss.plot(N_vals, Tss_vals, 's-', color='#e9c46a', lw=2)

    # ── 1.4a: perfiles radiales (para el N más grande disponible) ──
    N_big = Ns[-1]
    ax_prof.set_title(f"1.4  Perfiles radiales  (N={N_big}, promedio sobre realizaciones)")
    ax_prof.set_xlabel("S  [m]"); ax_prof.set_ylabel("valor (normalizado)")
    # Promediar sobre realizaciones
    all_rho, all_v, all_J = [], [], []
    S_ref = None
    for r in by_N[N_big]:
        all_rho.append(r['density']); all_v.append(r['velocity']); all_J.append(r['flux'])
        S_ref = r['S']
    rho_mean = np.mean(all_rho, axis=0)
    v_mean   = np.mean(all_v,   axis=0)
    J_mean   = np.mean(all_J,   axis=0)
    # Normalizar para mostrar en misma escala
    def safe_norm(arr):
        m = np.max(arr)
        return arr / m if m > 0 else arr
    ax_prof.plot(S_ref, safe_norm(rho_mean), color='#58a6ff', lw=2, label=r"$\langle\rho_f^{in}\rangle$ (norm)")
    ax_prof.plot(S_ref, safe_norm(v_mean),   color='#3fb950', lw=2, label=r"$|\langle v_f^{in}\rangle|$ (norm)")
    ax_prof.plot(S_ref, safe_norm(J_mean),   color='#f0883e', lw=2.5, label=r"$J_{in}$ (norm)")
    ax_prof.axvline(r_inner, color='#8b949e', ls=':', lw=1, label=f"r_inner={r_inner}")
    ax_prof.legend(fontsize=9, labelcolor='#c9d1d9', facecolor='#161b22', edgecolor='#30363d')

    # ── 1.4b: Jin @ S≈2 vs N ──
    ax_Jin.set_title("1.4  Jin en S≈2 vs N")
    ax_Jin.set_xlabel("N"); ax_Jin.set_ylabel("Jin(S≈2)")
    Jin_at2 = []
    for N in Ns:
        vals = []
        for r in by_N[N]:
            # encontrar el shell más cercano a S=2
            idx = np.argmin(np.abs(r['S'] - 2.0))
            vals.append(r['flux'][idx])
        Jin_at2.append(np.mean(vals))
    ax_Jin.plot(N_vals, Jin_at2, 'D-', color='#da3633', lw=2)

    fig.suptitle("Sistema 1 — Análisis de Scanning Rate, Fu(t) y Perfiles Radiales",
                 color='#e6edf3', fontsize=14, y=0.98)

    plt.savefig(out_prefix + ".png", dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    print(f"Figura guardada: {out_prefix}.png")
    plt.show()


# ──────────────────────────────────────────────────────────────
# 5. Entrypoint
# ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Análisis puntos 1.2, 1.3, 1.4 — Sistema 1")
    parser.add_argument("files",        nargs='+', help="Archivos output.txt a procesar")
    parser.add_argument("--N",          type=int,   default=None,  help="Número de partículas (si es único)")
    parser.add_argument("--Ns",         type=str,   default=None,  help="N por archivo, comma-separated. Ej: 50,100,100")
    parser.add_argument("--r_outer",    type=float, default=40.0)
    parser.add_argument("--r_inner",    type=float, default=1.0)
    parser.add_argument("--radius",     type=float, default=1.0,   help="Radio de la partícula")
    parser.add_argument("--dS",         type=float, default=0.2,   help="Ancho de capa para perfil radial")
    parser.add_argument("--tol",        type=float, default=0.5,   help="Tolerancia para detección de colisiones")
    parser.add_argument("--out",        type=str,   default="sistema1_analisis")
    args = parser.parse_args()

    # Determinar N para cada archivo
    if args.Ns:
        Ns_list = [int(x) for x in args.Ns.split(',')]
        assert len(Ns_list) == len(args.files), "Debe haber un N por archivo"
    elif args.N:
        Ns_list = [args.N] * len(args.files)
    else:
        # Intentar leer N del propio archivo (primer frame)
        Ns_list = []
        for f in args.files:
            times, states = parse_output(f)
            Ns_list.append(states.shape[1] if len(states) > 0 else 0)

    results = []
    for filepath, N in zip(args.files, Ns_list):
        print(f"Procesando {filepath}  (N={N}) ...")
        times, states = parse_output(filepath)
        if len(times) == 0:
            print(f"  ⚠ Archivo vacío o sin datos, saltando.")
            continue

        particle_state, cfc_series, _, _ = reconstruct_states(
            times, states, args.r_outer, args.r_inner, args.radius, tol=args.tol)

        Fu = compute_Fu(particle_state)

        S, density, velocity, flux = compute_radial_profiles(
            times, states, particle_state, args.r_inner, args.r_outer, dS=args.dS)

        results.append({
            'N':         N,
            'times':     times,
            'cfc_series': cfc_series,
            'Fu':        Fu,
            'S':         S,
            'density':   density,
            'velocity':  velocity,
            'flux':      flux,
        })
        J, se = compute_J(times, cfc_series)
        print(f"  → J = {J:.4f} ± {se:.4f}   Fest ≈ {Fu[-10:].mean():.3f}")

    if not results:
        print("No se procesó ningún archivo.")
        return

    plot_all(results, r_outer=args.r_outer, r_inner=args.r_inner,
             particle_radius=args.radius, out_prefix=args.out)


if __name__ == "__main__":
    main()