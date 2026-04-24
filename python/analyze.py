"""
analyze.py  –  Puntos 1.1, 1.2, 1.3 y 1.4 del Sistema 1 (Scanning rate en recinto circular)

Formato esperado del output.txt:
    # t
    # x y vx vy state
    <tiempo>
    x1 y1 vx1 vy1 state1
    ...

Uso (análisis 1.2/1.3/1.4):
    python analyze.py output.txt [output2.txt ...] [--Ns 200,200] [--out plot]

Uso (timing 1.1):
    python analyze.py --timing runs/timing.txt [--out timing]
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

# ──────────────────────────────────────────────────────────────
# 1. Parser del output.txt
# ──────────────────────────────────────────────────────────────

def parse_output(filepath):
    """
    Devuelve:
        times          : array (T,)
        states         : array (T, N, 4)   x, y, vx, vy
        particle_states: array (T, N)      0=fresca, 1=usada  (-1 si no disponible)
    """
    times, frames, state_frames = [], [], []
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
        sf = []
        while i < len(lines):
            parts = lines[i].split()
            if len(parts) in (4, 5):
                frame.append([float(v) for v in parts[:4]])
                sf.append(int(parts[4]) if len(parts) == 5 else -1)
                i += 1
            else:
                break
        if frame:
            times.append(t)
            frames.append(frame)
            state_frames.append(sf)

    times          = np.array(times)
    states         = np.array(frames)       # (T, N, 4)
    particle_states = np.array(state_frames) # (T, N)
    return times, states, particle_states


# ──────────────────────────────────────────────────────────────
# 2. Reconstruir estados y detectar colisiones
# ──────────────────────────────────────────────────────────────

def reconstruct_states(times, states, particle_states_raw,
                       r_outer, r_inner, particle_radius, tol=0.5):
    """
    Usa el estado real de Java si está disponible (columna 5),
    si no, infiere desde posiciones.

    Devuelve:
        particle_state : (T, N)   0=fresca, 1=usada
        cfc_series     : (T,)     conteo acumulado fresca→usada
    """
    T, N, _ = states.shape
    has_state = (particle_states_raw.shape == (T, N) and
                 np.all(particle_states_raw >= 0))

    if has_state:
        particle_state = particle_states_raw.copy()

        # Cfc: contar transiciones 0→1 entre snapshots consecutivos
        cfc = 0
        cfc_series = np.zeros(T)
        prev = particle_state[0].copy()
        for t in range(1, T):
            cfc += int(np.sum((prev == 0) & (particle_state[t] == 1)))
            cfc_series[t] = cfc
            prev = particle_state[t].copy()
        return particle_state, cfc_series

    # Fallback: inferencia por posición
    x   = states[:, :, 0]
    y   = states[:, :, 1]
    vx  = states[:, :, 2]
    vy  = states[:, :, 3]
    dist = np.sqrt(x**2 + y**2)

    contact_inner = np.abs(dist - (r_inner + particle_radius)) < tol
    contact_outer = np.abs(dist - (r_outer - particle_radius)) < tol

    particle_state = np.zeros((T, N), dtype=int)
    current_state  = np.zeros(N, dtype=int)
    cfc = 0
    cfc_series = np.zeros(T)

    for t in range(T):
        for j in range(N):
            if contact_inner[t, j] and current_state[j] == 0:
                current_state[j] = 1
                cfc += 1
            elif contact_outer[t, j] and current_state[j] == 1:
                current_state[j] = 0
        particle_state[t] = current_state.copy()
        cfc_series[t] = cfc

    return particle_state, cfc_series


# ──────────────────────────────────────────────────────────────
# 3. Observables
# ──────────────────────────────────────────────────────────────

def compute_J(times, cfc_series):
    slope, _, _, _, se = stats.linregress(times, cfc_series)
    return slope, se


def compute_Fu(particle_state):
    return particle_state.mean(axis=1)


def compute_radial_profiles(times, states, particle_state, r_inner, r_outer, dS=0.2):
    """
    Perfiles radiales de partículas frescas apuntando al centro.

    Densidad: promedio temporal sobre TODOS los snapshots (incluye ceros).
    Velocidad: promedio ponderado por número de partículas.

    Devuelve: S_centers, avg_density, avg_velocity, flux
    """
    num_shells = int(np.ceil((r_outer - r_inner) / dS))
    shell_edges = r_inner + np.arange(num_shells + 1) * dS

    # Acumuladores
    sum_density    = np.zeros(num_shells)       # sum(n_in/area) sobre todos los snapshots
    sum_vel_num    = np.zeros(num_shells)       # sum de |v_radial| individual
    count_particles = np.zeros(num_shells, dtype=int)  # total partículas acumuladas

    x   = states[:, :, 0]
    y   = states[:, :, 1]
    vx  = states[:, :, 2]
    vy  = states[:, :, 3]
    dist = np.sqrt(x**2 + y**2)
    T = len(times)

    areas = np.pi * (shell_edges[1:]**2 - shell_edges[:-1]**2)

    for t in range(T):
        fresh_mask   = (particle_state[t] == 0)
        radial_vel   = (x[t] * vx[t] + y[t] * vy[t]) / (dist[t] + 1e-15)
        toward_center = (radial_vel < 0) & fresh_mask

        for k in range(num_shells):
            in_shell = (toward_center &
                        (dist[t] >= shell_edges[k]) &
                        (dist[t] <  shell_edges[k + 1]))
            n_in = int(np.sum(in_shell))
            sum_density[k] += n_in / areas[k]   # 0 si n_in=0, promediado luego por T
            if n_in > 0:
                sum_vel_num[k]    += float(np.sum(np.abs(radial_vel[in_shell])))
                count_particles[k] += n_in

    avg_density  = sum_density / T
    valid        = count_particles > 0
    avg_velocity = np.where(valid,
                            sum_vel_num / np.where(count_particles > 0, count_particles, 1),
                            0.0)
    flux = avg_density * avg_velocity

    S_centers = r_inner + (np.arange(num_shells) + 0.5) * dS
    return S_centers, avg_density, avg_velocity, flux


# ──────────────────────────────────────────────────────────────
# 4. Plot de timing (1.1)
# ──────────────────────────────────────────────────────────────

def plot_timing(timing_file, out_prefix="analisis_timing"):
    data = np.loadtxt(timing_file, comments='#')
    if data.ndim == 1:
        data = data.reshape(1, -1)
    Ns      = data[:, 0].astype(int)
    elapsed = data[:, 1] / 1000.0  # ms → s

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(Ns, elapsed, 'o-', color='#457b9d', lw=2, ms=7)
    ax.set_xlabel("N  (número de partículas)", fontsize=12)
    ax.set_ylabel("Tiempo de ejecución  [s]", fontsize=12)
    ax.set_title("1.1  Tiempo de ejecución vs N  (tf = 5 s)", fontsize=13)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    plt.savefig(out_prefix + ".png", dpi=150, bbox_inches='tight')
    print(f"Figura guardada: {out_prefix}.png")
    plt.show()


# ──────────────────────────────────────────────────────────────
# 5. Plot principal (1.2, 1.3, 1.4)
# ──────────────────────────────────────────────────────────────

COLORS = ['#e63946', '#457b9d', '#2a9d8f', '#e9c46a', '#f4a261']


def plot_all(results, r_outer=40, r_inner=1, particle_radius=1, out_prefix="plot"):
    Ns   = sorted(set(r['N'] for r in results))
    by_N = {N: [r for r in results if r['N'] == N] for N in Ns}

    fig = plt.figure(figsize=(18, 14))
    fig.patch.set_facecolor('white')
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.45, wspace=0.38)

    ax_cfc  = fig.add_subplot(gs[0, 0])
    ax_J    = fig.add_subplot(gs[0, 1])
    ax_Fu   = fig.add_subplot(gs[1, 0])
    ax_Fest = fig.add_subplot(gs[1, 1])
    ax_Tss  = fig.add_subplot(gs[1, 2])
    ax_prof = fig.add_subplot(gs[2, :2])
    ax_Jin  = fig.add_subplot(gs[2, 2])

    for ax in [ax_cfc, ax_J, ax_Fu, ax_Fest, ax_Tss, ax_prof, ax_Jin]:
        ax.set_facecolor('white')
        ax.tick_params(colors='black')
        ax.xaxis.label.set_color('black')
        ax.yaxis.label.set_color('black')
        ax.title.set_color('black')
        for spine in ax.spines.values():
            spine.set_edgecolor('black')

    # ── 1.2a: Cfc(t) ──
    ax_cfc.set_title("1.2  Cfc(t) — conteo acumulado fresca→usada")
    ax_cfc.set_xlabel("t  [s]"); ax_cfc.set_ylabel("Cfc(t)")
    for i, N in enumerate(Ns):
        r0 = by_N[N][0]
        c  = COLORS[i % len(COLORS)]
        ax_cfc.plot(r0['times'], r0['cfc_series'], color=c, lw=1.5, label=f"N={N}")
        J, _ = compute_J(r0['times'], r0['cfc_series'])
        t_line = np.linspace(r0['times'][0], r0['times'][-1], 100)
        ax_cfc.plot(t_line, J * t_line, color=c, ls='--', lw=1, alpha=0.6)
    ax_cfc.legend(fontsize=8)

    # ── 1.2b: <J>(N) con barra de error ──
    ax_J.set_title("1.2  Scanning rate  ⟨J⟩  vs  N")
    ax_J.set_xlabel("N"); ax_J.set_ylabel("J  [colisiones/s]")
    J_means, J_stds, N_vals = [], [], []
    for N in Ns:
        Js = [compute_J(r['times'], r['cfc_series'])[0] for r in by_N[N]]
        J_means.append(np.mean(Js)); J_stds.append(np.std(Js)); N_vals.append(N)
    ax_J.errorbar(N_vals, J_means, yerr=J_stds, fmt='o-',
                  color='#457b9d', ecolor='#e63946', capsize=5, lw=2)

    # ── 1.3a: Fu(t) ──
    ax_Fu.set_title("1.3  Fu(t) — fracción de partículas usadas")
    ax_Fu.set_xlabel("t  [s]"); ax_Fu.set_ylabel("Fu(t) = Nu/N")
    for i, N in enumerate(Ns):
        for j, r in enumerate(by_N[N]):
            ax_Fu.plot(r['times'], r['Fu'], color=COLORS[i % len(COLORS)],
                       lw=1.2, alpha=0.5 if j > 0 else 1.0,
                       label=f"N={N}" if j == 0 else None)
    ax_Fu.legend(fontsize=8)
    ax_Fu.set_ylim(-0.05, 1.05)

    # ── 1.3b: Fest vs N ──
    ax_Fest.set_title("1.3  Fest (valor estacionario) vs N")
    ax_Fest.set_xlabel("N"); ax_Fest.set_ylabel("Fest")
    ax_Fest.set_ylim(0, 0.5)
    Fest_vals = []
    for N in Ns:
        fests = [np.mean(r['Fu'][max(1, int(len(r['Fu']) * 0.7)):]) for r in by_N[N]]
        Fest_vals.append(np.mean(fests))
    ax_Fest.plot(N_vals, Fest_vals, 'o-', color='#2a9d8f', lw=2)

    # ── 1.3c: T_ss vs N ──
    ax_Tss.set_title("1.3  T_estacionario vs N")
    ax_Tss.set_xlabel("N"); ax_Tss.set_ylabel("T_ss  [s]")
    ax_Tss.set_ylim(bottom=0)
    Tss_vals = []
    for N in Ns:
        tsss = []
        for r in by_N[N]:
            fu   = r['Fu']
            t    = r['times']
            fest = np.mean(fu[int(len(fu) * 0.7):])
            idx  = np.where(fu >= 0.9 * fest)[0]
            tsss.append(t[idx[0]] if len(idx) > 0 else t[-1])
        Tss_vals.append(np.mean(tsss))
    ax_Tss.plot(N_vals, Tss_vals, 's-', color='#e9c46a', lw=2)

    # ── 1.4a: perfiles radiales ──
    N_big = Ns[-1]
    ax_prof.set_title(f"1.4  Perfiles radiales (N={N_big}, promedio sobre realizaciones)")
    ax_prof.set_xlabel("S  [m]"); ax_prof.set_ylabel("valor (normalizado)")
    all_rho, all_v, all_J = [], [], []
    S_ref = None
    for r in by_N[N_big]:
        all_rho.append(r['density']); all_v.append(r['velocity']); all_J.append(r['flux'])
        S_ref = r['S']
    rho_mean = np.mean(all_rho, axis=0)
    v_mean   = np.mean(all_v,   axis=0)
    J_mean   = np.mean(all_J,   axis=0)

    def safe_norm(arr):
        m = np.max(arr)
        return arr / m if m > 0 else arr

    ax_prof.plot(S_ref, safe_norm(rho_mean), color='#457b9d', lw=2,
                 label=r"$\langle\rho_f^{in}\rangle$ (norm)")
    ax_prof.plot(S_ref, safe_norm(v_mean),   color='#2a9d8f', lw=2,
                 label=r"$|\langle v_f^{in}\rangle|$ (norm)")
    ax_prof.plot(S_ref, safe_norm(J_mean),   color='#e63946', lw=2.5,
                 label=r"$J_{in}$ (norm)")
    ax_prof.axvline(r_inner + particle_radius, color='gray', ls=':', lw=1,
                    label=f"S_min = {r_inner + particle_radius}")
    ax_prof.legend(fontsize=9)

    # ── 1.4b: Jin @ S≈2 vs N ──
    ax_Jin.set_title("1.4  Jin, ρ, v  en S≈2  vs  N")
    ax_Jin.set_xlabel("N")
    Jin_at2, rho_at2, v_at2 = [], [], []
    for N in Ns:
        j_vals, r_vals, v_vals = [], [], []
        for r in by_N[N]:
            idx = np.argmin(np.abs(r['S'] - 2.0))
            j_vals.append(r['flux'][idx])
            r_vals.append(r['density'][idx])
            v_vals.append(r['velocity'][idx])
        Jin_at2.append(np.mean(j_vals))
        rho_at2.append(np.mean(r_vals))
        v_at2.append(np.mean(v_vals))

    ax_Jin.plot(N_vals, Jin_at2, 'D-',  color='#e63946', lw=2, label=r"$J_{in}$")
    ax_Jin.plot(N_vals, rho_at2, 'o--', color='#457b9d', lw=1.5, label=r"$\langle\rho\rangle$")
    ax_Jin.plot(N_vals, v_at2,   's--', color='#2a9d8f', lw=1.5, label=r"$|\langle v\rangle|$")
    ax_Jin.legend(fontsize=8)

    fig.suptitle("Sistema 1 — Scanning Rate, Fu(t) y Perfiles Radiales",
                 fontsize=14, y=0.98)

    plt.savefig(out_prefix + ".png", dpi=150, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    print(f"Figura guardada: {out_prefix}.png")
    plt.show()


# ──────────────────────────────────────────────────────────────
# 6. Entrypoint
# ──────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Análisis puntos 1.1-1.4 — Sistema 1")
    parser.add_argument("files",     nargs='*', help="Archivos output.txt")
    parser.add_argument("--timing",  type=str,  default=None,
                        help="Modo 1.1: ruta al archivo timing.txt")
    parser.add_argument("--N",       type=int,   default=None)
    parser.add_argument("--Ns",      type=str,   default=None,
                        help="N por archivo, comma-separated. Ej: 50,100,100")
    parser.add_argument("--r_outer", type=float, default=40.0)
    parser.add_argument("--r_inner", type=float, default=1.0)
    parser.add_argument("--radius",  type=float, default=1.0)
    parser.add_argument("--dS",      type=float, default=0.2)
    parser.add_argument("--tol",     type=float, default=0.5)
    parser.add_argument("--out",     type=str,   default="sistema1_analisis")
    args = parser.parse_args()

    if args.timing:
        plot_timing(args.timing, out_prefix=args.out)
        return

    if not args.files:
        parser.print_help()
        return

    if args.Ns:
        Ns_list = [int(x) for x in args.Ns.split(',')]
        assert len(Ns_list) == len(args.files)
    elif args.N:
        Ns_list = [args.N] * len(args.files)
    else:
        Ns_list = []
        for f in args.files:
            times, states, _ = parse_output(f)
            Ns_list.append(states.shape[1] if len(states) > 0 else 0)

    results = []
    for filepath, N in zip(args.files, Ns_list):
        print(f"Procesando {filepath}  (N={N}) ...")
        times, states, particle_states_raw = parse_output(filepath)
        if len(times) == 0:
            print("  Archivo vacío, saltando.")
            continue

        particle_state, cfc_series = reconstruct_states(
            times, states, particle_states_raw,
            args.r_outer, args.r_inner, args.radius, tol=args.tol)

        Fu = compute_Fu(particle_state)

        S, density, velocity, flux = compute_radial_profiles(
            times, states, particle_state, args.r_inner, args.r_outer, dS=args.dS)

        results.append({
            'N':          N,
            'times':      times,
            'cfc_series': cfc_series,
            'Fu':         Fu,
            'S':          S,
            'density':    density,
            'velocity':   velocity,
            'flux':       flux,
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
