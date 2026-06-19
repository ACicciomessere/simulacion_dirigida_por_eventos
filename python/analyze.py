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
        states : array (T, N, 5)       x, y, vx, vy  para cada snapshot
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
            if len(parts) in (4, 5):
                frame.append([float(v) for v in parts])
                i += 1
            else:
                break
        if frame:
            times.append(t)
            frames.append(frame)

    times  = np.array(times)
    states = np.array(frames)   # (T, N, 5)
    return times, states


# ──────────────────────────────────────────────────────────────
# 2. Reconstruir estados fresca/usada y detectar colisiones
# ──────────────────────────────────────────────────────────────

# def reconstruct_states(times, states, r_outer, r_inner, particle_radius, tol=0.5):
#     """
#     Para cada snapshot, infiere si cada partícula es 'fresca' (0) o 'usada' (1).

#     Reglas (igual que la simulación):
#       - Arranca todo el mundo fresca.
#       - Si en cualquier snapshot una partícula está muy cerca del obstáculo interno
#         (dist_to_center ≈ r_inner + particle_radius), marcamos ese evento como
#         colisión con el centro → la partícula pasa a 'usada'.
#       - Si una partícula usada llega muy cerca de la pared externa
#         (dist_to_center ≈ r_outer - particle_radius), vuelve a ser 'fresca'.

#     Devuelve:
#         particle_state  : array (T, N)   0 = fresca, 1 = usada  en cada instante
#         cfc_series      : array (T,)     conteo acumulado de transiciones fresh→used
#         hit_center_mask : array (T, N)   True si en ese snapshot hubo colisión con centro
#         hit_outer_mask  : array (T, N)   True si en ese snapshot hubo colisión con pared
#     """
#     T, N, _ = states.shape
#     x  = states[:, :, 0]   # (T, N)
#     y  = states[:, :, 1]
#     vx = states[:, :, 2]
#     vy = states[:, :, 3]

#     dist = np.sqrt(x**2 + y**2)   # distancia al centro en cada snapshot

#     # Detectar snapshots donde la partícula está "en contacto" con cada pared
#     contact_inner = np.abs(dist - (r_inner + particle_radius)) < tol
#     contact_outer = np.abs(dist - (r_outer - particle_radius)) < tol

#     # Reconstruir estado snapshot a snapshot
#     particle_state = np.zeros((T, N), dtype=int)  # 0=fresca, 1=usada
#     current_state  = np.zeros(N, dtype=int)        # estado actual de cada partícula
#     cfc = 0
#     cfc_series = np.zeros(T)

#     for t in range(T):
#         for j in range(N):
#             if contact_inner[t, j] and current_state[j] == 0:
#                 current_state[j] = 1   # fresca → usada
#                 cfc += 1
#             elif contact_outer[t, j] and current_state[j] == 1:
#                 current_state[j] = 0   # usada → fresca
#         particle_state[t] = current_state.copy()
#         cfc_series[t] = cfc

#     return particle_state, cfc_series, contact_inner, contact_outer


def reconstruct_states(times, states, r_outer=40.0, r_inner=1.0, particle_radius=1.0, tol=0.5):
    """
    Si el output trae la columna 'fresh' (5ª columna), la usamos directamente
    (1=fresca, 0=usada, según Java). Si no, reconstruimos por proximidad a las
    paredes — esto soporta corridas viejas sin la columna 'fresh'.
    Devuelve (particle_state, cfc_series, None, None) para mantener la firma.
    """
    T, N, ncols = states.shape

    if ncols >= 5:
        # Columna 4 es el flag 'fresh' que escribe la simulación (1=fresca → 0=usada).
        particle_state = (1 - states[:, :, 4]).astype(int)
    else:
        # Output legado de 4 columnas: inferir estado por proximidad a las paredes.
        x  = states[:, :, 0]
        y  = states[:, :, 1]
        dist = np.sqrt(x**2 + y**2)
        contact_inner = np.abs(dist - (r_inner + particle_radius)) < tol
        contact_outer = np.abs(dist - (r_outer - particle_radius)) < tol

        particle_state = np.zeros((T, N), dtype=int)
        current = np.zeros(N, dtype=int)
        for t in range(T):
            for j in range(N):
                if contact_inner[t, j] and current[j] == 0:
                    current[j] = 1
                elif contact_outer[t, j] and current[j] == 1:
                    current[j] = 0
            particle_state[t] = current

    # Cfc: transiciones fresca→usada acumuladas en el tiempo.
    cfc_series = np.zeros(T)
    cfc = 0
    for t in range(1, T):
        newly_used = (particle_state[t] == 1) & (particle_state[t-1] == 0)
        cfc += int(np.sum(newly_used))
        cfc_series[t] = cfc

    return particle_state, cfc_series, None, None

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

def fit_relaxation(times, Fu, tail_frac=0.2):
    """
    Ajusta la curva de relajación  Fu(t) = Fest · (1 − exp(−t / τ)).

    Fu(t) crece desde 0 y satura en Fest. Para los N grandes la corrida puede
    terminar ANTES de llegar al estacionario (T_ss > t_final), por lo que el
    promedio de la cola subestima Fest; el ajuste extrapola al valor asintótico
    real. Todo el transitorio entra al fit, así que es robusto a las
    fluctuaciones del estacionario.

    Devuelve (Fest, tau). Si el ajuste no converge, cae al promedio de la cola
    para Fest y a un τ estimado por cruce.
    """
    times = np.asarray(times, dtype=float)
    Fu    = np.asarray(Fu,    dtype=float)
    n = len(Fu)
    Fest0 = np.mean(Fu[-max(1, int(n * tail_frac)):]) if n else 0.0

    if n < 3:
        return float(Fest0), float(times[-1]) if n else 0.0

    def model(t, Fest, tau):
        return Fest * (1.0 - np.exp(-t / tau))

    try:
        from scipy.optimize import curve_fit
        p0 = [Fest0 if Fest0 > 0 else 0.1, max(times[-1] * 0.2, 1e-3)]
        popt, _ = curve_fit(model, times, Fu, p0=p0, maxfev=10000)
        return float(popt[0]), float(abs(popt[1]))
    except Exception:
        win = max(1, int(n * 0.05))
        smooth = np.convolve(Fu, np.ones(win) / win, mode='same') if win > 1 else Fu
        idx = np.where(smooth >= 0.632 * Fest0)[0]   # 1 - 1/e
        tau = float(times[idx[0]]) if len(idx) else float(times[-1])
        return float(Fest0), tau


def steady_state_Fest(times, Fu):
    """
    Fest (fracción usada en equilibrio) = promedio temporal de Fu(t) en el
    régimen estacionario.

    Fu(t) sube desde 0 y satura en una meseta. El valor de equilibrio es,
    por definición, el promedio de Fu sobre esa meseta. Para aislarla:
      1. Se estima el tiempo de relajación τ con el ajuste de relajación.
      2. Se descarta el transitorio (~3·τ, donde Fu ya alcanzó el 95% de la
         meseta) y se promedia Fu sobre el resto de la corrida.
      3. Se conserva siempre al menos la última mitad del run (min con T/2),
         de modo que aunque τ sea grande el promedio use suficientes muestras.

    Este estimador es directo y de baja varianza: a diferencia de extrapolar
    el ajuste exponencial a t→∞ (donde Fest y τ están fuertemente
    correlacionados y el resultado es muy ruidoso), aquí se promedia el dato
    crudo ya estacionario, que es lo que pide la consigna.
    """
    times = np.asarray(times, dtype=float)
    Fu    = np.asarray(Fu,    dtype=float)
    if len(Fu) == 0:
        return 0.0
    _, tau = fit_relaxation(times, Fu)
    T = times[-1]
    t0 = min(3.0 * tau, 0.5 * T)          # inicio del régimen estacionario
    mask = times >= t0
    return float(np.mean(Fu[mask])) if mask.any() else float(Fu[-1])


def compute_Tss(times, Fu, frac=0.9):
    """
    Tiempo al estacionario: instante en que la curva ajustada alcanza `frac`
    (90%) de Fest, es decir  T_ss = τ · ln(1 / (1 − frac)).
    """
    _, tau = fit_relaxation(times, Fu)
    return float(tau * np.log(1.0 / (1.0 - frac)))

def compute_radial_profiles(times, states, particle_state, r_inner, r_outer,
                             dS=0.2):
    """
    Para cada snapshot: selecciona partículas frescas con velocidad radial
    apuntando al centro (Rj·vj < 0), las agrupa por capa S y acumula conteos.
    Devuelve, para cada capa:
        <rho_fin>(S) = (total_count_en_capa) / (n_frames * area_capa)
        <|v_fin|>(S) = (suma de |v_radial|) / (total_count_en_capa)
        Jin(S)       = <rho_fin>(S) * <|v_fin|>(S)
    Misma definición que TP4: el promedio temporal incluye los snapshots
    vacíos (no se excluyen), de modo que la densidad refleja la fracción
    real del tiempo que la capa está ocupada.
    """
    num_shells = int(np.ceil((r_outer - r_inner) / dS))
    shell_edges = r_inner + np.arange(num_shells + 1) * dS
    areas       = np.pi * (shell_edges[1:]**2 - shell_edges[:-1]**2)

    counts   = np.zeros(num_shells)         # total partículas-fresh-in en cada capa, sumado en t
    vrad_sum = np.zeros(num_shells)         # suma de |v_radial| sobre todas esas partículas

    x  = states[:, :, 0]
    y  = states[:, :, 1]
    vx = states[:, :, 2]
    vy = states[:, :, 3]
    dist = np.sqrt(x**2 + y**2)

    n_frames = len(times)
    for t in range(n_frames):
        fresh_mask    = (particle_state[t] == 0)
        radial_vel    = (x[t] * vx[t] + y[t] * vy[t]) / (dist[t] + 1e-15)
        toward_center = (radial_vel < 0) & fresh_mask
        d_in          = dist[t][toward_center]
        v_in          = np.abs(radial_vel[toward_center])

        for k in range(num_shells):
            in_shell = (d_in >= shell_edges[k]) & (d_in < shell_edges[k+1])
            counts  [k] += in_shell.sum()
            vrad_sum[k] += v_in[in_shell].sum()

    avg_density  = counts / max(n_frames, 1) / areas
    with np.errstate(invalid="ignore", divide="ignore"):
        avg_velocity = np.where(counts > 0, vrad_sum / np.maximum(counts, 1), 0.0)
    flux = avg_density * avg_velocity

    S_centers = r_inner + (np.arange(num_shells) + 0.5) * dS
    return S_centers, avg_density, avg_velocity, flux


# ──────────────────────────────────────────────────────────────
# 4. Plots
# ──────────────────────────────────────────────────────────────

COLORS = ['#e63946', '#457b9d', '#2a9d8f', '#e9c46a', '#f4a261']


def plot_1_4_separate(results, r_inner=1.0, particle_radius=1.0, out_prefix="plot"):
    """
    Genera UN PNG por panel del item 1.4 (perfiles radiales): rho, |v|, Jin.
    Promedia las realizaciones por N y usa el mismo estilo/colores que plot_all
    (leyenda en el gráfico, paleta COLORS, línea vertical en S_min).

    Archivos generados:
        {out_prefix}_1_4_rho.png
        {out_prefix}_1_4_v.png
        {out_prefix}_1_4_jin.png
    """
    Ns = sorted(set(r['N'] for r in results))
    if not Ns:
        return
    by_N = {N: [r for r in results if r['N'] == N] for N in Ns}

    panels = [
        ("rho", r"Densidad $\langle\rho_{in}\rangle$(S)",  r"$\rho$ (m$^{-2}$)",            'density'),
        ("v",   r"Velocidad radial $\langle v_{in}\rangle$(S)", "|v| (m/s)",                'velocity'),
        ("jin", r"Flujo $J_{in}$(S)",                      r"$J_{in}$ (m$^{-2}$ s$^{-1}$)", 'flux'),
    ]

    cmap   = plt.cm.plasma
    colors = [cmap(i / max(len(Ns) - 1, 1)) for i in range(len(Ns))]

    for tag, title, ylabel, key in panels:
        fig, ax = plt.subplots(figsize=(8, 5))

        for color, N in zip(colors, Ns):
            S_ref  = by_N[N][0]['S']
            stack  = np.array([r[key] for r in by_N[N]])
            y_mean = stack.mean(axis=0)
            ax.plot(S_ref, y_mean, color=color, lw=1.2)

        ax.set_title(title)
        ax.set_xlabel("S (m)")
        ax.set_ylabel(ylabel)
        ax.grid(True, alpha=0.3)

        # Colorbar gradual en lugar de leyenda (requerimiento de la consigna).
        sm = plt.cm.ScalarMappable(cmap=cmap,
                                   norm=plt.Normalize(vmin=min(Ns), vmax=max(Ns)))
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label="N")

        out = f"{out_prefix}_1_4_{tag}.png"
        fig.tight_layout()
        fig.savefig(out, dpi=150, bbox_inches='tight')
        plt.close(fig)
        print(f"Figura guardada: {out}")


def plot_all(results, r_outer=40, r_inner=1, particle_radius=1, out_prefix="plot"):
    """
    results: list of dicts, one per (N, realization):
        { 'N': int, 'times': ..., 'cfc_series': ..., 'Fu': ...,
          'S': ..., 'density': ..., 'velocity': ..., 'flux': ... }
    """
    Ns = sorted(set(r['N'] for r in results))

    # ── Agrupar por N ──
    by_N = {N: [r for r in results if r['N'] == N] for N in Ns}

    fig = plt.figure(figsize=(20, 20))
    fig.patch.set_facecolor('white')
    gs = gridspec.GridSpec(3, 3, figure=fig, hspace=0.85, wspace=0.7)

    ax_cfc  = fig.add_subplot(gs[0, 0])   # 1.2a: Cfc(t)
    ax_J    = fig.add_subplot(gs[0, 1])   # 1.2b: <J>(N)
    ax_Fu   = fig.add_subplot(gs[1, 0])   # 1.3a: Fu(t)
    ax_Fest = fig.add_subplot(gs[1, 1])   # 1.3b: Fest vs N
    ax_Tss  = fig.add_subplot(gs[1, 2])   # 1.3c: T_ss vs N
    prof_grid = gs[2, 0].subgridspec(3, 1, hspace=0.6)
    ax_rho   = fig.add_subplot(prof_grid[0, 0])  # 1.4a: rho(S)
    ax_vel   = fig.add_subplot(prof_grid[1, 0])  # 1.4a: |v|(S)
    ax_prof  = fig.add_subplot(prof_grid[2, 0])  # 1.4a: Jin(S)
    ax_Jin   = fig.add_subplot(gs[2, 1:])       # 1.4b: Jin @ S≈2 vs N (zoom)

    for ax in [ax_cfc, ax_J, ax_Fu, ax_Fest, ax_Tss, ax_rho, ax_vel, ax_prof, ax_Jin]:
        ax.set_facecolor('white')
        ax.tick_params(colors='black')
        ax.xaxis.label.set_color('black')
        ax.yaxis.label.set_color('black')
        ax.title.set_color('black')
        for spine in ax.spines.values():
            spine.set_edgecolor('black')

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
    ax_cfc.legend(fontsize=8, labelcolor='black', facecolor='white', edgecolor='black',
                  bbox_to_anchor=(1.05, 1), loc='upper left')

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
    ax_Fu.set_title("1.3a  Fu(t) — fracción de partículas usadas", fontweight='bold')
    ax_Fu.set_xlabel("t  [s]", fontsize=10); ax_Fu.set_ylabel("Fu(t) = Nu/N", fontsize=10)
    
    # Una sola curva por N: el promedio de las realizaciones (valores concretos,
    # sin la banda de dispersión que generaba superponer todas las corridas).
    cmap_fu = plt.cm.turbo
    colors_fu = [cmap_fu(i / max(len(Ns) - 1, 1)) for i in range(len(Ns))]
    for color, N in zip(colors_fu, Ns):
        runs_N = by_N[N]
        # Grilla temporal común (las corridas pueden diferir en longitud);
        # se interpola cada Fu y se promedia punto a punto.
        t_end = min(r['times'][-1] for r in runs_N)
        t_grid = np.linspace(0.0, t_end, 400)
        fu_stack = np.array([np.interp(t_grid, r['times'], r['Fu']) for r in runs_N])
        ax_Fu.plot(t_grid, fu_stack.mean(axis=0), color=color, lw=2, label=f"N={N}")

    ax_Fu.legend(fontsize=9, labelcolor='black', facecolor='white', edgecolor='#444444',
                bbox_to_anchor=(1.05, 1), loc='upper left', framealpha=0.95)
    ax_Fu.set_ylim(-0.05, 0.25)

    ax_Fu.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)

    # ── 1.3b: Fest vs N ──
    ax_Fest.set_title("1.3b  Fest (equilibrio) vs N", fontweight='bold')
    ax_Fest.set_xlabel("N", fontsize=10); ax_Fest.set_ylabel(r"$F_{est}$", fontsize=10)
    Fest_vals, Fest_errs = [], []
    for N in Ns:
        # Fest = promedio de Fu(t) en el régimen estacionario (tras descartar
        # el transitorio). Estimador directo y de baja varianza, a diferencia
        # de extrapolar el ajuste exponencial que aplanaba N=500/600/700.
        fests = [steady_state_Fest(r['times'], r['Fu']) for r in by_N[N]]
        Fest_vals.append(np.mean(fests))
        Fest_errs.append(np.std(fests))
    ax_Fest.errorbar(N_vals, Fest_vals, yerr=Fest_errs, fmt='o-',
                     color='#1f77b4', ecolor='#1f77b4', lw=2, markersize=7,
                     capsize=4, capthick=1.2)
    ax_Fest.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)

    # ── 1.3c: T_ss vs N (tiempo al estacionario) ──
    ax_Tss.set_title("1.3c  T_estacionario vs N", fontweight='bold')
    ax_Tss.set_xlabel("N", fontsize=10); ax_Tss.set_ylabel("T_ss  [s]", fontsize=10)
    Tss_vals = []
    for N in Ns:
        tsss = [compute_Tss(r['times'], r['Fu']) for r in by_N[N]]
        Tss_vals.append(np.mean(tsss))
    ax_Tss.plot(N_vals, Tss_vals, 's-', color='#ffb300', lw=2.5, markersize=8,
               markeredgecolor='#ff8f00', markeredgewidth=1.5)
    ax_Tss.grid(True, alpha=0.25, linestyle=':', linewidth=0.8)

    
    
    # ── 1.4a: perfiles radiales crudos (todos los N) ──
    ax_rho.set_title("1.4a  Perfil de densidad  <rho_in>(S)")
    ax_rho.set_ylabel(r"$\langle \rho_{in} \rangle$")

    ax_vel.set_title("1.4a  Perfil de velocidad  <|v_in|>(S)")
    ax_vel.set_ylabel(r"$\langle |v_{in}| \rangle$")
    
    ax_prof.set_title("1.4a  Perfil radial de flujo  J_in(S)")
    ax_prof.set_xlabel("S  [m]"); ax_prof.set_ylabel(r"$J_{in}(S)$")

    for i, N in enumerate(Ns):
        c = COLORS[i % len(COLORS)]
        all_rho = [r['density'] for r in by_N[N]]
        all_vel = [r['velocity'] for r in by_N[N]]
        all_J = [r['flux'] for r in by_N[N]]
        S_ref = by_N[N][0]['S']
        rho_mean = np.mean(all_rho, axis=0)
        vel_mean = np.mean(all_vel, axis=0)
        J_mean = np.mean(all_J, axis=0)
        ax_rho.plot(S_ref, rho_mean, color=c, lw=2, label=f"N={N}")
        ax_vel.plot(S_ref, vel_mean, color=c, lw=2, label=f"N={N}")
        ax_prof.plot(S_ref, J_mean, color=c, lw=2, label=f"N={N}")

    for ax in [ax_rho, ax_vel, ax_prof]:
        ax.axvline(r_inner + particle_radius, color='gray', ls=':', lw=1,
                    label=f"S_min = {r_inner + particle_radius}")
        ax.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc='upper left')
        ax.grid(True, alpha=0.2, linestyle=':', linewidth=0.8)

    # ── 1.4b: Jin @ S≈2 vs N ──
    # ax_Jin.set_title("1.4  Jin, ρ, v  en S≈2  vs  N")
    # ax_Jin.set_xlabel("N")
    # Jin_at2, rho_at2, v_at2 = [], [], []
    # for N in Ns:
    #     j_vals, r_vals, v_vals = [], [], []
    #     for r in by_N[N]:
    #         idx = np.argmin(np.abs(r['S'] - 2.0))
    #         j_vals.append(r['flux'][idx])
    #         r_vals.append(r['density'][idx])
    #         v_vals.append(r['velocity'][idx])
    #     Jin_at2.append(np.mean(j_vals))
    #     rho_at2.append(np.mean(r_vals))
    #     v_at2.append(np.mean(v_vals))

    # ax_Jin.plot(N_vals, Jin_at2, 'D-',  color='#e63946', lw=2, label=r"$J_{in}$")
    # ax_Jin.plot(N_vals, rho_at2, 'o--', color='#457b9d', lw=1.5, label=r"$\langle\rho\rangle$")
    # ax_Jin.plot(N_vals, v_at2,   's--', color='#2a9d8f', lw=1.5, label=r"$|\langle v\rangle|$")
    # ax_Jin.legend(fontsize=8, bbox_to_anchor=(1.05, 1), loc='upper left')

    fig.suptitle("Sistema 1 — Scanning Rate, Fu(t) y Perfiles Radiales",
                 fontsize=14, y=0.98)

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

    # Item 1.4: un PNG por panel (rho, |v|, Jin) para presentación/comparación con TP4.
    plot_1_4_separate(results, r_inner=args.r_inner,
                       particle_radius=args.radius, out_prefix=args.out)


if __name__ == "__main__":
    main()