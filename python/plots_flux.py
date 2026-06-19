"""
plots_flux.py — Gráficos finales del Sistema 1 (flujo de partículas frescas-in)

Genera dos figuras pedidas en la consigna:

  (A) Perfiles radiales en función de S, promediados sobre realizaciones y
      (opcionalmente) sobre los N elegidos:
          <rho_in>(S),  |<v_in>|(S)  y  Jin(S) = <rho_in>(S)·|<v_in>|(S)
      Las 3 curvas van en la misma figura usando 3 ejes Y (cada observable
      tiene su propia escala).  ->  {out}_perfil_S.png

  (B) Para la capa cercana a S=2, los observables en función de N:
        - Densidad y velocidad en un mismo gráfico con doble eje Y.
              ->  {out}_vsN_rho_v.png
        - Flujo Jin en un gráfico aparte.
              ->  {out}_vsN_jin.png

Reutiliza el parser y el cálculo de perfiles de analyze.py (sin duplicar lógica).

Uso:
    python plots_flux.py runs/N200/run_1/output.txt ... --Ns 200,200,...
    python plots_flux.py runs/N*/run_*/output.txt --Ns <csv> --out outputs/flux

Argumentos clave:
    --S0          capa objetivo para la figura (B)            (default 2.0)
    --profile_N   N a usar en la figura (A); "all" promedia todos (default all)
    --r_inner / --r_outer / --radius / --dS  igual que analyze.py
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt

from analyze import parse_output, reconstruct_states, compute_radial_profiles


# ──────────────────────────────────────────────────────────────
# Construcción de resultados (un dict por archivo/realización)
# ──────────────────────────────────────────────────────────────

def build_results(files, Ns_list, r_inner, r_outer, radius, dS, tol):
    results = []
    for filepath, N in zip(files, Ns_list):
        times, states = parse_output(filepath)
        if len(times) == 0:
            print(f"  ⚠ {filepath} vacío, se saltea.")
            continue
        particle_state, _, _, _ = reconstruct_states(
            times, states, r_outer, r_inner, radius, tol=tol)
        S, density, velocity, flux = compute_radial_profiles(
            times, states, particle_state, r_inner, r_outer, dS=dS)
        results.append({'N': N, 'S': S, 'density': density,
                        'velocity': velocity, 'flux': flux, 'dS': dS})
        print(f"  ✓ {filepath} (N={N})")
    return results


def band_observables(r, S_lo, S_hi):
    """
    Agrega los perfiles por capa de una corrida sobre la banda S_lo ≤ S < S_hi.

    El área de cada capa (centro S, ancho dS) es exactamente 2π·S·dS, así que:
      <rho>_banda = Σ(rho_k·a_k)/Σa_k            (promedio pesado por área)
      |<v>|_banda = Σ(v_k·rho_k·a_k)/Σ(rho_k·a_k) (pesado por densidad·área)
      Jin_banda   = <rho>_banda · |<v>|_banda
    Estos pesos reconstruyen el agregado exacto (Σcounts/Σarea, etc.) a partir
    de los arrays por capa, sin necesidad de los conteos crudos.
    """
    S   = r['S']
    rho = r['density']
    vel = r['velocity']
    a   = 2 * np.pi * S * r['dS']            # área de cada capa
    mask = (S >= S_lo) & (S < S_hi)
    aw = a[mask]
    rho_b = rho[mask]
    vel_b = vel[mask]

    area_tot = aw.sum()
    band_rho = float(np.sum(rho_b * aw) / area_tot) if area_tot > 0 else 0.0
    w_dens   = rho_b * aw
    wsum     = w_dens.sum()
    band_vel = float(np.sum(vel_b * w_dens) / wsum) if wsum > 0 else 0.0
    band_flux = band_rho * band_vel
    return band_rho, band_vel, band_flux


def average_by_N(results):
    """Promedia las realizaciones de cada N. Devuelve dict N -> medias y std."""
    Ns = sorted(set(r['N'] for r in results))
    out = {}
    for N in Ns:
        group = [r for r in results if r['N'] == N]
        S = group[0]['S']
        rho = np.array([r['density'] for r in group])
        vel = np.array([r['velocity'] for r in group])
        flx = np.array([r['flux'] for r in group])
        out[N] = {
            'S':       S,
            'rho':     rho.mean(axis=0), 'rho_std': rho.std(axis=0),
            'vel':     vel.mean(axis=0), 'vel_std': vel.std(axis=0),
            'flux':    flx.mean(axis=0), 'flux_std': flx.std(axis=0),
            'n_runs':  len(group),
        }
    return out


# ──────────────────────────────────────────────────────────────
# (A) Perfil radial: rho, |v| y Jin en función de S (3 ejes Y)
# ──────────────────────────────────────────────────────────────

def plot_profile_S(by_N, profile_N, out_prefix):
    """
    3 curvas en una misma figura (con 3 ejes Y) en función de S.
    profile_N: un N concreto, o 'all' para promediar todos los N disponibles.
    """
    Ns = sorted(by_N.keys())
    if profile_N == 'all':
        S   = by_N[Ns[0]]['S']
        rho = np.mean([by_N[N]['rho']  for N in Ns], axis=0)
        vel = np.mean([by_N[N]['vel']  for N in Ns], axis=0)
        flux = rho * vel        # Jin = <rho_in>·|<v_in>|
        label_N = f"promedio N={Ns[0]}–{Ns[-1]}"
    else:
        N = int(profile_N)
        S, rho, vel = by_N[N]['S'], by_N[N]['rho'], by_N[N]['vel']
        flux = rho * vel
        label_N = f"N={N}"

    fig, ax_rho = plt.subplots(figsize=(9, 5.5))
    ax_vel = ax_rho.twinx()
    ax_flx = ax_rho.twinx()
    ax_flx.spines['right'].set_position(('outward', 60))  # tercer eje desplazado

    c_rho, c_vel, c_flx = '#457b9d', '#2a9d8f', '#e63946'

    l1, = ax_rho.plot(S, rho,  color=c_rho, lw=2,
                      label=r"$\langle\rho_{in}\rangle$(S)")
    l2, = ax_vel.plot(S, vel,  color=c_vel, lw=2, ls='--',
                      label=r"$|\langle v_{in}\rangle|$(S)")
    l3, = ax_flx.plot(S, flux, color=c_flx, lw=2.2,
                      label=r"$J_{in}$(S) = $\langle\rho_{in}\rangle\,|\langle v_{in}\rangle|$")

    ax_rho.set_xlabel("S  [m]")
    ax_rho.set_ylabel(r"$\langle\rho_{in}\rangle$  [m$^{-2}$]", color=c_rho)
    ax_vel.set_ylabel(r"$|\langle v_{in}\rangle|$  [m/s]",      color=c_vel)
    ax_flx.set_ylabel(r"$J_{in}$  [m$^{-2}$s$^{-1}$]",          color=c_flx)
    for ax, c in [(ax_rho, c_rho), (ax_vel, c_vel), (ax_flx, c_flx)]:
        ax.tick_params(axis='y', colors=c)

    ax_rho.set_title(f"Perfiles radiales frescas-in vs S  ({label_N})")
    ax_rho.grid(True, alpha=0.25, ls=':')
    ax_rho.legend(handles=[l1, l2, l3], loc='upper right', fontsize=9)

    out = f"{out_prefix}_perfil_S.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Figura guardada: {out}")


def plot_flux_detail(by_N, S_max, out_prefix):
    """
    Detalle del flujo entrante Jin(S) en 0 ≤ S ≤ S_max, una curva por N
    (promedio de realizaciones), coloreadas con un colormap y leyenda con N.
    """
    Ns = sorted(by_N.keys())
    cmap = plt.cm.viridis
    colors = [cmap(i / max(len(Ns) - 1, 1)) for i in range(len(Ns))]

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for color, N in zip(colors, Ns):
        S, flux = by_N[N]['S'], by_N[N]['flux']
        m = S <= S_max
        ax.plot(S[m], flux[m], '-o', color=color, lw=1.5, ms=3.5, label=f"N={N}")

    ax.set_xlim(0, S_max)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("S  [m]")
    ax.set_ylabel(r"$J_{in}(S)$  [1/(m·s)]")
    ax.set_title(f"Flujo entrante: detalle 0 ≤ S ≤ {S_max:g} m")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, loc='upper left')

    out = f"{out_prefix}_flux_detalle_S.png"
    fig.tight_layout()
    fig.savefig(out, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Figura guardada: {out}")


# ──────────────────────────────────────────────────────────────
# (B) Observables en la banda S_lo ≤ S ≤ S_hi en función de N
# ──────────────────────────────────────────────────────────────

def band_by_N(results, S_lo, S_hi):
    """
    Para cada N agrega cada corrida sobre la banda [S_lo, S_hi] y luego promedia
    las realizaciones. Devuelve N_vals y, por observable, (mean, std) entre runs.
    """
    Ns = sorted(set(r['N'] for r in results))
    rho_m, rho_e, vel_m, vel_e, flux_m, flux_e = ([] for _ in range(6))
    for N in Ns:
        group = [r for r in results if r['N'] == N]
        bands = np.array([band_observables(r, S_lo, S_hi) for r in group])  # (runs, 3)
        rho_m.append(bands[:, 0].mean()); rho_e.append(bands[:, 0].std())
        vel_m.append(bands[:, 1].mean()); vel_e.append(bands[:, 1].std())
        flux_m.append(bands[:, 2].mean()); flux_e.append(bands[:, 2].std())
    return (np.array(Ns),
            np.array(rho_m), np.array(rho_e),
            np.array(vel_m), np.array(vel_e),
            np.array(flux_m), np.array(flux_e))


def plot_vsN(results, S_lo, S_hi, out_prefix):
    N_vals, rho, rho_e, vel, vel_e, flux, flux_e = band_by_N(results, S_lo, S_hi)
    band_lbl = f"{S_lo:g}–{S_hi:g} m"

    # ── Densidad y velocidad: misma figura, doble eje Y ──
    fig, ax_rho = plt.subplots(figsize=(8, 5))
    ax_vel = ax_rho.twinx()
    c_rho, c_vel = '#1f77b4', '#ff7f0e'

    l1 = ax_rho.errorbar(N_vals, rho, yerr=rho_e, fmt='o-', color=c_rho,
                         capsize=4, lw=2, label=r"$\langle\rho_{in}^{f}\rangle$")
    l2 = ax_vel.errorbar(N_vals, vel, yerr=vel_e, fmt='s-', color=c_vel,
                         capsize=4, lw=2, label=r"$|\langle v_{in}^{f}\rangle|$")

    ax_rho.set_xlabel("N")
    ax_rho.set_ylabel(rf"$\langle\rho_{{in}}^{{f}}\rangle_{{{band_lbl}}}$  [1/m$^2$]", color=c_rho)
    ax_vel.set_ylabel(rf"$|\langle v_{{in}}^{{f}}\rangle|_{{{band_lbl}}}$  [m/s]",      color=c_vel)
    ax_rho.tick_params(axis='y', colors=c_rho)
    ax_vel.tick_params(axis='y', colors=c_vel)
    ax_rho.set_title(f"Banda {band_lbl}  vs  N: densidad y velocidad")
    ax_rho.grid(True, alpha=0.25, ls=':')
    ax_rho.legend(handles=[l1, l2], loc='best', fontsize=9)

    out1 = f"{out_prefix}_vsN_rho_v.png"
    fig.tight_layout()
    fig.savefig(out1, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Figura guardada: {out1}")

    # ── Flujo Jin: figura aparte ──
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.errorbar(N_vals, flux, yerr=flux_e, fmt='D-', color='#e63946',
                capsize=4, lw=2)
    ax.set_xlabel("N")
    ax.set_ylabel(rf"$J_{{in}}\,|_{{{band_lbl}}}$  [m$^{{-2}}$ s$^{{-1}}$]")
    ax.set_title(f"Flujo $J_{{in}}$ en la banda {band_lbl}  vs  N")
    ax.grid(True, alpha=0.25, ls=':')

    out2 = f"{out_prefix}_vsN_jin.png"
    fig.tight_layout()
    fig.savefig(out2, dpi=150, bbox_inches='tight')
    plt.close(fig)
    print(f"Figura guardada: {out2}")


# ──────────────────────────────────────────────────────────────
# Entrypoint
# ──────────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(description="Gráficos de flujo frescas-in (Sistema 1)")
    p.add_argument("files", nargs='+', help="Archivos output.txt")
    p.add_argument("--N",   type=int, default=None, help="N único para todos los archivos")
    p.add_argument("--Ns",  type=str, default=None, help="N por archivo, comma-separated")
    p.add_argument("--r_outer",   type=float, default=40.0)
    p.add_argument("--r_inner",   type=float, default=1.0)
    p.add_argument("--radius",    type=float, default=1.0)
    p.add_argument("--dS",        type=float, default=0.2)
    p.add_argument("--tol",       type=float, default=0.5)
    p.add_argument("--S_lo",      type=float, default=2.0, help="Borde inferior de la banda para figura vs N")
    p.add_argument("--S_hi",      type=float, default=3.0, help="Borde superior de la banda para figura vs N")
    p.add_argument("--profile_N", type=str,   default="all",
                   help="N para el perfil vs S ('all' = promedio de todos)")
    p.add_argument("--S_detail",  type=float, default=5.0,
                   help="S máximo para el detalle de Jin(S) por N")
    p.add_argument("--out",       type=str,   default="outputs/flux")
    args = p.parse_args()

    if args.Ns:
        Ns_list = [int(x) for x in args.Ns.split(',')]
        assert len(Ns_list) == len(args.files), "Debe haber un N por archivo"
    elif args.N:
        Ns_list = [args.N] * len(args.files)
    else:
        Ns_list = []
        for f in args.files:
            t, s = parse_output(f)
            Ns_list.append(s.shape[1] if len(s) > 0 else 0)

    print(f"[INFO] {len(args.files)} archivos | banda=[{args.S_lo},{args.S_hi}] | profile_N={args.profile_N}")
    results = build_results(args.files, Ns_list, args.r_inner, args.r_outer,
                            args.radius, args.dS, args.tol)
    if not results:
        print("No se procesó ningún archivo.")
        return

    by_N = average_by_N(results)
    plot_profile_S(by_N, args.profile_N, args.out)
    plot_flux_detail(by_N, args.S_detail, args.out)
    plot_vsN(results, args.S_lo, args.S_hi, args.out)


if __name__ == "__main__":
    main()
