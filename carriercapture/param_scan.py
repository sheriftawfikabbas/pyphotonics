"""Parameter scan helpers (ported from ParamScan.jl)."""

from __future__ import annotations

import math
from pathlib import Path
from typing import Tuple

import numpy as np

from .constants import AMU_EV_C2, HBARC_EV_M
from .potential import Potential, harmonic, morse, solve_pot
from .capture_rate import ConfCoord, calc_overlap, calc_capt_coeff


def get_harmonic_q_m(hw: float, E0: float) -> float:
    """ΔQ for activationless Marcus regime for two harmonic PES."""
    a = AMU_EV_C2 / 2.0 * (hw / HBARC_EV_M / 1e10) ** 2
    return math.sqrt(E0 / a)


def get_morse_q_m(a: float, b: float, E0: float) -> float:
    """ΔQ for activationless Marcus regime for two Morse PES."""
    # Julia: Q_m = (1/b)*log(1-(E0/a)^0.5)
    return (1.0 / b) * math.log(1.0 - math.sqrt(E0 / a))


def fit_harmonic_params(hw_i: float, hw_f: float, dQ: float, dE: float) -> Tuple[Potential, Potential]:
    """Construct and solve two harmonic potentials (excited/relaxed)."""
    nev_excited = 180
    nf_min = (1.0 / hw_f) * (((nev_excited + 0.5) * hw_i) + 2.0 * dE) - 0.5
    nev_relaxed = int(math.floor(nf_min + 1.0))

    Q = np.linspace(-20.0 - dQ, 20.0 + dQ, 5000)

    potf = Potential(name="Relaxed state")
    potf.Q0 = dQ
    potf.E0 = 0.0
    potf.nev = nev_relaxed
    potf.func = lambda x: harmonic(x, hw_f, E0=potf.E0, Q0=potf.Q0)
    potf.Q = Q
    potf.E = potf.func(Q)
    solve_pot(potf)

    poti = Potential(name="Excited state")
    poti.Q0 = 0.0
    poti.E0 = dE
    poti.nev = nev_excited
    poti.func = lambda x: harmonic(x, hw_i, E0=poti.E0, Q0=poti.Q0)
    poti.Q = Q
    poti.E = poti.func(Q)
    solve_pot(poti)

    return poti, potf


def fit_morse_params(
    a_i: float, a_f: float, b_i: float, b_f: float, Q0: float, E0: float,
    *, output_dir: str | Path | None = None,
) -> Tuple[Potential, Potential]:
    """Construct and solve two Morse potentials (excited/relaxed).

    Set ``output_dir`` to save potential curves and all eigenvalues as CSV,
    full arrays (including wavefunctions) and input parameters as a compressed
    NumPy archive, and a PNG of the potentials and first five energy levels.
    The directory is created if needed; repeated runs replace its result files.
    With ``output_dir=None`` (the default), no files are written.
    Returns ``(poti, potf)`` in either case. Energies are in eV and Q is in
    sqrt(amu) Angstrom.
    """
    nev_excited = 180
    nev_relaxed = 800

    Q = np.linspace(-20.0 - Q0, 20.0 + Q0, 5000)

    poti = Potential(name="Excited state")
    poti.Q0 = 0.0
    poti.E0 = E0
    poti.nev = nev_excited
    poti.func = lambda x: morse(x, np.array([a_i, b_i]), E0=poti.E0, Q0=poti.Q0)
    poti.Q = Q
    poti.E = poti.func(Q)
    solve_pot(poti)

    potf = Potential(name="Relaxed state")
    potf.Q0 = Q0
    potf.E0 = 0.0
    potf.nev = nev_relaxed
    potf.func = lambda x: morse(x, np.array([a_f, b_f]), E0=potf.E0, Q0=potf.Q0)
    potf.Q = Q
    potf.E = potf.func(Q)
    solve_pot(potf)

    if output_dir is not None:
        _save_morse_results(
            poti, potf, Path(output_dir),
            dict(a_i=a_i, a_f=a_f, b_i=b_i, b_f=b_f, Q0=Q0, E0=E0),
        )

    return poti, potf


def _save_morse_results(poti: Potential, potf: Potential, output_dir: Path, params: dict) -> None:
    # Use a noninteractive canvas without changing the caller's pyplot backend.
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    output_dir.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        output_dir / "potentials.csv", np.column_stack((poti.Q, poti.E, potf.E)),
        delimiter=",", header="Q_sqrt_amu_A,initial_energy_eV,final_energy_eV", comments="",
    )
    for state, pot in (("initial", poti), ("final", potf)):
        np.savetxt(
            output_dir / f"eigenvalues_{state}.csv",
            np.column_stack((np.arange(pot.eps.size), pot.eps)),
            delimiter=",", header="state_index,energy_eV", comments="", fmt=("%d", "%.18e"),
        )
    np.savez_compressed(
        output_dir / "morse_results.npz",
        Q=poti.Q, E_i=poti.E, E_f=potf.E,
        eps_i=poti.eps, eps_f=potf.eps, chi_i=poti.chi, chi_f=potf.chi,
        **params,
    )

    fig = Figure(figsize=(11, 4.5), layout="constrained")
    FigureCanvasAgg(fig)
    ax_pot, ax_levels = fig.subplots(1, 2)
    for pot, color in ((poti, "tab:blue"), (potf, "tab:orange")):
        ax_pot.plot(pot.Q, pot.E, color=color, label=pot.name)
        levels = pot.eps[:5]
        ax_levels.plot(np.arange(levels.size), levels, "o-", color=color, label=pot.name)

    # Focus on the wells: the full solver grid includes an exponentially high
    # wall, which would make the low-energy structure invisible with autoscale.
    width_i = 1.0 / max(abs(params["b_i"]), 1e-12)
    width_f = 1.0 / max(abs(params["b_f"]), 1e-12)
    ax_pot.set_xlim(
        max(poti.Q.min(), min(poti.Q0 - 3 * width_i, potf.Q0 - 3 * width_f)),
        min(poti.Q.max(), max(poti.Q0 + 3 * width_i, potf.Q0 + 3 * width_f)),
    )
    low = min(poti.E0, potf.E0)
    high = max(poti.E0 + params["a_i"], potf.E0 + params["a_f"])
    margin = 0.1 * max(high - low, 0.01)
    ax_pot.set_ylim(low - margin, high + margin)
    ax_pot.set(title="Morse potentials", xlabel=r"Q ($\sqrt{\mathrm{amu}}$ Å)", ylabel="Energy (eV)")
    ax_levels.set(title="First five eigenvalues", xlabel="State index (n)", ylabel="Energy (eV)")
    ax_levels.set_xticks(range(5))
    for ax in (ax_pot, ax_levels):
        ax.legend()
        ax.grid(alpha=0.25)
    fig.savefig(output_dir / "morse_results.png", dpi=200)


def get_harmonic_capture(hw_i: float, hw_f: float, dQ: float, dE: float, W: float, V: float, Tmin: float = 100.0, Tmax: float = 800.0, NT: int = 100):
    poti, potf = fit_harmonic_params(hw_i, hw_f, dQ, dE)
    cc = ConfCoord(V1=poti, V2=potf, W=W, g=1)
    # Choose Q0 at midpoint as in Julia workflow for coupling matrix element reference
    Q0mid = float(poti.Q[len(poti.Q)//2])
    calc_overlap(cc, Q0=Q0mid)
    temps = np.linspace(Tmin, Tmax, NT)
    calc_capt_coeff(cc, V=V, temperature=temps)
    return cc


def get_morse_capture(a_i: float, a_f: float, b_i: float, b_f: float, Q0: float, E0: float, W: float, V: float, Tmin: float = 100.0, Tmax: float = 800.0, NT: int = 100):
    poti, potf = fit_morse_params(a_i, a_f, b_i, b_f, Q0, E0)
    cc = ConfCoord(V1=poti, V2=potf, W=W, g=1)
    Q0mid = float(poti.Q[len(poti.Q)//2])
    calc_overlap(cc, Q0=Q0mid)
    temps = np.linspace(Tmin, Tmax, NT)
    calc_capt_coeff(cc, V=V, temperature=temps)
    return cc
