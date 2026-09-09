from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, Iterable, Optional, Tuple

import numpy as np
import scipy.interpolate as si
import scipy.optimize as so

from .constants import AMU_EV_C2, HBARC_EV_M, BOLTZ_EV_K
from . import brooglie


def _as_1d(a) -> np.ndarray:
    arr = np.asarray(a, dtype=float)
    return arr.reshape(-1)


@dataclass
class Potential:
    """1D Potential Energy Surface container (ported from Potential.jl).

    QE_data is stored as two 1D arrays (Q_data, E_data) of equal length.
    After fitting, the fitted function is stored in `func` and the evaluation grid
    (Q, E) is stored as arrays.

    Eigenvalues/eigenvectors after `solve_pot`:
      - eps : shape (nev,)
      - chi : shape (nev, len(Q))  (each row is an eigenfunction on the Q grid)
    """

    name: str = ""

    Q_data: np.ndarray = field(default_factory=lambda: np.array([0.0]))
    E_data: np.ndarray = field(default_factory=lambda: np.array([0.0]))

    E0: float = 0.0
    Q0: float = 0.0

    func_type: str = "harmonic_fittable"
    func: Callable[[np.ndarray], np.ndarray] = field(default_factory=lambda: (lambda x: np.zeros_like(np.asarray(x, dtype=float))))

    params: Dict[str, Any] = field(default_factory=dict)

    Q: np.ndarray = field(default_factory=lambda: np.array([]))
    E: np.ndarray = field(default_factory=lambda: np.array([]))

    nev: int = 0
    eps: np.ndarray = field(default_factory=lambda: np.array([]))
    chi: np.ndarray = field(default_factory=lambda: np.zeros((0, 0)))

    T: float = 293.15
    T_weight: bool = False

    def copy(self) -> "Potential":
        return Potential(
            name=self.name,
            Q_data=self.Q_data.copy(),
            E_data=self.E_data.copy(),
            E0=float(self.E0),
            Q0=float(self.Q0),
            func_type=self.func_type,
            func=self.func,
            params=dict(self.params),
            Q=self.Q.copy(),
            E=self.E.copy(),
            nev=int(self.nev),
            eps=self.eps.copy(),
            chi=self.chi.copy(),
            T=float(self.T),
            T_weight=bool(self.T_weight),
        )


# --- Basic potential functions (Julia parity) ---------------------------------

def harmonic(x, hw: float, *, E0: float, Q0: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    a = AMU_EV_C2 / 2.0 * (hw / HBARC_EV_M / 1e10) ** 2
    return a * (x - Q0) ** 2 + E0


def harmonic_fittable(x, coeff, *, E0: float, Q0: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    coeff = _as_1d(coeff)
    return (0.0 * x) + E0 + coeff[0] * (x - Q0) ** 2


def polyfunc(x, coeffs, *, E0: float, Q0: float, poly_order: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    coeffs = _as_1d(coeffs)
    y = (0.0 * x) + E0
    # Julia loops i=2:poly_order+1, uses coeffs[i]
    for i in range(2, poly_order + 2):
        y = y + coeffs[i - 1] * (x - Q0) ** (i - 1)
    return y


def morse(x, coeffs, *, E0: float, Q0: float) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    coeffs = _as_1d(coeffs)
    A = coeffs[0]
    a = coeffs[1]
    return (0.0 * x) + E0 + A * (1.0 - np.exp(-a * (x - Q0))) ** 2


def morse_poly(x, coeffs, *, E0: float, Q0: float, poly_order) -> np.ndarray:
    """Morse potential with polynomial corrections (ported from Potential.jl)."""
    x = np.asarray(x, dtype=float)
    coeffs = _as_1d(coeffs)

    if isinstance(poly_order, (int, np.integer)):
        orders = [int(poly_order)]
    else:
        # accept "2 4" etc
        orders = [int(t) for t in str(poly_order).split()]

    A = abs(coeffs[0])
    a = coeffs[1]
    r0 = coeffs[2]

    def morse_core(xx):
        return A * (np.exp(-2.0 * a * (xx - Q0 - r0)) - 2.0 * np.exp(-a * (xx - Q0 - r0)))

    # polynomial coefficients vector up to max order
    maxo = max(orders)
    poly_coeff = np.zeros(maxo + 1)
    for idx, o in enumerate(orders, start=0):
        poly_coeff[o] = coeffs[3 + idx]
    poly_coeff[-1] = abs(poly_coeff[-1])

    # Find r1 as the real root of derivative closest to zero
    # tot_derv = polyder(poly,1) - A*(-2a)*(exp(2a*r0) - exp(a*r0))
    # In numpy polynomial convention: np.poly1d expects highest-first; we'll use np.polynomial.Polynomial
    P = np.polynomial.Polynomial(poly_coeff)  # coeff[0] + coeff[1] x + ...
    dP = P.deriv(1)
    const = -A * (-2.0 * a) * (math.exp(2.0 * a * r0) - math.exp(a * r0))
    # Solve dP(x) - const = 0
    dP2 = dP - const
    roots = dP2.roots()
    real_roots = roots[np.isclose(roots.imag, 0.0)].real
    if real_roots.size == 0:
        r1 = 0.0
    else:
        r1 = real_roots[np.argmin(np.abs(real_roots))]

    return morse_core(x) + P(x - Q0 - r1) + E0 - morse_core(Q0) - P(-r1)


# --- Fitting helpers ----------------------------------------------------------

def get_spline(Qs, Es, *, weight=None, smoothness: float = -1.0, order: int = 2) -> Callable[[np.ndarray], np.ndarray]:
    Qs = _as_1d(Qs)
    Es = _as_1d(Es)
    if Qs[0] > Qs[-1]:
        Qs = Qs[::-1]
        Es = Es[::-1]
        if weight is not None:
            weight = _as_1d(weight)[::-1]

    # Julia uses Dierckx.Spline1D with bc="extrapolate"
    # UnivariateSpline extrapolates if ext=0.
    s = None
    if smoothness is not None and smoothness >= 0:
        s = float(smoothness)
    else:
        s = 0.0  # interpolate by default

    w = None if weight is None else _as_1d(weight)
    spl = si.UnivariateSpline(Qs, Es, w=w, k=int(order), s=s, ext=0)
    return lambda x: spl(_as_1d(x))


def get_bspline(Qs, Es) -> Callable[[np.ndarray], np.ndarray]:
    Qs = _as_1d(Qs)
    Es = _as_1d(Es)
    if Qs[0] > Qs[-1]:
        Qs = Qs[::-1]
        Es = Es[::-1]
    # Quadratic spline (k=2) similar to Julia BSpline Quadratic
    spl = si.make_interp_spline(Qs, Es, k=2, bc_type="natural")
    return lambda x: np.asarray(spl(_as_1d(x)), dtype=float)


# --- Core operations ----------------------------------------------------------

def pot_from_file(filename: str, resolution: int = 3001) -> Potential:
    """Parse a two-column file (Q, E) ignoring comment lines starting with '#'."""
    Q_dat = []
    E_dat = []
    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            Q_dat.append(float(parts[0]))
            E_dat.append(float(parts[1]))
    Q_dat = np.asarray(Q_dat, dtype=float)
    E_dat = np.asarray(E_dat, dtype=float)

    # build interpolation grid
    qmin, qmax = float(Q_dat.min()), float(Q_dat.max())
    Q = np.linspace(qmin, qmax, int(resolution))
    # include any original Q points not exactly in linspace
    Q = np.unique(np.concatenate([Q, Q_dat]))
    Q.sort()

    lin = si.interp1d(Q_dat, E_dat, kind="linear", fill_value="extrapolate", assume_sorted=False)
    E = np.asarray(lin(Q), dtype=float)

    pot = Potential()
    pot.Q_data = Q_dat
    pot.E_data = E_dat
    pot.Q = Q
    pot.E = E

    # min point
    min_idx = int(np.argmin(E_dat))
    pot.Q0 = float(Q_dat[min_idx])
    pot.E0 = float(E_dat[min_idx])
    return pot


def _extremum_index(E: np.ndarray, which: str) -> int:
    """Index of the (single) minimum/maximum of E.

    Grid points immediately adjacent to the true extremum can tie with it to
    floating-point precision (e.g. a symmetric potential sampled on a grid
    that straddles the vertex), so a plain `np.isclose(E, E.min())` mask can
    contain several *contiguous* indices for what is really one extremum.
    Only treat this as ambiguous ("several minima/maxima") when the close
    indices fall into more than one contiguous run.
    """
    target = E.min() if which == "min" else E.max()
    idx = np.where(np.isclose(E, target))[0]
    if idx.size == 0:
        idx = np.array([int(np.argmin(E) if which == "min" else np.argmax(E))])

    groups = [[idx[0]]]
    for i in idx[1:]:
        if i == groups[-1][-1] + 1:
            groups[-1].append(i)
        else:
            groups.append([i])

    if len(groups) != 1:
        raise ValueError(f"The fitted potential has several {'minima' if which == 'min' else 'maxima'}.")

    lo, hi = groups[0][0], groups[0][-1]
    sub = E[lo : hi + 1]
    offset = int(np.argmin(sub) if which == "min" else np.argmax(sub))
    return lo + offset


def filter_sample_points(pot: Potential, thermal_energy: Optional[float] = None) -> None:
    """Filter QE_data points to the connected 'island' around the minimum below threshold."""
    if thermal_energy is None:
        thermal_energy = pot.T * BOLTZ_EV_K

    Q = pot.Q_data
    E = pot.E_data

    min_ind = _extremum_index(E, "min")

    below = E <= (thermal_energy + pot.E0)

    island = []
    current = []
    for idx, ok in enumerate(below):
        if not ok:
            if current:
                if min_ind in current:
                    island = current
                    break
                current = []
        else:
            current.append(idx)

    if not island and current and (min_ind in current):
        island = current

    if not island:
        # fall back to keep all points
        island = list(range(len(Q)))

    pot.Q_data = Q[island]
    pot.E_data = E[island]


def fit_pot(pot: Potential) -> None:
    """Fit pot.func_type to (Q_data, E_data) and populate pot.func, pot.E."""
    if pot.Q.size == 0:
        # default to using Q_data as the evaluation grid
        pot.Q = pot.Q_data.copy()

    # energy cutoff window as in Julia
    E_CUT = 2.0
    e_cut_mask = pot.E_data < (E_CUT + pot.E0)

    Qd = pot.Q_data[e_cut_mask]
    Ed = pot.E_data[e_cut_mask]

    if pot.func_type == "bspline":
        f = get_bspline(pot.Q_data, pot.E_data)
        pot.func = lambda x: f(x)
        pot.E = pot.func(pot.Q)
        return

    if pot.func_type == "spline":
        weight = pot.params.get("weight", None)
        if isinstance(weight, str):
            weight = np.array([float(t) for t in weight.split()], dtype=float)
        smoothness = pot.params.get("smoothness", 0.0)
        order = pot.params.get("order", 2)
        f = get_spline(pot.Q_data, pot.E_data, weight=weight, smoothness=float(smoothness), order=int(order))
        pot.func = lambda x: f(x)
        pot.E = pot.func(pot.Q)
        return

    if pot.func_type == "harmonic":
        hw = float(pot.params["hw"])
        pot.func = lambda x: harmonic(x, hw, E0=pot.E0, Q0=pot.Q0)
        pot.E = pot.func(pot.Q)
        return

    # fittable types: polyfunc, morse_poly, morse, harmonic_fittable
    if pot.func_type == "polyfunc":
        poly_order = int(pot.params.get("poly_order", 4))
        p0 = np.asarray(pot.params.get("p0", np.zeros(poly_order + 1)), dtype=float)
        def model(x, *p):
            # p is coeffs vector
            coeffs = np.asarray(p, dtype=float)
            return polyfunc(x, coeffs, E0=pot.E0, Q0=pot.Q0, poly_order=poly_order)
        p0 = p0.reshape(-1)
        if p0.size < poly_order + 1:
            p0 = np.pad(p0, (0, poly_order + 1 - p0.size))
    elif pot.func_type == "morse_poly":
        poly_order = pot.params.get("poly_order", 4)
        # initial parameter length = 3 + len(orders)
        if isinstance(poly_order, (int, np.integer)):
            n_orders = 1
        else:
            n_orders = len(str(poly_order).split())
        # A=0 or a=0 is a degenerate stationary point of the Morse term (the
        # term and its gradient both vanish), so curve_fit can get stuck
        # there if not given a nonzero starting guess.
        A0 = max(float(Ed.max() - Ed.min()), 1e-3)
        a0 = 1.0 / max(float(Qd.max() - Qd.min()), 1e-3)
        default_p0 = np.concatenate([[A0, a0, 0.0], np.zeros(n_orders)])
        p0 = np.asarray(pot.params.get("p0", default_p0), dtype=float).reshape(-1)
        def model(x, *p):
            coeffs = np.asarray(p, dtype=float)
            return morse_poly(x, coeffs, E0=pot.E0, Q0=pot.Q0, poly_order=poly_order)
    elif pot.func_type == "morse":
        A0 = max(float(Ed.max() - Ed.min()), 1e-3)
        a0 = 1.0 / max(float(Qd.max() - Qd.min()), 1e-3)
        p0 = np.asarray(pot.params.get("p0", [A0, a0]), dtype=float).reshape(-1)
        def model(x, *p):
            coeffs = np.asarray(p, dtype=float)
            return morse(x, coeffs, E0=pot.E0, Q0=pot.Q0)
    elif pot.func_type == "harmonic_fittable":
        p0 = np.asarray(pot.params.get("p0", np.zeros(1)), dtype=float).reshape(-1)
        def model(x, *p):
            coeffs = np.asarray(p, dtype=float)
            return harmonic_fittable(x, coeffs, E0=pot.E0, Q0=pot.Q0)
    else:
        raise ValueError(f"Unknown func_type: {pot.func_type}")

    sigma = None
    if pot.T_weight:
        weights = np.exp(-(Ed - Ed.min()) / (pot.T * BOLTZ_EV_K))
        # curve_fit uses sigma; weights in Julia correspond to w in least squares.
        # We approximate by sigma = 1/sqrt(w).
        sigma = 1.0 / np.sqrt(np.maximum(weights, 1e-300))

    popt, _pcov = so.curve_fit(model, Qd, Ed, p0=p0, sigma=sigma, absolute_sigma=False, maxfev=20000)

    pot.func = lambda x: model(_as_1d(x), *popt)
    pot.E = pot.func(pot.Q)
    pot.params["fit_param"] = np.asarray(popt, dtype=float)


def solve1d_ev_amu(func: Callable[[np.ndarray], np.ndarray], Q: np.ndarray, *, nev: int = 30, maxiter: Optional[int] = None):
    """Solve 1D Schrödinger equation using the Brooglie backend.

    Mirrors Potential.jl:
        factor = (1/amu) * (ħc*1E10)^2
        Brooglie.solve(x -> func.(x*sqrt(factor)); a=Qi/sqrt(factor), b=Qf/sqrt(factor), m=1)
        return eps, chi / factor^(1/4)
    """
    Q = _as_1d(Q)
    NQ = Q.size
    Qi, Qf = float(Q.min()), float(Q.max())

    if maxiter is None:
        maxiter = int(nev * NQ)

    factor = (1.0 / AMU_EV_C2) * (HBARC_EV_M * 1e10) ** 2
    sqrtf = math.sqrt(factor)

    # V(x) with x in scaled coordinates
    def V_scaled(x):
        x = np.asarray(x, dtype=float)
        return func(x * sqrtf)

    eps, wfs = brooglie.solve(V_scaled, N=NQ, a=Qi / sqrtf, b=Qf / sqrtf, m=1.0, nev=nev, maxiter=maxiter)

    # Convert list of (NQ,) arrays into (nev, NQ) and scale as in Julia
    chi = np.vstack([wf.reshape(-1) for wf in wfs]) / (factor ** 0.25)
    return np.asarray(eps, dtype=float), chi


def solve_pot(pot: Potential) -> None:
    if pot.nev <= 0:
        raise ValueError("pot.nev must be set to a positive integer before solve_pot")
    pot.eps, pot.chi = solve1d_ev_amu(pot.func, pot.Q, nev=int(pot.nev))


def find_crossing(pot_1: Potential, pot_2: Potential) -> Tuple[float, float]:
    """Find Q where pot_1(Q) == pot_2(Q)."""
    Q = pot_1.Q
    if Q.size == 0:
        raise ValueError("pot_1.Q is empty; fit or set an evaluation grid first")

    def diff(x):
        return float(pot_1.func(np.array([x]))[0] - pot_2.func(np.array([x]))[0])

    # Try to bracket using sign changes on the grid
    vals = np.array([diff(float(q)) for q in Q], dtype=float)
    sgn = np.sign(vals)
    idx = np.where(sgn[:-1] * sgn[1:] <= 0)[0]
    if idx.size > 0:
        i = int(idx[0])
        a, b = float(Q[i]), float(Q[i + 1])
        res = so.root_scalar(diff, bracket=(a, b), method="brentq")
        q0 = float(res.root)
    else:
        # secant from midpoint
        mid = float(Q[len(Q) // 2])
        res = so.root_scalar(diff, x0=mid, x1=mid + 1e-3, method="secant")
        q0 = float(res.root)

    e0 = float(pot_1.func(np.array([q0]))[0])
    return q0, e0


def cleave_pot(pot: Potential, *, discard_max: bool = False) -> Tuple[Potential, Potential]:
    """Split a potential into two by its single apex (global maximum).

    Returns (left_pot, right_pot). If discard_max=True, the maximum point is excluded.
    """
    Q = pot.Q_data
    E = pot.E_data
    imax = _extremum_index(E, "max")

    if discard_max:
        left_idx = np.arange(0, imax)
        right_idx = np.arange(imax + 1, len(Q))
    else:
        left_idx = np.arange(0, imax + 1)
        right_idx = np.arange(imax, len(Q))

    p1 = pot.copy()
    p1.Q_data = Q[left_idx].copy()
    p1.E_data = E[left_idx].copy()
    p2 = pot.copy()
    p2.Q_data = Q[right_idx].copy()
    p2.E_data = E[right_idx].copy()
    return p1, p2
