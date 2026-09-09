"""A small sparse-grid Schrödinger solver (ported from Brooglie.jl).

Notes:
  - Energies are returned in the same units as the potential function V.
  - This implementation supports arbitrary dimension D, but the CarrierCapture
    workflow uses it in 1D.
"""

from __future__ import annotations

import inspect
import itertools
import math
from typing import Callable, Iterable, List, Sequence, Tuple

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla

H2EV = 27.21138602  # kept for parity with Julia, not used internally


def _number_of_arguments(f: Callable) -> int:
    """Infer the number of positional arguments a callable expects.

    For plain Python functions/lambdas with a standard signature this works well.
    If inference fails, we default to 1 (CarrierCapture uses 1D).
    """
    try:
        sig = inspect.signature(f)
        params = [
            p for p in sig.parameters.values()
            if p.kind in (p.POSITIONAL_ONLY, p.POSITIONAL_OR_KEYWORD)
        ]
        # If *args present, treat as unknown; default to 1
        if any(p.kind == p.VAR_POSITIONAL for p in sig.parameters.values()):
            return 1
        return len(params)
    except Exception:
        return 1


def _build_skeleton(D: int, N: int) -> sp.csr_matrix:
    """Build the Laplacian-like skeleton matrix used by Brooglie.buildH."""

    n = N ** D
    diags = []
    offsets = []

    for i in range(1, D + 1):
        step = N ** (i - 1)
        # motif = [-ones(step*(N-1)); zeros(step)]
        motif = np.concatenate([ -np.ones(step * (N - 1)), np.zeros(step) ])
        # cycle motif and take n - step elements
        tiled = np.tile(motif, int(math.ceil((n - step) / motif.size) + 1))
        diag = tiled[: (n - step)]
        diags.append(diag)
        offsets.append(step)

    data = []
    offs = []
    for off, d in zip(offsets, diags):
        data.append(d)
        offs.append(off)

    S = sp.diags(data, offs, shape=(n, n), format="csr")
    return S + S.T


def combinations(els: Sequence[float], N: int) -> List[List[float]]:
    """All length-N combinations with repetition of `els` (lexicographic)."""
    if N < 0:
        raise ValueError("N must be >= 0")
    if N == 0:
        return []
    return [list(t) for t in itertools.product(els, repeat=N)]


def build_h(V: Callable, *, N: int = 20, a: float = -1.0, b: float = 1.0, m: float = 1.0) -> sp.csr_matrix:
    """Hamiltonian matrix for a particle in a D-dimensional box [a,b]^D.

    Discretization uses an N-point grid per coordinate.
    """
    D = _number_of_arguments(V)
    xs_1d = np.linspace(a, b, N)
    grid = list(itertools.product(xs_1d, repeat=D))
    h = (b - a) / N

    def theta(zeta: Tuple[float, ...]) -> float:
        return float(V(*zeta)) * 2.0 * m * h * h

    diag = np.array([2.0 * D + theta(z) for z in grid], dtype=float)
    H = _build_skeleton(D, N) + sp.diags(diag, 0, format="csr")
    return H


def _integrate(phi: np.ndarray, L: float) -> float:
    sizes = phi.shape
    if len(set(sizes)) != 1:
        raise ValueError("Array has different widths in different dimensions.")
    D = phi.ndim
    N = sizes[0]
    dA = (L / N) ** D
    return float(np.sum(phi) * dA)


def _normalize_wf(phi: np.ndarray, L: float) -> np.ndarray:
    return phi / math.sqrt(_integrate(phi * phi, L))


def solve(
    V: Callable,
    *,
    N: int = 500,
    a: float = -1.0,
    b: float = 1.0,
    m: float = 1.0,
    nev: int | None = None,
    maxiter: int = 1000,
) -> Tuple[np.ndarray, List[np.ndarray]]:
    """Solve the stationary Schrödinger equation on a grid.

    Returns (energies, eigenfunctions), where:
      - energies has shape (nev,)
      - eigenfunctions is a list of D-dimensional arrays of shape (N,)*D
    """
    D = _number_of_arguments(V)
    H = build_h(V, N=N, a=a, b=b, m=m)

    n = H.shape[0]
    if nev is None:
        nev = max(1, N // 20)
    nev = int(nev)
    if nev >= n:
        nev = n - 1

    try:
        if D == 1:
            # The 1D Hamiltonian is tridiagonal. Direct bisection handles steep
            # Morse walls without the slow convergence of smallest-value ARPACK.
            # Use an absolute tolerance: the default scales with the largest
            # diagonal entry and can obscure the low-energy states entirely.
            vals, vecs = la.eigh_tridiagonal(
                H.diagonal(), H.diagonal(1), select="i",
                select_range=(0, nev - 1), lapack_driver="stebz", tol=1e-12,
            )
        else:
            vals, vecs = spla.eigsh(H, k=nev, which="SA", maxiter=maxiter)
    except Exception as ex:
        raise RuntimeError(f"Eigen-solve failed: {ex}") from ex

    h = (b - a) / N
    energies = vals / (2.0 * m * h * h)

    # reshape eigenvectors into D-d arrays and normalize
    wfs: List[np.ndarray] = []
    for i in range(nev):
        wf = vecs[:, i].reshape((N,) * D)
        wfs.append(_normalize_wf(wf, b - a))

    return energies, wfs
