"""Tests for the carriercapture package.

Kept fast (small grids / few eigenstates) so the suite runs in seconds while
still exercising the full pipeline: Brooglie eigensolver -> Potential fitting
-> ConfCoord capture-coefficient calculation -> TransferCoord rates.
"""

import math

import numpy as np
import pytest

from carriercapture.brooglie import solve as brooglie_solve
from carriercapture.potential import (
    Potential,
    harmonic,
    morse,
    fit_pot,
    filter_sample_points,
    find_crossing,
    cleave_pot,
    solve_pot,
)
from carriercapture.capture_rate import ConfCoord, calc_overlap, calc_capt_coeff
from carriercapture.transfer_coord import (
    TransferCoord,
    get_coupling,
    get_reorg_energy,
    get_activation_energy,
    get_transfer_rate,
)


def test_brooglie_harmonic_oscillator_spectrum():
    # Quantum harmonic oscillator (m=1, hbar=1): E_n = n + 1/2.
    energies, wfs = brooglie_solve(lambda x: 0.5 * x**2, N=60, a=-6, b=6, m=1.0, nev=5)
    expected = np.array([0.5, 1.5, 2.5, 3.5, 4.5])
    assert np.allclose(energies, expected, atol=0.05)
    assert len(wfs) == 5
    assert wfs[0].shape == (60,)


def make_harmonic_potential(name, hw, E0, Q0, nev, npts=300, half_width=5.0):
    pot = Potential(name=name)
    pot.Q0 = Q0
    pot.E0 = E0
    pot.nev = nev
    pot.func = lambda x: harmonic(x, hw, E0=pot.E0, Q0=pot.Q0)
    pot.Q = np.linspace(Q0 - half_width, Q0 + half_width, npts)
    pot.E = pot.func(pot.Q)
    return pot


def test_solve_pot_and_capture_rate_pipeline():
    potf = make_harmonic_potential("ground", hw=0.03, E0=0.0, Q0=0.0, nev=40)
    solve_pot(potf)
    poti = make_harmonic_potential("excited", hw=0.03, E0=1.0, Q0=1.0, nev=40)
    solve_pot(poti)

    assert np.all(np.diff(potf.eps) > 0)  # eigenvalues sorted ascending
    assert potf.chi.shape == (40, potf.Q.size)

    cc = ConfCoord(V1=poti, V2=potf, W=1e-3, g=1)
    calc_overlap(cc, Q0=0.0)
    assert cc.overlap_matrix.shape == (40, 40)

    temps = np.linspace(100, 800, 10)
    calc_capt_coeff(cc, V=1e-21, temperature=temps)

    assert cc.capt_coeff.shape == (10,)
    assert np.all(np.isfinite(cc.capt_coeff))
    assert np.all(cc.capt_coeff >= 0)


def test_calc_capt_coeff_rejects_insufficient_states():
    # Too few eigenstates means the highest one is still significantly
    # occupied at T_max, which the sanity check should catch.
    potf = make_harmonic_potential("ground", hw=0.03, E0=0.0, Q0=0.0, nev=5)
    solve_pot(potf)
    poti = make_harmonic_potential("excited", hw=0.03, E0=1.0, Q0=1.0, nev=5)
    solve_pot(poti)

    cc = ConfCoord(V1=poti, V2=potf, W=1e-3, g=1)
    calc_overlap(cc, Q0=0.0)
    with pytest.raises(AssertionError):
        calc_capt_coeff(cc, V=1e-21, temperature=np.linspace(100, 800, 10))


def test_fit_pot_harmonic_fittable_recovers_curvature():
    Q_data = np.linspace(-2, 2, 50)
    true_a = 0.017
    E_data = true_a * Q_data**2

    pot = Potential()
    pot.Q_data = Q_data
    pot.E_data = E_data
    pot.E0 = float(E_data.min())
    pot.Q0 = float(Q_data[np.argmin(E_data)])
    pot.func_type = "harmonic_fittable"
    pot.Q = np.linspace(-2, 2, 200)

    fit_pot(pot)

    # pot.E0/pot.Q0 are pinned to the nearest grid point rather than the
    # exact analytic vertex, so a small residual is expected.
    assert pot.params["fit_param"][0] == pytest.approx(true_a, rel=1e-2)


def test_fit_pot_morse_converges_from_default_initial_guess():
    # Regression test: fit_pot's default initial guess for the Morse model
    # used to be all zeros, which is a degenerate stationary point of the
    # Morse functional (both the term and its gradient vanish at a=0), so
    # curve_fit got stuck there and silently returned [0, 0].
    # fit_pot only fits points within 2 eV of the minimum (E_CUT), so keep
    # the whole test domain inside that window rather than the full,
    # steeply-diverging tail of the Morse curve.
    Q_data = np.linspace(0.2, 3, 60)
    true_A, true_a = 0.5, 1.2
    E_data = true_A * (1 - np.exp(-true_a * (Q_data - 1))) ** 2
    assert E_data.max() < 2.0

    pot = Potential()
    pot.Q_data = Q_data
    pot.E_data = E_data
    pot.E0 = float(E_data.min())
    pot.Q0 = float(Q_data[np.argmin(E_data)])
    pot.func_type = "morse"
    pot.Q = Q_data

    fit_pot(pot)

    fitted_A, fitted_a = pot.params["fit_param"]
    assert fitted_A == pytest.approx(true_A, rel=0.1)
    assert fitted_a == pytest.approx(true_a, rel=0.1)
    assert np.max(np.abs(pot.func(Q_data) - E_data)) < 0.05


def test_fit_pot_morse_poly_fits_data():
    Q_data = np.linspace(-1, 3, 60)
    E_data = 0.5 * (1 - np.exp(-1.2 * (Q_data - 1))) ** 2 + 0.01 * (Q_data - 1) ** 2

    pot = Potential()
    pot.Q_data = Q_data
    pot.E_data = E_data
    pot.E0 = float(E_data.min())
    pot.Q0 = float(Q_data[np.argmin(E_data)])
    pot.func_type = "morse_poly"
    pot.params = {"poly_order": 2}
    pot.Q = Q_data

    fit_pot(pot)

    relative_error = np.abs(pot.func(Q_data) - E_data) / np.maximum(np.abs(E_data), 1e-2)
    assert np.max(relative_error) < 0.02


def test_filter_sample_points_handles_tied_grid_minimum():
    # Regression test: with an even-length symmetric grid, the two points
    # nearest the vertex of a parabola can tie to floating-point precision,
    # which used to make filter_sample_points falsely raise "several minima"
    # for a perfectly ordinary single-well potential.
    Q_data = np.linspace(-2, 2, 50)
    E_data = 0.02 * Q_data**2

    pot = Potential()
    pot.Q_data = Q_data
    pot.E_data = E_data
    pot.E0 = float(E_data.min())
    pot.T = 300

    filter_sample_points(pot)

    assert pot.Q_data.size > 0
    assert pot.Q_data.size <= Q_data.size


def test_filter_sample_points_still_rejects_genuine_double_well():
    Q_data = np.linspace(-5, 5, 101)
    E_data = np.minimum(0.02 * (Q_data - 2) ** 2, 0.02 * (Q_data + 2) ** 2)

    pot = Potential()
    pot.Q_data = Q_data
    pot.E_data = E_data
    pot.E0 = float(E_data.min())
    pot.T = 300

    with pytest.raises(ValueError):
        filter_sample_points(pot)


def test_cleave_pot_handles_tied_grid_maximum():
    Q_data = np.linspace(-2, 2, 50)
    E_data = -0.02 * Q_data**2

    pot = Potential()
    pot.Q_data = Q_data
    pot.E_data = E_data

    p1, p2 = cleave_pot(pot)

    assert p1.Q_data.size + p2.Q_data.size == Q_data.size + 1  # apex shared


def test_find_crossing_and_transfer_rate():
    Q = np.linspace(-2, 2, 200)

    p1 = Potential()
    p1.Q0, p1.E0, p1.Q = -0.5, 0.0, Q
    p1.func = lambda x: 0.02 * (np.asarray(x) + 0.5) ** 2
    p1.E = p1.func(Q)

    p2 = Potential()
    p2.Q0, p2.E0, p2.Q = 0.5, 0.0, Q
    p2.func = lambda x: 0.02 * (np.asarray(x) - 0.5) ** 2
    p2.E = p2.func(Q)

    q0, e0 = find_crossing(p1, p2)
    assert q0 == pytest.approx(0.0, abs=1e-6)
    assert e0 == pytest.approx(0.005, abs=1e-6)

    pa = Potential()
    pa.T = 300
    pa.func = lambda x: 0.0 * np.asarray(x) - 0.01

    tc = TransferCoord(pot_d_i=p1, pot_d_f=p2, pot_a=pa)
    assert get_coupling(tc) == pytest.approx(0.015, abs=1e-6)
    assert get_reorg_energy(tc) == pytest.approx(0.02, abs=1e-6)
    assert get_activation_energy(tc) == pytest.approx(0.005, abs=1e-6)
    assert get_transfer_rate(tc) > 0
    assert math.isfinite(get_transfer_rate(tc))
