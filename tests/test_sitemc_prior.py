#!/usr/bin/env python3
"""Gates for the arrival-prior machinery (numpy analysis-side)."""
import math

import pytest

np = pytest.importorskip("numpy")

from deorbit.sitemc import arrival_prior as ap


def test_torus_tisserand():
    """Every torus hit's |U| must satisfy the Tisserand relation for its
    (a,e,i) — the whole encounter-kinematics chain in one identity."""
    U, w, aei = ap.torus_sample(2_000_000, ap.model_density("flat"), seed=3)
    assert U.shape[1] > 1000
    Umag = np.linalg.norm(U, axis=0)
    Ut = ap.tisserand_U(*aei)
    err = np.abs(Umag - Ut) / np.maximum(Ut, 0.3)
    assert np.median(err) < 0.03


def test_low_U_hits_are_earth_like():
    U, w, aei = ap.torus_sample(2_000_000, ap.model_density("flat"), seed=4)
    Umag = np.linalg.norm(U, axis=0)
    sel = Umag < 5.0
    a, e, i = aei[0][sel], aei[1][sel], aei[2][sel]
    # Tisserand at low U pins sqrt(p)*cos i = (3 - 1/a - U^2)/2, which for
    # the box a in [0.8,1.3] and U<5 km/s spans [0.86, 1.12] (NOT a near 1)
    spci = np.sqrt(a * (1 - e * e)) * np.cos(np.radians(i))
    assert spci.min() > 0.84 and spci.max() < 1.13
    assert i.max() <= 15.0


def test_earth_frame():
    """Analytic-ephemeris Earth frame: right speed, near-orthogonal basis,
    z-axis near the ecliptic pole (23.4 deg from equatorial z)."""
    for jd in (2451545.0, 2455000.25, 2461000.7):
        xh, yh, zh, v = ap.earth_frame(jd)
        assert 29.0e3 < v < 30.6e3
        assert abs(sum(a * b for a, b in zip(xh, yh))) < 0.02
        # ecliptic pole in equatorial frame: (0, -sin eps, cos eps)
        eps = math.radians(23.43929111)
        dot = -zh[1] * math.sin(eps) + zh[2] * math.cos(eps)
        assert dot > 0.999


def test_direction_structure():
    """Low-U radiants are NOT a thin ecliptic cone: vertical-dominated
    (U_z) encounters are a real sub-population (the i-dominated branch)."""
    U, w, aei = ap.torus_sample(4_000_000, ap.model_density("tilt"), seed=5)
    Umag = np.linalg.norm(U, axis=0)
    sel = (Umag > 0.3) & (Umag < 5.0)
    frac_vertical = float((np.abs(U[2][sel]) / Umag[sel] > 0.7).mean())
    assert 0.05 < frac_vertical < 0.9
