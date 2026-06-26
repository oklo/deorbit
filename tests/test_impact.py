"""Bounce-off scout regime classifier (pure stdlib). Locks the regime topology +
the two anchor cases that test whether the (theta, v/v_esc) framing is right."""
from deorbit.impact import scout


def test_anchor_regimes():
    # Cayambe (Earth, v=7.4 km/s, theta=2 deg): intact suborbital ricochet
    assert scout.regime(2.0, 7400.0 / scout.V_ESC["earth"]) == "I"
    # Orcus (Mars, v=10 km/s, theta=30 deg): conventional crater
    assert scout.regime(30.0, 10000.0 / scout.V_ESC["mars"], scout.V_ESC["mars"]) == "C"


def test_escape_requires_grazing_and_high_speed():
    assert scout.regime(2.0, 0.6) in ("I", "D")   # sub-escape grazing cannot leave
    assert scout.regime(2.0, 2.2) == "E"           # grazing + v >> v_esc -> escape


def test_steep_is_always_crater():
    assert scout.regime(40.0, 0.5) == "C"
    assert scout.regime(40.0, 2.2) == "C"


def test_intact_to_dispersed_at_melt_threshold():
    # increasing v/v_esc at fixed grazing angle: intact -> dispersed (contact melts iron)
    assert scout.regime(2.0, 0.5) == "I"
    assert scout.regime(2.0, 1.1) == "D"
