"""Tests for foamgci.gci.

Anchor case: **Celik et al. (2008), J. Fluids Eng. 130, 078001**,
"Procedure for estimation and reporting of uncertainty due to
discretization in CFD applications", Table 1 (first column —
reattachment length on the Celik & Karatekin 1997 backward-facing-step
non-uniform structured grids of N = 18000 / 8000 / 4500 cells).

The paper uses NON-UNIFORM refinement ratios:
    r_21 = 1.5,   r_32 = 1.333
    φ_1 = 6.063,  φ_2 = 5.972,  φ_3 = 5.863     (φ_k = reattachment length on grid k)
Expected results reproduced in the paper's Table 1:
    p̂ ≈ 1.53
    φ_exact (extrapolated) ≈ 6.1685
    ea_21 (approx. rel. err.) ≈ 1.5 %
    eext_21 (extrapolated rel. err.) ≈ 1.7 %
    GCI_21^fine ≈ 2.2 %     (with safety factor 1.25)

Cell sizes are reconstructed from N_i and r_ij: in 2D, h ∝ 1/√N, and the
absolute scale is irrelevant — only the ratios matter — so we use
h_3 = 1.0 and walk inward.
"""
from __future__ import annotations

import math

import pytest

from foamgci.gci import (
    apparent_order,
    richardson_extrapolation,
    roache_gci,
    gci_over_hierarchy,
)


# ---- Celik 2008 Table 1, column 1 ---------------------------------------

CELIK_PHIs = dict(phi_fine=6.063, phi_medium=5.972, phi_coarse=5.863)
# r_21 = h_medium / h_fine = 1.5, r_32 = h_coarse / h_medium = 1.333
# Walk inward from h_coarse = 1.0:
CELIK_HS = dict(
    h_coarse=1.0,
    h_medium=1.0 / 1.333,
    h_fine=(1.0 / 1.333) / 1.5,
)


def test_celik_apparent_order() -> None:
    p = apparent_order(
        phi1=CELIK_PHIs["phi_fine"],
        phi2=CELIK_PHIs["phi_medium"],
        phi3=CELIK_PHIs["phi_coarse"],
        r21=1.5, r32=1.333,
    )
    # Paper Table 1, column 1 reports p ≈ 1.53.
    assert 1.50 < p < 1.56


def test_celik_richardson_extrapolation() -> None:
    p = apparent_order(
        CELIK_PHIs["phi_fine"], CELIK_PHIs["phi_medium"],
        CELIK_PHIs["phi_coarse"], 1.5, 1.333,
    )
    phi_e = richardson_extrapolation(
        CELIK_PHIs["phi_fine"], CELIK_PHIs["phi_medium"], 1.5, p
    )
    # Paper: φ_ext ≈ 6.1685
    assert math.isclose(phi_e, 6.1685, rel_tol=0.002)


def test_celik_full_gci() -> None:
    g = roache_gci(
        phi_coarse=CELIK_PHIs["phi_coarse"],
        phi_medium=CELIK_PHIs["phi_medium"],
        phi_fine=CELIK_PHIs["phi_fine"],
        h_coarse=CELIK_HS["h_coarse"],
        h_medium=CELIK_HS["h_medium"],
        h_fine=CELIK_HS["h_fine"],
    )
    assert 1.50 < g.p_apparent < 1.56
    # Paper Table 1: GCI_21^fine = 2.2% with Fs=1.25.
    assert math.isclose(g.gci_fine_21_pct, 2.2, rel_tol=0.05)
    # Approximate relative error e_a21 ≈ 1.5%.
    assert math.isclose(g.e21_relative * 100.0, 1.5, rel_tol=0.05)


# ---- Synthetic exact-second-order benchmark ------------------------------

def test_exact_second_order_recovers_p_equals_2() -> None:
    """If φ(h) = φ_exact + C h^2 with C constant, p̂ should be 2."""
    phi_exact = 10.0
    C = 0.1
    h1, h2, h3 = 1.0, 2.0, 4.0
    phi1 = phi_exact + C * h1 ** 2
    phi2 = phi_exact + C * h2 ** 2
    phi3 = phi_exact + C * h3 ** 2
    p = apparent_order(phi1, phi2, phi3, 2.0, 2.0)
    assert math.isclose(p, 2.0, abs_tol=1e-6)
    phi_e = richardson_extrapolation(phi1, phi2, 2.0, p)
    assert math.isclose(phi_e, phi_exact, abs_tol=1e-9)


def test_asymptotic_ratio_is_one_for_pure_order_p() -> None:
    """When the three grids ARE in the asymptotic range and φ(h) = φ_e + C h^p,
    the GCIs of overlapping triplets should be in ratio r_21^p, i.e.
    R_asymptotic = 1.0 exactly."""
    phi_exact = 5.0
    C = 0.05
    p_true = 2.0
    # 4 grids, r = 2
    hs = [4.0, 2.0, 1.0, 0.5]
    phis = [phi_exact + C * h ** p_true for h in hs]
    results = gci_over_hierarchy(phis, hs, labels=["C", "M", "F", "XF"])
    assert len(results) == 2
    for g in results:
        assert math.isclose(g.p_apparent, p_true, abs_tol=1e-6)
        assert math.isclose(g.asymptotic_ratio, 1.0, abs_tol=1e-6)


# ---- Defensive --------------------------------------------------------------

def test_gci_rejects_bad_ordering() -> None:
    with pytest.raises(ValueError, match="strictly decreasing|h_coarse > h_medium"):
        roache_gci(
            phi_coarse=1.0, phi_medium=1.0, phi_fine=1.0,
            h_coarse=1.0, h_medium=2.0, h_fine=4.0,
        )


def test_gci_rejects_stalled_refinement() -> None:
    with pytest.raises(ValueError, match="exact|stalled"):
        roache_gci(
            phi_coarse=2.0, phi_medium=2.0, phi_fine=2.0,
            h_coarse=4.0, h_medium=2.0, h_fine=1.0,
        )


def test_hierarchy_requires_three_grids() -> None:
    with pytest.raises(ValueError, match="at least 3"):
        gci_over_hierarchy([1.0, 2.0], [2.0, 1.0])
