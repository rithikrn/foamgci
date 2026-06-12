"""foamgci.gci — Roache Grid Convergence Index.

Implements the GCI formulation of

    Roache, P.J. (1994), "Perspective: A method for uniform reporting of grid
    refinement studies", J. Fluids Eng., 116, 405–413.

with the four-grid extension and asymptotic-range diagnostic described in

    Celik, I.B., Ghia, U., Roache, P.J., Freitas, C.J., Coleman, H., Raad,
    P.E. (2008), "Procedure for estimation and reporting of uncertainty due
    to discretization in CFD applications", J. Fluids Eng., 130, 078001.

A *triplet* consists of three grids with cell sizes :math:`h_1 < h_2 < h_3`
(fine, medium, coarse) and a converged scalar quantity of interest
:math:`\\phi_i` on each. Define the refinement ratios

.. math::
    r_{21} = h_2 / h_1, \\qquad r_{32} = h_3 / h_2.

For a refinement study with constant ratio :math:`r_{21}=r_{32}=r`, the
**apparent order** :math:`\\hat p` is

.. math::
    \\hat p = \\frac{1}{\\ln r}\\,
        \\Bigl|\\ln\\!\\bigl|(\\phi_3 - \\phi_2)/(\\phi_2 - \\phi_1)\\bigr|
        \\Bigr|.

Celik's iterative form handles non-constant :math:`r` and is implemented
in :func:`apparent_order`. The **Richardson-extrapolated** exact-grid
value is

.. math::
    \\phi_{\\text{exact}} \\approx \\phi_1 + \\frac{\\phi_1 - \\phi_2}{r_{21}^{\\hat p} - 1}.

The **GCI on the fine grid** with the standard safety factor :math:`F_s=1.25`
is

.. math::
    \\text{GCI}_{21}^{\\text{fine}} = F_s \\frac{|\\phi_1 - \\phi_2|}{|\\phi_1| (r_{21}^{\\hat p} - 1)}.

A pair of overlapping triplets (coarse–medium–fine, medium–fine–extra-fine)
should give GCIs that are in the **asymptotic ratio**

.. math::
    R = \\frac{r_{21}^{\\hat p} \\cdot \\text{GCI}^{\\text{fine}}_{21}}{\\text{GCI}^{\\text{medium}}_{32}}
        \\approx 1.

Values far from unity indicate that the grids are not in the asymptotic
range and the GCI uncertainty band should not be trusted at face value.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np


# ---------------------------------------------------------------------------
# Apparent order
# ---------------------------------------------------------------------------

def apparent_order(
    phi1: float, phi2: float, phi3: float,
    r21: float, r32: float,
    p_min: float = 0.1, p_max: float = 6.0, tol: float = 1e-8,
    max_iter: int = 100,
) -> float:
    """Iteratively solve for Celik et al. (2008) apparent order :math:`\\hat p`.

    The fixed-point iteration is

    .. math::
        p_{k+1} = \\frac{1}{\\ln r_{21}}\\,
        \\Bigl|\\ln|\\epsilon_{32}/\\epsilon_{21}| + q(p_k)\\Bigr|, \\\\
        q(p) = \\ln\\!\\Bigl(\\frac{r_{21}^{p} - s}{r_{32}^{p} - s}\\Bigr),
        \\qquad s = \\operatorname{sign}(\\epsilon_{32}/\\epsilon_{21}),

    where :math:`\\epsilon_{32}=\\phi_3-\\phi_2`,
    :math:`\\epsilon_{21}=\\phi_2-\\phi_1`. For a uniform refinement
    :math:`r_{21}=r_{32}\\Rightarrow q\\equiv 0` and the iteration
    collapses to the closed-form expression.

    Parameters
    ----------
    phi1, phi2, phi3 : float
        Quantity of interest on fine, medium, coarse grid (in that order).
    r21, r32 : float
        Refinement ratios :math:`h_2/h_1`, :math:`h_3/h_2`. Both > 1.
    p_min, p_max : float
        Apparent order is clipped to ``[p_min, p_max]`` after each step
        for numerical stability (Celik recommendation).

    Returns
    -------
    p : float
        Apparent order :math:`\\hat p`.

    Raises
    ------
    ValueError
        If consecutive differences are zero or refinement ratios ≤ 1.
    """
    if r21 <= 1.0 or r32 <= 1.0:
        raise ValueError("Refinement ratios must be > 1.")
    eps32 = phi3 - phi2
    eps21 = phi2 - phi1
    if eps21 == 0.0:
        raise ValueError("phi2 - phi1 = 0; refinement is exact (or stalled).")
    if eps32 == 0.0:
        raise ValueError("phi3 - phi2 = 0; coarsening is exact (or stalled).")

    s = float(np.sign(eps32 / eps21))
    p = 2.0  # benign starting guess
    for _ in range(max_iter):
        try:
            num = r21 ** p - s
            den = r32 ** p - s
            if num == 0.0 or den == 0.0:
                # Avoid log(0); fall back to constant-r limit.
                q = 0.0
            else:
                q = np.log(num / den)
        except (OverflowError, FloatingPointError):
            q = 0.0
        p_new = (1.0 / np.log(r21)) * np.abs(np.log(abs(eps32 / eps21)) + q)
        p_new = float(np.clip(p_new, p_min, p_max))
        if abs(p_new - p) < tol:
            return p_new
        p = p_new
    return p


# ---------------------------------------------------------------------------
# Richardson + GCI
# ---------------------------------------------------------------------------

def richardson_extrapolation(
    phi1: float, phi2: float, r21: float, p: float
) -> float:
    """Exact-grid value via Richardson extrapolation."""
    denom = r21 ** p - 1.0
    if denom == 0.0:
        return float(phi1)
    return float(phi1 + (phi1 - phi2) / denom)


@dataclass
class GCIResult:
    """Result of Roache GCI on a triplet (coarse → medium → fine)."""

    label_coarse: str
    label_medium: str
    label_fine: str
    h_coarse: float
    h_medium: float
    h_fine: float
    phi_coarse: float
    phi_medium: float
    phi_fine: float
    r21: float
    r32: float
    p_apparent: float
    phi_exact: float
    gci_fine_21_pct: float      # GCI_21 on the fine grid (%)
    gci_medium_32_pct: float    # GCI_32 on the medium grid (%)
    asymptotic_ratio: float     # ≈ 1 for asymptotic range
    Fs: float                   # safety factor (1.25 for ≥ 3 grids)
    e21_relative: float         # |φ1-φ2|/|φ1|
    e32_relative: float         # |φ2-φ3|/|φ2|
    regime: str = "monotonic"   # monotonic | oscillatory | divergent | degenerate
        note: str = ""

    def __repr__(self) -> str:  # pragma: no cover — printout
        return (
            f"GCIResult(triplet={self.label_coarse}/"
            f"{self.label_medium}/{self.label_fine}, "
            f"p̂={self.p_apparent:.3f}, "
            f"φ_exact={self.phi_exact:.6g}, "
            f"GCI_fine={self.gci_fine_21_pct:.3f}%, "
            f"R_asym={self.asymptotic_ratio:.3f})"
        )


def roache_gci(
    phi_coarse: float, phi_medium: float, phi_fine: float,
    h_coarse: float, h_medium: float, h_fine: float,
    label_coarse: str = "coarse",
    label_medium: str = "medium",
    label_fine: str = "fine",
    Fs: float = 1.25,
) -> GCIResult:
    """Roache GCI on a single refinement triplet.

    Parameters
    ----------
    phi_coarse, phi_medium, phi_fine : float
        Quantity of interest on each grid.
    h_coarse, h_medium, h_fine : float
        Representative cell size on each grid (units irrelevant — only
        ratios enter).
    label_* : str
        Used only for human-readable output.
    Fs : float
        Safety factor. Celik et al. (2008) recommend ``1.25`` when three
        or more grids are available, ``3.0`` otherwise.

    Returns
    -------
    GCIResult
    """
    if not (h_coarse > h_medium > h_fine > 0):
        raise ValueError(
            "Cell sizes must satisfy h_coarse > h_medium > h_fine > 0, "
            f"got h_c={h_coarse}, h_m={h_medium}, h_f={h_fine}."
        )
    r21 = h_medium / h_fine
    r32 = h_coarse / h_medium

    d21 = phi_fine - phi_medium
    d32 = phi_medium - phi_coarse
    eps21 = abs(d21)
    eps32 = abs(d32)
    e21_rel = eps21 / max(abs(phi_fine), 1e-300)
    e32_rel = eps32 / max(abs(phi_medium), 1e-300)

    if eps21 < 1e-14 or eps32 < 1e-14:
        regime, note = "degenerate", "vanishing residual; observed order undefined"
    else:
        R_conv = d21 / d32
        if R_conv < 0.0:
            regime = "oscillatory"
            note = f"R={R_conv:.3f} < 0: sign change between successive differences"
        elif R_conv > 1.0:
            regime = "divergent"
            note = f"R={R_conv:.3f} > 1: residual grows under refinement (pre-asymptotic)"
        else:
            regime, note = "monotonic", "monotonic convergence; GCI valid"

    if regime != "monotonic":
        nan = float("nan")
        return GCIResult(
            label_coarse=label_coarse, label_medium=label_medium,
            label_fine=label_fine,
            h_coarse=h_coarse, h_medium=h_medium, h_fine=h_fine,
            phi_coarse=phi_coarse, phi_medium=phi_medium, phi_fine=phi_fine,
            r21=r21, r32=r32,
            p_apparent=nan, phi_exact=nan,
            gci_fine_21_pct=nan, gci_medium_32_pct=nan, asymptotic_ratio=nan,
            Fs=Fs, e21_relative=e21_rel, e32_relative=e32_rel,
            regime=regime, note=note,
        )

    p = apparent_order(phi_fine, phi_medium, phi_coarse, r21, r32)
    phi_exact = richardson_extrapolation(phi_fine, phi_medium, r21, p)
    gci21 = Fs * e21_rel / (r21 ** p - 1.0)
    gci32 = Fs * e32_rel / (r32 ** p - 1.0)
    R = r21 ** p * eps21 / max(eps32, 1e-300)

     return GCIResult(
        label_coarse=label_coarse,
        label_medium=label_medium,
        label_fine=label_fine,
        h_coarse=h_coarse, h_medium=h_medium, h_fine=h_fine,
        phi_coarse=phi_coarse, phi_medium=phi_medium, phi_fine=phi_fine,
        r21=r21, r32=r32,
        p_apparent=p,
        phi_exact=phi_exact,
        gci_fine_21_pct=100.0 * gci21,
        gci_medium_32_pct=100.0 * gci32,
        asymptotic_ratio=R,
        Fs=Fs,
        e21_relative=e21_rel,
        e32_relative=e32_rel,
        regime="monotonic",
        note=note,
    )


def gci_over_hierarchy(
    phis: Sequence[float],
    hs: Sequence[float],
    labels: Sequence[str] | None = None,
    Fs: float = 1.25,
) -> list[GCIResult]:
    """Compute GCI on every consecutive triplet of a refinement hierarchy.

    Pass the grids **coarse-to-fine** (``hs`` strictly decreasing).
    For 4 grids :math:`(C, M, F, XF)` this returns two GCIs:
    one on the ``(C, M, F)`` triplet and one on the ``(M, F, XF)`` triplet.
    """
    if len(phis) != len(hs):
        raise ValueError("phis and hs must have the same length.")
    n = len(phis)
    if n < 3:
        raise ValueError("Need at least 3 grids.")
    if labels is None:
        labels = [f"grid{i}" for i in range(n)]
    if any(hs[i] <= hs[i + 1] for i in range(n - 1)):
        raise ValueError("hs must be strictly decreasing (coarse → fine).")

    out: list[GCIResult] = []
    for i in range(n - 2):
        out.append(
            roache_gci(
                phi_coarse=phis[i], phi_medium=phis[i + 1], phi_fine=phis[i + 2],
                h_coarse=hs[i], h_medium=hs[i + 1], h_fine=hs[i + 2],
                label_coarse=labels[i], label_medium=labels[i + 1],
                label_fine=labels[i + 2], Fs=Fs,
            )
        )
    return out
