"""
gci.py
------
Roache Grid Convergence Index per the formulation standardised in

  Celik IB, Ghia U, Roache PJ, Freitas CJ, Coleman H, Raad PE.
  Procedure for estimation and reporting of uncertainty due to
  discretization in CFD applications.  ASME J. Fluids Eng. 130 (2008).

Convention used throughout (Celik et al.):
    subscript 1 = finest, 3 = coarsest;  r = h_coarse / h_fine > 1.
    eps21 = phi_2 - phi_1     (medium minus fine)
    eps32 = phi_3 - phi_2     (coarse minus medium)
    R     = eps21 / eps32

Convergence regime (Celik et al. sec. 2.3, eq. 8):
    0 < R < 1  : monotonic convergence  (asymptotic; GCI is meaningful)
    R < 0      : oscillatory convergence
    R > 1      : divergence  (NOT in asymptotic range; GCI undefined)
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from math import log
from typing import Literal


Regime = Literal["monotonic", "oscillatory", "divergent", "degenerate"]


@dataclass
class GCIResult:
    label: str
    phi1: float
    phi2: float
    phi3: float
    h1: float
    h2: float
    h3: float
    r: float
    eps21: float
    eps32: float
    R: float
    regime: Regime
    p_obs: float
    phi_ext: float
    gci_fine_pct: float
    gci_med_pct: float
    asymptotic_ratio: float
    note: str

    def as_dict(self) -> dict:
        return asdict(self)


def _classify(R: float) -> Regime:
    if not (R == R):
        return "degenerate"
    if R < 0:
        return "oscillatory"
    if R > 1.0:
        return "divergent"
    return "monotonic"


def gci_triplet(phi_fine: float, phi_medium: float, phi_coarse: float,
                h_fine: float, h_medium: float, h_coarse: float,
                label: str = "",
                Fs: float = 1.25) -> GCIResult:
    phi1, phi2, phi3 = phi_fine, phi_medium, phi_coarse
    h1, h2, h3 = h_fine, h_medium, h_coarse
    r = h2 / h1
    eps21 = phi2 - phi1
    eps32 = phi3 - phi2

    if abs(eps21) < 1e-14 or abs(eps32) < 1e-14:
        return GCIResult(
            label, phi1, phi2, phi3, h1, h2, h3, r,
            eps21, eps32, float("nan"),
            "degenerate",
            float("nan"), float("nan"), float("nan"), float("nan"),
            float("nan"),
            "vanishing residual; observed order undefined",
        )

    R = eps21 / eps32
    regime = _classify(R)

    if regime != "monotonic":
        if regime == "divergent":
            msg = (f"R = {R:.3f} > 1: residual grows with refinement; "
                   "this triplet is pre-asymptotic")
        elif regime == "oscillatory":
            msg = (f"R = {R:.3f} < 0: sign change between successive "
                   "differences; not monotonic")
        else:
            msg = "degenerate"
        return GCIResult(
            label, phi1, phi2, phi3, h1, h2, h3, r,
            eps21, eps32, R, regime,
            float("nan"), float("nan"), float("nan"), float("nan"),
            float("nan"),
            msg,
        )

    p_obs = log(abs(eps32 / eps21)) / log(r)
    denom = r ** p_obs - 1.0
    phi_ext = (r ** p_obs * phi1 - phi2) / denom

    gci_fine = Fs * abs(eps21 / phi1) / denom
    gci_med  = Fs * abs(eps32 / phi2) / denom
    AR = gci_med / (r ** p_obs * gci_fine)

    if abs(AR - 1.0) < 0.05:
        note = "monotonic and well within asymptotic range (AR -> 1)"
    elif abs(AR - 1.0) < 0.15:
        note = "monotonic, approximately asymptotic"
    else:
        note = f"monotonic but AR = {AR:.3f} deviates from unity"

    return GCIResult(
        label, phi1, phi2, phi3, h1, h2, h3, r,
        eps21, eps32, R, "monotonic",
        p_obs, phi_ext,
        gci_fine * 100.0, gci_med * 100.0,
        AR, note,
    )
