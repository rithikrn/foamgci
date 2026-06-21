"""Exact oblique-shock relations — the analytical reference for this case.

For an inviscid, calorically-perfect gas the state behind an attached
oblique shock is fixed exactly by the upstream Mach number ``M1``, the
flow-deflection (wedge) angle ``theta``, and ``gamma``. These relations are
the gold-standard verification reference: they are exact, closed-form, and
independent of mesh, solver version, and time step (NACA Report 1135, 1953;
Anderson, *Modern Compressible Flow*, Ch. 4).

This case uses the *ramp-surface (wall) static-pressure ratio* ``p2/p1`` as
the GCI quantity of interest. For inviscid flow the ramp-surface pressure
equals the post-shock static pressure ``p2`` everywhere on the straight
ramp, so the area-averaged wall pressure written by the ``surfaceFieldValue``
function object should converge to ``p2/p1`` under mesh refinement.

Dependency-light on purpose: NumPy only, weak-shock root by bracketed
bisection (no SciPy), mirroring the foamgci minimum-dependency philosophy.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def _theta_of_beta(beta: float, M1: float, gamma: float) -> float:
    """theta-beta-M relation: flow deflection produced by shock angle beta."""
    num = M1 ** 2 * np.sin(beta) ** 2 - 1.0
    den = M1 ** 2 * (gamma + np.cos(2.0 * beta)) + 2.0
    return float(np.arctan2(2.0 * (num / np.tan(beta)), den))


def shock_angle(M1: float, theta_deg: float, gamma: float = 1.4) -> float:
    """Weak-shock solution beta (radians) of the theta-beta-M relation.

    Bracketed bisection between the Mach angle and the detachment angle.
    Raises ValueError if the deflection exceeds the maximum for an attached
    shock (i.e. the shock would detach).
    """
    theta = np.radians(theta_deg)
    mach_angle = np.arcsin(1.0 / M1)
    betas = np.linspace(mach_angle + 1e-6, np.radians(89.999), 20000)
    diffs = np.array([_theta_of_beta(b, M1, gamma) for b in betas]) - theta
    sign_changes = np.where(np.sign(diffs[:-1]) != np.sign(diffs[1:]))[0]
    if sign_changes.size == 0:
        raise ValueError(
            f"No attached oblique shock for M1={M1}, theta={theta_deg} deg, "
            f"gamma={gamma} (deflection exceeds detachment)."
        )
    # The first sign change (smallest beta) is the weak-shock branch.
    lo, hi = betas[sign_changes[0]], betas[sign_changes[0] + 1]
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if (_theta_of_beta(lo, M1, gamma) - theta) * (
            _theta_of_beta(mid, M1, gamma) - theta
        ) <= 0.0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi)


@dataclass(frozen=True)
class ObliqueShockState:
    M1: float
    theta_deg: float
    gamma: float
    beta_deg: float       # shock angle
    Mn1: float            # upstream normal Mach
    p2_p1: float          # static-pressure ratio  <-- WALL-PRESSURE QoI ref
    rho2_rho1: float
    T2_T1: float
    M2: float


def oblique_shock(M1: float, theta_deg: float, gamma: float = 1.4) -> ObliqueShockState:
    """Full downstream state for an attached weak oblique shock."""
    beta = shock_angle(M1, theta_deg, gamma)
    th = np.radians(theta_deg)
    Mn1 = M1 * np.sin(beta)
    p2_p1 = 1.0 + 2.0 * gamma / (gamma + 1.0) * (Mn1 ** 2 - 1.0)
    rho2_rho1 = (gamma + 1.0) * Mn1 ** 2 / ((gamma - 1.0) * Mn1 ** 2 + 2.0)
    T2_T1 = p2_p1 / rho2_rho1
    Mn2 = np.sqrt(
        (1.0 + 0.5 * (gamma - 1.0) * Mn1 ** 2)
        / (gamma * Mn1 ** 2 - 0.5 * (gamma - 1.0))
    )
    M2 = Mn2 / np.sin(beta - th)
    return ObliqueShockState(
        M1=float(M1),
        theta_deg=float(theta_deg),
        gamma=float(gamma),
        beta_deg=float(np.degrees(beta)),
        Mn1=float(Mn1),
        p2_p1=float(p2_p1),
        rho2_rho1=float(rho2_rho1),
        T2_T1=float(T2_T1),
        M2=float(M2),
    )


if __name__ == "__main__":
    s = oblique_shock(5.0, 15.0, 1.4)
    print(f"M1={s.M1}  theta={s.theta_deg} deg  gamma={s.gamma}")
    print(f"  shock angle beta = {s.beta_deg:.5f} deg")
    print(f"  p2/p1            = {s.p2_p1:.6f}   <- wall-pressure QoI reference")
    print(f"  rho2/rho1        = {s.rho2_rho1:.6f}")
    print(f"  T2/T1            = {s.T2_T1:.6f}")
    print(f"  M2               = {s.M2:.6f}")
