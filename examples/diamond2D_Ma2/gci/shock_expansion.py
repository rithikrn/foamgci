"""Closed-form shock-expansion reference for a symmetric diamond airfoil.

Pure-numpy, bisection only (no scipy), matching the style of the wedge case's
``oblique_shock.py``.  For a diamond of half-angle ``eps`` at zero incidence in
a Mach ``M_inf`` stream:

* the front facet turns the flow by ``eps`` through an attached oblique shock,
  giving the surface pressure ``p_front/p_inf`` from the Rankine-Hugoniot jump;
* the rear facet turns the flow back by ``2*eps`` through a Prandtl-Meyer
  expansion, giving ``p_rear/p_inf`` isentropically;
* the wave-drag coefficient follows from the facet pressures,

      Cd = (p_front/p_inf - p_rear/p_inf) * (t/c) / (0.5 * gamma * M_inf^2),

  with thickness ratio ``t/c = tan(eps)``.

These are the references the grid-converged QoIs are checked against.
"""

from __future__ import annotations
import math
from dataclasses import dataclass


def _bisect(f, lo, hi, tol=1e-12, itmax=200):
    flo, fhi = f(lo), f(hi)
    if flo * fhi > 0:
        raise ValueError("root not bracketed")
    for _ in range(itmax):
        mid = 0.5 * (lo + hi)
        fm = f(mid)
        if abs(fm) < tol or 0.5 * (hi - lo) < tol:
            return mid
        if flo * fm < 0:
            hi, fhi = mid, fm
        else:
            lo, flo = mid, fm
    return 0.5 * (lo + hi)


def theta_beta_M(M, theta, gamma=1.4):
    """Weak-shock wave angle beta (rad) for deflection ``theta`` (rad).

    The theta-beta-M residual is positive at the Mach angle and at 90 deg and
    dips negative only between the weak and strong roots, so a fixed bracket
    of the whole range is not sign-changing.  Scan upward from the Mach angle
    in small steps for the first sign change (the weak root) and bisect there.
    """
    def resid(beta):
        return (
            math.tan(theta)
            - 2.0 / math.tan(beta)
            * (M ** 2 * math.sin(beta) ** 2 - 1.0)
            / (M ** 2 * (gamma + math.cos(2.0 * beta)) + 2.0)
        )
    mu = math.asin(1.0 / M)                 # Mach angle (lower limit)
    lo = mu + 1e-6
    flo = resid(lo)
    step = math.radians(0.2)
    b = lo + step
    hi = math.radians(89.5)
    while b < hi:
        fb = resid(b)
        if flo * fb <= 0.0:
            return _bisect(resid, lo, b)
        lo, flo, b = b, fb, b + step
    raise ValueError("no weak-shock root: deflection exceeds detachment angle")


def oblique_shock(M1, theta, gamma=1.4):
    """Oblique-shock jumps for deflection ``theta`` (rad). Returns a dict."""
    beta = theta_beta_M(M1, theta, gamma)
    M1n = M1 * math.sin(beta)
    p2_p1 = 1.0 + 2.0 * gamma / (gamma + 1.0) * (M1n ** 2 - 1.0)
    rho2_rho1 = (gamma + 1.0) * M1n ** 2 / ((gamma - 1.0) * M1n ** 2 + 2.0)
    T2_T1 = p2_p1 / rho2_rho1
    M2n = math.sqrt(
        (1.0 + 0.5 * (gamma - 1.0) * M1n ** 2)
        / (gamma * M1n ** 2 - 0.5 * (gamma - 1.0))
    )
    M2 = M2n / math.sin(beta - theta)
    return dict(beta=beta, M1n=M1n, p2_p1=p2_p1, rho2_rho1=rho2_rho1,
                T2_T1=T2_T1, M2=M2)


def prandtl_meyer(M, gamma=1.4):
    """Prandtl-Meyer angle nu(M) in radians."""
    gp, gm = gamma + 1.0, gamma - 1.0
    return (math.sqrt(gp / gm)
            * math.atan(math.sqrt(gm / gp * (M ** 2 - 1.0)))
            - math.atan(math.sqrt(M ** 2 - 1.0)))


def inverse_prandtl_meyer(nu_target, gamma=1.4):
    """Mach number with Prandtl-Meyer angle ``nu_target`` (rad)."""
    return _bisect(lambda M: prandtl_meyer(M, gamma) - nu_target, 1.0 + 1e-9, 60.0)


def _p0_over_p(M, gamma=1.4):
    return (1.0 + 0.5 * (gamma - 1.0) * M ** 2) ** (gamma / (gamma - 1.0))


@dataclass
class DiamondReference:
    M_inf: float
    eps_deg: float
    gamma: float
    beta_deg: float
    p_front: float          # /p_inf
    p_rear: float           # /p_inf
    M2: float               # behind the front-facet shock
    M3: float               # behind the rear-facet expansion
    t_over_c: float
    Cd_shock_expansion: float
    Cd_linear: float        # Ackeret 4 eps^2 / sqrt(M^2 - 1)


def diamond_reference(M_inf=2.0, eps_deg=10.0, gamma=1.4) -> DiamondReference:
    eps = math.radians(eps_deg)
    tc = math.tan(eps)
    # front facet: oblique shock, deflection eps
    sh = oblique_shock(M_inf, eps, gamma)
    p_front = sh["p2_p1"]
    M2 = sh["M2"]
    # rear facet: Prandtl-Meyer expansion through 2*eps from M2 to M3
    nu3 = prandtl_meyer(M2, gamma) + 2.0 * eps
    M3 = inverse_prandtl_meyer(nu3, gamma)
    # isentropic pressure from state 2 to state 3 (stagnation pressure conserved)
    p_rear = p_front * _p0_over_p(M2, gamma) / _p0_over_p(M3, gamma)
    Cd_se = (p_front - p_rear) * tc / (0.5 * gamma * M_inf ** 2)
    Cd_lin = 4.0 * eps ** 2 / math.sqrt(M_inf ** 2 - 1.0)
    return DiamondReference(
        M_inf=M_inf, eps_deg=eps_deg, gamma=gamma,
        beta_deg=math.degrees(sh["beta"]),
        p_front=p_front, p_rear=p_rear, M2=M2, M3=M3, t_over_c=tc,
        Cd_shock_expansion=Cd_se, Cd_linear=Cd_lin,
    )


if __name__ == "__main__":
    r = diamond_reference()
    print(f"M_inf={r.M_inf}, eps={r.eps_deg} deg, gamma={r.gamma}")
    print(f"beta            = {r.beta_deg:.4f} deg")
    print(f"p_front/p_inf   = {r.p_front:.4f}")
    print(f"p_rear/p_inf    = {r.p_rear:.4f}")
    print(f"M2, M3          = {r.M2:.4f}, {r.M3:.4f}")
    print(f"t/c             = {r.t_over_c:.4f}")
    print(f"Cd (shock-exp)  = {r.Cd_shock_expansion:.5f}")
    print(f"Cd (linear)     = {r.Cd_linear:.5f}")
