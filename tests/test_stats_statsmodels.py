"""Optional cross-validation of foamgci.stats.kpss_test against statsmodels.

Skipped when statsmodels is not installed (it is NOT a dependency of
foamgci; this test exists so CI can verify the 'matches statsmodels'
claim in the kpss_test docstring instead of asserting it on faith).
"""
from __future__ import annotations

import numpy as np
import pytest

sm = pytest.importorskip("statsmodels.tsa.stattools")

from foamgci.stats import kpss_test  # noqa: E402


@pytest.mark.parametrize("regression", ["c", "ct"])
def test_kpss_matches_statsmodels(regression: str) -> None:
    import warnings
    rng = np.random.default_rng(7)
    series = {
        "white": rng.standard_normal(800),
        "trend": np.linspace(0, 3, 800) + 0.3 * rng.standard_normal(800),
        "rw": np.cumsum(rng.standard_normal(800)),
    }
    for x in series.values():
        mine = kpss_test(x, regression=regression)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            stat, _, _, _ = sm.kpss(x, regression=regression,
                                    nlags=mine["lag"])
        assert abs(mine["statistic"] - stat) < 1e-12
