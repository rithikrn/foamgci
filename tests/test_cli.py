"""Smoke test for the ``python -m foamgci`` CLI."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np


def _make_minmax(p: Path, mean: float, sigma: float = 0.01,
                 n: int = 500, t_end: float = 10.0, seed: int = 0) -> None:
    rng = np.random.default_rng(seed)
    t = np.linspace(0.0, t_end, n)
    s = mean + sigma * rng.standard_normal(n)
    lines = ["# Time  field  min  location(min)  proc  max  location(max)  proc"]
    for ti, si in zip(t, s):
        lines.append(f"{ti:.6f}  p  {si-1:.6f}  (0 0 0)  0  {si:.6f}  (0.6 0 0)  0")
    p.write_text("\n".join(lines) + "\n")


def test_cli_report_runs(tmp_path: Path) -> None:
    files = {}
    means = {"coarse": 11.99, "medium": 12.02, "fine": 12.05, "extra-fine": 12.058}
    hs = {"coarse": 0.025, "medium": 0.0125, "fine": 0.00625, "extra-fine": 0.003125}
    for i, (label, m) in enumerate(means.items()):
        f = tmp_path / f"{label}.dat"
        _make_minmax(f, mean=m, seed=200 + i)  # deterministic; not hash(label)
        files[label] = f
    out_text = tmp_path / "rep.txt"
    out_tex = tmp_path / "tab.tex"

    args = [sys.executable, "-m", "foamgci", "report",
            "--field", "p", "--quantity", "max",
            "--window", "3", "10",
            "--reference", "rayleigh-pitot", "--mach", "3", "--gamma", "1.4",
            "--text", str(out_text),
            "--latex", str(out_tex)]
    for label, f in files.items():
        args.extend(["--case", f"{label}:{f}:{hs[label]}:1000"])

    result = subprocess.run(args, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert "Roache GCI" in result.stdout
    assert out_text.exists()
    assert out_tex.exists()
    tex = out_tex.read_text()
    assert r"\bottomrule" in tex
