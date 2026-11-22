#!/usr/bin/env python
"""Generate a small set of synthetic McStas-format output files for testing
the plotter:  python make_example_data.py [output_folder]"""

import sys
from pathlib import Path

import numpy as np

SIM_HEADER = """# Format: McCode with text headers
# URL: http://www.mccode.org
# Creator: McStas 3.4 - synthetic example
# Instrument: Example_SANS.instr
# Ncount: 10000000
# Trace: no
# Gravitation: no
# Seed: 1234
# Directory: {directory}
# Param: lambda_min={lmin}
# Param: lambda_max=13
"""


def write_1d(path, comp, title, xvar, xlabel, xlim, x, y, var_name="I"):
    rng = np.random.default_rng(abs(hash(str(path))) % 2**32)
    n = np.maximum(np.round(y / y.max() * 5000), 0)
    err = np.where(n > 0, y / np.sqrt(np.maximum(n, 1)), 0.0)
    y_noisy = np.maximum(y + rng.normal(0, 1, y.size) * err, 0)
    header = SIM_HEADER.format(directory=path.parent, lmin=xlim[2]) + f"""# Date: Sun Nov 22 12:00:00 2026 (1795262400)
# type: array_1d({x.size})
# Source: Example_SANS (Example_SANS.instr)
# component: {comp}
# position: 0 0 {xlim[3]}
# title: {title}
# Ncount: 10000000
# filename: {path.name}
# statistics: X0=0; dX=0;
# signal: Min={y_noisy.min():.6g}; Max={y_noisy.max():.6g}; Mean={y_noisy.mean():.6g};
# values: {y_noisy.sum():.6g} {np.sqrt((err**2).sum()):.6g} {int(n.sum())}
# xvar: {xvar}
# yvar: (I,I_err)
# xlabel: {xlabel}
# ylabel: Intensity
# xlimits: {xlim[0]} {xlim[1]}
# variables: {xvar} {var_name} {var_name}_err N
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(header)
        for row in zip(x, y_noisy, err, n):
            fh.write(f"{row[0]:.6g} {row[1]:.6g} {row[2]:.6g} {int(row[3])}\n")


def maxwell(lam, T=60.0, lmin=0.9):
    a = 949.0 / T  # Maxwellian in wavelength, arbitrary scale
    s = (a**2 / lam**5) * np.exp(-a / lam**2)
    return 1e3 * s * (lam > lmin)


def make(out: Path):
    L = 25.0  # m, moderator-detector distance
    for i, lmin in enumerate([0.9, 1.75, 3.0, 4.5]):
        d = out / "scan" / str(i)
        lam = np.linspace(0.5, 15, 100)
        write_1d(d / "lambda_mon.L", "lam_mon", "Wavelength monitor", "L",
                 "Wavelength [AA]", (0.5, 15, lmin, 24.9), lam,
                 maxwell(lam, lmin=lmin) + 1e-3)
        t = np.linspace(0, 100000, 200)   # microseconds
        lam_t = t / (252.78 * L)
        write_1d(d / "tof_mon.t", "tof_mon", "Time-of-flight monitor", "t",
                 "Time-of-flight [\\gms]", (0, 100000, lmin, 25.0), t,
                 np.nan_to_num(maxwell(np.maximum(lam_t, 1e-3), lmin=lmin)))

    # Monitor_nD-style TOF output (time in seconds) in the top-level folder
    t_s = np.linspace(0, 0.1, 150)
    lam_t = t_s * 1e6 / (252.78 * 12.0)
    write_1d(out / "Mon_nD_tof_1795262400.t", "MonND", "Monitor_nD TOF",
             "t", "TOF [s]", (0, 0.1, 0.9, 12.0), t_s,
             maxwell(np.maximum(lam_t, 1e-3)))
    # An energy monitor (shows up as 'Other 1D')
    E = np.linspace(0, 50, 80)
    write_1d(out / "energy.E", "E_mon", "Energy monitor", "E",
             "Energy [meV]", (0, 50, 0.9, 10.0), E,
             E * np.exp(-E / 5.0) + 1e-3)

    # A 2D PSD file - must be ignored by the plotter
    psd = out / "psd.dat"
    psd.write_text(SIM_HEADER.format(directory=out, lmin=0.9)
                   + "# type: array_2d(4, 4)\n# component: psd\n"
                   + "# xvar: X\n# yvar: Y\n# Data [psd/psd.dat] I:\n"
                   + "\n".join(" ".join("1" for _ in range(4))
                               for _ in range(4)) + "\n")
    # McStas also writes mccode.sim - must be ignored
    (out / "mccode.sim").write_text("begin instrument: Example_SANS\nend\n")
    print(f"Example data written to {out.resolve()}")


if __name__ == "__main__":
    make(Path(sys.argv[1] if len(sys.argv) > 1 else "example_data"))
