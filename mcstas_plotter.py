#!/usr/bin/env python
"""
McStas 1D Monitor Plotter
=========================

A Qt (PySide6) + matplotlib application for browsing, comparing,
manipulating and plotting the text output files written by McStas / McCode
1D monitors, in particular time-of-flight (TOF_monitor, Monitor_nD "t") and
wavelength (L_monitor, Monitor_nD "lambda") monitors.

Features
--------
* Load one or more data folders (A, B, C ...) so that different simulations
  can be compared; optionally include sub-folders (e.g. scan points).
* All McStas 1D monitor files found are listed with folder tag, type
  (TOF / Wavelength / Other 1D), component, bins and axis label.
* Plot selections in a new figure (tab) or overplot on the current figure,
  with or without error bars (per curve).
* Data operations: A + B, A - B, A x B, A / B between two datasets whose
  x axes agree to better than 0.1 % and have the same number of points, or
  with a constant (offset / scale). Uncertainties are propagated assuming
  uncorrelated errors. Results appear in the "Derived data" list and can be
  plotted, used in further operations and saved.
* Save any dataset in McStas 1D text format (re-readable by McStas tools and
  by this program) or as simple 3-column text (x, y, y_err).
* Manual / automatic axis limits, linear or log axes, line / marker / step
  styles, normalisation, legend, grid, matplotlib toolbar (zoom, save ...).

Usage
-----
    python mcstas_plotter.py [data_folder ...]
"""

from __future__ import annotations

import os
import re
import string
import sys
import time
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np

# Import the Qt binding *before* matplotlib so that matplotlib's QtAgg
# backend picks PySide6.
from PySide6 import QtCore, QtGui, QtWidgets
from PySide6.QtCore import Qt

import matplotlib

matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import (  # noqa: E402
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure  # noqa: E402

APP_NAME = "McStas 1D Monitor Plotter"
__version__ = "1.1.0"

KIND_TOF = "TOF"
KIND_WL = "Wavelength"
KIND_OTHER = "Other 1D"

# Relative tolerance for x-axis agreement in data operations (0.1 %)
X_TOLERANCE = 1e-3

# File types that are never McStas text monitor files - skipped for speed.
SKIP_SUFFIXES = {
    ".sim", ".instr", ".comp", ".c", ".h", ".o", ".out", ".exe", ".so",
    ".dll", ".py", ".pyc", ".png", ".jpg", ".jpeg", ".gif", ".pdf", ".svg",
    ".eps", ".ps", ".h5", ".hdf", ".hdf5", ".nxs", ".nx", ".zip", ".gz",
    ".tar", ".json", ".xml", ".html", ".log",
}
MAX_FILE_SIZE = 200 * 1024 * 1024  # bytes
MAX_FILES_SCANNED = 20000

HEADER_RE = re.compile(r"^#\s*([^:]+?)\s*:\s?(.*)$")
NBINS_RE = re.compile(r"array_1d\s*\(\s*(\d+)\s*\)", re.IGNORECASE)

STYLES = ["Line", "Markers", "Line + markers", "Steps"]

# operation code -> (symbol, description)
OPERATIONS = {
    "+": ("+", "add"),
    "-": ("−", "subtract"),
    "*": ("×", "multiply / scale"),
    "/": ("÷", "divide"),
}


# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------
def natural_key(text: str):
    """Sort key that orders 'scan/2' before 'scan/10'."""
    return [int(tok) if tok.isdigit() else tok.lower()
            for tok in re.split(r"(\d+)", text)]


def pretty_label(label: str) -> str:
    """Convert McStas label conventions to readable unicode text."""
    if not label:
        return ""
    out = label.strip().strip('"')
    replacements = [
        (r"\gms", "µs"), (r"\gm", "µ"), (r"\gl", "λ"), (r"\gD", "Δ"),
        (r"\gq", "θ"), (r"\gf", "φ"), (r"\gW", "Ω"),
        ("[AA]", "[Å]"), ("[Angs]", "[Å]"), ("[Angstrom]", "[Å]"),
        ("[AA^-1]", "[Å⁻¹]"), ("[Angs-1]", "[Å⁻¹]"),
    ]
    for old, new in replacements:
        out = out.replace(old, new)
    # A lone '$' would trigger matplotlib mathtext parsing errors.
    if out.count("$") % 2 == 1:
        out = out.replace("$", r"\$")
    return out


def classify(meta: dict) -> str:
    """Decide whether a 1D monitor is TOF, wavelength or something else."""
    xvar = meta.get("xvar", "").strip().lower()
    xlabel = meta.get("xlabel", "").lower()
    title = meta.get("title", "").lower()
    if (xvar in ("t", "tof", "time")
            or "time" in xlabel or "tof" in xlabel
            or "time-of-flight" in title or "tof" in title):
        return KIND_TOF
    if (xvar in ("l", "lambda", "wl", "wavelength")
            or "wavelength" in xlabel or "lambda" in xlabel
            or "[aa]" in xlabel or "[angs" in xlabel or "\\gl" in xlabel
            or "wavelength" in title):
        return KIND_WL
    return KIND_OTHER


_ASCII_MAP = {"÷": "/", "×": "*", "−": "-", "Å": "AA", "µ": "u",
              "λ": "lambda", "⁻¹": "^-1", "Δ": "D", "θ": "theta",
              "φ": "phi", "Ω": "Omega"}


def to_ascii(text: str) -> str:
    """Plain-ASCII version of a header string (for portable files)."""
    for k, v in _ASCII_MAP.items():
        text = text.replace(k, v)
    return text.encode("ascii", "replace").decode("ascii")


def safe_filename(name: str) -> str:
    out = re.sub(r"[^\w.\-+]+", "_", name).strip("._")
    return out or "dataset"


# --------------------------------------------------------------------------
# Data model
# --------------------------------------------------------------------------
@dataclass
class MonitorFile:
    """A McStas 1D monitor file found on disk (header only)."""
    path: Path
    rel: str
    tag: str
    meta: dict
    header_text: str
    kind: str
    nbins: int | None

    @property
    def component(self) -> str:
        return self.meta.get("component", "")

    @property
    def title(self) -> str:
        return self.meta.get("title", "")

    @property
    def xlabel(self) -> str:
        return pretty_label(self.meta.get("xlabel", "")) or "x"

    @property
    def ylabel(self) -> str:
        return pretty_label(self.meta.get("ylabel", "")) or "Intensity"

    @property
    def params(self) -> str:
        return self.meta.get("param", "")


@dataclass
class Dataset:
    """A loaded 1D dataset: from a file or derived by an operation."""
    name: str
    x: np.ndarray
    y: np.ndarray
    e: np.ndarray
    n: np.ndarray | None
    xlabel: str
    ylabel: str
    kind: str
    source: str
    meta: dict = field(default_factory=dict)
    derived: bool = False
    description: str = ""

    @property
    def npts(self) -> int:
        return int(self.x.size)


def parse_monitor_header(path: Path, root: Path,
                         tag: str = "A") -> MonitorFile | None:
    """Read the '#' header of a file; return a MonitorFile if it is a
    McStas/McCode 1D monitor, otherwise None."""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            first = fh.readline()
            if not first.startswith("#"):
                return None
            lines = [first]
            for _ in range(1000):
                line = fh.readline()
                if not line or not line.startswith("#"):
                    break
                lines.append(line)
    except (OSError, UnicodeError):
        return None

    meta: dict[str, str] = {}
    params: list[str] = []
    for line in lines:
        m = HEADER_RE.match(line.rstrip("\r\n"))
        if not m:
            continue
        key = m.group(1).strip().lower()
        value = m.group(2).strip()
        if key == "param":
            params.append(value)
        else:
            meta[key] = value
    if params:
        meta["param"] = "; ".join(params)

    mtype = meta.get("type", "").strip().lower()
    if not mtype.startswith("array_1d"):
        return None

    nb = NBINS_RE.search(mtype)
    try:
        rel = str(path.relative_to(root))
    except ValueError:
        rel = path.name
    return MonitorFile(
        path=path,
        rel=rel.replace(os.sep, "/"),
        tag=tag,
        meta=meta,
        header_text="".join(lines),
        kind=classify(meta),
        nbins=int(nb.group(1)) if nb else None,
    )


def scan_folder(folder: Path, recursive: bool, tag: str = "A",
                progress_cb=None) -> tuple[list[MonitorFile], bool]:
    """Find McStas 1D monitor files. Returns (files, truncated)."""
    results: list[MonitorFile] = []
    n_checked = 0
    if recursive:
        walker = os.walk(folder)
    else:
        try:
            names = [e.name for e in os.scandir(folder) if e.is_file()]
        except OSError:
            names = []
        walker = [(str(folder), [], names)]

    for root, dirs, fnames in walker:
        dirs[:] = sorted((d for d in dirs if not d.startswith(".")),
                         key=natural_key)
        for fn in sorted(fnames, key=natural_key):
            if fn.startswith(".") or Path(fn).suffix.lower() in SKIP_SUFFIXES:
                continue
            p = Path(root) / fn
            try:
                if p.stat().st_size > MAX_FILE_SIZE:
                    continue
            except OSError:
                continue
            n_checked += 1
            if n_checked > MAX_FILES_SCANNED:
                return results, True
            if progress_cb and n_checked % 200 == 0:
                progress_cb()
            mf = parse_monitor_header(p, folder, tag)
            if mf is not None:
                results.append(mf)
    return results, False


def load_monitor_data(path: Path):
    """Return x, I, I_err, N arrays from a McStas 1D monitor file."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            arr = np.loadtxt(path, comments="#", ndmin=2)
        except ValueError:
            arr = np.atleast_2d(np.genfromtxt(path, comments="#"))
    if arr.size == 0:
        raise ValueError("file contains no data rows")
    ncol = arr.shape[1]
    if ncol < 2:
        raise ValueError(f"expected at least 2 columns, found {ncol}")
    x = arr[:, 0]
    y = arr[:, 1]
    e = np.abs(arr[:, 2]) if ncol > 2 else np.zeros_like(y)
    n = arr[:, 3] if ncol > 3 else None
    return x, y, e, n


def dataset_from_file(mf: MonitorFile, name: str) -> Dataset:
    x, y, e, n = load_monitor_data(mf.path)
    return Dataset(name=name, x=x, y=y, e=e, n=n, xlabel=mf.xlabel,
                   ylabel=mf.ylabel, kind=mf.kind, source=str(mf.path),
                   meta=dict(mf.meta), description=mf.header_text)


# --------------------------------------------------------------------------
# Data operations
# --------------------------------------------------------------------------
def check_x_compat(xa: np.ndarray, xb: np.ndarray,
                   tol: float = X_TOLERANCE):
    """Check two x axes have the same length and agree point by point to
    within `tol` (relative). Returns (ok, message, max_relative_deviation).

    The deviation of each point is |xa - xb| / max(|xa|, |xb|); for points
    close to zero (|x| < tol * x-range) the full x-range is used as the
    reference instead, so that a bin centred on 0 does not fail."""
    if xa.size != xb.size:
        return (False, f"different number of points ({xa.size} vs "
                f"{xb.size})", None)
    if xa.size == 0:
        return False, "empty dataset", None
    span = max(float(np.ptp(xa)), float(np.ptp(xb)))
    ref = np.maximum(np.abs(xa), np.abs(xb))
    ref = np.where(ref < tol * span, span, ref)
    ref = np.where(ref == 0, 1.0, ref)
    dev = np.abs(xa - xb) / ref
    maxdev = float(np.nanmax(dev))
    if not np.isfinite(maxdev) or maxdev > tol:
        i = int(np.nanargmax(dev))
        return (False, f"x values differ by {100 * maxdev:.3g} % at point "
                f"{i} ({xa[i]:.6g} vs {xb[i]:.6g}); limit is "
                f"{100 * tol:g} %", maxdev)
    return (True, f"x axes match: {xa.size} points, max deviation "
            f"{100 * maxdev:.2g} % (limit {100 * tol:g} %)", maxdev)


def apply_operation(a: Dataset, op: str, b: Dataset | None = None,
                    const: float | None = None, name: str = ""):
    """Return (result Dataset, number of undefined points).

    op is one of '+', '-', '*', '/'. Either b (a Dataset) or const (a
    float) must be given. Errors are propagated assuming the two operands
    are uncorrelated; constants are treated as exact."""
    if op not in OPERATIONS:
        raise ValueError(f"unknown operation {op!r}")
    sym = OPERATIONS[op][0]
    ya = a.y.astype(float)
    ea = a.e.astype(float)
    n_bad = 0
    with np.errstate(all="ignore"):
        if b is not None:
            ok, msg, _ = check_x_compat(a.x, b.x)
            if not ok:
                raise ValueError(msg)
            yb = b.y.astype(float)
            eb = b.e.astype(float)
            if op == "+":
                y, e = ya + yb, np.hypot(ea, eb)
            elif op == "-":
                y, e = ya - yb, np.hypot(ea, eb)
            elif op == "*":
                y, e = ya * yb, np.hypot(yb * ea, ya * eb)
            else:
                bad = yb == 0
                n_bad = int(bad.sum())
                y = ya / yb
                e = np.hypot(ea / yb, ya * eb / yb ** 2)
                y[bad] = np.nan
                e[bad] = np.nan
            if op in "+-" and a.n is not None and b.n is not None:
                n = a.n + b.n
            else:
                n = None if a.n is None else a.n.copy()
            expr = f"({a.name}) {sym} ({b.name})"
            if op == "/":
                ylabel = "Ratio"
            elif op == "*":
                ylabel = "Product"
            else:
                ylabel = a.ylabel
            operands = [a, b]
        else:
            if const is None or not np.isfinite(const):
                raise ValueError("a finite constant is required")
            c = float(const)
            if op == "+":
                y, e = ya + c, ea.copy()
            elif op == "-":
                y, e = ya - c, ea.copy()
            elif op == "*":
                y, e = ya * c, ea * abs(c)
            else:
                if c == 0:
                    raise ValueError("division by zero")
                y, e = ya / c, ea / abs(c)
            n = None if a.n is None else a.n.copy()
            expr = f"({a.name}) {sym} {c:g}"
            ylabel = a.ylabel
            operands = [a]

    meta = dict(a.meta)
    meta["title"] = expr
    meta["ylabel"] = ylabel
    desc = [f"Derived dataset: {name or expr}",
            f"Operation: {expr}",
            f"Created: {datetime.now():%Y-%m-%d %H:%M:%S}", ""]
    for lbl, ds in zip("AB", operands):
        desc.append(f"Operand {lbl}: {ds.name}")
        desc.append(f"    source: {ds.source}")
    if const is not None and b is None:
        desc.append(f"Constant: {const:g}")
    desc += ["",
             f"Points: {a.npts}   x range: {np.nanmin(a.x):.6g} .. "
             f"{np.nanmax(a.x):.6g}",
             "Errors propagated assuming uncorrelated operands."]
    if n_bad:
        desc.append(f"{n_bad} point(s) undefined (division by zero) "
                    "set to NaN.")
    result = Dataset(name=name or expr, x=a.x.copy(), y=y, e=e, n=n,
                     xlabel=a.xlabel, ylabel=ylabel, kind=a.kind,
                     source=f"derived: {expr}", meta=meta, derived=True,
                     description="\n".join(desc))
    return result, n_bad


# --------------------------------------------------------------------------
# Writers
# --------------------------------------------------------------------------
def write_mcstas_1d(path: Path, ds: Dataset):
    """Write a dataset in McStas/McCode 1D text format."""
    x = np.asarray(ds.x, float)
    y = np.asarray(ds.y, float)
    e = np.asarray(ds.e, float)
    n = (np.asarray(ds.n, float) if ds.n is not None
         else np.zeros_like(y))
    m = dict(ds.meta)
    ok = np.isfinite(y)
    s_i = float(y[ok].sum()) if ok.any() else 0.0
    if s_i != 0:
        x0 = float((x[ok] * y[ok]).sum() / s_i)
        dx = float(np.sqrt(abs((x[ok] ** 2 * y[ok]).sum() / s_i - x0 ** 2)))
    else:
        x0 = dx = 0.0
    ymin = float(np.nanmin(y)) if ok.any() else 0.0
    ymax = float(np.nanmax(y)) if ok.any() else 0.0
    ymean = float(np.nanmean(y)) if ok.any() else 0.0
    err_tot = float(np.sqrt(np.nansum(e ** 2)))
    xvar = m.get("xvar") or {KIND_TOF: "t", KIND_WL: "L"}.get(ds.kind, "x")
    xlabel = m.get("xlabel") or ds.xlabel
    ylabel = m.get("ylabel") or ds.ylabel
    now = time.time()
    lines = [
        "Format: McCode with text headers",
        "URL: http://www.mccode.org",
        f"Creator: {APP_NAME} {__version__}"
        + (" (derived data)" if ds.derived else " (re-saved data)"),
        f"Instrument: {m.get('instrument', '')}",
        f"Ncount: {m.get('ncount', '0')}",
    ]
    if ds.derived:
        lines.append(f"Derived: {m.get('title', ds.name)}")
    else:
        lines.append(f"Original: {ds.source}")
    if m.get("param"):
        for p in m["param"].split("; "):
            lines.append(f"Param: {p}")
    lines += [
        f"Date: {datetime.fromtimestamp(now):%a %b %d %H:%M:%S %Y} "
        f"({int(now)})",
        f"type: array_1d({x.size})",
        f"Source: {m.get('source', ds.source)}",
        f"component: {m.get('component', 'derived') or 'derived'}",
        f"position: {m.get('position', '0 0 0')}",
        f"title: {ds.name}",
        f"Ncount: {m.get('ncount', '0')}",
        f"filename: {path.name}",
        f"statistics: X0={x0:.6g}; dX={dx:.6g};",
        f"signal: Min={ymin:.6g}; Max={ymax:.6g}; Mean={ymean:.6g};",
        f"values: {s_i:.6g} {err_tot:.6g} {int(np.nansum(n))}",
        f"xvar: {xvar}",
        "yvar: (I,I_err)",
        f"xlabel: {xlabel}",
        f"ylabel: {ylabel}",
        f"xlimits: {np.nanmin(x):.6g} {np.nanmax(x):.6g}",
        f"variables: {xvar} I I_err N",
    ]
    header = "\n".join("# " + to_ascii(ln) for ln in lines)
    data = np.column_stack([x, y, e, n])
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(header + "\n")
        np.savetxt(fh, data, fmt=["%.10g", "%.10g", "%.10g", "%.10g"])


def write_three_column(path: Path, ds: Dataset):
    """Write a dataset as plain text: x  y  y_err (with '#' comments)."""
    data = np.column_stack([ds.x, ds.y, ds.e])
    xlabel = ds.meta.get("xlabel") or ds.xlabel
    ylabel = ds.meta.get("ylabel") or ds.ylabel
    header = to_ascii(f"{ds.name}\n"
                      f"source: {ds.source}\n"
                      f"x: {xlabel}   y: {ylabel}\n"
                      f"columns: x y y_err")
    np.savetxt(path, data, fmt="%.10g", header=header, comments="# ",
               encoding="ascii")


# --------------------------------------------------------------------------
# Plot tab
# --------------------------------------------------------------------------
@dataclass
class Curve:
    label: str
    x: np.ndarray
    y: np.ndarray
    e: np.ndarray
    xlabel: str
    ylabel: str
    kind: str
    source: str
    color: str = "C0"
    errorbars: bool = True
    visible: bool = True

    @classmethod
    def from_dataset(cls, ds: Dataset, errorbars: bool) -> "Curve":
        return cls(label=ds.name, x=ds.x, y=ds.y, e=ds.e, xlabel=ds.xlabel,
                   ylabel=ds.ylabel, kind=ds.kind, source=ds.source,
                   errorbars=errorbars)


@dataclass
class AxisState:
    xscale: str = "linear"
    yscale: str = "linear"
    xlim: tuple | None = None
    ylim: tuple | None = None
    style: str = "Line"
    normalise: bool = False
    legend: bool = True
    grid: bool = True


class PlotTab(QtWidgets.QWidget):
    """One figure: a matplotlib canvas + toolbar and the curves on it."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.figure = Figure(figsize=(7, 5), layout="constrained")
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                  QtWidgets.QSizePolicy.Expanding)
        self.toolbar = NavigationToolbar(self.canvas, self)
        lay = QtWidgets.QVBoxLayout(self)
        lay.setContentsMargins(2, 2, 2, 2)
        lay.addWidget(self.toolbar)
        lay.addWidget(self.canvas)
        self.ax = self.figure.add_subplot(111)
        self.curves: list[Curve] = []
        self.state = AxisState()
        self._color_index = 0
        self.redraw()

    def _next_color(self) -> str:
        cycle = matplotlib.rcParams["axes.prop_cycle"].by_key().get(
            "color", ["C0"])
        c = cycle[self._color_index % len(cycle)]
        self._color_index += 1
        return c

    def add_curve(self, curve: Curve):
        curve.color = self._next_color()
        self.curves.append(curve)

    def kinds(self) -> set[str]:
        return {c.kind for c in self.curves}

    def _draw_curve(self, c: Curve):
        ax = self.ax
        st = self.state
        x, y, e = c.x, c.y.astype(float), c.e.astype(float)
        if st.normalise:
            finite = y[np.isfinite(y)]
            m = np.max(np.abs(finite)) if finite.size else 0.0
            if m > 0:
                y = y / m
                e = e / m
        e = np.nan_to_num(e, nan=0.0, posinf=0.0, neginf=0.0)

        yerr = e
        if c.errorbars and st.yscale == "log":
            # Keep the lower bar above zero so it can be drawn on a log axis
            with np.errstate(invalid="ignore"):
                lower = np.where(y - e > 0, e, 0.95 * y)
                lower = np.where(y > 0, lower, 0.0)
            yerr = np.vstack([np.nan_to_num(lower), e])

        kw = dict(color=c.color, label=c.label)
        if st.style == "Steps":
            (handle,) = ax.step(x, y, where="mid", lw=1.2, **kw)
            if c.errorbars:
                ax.errorbar(x, y, yerr=yerr, fmt="none", ecolor=c.color,
                            elinewidth=0.8, capsize=0, label="_nolegend_")
        else:
            fmt = {"Line": "-", "Markers": "o",
                   "Line + markers": "-o"}.get(st.style, "-")
            if c.errorbars:
                handle = ax.errorbar(x, y, yerr=yerr, fmt=fmt, ms=3, lw=1.2,
                                     elinewidth=0.8, capsize=2, **kw)
            else:
                (handle,) = ax.plot(x, y, fmt, ms=3, lw=1.2, **kw)
        return handle

    def redraw(self):
        ax = self.ax
        st = self.state
        ax.clear()
        shown = [c for c in self.curves if c.visible]
        handles = [self._draw_curve(c) for c in shown]

        xlabels = list(dict.fromkeys(c.xlabel for c in shown))
        ylabels = list(dict.fromkeys(c.ylabel for c in shown))
        ax.set_xlabel(" / ".join(xlabels) if xlabels else "x")
        if st.normalise:
            ax.set_ylabel("Intensity (normalised to max)")
        else:
            ax.set_ylabel(" / ".join(ylabels) if ylabels else "Intensity")

        ax.set_xscale("log", nonpositive="mask") if st.xscale == "log" \
            else ax.set_xscale("linear")
        ax.set_yscale("log", nonpositive="mask") if st.yscale == "log" \
            else ax.set_yscale("linear")

        if shown:
            ax.relim()
            ax.autoscale_view()
        if st.xlim is not None:
            ax.set_xlim(*st.xlim)
        if st.ylim is not None:
            ax.set_ylim(*st.ylim)

        if st.grid:
            ax.grid(True, which="major", alpha=0.35)
            if "log" in (st.xscale, st.yscale):
                ax.grid(True, which="minor", alpha=0.12)
        else:
            ax.grid(False)

        if st.legend and shown:
            ax.legend(handles, [c.label for c in shown], fontsize="small",
                      loc="best", framealpha=0.85)

        if not self.curves:
            ax.text(0.5, 0.5,
                    "Select data on the left and click\n"
                    "'Plot in new figure' or 'Overplot on current'",
                    transform=ax.transAxes, ha="center", va="center",
                    color="0.55")
        self.canvas.draw_idle()


# --------------------------------------------------------------------------
# Tree items
# --------------------------------------------------------------------------
class FileItem(QtWidgets.QTreeWidgetItem):
    COL_SET, COL_FILE, COL_KIND, COL_COMP, COL_BINS, COL_X, COL_TITLE = \
        range(7)

    def __init__(self, mf: MonitorFile):
        super().__init__([mf.tag, mf.rel, mf.kind, mf.component,
                          "" if mf.nbins is None else str(mf.nbins),
                          mf.xlabel, mf.title])
        self.mf = mf
        tip = f"{mf.path}"
        if mf.params:
            tip += f"\nParameters: {mf.params}"
        for col in range(7):
            self.setToolTip(col, tip)
        self.setTextAlignment(self.COL_BINS, Qt.AlignRight | Qt.AlignVCenter)
        self.setTextAlignment(self.COL_SET, Qt.AlignCenter)

    def __lt__(self, other):
        tree = self.treeWidget()
        col = tree.sortColumn() if tree else 0
        if col == self.COL_BINS:
            return (self.mf.nbins or 0) < (other.mf.nbins or 0)
        return natural_key(self.text(col)) < natural_key(other.text(col))


class DerivedItem(QtWidgets.QTreeWidgetItem):
    def __init__(self, ds: Dataset):
        super().__init__([ds.name, ds.meta.get("title", ""), str(ds.npts),
                          ds.xlabel])
        self.ds = ds
        self.setToolTip(0, ds.description)
        self.setToolTip(1, ds.description)


# --------------------------------------------------------------------------
# Operations dialog
# --------------------------------------------------------------------------
class OperationsDialog(QtWidgets.QDialog):
    """Build a derived dataset from A (op) B or A (op) constant."""

    AFTER = ["Plot in new figure", "Overplot on current figure",
             "Don't plot"]

    def __init__(self, win: "MainWindow", pre_a=None, pre_b=None):
        super().__init__(win)
        self.win = win
        self.setWindowTitle("Data operations")
        self.setMinimumWidth(640)
        self._cache: dict[int, Dataset] = {}
        self._name_auto = True

        form = QtWidgets.QFormLayout(self)
        self.a_combo = QtWidgets.QComboBox()
        self.a_combo.setSizeAdjustPolicy(
            QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        form.addRow("Dataset A:", self.a_combo)

        self.op_combo = QtWidgets.QComboBox()
        for code, (sym, desc) in OPERATIONS.items():
            self.op_combo.addItem(f"{sym}   ({desc})", code)
        form.addRow("Operation:", self.op_combo)

        b_widget = QtWidgets.QWidget()
        bl = QtWidgets.QGridLayout(b_widget)
        bl.setContentsMargins(0, 0, 0, 0)
        self.b_ds_radio = QtWidgets.QRadioButton("Dataset B")
        self.b_const_radio = QtWidgets.QRadioButton("Constant")
        self.b_ds_radio.setChecked(True)
        self.b_combo = QtWidgets.QComboBox()
        self.b_combo.setSizeAdjustPolicy(
            QtWidgets.QComboBox.AdjustToMinimumContentsLengthWithIcon)
        self.const_edit = QtWidgets.QLineEdit("1")
        self.const_edit.setToolTip(
            "Exact constant, e.g. 2, 0.5, 1e-3.\n"
            "Use × to scale, + / − to offset.")
        bl.addWidget(self.b_ds_radio, 0, 0)
        bl.addWidget(self.b_combo, 0, 1)
        bl.addWidget(self.b_const_radio, 1, 0)
        bl.addWidget(self.const_edit, 1, 1)
        bl.setColumnStretch(1, 1)
        form.addRow("Operand:", b_widget)

        self.name_edit = QtWidgets.QLineEdit()
        self.name_edit.textEdited.connect(self._on_name_edited)
        form.addRow("Result name:", self.name_edit)

        self.status = QtWidgets.QLabel()
        self.status.setWordWrap(True)
        self.status.setTextFormat(Qt.RichText)
        self.status.setMinimumHeight(48)
        form.addRow("Check:", self.status)

        self.after_combo = QtWidgets.QComboBox()
        self.after_combo.addItems(self.AFTER)
        self.after_combo.setCurrentIndex(int(
            win.settings.value("ops_after", 0, type=int)))
        form.addRow("After creating:", self.after_combo)

        note = QtWidgets.QLabel(
            f"Datasets must have the same number of points and x values "
            f"agreeing to within {100 * X_TOLERANCE:g} %. Errors are "
            "propagated assuming uncorrelated operands; constants are "
            "exact. The result takes the x values of A.")
        note.setWordWrap(True)
        note.setStyleSheet("color: gray; font-size: 11px;")
        form.addRow(note)

        bb = QtWidgets.QDialogButtonBox()
        self.create_btn = bb.addButton(
            "Create", QtWidgets.QDialogButtonBox.AcceptRole)
        close_btn = bb.addButton(QtWidgets.QDialogButtonBox.Close)
        self.create_btn.clicked.connect(self.create)
        close_btn.clicked.connect(self.reject)
        form.addRow(bb)

        self._load_candidates(pre_a, pre_b)
        if pre_b is None and pre_a is not None:
            self.b_const_radio.setChecked(True)

        for w in (self.a_combo, self.b_combo, self.op_combo):
            w.currentIndexChanged.connect(self.update_status)
        self.b_ds_radio.toggled.connect(self.update_status)
        self.const_edit.textChanged.connect(self.update_status)
        self.update_status()

    # -- candidates -----------------------------------------------------
    def _load_candidates(self, pre_a=None, pre_b=None):
        """pre_a / pre_b: MonitorFile or Dataset objects to preselect."""
        self.cands = self.win.operand_candidates()
        for combo, pre in ((self.a_combo, pre_a), (self.b_combo, pre_b)):
            old = combo.currentIndex()
            combo.blockSignals(True)
            combo.clear()
            for i, (label, _obj) in enumerate(self.cands):
                combo.addItem(label, i)
            idx = 0
            if pre is not None:
                for i, (_l, obj) in enumerate(self.cands):
                    if obj is pre:
                        idx = i
                        break
            elif 0 <= old < combo.count():
                idx = old
            combo.setCurrentIndex(idx)
            combo.blockSignals(False)
        if pre_a is not None and pre_b is None and len(self.cands) > 1:
            # sensible default for B: the next candidate
            self.b_combo.setCurrentIndex(
                (self.a_combo.currentIndex() + 1) % len(self.cands))

    def _dataset(self, combo) -> Dataset:
        i = combo.currentData()
        if i is None:
            raise ValueError("no data available")
        if i not in self._cache:
            self._cache[i] = self.win.dataset_for(self.cands[i][1])
        return self._cache[i]

    # -- UI logic ---------------------------------------------------------
    def _on_name_edited(self, _text):
        self._name_auto = not self.name_edit.text().strip()

    def _set_status(self, html, ok):
        self.status.setText(html)
        self.create_btn.setEnabled(ok)

    def _inputs(self):
        op = self.op_combo.currentData()
        a = self._dataset(self.a_combo)
        if self.b_ds_radio.isChecked():
            return a, op, self._dataset(self.b_combo), None
        txt = self.const_edit.text().strip()
        try:
            c = float(txt)
        except ValueError:
            raise ValueError(f"'{txt}' is not a number") from None
        if not np.isfinite(c):
            raise ValueError("constant must be finite")
        if op == "/" and c == 0:
            raise ValueError("division by zero")
        return a, op, None, c

    def update_status(self, *_):
        use_b = self.b_ds_radio.isChecked()
        self.b_combo.setEnabled(use_b)
        self.const_edit.setEnabled(not use_b)
        if not self.cands:
            self._set_status("<span style='color:#b00'>No datasets "
                             "available – load a folder first.</span>", False)
            return
        try:
            a, op, b, c = self._inputs()
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"<span style='color:#b00'>✗ {exc}</span>",
                             False)
            return
        sym = OPERATIONS[op][0]
        if b is not None:
            expr = f"({a.name}) {sym} ({b.name})"
            ok, msg, _ = check_x_compat(a.x, b.x)
            if not ok:
                self._set_status(f"<span style='color:#b00'>✗ {msg}</span>",
                                 False)
                return
            html = f"<span style='color:#080'>✓ {msg}</span>"
            if a.xlabel != b.xlabel or a.kind != b.kind:
                html += ("<br><span style='color:#a60'>⚠ x-axis labels "
                         f"differ ('{a.xlabel}' vs '{b.xlabel}') – check "
                         "the units are the same.</span>")
            if op == "/" and np.any(b.y == 0):
                html += ("<br><span style='color:#a60'>⚠ B is zero at "
                         f"{int(np.sum(b.y == 0))} point(s); the result "
                         "will be undefined (NaN) there.</span>")
        else:
            expr = f"({a.name}) {sym} {c:g}"
            html = (f"<span style='color:#080'>✓ {a.npts} points, constant "
                    f"{c:g}</span>")
        if self._name_auto:
            self.name_edit.setText(expr)
        self._set_status(html, True)

    def create(self):
        try:
            a, op, b, c = self._inputs()
            name = self.name_edit.text().strip() or None
            ds, n_bad = apply_operation(a, op, b=b, const=c,
                                        name=name or "")
        except Exception as exc:  # noqa: BLE001
            QtWidgets.QMessageBox.warning(self, "Operation failed", str(exc))
            return
        self.win.settings.setValue("ops_after",
                                   self.after_combo.currentIndex())
        ds = self.win.add_derived(ds, self.after_combo.currentText())
        self._cache.clear()
        a_i, b_i = self.a_combo.currentIndex(), self.b_combo.currentIndex()
        self._load_candidates()
        self.a_combo.setCurrentIndex(a_i)
        self.b_combo.setCurrentIndex(b_i)
        self._name_auto = True
        self.update_status()
        extra = (f" ({n_bad} undefined point(s) set to NaN)"
                 if n_bad else "")
        self.status.setText(self.status.text()
                            + f"<br><b>Created '{ds.name}'{extra}.</b>")


# --------------------------------------------------------------------------
# Main window
# --------------------------------------------------------------------------
class MainWindow(QtWidgets.QMainWindow):
    TAB_FILES, TAB_DERIVED = 0, 1

    def __init__(self, folders: list[str] | None = None):
        super().__init__()
        self.setWindowTitle(APP_NAME)
        self.resize(1600, 920)
        self.settings = QtCore.QSettings("McStasTools", "McStas1DPlotter")
        self.folders: list[tuple[str, Path]] = []
        self.files: list[MonitorFile] = []
        self.derived: list[Dataset] = []
        self._plot_counter = 0
        self._syncing = False

        self._build_ui()
        self._build_menus()
        self.new_tab()

        if not folders:
            folders = self._stored_folders()
        added = False
        for f in folders:
            if Path(f).expanduser().is_dir():
                self.add_folder(f, rescan=False)
                added = True
        if added:
            self.rescan()
        else:
            self.statusBar().showMessage(
                "Open a McStas output folder (Ctrl+O); add more with "
                "'Add folder…' to compare simulations")

    def _stored_folders(self) -> list[str]:
        v = self.settings.value("folders", None)
        if v is None:
            last = self.settings.value("last_folder", "", type=str)
            return [last] if last else []
        if isinstance(v, str):
            return [v] if v else []
        return [str(x) for x in v]

    # ======================================================== UI building
    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        vlay = QtWidgets.QVBoxLayout(central)
        vlay.setContentsMargins(6, 6, 6, 6)

        # --- folder row ---
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Data folders:"))
        self.folder_combo = QtWidgets.QComboBox()
        self.folder_combo.setSizePolicy(QtWidgets.QSizePolicy.Expanding,
                                        QtWidgets.QSizePolicy.Fixed)
        self.folder_combo.setToolTip(
            "Loaded folders. Files are tagged with the letter of their "
            "folder (A, B, …)")
        row.addWidget(self.folder_combo, 1)
        open_btn = QtWidgets.QPushButton("Open folder…")
        open_btn.setToolTip("Replace all loaded folders with one folder "
                            "(Ctrl+O)")
        open_btn.clicked.connect(self.open_folder_dialog)
        add_btn = QtWidgets.QPushButton("Add folder…")
        add_btn.setToolTip("Load another folder alongside the current ones, "
                           "e.g. a second simulation to compare (Ctrl+Shift+O)")
        add_btn.clicked.connect(self.add_folder_dialog)
        rm_btn = QtWidgets.QPushButton("Remove")
        rm_btn.setToolTip("Remove the folder shown in the list")
        rm_btn.clicked.connect(self.remove_folder)
        for b in (open_btn, add_btn, rm_btn):
            row.addWidget(b)
        self.subdir_check = QtWidgets.QCheckBox("Include subfolders")
        self.subdir_check.setToolTip(
            "Also search sub-directories, e.g. the numbered point folders "
            "of a McStas scan")
        self.subdir_check.setChecked(
            self.settings.value("include_subfolders", True, type=bool))
        self.subdir_check.toggled.connect(self._on_subdir_toggled)
        row.addWidget(self.subdir_check)
        refresh = QtWidgets.QPushButton("Refresh")
        refresh.setToolTip("Re-scan all folders (F5)")
        refresh.clicked.connect(self.rescan)
        row.addWidget(refresh)
        vlay.addLayout(row)

        splitter = QtWidgets.QSplitter(Qt.Horizontal)
        vlay.addWidget(splitter, 1)

        # --- left: data browser ---
        left_split = QtWidgets.QSplitter(Qt.Vertical)
        data_box = QtWidgets.QGroupBox("Data")
        dl = QtWidgets.QVBoxLayout(data_box)
        self.left_tabs = QtWidgets.QTabWidget()
        dl.addWidget(self.left_tabs, 1)

        # files tab
        files_w = QtWidgets.QWidget()
        fl = QtWidgets.QVBoxLayout(files_w)
        fl.setContentsMargins(4, 4, 4, 4)
        frow = QtWidgets.QHBoxLayout()
        frow.addWidget(QtWidgets.QLabel("Show:"))
        self.tof_check = QtWidgets.QCheckBox("TOF")
        self.wl_check = QtWidgets.QCheckBox("Wavelength")
        self.other_check = QtWidgets.QCheckBox("Other 1D")
        for cb in (self.tof_check, self.wl_check, self.other_check):
            cb.setChecked(True)
            cb.toggled.connect(self.apply_filter)
            frow.addWidget(cb)
        frow.addStretch(1)
        fl.addLayout(frow)
        self.filter_edit = QtWidgets.QLineEdit()
        self.filter_edit.setPlaceholderText(
            "Filter by folder tag, file, component or title…  (e.g. 'B: tof')")
        self.filter_edit.setClearButtonEnabled(True)
        self.filter_edit.textChanged.connect(self.apply_filter)
        fl.addWidget(self.filter_edit)
        self.file_tree = QtWidgets.QTreeWidget()
        self.file_tree.setHeaderLabels(
            ["Set", "File", "Type", "Component", "Bins", "X axis", "Title"])
        self.file_tree.setRootIsDecorated(False)
        self.file_tree.setAlternatingRowColors(True)
        self.file_tree.setUniformRowHeights(True)
        self.file_tree.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection)
        self.file_tree.setSortingEnabled(True)
        self.file_tree.sortByColumn(0, Qt.AscendingOrder)
        self.file_tree.currentItemChanged.connect(self._on_data_current)
        self.file_tree.itemDoubleClicked.connect(lambda *_: self.overplot())
        self.file_tree.itemSelectionChanged.connect(self._update_buttons)
        fl.addWidget(self.file_tree, 1)
        self.left_tabs.addTab(files_w, "Files")

        # derived tab
        der_w = QtWidgets.QWidget()
        dvl = QtWidgets.QVBoxLayout(der_w)
        dvl.setContentsMargins(4, 4, 4, 4)
        self.derived_tree = QtWidgets.QTreeWidget()
        self.derived_tree.setHeaderLabels(
            ["Name", "Operation", "Points", "X axis"])
        self.derived_tree.setRootIsDecorated(False)
        self.derived_tree.setAlternatingRowColors(True)
        self.derived_tree.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection)
        self.derived_tree.currentItemChanged.connect(self._on_data_current)
        self.derived_tree.itemDoubleClicked.connect(
            lambda *_: self.overplot())
        self.derived_tree.itemSelectionChanged.connect(self._update_buttons)
        dvl.addWidget(self.derived_tree, 1)
        drow = QtWidgets.QHBoxLayout()
        ren = QtWidgets.QPushButton("Rename…")
        ren.clicked.connect(self.rename_derived)
        dele = QtWidgets.QPushButton("Delete")
        dele.clicked.connect(self.delete_derived)
        drow.addWidget(ren)
        drow.addWidget(dele)
        drow.addStretch(1)
        dvl.addLayout(drow)
        self.left_tabs.addTab(der_w, "Derived data (0)")
        self.left_tabs.currentChanged.connect(self._on_left_tab_changed)

        # action buttons (act on the active list)
        brow = QtWidgets.QHBoxLayout()
        self.err_check = QtWidgets.QCheckBox("Error bars")
        self.err_check.setToolTip(
            "Draw error bars for the curves you are about to plot.\n"
            "Can be changed per curve afterwards in 'Curves in current "
            "figure'.")
        self.err_check.setChecked(True)
        brow.addWidget(self.err_check)
        brow.addStretch(1)
        self.plot_new_btn = QtWidgets.QPushButton("Plot in new figure")
        self.plot_new_btn.setToolTip("Ctrl+N")
        self.plot_new_btn.clicked.connect(self.plot_new)
        self.overplot_btn = QtWidgets.QPushButton("Overplot on current")
        self.overplot_btn.setToolTip("Ctrl+Shift+N, or double-click")
        self.overplot_btn.clicked.connect(self.overplot)
        brow.addWidget(self.plot_new_btn)
        brow.addWidget(self.overplot_btn)
        dl.addLayout(brow)
        brow2 = QtWidgets.QHBoxLayout()
        self.ops_btn = QtWidgets.QPushButton("Data operations…")
        self.ops_btn.setToolTip(
            "Add, subtract, multiply, divide or scale datasets (Ctrl+M).\n"
            "The first two selected datasets are pre-filled as A and B.")
        self.ops_btn.clicked.connect(self.open_operations)
        self.save_btn = QtWidgets.QPushButton("Save data…")
        self.save_btn.setToolTip(
            "Save the selected dataset(s) in McStas 1D format or as "
            "3-column x, y, y_err text")
        self.save_btn.clicked.connect(self.save_selected)
        brow2.addWidget(self.ops_btn)
        brow2.addStretch(1)
        brow2.addWidget(self.save_btn)
        dl.addLayout(brow2)
        left_split.addWidget(data_box)

        info_box = QtWidgets.QGroupBox("Details")
        il = QtWidgets.QVBoxLayout(info_box)
        self.info_text = QtWidgets.QPlainTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
        self.info_text.setFont(
            QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont))
        il.addWidget(self.info_text)
        left_split.addWidget(info_box)
        left_split.setSizes([600, 200])
        splitter.addWidget(left_split)

        # --- centre: plot tabs ---
        self.tabs = QtWidgets.QTabWidget()
        self.tabs.setTabsClosable(True)
        self.tabs.setMovable(True)
        self.tabs.tabCloseRequested.connect(self.close_tab)
        self.tabs.currentChanged.connect(self._on_tab_changed)
        plus = QtWidgets.QToolButton()
        plus.setText(" + ")
        plus.setToolTip("New empty figure")
        plus.clicked.connect(lambda: self.new_tab())
        self.tabs.setCornerWidget(plus, Qt.TopRightCorner)
        splitter.addWidget(self.tabs)

        # --- right: controls ---
        splitter.addWidget(self._build_control_panel())
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        splitter.setStretchFactor(2, 0)
        splitter.setSizes([520, 760, 320])

        QtGui.QShortcut(QtGui.QKeySequence("F5"), self, self.rescan)
        QtGui.QShortcut(QtGui.QKeySequence("Ctrl+N"), self, self.plot_new)
        QtGui.QShortcut(QtGui.QKeySequence("Ctrl+Shift+N"), self,
                        self.overplot)
        QtGui.QShortcut(QtGui.QKeySequence("Ctrl+M"), self,
                        self.open_operations)
        self._update_buttons()

    def _build_control_panel(self) -> QtWidgets.QWidget:
        panel = QtWidgets.QWidget()
        pl = QtWidgets.QVBoxLayout(panel)
        pl.setContentsMargins(0, 0, 0, 0)

        axes_box = QtWidgets.QGroupBox("Axes (current figure)")
        g = QtWidgets.QGridLayout(axes_box)
        g.addWidget(QtWidgets.QLabel("Min"), 0, 1)
        g.addWidget(QtWidgets.QLabel("Max"), 0, 2)
        g.addWidget(QtWidgets.QLabel("Auto"), 0, 3)
        g.addWidget(QtWidgets.QLabel("Scale"), 0, 4)
        self.lim_edits = {}
        self.auto_checks = {}
        self.scale_combos = {}
        for r, axis in enumerate(("x", "y"), start=1):
            g.addWidget(QtWidgets.QLabel(axis.upper()), r, 0)
            lo = QtWidgets.QLineEdit()
            hi = QtWidgets.QLineEdit()
            for w in (lo, hi):
                w.setMinimumWidth(70)
                w.setPlaceholderText("auto")
                w.returnPressed.connect(self.apply_limits)
                w.setEnabled(False)
            g.addWidget(lo, r, 1)
            g.addWidget(hi, r, 2)
            auto = QtWidgets.QCheckBox()
            auto.setChecked(True)
            auto.toggled.connect(
                lambda checked, a=axis: self._on_auto_toggled(a, checked))
            g.addWidget(auto, r, 3, alignment=Qt.AlignCenter)
            combo = QtWidgets.QComboBox()
            combo.addItems(["linear", "log"])
            combo.currentTextChanged.connect(
                lambda text, a=axis: self._on_scale_changed(a, text))
            g.addWidget(combo, r, 4)
            self.lim_edits[axis] = (lo, hi)
            self.auto_checks[axis] = auto
            self.scale_combos[axis] = combo
        brow = QtWidgets.QHBoxLayout()
        apply_btn = QtWidgets.QPushButton("Apply limits")
        apply_btn.clicked.connect(self.apply_limits)
        from_btn = QtWidgets.QPushButton("From plot")
        from_btn.setToolTip(
            "Copy the limits currently shown (e.g. after zooming with the "
            "toolbar) into the fields and keep them fixed")
        from_btn.clicked.connect(self.limits_from_plot)
        auto_btn = QtWidgets.QPushButton("Autoscale")
        auto_btn.clicked.connect(self.autoscale)
        for b in (apply_btn, from_btn, auto_btn):
            brow.addWidget(b)
        g.addLayout(brow, 3, 0, 1, 5)
        pl.addWidget(axes_box)

        disp_box = QtWidgets.QGroupBox("Display (current figure)")
        f = QtWidgets.QFormLayout(disp_box)
        self.style_combo = QtWidgets.QComboBox()
        self.style_combo.addItems(STYLES)
        self.style_combo.currentTextChanged.connect(self._on_display_changed)
        f.addRow("Style:", self.style_combo)
        self.norm_check = QtWidgets.QCheckBox("Normalise each curve to max")
        self.legend_check = QtWidgets.QCheckBox("Legend")
        self.grid_check = QtWidgets.QCheckBox("Grid")
        for cb in (self.norm_check, self.legend_check, self.grid_check):
            cb.toggled.connect(self._on_display_changed)
            f.addRow(cb)
        pl.addWidget(disp_box)

        curves_box = QtWidgets.QGroupBox("Curves in current figure")
        cl = QtWidgets.QVBoxLayout(curves_box)
        hint = QtWidgets.QLabel(
            "Tick to show/hide · tick 'Errors' for error bars · "
            "double-click a label to rename")
        hint.setWordWrap(True)
        hint.setStyleSheet("color: gray; font-size: 11px;")
        cl.addWidget(hint)
        self.curve_tree = QtWidgets.QTreeWidget()
        self.curve_tree.setHeaderLabels(["Curve", "Errors"])
        self.curve_tree.setRootIsDecorated(False)
        self.curve_tree.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection)
        self.curve_tree.setEditTriggers(
            QtWidgets.QAbstractItemView.DoubleClicked
            | QtWidgets.QAbstractItemView.EditKeyPressed)
        hdr = self.curve_tree.header()
        hdr.setStretchLastSection(False)
        hdr.setSectionResizeMode(0, QtWidgets.QHeaderView.Stretch)
        hdr.setSectionResizeMode(1, QtWidgets.QHeaderView.ResizeToContents)
        self.curve_tree.itemChanged.connect(self._on_curve_item_changed)
        cl.addWidget(self.curve_tree, 1)
        r1 = QtWidgets.QHBoxLayout()
        b_on = QtWidgets.QPushButton("All errors on")
        b_on.clicked.connect(lambda: self.set_all_errorbars(True))
        b_off = QtWidgets.QPushButton("All errors off")
        b_off.clicked.connect(lambda: self.set_all_errorbars(False))
        r1.addWidget(b_on)
        r1.addWidget(b_off)
        cl.addLayout(r1)
        r2 = QtWidgets.QHBoxLayout()
        b_rm = QtWidgets.QPushButton("Remove selected")
        b_rm.clicked.connect(self.remove_selected_curves)
        b_clr = QtWidgets.QPushButton("Clear figure")
        b_clr.clicked.connect(self.clear_plot)
        r2.addWidget(b_rm)
        r2.addWidget(b_clr)
        cl.addLayout(r2)
        pl.addWidget(curves_box, 1)

        scroll = QtWidgets.QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QtWidgets.QFrame.NoFrame)
        scroll.setWidget(panel)
        scroll.setMinimumWidth(300)
        self.control_panel = panel
        return scroll

    def _build_menus(self):
        mb = self.menuBar()
        fm = mb.addMenu("&File")
        a = fm.addAction("Open data folder…", self.open_folder_dialog)
        a.setShortcut(QtGui.QKeySequence.Open)
        a = fm.addAction("Add data folder…", self.add_folder_dialog)
        a.setShortcut(QtGui.QKeySequence("Ctrl+Shift+O"))
        fm.addAction("Refresh file lists", self.rescan)
        fm.addSeparator()
        fm.addAction("Save selected data…", self.save_selected)
        a = fm.addAction("Save current figure…", self.save_figure)
        a.setShortcut(QtGui.QKeySequence.Save)
        fm.addSeparator()
        a = fm.addAction("Quit", self.close)
        a.setShortcut(QtGui.QKeySequence.Quit)

        dm = mb.addMenu("&Data")
        dm.addAction("Data operations…", self.open_operations)
        dm.addAction("Save selected data…", self.save_selected)

        pm = mb.addMenu("&Plot")
        pm.addAction("Plot selection in new figure", self.plot_new)
        pm.addAction("Overplot selection on current figure", self.overplot)
        pm.addSeparator()
        pm.addAction("New empty figure", lambda: self.new_tab())
        pm.addAction("Clear current figure", self.clear_plot)

        hm = mb.addMenu("&Help")
        hm.addAction("Usage", self.show_help)
        hm.addAction("About", lambda: QtWidgets.QMessageBox.about(
            self, "About",
            f"<b>{APP_NAME}</b> v{__version__}<br>"
            "Plots, compares and combines McStas/McCode 1D monitor text "
            "files (TOF, wavelength and other 1D monitors).<br><br>"
            f"Qt (PySide6) {QtCore.qVersion()} · "
            f"matplotlib {matplotlib.__version__} · numpy {np.__version__}"))

    # ======================================================== folders
    def _next_tag(self) -> str:
        used = {t for t, _ in self.folders}
        for t in string.ascii_uppercase:
            if t not in used:
                return t
        i = 27
        while f"F{i}" in used:
            i += 1
        return f"F{i}"

    def _start_dir(self) -> str:
        if self.folders:
            return str(self.folders[-1][1].parent)
        return str(Path.home())

    def open_folder_dialog(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Open McStas data folder", self._start_dir())
        if d:
            self.folders.clear()
            self.add_folder(d)

    def add_folder_dialog(self):
        d = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Add McStas data folder for comparison", self._start_dir())
        if d:
            self.add_folder(d)

    def add_folder(self, folder: str, rescan: bool = True):
        p = Path(folder).expanduser()
        if not p.is_dir():
            QtWidgets.QMessageBox.warning(
                self, "Not a folder", f"'{folder}' is not a directory.")
            return
        p = p.resolve()
        if any(p == q for _, q in self.folders):
            self.statusBar().showMessage(f"{p} is already loaded", 5000)
            return
        self.folders.append((self._next_tag(), p))
        self._update_folder_combo()
        if rescan:
            self.rescan()

    def remove_folder(self):
        i = self.folder_combo.currentIndex()
        if 0 <= i < len(self.folders):
            self.folders.pop(i)
            self._update_folder_combo()
            self.rescan()

    def _update_folder_combo(self):
        self.folder_combo.clear()
        for tag, p in self.folders:
            self.folder_combo.addItem(f"{tag}:  {p}")
            self.folder_combo.setItemData(self.folder_combo.count() - 1,
                                          str(p), Qt.ToolTipRole)
        self.folder_combo.setCurrentIndex(self.folder_combo.count() - 1)
        self.settings.setValue("folders", [str(p) for _, p in self.folders])

    def _on_subdir_toggled(self, checked):
        self.settings.setValue("include_subfolders", checked)
        self.rescan()

    def rescan(self):
        self.files = []
        msgs = []
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            for tag, folder in self.folders:
                if not folder.is_dir():
                    msgs.append(f"{tag}: folder missing")
                    continue
                found, truncated = scan_folder(
                    folder, self.subdir_check.isChecked(), tag,
                    progress_cb=QtWidgets.QApplication.processEvents)
                self.files.extend(found)
                counts = {k: sum(f.kind == k for f in found)
                          for k in (KIND_TOF, KIND_WL, KIND_OTHER)}
                m = (f"{tag}: {len(found)} files (TOF {counts[KIND_TOF]}, "
                     f"λ {counts[KIND_WL]}, other {counts[KIND_OTHER]})")
                if truncated:
                    m += " – scan truncated, choose a more specific folder"
                msgs.append(m)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self._populate_tree()
        self.statusBar().showMessage(
            "   |   ".join(msgs) if msgs else "No folders loaded")

    def file_label(self, mf: MonitorFile) -> str:
        """Label used for a file in legends and operations."""
        if len(self.folders) > 1:
            return f"{mf.tag}:{mf.rel}"
        return mf.rel

    def _populate_tree(self):
        self.file_tree.setSortingEnabled(False)
        self.file_tree.clear()
        for mf in self.files:
            self.file_tree.addTopLevelItem(FileItem(mf))
        self.file_tree.setSortingEnabled(True)
        for col in range(6):
            self.file_tree.resizeColumnToContents(col)
        self.file_tree.setColumnWidth(
            FileItem.COL_FILE,
            min(self.file_tree.columnWidth(FileItem.COL_FILE), 240))
        self.apply_filter()
        if self.left_tabs.currentIndex() == self.TAB_FILES:
            self.info_text.clear()
        self._update_buttons()

    def apply_filter(self):
        show = {KIND_TOF: self.tof_check.isChecked(),
                KIND_WL: self.wl_check.isChecked(),
                KIND_OTHER: self.other_check.isChecked()}
        text = self.filter_edit.text().strip().lower()
        for i in range(self.file_tree.topLevelItemCount()):
            it = self.file_tree.topLevelItem(i)
            mf = it.mf
            ok = show.get(mf.kind, True)
            if ok and text:
                hay = " ".join([f"{mf.tag}:{mf.rel}", mf.component, mf.title,
                                mf.xlabel, mf.params]).lower()
                ok = all(tok in hay for tok in text.split())
            it.setHidden(not ok)
            if not ok:
                it.setSelected(False)
        self._update_buttons()

    # ======================================================== selection
    def _on_left_tab_changed(self, _i):
        tree = (self.file_tree if self.left_tabs.currentIndex()
                == self.TAB_FILES else self.derived_tree)
        self._on_data_current(tree.currentItem(), None)
        self._update_buttons()

    def _on_data_current(self, current, _previous):
        if current is None:
            self.info_text.clear()
        elif isinstance(current, FileItem):
            self.info_text.setPlainText(current.mf.header_text)
        elif isinstance(current, DerivedItem):
            self.info_text.setPlainText(current.ds.description)

    def _update_buttons(self):
        has = bool(self.selected_objects())
        for b in (self.plot_new_btn, self.overplot_btn, self.save_btn):
            b.setEnabled(has)

    def selected_objects(self) -> list:
        """Selected MonitorFile or Dataset objects in the active list."""
        if self.left_tabs.currentIndex() == self.TAB_FILES:
            items = [it for it in self.file_tree.selectedItems()
                     if not it.isHidden()]
            items.sort(key=self.file_tree.indexOfTopLevelItem)
            return [it.mf for it in items]
        items = list(self.derived_tree.selectedItems())
        items.sort(key=self.derived_tree.indexOfTopLevelItem)
        return [it.ds for it in items]

    def dataset_for(self, obj) -> Dataset:
        if isinstance(obj, Dataset):
            return obj
        return dataset_from_file(obj, self.file_label(obj))

    def collect_selected(self) -> tuple[list[Dataset], list[str]]:
        datasets, errors = [], []
        for obj in self.selected_objects():
            try:
                datasets.append(self.dataset_for(obj))
            except Exception as exc:  # noqa: BLE001
                name = obj.rel if isinstance(obj, MonitorFile) else obj.name
                errors.append(f"{name}: {exc}")
        return datasets, errors

    def operand_candidates(self) -> list[tuple[str, object]]:
        """(display text, object) for every file and derived dataset."""
        out = []
        for mf in sorted(self.files,
                         key=lambda m: (m.tag, natural_key(m.rel))):
            npts = f"{mf.nbins} pts" if mf.nbins else "? pts"
            out.append((f"{self.file_label(mf)}   [{mf.kind}, {npts}]", mf))
        for ds in self.derived:
            out.append((f"{ds.name}   [derived, {ds.npts} pts]", ds))
        return out

    # ======================================================== plotting
    def current_tab(self) -> PlotTab | None:
        w = self.tabs.currentWidget()
        return w if isinstance(w, PlotTab) else None

    def new_tab(self, title: str | None = None) -> PlotTab:
        self._plot_counter += 1
        tab = PlotTab()
        idx = self.tabs.addTab(tab, title or f"Figure {self._plot_counter}")
        self.tabs.setCurrentIndex(idx)
        return tab

    def close_tab(self, index: int):
        w = self.tabs.widget(index)
        self.tabs.removeTab(index)
        if w is not None:
            w.figure.clear()
            w.deleteLater()
        if self.tabs.count() == 0:
            self.new_tab()

    def plot_new(self, datasets: list[Dataset] | None = None):
        if not datasets:
            datasets, errors = self.collect_selected()
            self._report_errors(errors)
        if not datasets:
            self.statusBar().showMessage("Nothing selected", 4000)
            return
        tab = self.current_tab()
        if tab is None or tab.curves:   # re-use an empty figure
            tab = self.new_tab()
        self._add_datasets(tab, datasets)

    def overplot(self, datasets: list[Dataset] | None = None):
        if not datasets:
            datasets, errors = self.collect_selected()
            self._report_errors(errors)
        if not datasets:
            self.statusBar().showMessage("Nothing selected", 4000)
            return
        tab = self.current_tab() or self.new_tab()
        existing = tab.kinds()
        new = {d.kind for d in datasets}
        if existing and (new | existing) != existing:
            r = QtWidgets.QMessageBox.question(
                self, "Different x-axis quantities",
                "The current figure contains "
                f"{', '.join(sorted(existing))} data and you are adding "
                f"{', '.join(sorted(new))} data.\n\n"
                "The x axes probably have different units. Overplot anyway?",
                QtWidgets.QMessageBox.Yes | QtWidgets.QMessageBox.No,
                QtWidgets.QMessageBox.No)
            if r != QtWidgets.QMessageBox.Yes:
                return
        self._add_datasets(tab, datasets)

    def _add_datasets(self, tab: PlotTab, datasets: list[Dataset]):
        err_on = self.err_check.isChecked()
        for ds in datasets:
            tab.add_curve(Curve.from_dataset(ds, err_on))
        tab.redraw()
        self._refresh_curve_list()
        self.statusBar().showMessage(
            f"Plotted {len(datasets)} curve(s) in "
            f"'{self.tabs.tabText(self.tabs.indexOf(tab))}'", 5000)

    def _report_errors(self, errors):
        if errors:
            QtWidgets.QMessageBox.warning(
                self, "Some data could not be read", "\n".join(errors))

    # ======================================================== derived data
    def open_operations(self):
        objs = self.selected_objects()
        pre_a = objs[0] if objs else None
        pre_b = objs[1] if len(objs) > 1 else None
        if not self.files and not self.derived:
            QtWidgets.QMessageBox.information(
                self, "No data", "Load a data folder first.")
            return
        dlg = OperationsDialog(self, pre_a, pre_b)
        dlg.exec()

    def _unique_name(self, name: str) -> str:
        names = {d.name for d in self.derived}
        if name not in names:
            return name
        i = 2
        while f"{name} ({i})" in names:
            i += 1
        return f"{name} ({i})"

    def add_derived(self, ds: Dataset, after: str = "Don't plot") -> Dataset:
        ds.name = self._unique_name(ds.name)
        self.derived.append(ds)
        item = DerivedItem(ds)
        self.derived_tree.addTopLevelItem(item)
        for col in range(4):
            self.derived_tree.resizeColumnToContents(col)
        self.derived_tree.setColumnWidth(
            0, min(self.derived_tree.columnWidth(0), 260))
        self._update_derived_tab_title()
        self.left_tabs.setCurrentIndex(self.TAB_DERIVED)
        self.derived_tree.clearSelection()
        self.derived_tree.setCurrentItem(item)
        item.setSelected(True)
        if after.startswith("Plot in new"):
            self.plot_new([ds])
        elif after.startswith("Overplot"):
            self.overplot([ds])
        return ds

    def _update_derived_tab_title(self):
        self.left_tabs.setTabText(self.TAB_DERIVED,
                                  f"Derived data ({len(self.derived)})")

    def rename_derived(self):
        it = self.derived_tree.currentItem()
        if it is None:
            return
        name, ok = QtWidgets.QInputDialog.getText(
            self, "Rename dataset", "New name:", text=it.ds.name)
        name = name.strip()
        if ok and name and name != it.ds.name:
            it.ds.name = self._unique_name(name)
            it.setText(0, it.ds.name)

    def delete_derived(self):
        items = self.derived_tree.selectedItems()
        if not items:
            return
        r = QtWidgets.QMessageBox.question(
            self, "Delete derived data",
            f"Delete {len(items)} derived dataset(s)? Curves already "
            "plotted are not affected.")
        if r != QtWidgets.QMessageBox.Yes:
            return
        for it in items:
            self.derived.remove(it.ds)
            self.derived_tree.takeTopLevelItem(
                self.derived_tree.indexOfTopLevelItem(it))
        self._update_derived_tab_title()
        self._update_buttons()

    # ======================================================== saving
    FMT_MCSTAS = "McStas 1D format (*.dat)"
    FMT_3COL = "3-column x y y_err (*.txt)"

    def save_selected(self):
        datasets, errors = self.collect_selected()
        self._report_errors(errors)
        if not datasets:
            self.statusBar().showMessage("Select data to save first", 4000)
            return
        start = Path(self.settings.value(
            "save_dir", str(self.folders[0][1]) if self.folders
            else str(Path.home()), type=str))
        if len(datasets) == 1:
            ds = datasets[0]
            last = self.settings.value("save_filter", self.FMT_MCSTAS,
                                       type=str)
            ext = ".txt" if last == self.FMT_3COL else ".dat"
            name, flt = QtWidgets.QFileDialog.getSaveFileName(
                self, f"Save '{ds.name}'",
                str(start / (safe_filename(ds.name) + ext)),
                f"{self.FMT_MCSTAS};;{self.FMT_3COL}", last)
            if not name:
                return
            path = Path(name)
            fmt = "3col" if flt == self.FMT_3COL else "mcstas"
            if not path.suffix:
                path = path.with_suffix(".txt" if fmt == "3col" else ".dat")
            self.settings.setValue("save_filter", flt)
            jobs = [(ds, path)]
        else:
            fmt_txt, ok = QtWidgets.QInputDialog.getItem(
                self, "Save data", f"Format for {len(datasets)} datasets:",
                [self.FMT_MCSTAS, self.FMT_3COL], 0, False)
            if not ok:
                return
            d = QtWidgets.QFileDialog.getExistingDirectory(
                self, "Folder to save into", str(start))
            if not d:
                return
            fmt = "3col" if fmt_txt == self.FMT_3COL else "mcstas"
            ext = ".txt" if fmt == "3col" else ".dat"
            jobs, used = [], set()
            for ds in datasets:
                stem = safe_filename(ds.name)
                cand, i = stem, 2
                while cand in used:
                    cand = f"{stem}_{i}"
                    i += 1
                used.add(cand)
                jobs.append((ds, Path(d) / (cand + ext)))
            exist = [p for _, p in jobs if p.exists()]
            if exist:
                r = QtWidgets.QMessageBox.question(
                    self, "Overwrite?",
                    f"{len(exist)} file(s) already exist, e.g.\n"
                    f"{exist[0]}\n\nOverwrite?")
                if r != QtWidgets.QMessageBox.Yes:
                    return
        written, errs = [], []
        for ds, path in jobs:
            try:
                if fmt == "3col":
                    write_three_column(path, ds)
                else:
                    write_mcstas_1d(path, ds)
                written.append(path)
            except Exception as exc:  # noqa: BLE001
                errs.append(f"{path}: {exc}")
        if written:
            self.settings.setValue("save_dir", str(written[0].parent))
            self.statusBar().showMessage(
                f"Saved {len(written)} file(s) to {written[0].parent}", 8000)
        self._report_errors(errs)

    # ======================================================== tab <-> panel
    def _on_tab_changed(self, _index):
        self._sync_panel_from_tab()
        self._refresh_curve_list()

    def _sync_panel_from_tab(self):
        tab = self.current_tab()
        self.control_panel.setEnabled(tab is not None)
        if tab is None:
            return
        st = tab.state
        self._syncing = True
        try:
            for axis in ("x", "y"):
                lim = st.xlim if axis == "x" else st.ylim
                lo, hi = self.lim_edits[axis]
                self.auto_checks[axis].setChecked(lim is None)
                lo.setEnabled(lim is not None)
                hi.setEnabled(lim is not None)
                lo.setText("" if lim is None else f"{lim[0]:.6g}")
                hi.setText("" if lim is None else f"{lim[1]:.6g}")
                self.scale_combos[axis].setCurrentText(
                    st.xscale if axis == "x" else st.yscale)
            self.style_combo.setCurrentText(st.style)
            self.norm_check.setChecked(st.normalise)
            self.legend_check.setChecked(st.legend)
            self.grid_check.setChecked(st.grid)
        finally:
            self._syncing = False

    def _refresh_curve_list(self):
        tab = self.current_tab()
        self.curve_tree.blockSignals(True)
        self.curve_tree.clear()
        if tab is not None:
            for i, c in enumerate(tab.curves):
                it = QtWidgets.QTreeWidgetItem([c.label, ""])
                it.setData(0, Qt.UserRole, i)
                it.setFlags(it.flags() | Qt.ItemIsUserCheckable
                            | Qt.ItemIsEditable)
                it.setCheckState(0, Qt.Checked if c.visible else Qt.Unchecked)
                it.setCheckState(1, Qt.Checked if c.errorbars
                                 else Qt.Unchecked)
                pix = QtGui.QPixmap(14, 14)
                pix.fill(QtGui.QColor(matplotlib.colors.to_hex(c.color)))
                it.setIcon(0, QtGui.QIcon(pix))
                it.setToolTip(0, c.source)
                self.curve_tree.addTopLevelItem(it)
        self.curve_tree.blockSignals(False)

    def _on_curve_item_changed(self, item, column):
        tab = self.current_tab()
        if tab is None:
            return
        idx = item.data(0, Qt.UserRole)
        if idx is None or idx >= len(tab.curves):
            return
        c = tab.curves[idx]
        if column == 0:
            c.visible = item.checkState(0) == Qt.Checked
            label = item.text(0).strip()
            if label:
                c.label = label
        elif column == 1:
            c.errorbars = item.checkState(1) == Qt.Checked
        tab.redraw()

    def set_all_errorbars(self, on: bool):
        tab = self.current_tab()
        if tab is None:
            return
        for c in tab.curves:
            c.errorbars = on
        tab.redraw()
        self._refresh_curve_list()

    def remove_selected_curves(self):
        tab = self.current_tab()
        if tab is None:
            return
        idxs = {it.data(0, Qt.UserRole)
                for it in self.curve_tree.selectedItems()}
        if not idxs:
            self.statusBar().showMessage("Select curve(s) in the list first",
                                         4000)
            return
        tab.curves = [c for i, c in enumerate(tab.curves) if i not in idxs]
        tab.redraw()
        self._refresh_curve_list()

    def clear_plot(self):
        tab = self.current_tab()
        if tab is None:
            return
        tab.curves.clear()
        tab._color_index = 0
        tab.redraw()
        self._refresh_curve_list()

    # ======================================================== axes controls
    def _on_auto_toggled(self, axis: str, checked: bool):
        lo, hi = self.lim_edits[axis]
        lo.setEnabled(not checked)
        hi.setEnabled(not checked)
        if self._syncing:
            return
        tab = self.current_tab()
        if tab is None:
            return
        if checked:
            setattr(tab.state, f"{axis}lim", None)
            lo.clear()
            hi.clear()
            tab.redraw()
        else:
            cur = tab.ax.get_xlim() if axis == "x" else tab.ax.get_ylim()
            if not lo.text():
                lo.setText(f"{cur[0]:.6g}")
            if not hi.text():
                hi.setText(f"{cur[1]:.6g}")
            lo.setFocus()

    def _on_scale_changed(self, axis: str, text: str):
        if self._syncing:
            return
        tab = self.current_tab()
        if tab is None:
            return
        setattr(tab.state, f"{axis}scale", text)
        lim = getattr(tab.state, f"{axis}lim")
        if text == "log" and lim is not None and lim[0] <= 0:
            setattr(tab.state, f"{axis}lim", None)
            self._sync_panel_from_tab()
            self.statusBar().showMessage(
                f"{axis.upper()} limits reset to auto: log scale needs "
                "positive limits", 6000)
        tab.redraw()

    def apply_limits(self):
        tab = self.current_tab()
        if tab is None:
            return
        new = {}
        for axis in ("x", "y"):
            if self.auto_checks[axis].isChecked():
                new[axis] = None
                continue
            lo, hi = self.lim_edits[axis]
            try:
                vmin = float(lo.text())
                vmax = float(hi.text())
            except ValueError:
                QtWidgets.QMessageBox.warning(
                    self, "Invalid limits",
                    f"{axis.upper()} limits must be numbers "
                    "(e.g. 1000, 2.5e4).")
                return
            if not vmin < vmax:
                QtWidgets.QMessageBox.warning(
                    self, "Invalid limits",
                    f"{axis.upper()} min must be smaller than max.")
                return
            if getattr(tab.state, f"{axis}scale") == "log" and vmin <= 0:
                QtWidgets.QMessageBox.warning(
                    self, "Invalid limits",
                    f"{axis.upper()} axis is logarithmic: limits must be > 0.")
                return
            new[axis] = (vmin, vmax)
        tab.state.xlim = new["x"]
        tab.state.ylim = new["y"]
        tab.redraw()

    def limits_from_plot(self):
        tab = self.current_tab()
        if tab is None:
            return
        tab.state.xlim = tuple(tab.ax.get_xlim())
        tab.state.ylim = tuple(tab.ax.get_ylim())
        self._sync_panel_from_tab()
        tab.redraw()

    def autoscale(self):
        tab = self.current_tab()
        if tab is None:
            return
        tab.state.xlim = None
        tab.state.ylim = None
        self._sync_panel_from_tab()
        tab.redraw()

    def _on_display_changed(self, *_):
        if self._syncing:
            return
        tab = self.current_tab()
        if tab is None:
            return
        st = tab.state
        st.style = self.style_combo.currentText()
        st.normalise = self.norm_check.isChecked()
        st.legend = self.legend_check.isChecked()
        st.grid = self.grid_check.isChecked()
        if st.normalise and st.ylim is not None:
            st.ylim = None
            self._sync_panel_from_tab()
        tab.redraw()

    # ======================================================== misc
    def save_figure(self):
        tab = self.current_tab()
        if tab is None:
            return
        start = self.folders[0][1] if self.folders else Path.home()
        name, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save figure", str(start / "mcstas_plot.png"),
            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg);;EPS (*.eps)")
        if name:
            tab.figure.savefig(name, dpi=200)
            self.statusBar().showMessage(f"Saved {name}", 5000)

    def show_help(self):
        QtWidgets.QMessageBox.information(
            self, "Usage",
            "<b>1.</b> <i>Open folder…</i> (Ctrl+O) loads a McStas output "
            "folder; <i>Add folder…</i> (Ctrl+Shift+O) loads more folders to "
            "compare. Files are tagged A, B, C… by folder.<br>"
            "<b>2.</b> Select files (Ctrl/Shift-click) and click <i>Plot in "
            "new figure</i> or <i>Overplot on current</i>. The <i>Error "
            "bars</i> box applies to newly added curves.<br>"
            "<b>3.</b> <i>Data operations…</i> (Ctrl+M) computes A+B, A−B, "
            "A×B, A÷B or A with a constant. Datasets must have the same "
            f"number of points and x values within {100 * X_TOLERANCE:g} %. "
            "The first two selected datasets are pre-filled as A and B. "
            "Results go to the <i>Derived data</i> tab, where they can be "
            "plotted, combined further, renamed or deleted.<br>"
            "<b>4.</b> <i>Save data…</i> writes the selected dataset(s) in "
            "McStas 1D format or as 3-column x, y, y_err text.<br>"
            "<b>5.</b> In <i>Curves in current figure</i>: tick to show/hide, "
            "tick <i>Errors</i> per curve, double-click to rename.<br>"
            "<b>6.</b> Set axis limits (untick <i>Auto</i>, type values, "
            "Enter) and linear/log scales. <i>From plot</i> fixes the "
            "toolbar zoom.<br>"
            "<b>7.</b> Save figures with the toolbar disk icon or Ctrl+S.")

    def closeEvent(self, event):
        self.settings.sync()
        super().closeEvent(event)


def main():
    QtWidgets.QApplication.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName(APP_NAME)
    win = MainWindow(sys.argv[1:] or None)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
