#!/usr/bin/env python
"""
McStas 1D Monitor Plotter
=========================

A Qt (PySide6) + matplotlib application for browsing and plotting the
text output files written by McStas / McCode 1D monitors, in particular
time-of-flight (TOF_monitor, Monitor_nD "t") and wavelength
(L_monitor, Monitor_nD "lambda") monitors.

Features
--------
* Choose a data folder (optionally including sub-folders, e.g. scan points
  0/, 1/, 2/ ...). All McStas 1D monitor files found are listed with their
  type (TOF / Wavelength / Other 1D), component name, number of bins and
  axis label.
* Plot the selected file(s) in a new figure (tab) or overplot them on the
  current figure.
* Error bars on/off per curve; curves can be hidden, renamed (double-click
  the label) or removed.
* Axis limits (manual or automatic) and linear/log scale for each axis.
* Line, marker, line+marker or histogram-step drawing, optional
  normalisation to the maximum, legend and grid toggles.
* Standard matplotlib toolbar (zoom, pan, save figure as PNG/PDF/SVG...).

Usage
-----
    python mcstas_plotter.py [data_folder]
"""

from __future__ import annotations

import os
import re
import sys
import warnings
from dataclasses import dataclass, field
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
__version__ = "1.0.0"

KIND_TOF = "TOF"
KIND_WL = "Wavelength"
KIND_OTHER = "Other 1D"

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


# --------------------------------------------------------------------------
# File model
# --------------------------------------------------------------------------
@dataclass
class MonitorFile:
    path: Path
    rel: str
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


def parse_monitor_header(path: Path, root: Path) -> MonitorFile | None:
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
        meta=meta,
        header_text="".join(lines),
        kind=classify(meta),
        nbins=int(nb.group(1)) if nb else None,
    )


def scan_folder(folder: Path, recursive: bool,
                progress_cb=None) -> tuple[list[MonitorFile], bool]:
    """Find McStas 1D monitor files. Returns (files, truncated)."""
    results: list[MonitorFile] = []
    n_checked = 0
    truncated = False
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
                truncated = True
                return results, truncated
            if progress_cb and n_checked % 200 == 0:
                progress_cb()
            mf = parse_monitor_header(p, folder)
            if mf is not None:
                results.append(mf)
    return results, truncated


def load_monitor_data(mf: MonitorFile):
    """Return x, I, I_err, N arrays from a McStas 1D monitor file."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            arr = np.loadtxt(mf.path, comments="#", ndmin=2)
        except ValueError:
            arr = np.genfromtxt(mf.path, comments="#")
            arr = np.atleast_2d(arr)
    if arr.size == 0:
        raise ValueError("file contains no data rows")
    ncol = arr.shape[1]
    if ncol < 2:
        raise ValueError(f"expected at least 2 columns, found {ncol}")
    x = arr[:, 0]
    y = arr[:, 1]
    e = arr[:, 2] if ncol > 2 else np.zeros_like(y)
    n = arr[:, 3] if ncol > 3 else None
    return x, y, np.abs(e), n


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

    # ---- curves ---------------------------------------------------------
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

    # ---- drawing --------------------------------------------------------
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
            lower = np.where(y - e > 0, e, 0.95 * y)
            lower = np.where(y > 0, lower, 0.0)
            yerr = np.vstack([lower, e])

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

        if st.xscale == "log":
            ax.set_xscale("log", nonpositive="mask")
        else:
            ax.set_xscale("linear")
        if st.yscale == "log":
            ax.set_yscale("log", nonpositive="mask")
        else:
            ax.set_yscale("linear")

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
                    "Select monitor files on the left and click\n"
                    "'Plot in new figure' or 'Overplot on current'",
                    transform=ax.transAxes, ha="center", va="center",
                    color="0.55")
        self.canvas.draw_idle()


# --------------------------------------------------------------------------
# Tree items
# --------------------------------------------------------------------------
class FileItem(QtWidgets.QTreeWidgetItem):
    COL_FILE, COL_KIND, COL_COMP, COL_BINS, COL_X, COL_TITLE = range(6)

    def __init__(self, mf: MonitorFile):
        super().__init__([mf.rel, mf.kind, mf.component,
                          "" if mf.nbins is None else str(mf.nbins),
                          mf.xlabel, mf.title])
        self.mf = mf
        tip = f"{mf.path}"
        if mf.params:
            tip += f"\nParameters: {mf.params}"
        for col in range(6):
            self.setToolTip(col, tip)
        self.setTextAlignment(self.COL_BINS, Qt.AlignRight | Qt.AlignVCenter)

    def __lt__(self, other):
        tree = self.treeWidget()
        col = tree.sortColumn() if tree else 0
        if col == self.COL_BINS:
            a = self.mf.nbins or 0
            b = other.mf.nbins or 0
            return a < b
        return natural_key(self.text(col)) < natural_key(other.text(col))


# --------------------------------------------------------------------------
# Main window
# --------------------------------------------------------------------------
class MainWindow(QtWidgets.QMainWindow):
    def __init__(self, folder: str | None = None):
        super().__init__()
        self.setWindowTitle(APP_NAME)
        self.resize(1550, 900)
        self.settings = QtCore.QSettings("McStasTools", "McStas1DPlotter")
        self.folder: Path | None = None
        self.files: list[MonitorFile] = []
        self._plot_counter = 0
        self._syncing = False

        self._build_ui()
        self._build_menus()
        self.new_tab()

        start = folder or self.settings.value("last_folder", "", type=str)
        if start and Path(start).is_dir():
            self.set_folder(start)
        else:
            self.statusBar().showMessage(
                "Choose a McStas output folder with 'Browse…' (Ctrl+O)")

    # ======================================================== UI building
    def _build_ui(self):
        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        vlay = QtWidgets.QVBoxLayout(central)
        vlay.setContentsMargins(6, 6, 6, 6)

        # --- folder row ---
        row = QtWidgets.QHBoxLayout()
        row.addWidget(QtWidgets.QLabel("Data folder:"))
        self.folder_edit = QtWidgets.QLineEdit()
        self.folder_edit.setPlaceholderText(
            "Path to a McStas output directory (press Enter to load)")
        self.folder_edit.returnPressed.connect(
            lambda: self.set_folder(self.folder_edit.text()))
        row.addWidget(self.folder_edit, 1)
        browse = QtWidgets.QPushButton("Browse…")
        browse.clicked.connect(self.choose_folder)
        row.addWidget(browse)
        self.subdir_check = QtWidgets.QCheckBox("Include subfolders")
        self.subdir_check.setToolTip(
            "Also search sub-directories, e.g. the numbered point folders "
            "of a McStas scan")
        self.subdir_check.setChecked(
            self.settings.value("include_subfolders", True, type=bool))
        self.subdir_check.toggled.connect(self._on_subdir_toggled)
        row.addWidget(self.subdir_check)
        refresh = QtWidgets.QPushButton("Refresh")
        refresh.setToolTip("Re-scan the folder (F5)")
        refresh.clicked.connect(self.rescan)
        row.addWidget(refresh)
        vlay.addLayout(row)

        splitter = QtWidgets.QSplitter(Qt.Horizontal)
        vlay.addWidget(splitter, 1)

        # --- left: file browser ---
        left_split = QtWidgets.QSplitter(Qt.Vertical)
        files_box = QtWidgets.QGroupBox("Available 1D monitor files")
        fl = QtWidgets.QVBoxLayout(files_box)

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
            "Filter by file, component or title…")
        self.filter_edit.setClearButtonEnabled(True)
        self.filter_edit.textChanged.connect(self.apply_filter)
        fl.addWidget(self.filter_edit)

        self.file_tree = QtWidgets.QTreeWidget()
        self.file_tree.setHeaderLabels(
            ["File", "Type", "Component", "Bins", "X axis", "Title"])
        self.file_tree.setRootIsDecorated(False)
        self.file_tree.setAlternatingRowColors(True)
        self.file_tree.setUniformRowHeights(True)
        self.file_tree.setSelectionMode(
            QtWidgets.QAbstractItemView.ExtendedSelection)
        self.file_tree.setSortingEnabled(True)
        self.file_tree.sortByColumn(0, Qt.AscendingOrder)
        self.file_tree.currentItemChanged.connect(self._on_file_current)
        self.file_tree.itemDoubleClicked.connect(lambda *_: self.overplot())
        self.file_tree.itemSelectionChanged.connect(self._update_buttons)
        fl.addWidget(self.file_tree, 1)

        brow = QtWidgets.QHBoxLayout()
        self.err_check = QtWidgets.QCheckBox("Error bars")
        self.err_check.setToolTip(
            "Draw error bars for the curves you are about to plot.\n"
            "Can be changed per curve afterwards in 'Curves in current plot'.")
        self.err_check.setChecked(True)
        brow.addWidget(self.err_check)
        brow.addStretch(1)
        self.plot_new_btn = QtWidgets.QPushButton("Plot in new figure")
        self.plot_new_btn.setToolTip("Ctrl+N")
        self.plot_new_btn.clicked.connect(self.plot_new)
        self.overplot_btn = QtWidgets.QPushButton("Overplot on current")
        self.overplot_btn.setToolTip("Ctrl+Shift+N, or double-click a file")
        self.overplot_btn.clicked.connect(self.overplot)
        brow.addWidget(self.plot_new_btn)
        brow.addWidget(self.overplot_btn)
        fl.addLayout(brow)
        left_split.addWidget(files_box)

        info_box = QtWidgets.QGroupBox("File header")
        il = QtWidgets.QVBoxLayout(info_box)
        self.info_text = QtWidgets.QPlainTextEdit()
        self.info_text.setReadOnly(True)
        self.info_text.setLineWrapMode(QtWidgets.QPlainTextEdit.NoWrap)
        self.info_text.setFont(
            QtGui.QFontDatabase.systemFont(QtGui.QFontDatabase.FixedFont))
        il.addWidget(self.info_text)
        left_split.addWidget(info_box)
        left_split.setSizes([560, 220])
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
        splitter.setSizes([480, 760, 320])

        QtGui.QShortcut(QtGui.QKeySequence("F5"), self, self.rescan)
        QtGui.QShortcut(QtGui.QKeySequence("Ctrl+N"), self, self.plot_new)
        QtGui.QShortcut(QtGui.QKeySequence("Ctrl+Shift+N"), self,
                        self.overplot)
        self._update_buttons()

    def _build_control_panel(self) -> QtWidgets.QWidget:
        panel = QtWidgets.QWidget()
        pl = QtWidgets.QVBoxLayout(panel)
        pl.setContentsMargins(0, 0, 0, 0)

        # --- Axes ---
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
            lo.setEnabled(False)
            hi.setEnabled(False)
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
        brow.addWidget(apply_btn)
        brow.addWidget(from_btn)
        brow.addWidget(auto_btn)
        g.addLayout(brow, 3, 0, 1, 5)
        pl.addWidget(axes_box)

        # --- Display ---
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

        # --- Curves ---
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
        a = fm.addAction("Open data folder…", self.choose_folder)
        a.setShortcut(QtGui.QKeySequence.Open)
        fm.addAction("Refresh file list", self.rescan)
        fm.addSeparator()
        a = fm.addAction("Save current figure…", self.save_figure)
        a.setShortcut(QtGui.QKeySequence.Save)
        fm.addSeparator()
        a = fm.addAction("Quit", self.close)
        a.setShortcut(QtGui.QKeySequence.Quit)

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
            "Plots McStas/McCode 1D monitor text files "
            "(TOF, wavelength and other 1D monitors).<br><br>"
            f"Qt (PySide6) {QtCore.qVersion()} · "
            f"matplotlib {matplotlib.__version__} · numpy {np.__version__}"))

    # ======================================================== folder / files
    def choose_folder(self):
        start = str(self.folder) if self.folder else str(Path.home())
        d = QtWidgets.QFileDialog.getExistingDirectory(
            self, "Select McStas data folder", start)
        if d:
            self.set_folder(d)

    def set_folder(self, folder: str):
        p = Path(folder).expanduser()
        if not p.is_dir():
            QtWidgets.QMessageBox.warning(
                self, "Not a folder", f"'{folder}' is not a directory.")
            return
        self.folder = p.resolve()
        self.folder_edit.setText(str(self.folder))
        self.settings.setValue("last_folder", str(self.folder))
        self.rescan()

    def _on_subdir_toggled(self, checked):
        self.settings.setValue("include_subfolders", checked)
        self.rescan()

    def rescan(self):
        if self.folder is None:
            return
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            self.files, truncated = scan_folder(
                self.folder, self.subdir_check.isChecked(),
                progress_cb=QtWidgets.QApplication.processEvents)
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self._populate_tree()
        n = len(self.files)
        counts = {k: sum(f.kind == k for f in self.files)
                  for k in (KIND_TOF, KIND_WL, KIND_OTHER)}
        msg = (f"{n} 1D monitor file(s) found in {self.folder}  "
               f"[TOF: {counts[KIND_TOF]}, wavelength: {counts[KIND_WL]}, "
               f"other: {counts[KIND_OTHER]}]")
        if truncated:
            msg += (f"  — stopped after checking {MAX_FILES_SCANNED} files; "
                    "choose a more specific folder")
        self.statusBar().showMessage(msg)

    def _populate_tree(self):
        self.file_tree.setSortingEnabled(False)
        self.file_tree.clear()
        for mf in self.files:
            self.file_tree.addTopLevelItem(FileItem(mf))
        self.file_tree.setSortingEnabled(True)
        for col in range(5):
            self.file_tree.resizeColumnToContents(col)
        self.file_tree.setColumnWidth(
            0, min(self.file_tree.columnWidth(0), 260))
        self.apply_filter()
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
                hay = " ".join([mf.rel, mf.component, mf.title,
                                mf.xlabel, mf.params]).lower()
                ok = all(tok in hay for tok in text.split())
            it.setHidden(not ok)
            if not ok:
                it.setSelected(False)
        self._update_buttons()

    def _on_file_current(self, current, _previous):
        if current is None:
            self.info_text.clear()
            return
        self.info_text.setPlainText(current.mf.header_text)

    def _update_buttons(self):
        has = bool(self.selected_files())
        self.plot_new_btn.setEnabled(has)
        self.overplot_btn.setEnabled(has)

    def selected_files(self) -> list[MonitorFile]:
        items = [it for it in self.file_tree.selectedItems()
                 if not it.isHidden()]
        items.sort(key=self.file_tree.indexOfTopLevelItem)
        return [it.mf for it in items]

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

    def plot_new(self):
        files = self.selected_files()
        if not files:
            self.statusBar().showMessage("No files selected", 4000)
            return
        tab = self.current_tab()
        # Re-use the current figure if it is still empty
        if tab is None or tab.curves:
            tab = self.new_tab()
        self._add_files(tab, files)

    def overplot(self):
        files = self.selected_files()
        if not files:
            self.statusBar().showMessage("No files selected", 4000)
            return
        tab = self.current_tab() or self.new_tab()
        existing = tab.kinds()
        new = {f.kind for f in files}
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
        self._add_files(tab, files)

    def _add_files(self, tab: PlotTab, files: list[MonitorFile]):
        errors = []
        added = 0
        err_on = self.err_check.isChecked()
        QtWidgets.QApplication.setOverrideCursor(Qt.WaitCursor)
        try:
            for mf in files:
                try:
                    x, y, e, _n = load_monitor_data(mf)
                except Exception as exc:  # noqa: BLE001
                    errors.append(f"{mf.rel}: {exc}")
                    continue
                tab.add_curve(Curve(label=mf.rel, x=x, y=y, e=e,
                                    xlabel=mf.xlabel, ylabel=mf.ylabel,
                                    kind=mf.kind, source=str(mf.path),
                                    errorbars=err_on))
                added += 1
            tab.redraw()
        finally:
            QtWidgets.QApplication.restoreOverrideCursor()
        self._refresh_curve_list()
        if added:
            self.statusBar().showMessage(
                f"Plotted {added} curve(s) in "
                f"'{self.tabs.tabText(self.tabs.indexOf(tab))}'", 5000)
        if errors:
            QtWidgets.QMessageBox.warning(
                self, "Some files could not be read", "\n".join(errors))

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
            # Pre-fill with the limits currently displayed
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
        name, _ = QtWidgets.QFileDialog.getSaveFileName(
            self, "Save figure",
            str((self.folder or Path.home()) / "mcstas_plot.png"),
            "PNG (*.png);;PDF (*.pdf);;SVG (*.svg);;EPS (*.eps)")
        if name:
            tab.figure.savefig(name, dpi=200)
            self.statusBar().showMessage(f"Saved {name}", 5000)

    def show_help(self):
        QtWidgets.QMessageBox.information(
            self, "Usage",
            "<b>1.</b> Choose a McStas output folder (Browse… or Ctrl+O). "
            "Tick <i>Include subfolders</i> for scans.<br>"
            "<b>2.</b> Select one or more files (Ctrl/Shift-click). "
            "Use the TOF / Wavelength / Other checkboxes and the filter box "
            "to narrow the list.<br>"
            "<b>3.</b> Click <i>Plot in new figure</i> or <i>Overplot on "
            "current</i> (double-click a file also overplots). The "
            "<i>Error bars</i> box sets error bars for newly added curves."
            "<br><b>4.</b> In <i>Curves in current figure</i>: tick to "
            "show/hide, tick <i>Errors</i> to toggle error bars per curve, "
            "double-click to rename.<br>"
            "<b>5.</b> Set axis limits (untick <i>Auto</i>, type values, "
            "press Enter or <i>Apply limits</i>) and choose linear/log "
            "scales. <i>From plot</i> fixes the current toolbar zoom.<br>"
            "<b>6.</b> Save with the toolbar disk icon or Ctrl+S.")

    def closeEvent(self, event):
        self.settings.sync()
        super().closeEvent(event)


def main():
    QtWidgets.QApplication.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough)
    app = QtWidgets.QApplication(sys.argv)
    app.setApplicationName(APP_NAME)
    folder = sys.argv[1] if len(sys.argv) > 1 else None
    win = MainWindow(folder)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
