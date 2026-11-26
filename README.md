# McStas 1D Monitor Plotter

A desktop application (Qt/PySide6 + matplotlib) to browse, compare,
combine and plot the text output files written by McStas/McCode
**1D monitors**: time-of-flight (`TOF_monitor`, `Monitor_nD` with `t`),
wavelength (`L_monitor`, `Monitor_nD` with `lambda`) and any other 1D
monitor (energy, position, …).

**New in 1.1**
* You can load several data folders (A, B, C…) at once to compare
  simulations.
* Data operations: A + B, A − B, A × B and A ÷ B between datasets, and
  offsetting or scaling by a constant. Uncertainties are propagated.
* Results go into a *Derived data* list, where they can be plotted like
  any other data and used in further operations.
* Any dataset can be saved in McStas 1D format or as 3-column
  `x y y_err` text.

## Contents

| File | Purpose |
|---|---|
| `mcstas_plotter.py` | The application |
| `environment.yml` | Conda environment definition (conda-forge only) |
| `install.sh` | Installer for Linux / macOS |
| `install.bat` | Installer for Windows |
| `make_example_data.py` | Writes synthetic McStas files for testing |
| `example_data/` | Two example "simulations", `sim_A` (60 K moderator) and `sim_B` (40 K). Each has a 4-point scan, a Monitor_nD file and an energy monitor, plus a 2D file and `mccode.sim` that are (correctly) ignored. `sim_B` also has a 200-bin wavelength monitor that does not match `sim_A`'s 100 bins, to demonstrate the x-axis check. |

---

## Installation

### 1. Install conda (skip if you already have it)

Install **Miniforge**, the conda-forge distribution of conda/mamba, from
<https://conda-forge.org/download/>. It is free to use in any organisation,
unlike the Anaconda `defaults` channel, which needs a commercial licence for
organisations with more than 200 staff. The environment file only uses
conda-forge.

* **Linux/macOS:**
  ```bash
  curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
  bash Miniforge3-$(uname)-$(uname -m).sh
  ```
  Accept the defaults and allow it to initialise your shell, then open a new
  terminal.
* **Windows:** run the `.exe` installer, then use **Miniforge Prompt** from
  the Start menu for the steps below.

### 2. Run the installer

Unpack `mcstas_plotter.zip` somewhere permanent (the launcher records this
path), then:

**Linux / macOS**
```bash
cd mcstas_plotter
chmod +x install.sh
./install.sh                # or: ./install.sh my-env-name
```

**Windows** (in Miniforge Prompt)
```bat
cd path\to\mcstas_plotter
install.bat                 REM or: install.bat my-env-name
```

The installer:
1. creates the conda environment `mcstas-plotter` (or updates it if it
   already exists), using `mamba` if available;
2. checks that numpy, matplotlib and PySide6 import correctly;
3. writes a launcher, `run_plotter.sh` or `run_plotter.bat`, that runs the
   app with the environment's Python, so you don't need to activate
   anything;
4. creates `example_data/` if it is missing.

### Manual installation (equivalent)

```bash
conda env create -f environment.yml          # once
conda activate mcstas-plotter
python mcstas_plotter.py [path/to/mcstas/output]
```
To update later: `conda env update -f environment.yml --prune`.
To uninstall: `conda env remove -n mcstas-plotter`, then delete the folder.

---

## Running

```bash
./run_plotter.sh                       # Linux/macOS
./run_plotter.sh /data/sims/larmor_run # open a folder directly
./run_plotter.sh /data/runA /data/runB # open two folders to compare
run_plotter.bat                        # Windows (or double-click it)
```
To try it out, run
`./run_plotter.sh example_data/sim_A example_data/sim_B`.

The loaded folders and the *Include subfolders* setting are remembered
between sessions.

---

## Using the application

**Folders**
* **Open folder… (Ctrl+O)** loads one folder and replaces any folders
  already loaded.
* **Add folder… (Ctrl+Shift+O)** loads another folder alongside the
  current ones, e.g. a second simulation to compare. Each folder gets a
  tag (A, B, C, …), shown in the **Set** column.
* When more than one folder is loaded, files are labelled `A:scan/0/tof_mon.t`,
  `B:scan/0/tof_mon.t` and so on. The same labels are used in legends and
  data operations.
* **Remove** drops the folder currently shown in the drop-down.
* **Refresh / F5** re-scans all folders, e.g. while a simulation is
  still running.
* **Include subfolders** also searches sub-directories, such as the
  numbered point folders `0/ 1/ 2/ …` of a McStas scan (`-N` option).
* Every file is checked by its header, not its extension. Any file
  with `# type: array_1d(N)` is listed. 2D monitors, `mccode.sim`,
  `mccode.dat` scan summaries and non-McStas files are ignored.
**Left panel – Files tab**
* Columns: set (folder tag), file (relative path), **type** (TOF / Wavelength / Other 1D,
  worked out from `xvar`, `xlabel` and `title`), component, number of
  bins, x-axis label and title. Click a header to sort.
* The **Show** check-boxes filter by type. The text box filters by folder
  tag, file name, component, title or scan parameter, e.g. `B: tof` or
  `tof 3/`.
* Selecting a file shows its full header in **Details**. Hover over a
  row to see its full path and scan parameters.
* Select one or more files (Ctrl/Shift-click). The buttons below the list
  act on whatever is selected in the active tab, Files or Derived data:
  * **Plot in new figure** (Ctrl+N) opens a new figure tab. An empty
    current figure is reused.
  * **Overplot on current** (Ctrl+Shift+N, or double-click a file) adds
    the curves to the current figure. You are warned if you mix
    TOF and wavelength data.
  * **Error bars** sets whether the curves you are adding get error bars,
    so one figure can hold curves with and without them.
  * **Data operations… (Ctrl+M)** opens the operations dialog (see below).
  * **Save data…** writes the selected dataset(s) to disk (see below).

**Data operations**

Click **Data operations…** or press Ctrl+M. The first two datasets
selected in the list are pre-filled as A and B.

| Operation | with dataset B | with a constant c |
|---|---|---|
| + | y = a + b,  σ = √(σa² + σb²) | offset: y = a + c, σ = σa |
| − | y = a − b,  σ = √(σa² + σb²) | y = a − c, σ = σa |
| × | y = a·b,  σ = √((b σa)² + (a σb)²) | **scale**: y = c·a, σ = \|c\| σa |
| ÷ | y = a/b,  σ = \|y\|·√((σa/a)² + (σb/b)²) | y = a/c, σ = σa/\|c\| |

* **Compatibility check.** Datasets can only be combined if they
  have the **same number of points** and their **x values agree to
  better than 0.1 %** point by point.
  * Deviation is measured as |xa − xb| / max(|xa|, |xb|).
  * For points close to x = 0, the full x range is used as the reference
    instead.
  * The dialog shows the result live, and **Create** is disabled with
    the reason if the check fails.
  * You are also warned if the x-axis labels differ (e.g. µs vs s), or
    if B has zeros when dividing. Points where B is zero become NaN and
    show as gaps in the plot.
* **Assumptions.** Errors assume A and B are uncorrelated (independent
  simulations), and constants are treated as exact. The result takes
  A's x values and A's x label. The y label becomes *Ratio* or
  *Product* for ÷ and × of two datasets.
* **Naming.** The result is named automatically from the expression,
  e.g. `(A:scan/0/lambda_mon.L) ÷ (B:scan/0/lambda_mon.L)`. You can type
  your own name in the dialog.
* **After creating.** Choose *Plot in new figure*, *Overplot on current
  figure* or *Don't plot*. The dialog stays open, so you can chain
  operations, e.g. divide and then scale the result.

**Left panel – Derived data tab**
* Lists every result of the operations dialog. Selecting one shows its
  operation, its operands and their source files in **Details**.
* Derived datasets can be plotted, overplotted (double-click works too),
  used as A or B in further operations, saved, renamed or deleted.
* They exist only for the current session, so **save the ones you want
  to keep**.

**Saving data**

Select one or more datasets, either files or derived data, and click
**Save data…**. With one dataset you choose the file name and format.
With several you choose a format and a folder, and file names come from
the dataset names. Two formats are available:
* **McStas 1D format (`.dat`).** This is the standard McCode text header
  plus `x I I_err N` columns. It includes `xvar`/`xlabel`, statistics and
  a `# Derived:` line recording the operation. The file can be read by
  McStas tools (mcplot) and by this program: add its folder and it
  appears in the file list.
* **3-column text (`.txt`).** This is `x y y_err`, with a short `#`
  comment header giving the name, source and axis labels. It can be read
  directly by `numpy.loadtxt`, gnuplot, Excel and so on.

Headers are written in plain ASCII (Å → AA, ÷ → /, …) for portability.
Undefined points are written as `nan`.

**Centre – figures**
* Each figure is a tab. Use **+** for an empty one and **×** to close.
* The matplotlib toolbar gives zoom, pan, home, subplot spacing, and
  **save** (PNG/PDF/SVG/EPS). **File → Save current figure (Ctrl+S)** does
  the same at 200 dpi.

**Right panel – controls for the current figure**
* **Axes**
  * Untick **Auto** for X or Y and type **Min/Max**, then press Enter or
    **Apply limits**.
  * The **Scale** drop-down switches each axis between *linear* and
    *log*. On log axes, zero or negative points are masked. A lower error
    bar that would cross zero is clipped so it can still be drawn.
  * **From plot** fixes the limits currently shown, e.g. after a toolbar
    zoom, so they survive later changes such as toggling error bars.
  * **Autoscale** returns both axes to automatic limits.
* **Display**
  * Style: Line, Markers, Line + markers or Steps (histogram).
  * Normalise each curve to its maximum, which is useful for comparing
    spectral shapes.
  * Turn the legend and grid on or off.
* **Curves in current figure**
  * The tick box beside a name shows or hides that curve.
  * The **Errors** tick box turns error bars on or off for that curve.
  * Double-click a name to rename it in the legend.
  * **All errors on/off**, **Remove selected** and **Clear figure** act
    on the curves in the list.

Each figure keeps its own axis, scale and display settings, and the panel
updates when you switch tabs.

---

## Supported file format

The McCode text format written by McStas 2.x and 3.x, for example:

```
# Format: McCode with text headers
# type: array_1d(200)
# component: tof_mon
# title: Time-of-flight monitor
# xvar: t
# xlabel: Time-of-flight [\gms]
# ylabel: Intensity
# variables: t I I_err N
0.5   1.23e-08  3.1e-10  123
...
```
Column 1 is x, column 2 is the intensity I, column 3 is the error on I,
and column 4 (N, optional) is the event count. McStas label codes are
converted for display (`\gms`→µs, `\gl`→λ, `[AA]`/`[Angs]`→Å).
`Monitor_nD` TOF output is in **seconds**, while `TOF_monitor` output is
in **µs**. The axis label shows which one you have.

---

## Troubleshooting

| Problem | Fix |
|---|---|
| Linux: `qt.qpa.plugin: Could not load the Qt platform plugin "xcb"` | Qt ≥ 6.5 needs `libxcb-cursor0`. Run `sudo apt install libxcb-cursor0` on Debian/Ubuntu, or `sudo dnf install xcb-util-cursor` on RHEL/Rocky/Fedora. |
| Running over SSH | Use `ssh -X`/`-Y` or a NoMachine/VNC session. For quick tests, `QT_QPA_PLATFORM=offscreen` runs the app without a display. |
| Wayland glitches | Start with `QT_QPA_PLATFORM=xcb ./run_plotter.sh`. |
| `conda: command not found` (Windows) | Use **Miniforge Prompt**, or run `conda init cmd.exe` / `conda init powershell` once. |
| Environment solve is slow | Install `mamba` (included in Miniforge). The installer uses it automatically. |
| A monitor file isn't listed | Check that its header contains `# type: array_1d(...)`. Untick no type filters and clear the text filter. Tick *Include subfolders* if the file is in a sub-directory. |
| A file is listed as "Other 1D" instead of TOF/Wavelength | Classification is by `xvar`/`xlabel`/`title`. It still plots normally. |
| Moved the folder after installing | Re-run the installer to regenerate the launcher. |
