# McStas 1D Monitor Plotter

A desktop application (Qt/PySide6 + matplotlib) to browse and plot the text
output files written by McStas/McCode **1D monitors**: time-of-flight
(`TOF_monitor`, `Monitor_nD` with `t`), wavelength (`L_monitor`,
`Monitor_nD` with `lambda`) and any other 1D monitor (energy, position, …).

## Contents

| File | Purpose |
|---|---|
| `mcstas_plotter.py` | The application |
| `environment.yml` | Conda environment definition (conda-forge only) |
| `install.sh` | Installer for Linux / macOS |
| `install.bat` | Installer for Windows |
| `make_example_data.py` | Writes synthetic McStas files for testing |
| `example_data/` | Example output: a 4-point scan, a Monitor_nD file, an energy monitor, plus a 2D file and `mccode.sim` that are (correctly) ignored |

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
run_plotter.bat                        # Windows (or double-click it)
```
To try it out, open `example_data` and tick **Include subfolders**.

The last folder you used and the *Include subfolders* setting are remembered
between sessions.

---

## Using the application

**Left panel – files**
* **Browse… / Ctrl+O** selects a data folder. You can also type a path and
  press Enter. **Refresh / F5** re-scans it, e.g. while a simulation is still
  running.
* **Include subfolders** also searches sub-directories, such as the
  numbered point folders `0/ 1/ 2/ …` of a McStas scan (`-N` option).
* Every file is checked by its header, not its extension. Any file
  with `# type: array_1d(N)` is listed. 2D monitors, `mccode.sim`,
  `mccode.dat` scan summaries and non-McStas files are ignored.
* Columns: file (relative path), **type** (TOF / Wavelength / Other 1D,
  worked out from `xvar`, `xlabel` and `title`), component, number of
  bins, x-axis label and title. Click a header to sort.
* The **Show** check-boxes filter by type. The text box filters by file
  name, component, title or scan parameter, e.g. `tof 3/`.
* Selecting a file shows its full header in **File header**. Hover over a
  row to see its full path and scan parameters.
* Select one or more files (Ctrl/Shift-click), then:
  * **Plot in new figure** (Ctrl+N) opens a new figure tab. An empty
    current figure is reused.
  * **Overplot on current** (Ctrl+Shift+N, or double-click a file) adds
    the curves to the current figure. You are warned if you mix
    TOF and wavelength data.
  * **Error bars** sets whether the curves you are adding get error bars,
    so one figure can hold curves with and without them.

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
