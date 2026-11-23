# PopBio Studio — User Guide (`v1.2.0`)

PopBio Studio is released under the **GNU General Public License v3.0 or later (`GPL-3.0-or-later`)** as part of the CRV Science open-source observatory network ([`https://crvscience.com/apps`](https://crvscience.com/apps)).

## 1. Navigation & Global Hotkeys

All single-key workspace shortcuts check `!ctx.wants_keyboard_input()`, so typing in the alignment editor, Front Menu search bar, or SQLite checkpoint slot field behaves normally without triggering simulation commands.

| Control | Scope | Action |
| --- | --- | --- |
| `1` – `4` | Workspace | Select Tabs `1..4` on `V1..T6` (`Visualization`, `Graphs`, `Phase/Diagnostics`, `Equations/CSV`) or Views `E1..E4` on empirical datasets |
| `Space` | Workspace | Play / Pause simulation (`V1..T6`) or Wright-Fisher haplotype evolution (`E2`) |
| `.` (Period) | Workspace | Pause and advance exactly 1 time step or 1 `E2` generation |
| `N` | Dynamical (`V1..T6`) | Advance step-by-step until next biological growth phase transition |
| `R` | Workspace | Pause and reset active module or empirical haplotype state to baseline |
| `F1` | Global | Toggle Shortcuts & Citations window |
| `F5` / `F9` | Global | Quicksave / Quickload all 19 modules to SQLite slot `"Quicksave"` |
| `Esc` | Global | Close `About` -> `Help`, then return to Front Menu |

## 2. 3D Adaptive Landscapes & Spatial Arenas

- **Sewall Wright 3D Adaptive Landscape (`T5`)**:
  - **Orbit**: Left-drag on the 3D viewport (or adjust `Azimuth` and `Elevation` sliders).
  - **Zoom**: Scroll the mouse wheel over the 3D surface (`0.55x..2.40x`) or use the `Zoom` slider.
  - **Reset 3D View**: Click `Reset 3D View` (in either the scenario bar or viewport footer) to restore `(azimuth = -38 deg, elevation = 30 deg, zoom = 1.0x)`.
  - **Presets**: Compare `Valley Crossing (N_e=20)` (drift-assisted saddle traversal across Peaks 1 & 2) vs. `Peak Trapped (N_e=350)`.
- **Spatial Patch Hover Telemetry (`I2`, `M3`, `T1`)**:
  - Hovering any orange patch in **Huffaker Mites (`I2`)** reveals exact prey (`Eotetranychus`) density vs. carrying capacity `K`, predator density, and Vaseline barrier fraction.
  - Hovering any agar cell in **Phage T4 Plaque (`M3`)** shows local `E. coli` bacterial lawn percentage and virion density.
  - Hovering any habitat node in **Levins Metapopulation (`T1`)** displays patch index, occupancy/destruction status, and connected degree.

## 3. Empirical PopGen, Bayesian Skyline & Sliding-Window Workflows (`E1`–`E4`)

1. **Loading Data**: Select one of the built-in FASTA benchmarks in the left panel or click **Load FASTA / PHYLIP / NEXUS File...** (tag subpopulations in FASTA headers using `>SampleID|DemeName` or `>SampleID [Pop=DemeName]`).
2. **Switching Views**: Use `1`–`4` or the top sub-tabs (`1. Global PopGen`, `2. Haplotype Network`, `3. Bayesian Skyline`, `4. Sliding-Window Scan`).
3. **Interactive Log-Y & Adaptive Ticks**:
   - In **`E3` Bayesian Skyline**, toggle **Log10 Ne(t) Y-Axis** to compare multi-order-of-magnitude bottlenecks/expansions cleanly; linear `N_e` up to `10^6` and tiny `pi ~ 5e-4` render cleanly in scientific notation (`1.00e6`, `5.00e-4`).
   - In dynamical time-series tabs, click the **Log-Y** checkbox in the plot legend strip for positive series.

## 4. Headless CLI & Portable Worksets

- Export portable SHA-256 checksummed `.popworkset.json` files or SQLite slots (`popbio_worksets.sqlite3`) for reproducible lab handouts.
- Run automated analyses headlessly:
  ```powershell
  cargo run --release -- simulate --module v1-hare-lynx --steps 600 --out hare_lynx.csv
  cargo run --release -- analyze-fasta examples/chromosomal_sweep_balancing_180bp.fasta --bsp-groups 4 --bsp-mcmc-steps 3000 --out-stats popgen.json --out-skyline-csv skyline.csv --out-windows-csv windows.csv --out-svg mj_network.svg
  ```
