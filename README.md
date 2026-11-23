# Geothermal Data Analysis: natural-state interpretation

This repository is for analysing geothermal reservoir data with AI agents, using Claude Code.

Its current process is to take downhole data and well geometry from a field and interpret it into **natural-state calibration targets** and a **conceptual model**. The intended user is a reservoir engineer, or another agent, who is calibrating a numerical reservoir model.

`CLAUDE.md` holds the binding rules for agents. This README explains how the repository is laid out, how to run the analysis, and **what format input data must have** for the process to work.

---

## 1. Repository layout

```
GDA Files/
├── CLAUDE.md          Rules for agents (scope, data vs reference, required checks, deliverables)
├── README.md          Repository overview
├── Baseline/          Input data. READ-ONLY: never edited by the analysis
│   └── <Project>/     e.g. Teaching_Model_3D
│       ├── Data/      Field data: the factual inputs
│       ├── Models/    Reference model: geometry is used, everything else is reference only
│       └── Parameters/ Production history spreadsheets (not used for natural state)
├── Scripts/           Analysis scripts (natural_state_analysis.py)
└── Output/            All generated tables, figures, reports and logs
    └── figures/
```

Git history keeps earlier runs, for example one each by Opus, Sonnet and Haiku. A new run starts from `Baseline/` and does not reuse earlier outputs.

## 2. Running the analysis

Run this from the repository root:

```powershell
python Scripts/natural_state_analysis.py
```

**Requirements:** Python 3 with numpy and matplotlib.

**Not required:** pandas, scipy and iapws. The script implements IAPWS-IF97 itself:
- Region 4 for saturation;
- Region 1 for liquid density.

**What the script does:** it regenerates everything listed in §5. Console output is saved to `Output/run_console_log.txt`.

**What the agent does as well:**
- writes `Output/NS_data_interpretation.md`;
- writes `Output/process_log.md`, containing the prompt and all steps taken.

## 3. Which inputs are used

| Role | Files | How they are used |
|---|---|---|
| **Data** (factual) | `Data/Data_json/WELL-*_data.json`, `Data/Mapping/Mapping.json`, `Data/Groups_data.json` | Primary evidence |
| **Cross-check** | `Data/Data_dat/*.dat` | Must match the JSON exactly |
| **Geometry** | `Models/<model>/g<model>.dat` (the `VERTICES`, `GRID`, `LAYERS`, `SURFA` and `WELLS` sections), `Models/<model>/g<model>_wells.json` | Mapping wells to model blocks |
| **Reference only** | Rock types, boundary conditions, natural-state `.dat` and `.incon`, `current*Setup.json`, `Data/Plotting/`, `Parameters/` | Describes the starting model only. Never used as evidence |
| **Timing check only** | The first value of each transient pressure series | Tests whether the downhole surveys are undisturbed |

## 4. Required input data formats

A new dataset must follow these formats if the existing script is to run on it unchanged. Files go under `Baseline/<Project>/`.

### 4.1 Well file: `Data/Data_json/WELL-<n>_data.json`

There is one file per well. Here is the top-level structure, with fields the natural-state analysis does not need left out:

```jsonc
{
  "Well_Name": "WELL-2",                  // must match Groups_data.json
  "AUT2_Name": "WELL2",                   // well name in the model (g<model>.dat WELLS section)
  "Feed_Zones": {
    "Units": {"Top": "m", "Bottom": "m", "Quality": "Text", "Source": "Text"},
    "Data": ["[-250, -270, 'LOW']", "[-720, -770, 'HIGH']"]
  },
  "Well_Track": {
    "Units": {"x": "m", "y": "m", "z": "masl"},
    "Data": ["[2819527.0, 6335494.0, 225.0]", "[2819527.0, 6335494.0, -2348.0]"]
  },
  "Well_Design": {
    "Units": {"Top": "mRL", "Bottom": "mRL", "Diameter": "\"", "Slotted": "Boolean", "Spec": "Text"},
    "Data": ["[225, -100, '13-3/8', False, 'Casing 13-3/8\"']",
             "[-50, -2348, '8-5/8', True, 'Liner 8-5/8\"']"]
  },
  "Downhole_Temperature": {
    "29/07/2015": {                        // survey date, dd/mm/yyyy
      "Units": {"Date": "dd/mm/yyyy", "Elevation": "mRL", "Temperature": "degC"},
      "Source": "", "Comment": "",
      "StableFlag": true,
      "Data": {"225.0": 16.5, "210.0": 19.3}   // elevation (string key) -> value
    }
  },
  "Downhole_Pressure": { "29/07/2015": { "...": "same structure, value in bar" } },

  // Optional. Used only for the natural-state timing check:
  "Transient_Pressure": {
    "-750.0": {                            // gauge elevation, mRL
      "Units": {"Date": "dd/mm/yyyy", "Pressure": "bar"},
      "Data": {"01/01/2014": 98.3}
    }
  }
}
```

Formatting rules the parser relies on:

- **Row format.** Every row of `Feed_Zones`, `Well_Track` and `Well_Design` is a **string containing a Python literal**, so it may use single quotes and `True`/`False`. Parse it with `ast.literal_eval`, not with a JSON parser.
- **Elevations, not depths.** Every elevation is in **mRL (metres relative to sea level, positive upwards)**. This includes feed-zone `Top` and `Bottom`, even though they are labelled "m". Measured depths are not supported.
- **Feed-zone order.** A feed zone is `[top, bottom, quality]` with top > bottom. Quality is one of `LOW`, `MEDIUM`, `HIGH` or `Synthetic`.
- **Well track.**
  - Points are `[x, y, z]` in the model's projected coordinate system (NZGM here), listed from the wellhead down.
  - A vertical well needs only two points.
  - A deviated well lists every survey station, and positions between stations are linearly interpolated.
- **Casing and liner.** A casing or liner row is `[top, bottom, diameter, slotted, spec]`. `slotted = True` marks the open or permeable interval.
- **Survey data.**
  - Survey `Data` maps elevation to value, with the elevation as a string key.
  - T and P surveys should share the same elevations.
  - Units are °C and bar. State whether pressure is absolute or gauge in `Comment`; the current dataset appears to be absolute.
- **Stable surveys.** `StableFlag: true` marks a survey taken after heat-up. Only stable surveys are used as natural-state targets.
- **Date format.** Dates are `dd/mm/yyyy` throughout.

### 4.2 Plain-text copies: `Data/Data_dat/*.dat`

These are two whitespace-separated columns in scientific notation, with no header.

| File pattern | Column 1 | Column 2 |
|---|---|---|
| `<Well>_downhole_temperature.dat` | Temperature (°C) | Elevation (mRL) |
| `<Well>_downhole_pressure.dat` | Pressure (bar) | Elevation (mRL) |
| `<Well>_<elev>_transient_pressure.dat` | Time (years since 01/01/2014) | Pressure (bar) |
| `<Well>_production_mass_flow.dat`, `_injection_mass_flow.dat` | Time (years) | Mass flow (kg/s; production is negative) |
| `<Well>_transient_enthalpy.dat`, `_injection_enthalpy.dat` | Time (years) | Enthalpy (kJ/kg) |

- **Well names.** `<Well>` is the `AUT2_Name` in title case (`Well2`). The model's 5-character name for WELL-10 is `WEL10`, so its files are named `Wel10_*`.
- **Gauge elevation.** For transient pressure, the gauge elevation is in the filename, for example `Well2_-750.0_transient_pressure.dat`.

### 4.3 Groups: `Data/Groups_data.json`

```json
{"Areas": {"<Project>": {"All": ["WELL-1", "WELL-2", "..."]}}}
```

### 4.4 Map overlays: `Data/Mapping/Mapping.json`

This is one object per feature, keyed by name (lake, rivers, faults):

```jsonc
{"Fault 1": {"Data": [[x, y], [x, y], ...], "type": "...", "c": "...", "ls": "...", "lw": 1, "legend": "..."}}
```

`Data` is a polyline (or polygon) of `[x, y]` points, in the same coordinate system as the well tracks. The other keys are matplotlib styling.

### 4.5 Model geometry: `Models/<model>/g<model>.dat` (AUTOUGH2 / MULgraph format)

| Section | Format | Used for |
|---|---|---|
| `VERTICES` | `<3-char id><x><y>`, fixed width | Column polygons |
| `GRID` | `<column id><number of vertices>`, then the vertex ids | Column polygons |
| `LAYERS` | `<index> <bottom elevation> <centre elevation>` (mRL); layer 0 is the top reference | Layer intervals and centres |
| `SURFA` | `<column id> <surface elevation>` | Model surface, compared with the wellhead |
| `WELLS` | `<AUT2_Name> x y z`, one line per track point | Cross-check against `Well_Track` |

**Block names** are `<column id><layer index>`, for example `arz38`.

`g<model>_wells.json` lists each well's block intersections as `{"block": "...", "interval": [top, bottom]}`. Use it to cross-check the blocks the script works out from the track.

## 5. Outputs

`CLAUDE.md` §5 defines the standard output names, so that runs can be compared:

| File | Format |
|---|---|
| `NS_data_interpretation.md` | Report: conclusions, data checks, T and P structure, conceptual model, calibration set-up, data issues |
| `NS_calibration_targets.csv` | One row per survey point. Columns: `well, x, y, elevation_mRL, T_degC, P_bar, Tsat_degC, boiling_margin_degC, column, layer, block, in_feedzone, behind_casing, above_model_surface, ns_confidence, weight, weight_reason` |
| `feedzones.csv` | One row per feed zone: survey T, P and boiling margin, model block, quality |
| `T_P_at_layer_centres.csv` | T and P for each well at each model layer centre (interpolated, never extrapolated) |
| `well_summary.csv` | One row per well: Tmax and its elevation, inversion, gradients, isotherm depths, boiling margin, water level, regime |
| `pressure_consistency_check.csv` | Survey pressure vs the first transient pressure at the gauge elevation |
| `hydrostatic_check.csv` | Tests of column density over long single-phase intervals |
| `well_geometry_check.csv` | Wellhead vs model surface, deviation, distances to faults and rivers |
| `figures/*.png` | Map, profiles, maps at fixed elevations, cross-sections, P–T state plot |
| `run_console_log.txt`, `process_log.md` | Run log, and the prompt plus the steps taken |

CSV conventions:
- comma-separated, with a header row;
- SI-style units are written into the column names (`_mRL`, `_degC`, `_bar`);
- blank cells mean no data, never zero.

## 6. Known caveats in the Teaching_Model_3D dataset

`CLAUDE.md` §3 has the full list. The main ones:

- **Survey timing.** All downhole surveys are dated 29/07/2015, after production and injection started on 01/01/2014. Whether each survey is natural state has to be tested per well.
- **Model surface.** The model surface (`SURFA`) is up to about 100 m above some wellheads.
- **WELL-4.** Its casing record extends below its track.
- **Late survey starts.** The surveys for WELL-4 and WELL-10 start below the wellhead.
- **WELL-10 feed zone.** It is tagged `Synthetic`.
- **Plotting files.** The JSON files in `Data/Plotting/` contain absolute paths from another machine. Do not use them.
