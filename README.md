# Geothermal Data Analysis: natural-state interpretation

This repository analyses geothermal well data with AI agents (Claude Code). It turns downhole surveys and well geometry into **natural-state calibration targets** and a **conceptual model**, for a reservoir engineer or agent calibrating a numerical reservoir model.

- **`CLAUDE.md`** holds the binding rules for agents. It is dataset-agnostic: it names *roles* of input, not files.
- **This README** is for the people who prepare data. It sets out the **required file structure and formats** for `Baseline/`. Agents are told to ignore it, so the structure must be right in the data itself.

**Why the structure matters.**
- Every run starts by listing `Baseline/` and assigning each file a role (`CLAUDE.md` §2.1). When the data follow this standard, that step is mechanical and the agent's effort goes into interpretation.
- When they don't, each run has to reverse-engineer the layout, and runs diverge. In past runs, most of the differences between models came from parsing, not from interpretation.

---

## 1. Repository layout

```
<repository root>/
├── CLAUDE.md                     Rules for agents
├── README.md                     This file: data standard
├── Baseline/                     Input data. READ-ONLY during analysis
│   └── <Project>/                One project (field + model) per folder; see §2
├── Scripts/
│   ├── check_deliverables.py     Acceptance check (provided; agents must not edit it)
│   └── natural_state_analysis.py Written by the agent in each run
├── Output/                       Everything a run generates
│   └── figures/
└── Logs/                         Earlier reports and notes. Ignored by agents
```

- **One project at a time.** Keep a single `<Project>/` under `Baseline/` per analysis. If there are several, the inventory step has to decide which one is meant.
- **Earlier runs** are kept in git history and on branches. A new run never reuses them.

## 2. Required `Baseline/` structure

```
Baseline/<Project>/
├── Data/
│   ├── Groups_data.json                    Well list                      → well_groups
│   ├── Data_json/
│   │   └── <Well>_data.json                One per well                   → well_data_json
│   ├── Data_dat/
│   │   ├── <Well>_downhole_temperature.dat                                → well_data_dat
│   │   ├── <Well>_downhole_pressure.dat
│   │   ├── <Well>_<elev>_transient_pressure.dat      (optional)
│   │   ├── <Well>_production_mass_flow.dat           (optional)
│   │   ├── <Well>_transient_enthalpy.dat             (optional)
│   │   ├── <Well>_injection_mass_flow.dat            (optional)
│   │   └── <Well>_injection_enthalpy.dat             (optional)
│   └── Mapping/
│       └── Mapping.json                    Faults, rivers, lakes          → mapping
├── Models/
│   └── <model>/
│       ├── g<model>.dat                    Grid (MULgraph format)         → geometry_dat
│       ├── g<model>_wells.json             Well-to-block intersections    → well_blocks_json
│       ├── g<model>_dict_rocktype.json     Rock types, permeabilities     → reference_rock
│       ├── g<model>_boundary_conditions.json                              → reference_bc
│       ├── NaturalState/                   Model input, results, incon    → reference_results
│       ├── Production/                                                    → out_of_scope
│       └── FutureScenarios/                                               → out_of_scope
└── Parameters/                             Production histories, scenarios → out_of_scope
```

The arrows give the role each item should receive in the agent's `Output/input_inventory.csv` (`CLAUDE.md` §2.2).

**Naming rules**

| Item | Rule |
|---|---|
| `<Project>` | Any folder name without spaces |
| `<Well>` | **One name per well, used everywhere:** in `Groups_data.json`, in the JSON `Well_Name`, in the JSON and dat file names, and as the key in `g<model>_wells.json`. Do not abbreviate or re-case it in file names |
| `<elev>` | Gauge elevation in mRL, with one decimal place, for example `-750.0` |
| `<model>` | The model name. The grid file and well-to-block file must share it |

**Keep `Baseline/` clean.** Do not add plotting files, scratch copies, OS files (`.DS_Store`) or alternative versions of the same data. Every extra file must be classified in each run. Data that is not part of this standard belongs in a separately named sub-folder (for example `Data/Extra/`), with a `README.txt` saying what it is.

## 3. Required file formats

These formats are what `CLAUDE.md` §3 tells agents to expect. **"Required"** means the analysis or the acceptance check depends on it. **"Tolerated"** means agents are told to cope with it, but it should be fixed in new data.

### 3.1 Well file: `Data/Data_json/<Well>_data.json`

```jsonc
{
  "Well_Name": "<Well>",
  "AUT2_Name": "<name of the well in g<model>.dat>",
  "Feed_Zones": {
    "Units": {"Top": "mRL", "Bottom": "mRL", "Quality": "Text", "Source": "Text"},
    "Data": [[<top>, <bottom>, "<QUALITY>"], ...]
  },
  "Well_Track": {
    "Units": {"x": "m", "y": "m", "z": "mRL"},
    "Data": [[<x>, <y>, <z>], ...]                    // wellhead first, downwards
  },
  "Well_Design": {
    "Units": {"Top": "mRL", "Bottom": "mRL", "Diameter": "\"", "Slotted": "Boolean", "Spec": "Text"},
    "Source": "<where the design data came from>",
    "Data": [[<top>, <bottom>, "<diameter>", <slotted>, "<spec>"], ...]
  },
  "Downhole_Temperature": {
    "<dd/mm/yyyy>": {
      "Units": {"Date": "dd/mm/yyyy", "Elevation": "mRL", "Temperature": "degC"},
      "Source": "<logging tool / contractor>", "Comment": "<conditions: shut-in time, static/flowing>",
      "StableFlag": true,
      "Data": {"<elevation>": <value>, ...}
    }
  },
  "Downhole_Pressure": { "<dd/mm/yyyy>": { "...": "same structure; value in bar absolute" } },

  // Optional: used only for the natural-state timing test
  "Transient_Pressure": {
    "<gauge elevation mRL>": {"Units": {"Date": "dd/mm/yyyy", "Pressure": "bar"}, "Data": {"<dd/mm/yyyy>": <value>, ...}}
  },
  "Production_Mass_Flow": {"Units": {...}, "Data": {"<dd/mm/yyyy>": <kg/s>, ...}},
  "Transient_Enthalpy":   {"Units": {...}, "Data": {"<dd/mm/yyyy>": <kJ/kg>, ...}},
  "Injection_Mass_Flow":  {"Units": {...}, "Data": {"<dd/mm/yyyy>": <kg/s>, ...}},
  "Injection_Enthalpy":   {"Units": {...}, "Data": {"<dd/mm/yyyy>": <kJ/kg>, ...}}
}
```

| Rule | Status |
|---|---|
| **Elevations, not depths.** All elevations are mRL (metres relative to sea level, positive up), and must be labelled `mRL` | Required. A label of "m" on a value that is really mRL is tolerated |
| **Feed zone** `[top, bottom, quality]`, with top > bottom, lying inside the track and inside a `slotted = true` interval | Required |
| **Quality tags** are a graded set, `HIGH`, `MEDIUM`, `LOW`. Any other tag (for example for a generated value) must be explained in the well's `Source` | Required |
| **Well track** in the same projected coordinate system as the model grid and `Mapping.json`. Vertical wells need only two points; deviated wells list every survey station | Required |
| **Casing/liner rows** `[top, bottom, diameter, slotted, spec]`. `slotted = true` marks the open interval, and the production casing shoe must be identifiable | Required |
| **T and P surveys** use the same elevations and are listed per survey date. Pressure is in **bar absolute** | Required |
| **`StableFlag: true`** only on surveys taken after heat-up | Required |
| **Survey dates** are the actual logging date. A placeholder or file date is not acceptable: the natural-state test compares surveys with the gauge records | Required |
| **`Source` and `Comment`** say where the data came from, and whether they are measured, processed or synthetic | Required |
| Rows stored as **strings containing Python literals** (`"[-250, -270, 'LOW']"`) rather than JSON arrays | Tolerated. Prefer real JSON arrays |
| Dates as `dd/mm/yyyy` | Required |

### 3.2 Plain-column copies: `Data/Data_dat/*.dat`

These are two whitespace-separated numeric columns, with no header. They must equal the JSON values exactly; they are used only to cross-check them.

| File | Column 1 | Column 2 |
|---|---|---|
| `<Well>_downhole_temperature.dat` | Temperature (°C) | Elevation (mRL) |
| `<Well>_downhole_pressure.dat` | Pressure (bar abs.) | Elevation (mRL) |
| `<Well>_<elev>_transient_pressure.dat` | Time (years since the reference date) | Pressure (bar abs.) |
| `<Well>_production_mass_flow.dat`, `<Well>_injection_mass_flow.dat` | Time (years since the reference date) | Mass flow (kg/s) |
| `<Well>_transient_enthalpy.dat`, `<Well>_injection_enthalpy.dat` | Time (years since the reference date) | Enthalpy (kJ/kg) |

- **Reference date.** State the time-axis reference date in a `Data/Data_dat/README.txt`. Otherwise every run has to infer it by matching against the JSON.
- **Coverage.** Every well in `Groups_data.json` needs at least the temperature and pressure files.

### 3.3 Well list: `Data/Groups_data.json`

```json
{"Areas": {"<Project>": {"All": ["<Well>", "<Well>", "..."]}}}
```

Every listed well must have a `<Well>_data.json`, and there must be no JSON file for a well that is not listed.

### 3.4 Map features: `Data/Mapping/Mapping.json`

```jsonc
{"<Feature name>": {"Data": [[x, y], ...], "type": "...", "legend": "...", "...": "styling"}}
```

- **Coordinates.** Points are in the same coordinate system as the well tracks.
- **Single features.** Faults are single polylines. A river or shoreline may be a list of segments (`[[[x, y], ...], [[x, y], ...]]`), but each feature must be one entry.
- **Feature names.** Use names that say what the feature is (`Fault <id>`, `River <name>`, `Lake <name>`). The agent reports distances per feature, using these names.

### 3.5 Model geometry: `Models/<model>/g<model>.dat` (MULgraph format)

| Section | Format | Used for |
|---|---|---|
| `VERTICES` | vertex id, x, y | Column polygons |
| `GRID` | column id, number of vertices, then the vertex ids | Column polygons |
| `LAYERS` | `index  bottom_elevation  centre_elevation` (mRL); layer 0 is the top surface | Layer intervals and centres |
| `SURFA` | column id, surface elevation | Model surface vs wellhead |
| `WELLS` | well name, then x y z per track point | Cross-check against the JSON track |

- **Block names** are the column id followed by the layer number.
- **Wellheads.** The model surface should be at or below each wellhead. Where it is above one, either fix `SURFA` or expect the agent to report it as a data issue.

### 3.6 Well-to-block file: `Models/<model>/g<model>_wells.json`

```jsonc
{"<Well>": {
   "Blocks":     [{"block": "<block>", "interval": [<top mRL>, <bottom mRL>]}, ...],
   "Feed_Zones": {"<feed-zone name>": {"<block>": {...}}, ...},     // in the same order as the well file
   "Interp_Weights": {...}
}}
```

**This file is the authority for block assignment.** The acceptance check compares every survey point and feed zone against it (checks G1 and G2). It must be regenerated whenever the grid or a well track changes.

### 3.7 Reference-model files

Rock types, boundary conditions, natural-state input and results, and set-up files are **reference only**. Agents may describe them, labelled "[reference model]", but never use them as evidence.
- Keep them under `Models/<model>/` so that they cannot be mistaken for field data.
- Set-up files with absolute paths from another machine are acceptable, but they must not be needed to run the analysis.

## 4. Preparing a new dataset: checklist

1. Create `Baseline/<Project>/` with exactly the folders in §2, and nothing else at top level.
2. Use one well name everywhere (§2 naming rules), including in the dat file names.
3. Store every elevation as mRL and label it `mRL`.
4. Give real survey dates. Fill `Source` and `Comment`, including whether the data are measured or synthetic.
5. Use `HIGH`/`MEDIUM`/`LOW` quality tags, or document any other tag.
6. Make sure the dat files equal the JSON, and record the time-axis reference date in `Data/Data_dat/README.txt`.
7. Regenerate `g<model>_wells.json` from the current grid and tracks.
8. Remove plotting files, OS files and duplicates, or move them into a documented sub-folder.
9. Commit `Baseline/` before the first run. Runs must never modify it.

## 5. Running an analysis

**Requirements:** Python 3 with numpy, matplotlib and pandas. No IF97 package is needed: the analysis implements IAPWS-IF97 itself, and the acceptance check has its own independent implementation.

An agent run, following `CLAUDE.md`:
1. Lists `Baseline/` and writes `Output/input_inventory.csv`.
2. Writes `Scripts/natural_state_analysis.py`, then runs it from the repository root:

   ```powershell
   python Scripts/natural_state_analysis.py
   ```

3. Runs the acceptance check, which must finish with no FAIL:

   ```powershell
   python Scripts/check_deliverables.py --rerun
   ```

4. Writes `Output/process_log.md`, including the acceptance-check output.

The required output files and their exact columns are defined in `CLAUDE.md` §5, and are not repeated here, so that there is one source of truth. The main report is `Output/NS_data_interpretation.md`.

## 6. Current dataset: departures from this standard

The dataset now in `Baseline/` predates this standard. The agent rules tolerate the following points, but each one costs effort in every run and should be fixed when the data are next revised.

| Departure | Standard (§) |
|---|---|
| One well's dat files use a shortened, re-cased name (`Wel10_*` instead of the JSON well name) | 2 |
| Feed-zone `Top`/`Bottom` are labelled "m" although they are mRL; the track `z` is labelled "masl" | 3.1 |
| JSON rows are stored as strings containing Python literals | 3.1 |
| Every survey carries the same date, which conflicts with the gauge records | 3.1 |
| One feed zone has the quality tag `Synthetic`; every `Well_Design.Source` reads "Synthetic data" without further explanation | 3.1 |
| One well's design lists a liner deeper than its track | 3.1 |
| The dat time-axis reference date is not documented | 3.2 |
| Several wellheads lie below the model surface | 3.5 |
| Extra, unstandardised items: `Data/Area_data.json`, `Data/Plotting/`, an empty `Data/Make_up_json/`, `Models/<model>/g<model>.auxjson`, `currentGeomSetup.json`, and a `.DS_Store` under `Parameters/` | 2 |
