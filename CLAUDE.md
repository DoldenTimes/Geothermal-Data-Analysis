# CLAUDE.md: natural-state data interpretation for geothermal model calibration

## 1. Ground rules

- **Ignore `README.md`.**
- **`Baseline/` is read-only.** Read all input data from it. Never edit, move, rename or write anything inside it.
- **Ignore `Logs/`.**
- **Scripts go in `Scripts/`.** Data outputs, figures and reports go in `Output/`, with figures in `Output/figures/`.
- **Scope is the natural state only.** Production, injection and future scenarios are out of scope.
- **Audience.** Write every output for a reservoir engineer, or another agent, who will calibrate the natural-state model. Each finding should end in something they can act on: a target, a weight, a structural feature the model needs, or a data problem to resolve.
- **Process log.** Write `Output/process_log.md`. It must contain:
  - the user's prompt, verbatim;
  - a summary of this file's rules as they were applied;
  - every step taken, including errors found and how they were fixed;
  - any departure from these instructions, with the reason.
- **Previous runs.** Earlier results may be deleted from the working tree but still exist in git history. Do not restore or reuse them as inputs. Start from `Baseline/`.

## 2. What counts as data and what counts as reference

| Use as data (factual) | Use for geometry only | Reference only, never evidence |
|---|---|---|
| `Data/Data_json/WELL-*_data.json`: tracks, casing/liner design, feed zones, downhole T and P | `Models/TE14515/gTE14515.dat`: the `VERTICES`, `GRID`, `LAYERS`, `SURFA` and `WELLS` sections only | `gTE14515_dict_rocktype.json` (permeabilities, rock types) |
| `Data/Data_dat/*.dat`: the same T/P as plain columns. Use these only to cross-check the JSON | `Models/TE14515/gTE14515_wells.json`: well-to-block intersections | `gTE14515_boundary_conditions.json` (heat flux, upflows, rain, lake) |
| `Data/Mapping/Mapping.json`: lake, river and fault polylines | | `NaturalState/TE14515_NS_000.dat` and `.incon` (model results and initial conditions) |
| `Data/Groups_data.json`: the list of wells | | `current*Setup.json`, `Parameters/`, `Data/Plotting/` (the plotting files hold hard-coded paths from another machine) |

- **Model properties are not evidence.** Never present a property of the reference model, such as permeability, heat flux, upflow location or rock conductivity, as if the data showed it.
- **Describing the starting model is allowed.** You may describe the reference model so the engineer can see how far the data are from it. Label that content clearly as "reference model".
- **State any assumption you need.** If you need a parameter that is not in the data, for example rock thermal conductivity to turn a gradient into a heat flux, state the assumed value and where it came from.


## 3. Required analysis steps

### 3.1 Data integrity and timing
1. **Cross-check sources.** Compare the JSON with the dat files for every T and P point. Check that every survey point and feed zone lies within the well track and within the slotted or open interval.
2. **Test whether the surveys represent the natural state.**
   - Compare each survey pressure, interpolated to the transient-gauge elevation, with the **first** transient pressure value (WELL-2, 3 and 9 start in 2014; WELL-10's series starts in 2020).
   - Note how the injector profiles might be affected by injection.
   - Give a per-well verdict: natural state, probably natural state, or possibly disturbed. Give each verdict a confidence level, and carry it through to the target weights.

### 3.2 Wellbore physics: read the profiles correctly
- **Pressure.** Measured pressure represents the formation only at feed zones. Elsewhere it follows the fluid column inside the well, which pivots at the dominant feed zone. Use pressure at feed zones, or the pivot point, as the formation pressure target. Do not treat the shallow wellbore column, for example about 1 bar at the wellhead, as reservoir pressure.
- **Temperature.** Between feed zones in open hole, internal flow in the well can distort temperature, for example by producing isothermal sections. Behind casing, temperatures are usually closer to the conductive formation temperature. Before interpreting any feature as a reservoir feature, consider whether the well itself could have caused it.
- **Inversions and isothermal sections.** Interpret temperature inversions, isothermal sections and boiling-point-for-depth sections physically: lateral outflow, cold downflow, or two-phase conditions. Locate each one against the casing shoe and the feed zones.

### 3.3 Thermodynamics
- **Use IAPWS-IF97.** Use Region 4 for saturation temperature and pressure, and Region 1 for liquid density. **Do not use the Kell polynomial, or any fit only valid below 150 °C,** at reservoir temperatures.
- **Validate your implementation.** Check it against published reference values, for example Tsat(1.01325 bar) = 99.97 °C and Tsat(100 bar) = 311.0 °C. Make the script fail if any error exceeds 0.1 K.
- **Hydrostatic checks.**
  - Only test intervals that are single-phase liquid. Do not apply a liquid-density test where T is within a few °C of Tsat(P).
  - Use intervals hundreds of metres long. Over 10–20 m, the 0.1 bar rounding causes errors of tens of kg/m³.
  - Separate any field-wide gradient excess (seen even in cold wells) from the excess in individual wells. Report the field-wide excess as a possible datum or convention issue, not as physical overpressure.
- **Boiling margin.** Report Tsat(P) − T at every point, and flag every interval where it is below about 5 °C.

### 3.4 Interpretation for calibration
Produce at least the following:
- **Per-well characteristics:**
  - maximum temperature and its elevation;
  - temperature at total depth;
  - size and depth range of any inversion;
  - shallow (cap) and deep temperature gradients;
  - elevations of the 100 °C and 200 °C isotherms;
  - minimum boiling margin;
  - estimated water level or pivot pressure;
  - regime: upflow, outflow, conductive, recharge/downflow, or two-phase.
- **Lateral structure.** T and P at common elevations, and at model layer centres. Use interpolation only and never extrapolate. Show this as maps at several elevations, plus at least one cross-section through the hottest wells.
- **Conceptual model, built from the data alone.** Show:
  - where the upflow is;
  - the paths, depths and directions of outflow;
  - recharge and downflow areas;
  - the extent and thickness of the cap;
  - any two-phase zone;
  - how deep pressures are connected across the field.
  Relate each feature to the faults and surface water in `Mapping.json`. For every claim, state which observation supports it.
- **Feed zones.** Give survey T, P and boiling margin at each feed zone, its model block (from the track position), and its quality tag.

## 4. Standard deliverables

Use these exact names, so that runs can be compared with one another. Additional files are welcome.

| File | Contents |
|---|---|
| `Scripts/natural_state_analysis.py` | One script that regenerates every output when run from the repository root with `python Scripts/natural_state_analysis.py` |
| `Output/NS_data_interpretation.md` | The main report (structure in §6) |
| `Output/NS_calibration_targets.csv` | One row per survey point. Columns: `well, x, y, elevation_mRL, T_degC, P_bar, Tsat_degC, boiling_margin_degC, column, layer, block, in_feedzone, behind_casing, above_model_surface, ns_confidence, weight, weight_reason` |
| `Output/feedzones.csv` | Feed-zone targets with model blocks |
| `Output/T_P_at_layer_centres.csv` | T and P for each well at each model layer centre. No extrapolation |
| `Output/well_summary.csv` | The per-well characteristics listed in §4.4 |
| `Output/pressure_consistency_check.csv` | Survey pressure vs first transient pressure (§4.1) |
| `Output/hydrostatic_check.csv` | Column density tests (§4.3) |
| `Output/well_geometry_check.csv` | Wellhead vs model surface, deviation, distances to faults and rivers |
| `Output/figures/` | Well map with grid, faults and water bodies; T and P profiles with the Tsat curve, casing and feed zones; T and P maps by elevation; cross-sections; a P–T state plot against the saturation curve |
| `Output/run_console_log.txt` | Console output from the final run |
| `Output/process_log.md` | Prompt and steps (§1) |

**Weighting.** Choose a scheme and justify it in the report. The previous default, which can be adjusted:

| Location of the point | Weight |
|---|---|
| Feed zone | 1.0 |
| Open hole | 0.8 |
| Behind casing | 0.5 |
| Near the wellhead / top layer | 0.3 |
| Above the model surface | 0 |

Multiply the weight by the natural-state confidence from §4.1.

## 6. Report structure (`NS_data_interpretation.md`)

1. **Key conclusions for calibration.** Numbered. Each one names the observation that supports it.
2. **The data.** What exists and how it was checked. Includes the natural-state timing verdict for each well.
3. **Temperature structure.** Includes a well-by-well table.
4. **Pressure structure.** Includes a table at common elevations.
5. **Conceptual model.** A sketch or diagram, plus a description.
6. **Recommended calibration set-up.** Targets and weights. Where heat and mass need to enter and leave. The permeability structure the data imply. Boundary requirements.
7. **Data quality issues and open questions.** Rank them by how much they affect calibration.
8. **Output file index.**

## 7. Quality bar

- **Every number in the report comes from an output file** written by the script. Check the report against those tables before finishing.
- **Look at every figure** after it is generated, to confirm it renders and is readable.
- **Use plain units consistently** throughout: mRL, bar, °C, kJ/kg, kg/s.
- **Label the strength of each claim.** Distinguish observations ("WELL-3 reaches 261 °C at −450 mRL") from interpretations ("lateral outflow").
- **Never present an anomaly as fact before checking it.** If something looks physically impossible, such as a gradient steeper than a column of water, first check for a datum, unit or wellbore explanation, and report what you found.
