# Process log

## 1. The user's prompt (verbatim)

```
Analyse the data
```

## 2. Rules of `CLAUDE.md` as they were applied

- `README.md` and `Logs/` were not read. `Baseline/` was only read (no write, move or rename). Earlier results in git history and on other branches were not restored or used; the working tree had no `Output/` results and everything was generated again from `Baseline/`.
- Scripts are in `Scripts/` (`natural_state_analysis.py` written by me; `check_deliverables.py` provided and not edited). Tables and the report are in `Output/`, figures in `Output/figures/`. `Output/process_log.md` is the only hand-written output.
- Scope: natural state only. Production and injection records (the flow, enthalpy and transient-pressure series inside the well JSON files) were used only for the timing test of section 4.1 and to derive the injectate temperature. The production spreadsheets, forecast set-ups and the field-total file `Area_data.json` are `out_of_scope` and were not read.
- Inventory first (section 2.1): `Output/input_inventory.csv` is written by the script, every file under `Baseline/` is covered (22 rows). The roles that could not be filled: none; `Data/Make_up_json` and two `Reports` folders are empty directories (listed as `out_of_scope`). `Plotting` and `.auxjson` hold hard-coded paths from other machines (listed as `reference_other`).
- Roles (section 2.2): reference model files (rock types, boundary conditions, incons, set-up, plotting, auxjson) are never evidence. The report describes the reference model only in a paragraph labelled "[reference model]". The plain-column files were used only to cross-check the JSON files. The grid file was used for geometry only, the well-to-block file as the authority for blocks.
- Assumptions: the only parameter that is not in the data is the rock thermal conductivity (2.0 W/m/K) used to convert gradients to heat flux; it is labelled as an assumption in the report and listed in `analysis_parameters.csv`.
- Provenance: every Source, Comment, StableFlag and quality field was read and reported (`provenance_labels.csv`, report section 2c).
- Wellbore physics (4.2), thermodynamics (4.3: IF97 Regions 4 and 1, validation that stops the script above 0.1 K, windows at least 300 m, bands cased/open, cold-well component, boiling margins), interpretation (4.4) and deliverables (section 5) are implemented in the script; the weight tables and the quality-tag mapping follow section 5 without change.
- Report structure follows section 6 (eight sections, claims labelled observation / interpretation / assumption / reference model, plain units).

## 3. Steps taken

1. Read `CLAUDE.md` and `Scripts/check_deliverables.py` (the check also fixes some conventions: the report-number match, a deterministic re-run, the inventory-driven geometry paths).
2. Listed `Baseline/` (75 files) and inspected every kind of file: well JSON (keys, track, design, feed zones, surveys, records), plain-column files, grid file (sections `VERTICES`, `GRID`, `CONNECTIONS`, `LAYERS`, `SURFA`, `WELLS`), well-to-block file, mapping, reference files (only to classify them).
3. Confirmed or rejected each statement of section 3 against the files (results in section 5 below).
4. Wrote the script in sections (IF97, inventory, readers, per-well analysis, field-wide analysis, tables, figures, report) and tested each part against the data before going on.
5. Ran the script, inspected the tables, looked at every figure, fixed layout problems (fault labels over wells, legends over data, clipped labels), ran the acceptance check, fixed what it reported, reviewed the report against the tables, re-ran.

### Errors found and how they were fixed

- *Shell tool:* two multi-line heredocs were rejected by the tool; the files were written with the file tool instead (no effect on results).
- *Crash in the reader:* the helper that sorts JSON files by the number in their name failed on `Mapping.json` (no number). Replaced by a natural sort key.
- *Plain-column time reference:* the first rule (first JSON date minus first time, rounded to a day) gave 31/12/2013 for ten files because the series use 1/12-year steps. Now the candidates within a few days are scored and 1 January is preferred when its error is within 0.003 yr of the best; all 18 time-series files give 01/01/2014.
- *Deep gradient of W2:* the interval below the deepest feed zone holds two survey points; a three-point minimum left it blank. Two points now give the (exact) least-squares slope; a note in the report says so.
- *Injector rule:* my first rule compared the profile maximum with the injectate temperature; W8's maximum (93.1 degC) is at its bottom and not at the feed zone, so W8 was given low confidence. The rule now compares the feed-zone temperatures (63.9 degC for W8) and looks for an isothermal section near the injectate temperature; W8 becomes medium, W5 stays low. This change was made after seeing the result, for the physical reason just given.
- *Acceptance check quirk:* the check reads "8,5861.0" in a CSV as one number (comma as thousands separator). The report-number test then could not find some values. Fixed in the script by writing integer columns last in `hydrostatic_summary.csv`, writing `weights_by_well.csv`, registering derived numbers in `report_derived_numbers.csv` and avoiding decimal section numbers in the report headings (4a, 4b ...). The check itself was not changed.
- *Pandas deprecation warning* (`reset_index(inplace=True)`): removed.

## 4. Independent review: conclusions re-read against the tables, tried to disprove

| Item | Attempt to disprove | Result / change |
|---|---|---|
| Surveys are natural state | Is the 0.1 bar match of W2 and W3 trivial? The gauge series start on the first production day with a flow already running. The survey is 9.1 to 14.7 bar above the gauge interpolated to the labelled date. | Kept ("natural state" high for W2, W3) as the instructions require; the date label is a data question (report 2d, issue 2). |
| W9 | Survey 2.05 bar above its first gauge value: a match within precision fails. | Verdict "probably natural state", medium, with the reason; not "high". |
| Feed-zone pressure trend | The trend is built from only W2 and W3. | Reported with the number of trend wells; limit 3 bar fixed before the results were seen; W1 is at 2.92 bar and is flagged as borderline (issue 10). |
| Boiling at W10 | Is the 2.6 degC margin real? Measured pressures are above the IF97 column, which raises Tsat. | The first draft told the user not to force boiling onto W2 and W3, which is the wrong direction (too-high pressures hide boiling). Corrected: margins of W2 and W3 are upper bounds until the pressure convention is known (conclusion 6). |
| Hot core and upflow | Could the hot wells be an outflow only? | Their temperatures are the highest and shallowest, the 200 degC isotherm is shallowest there and W2 is still heating at its bottom; the report states the alternatives (upflow on Fault 1 or between faults) and the deciding observations. |
| Outflow directions | Tmax decreases with distance, but a decrease with distance also fits conduction away from a source. | Labelled interpretation; the alternative (cold water from below) and the deciding observation are in section 5. W1 has no inversion and is not called outflow. |
| W8 as cold downflow | Its shallow gradient (84.6 degC/km) is conductive-like. | Wording changed: cooling only below -550 mRL, "lateral cool inflow as well as downflow"; W4 and W8 are described separately. |
| Hydrostatic excess | Could a convention explain 1.55? | Datum, gravity, two-phase and rounding each checked in section 4c; the hot cased excess is not explained; the report says so and does not call it overpressure. |
| Tsat above the critical pressure | List of wells exceeding 220.64 bar was wrong in the draft (W2 and W6 missing). | Now computed from the data. |
| Feed-zone elevation ranges (6c) | Draft listed W6 in both shallow and deep groups by hand. | Now computed from the feed-zone table. |
| Heat flux | The draft said the margins have weaker heat input than the reference value, but W8's assumed-k flux (169 mW/m2) is above it. | Corrected: the margin values bracket the reference value. |
| Section 4c claim about gradients | "Largest in wells with the steepest gradients" | Checked: the cased ratio rises with the shallow gradient for eight wells with a cased window (W8 and W1 swap order); stated as "broadly". |

## 5. Checks of the statements of section 3 of `CLAUDE.md`

- Grid file: MULgraph, `LAYERS` lines are `index bottom centre`, layer 0 gives the top surface (1000 mRL), layers tile the grid (checked: no gaps). **Departure:** the vertex lines have the two coordinates run together (`aab2812872.706327567.88`), so whitespace parsing is impossible there; they were read as a 3-character name followed by two numbers with two decimals. The `GRID` header lines are `name` plus a count (`aab0 4`), also read as given.
- Well-to-block file: blocks along the track, feed-zone blocks in order, interpolation weights: confirmed. Own point-in-polygon and layer search agrees with it at all 266 survey points.
- Well JSON: list entries are strings, read with a literal evaluator; feed-zone top and bottom are labelled "m" but are mRL (they match the track in mRL); surveys keyed by date then elevation string: confirmed.
- Plain-column files: survey files are `value elevation`, transient files are `years value` (reference date 01/01/2014 in all files, including W10), production flows have a negative sign in the plain files, injection flows are positive. Irregular names: the five `Wel10_*` files.
- Mapping: the river entries are ragged lists of segments, faults are single polylines. 'Lake Tarawera' is a single closed outline (not ragged), which the reader handles.

## 6. Departures from the instructions (and why)

1. Vertex lines of the grid file are parsed with a fixed pattern, not by whitespace (impossible: no separator).
2. Hydrostatic windows never straddle the casing shoe. The instructions say "split at the casing shoe where the 300 m minimum allows"; where it does not allow, I do not test the short side (4 pieces: W2 cased 290 m, W4 open 250 m, W10 cased 160 m and the empty W10 open piece), so every window has one of the two band labels. They are listed in `run_console_log.txt`.
3. Deep-gradient fits with two points are accepted (the least-squares line through two points); the number of points is in `well_summary.csv` (`deep_fit_n`).
4. "Top block of its column" is read as the model layer that contains the model surface of that column. In this dataset no survey point lies in that block more than 30 m below its wellhead (checked), so the two conditions give the same 0.3 weights; W4 and W10 surveys start 110 m and 50 m below the wellhead and have no such points.
5. The thresholds of the natural-state tests that the instructions leave open were set by me before the results were seen: match within 0.15 bar (reading precision 0.1 bar), feed-zone trend limit 3 bar, injector rule as described above. They are in `analysis_parameters.csv`.
6. Isothermal section = at least 300 m long with a temperature range of at most 5 degC (the instructions give only the length).
7. Regime classes follow a rule that I wrote (report section 3); the instructions give the class names only.
8. Additional outputs beyond section 5 (all generated by the script): `dat_json_crosscheck.csv`, `integrity_checks.csv`, `provenance_labels.csv`, `pressure_trend_check.csv`, `hydrostatic_summary.csv`, `isothermal_sections.csv`, `boiling_margin_flags.csv`, `TP_at_common_elevations.csv`, `temperature_at_common_elevations.csv`, `pressure_at_common_elevations.csv`, `weights_by_well.csv`, `if97_validation.csv`, `analysis_parameters.csv`, `reference_model_summary.csv`, `report_derived_numbers.csv`, a pressure section figure.
9. The report cites a few round qualitative numbers by hand where the data support them (for example "within a kilometre" for river distances); all of them were checked against the tables.

## 7. Acceptance check (`python Scripts/check_deliverables.py --rerun`)

Final run: `python Scripts/natural_state_analysis.py` followed by `python Scripts/check_deliverables.py --rerun` (the re-run regenerates every table and the report in a temporary copy and compares them; the check exited with code 0). Full output (`Output/acceptance_check.txt`):

```
PASS  I1   22 inventory rows, all roles valid and paths present
PASS  I1   every Baseline file is covered by the inventory
PASS  F1   all 12 required files present
PASS  F1   all required figure types present
PASS  S1   input_inventory.csv: columns OK
PASS  S1   NS_calibration_targets.csv: columns OK
PASS  S1   feedzones.csv: columns OK
PASS  S1   T_P_at_layer_centres.csv: columns OK
PASS  S1   well_summary.csv: columns OK
PASS  S1   pressure_consistency_check.csv: columns OK
PASS  S1   hydrostatic_check.csv: columns OK
PASS  S1   well_geometry_check.csv: columns OK
PASS  S2   one pressure-consistency row per well
PASS  S2   hydrostatic band labels OK
PASS  T1   Tsat agrees with an independent IF97 Region 4 within 0.1 K at every point
PASS  T2   boiling margins consistent
PASS  W1   weights valid
PASS  W2   ns_confidence takes 3 values across wells
PASS  R1   0/892 report numbers (0%) not found in any output table; first: []
PASS  R2   analysis script writes the report
PASS  P1   process log has prompt, departures and acceptance check
PASS  G1   all 266 survey-point blocks match the well-to-block file
PASS  G2   all 18 feed-zone blocks match
PASS  G3   all 30 layer centres match the grid file
PASS  G3   no layer-centre row lies outside its well's survey range
PASS  R3   re-run reproduces every CSV and the report

0 FAIL, 0 WARN, 26 PASS
```

**WARN items: none.** Notes on checks that were close: W2 passes because `ns_confidence` takes three values (high, medium, low), so the confidence factor does change relative weights (explained in report section 2d). R1 passes with 0 of 892 report numbers unmatched; before the fixes listed in section 3 it showed 44 unmatched numbers (a WARN) caused by the way the check reads comma-separated integers, and 0 FAIL throughout the later runs. The earlier FAIL on the missing process log and the crash of the process-log check were only because this file had not been written yet.

