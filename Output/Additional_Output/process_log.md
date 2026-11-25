# Process log: natural-state data interpretation

## 1. User prompt (verbatim)

> Analyse the data

## 2. CLAUDE.md rules as applied

- `README.md` and `Logs/` ignored. `Baseline/` read only; nothing written there. No earlier results restored from git history; everything regenerated from `Baseline/`.
- Scripts in `Scripts/`; outputs only in `Output/Calibration_Input/` (report, four target tables, figure index, standard figures F1-F7) and `Output/Additional_Output/` (everything else). Nothing loose in `Output/`.
- Scope: natural state only. Production/injection records were read only for the timing test (§4.1): first-flow dates, flow periods for F7, injectate enthalpy. Production-model spreadsheets were inventoried as `out_of_scope` and not read.
- Inventory first (§2.1), with roles from §2.2. Reference-model files were used only for a short, labelled `[reference model]` paragraph (boundary-condition defaults); no reference property is presented as evidence.
- Defaults in standard columns and figures: §4.3 window algorithm, §4.4 definitions, §4.1 rubric, §5.5 weights, §8 house style. Extra items are added as extra columns or files (e.g. `cap_base_mRL`, `regime_note`, `near_boiling`, `hydrostatic_summary.csv`, `parameters_used.csv`), never in place of the default.
- IF97 Region 1 and Region 4 implemented and checked against the IF97 verification tables; the script stops if any check fails.
- Every report number is generated from the output tables by the script; the report is written by the script.

## 3. Steps taken

1. **Inventory.** Listed all 75 files under `Baseline/Teaching_Model_3D`. Assigned roles by rule (20 inventory rows, directory rows where possible). Empty directories found: `Data/Make_up_json`, `NaturalState/Reports`, `Production/Reports`, `Parameters/FutureScenarios`. Ambiguous roles recorded here:
   - `Data/Area_data.json` (field-total steam flow and power) has no exact role; assigned `well_data_json` (flow records) and read only for the first-field-production date.
   - `NaturalState/TE14515_NS_000.dat` holds rocks, generators and parameters together; assigned `reference_other` (it could also be `reference_rock`/`reference_bc`). Not used as evidence.
   - No role was left unfilled.
2. **Format checks (§3).** Confirmed: well JSON list entries are Python-literal strings (parsed with `ast.literal_eval`); feed-zone top/bottom labelled "m" are elevations; surveys keyed by date then elevation string; plain-column survey files are `value elevation`; transient files are years since a reference date, established by matching as 2014-01-01 (365.25-day years; monthly JSON dates differ from the dat times by up to about 2 days, so matching uses a 3-day tolerance); production mass flow in `.dat` files has the opposite sign to the JSON; files for WELL-10 use the irregular prefix `Wel10`; mapping rivers are ragged multi-segment lists and the lake is a filled polygon.
   - **Format that does not hold:** the MULgraph VERTICES section has no whitespace between vertex name and coordinates (e.g. `aab2812872.706327567.88`), so it cannot be parsed by whitespace. It is parsed with a regex (name letters, then two numbers with exactly two decimals); every vertex line matches. The LAYERS, SURFA, GRID and WELLS sections are parsed by whitespace.
   - Every column has a surface entry, so the grid-top fallback was never needed.
3. **Cross-checks.** 38 dat files vs JSON: all values identical. Integrity checks (106): one issue — WELL-4 liner bottom -3251 mRL vs track TD -750 mRL. All survey points and feed zones within tracks; all feed zones within slotted liners.
4. **Natural-state test.** Gauge precision 0.1 bar (one decimal). WELL-2 and WELL-3 match their first readings (2014-01-01, the same day as first field production) within 0.1 bar → natural state/high. The reading taken on the first-production date was treated as pre-production (it is the initial reading). WELL-9 misses by +2.05 bar → medium (feed zones on trend). WELL-10 gauge starts 2020 → medium. Injectate: 335.0 kJ/kg → 80.0 °C (IF97 Region 1 inverse at 1.01325 bar). WELL-5 (feed zones 147 °C > 80 °C) → low. WELL-4 and WELL-8: injector rule judged at feed-zone temperature (27 and 64 °C, colder than injectate, no isothermal section near 80 °C) → evidence against injection effects → medium via trend test. WELL-1: feed zone 3.32 bar below the trend, over the 3.00 bar tolerance → "possibly disturbed"/low by the rubric, although it has no flow records (reported as a data question).
5. **Analysis tables**, then figures F1-F7, then the report.
6. **Errors found and fixed:**
   - IF97 Region 1 verification failed at 500 K (v off by 25 %): coefficient n33 had been typed as 0.18228094581404e-20; correct value is 0.18228094581404e-23. The built-in verification caught it before any result was produced.
   - `DataFrame.itertuples` drops underscore-prefixed column names → AttributeError; fixed by renaming before iterating.
   - First vertex regex was greedy and split coordinates wrongly (point-in-polygon found no column for any wellhead). Fixed with the two-decimal pattern; afterwards point-in-polygon agrees with the well-to-block file for every wellhead.
   - Section parser stopped at the first WELLS line starting with four capitals (`WELL1`), so the WELLS section read empty; fixed by stopping only at known section headers or a blank line.
   - Dat reference-date rounding gave two candidates (2013-12-31 / 2014-01-01); fixed by taking the most frequent candidate across all files.
   - Figure layout: overlapping labels in profiles, section titles and F6 cap labels; moved.
   - Checker regex merges an integer column with the next number (`2,1.0458` read as `21.0458`); moved `n_windows` to the last column of `hydrostatic_summary.csv` so the report's medians are matched.

## 4. Geometry self-check (§7)

| Item | Source compared | Result |
|---|---|---|
| Model surface | `model_surface_mRL` vs SURFA value of the wellhead column (column of the top block in the well-to-block file) | All 10 match (G4 PASS). 8 wellheads lie below the model surface, by up to 105 m (WELL-1); WELL-4 and WELL-10 are at the surface. |
| Blocks | Survey-point and feed-zone blocks vs well-to-block file | Taken directly from the file; 266/266 and 18/18 match (G1, G2). Wellhead column by point-in-polygon equals the block-file column for all wells. |
| Layer centres | `T_P_at_layer_centres.csv` vs LAYERS third column | All 30 layers match; no row outside a well's survey range (G3). |
| Wellhead elevation | First track point vs block-file top interval vs grid WELLS section | Equal for all 10 wells (`well_geometry_check.csv`, `integrity_checks.csv`). |

## 5. Independent review (attempts to disprove conclusions)

- Conclusion 6 (surveys pre-date production): checked that a disturbed survey would read lower, not higher; WELL-2/3 differences are -0.10 bar while their gauges fell 8.2-11.5 bar by the second reading; WELL-10's survey is 8.0 bar above its 2020 reading. Replaced a hand-written "about 9 bar" with the computed range.
- WELL-9 miss: replaced a hand-estimated "about 25 m" with the computed elevation at which the survey reaches the first gauge pressure (-651.2 mRL, 23.8 m above the gauge).
- Conclusion 1: "blocks ajf, arz" were columns, not blocks — wording corrected. A two-tag claim `[obs, IF97]` violated the tag rule — fixed.
- Conclusion 8 (super-hydrostatic cased gradients): checked datum (cannot change gradient), convention (would act at all depths; open-hole windows of the same wells are near the cold component), two-phase (lighter), rounding (≤ 0.0070 in ratio). Kept as a data question, not presented as physical fact.
- Section 4 pressure table: large positive deviations at depth in WELL-8 could be read as overpressure; added the explanation that they are a cold-column density effect.
- Regime of WELL-10 ("upflow") checked against the rule; the well is only 610 m deep, so a `regime_note` states the alternative (near-source outflow) without changing `regime`.

## 6. Departures from the instructions

- VERTICES section parsed by regex, not whitespace (no whitespace exists; §3). Documented above.
- Hydrostatic windows: eligible intervals are split at the production casing shoe before the window algorithm runs, so no window straddles the shoe (needed for the cased/open band labels). The interval between the last cased point and the first open-hole point is therefore not tested. A cased part shorter than 300 m (WELL-10) forms a single short window.
- F6 section uses line B (an F4 line from the hottest well WELL-10 east through WELL-3, WELL-5, WELL-1) rather than line A, because line A crosses only WELL-8, WELL-10 and WELL-2 and shows no outflow. Line A is the default F4 line and is drawn.
- Elevation axis fixed at 400 to -2900 mRL as required; the model surface west of the well field (up to about 560 mRL) is clipped in the sections.
- Injector test (§4.1) evaluated at feed-zone temperatures (where injectate would enter); WELL-8's deepest survey temperature (93 °C) exceeds the injectate, but its feed zone (64 °C) does not.
- Optional reference-model companion figure (§8.3) not produced.

## 7. Acceptance check (final run, `python Scripts/check_deliverables.py --rerun`)

```
PASS  I1   20 inventory rows, all roles valid and paths present
PASS  I1   every Baseline file is covered by the inventory
PASS  F1   all 14 required files present in their folders
PASS  F1   all standard figure types present
PASS  F2   figure index lists 23 files covering ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7']
PASS  F2   every Calibration_Input figure is indexed
PASS  F3   no standard file is misplaced
PASS  F3   Calibration_Input holds only standard files
PASS  S1   input_inventory.csv: columns OK
PASS  S1   NS_calibration_targets.csv: columns OK
PASS  S1   feedzones.csv: columns OK
PASS  S1   T_P_at_layer_centres.csv: columns OK
PASS  S1   well_summary.csv: columns OK
PASS  S1   pressure_consistency_check.csv: columns OK
PASS  S1   hydrostatic_check.csv: columns OK
PASS  S1   well_geometry_check.csv: columns OK
PASS  S1   T_P_at_elevation_bands.csv: columns OK
PASS  S1   figures_index.csv: columns OK
PASS  S2   one pressure-consistency row per well
PASS  S2   hydrostatic band labels OK
PASS  S3   UTF-8, integer feed zones, ISO dates
PASS  S3   LF line endings
PASS  T1   Tsat agrees with an independent IF97 Region 4 within 0.1 K at every point
PASS  T2   boiling margins consistent
PASS  D1   all 160 well_summary values match the §4.4 definitions
PASS  W1   weights valid
PASS  W2   ns_confidence takes 3 values across wells
PASS  W3   weights follow §5.5; verdict and confidence agree across files
PASS  R1   1/781 report numbers (0%) not found in any output table; first: ['2.2']
PASS  R2   analysis script writes the report
PASS  R4   no unfilled template text in the report
PASS  P1   process log has prompt, departures and acceptance check
PASS  G1   all 266 survey-point blocks match the well-to-block file
PASS  G2   all 18 feed-zone blocks match
PASS  G3   all 30 layer centres match the grid file
PASS  G3   no layer-centre row lies outside its well's survey range
PASS  G4   model_surface_mRL matches the surface section for all 10 wells
PASS  R3   re-run reproduces every CSV and the report

0 FAIL, 0 WARN, 38 PASS
```

No WARN to explain. The single unmatched report number in R1 ('2.2') is the section reference "2.2", not a value.
