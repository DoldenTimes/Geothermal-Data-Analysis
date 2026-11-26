# Process log: natural-state analysis of Teaching_Model_3D (TE14515)

## 1. User prompt (verbatim)

> Analyse the data

Interpreted as: run the natural-state workflow in `CLAUDE.md` on `Baseline/` and produce the standard deliverables.

## 2. Rules as applied (summary of CLAUDE.md)

- `Baseline/` read only. Scripts in `Scripts/`. Outputs in `Output/Calibration_Input/` (calibration set) and `Output/Additional_Output/` (everything else). Nothing loose in `Output/`.
- Natural state only. Production and injection records used only to test whether surveys represent natural state.
- Inventory first (`input_inventory.csv`, 20 rows covering every file under `Baseline/`).
- Verdicts from the gauge test and the feed-zone trend (tolerance = max(3 bar, 2 x RMS)); confidence carried into weights (high 1.0, medium 0.75, low 0.5).
- Wellbore physics: feed-zone pressures are the formation targets; isothermal sections in open hole only.
- IAPWS-IF97: Region 4 for saturation (own implementation, checked against IF97 verification values and the iapws package, max error 3.9e-7 K); Region 1 for liquid density (iapws package).
- Hydrostatic windows (300 m, close at feed zones or casing shoe), cased and open bands, cold-well component.
- Standard tables with the required leading columns; standard figures F1-F7 with the house style; interpolated companions in `Additional_Output/figures/` titled "[interpolated, not data]".
- Report claims tagged `[obs]`, `[interp]`, `[assumption]` or `[reference model]`; every report number traceable to an output table.
- Acceptance check `Scripts/check_deliverables.py --rerun` run after the final run; no FAIL left; every WARN explained (section 6).

## 3. Steps taken

1. **Orientation.** Read `CLAUDE.md`, `Scripts/check_deliverables.py` (the acceptance check, not edited), the baseline file list and the data formats: well JSON (track, design, feed zones, surveys, transient, flows), plain-column `.dat` files (value then elevation; transient as time in years), `gTE14515.dat` (MULgraph sections: VERTICES fixed width, GRID column polygons, LAYERS, SURFA), `gTE14515_wells.json` (blocks and feed-zone blocks), `Mapping.json` (faults, water bodies). `Logs/` and `README.md` ignored as instructed.
2. **Inventory.** Wrote `input_inventory.csv` (20 rows). Roles: well data JSON and DAT as data; mapping, groups, geometry, well blocks; reference rock, boundary conditions, results and other set-up files as reference only; production, future scenarios, `Area_data.json` and an empty folder as out of scope. `Make_up_json` is an empty directory (no files, noted here).
3. **Cross-check (§4.1 step 1).** All 24 survey and transient series in the `.dat` files match the JSON exactly (`data_crosscheck.csv`). Transient reference date established as 2014-01-01 (time in whole years: t = 6 gives 2020-01-01). Flow and enthalpy `.dat` files were not cross-checked (not part of step 1).
4. **Thermodynamics.** Region 4 implemented (IF97 eq. 31); validated against the three IF97 verification values, the CLAUDE values (99.97 degC at 1.01325 bar, 311.0 degC at 100 bar) and against the iapws package over 1-220 bar: max error 3.9e-7 K, well inside the 0.1 K limit. Injectate temperature from the injection enthalpy (about 335 kJ/kg, median) at 10 bar with iapws Region 1 (see departures).
5. **Targets.** One row per survey point (266 points). Block assignment from the well-to-block file (block whose interval contains the elevation; checked against it in G1). Column, layer, model surface (surface section of the point's column), above-surface flag, feed-zone and casing-shoe flags, open-hole isothermal sections (run of at least 300 m with range at most 5 degC), weights by the §5.5 rules in order.
6. **Feed zones and layers.** 18 feed zones (T, P, Tsat, margin at feed-zone mid elevation by linear interpolation within the survey) and 237 layer-centre rows (interpolated only within each survey range).
7. **Natural-state tests (§4.1).** Gauge test: survey P at gauge elevation minus first reading (precision 0.1 bar). Trend from the gauge-verified wells' feed zones (WELL-2, WELL-3): P = 32.76 + 0.08767 x depth, RMS 0.70 bar, tolerance 3.0 bar. Verdict rules 1 to 5 applied per well in order. Results: WELL-2, WELL-3 natural state (high); WELL-4, 6, 7, 8, 9, 10 probably natural state (medium); WELL-5 probably natural state at best (low; formation about 145 degC, hotter than the 80 degC injectate); WELL-1 possibly disturbed (low; feed zone 3.3 bar below trend, outside the 3 bar tolerance).
8. **Well summary (§4.4).** Computed by the script; the acceptance check recomputes all 160 values independently (D1 PASS).
9. **Hydrostatic windows (§4.3).** Window algorithm as written; densities by IF97 Region 1; implied density, rounding uncertainty; cold wells WELL-4 and WELL-8 (Tmax below 100 degC) give the cold-well component.
10. **Tables.** Pressure consistency, hydrostatic checks (per window and per well/band), well geometry (distances to mapped features), elevation bands (500 m, from the shallowest band to the deepest band with at least three wells), plus supporting tables: T and P at common elevations, injector and isothermal check, weight counts, conceptual-model key, analysis parameters.
11. **Figures.** F1 base map; F2 profile per well and per group panel; F3 band symbol maps (temperature and pressure); F4 section A-A' (through the two hottest wells, WELL-10 and WELL-2, extended to the model extent) for temperature and pressure deviation with data envelope and no-data hatching; F5 P-T state; F6 conceptual model (plan and section, numbered elements); F7 gauge check; interpolated companions (section temperature and pressure, band maps temperature and pressure) in `Additional_Output/figures/`.
12. **Report.** `NS_data_interpretation.md` generated by the script, every number read from an output table.
13. **Figure review.** Viewed after generation: F1, F4 temperature and pressure, F6, F2 WELL-10 profile, F2 group panels, F3 temperature bands, F5, F7, interpolated section (temperature), interpolated temperature band maps. Not individually inspected after generation: the per-well profile files for WELL-1 to WELL-9 (the group panel shows them), the pressure band maps (F3 pressure, interpolated pressure maps). Corrections made from review are in section 4.

## 4. Errors found and how they were fixed

- **Straddling windows.** A window closed at the shoe with a 200 m remnant below it: the remnant joined the window above (the spec), so WELL-4 (230 to -750 mRL) and WELL-2 (190 to -150 mRL) windows straddle the shoe. First version labelled them "open" by bottom elevation, which would have mixed them into the open band. Fixed: labelled "straddles shoe", excluded from both band medians and from the cold component. This produces WARN S2 (see section 6).
- **Cold-well component.** Changing the straddle rule left the cased component with one window (WELL-8) and the open component with two. Reported as such.
- **Figure layout.** Legends clipped or overlapping axis labels (F4, F6); offset labels overlapping (F4); profile text box overflowing (F2); fault and lake labels placed outside the map (F1). Fixed by layout changes and by placing labels at vertices inside the extent.
- **Report text (independent review of key conclusions against the tables).** Corrected: "in every well" (cased ratios exist only for wells with a cased window); "agree within a few percent" (WELL-2 open windows reach 1.12, so the wording now names the wells above 3 percent); a cold-recharge sentence tagged `[obs]` that was an interpretation; element numbers in section 6 (outflow elements 2 to 7, cold elements 8 and 9); wellhead conclusion numbered 8 in cross-references; "no pre-production gauge reading exists" contradicted the table, now stated as an assumption; cold-well component wording (WELL-8 alone); tolerance tag (interpretation, not assumption).
- **Numbers not traceable (R1).** Trend coefficients and the pressure difference in the conceptual key were not in any table; written to `analysis_parameters.csv` and `conceptual_model_key.csv`.
- **Template text (R4).** Two f-string prefixes missing in report lines; fixed.
- **Line endings (S3).** Console log written with CRLF by the Windows console; normalised to LF after the run.
- **Injector sentence.** Replaced a hard-coded placeholder with values read from `injector_and_isothermal_check.csv`.
- **Attempted install.** An attempt to run `pip install iapws` was blocked by the auto-mode classifier. Not retried. iapws was already installed and imports successfully, so nothing depended on the install.

## 5. Geometry self-check (§7)

- Model surface against the surface section: G4 PASS, all 10 wells.
- Blocks against the well-to-block file: G1 PASS (266 survey points), G2 PASS (18 feed zones).
- Layer centres against the layer section: G3 PASS (30 rows, no extrapolation).
- Wellhead elevation against the first track point: identical by construction (the first track point after sorting by elevation is the wellhead).
- Plain-column against JSON: all 24 T, P and transient series match (cross-check table).

## 6. Acceptance check, WARN explanations

- **S2 WARN, band labels other than cased and open ("straddles shoe").** Two windows (WELL-2 and WELL-4) have their top above the production casing shoe and their bottom below it, because the spec's remnant rule joins a short remnant to the window above. They fit neither band. They are kept in `hydrostatic_check.csv` with the label and are excluded from the band medians. Recorded as a departure in section 7.

## 7. Departures from the instructions, with reasons

1. **First reading on the own first-production date treated as pre-production.** For WELL-2, WELL-3 and WELL-9 the gauge series begins on the first flow date. Strictly, no reading predates production. These readings are treated as pre-production (reading date not after own first flow) because they are the earliest record and represent the drawdown-free state; without this, WELL-2 and WELL-3 would be downgraded to medium and the rule-1 test could not be applied to any well. This assumption drives the high verdicts of WELL-2 and WELL-3 and is flagged in the report (section 7).
2. **Gauge reading precision of 0.1 bar** is an assumption (one unit of the recorded values). Tagged `[assumption]` in the report.
3. **Region 1 liquid density and the injectate temperature** use the iapws package (an IF97 implementation). Region 4 is implemented in the script and validated. The injectate temperature uses 10 bar, assumed (subcooled liquid, weak pressure dependence). Tagged `[assumption]`.
4. **Straddling windows** labelled "straddles shoe" rather than forced into cased or open (see section 6).
5. **Section line.** One section (A-A') through the two hottest wells. The optional reference-model companion (§8.3) was not produced, since the reference model is not evidence and the spec makes it optional.
6. **Cap not drawn.** The data do not resolve a cap; F6 shows it as absent with the decisive observation stated in the report.
7. **Area_data.json** (field totals) classified out of scope and not read. First production in the field is taken from the well flow records (earliest positive flow record, 2014-01-01).
8. **Gauge-verified wells for the trend.** Taken as the wells passing rule 1, which are WELL-2 and WELL-3 only, so the trend rests on five feed zones.

## 8. Acceptance check output

(Final run of `python Scripts/check_deliverables.py --rerun` after this log was written. Copied from `acceptance_check.txt`.)

```
PASS  I1   20 inventory rows, all roles valid and paths present
PASS  I1   every Baseline file is covered by the inventory
PASS  F1   all 14 required files present in their folders
PASS  F1   all standard figure types present
PASS  F2   figure index lists 19 files covering ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7']
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
WARN  S2   hydrostatic band labels other than ['cased', 'open']: ['straddles shoe']
PASS  S3   UTF-8, integer feed zones, ISO dates
PASS  S3   LF line endings
PASS  T1   Tsat agrees with an independent IF97 Region 4 within 0.1 K at every point
PASS  T2   boiling margins consistent
PASS  D1   all 160 well_summary values match the §4.4 definitions
PASS  W1   weights valid
PASS  W2   ns_confidence takes 3 values across wells
PASS  W3   weights follow §5.5; verdict and confidence agree across files
PASS  R1   0/389 report numbers (0%) not found in any output table; first: []
PASS  R2   analysis script writes the report
PASS  R4   no unfilled template text in the report
PASS  P1   process log has prompt, departures and acceptance check
PASS  G1   all 266 survey-point blocks match the well-to-block file
PASS  G2   all 18 feed-zone blocks match
PASS  G3   all 30 layer centres match the grid file
PASS  G3   no layer-centre row lies outside its well's survey range
PASS  G4   model_surface_mRL matches the surface section for all 10 wells
PASS  R3   re-run reproduces every CSV and the report

0 FAIL, 1 WARN, 37 PASS
```
