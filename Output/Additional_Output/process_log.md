# Process log

## User prompt (verbatim)

> Analyse the data

## Rules from CLAUDE.md as applied

- `README.md` and `Logs/` ignored; `Baseline/` only read; scripts in `Scripts/`; outputs only under `Output/Calibration_Input/` (report, targets, standard figures, index) and `Output/Additional_Output/` (everything else). Previous results in git history were not restored or reused.
- Scope: natural state only. Production/injection records inside the well files were used only for the natural-state test (gauge vs survey, injectate temperature).
- Roles assigned to every file (`input_inventory.csv`); reference model files (rock, boundary conditions, natural-state results, set-up, plotting) were not read and not used as evidence; production/forecast files out of scope.
- IF97 Region 4 (Tsat) and Region 1 (density, enthalpy) implemented in the script and validated against the IF97 verification tables, Tsat(1.01325 bar) = 99.97 degC and Tsat(100 bar) = 311.0 degC; the script stops if an error exceeds 0.1 K (maximum Region 4 error 0.004 K).
- Defaults followed in the standard columns: §4.4 definitions (summary values recomputed by the checker from the targets table), §4.3 hydrostatic window algorithm with cased/open bands, §4.1 verdict rubric and tolerance, §5.5 weights, §8 figure layouts. Extras are in extra columns / Additional_Output.
- Report: every claim tagged [obs], [interp], [assumption] or [reference model]; every quoted derived number is also written to `report_numbers.csv`.

## Steps taken

1. Listed `Baseline/`: 10 well JSON files, 38 plain-column files, mapping (3 faults, Tarawera River, river, Lake Tarawera), groups file, MULgraph grid (398 columns, 70 layer lines), well-to-block file, reference model files, production/forecast files. Wrote the inventory from the script.
2. Parsed the grid file by section header (layers by whitespace; vertices by fixed width, since coordinates run together). Confirmed: layer 0 = 1000 mRL top, layers tile without gaps, surface section has 398 entries.
3. Parsed wells (literal evaluation of string-encoded lists). One survey date per well (2015-07-29), identical T and P elevations.
4. Cross-checked JSON against plain-column files. Found: `Wel10_*` file names (irregular); production flow is stored negative in the .dat files (generation-rate sign convention) - first comparison reported a false mismatch, fixed by comparing magnitudes and recording the convention; time axis is years since 2014-01-01.
5. Geometry: blocks from the well-to-block file are the authority; an independent point-in-polygon search over the grid polygons agrees at every survey point. Model surface taken from the surface section at the column of the top block.
6. Natural-state test: gauge series start 2014-01-01 (WELL-2, -3, -9) and 2020-01-01 (WELL-10). WELL-2 and WELL-3 miss the first gauge reading by 0.10 bar, equal to the reading precision (0.1 bar), so pass at the limit. WELL-9 misses by 2.05 bar. WELL-10 gauge starts after its first flow. Injectate temperature 77-78 degC from the 335 kJ/kg enthalpy.
7. Error found and fixed during the run: the first injector rule used a 15 degC window for an isothermal section "near" the injectate temperature, which labelled WELL-8 (isothermal section 10-14 degC colder than the injectate, which injection cannot produce) as injection-affected. The window was tightened to the 5 degC isothermal tolerance; WELL-8 now follows the trend rows (medium).
8. Weights, summary, layer centres, hydrostatic windows, elevation bands, geometry table, pressure connections, report. Figures F1-F7 and the interpolated companions were generated and each was looked at (F1, F2 profile and panels, F3 temperature, F4 section A, F4 B pressure interpolated, F5, F6, F7). Fixes after looking: fault labels moved away from wells, profile spacing, section title wrapping, F6 title padding, overlapping pressure labels offset.
9. Section lines: A through the two hottest wells (WELL-10, WELL-2), extended to the model extent, wells within 1000 m; B through the hottest well, perpendicular to A, wells within 600 m (stated in the report and `section_offsets.csv`).

## Geometry self-check

- Model surface vs surface section: matches for all 10 wells (G4 PASS).
- Blocks vs well-to-block file: all 266 survey points and 18 feed zones match (G1, G2 PASS).
- Layer centres vs layer section: all 30 distinct layers match; no extrapolation (G3 PASS).
- Wellhead elevation vs first track point: equal for all wells (`well_geometry_check.csv`). Wellheads lie below the model surface of their column in 8 wells (up to 105 m, WELL-1): recorded as a data question, not changed.

## Independent review (try to disprove) and changes made

- Conclusion 1 originally called every non-gauge well "probably natural state at best"; WELL-1 is "possibly disturbed", so the text now lists each verdict.
- The 'thin cap' and 'fed from the fault and the river' wording was not supported (WELL-4 is 3 km from the river); softened to a conductive cap and "possibly".
- WELL-10 'about 300 m near boiling' was wrong (187 m between the points with margin below 5 degC); now computed from the table.
- Water levels: the text said most are extrapolated; six are, four lie within the survey; the text now says so.
- WELL-1 'possibly disturbed' rests on a 0.3 bar exceedance of the 3 bar tolerance; kept because the rubric is binding, flagged as marginal in the report.
- Shallow cased-interval pressure gradients exceed that of water in the hot wells; no datum, convention, two-phase or rounding explanation fits; reported as unexplained and excluded from pressure targets.

## Departures from the instructions

- Hydrostatic window markers: a feed zone is "reached" at its top elevation; windows that run past the shoe because the remnant is under 300 m are labelled by the larger share of their length and flagged in `straddles_shoe` (rule not specified in CLAUDE.md).
- Gauge "pre-production" means the first gauge date is not later than the well's first flow record and the field's first production (same-day readings count); reading precision taken as 0.1 bar (stated resolution of survey and gauge values).
- Top block of a column = the layer whose interval contains the model surface (CLAUDE.md does not define it).
- Optional reference-model companion figure (§8.3) not produced: no reference model file was read.
- The F4 sections plot x-limits are cropped to the wells +/- 1.2 km; the lines themselves extend to the model extent.

## Acceptance check

Command: `python Scripts/check_deliverables.py --rerun`

```
PASS  I1   18 inventory rows, all roles valid and paths present
PASS  I1   every Baseline file is covered by the inventory
PASS  F1   all 14 required files present in their folders
PASS  F1   all standard figure types present
PASS  F2   figure index lists 21 files covering ['F1', 'F2', 'F3', 'F4', 'F5', 'F6', 'F7']
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
PASS  R1   25/665 report numbers (4%) not found in any output table; first: ['2.2', '72.7', '266.1', '241.1', '243.3', '152.9', '154.5', '218.8', '163.5', '153.9', '61.9', '195.2']
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

WARN explanation: none outstanding in the final run (the only WARN in earlier runs was P1 while this log was a placeholder).
