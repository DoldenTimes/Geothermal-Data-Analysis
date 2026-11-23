# Process log: natural-state analysis

## Prompt

User request (Claude Code, repository `C:\Users\Dold718\GDA Files`):

> Analyse the data

Governing instructions from `CLAUDE.md` (applied to this run):
- Initial data from `Baseline/`; no edits to `Baseline/`.
- Use the provided downhole data and geometry. The model files are reference only, not factual data.
- Scripts go in `Scripts/`; data outputs go in `Output/`.
- Include the prompt and all steps taken in a separate file in `Output/` (this file).
- Purpose: interpret the data so an engineer or agent can calibrate a geothermal reservoir model.
- Natural state only. Production and future states are out of scope for now.

## Steps taken

1. **Survey the repository.** Listed `Baseline/`, `Output/`, `Scripts/`. `Output/` and `Scripts/` were empty in the working tree (the earlier results were deleted, per the git status at session start; they remain in git history at commit `4c18c22` and were not restored or used).
2. **Identify inputs.** Data used: `Baseline/Teaching_Model_3D/Data/Data_json/WELL-n_data.json` (geometry, feed zones, casing design, downhole T and P, flow records) and `Data_dat/*.dat` (same P/T, as a cross-check). `Mapping/Mapping.json` provides the lake and fault overlays for the map. `Models/` and `Parameters/` were not used as data.
3. **Inspect file structure.** Confirmed the JSON schema, the dat column order (value, then elevation), and the date and flow coverage per well. Found one survey date (29/07/2015) for all downhole P/T, all flagged stable. Found flow records from 2014 for WELL-2, 3, 4, 5, 8, 9, and from 2019 for WELL-10.
4. **Check available tools.** Python with numpy and matplotlib only. pandas, scipy, iapws and CoolProp are not installed, so the saturation curve and water density were implemented directly.
5. **Write `Scripts/natural_state_analysis.py`.** The script:
   - loads geometry and downhole P/T for all wells;
   - cross-checks the dat files against the JSON;
   - checks feed zones and measurements against the well tracks and liners;
   - builds a long table of T, P, boiling-curve temperature (Tsat), Tsat−T margin, feed-zone flags and geothermal gradient;
   - tests hydrostatic consistency using Kell water density, integrated from the shallowest measurement;
   - reduces pressures to datum 0 masl, for cross-well comparison;
   - summarises feed zones and interpolates T and P at common elevations (no extrapolation);
   - writes CSV, JSON and PNG outputs.
6. **Validate the saturation function.** Compared the Wagner-type backward equation (IAPWS-IF97 Region 4) against reference saturation temperatures at 1.01325, 1, 10 and 100 bar. The script fails if any error exceeds 0.5 K.
7. **First run and corrections.** Three issues were found and fixed during the run:
   - The first validation table used wrong reference values (I had assigned 179.88 °C to 1 bar and 311 °C to 10 bar). The function was correct. Corrected the references to 1.01325 bar → 99.97 °C, 1 bar → 99.61 °C, 10 bar → 179.88 °C, 100 bar → 311.00 °C.
   - The hydrostatic summary multiplied a percentage by 100 again, giving values in the thousands of percent. Fixed, and added a whole-column comparison (measured mean density vs Kell mean density over the same interval).
   - The temperature gradient sign was reversed in the long table and feed-zone summary. Fixed so that positive means warmer with depth.
   - The mapping overlay failed on ragged polylines. Fixed by handling single and multi-segment items.
8. **Final run.** Script ran to completion. Console output saved to `Output/run_console_log.txt`.
9. **Check outputs.** Read `hydrostatic_check.csv`, `pressure_datum_reduction.csv`, `feedzones_summary.csv`, `T_P_at_common_elevations.csv` and `well_geometry_check.csv` against the console log. Spot-checked WELL-2 row by row and checked the common-elevation table against the long table. Viewed figures 02 (profiles) and 05 (maps by elevation) to confirm they render correctly.
10. **Written interpretation.** `Output/NS_data_interpretation.md` summarises the findings, caveats and calibration targets. Numbers in it were checked against the output tables; three wording errors were corrected (feed-zone count, WELL-10 near-boiling depth range, and the description of where the hot and cold wells sit).

## Key decisions and assumptions

- **Natural state** is taken to be the 29/07/2015 downhole surveys, which are the only P/T observations in the dataset. Wells with flow records from 2014 may already be affected by production or injection; this is flagged in the interpretation and not corrected for.
- **Water density**: Kell (1975) polynomial. It is fitted for 0–150 °C, and the data reach 270 °C. The hydrostatic discrepancies found (2–15% in most wells, 38% in WELL-10) are larger than the expected error from this, so the conclusion is considered robust, but an IAPWS-IF97 density should be used to confirm it.
- **Pressure resolution** is 0.1 bar. Densities computed over short segments are noisy, so the column-scale comparison is the main test.
- **Maps** are scatter plots at well locations, not interpolated surfaces, because there are only 10 wells.

## Files produced (in `Output/`)

- `NS_data_interpretation.md`: interpretation and calibration targets
- `process_log.md`: this file
- `NS_summary_stats.json`: per-well summary, including hydrostatic statistics
- `well_geometry_check.csv`
- `dat_vs_json_check.csv`
- `downhole_profiles_long.csv`: every measured point with T, P, Tsat, margin, gradient and feed-zone flag
- `hydrostatic_check.csv`
- `pressure_datum_reduction.csv`
- `well_pressure_consistency.csv`
- `feedzones_summary.csv`
- `T_P_at_common_elevations.csv`
- `run_console_log.txt`
- `figures/01_well_map.png`, `figures/02_T_P_profiles.png`, `figures/03_PT_state.png`, `figures/04_datum_pressure_consistency.png`, `figures/05_T_P_maps_by_elevation.png`

Script: `Scripts/natural_state_analysis.py` (re-run with `python Scripts/natural_state_analysis.py` from the repository root).

## Not done

- No commit was made. The working tree still shows the earlier deleted `Output/` and `Scripts/` files from the session-start git status; this run did not restore them.
- No production or future-state analysis, per CLAUDE.md.
- The transient pressure, flow and enthalpy records were not analysed, since they are outside the natural-state scope.
