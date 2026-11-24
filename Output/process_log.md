# Process log: natural-state data interpretation

## User prompt (verbatim)

> Analyse the data

## CLAUDE.md rules as applied

- `README.md`, `Logs/` and the git history of earlier runs were ignored. Nothing was restored from previous runs. All inputs were read from `Baseline/Teaching_Model_3D/`, which was not modified.
- Script: `Scripts/natural_state_analysis.py`. Outputs: `Output/`, figures in `Output/figures/`.
- Scope: natural state only. Production and injection histories were used only to test whether the surveys are pre-production (§3.1).
- **Data (evidence):** `Data_json/WELL-*_data.json`; `Data_dat/*.dat`, used only for cross-checks; `Mapping.json`; `Groups_data.json`.
- **Geometry only:** the VERTICES, GRID, LAYERS, SURFA and WELLS sections of `gTE14515.dat`, and `gTE14515_wells.json`, which was used to cross-check feed-zone blocks.
- **Reference model only:** `gTE14515_boundary_conditions.json` and the GENER block of `NaturalState/TE14515_NS_000.dat`. They are described in the report and labelled **[reference model]**; nothing from them is used as evidence. `gTE14515_dict_rocktype.json`, `.incon`, `current*Setup.json` (except reading the NS setup file once to understand the folder), `Parameters/` and `Data/Plotting/` were not used.
- **Assumption stated:** rock thermal conductivity k = 2.0 W/m/K, used only to convert gradients into conductive heat flux.
- **Thermodynamics:** IAPWS-IF97 Region 4 (Tsat, Psat) and Region 1 (density, enthalpy), implemented in the script. They are validated against the IF97 verification tables and the CLAUDE.md values; the script stops if any Tsat error exceeds 0.1 K. No Kell polynomial or low-temperature fit is used.
- **Hydrostatic tests:** single-phase liquid only (P ≥ 1.5 bar, boiling margin ≥ 10 °C), windows ≥ 300 m. The field-wide component (cold wells) is separated from the per-well excess, and the analysis is split into depth bands.
- **Interpolation only:** no extrapolation in the maps, sections, common-elevation tables or layer-centre tables.
- **Weights:** CLAUDE.md defaults (1.0 / 0.8 / 0.5 / 0.3 / 0), multiplied by the natural-state confidence (high 1.0, medium 0.8). One class was added (see departures).
- **Quality bar:** every number in the report is produced by the script from the output tables. Every figure was opened and checked after generation.

## Steps

1. Inspected the structure of every Baseline input: well JSON keys, dat file formats (the T and P dat files are `value elevation`; the transient-P dat files are `years-since-2014 value`), Mapping.json (rivers are multi-segment "multiplot" lists), and the geometry file sections (fixed-width VERTICES; GRID column headers followed by vertex names; LAYERS `name bottom centre`; SURFA per column).
2. Printed every T/P profile with Tsat and density for a first look. Found that the WELL-1 wellhead pressure is 10.8 bar at 15 mRL (artesian), that WELL-10 is near boiling, and that pressures above 220.6 bar have no Tsat (the boiling margin is set to NaN and flagged).
3. Looked at the GENER block of the reference model: only basal HEAT and near-surface rain MASS generators, with no deep upflow. This is reported as [reference model] only.
4. Wrote `Scripts/natural_state_analysis.py`: IF97 and its validation, loaders, geometry (point-in-polygon column search, layer lookup, top-block handling), integrity checks, natural-state test, per-point target table, geometry check, feed zones, hydrostatic test, per-well summary, layer-centre and common-elevation tables, figures, and the report writer.
5. **Error:** `Mapping.json` rivers are ragged lists of segments, so `np.array` failed. **Fixed:** each river is stored as a list of segment arrays.
6. **Finding checked before reporting:** the shallow hydrostatic ratios reach 1.54 (implied density 1363 kg/m³), which is physically impossible for a static water column. Checked: datum offset (cannot change a gradient); g or unit convention (would act at every depth, but the WELL-5/6/7 deep windows are within 0.3 % of IF97); two-phase fluid (would lower the gradient, not raise it); 0.1 bar rounding (below 0.5 % over 300 m); elevation reference (the wells are vertical except WELL-8/9, whose elevations are already TVD). Conclusion: the anomaly is confined to above about −450 mRL and grows with temperature. It is reported as a data issue with an interpretation of low confidence, and pressure weight is applied at feed zones only. A depth-band split was added to the hydrostatic analysis (`hydrostatic_band_summary.csv`).
7. **Refined rules after the first run:**
   - The deep gradient was first computed below the deepest feed zone and picked up the inversion, giving negative values for WELL-3. It now uses the static interval below both the deepest feed zone and the inversion base. If that interval is too short it falls back to below the inversion base, and if there is no inversion to the deepest 500 m. It is not reported for the near-boiling WELL-10.
   - The regime rule labelled WELL-7's 46 °C inversion "weak". It now separates outflow (≥ 20 °C inversion) from weak distal outflow (5–20 °C).
8. Viewed every figure. Fixes: fault labels were drawn off-map, so they are now placed on the visible trace; well traces sat on the section-plot edge, so x-margins were added.
9. Added the IF97 Region 1 enthalpy, validated against the IF97 table, for the upflow enthalpy and the feed-zone enthalpies.
10. Wrote the report generator inside the main script, so that a single script produces every output. Reviewed the generated report against the tables and corrected statements the data did not support:
    - "inversion maxima at the shallowest feed zone" is now the computed distance per well; for WELL-5 the maximum is behind casing, 650 m from a feed zone;
    - "the 200 °C isotherm is behind casing in all hot wells" now names the exception, WELL-9;
    - "within 2 bar of the median from −500 to −1500 mRL" is now a computed count per elevation;
    - "outflow crosses Fault 2 at WELL-9/7" is now "reaches WELL-5, east of Fault 2";
    - "WELL-10 next to the F1×F3 intersection" now gives its distances to each fault (758 m from the intersection);
    - fixed numeric ranges (cold-water temperatures, the range of the survey-date drawdown) are now computed values.
11. **Shell issues:** one bash heredoc failed to parse. The patch was written to a scratchpad Python file and run from there. A heredoc patch also left a literal newline inside a `print` string, which was fixed.
12. **Blocked action:** a `rm -rf Output/*` to clear old outputs before the final run was refused by the safety check. Instead, the final run overwrites every output. The output files from earlier runs that this run does not produce had already been deleted from the working tree at the start of the session (git status), so none remain.
13. Final run: `python Scripts/natural_state_analysis.py`, exit code 0. The console output is in `run_console_log.txt`.

## Departures from CLAUDE.md, with reasons

- **Extra weight class:** "open hole, possible internal flow", weight 0.4. Four WELL-2 points between −1100 and −1900 mRL are isothermal between feed zones, the wellbore-flow signature that §3.2 asks to be considered.
- **Extra columns in `NS_calibration_targets.csv`**, after the required columns, in the required order: `location_class`, `weight_P` (pressure weight, non-zero only at feed zones), `rho_IF97_kgm3` and `flags`. CLAUDE.md asks for pressure targets only at feed zones, which a single `weight` column cannot express for both T and P.
- **"Near wellhead / top layer"** is defined as within 30 m of the wellhead, or in the top model block of the column.
- **Field-wide gradient excess:** CLAUDE.md expects a single field-wide excess. The data show that the excess is not uniform with depth. The cold-well median is reported per depth band, together with the full-run value (2.5 %, labelled a possible datum or convention issue).
- **`process_log.md` is written by hand**, not by the script, because it records the process. Every other output is regenerated by the script.
- **Section references:** CLAUDE.md §4 and §5 refer to "§4.1", "§4.3" and "§4.4", which match the content of §3.1, §3.3 and §3.4. They were read that way.
- **Additional outputs** (allowed): `data_integrity_check.csv`, `T_P_at_common_elevations.csv`, `hydrostatic_band_summary.csv`, `key_numbers.json`, and the extra figures (all-well profiles, hydrostatic ratio, feed-zone pressures, section C).
