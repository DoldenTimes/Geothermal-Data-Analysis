# Process log

## Prompt
User request (verbatim): **"Analyse the data"**
Governing instructions: `CLAUDE.md` (README.md ignored; read only from `Baseline/`, no edits there; scripts → `Scripts/`; outputs → `Output/`; natural state only; supplied model is reference only; include prompt and all steps in a separate Output file).

## Steps
1. Listed the project: `Baseline/Teaching_Model_3D` contains `Data` (well JSON/DAT, groups, mapping), `Models/TE14515` (geometry `gTE14515.dat`, `_wells.json`, `_boundary_conditions.json`, rock-type dictionary, natural-state model and incon, FutureScenarios/Production setup) and `Parameters/ProductionHistory` (xlsx). Git status showed an earlier set of Output/Scripts files deleted in the working tree; they were re-created from scratch rather than restored.
2. Inspected formats: per-well JSON (feed zones, track, design, one downhole T and P survey dated 29/07/2015, production/injection rate, enthalpy, transient pressure); geometry file sections (VERTICES fixed-width 3/10/10 chars, GRID column→vertex lists, LAYERS bottom/centre, SURFA surface elevation, WELLS tracks); well→block mapping with feed-zone weights; boundary conditions (atmosphere, rain, heat flux, lake columns); GIS polylines (lake, river, faults).
3. Checked the environment: only numpy and matplotlib available (no pandas/scipy/iapws). Implemented IAPWS-IF97 region-4 saturation (Psat/Tsat), a steam-table interpolation for saturated-liquid enthalpy and a saturated-liquid density with a compressibility term.
4. Wrote `Scripts/natural_state_analysis.py` (fixes during development: column centroids from GRID rather than vertex lines; GIS entries with multiple polylines; run without `-I` because that flag hides user site-packages).
5. Ran it. It computes per-well metrics (Tmax, inversion, gradients, T=100/150/200 °C elevations, boiling margin, estimated water level), feed-zone T/P, T-integrated hydrostatic pressure and residual, survey-vs-first-transient pressure, initial enthalpy → liquid-T comparison, well–block T/P targets, T/P at common elevations, collar vs model surface, and distances to faults/rivers/lake. Produced 7 figure files and 8 CSVs.
6. Reviewed key figures (map, per-well panels) and cross-checked a suspicious result: the hydrostatic residual is positive in every well, including cold W4/W8, so the implied pressure gradients were inspected (≈2–3 % steeper than cold-water hydrostatic, partly rounding noise). Recorded as a caveat rather than a defect.
7. Wrote `NS_data_interpretation.md` summarising data, natural-state structure, targets, and caveats.

## Limitations
* One survey per well; no heat-up series, so formation vs. flowing temperature cannot be distinguished.
* Water levels are extrapolated from the survey top; W4/W10 surveys start below the collar.
* Liquid-density approximation for the hydrostatic check is accurate to ~1 %.
* The supplied natural-state model and lithology/rock-type properties were not used as data.
* Interpretations of upflow location and cap geometry are hypotheses from 10 vertical-ish wells.
