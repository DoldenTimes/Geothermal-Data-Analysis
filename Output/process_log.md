# Process Log: Natural-State Data Analysis

**Date:** 2026-11-23
**Prompt:** "Analyse the data"

**Instructions followed (from CLAUDE.md):**
- Take data from Baseline and make no edits to it.
- Use the downwell data and geometry. The model is a reference only.
- Put scripts in `Scripts/` and outputs in `Output/`.
- Analyse the natural state only.
- Write the results for calibration use.
- Record the prompt and the steps taken.

## Steps

1. **Surveyed the repository.**
   - Listed the Baseline contents: 10 well JSON files, the `.dat` profile exports, mapping, the mulgrid geometry `gTE14515.dat`, and the reference model files.
   - Noted that `gTE14515.dat` and `gTE14515_wells.json` had uncommitted changes (a WELLS section added). They were later committed by the user as 7b0992b "Fixed well geometry". This analysis did not modify them.
2. **Inspected the data formats.**
   - Each `WELL-*_data.json` holds the well track, casing design, feed zones (with quality) and one stable T/P survey (29/07/2015). Producers also carry transient P and enthalpy.
   - The `.dat` files duplicate the JSON profiles.
   - The mulgrid holds VERTICES, GRID (columns), LAYERS, SURFA and WELLS sections.
3. **Chose the scope.**
   - Used only the downwell T/P, feed zones, tracks and casing, plus the geometry (columns, layers, surface) and map features.
   - Did not use the reference model's rock types, boundary conditions or NS results.
   - Used production data only to check whether the 2015 profiles are pre-production (natural state).
4. **Wrote `Scripts/natural_state_analysis.py`** (numpy and matplotlib only). It:
   - parses the mulgrid geometry and finds the model column (point-in-polygon) and layer for every datum, following the deviated tracks of WELL-8 and WELL-9;
   - computes Tsat(P) with IAPWS-IF97 region 4 to get the boiling margin at each depth;
   - classifies each profile (reversal, boiling, conductive, cold downflow);
   - estimates the deep gradient and water level for each well;
   - compares each interval's dP/dz with ρ_liquid(T)·g to flag vertical flow;
   - compares profile pressures with the first transient pressures (the natural-state check);
   - interpolates T and P to the model layer centres;
   - assigns suggested calibration weights;
   - produces the maps, profiles, per-well panels, elevation-slice T/P maps and two indicative cross-sections.
5. **Ran the script and fixed two problems.**
   - The map river features are multi-part polylines, so their handling was fixed.
   - The first regime classification was too coarse, so it was refined and re-run.
   - Ylims and labels were adjusted after reviewing the figures.
6. **Reviewed the figures and tables** and wrote the interpretation in `Output/NS_data_interpretation.md`. It covers:
   - the upflow location;
   - the boiling zone;
   - the outflow aquifer;
   - recharge;
   - the pressure field;
   - recommended targets and model structure;
   - data issues.

## Reproduce

Run from the repository root:

```
python Scripts/natural_state_analysis.py
```
