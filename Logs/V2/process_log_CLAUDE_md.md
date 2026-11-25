# Process log: writing an improved CLAUDE.md

## Prompt

> Given the already provided instructions and the example data provided, write an improved CLAUDE.md file into Output

## Steps

1. **Read the current root `CLAUDE.md`** and listed the repository. `Output/` and `Scripts/` had been deleted from the working tree but are still in git history. `Baseline/` was only read, never edited.
2. **Inspected the data formats in `Baseline/Teaching_Model_3D/Data`:**
   - the JSON schema for each well (feed zones, track, design, downhole T and P, flow and transient records);
   - the column order in the `.dat` files;
   - `Groups_data.json`, `Mapping.json` and the plotting JSON files.
3. **Checked when each well's records start:**
   - all downhole surveys are dated 29/07/2015;
   - production and injection start on 01/01/2014 for WELL-2, 3, 4, 5, 8 and 9;
   - WELL-10's records start in 12/2019 (production) and 01/2020 (transient pressure).
4. **Inspected the model geometry and setup files:** the sections of `gTE14515.dat`, `gTE14515_wells.json` and `currentGeomSetup.json`.
5. **Checked the Python environment:** numpy and matplotlib were available; pandas, scipy and iapws were not.
6. **Compared the three earlier analysis runs in git history** (3352c8a, 4c18c22 and fe15a96) to find where they diverged or went wrong:
   - **Natural state.** The runs disagreed on whether the 2015 surveys are natural state, because only two runs compared them with the first transient pressures.
   - **Data vs reference.** The boundary-condition file was used as data in one run.
   - **Water density.** One run used the Kell density above its 150 °C validity range, and applied a liquid hydrostatic test to the near-boiling WELL-10.
   - **Wellbore pressure.** No run separated the pressure in the wellbore column from the formation pressure at feed zones.
   - **File names and coverage.** Output file names and coverage differed between runs, so the runs could not be compared.
7. **Wrote `Output/CLAUDE.md`.** It keeps every original rule and adds:
   - a table classifying each file as data, geometry or reference;
   - known quirks of the dataset;
   - a required test of whether the surveys are natural state;
   - guidance on wellbore physics and thermodynamics (IAPWS-IF97);
   - a fixed set of deliverables, with the columns for the targets table;
   - a report structure;
   - a quality bar.
