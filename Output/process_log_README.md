# Process log: writing Output/README.md

## Prompt

> Ignoring the instruction to ignore the README file, write a new README file into Output with relevant details and required data formatting for the current analysis process

## Rules applied

- **README.** The CLAUDE.md rule "Ignore `README.md`" was overridden by the user for this task only. I read the existing root README, which is a two-line description of the repository, and left it unchanged.
- **Baseline.** `Baseline/` was only read.
- **Output location.** The new file was written to `Output/`.

## Steps

1. **Read the root `README.md`** and listed `Output/`.
2. **Confirmed the formats not checked in the previous task:**
   - **Plain-text files.** The `.dat` flow and enthalpy files have time in years since 01/01/2014, and production flow is stored as negative.
   - **Model geometry.** In `gTE14515.dat`, the `LAYERS` rows are `index, bottom, centre`; the `SURFA` rows are `column, elevation`; the `GRID` rows are `column, vertex count`, followed by the vertex ids.
   - **Well JSON.** The JSON structure of the pressure and transient-pressure blocks was checked.
   - **Mapping.json.** The keys in `Mapping.json` were listed.
3. **Wrote `Output/README.md`.** It covers:
   - the repository layout;
   - how to run the analysis;
   - how each input is classed (data, geometry, reference);
   - the required input formats, including the parsing rules for rows stored as Python-literal strings;
   - the standard output formats, consistent with `CLAUDE.md` §5;
   - the dataset caveats.

## Departure from the instructions

`CLAUDE.md` names the process log `Output/process_log.md`. I used a separate file name for this log, so that it does not overwrite the log of an analysis run.
