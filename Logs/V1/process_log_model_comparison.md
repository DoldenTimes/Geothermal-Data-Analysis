# Process log: comparing the Opus1, Sonnet1 and Haiku1 runs

## Prompt

> Accessing the previous runs using various models in Opus1, Sonnet1 and Haiku1, analyse the differences between model performances and write a report on this into Output

## Rules applied (`CLAUDE.md`)

- `Baseline/` was only read.
- The check script is in `Scripts/`, and the report and check output are in `Output/`.
- The runs' outputs were used only as material to assess, never as inputs to a new analysis.

## Steps

1. **Located the runs.** Opus1, Sonnet1 and Haiku1 are not folders. They are remote branches, `origin/Opus1`, `origin/Sonnet1` and `origin/Haiku1`, pointing at commits `3352c8a`, `4c18c22` and `fe15a96` on `main`.
2. **Checked that the comparison was fair:**
   - `CLAUDE.md` has the same md5 at all three commits;
   - `Baseline/` has no differences between the commits;
   - all three process logs give the same prompt;
   - the process logs show the same environment.
3. **Read each run's material:** the process log, the interpretation report, the CSV headers and row counts, and the list of figures.
4. **Re-ran each run's script.** I extracted each branch with `git archive` into a scratch folder (outside the repository), cleared its `Output/`, ran `python Scripts/natural_state_analysis.py`, and compared the regenerated files with the committed ones:
   - every CSV and JSON file matched exactly;
   - Haiku's console log differed only in the output path.
5. **Wrote an independent check, `Scripts/model_comparison_check.py`,** and saved its output to `Output/model_comparison_check.txt`. It computes, directly from `Baseline/`:
   - Tmax, T at total depth and inversion size;
   - the minimum boiling margin (IF97 Region 4, validated);
   - pressure at −250, −750 and −1500 mRL;
   - survey vs first transient pressure;
   - the ratio of each column's pressure gradient to saturated-liquid density;
   - the apparent density of the shallow column;
   - wellhead vs model surface, via point-in-polygon on the grid.
   I also worked out the Fault 1 / Fault 3 intersection by hand from `Mapping.json`.
6. **Fixed one bug in my check script.** The GRID header line is `aab0 4`: the vertex count is the last field, not the first. Fixed, then re-ran.
7. **Viewed the T/P profile figure** from each run.
8. **Wrote `Output/model_comparison_report.md`.** It covers:
   - a summary and ranking;
   - fairness;
   - accuracy against the independent check;
   - interpretive differences;
   - deliverables;
   - reproducibility;
   - compliance with the instructions;
   - limitations;
   - lessons mapped to the current `CLAUDE.md`.
9. **Corrected the section references.** The current `CLAUDE.md` numbers the deliverables as §4, not §5, so I updated the references in the report.

## Departures from the instructions

- **Process log name.** I wrote this log to a separate file so it does not overwrite the `process_log.md` of an analysis run.
- **Density check.** The check script uses saturated-liquid density, not the IF97 Region 1 density that `CLAUDE.md` §3.3 requires for analysis runs. It is only a coarse (±5 %) check of the runs' claims, and the report says so.
