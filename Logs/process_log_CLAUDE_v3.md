# Process log: dataset-agnostic CLAUDE.md v3 draft and acceptance check

## User prompt (verbatim)

> Repeat this but ensure the CLAUDE file is dataset agnostic. The information provided can be assumed to be the same but specific dates, values and file names cannot.

"This" is the previous request, "Draft both into Output": a CLAUDE.md v3 draft and an acceptance-check script. Those earlier drafts were moved by the user to `Logs/`. They were not read from there; their content was available in the session.

## Rules applied

- `Baseline/` was only read: by the checker during tests, and by one scan of its file names.
- `Logs/` was not read. Its files are referred to by name only, in the changes note.
- Outputs went to `Output/`, as the user asked.
- The live `CLAUDE.md` and `Scripts/` were not modified.

## Steps

1. Found that `Output/` was empty and the earlier drafts were in `Logs/`. Wrote fresh drafts into `Output/`.
2. Rewrote the CLAUDE.md draft around input **roles** instead of file names:
   - a new required `Output/input_inventory.csv`;
   - hydrostatic bands defined by each well's casing shoe instead of a fixed elevation;
   - generic distance columns for mapped features;
   - rank-based quality-tag factors;
   - generic figure names and examples.
3. Rewrote the checker so that it takes the grid and well-to-block paths from the inventory. Added check I1 (inventory validity and coverage) and check S2 (one pressure row per well; band labels). Feed-zone matching now uses order instead of literal key names.
4. Tested the checker on the three v2 runs (scratch copies, each with a test inventory), with `--rerun`. It catches every Haiku error. **Error found:** the report-number check R1 drew numbers from the previous `acceptance_check.txt` in the same folder, which matched the very numbers it had flagged. **Fixed:** that file is excluded from the pool. Re-tested.
5. Negative test: an inventory missing the well-to-block role, with one non-existent path, fails I1 and the geometry checks as intended. The test copy was deleted afterwards.
6. Scanned both drafts for every Baseline file and directory name, mapped feature name, well name, model name, year, date, block-name pattern and the old band elevation. The only hits were generic words ("data", "production", "river", "mapping").
7. Wrote `Output/CLAUDE_v3_changes.md`: the agnostic rewrite, the v2 failures each change addresses, the checker test results, and the limitations.

## Departures

- **The draft script is in `Output/`, not `Scripts/`**, because the user asked for the drafts in `Output`. Its intended location is `Scripts/check_deliverables.py`.
- **The changes note refers to v2 evidence** (run names and counts). It is meant for the user, not for a model, so only the two drafts themselves are dataset-agnostic.
