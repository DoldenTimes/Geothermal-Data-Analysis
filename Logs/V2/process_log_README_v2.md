# Process log: README rewrite (data structure standard)

## User prompt (verbatim)

> Rewrite the README file to specify the required file structure for efficient data analysis

## Rules applied

- `Baseline/` was only listed and read, to describe its current layout.
- `Logs/` was not read.
- `README.md` was rewritten at the user's request. `CLAUDE.md` tells agents to ignore it, so the README is aimed at the people who prepare data, and points to `CLAUDE.md` for the agent rules.

## Steps

1. Read the old README. It described this dataset's layout, including dataset-specific dates, values and quirks.
2. Confirmed that `CLAUDE.md` is now the v3 (role-based) version and that `Scripts/check_deliverables.py` is in place. The README was aligned with both:
   - the `CLAUDE.md` section numbers (§2.1, §2.2, §3, §5);
   - the role names;
   - the acceptance-check IDs (G1, G2);
   - pandas as a requirement (the checker imports it).
3. Listed `Baseline/` to cover the items the old README omitted: `Area_data.json`, the empty `Make_up_json/`, `Plotting/`, `.auxjson`, the set-up files, the `Production/`, `FutureScenarios/` and `Parameters/` contents, and `.DS_Store`.
4. Wrote the new README:
   - repository layout;
   - the required `Baseline/<Project>/` tree, with the CLAUDE.md role for each item;
   - naming rules (one well name everywhere);
   - required and tolerated formats for each file type, using placeholders only;
   - a preparation checklist;
   - how to run the analysis and the acceptance check;
   - a table of the current dataset's departures from the standard.
5. Checked each departure in §6 against `Baseline/`:
   - `Wel10_*` file names;
   - feed-zone units labelled "m" and the track `z` labelled "masl";
   - rows stored as string literals;
   - a single survey date for every well;
   - the `Synthetic` tag and source label;
   - a liner deeper than its track;
   - wellheads below the model surface;
   - the extra files.

## Departures

- **README §6 is dataset-specific**, because it lists the current data's departures from the standard. The standard itself (§1–5) uses placeholders only. Agents do not read the README, so the dataset-agnostic rule for `CLAUDE.md` is not affected.
