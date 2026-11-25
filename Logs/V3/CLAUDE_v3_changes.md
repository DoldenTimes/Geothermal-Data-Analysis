# CLAUDE.md v3 draft (dataset-agnostic): what changed and why

Drafts in `Output/`:
- `CLAUDE_v3_draft.md`, intended to replace `CLAUDE.md`;
- `check_deliverables_draft.py`, intended to become `Scripts/check_deliverables.py`.

They come from the comparison of the second-iteration runs (`Logs/model_comparison_report_2.md`). The live `CLAUDE.md` and `Scripts/` are unchanged.

**Dataset-agnostic.** Neither draft names a file, well, date, mapped feature or value from this dataset. They assume only that the same *kinds* of information are supplied: per-well JSON with plain-column copies, a mapping file, a well list, a grid file in MULgraph format, a well-to-block file, and reference-model files. They also assume that the information is stored in the same way. Nothing in them states a conclusion about any field.

## 1. How the drafts were made dataset-agnostic

| Previously dataset-specific | Now |
|---|---|
| Input file and directory names in the data/reference table | **Roles** (`well_data_json`, `geometry_dat`, `well_blocks_json`, `reference_bc`, …), described by content. The model must inventory `Baseline/` first and write `Output/input_inventory.csv` (`role, path, note`), with every Baseline file covered (§2.1) |
| Checker paths built from the known directory layout | The checker reads the grid and well-to-block paths **from the inventory**, and checks that the inventory is valid and complete (new check I1) |
| Grid file named by model; well-to-block file named by model | "Model grid file (MULgraph format)" and "well-to-block file" |
| A specific irregular dat file name | "File names do not always follow one pattern; match by content and name, and report irregular names" |
| Transient time column measured from a specific year | "Years since a reference date; establish the date by matching against the JSON series" |
| Hydrostatic band boundary at a fixed elevation (taken from this field's data) | Bands defined per well by its **production casing shoe**: `cased` and `open` |
| `dist_Fault1_m`, `dist_Fault2_m`, `dist_Fault3_m`, a named-river column | `nearest_fault, dist_nearest_fault_m, nearest_water_body, dist_nearest_water_body_m`, plus one `dist_<feature>_m` per mapped feature, using the feature names from the mapping file |
| Feed-zone quality factors for named tags | Mapped by rank (highest, intermediate, lowest or non-grade tag); the report states the mapping |
| `profile_WELL-*.png`, with a fixed count of 10 | `profile_<well>.png`, one per well in the targets table |
| Example claim quoting a real well, temperature and elevation | "well X reaches T °C at z mRL" |
| Feed zones matched to the block file by literal key names | Matched by order within each well; `feedzones.csv` must list feed zones in order |

**Kept, because they are format or physics, not dataset values:**
- the MULgraph layer-line layout;
- string-encoded JSON lists;
- feed-zone units labelled "m" that are mRL;
- the two-column dat layout;
- ragged river arrays;
- the IF97 reference values (99.97 °C at 1.01325 bar, 311.0 °C at 100 bar);
- the critical pressure (220.64 bar);
- the generic thresholds (300 m windows, 3 °C inversion, 100 °C cold-well cut-off, 30 m wellhead zone);
- the default weights.

**Checked mechanically.** Both drafts were scanned for every file and directory name under `Baseline/`, every mapped feature name, every well name, the model name, 4-digit years, dd/mm/yyyy dates, block-name patterns and the previous band elevation. The only hits were generic words that happen to be directory names ("data", "production", "river", "mapping").

## 2. Changes from v2, by the failure they address

| # | Change in v3 | Failure seen in the v2 runs | Where in v3 |
|---|---|---|---|
| 1 | Input inventory by role, written to a file and checked | Runs located inputs ad hoc; nothing recorded which files were read as data and which as reference | §2.1, check I1 |
| 2 | §3 **Input file formats** (layer-line layout, block names, well-to-block structure, string-encoded lists, feed-zone units, dat columns, ragged polylines) | Haiku read the layer lines as bottom/top, so every layer centre was wrong by up to 125 m, with spurious "layer gaps" | §3, check G3 |
| 3 | The well-to-block file is the **authority** for blocks; a disagreement means your own code is suspect first | Haiku put 33 points of one well in the wrong column and reported it as a model problem | §2.2, §3, §7, checks G1 and G2 |
| 4 | **Verdict rule:** the gauge comparison decides, not the date label. Explicit injector logic. Gauge-less wells compared at feed-zone elevations. Uniform confidence must be justified | Haiku downgraded wells whose survey matched the first gauge reading to 0.1 bar, used reversed injector logic, and gave every well the same confidence | §4.1, check W2 |
| 5 | Default **confidence factors** and a rank-based feed-zone quality factor | The three runs used three different factor sets, so their weights could not be compared | §5 |
| 6 | Required **`weight_P`**, non-zero only at feed zones | One run had no pressure weight; another weighted pressures away from feed zones | §5, check W1 |
| 7 | Weight class for open hole inside an isothermal section ≥ 300 m long | Two runs each invented their own class and weight | §5 |
| 8 | **Hydrostatic method fixed:** non-overlapping windows ≥ 300 m, `cased`/`open` bands, cold well Tmax < 100 °C | The field-wide excess came out as 2.5 %, 3.0 % and 8.0 %. Short, pooled windows missed that the deep column is hydrostatic | §4.3, check S2 |
| 9 | **Metric definitions:** shallow and deep gradient, inversion threshold, isotherm crossing | Shallow gradients differed by up to 197 °C/km between runs, and a "deep gradient" was reported for a well only 610 m deep | §4.4 |
| 10 | **Exact column names and order** for every standard CSV; figure-name patterns; one pressure-consistency row per well | Each run named the same quantity differently | §5, checks S1 and S2 |
| 11 | The script writes the report and the inventory; `process_log.md` is the only hand-written output | One report was hand-written while its log said otherwise. All three runs logged the process log as a departure | §1, §5, checks R2 and R3 |
| 12 | **Acceptance check**: no FAIL allowed; WARNs explained in the process log | No run had an enforced gate | §7 |
| 13 | **Provenance:** read and report source, comment and quality fields | Only one run reported a provenance label on the well designs | §2.2 |
| 14 | Tsat blank above the critical pressure; no clamped curve | One run plotted Tsat clamped at the critical temperature | §4.3 |
| 15 | Geometry anomalies treated like physical anomalies | Parsing errors were reported as data-quality issues | §7 |
| 16 | Independent review pass before finishing | The two stronger runs caught their own over-claims this way; the weakest did not | §7 |
| 17 | Section numbering fixed (1–7) | The v2 cross-references pointed at the wrong sections | throughout |

## 3. The acceptance check (`check_deliverables_draft.py`)

The script only reads `Baseline/`, and the only input paths it uses come from `Output/input_inventory.csv`. It writes `Output/acceptance_check.txt` and exits 1 on any FAIL. It has its own IF97 Region 4 implementation, validated against IF97 verification values at start-up.

| ID | Check |
|---|---|
| I1 | inventory: valid roles, paths exist, required roles present (exactly one grid file and one well-to-block file); WARN for Baseline files not covered |
| F1 | required files, and figure types (one profile per well) |
| S1 | required columns present, first and in order |
| S2 | one pressure-consistency row per well; hydrostatic band labels are `cased` or `open` |
| G1, G2 | survey-point and feed-zone blocks match the well-to-block file |
| G3 | layer centres match the grid file; no extrapolation beyond the survey range |
| T1, T2 | Tsat agrees with independent IF97 within 0.1 K; boiling margin = Tsat − T |
| W1, W2 | weights in [0, 1]; zero above the surface; `weight_P` only at feed zones; WARN if confidence is uniform |
| R1 | report numbers found in output tables (FAIL above 25 % unmatched, WARN above 5 %) |
| R2 | the analysis script writes the report |
| R3 | `--rerun`: a clean re-run reproduces every CSV and the report (one differing report line allowed for a timestamp) |
| P1 | the process log mentions the prompt, the departures and the acceptance check |

### Tested against the three v2 runs

A test inventory was added to each run's copy, because v2 did not ask for one. All three runs were checked with `--rerun`.

| Check | Opus2 | Sonnet2 | Haiku2 |
|---|---|---|---|
| I1 inventory | PASS; WARN (test inventory left 16 files uncovered) | same | same |
| G1 survey-point blocks | PASS | PASS | **FAIL** (37/266) |
| G3 layer centres | PASS | PASS | **FAIL** (30 layers, up to 125 m) |
| T1, T2 thermodynamics | PASS | PASS | PASS |
| W1 weights | PASS | **FAIL** (`weight_P` at 234 off-feed-zone points) | **FAIL** (no `weight_P`) |
| W2 confidence | PASS | PASS | **WARN** (uniform) |
| R1 report numbers | PASS (8/752 unmatched) | PASS (6/1139) | PASS (4/528) |
| R2 script writes report | PASS | PASS | **FAIL** |
| F1 conceptual sketch | PASS | PASS | **FAIL** |
| S1, S2, G2, R3 (v3 names, one row per well, `block` column, inventory regenerated) | FAIL (expected) | FAIL (expected) | FAIL (expected) |

The negative test passed: an inventory with a missing well-to-block role and a non-existent path fails I1 and skips the geometry checks with an explicit FAIL. Opus fails only on requirements that are new in v3. Haiku's v2 errors are all caught.

### Known limitations
- **R1 is weak.** It matches each report number against every number in the output tables, so it catches invented numbers but not misattributed ones.
- **G2 depends on order.** It relies on `feedzones.csv` listing feed zones in the same order as the well-to-block file.
- **The grid reader assumes MULgraph format.** A different grid format needs a new reader.
- **The check does not judge interpretation.** Verdict logic, regime labels and the conceptual model still need review. The comparison in `Logs/model_comparison_report_2.md` §3 can serve for that.

## 4. Before adopting
1. Copy `check_deliverables_draft.py` to `Scripts/check_deliverables.py`, and `CLAUDE_v3_draft.md` to `CLAUDE.md`.
2. Run each model two or three times under v3. Score the runs with the acceptance check, and with `Logs/compare_model_runs_2.py` after updating its column mappings to the v3 names, to separate model differences from run-to-run variance.
