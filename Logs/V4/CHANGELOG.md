# Changelog: CLAUDE.md and acceptance check

## v4 (2026-11-26)

v4 has three goals:
1. Remove the ambiguities that caused the remaining differences between models in iteration 3.
2. Standardise the figures, using the professional model report O'Sullivan et al. (2025), *Geothermics* 130, 103339, as a template.
3. Split the outputs so the calibrator gets a short, self-contained package.

Files changed: `CLAUDE.md` (v3, 220 lines → v4, 413 lines) and `Scripts/check_deliverables.py` (26 → 37 check lines; 15 → 22 check IDs).

### 1. Output folders (new)

Every output now goes in one of two folders. Nothing is written loose in `Output/`.

| `Output/Calibration_Input/`: what the calibrator needs | `Output/Additional_Output/`: supporting and excess information |
|---|---|
| `NS_data_interpretation.md` (main report; must stand on its own) | `input_inventory.csv` |
| `NS_calibration_targets.csv` | `pressure_consistency_check.csv` |
| `feedzones.csv` | `hydrostatic_check.csv` |
| `T_P_at_layer_centres.csv` | `well_geometry_check.csv` |
| `well_summary.csv` | `T_P_at_elevation_bands.csv` (new) |
| `figures_index.csv` (new) | `run_console_log.txt`, `acceptance_check.txt`, `process_log.md` |
| `figures/`: standard figures F1–F7 only | `figures/`: extra figures, including the optional reference-model companion |
| | any other diagnostic or intermediate file |

**Why:** the iteration-3 runs wrote 25–40 files into one folder. The engineer had to pick out the five tables and seven figures that matter. The checks, logs and diagnostics are still kept, but out of the way.

### 2. Definitions fixed (§4.2, §4.3, §4.4)

**General rule (§1).** A defined default is binding in the standard column. Alternative methods go in extra columns or extra figures.

| Item | v3 | v4 | Iteration-3 issue it fixes |
|---|---|---|---|
| Source for per-well values | not stated | rows of `NS_calibration_targets.csv`, top down | makes the values recomputable (D1) |
| Liquid-filled point | undefined | P ≥ 1.5 bar and boiling margin ≥ 0 (supercritical points count) | three different readings |
| Shallow gradient | "top 300 m of the liquid-filled survey" | shallowest liquid-filled point plus all liquid-filled points within 300 m below it; ≥ 3 points | spread of up to 180 °C/km |
| Deep gradient | "below" inversion base and deepest feed zone | at or below (inclusive); ≥ 3 points spanning ≥ 300 m; adds `deep_fit_n` and `deep_fit_span_m` | WELL-6: 7.1 vs 20.2 °C/km |
| Inversion | largest drop below running maximum | the same, with top and base defined explicitly | |
| Isotherm crossing | "first downward crossing" | first consecutive pair with T above < X ≤ T below; linear interpolation | |
| Minimum boiling margin | all points | liquid-filled points only | WELL-7: 97.1 vs 80.5 °C |
| Water level | "state the method" | line through the 3 shallowest liquid-filled P points, taken to 1.01325 bar; adds `water_level_extrapolation_m` | up to 18.4 m apart |
| Pivot | "state the rule" | highest quality tag, then thickest, then deepest | |
| Regime | list of names | ordered rule; adds `near_boiling` and `regime_note` columns | WELL-10 split between models |
| Isothermal section | "≥ 300 m" | consecutive points spanning ≥ 300 m within 5 °C, open hole only | 19, 29 and 22 points flagged |
| Hydrostatic windows | "split where 300 m allows" | top-down algorithm; a remnant under 300 m joins the window above | 24, 25 and 23 windows |

### 3. Natural-state confidence rubric (§4.1)

There is now a five-row rubric, with a tolerance of 3 bar or twice the RMS residual, whichever is larger. Verdict and confidence must agree in every file. **Why:** Haiku3 graded five wells lower than Sonnet3 and Opus3, which cut its total temperature weight by 13–14 %.

### 4. Geometry (§3, §7)

`model_surface_mRL` comes from the grid file's surface section, not from a layer top. A geometry self-check list was added to §7. **Why:** Haiku3 took the model surface from the layer top and got it wrong for 6 of 10 wells.

### 5. Weights (§5.5)

- The location rule is now ordered: above surface, near wellhead or top block, feed zone, behind casing, isothermal section, open hole. The first match wins. v3 listed the cases but gave no precedence.
- `weight_T` in `feedzones.csv` is defined.

### 6. Report hygiene (§6)

- **Claim tags:** exactly `[obs]`, `[interp]`, `[assumption]` or `[reference model]`.
- **Key conclusions:** each ends with an `**Action:**` line.
- **Template text:** none may be left unfilled.
- **Consistency:** values in the text must match the tables.
- **Conceptual model:** it now needs a numbered evidence table.
- **Figures:** each is embedded with a caption.

### 7. File format (§5.4)

Text files are UTF-8 with LF line endings, in a fixed row order. Feed zones are 1-based integers, and dates are ISO (`YYYY-MM-DD`). Only the report has a timestamp, in one line.

### 8. Standard figures (§8, new)

The v3 figure list (seven name patterns, little content specification) is replaced by seven specified figures and a shared house style. The layouts are adapted from O'Sullivan et al. (2025):

| ID | Figure | Adapted from (paper figure) | Main reason |
|---|---|---|---|
| F1 | Base map, with grid, faults, water and section lines | Figs 1, 3, 5 | one base layer under every plan view |
| F2 | Per-well profiles with a locator inset, plus group small multiples on shared axes | Figs 15–16 | wells comparable at a glance; the same layout as the calibration match plots |
| F3 | Elevation-band symbol maps (deviation from the field median, in fixed classes) | Fig. 17 | no interpolation between sparse wells; becomes the residual map later |
| F4 | Data section through the hottest wells, contoured only within the data | Fig. 14 (layout of Figs 9–12) | Opus3 contoured below its data |
| F5 | P–T state | kept from v3 | |
| F6 | Conceptual model, plan and section, fixed symbols, solid for observed and dashed for inferred, numbered | Figs 6, 18, 25 | Haiku3's sketch had no structure; makes `[obs]` and `[interp]` visible |
| F7 | Gauge check: transient against survey value, with production timing | Figs 20, 22 | makes the verdict that drives the weights checkable by eye |

**House style:**
- one colour per well;
- fixed temperature and elevation scales;
- fixed isotherms;
- data drawn as points, derived curves as lines;
- no colour where there are no data;
- reference-model content only in separate "[reference model]" figures in `Additional_Output/figures/`.

**Figures deliberately not adopted:**
- production and enthalpy histories, and time slices (out of scope);
- permeability, porosity, rock-type and upflow maps (reference model, not evidence);
- 3D alteration and resistivity views (no such data in the inputs).

### 9. Acceptance check (`Scripts/check_deliverables.py`)

| ID | Status | Change |
|---|---|---|
| all | changed | Paths moved to the two folders. The inventory and process log are read from `Additional_Output/`, and `acceptance_check.txt` is written there with LF endings |
| F1 | changed | Checks each file in its required folder. Figures are read from `Calibration_Input/figures/`. Adds `profile_panels_*.png`, and `gauge_check*.png` when any well has a gauge. Group panels no longer count towards the one-profile-per-well rule |
| F2 | new | `figures_index.csv` covers F1–F7, every listed file exists, `claim_type` values are valid; WARN for figures missing from the index |
| F3 | new | FAIL if a standard file is outside its folder; WARN for non-standard items in `Calibration_Input/` |
| S1 | changed | Adds columns for `T_P_at_elevation_bands.csv` and `figures_index.csv` |
| S3 | new | FAIL on non-UTF-8 files, non-integer or non-consecutive feed zones, or non-ISO dates; WARN on CRLF line endings |
| G4 | new | `model_surface_mRL` against the surface section at the wellhead column |
| D1 | new | Recomputes 14 `well_summary.csv` values per well, plus regime and pivot, from the targets using the §4.4 definitions |
| W3 | new | Weight = location × confidence factor under the ordered rule; `weight_P` is a quality factor × confidence; feed-zone weights; verdict and confidence agree across all four files |
| R1, R3 | changed | Search and compare files in both folders |
| R4 | new | No `{…}` or `nan` text in the report outside code |

Unchanged: I1, S2, G1–G3, T1, T2, W1, W2, R2, P1.

**Why CRLF is a WARN, not a FAIL.** This repository has `core.autocrlf=true`, so a correct LF file can come back as CRLF after a git checkout. CLAUDE.md still requires LF.

### 10. Verification of the new checker

The iteration-3 outputs were copied into the v4 folder layout in a scratch directory. They were not changed in the repository.

| Run | Result | What the new checks found |
|---|---|---|
| Opus3 | 5 FAIL | v4-only files missing (expected); D1: 5 values differ under the pinned definitions (water level for WELL-1, 4 and 10, WELL-2 deep gradient over a 98 m span, WELL-5 pivot); S3: non-ISO dates; G4 PASS |
| Haiku3 | 6 FAIL | **G4 FAIL for 6/10 wells**, matching the comparison report; D1: 23 differences |
| Sonnet3 | 6 FAIL | **R4 FAIL: five unfilled placeholders** in the report; the comparison report had found only one. D1: 9 differences |
| Opus3, patched to v4 | **0 FAIL, 1 WARN, 37 PASS** | Confirms that a compliant run passes. R1's WARN comes from diagnostic files the patch left out of the two folders. R3 tested with a stand-in script: PASS |

The D1 water level for WELL-1 (112.87 mRL) was checked by hand from the three shallowest points.

**Not yet tested:** no agent has run under v4. The pinned definitions may change some iteration-3 values, for example WELL-10's regime, and the figure specifications have not been tried. The figure content rules in §8 (data envelope, axes, symbology) cannot be checked from PNG files. They rely on the per-figure checklists and `figures_index.csv`.

### Other notes

- `Output/CLAUDE_v4_proposed.md` is the earlier draft, written without the folder split. The adopted `CLAUDE.md` supersedes it, and it can be deleted.
- The v3 checker is still in git history (commit `7975225`) if the earlier iterations need to be re-scored on the v3 yardstick.
