# Iteration 4: how Haiku, Sonnet and Opus performed, and what changed from iteration 3

*Compares the three iteration-4 runs (branches `Haiku4`, `Sonnet4`, `Opus4`, all 26 Nov 2026) with each other and with iteration 3 (`Haiku3`, `Sonnet3`, `Opus3`, 25 Nov). Iteration-3 figures are quoted from `Logs/V3/model_comparison_3/` and `Logs/V4/CHANGELOG.md`. Iteration-4 figures were computed directly from each branch's `Output/` and `Scripts/` (extracted to a scratch folder; nothing in the repository was changed). The three run folders were not re-run here. The acceptance results are those each run wrote to its own `acceptance_check.txt`.*

*Placement note: CLAUDE.md says nothing goes loose in `Output/`. You asked for this report in the output folder, so it sits at `Output/`. Move it to `Output/Additional_Output/` if you want to follow that rule strictly.*

---

## The short version

- **The three models now give the same answer on every number that can be compared.** All 160 per-well summary values (16 quantities × 10 wells) agree across the three runs, with a maximum spread of 0.05 on any quantity. In iteration 3, water level agreed on 0/10 wells, the shallow gradient on 5/10, the deep gradient on 1/6 and natural-state confidence on 5/10. Regime, verdict and confidence are identical on all 10 wells. v4's pinned definitions did what they were meant to do.
- **All three pass the acceptance check with no FAIL.** Opus and Sonnet score 38 PASS and 0 WARN. Haiku scores 37 PASS and 1 WARN (a hydrostatic band label). The checker is now stricter (38 checks, up from 26), and every run also passes the clean re-run (R3).
- **Haiku has one genuine error that the checker missed, and the differences left are in interpretation and figures.** Haiku gave three deep survey points the "within 30 m of wellhead" weight (0.225 instead of 0.375 or 0.6). Sonnet and Opus agree on all 266 point weights. Haiku also wrote the thinnest conceptual model: it says the cap is "not resolved", where Sonnet and Opus measure one.
- **One new ambiguity appeared in the spec: hydrostatic windows that straddle the casing shoe.** The three models handled it three different ways. The "cased" cold-well ratio is 1.18, 1.17 or 1.47 depending on the run, so it is not comparable between runs.

---

## 1. What changed from iteration 3 to iteration 4

Source: `Logs/V4/CHANGELOG.md` and the CLAUDE.md diff.

| | Iteration 3 | Iteration 4 |
|---|---|---|
| CLAUDE.md | 220 lines | 413 lines at first, then a small edit (see §6 caveat) |
| Checker | 26 check lines, 15 IDs | 38 check lines, 22 IDs. New: F2, F3, S3, G4 (surface), D1 (recomputes 14 summary values per well), W3, R4 |
| Output layout | 25–40 files in one folder | Two folders: `Calibration_Input/` (5 tables, report, figures index, figures) and `Additional_Output/` |
| Definitions | Several left to the model | Pinned: liquid-filled point, shallow and deep gradient fit sets, isotherm crossing, boiling margin population, water level, pivot, regime rule, isothermal section, hydrostatic window algorithm |
| NS confidence | Loose guidance | Five-row rubric, tolerance 3 bar or twice the RMS residual |
| Weights | Cases listed | Ordered rule (first match wins), `weight_T` defined |
| Geometry | Model surface not defined precisely | Taken from the grid file's SURFA section, with a self-check |
| Figures | Seven name patterns | Seven specified figures (F1–F7), a house style, per-figure checklists, data-envelope colour fill, and a required interpolated companion set |
| Report | Tags and structure | Exactly four claim tags, an `**Action:**` line on each key conclusion, no template text, a numbered conceptual-model table |

Every iteration-3 error that a v4 rule targeted is gone:

| Iteration-3 problem | v4 rule | Iteration 4 |
|---|---|---|
| Haiku's model surface wrong for 6/10 wells | G4 check and SURFA definition | 10/10 correct in all three runs |
| Sonnet's unfilled placeholder `{len(irregular)}` | R4 check | None in any run |
| Water level up to 18.4 m apart | pinned method | Spread ≤ 0.05 m |
| Shallow gradient up to 180 °C/km apart | pinned fit set | Spread ≤ 0.03 °C/km |
| Deep gradient up to 13 °C/km apart | pinned fit set | Spread ≤ 0.04 °C/km |
| WELL-10 regime split | ordered rule plus `near_boiling` | All three give "upflow" |
| Haiku graded 5 wells lower than the others | five-row rubric | All 10 verdicts and confidences identical |
| 19, 29 and 22 isothermal points | pinned definition | Weights agree (see §3) |

---

## 2. Compliance

| | Haiku4 | Sonnet4 | Opus4 |
|---|---|---|---|
| Checks (PASS / WARN / FAIL) | 37 / 1 / 0 | 38 / 0 / 0 | 38 / 0 / 0 |
| WARN | S2: a `straddles shoe` band label | none | none |
| D1 (14 values × 10 wells + regime and pivot) | all 160 match | all 160 match | all 160 match |
| G1 / G2 / G3 / G4 | pass (266 points, 18 feed zones, 30 layer centres, 10 surfaces) | same | same |
| R1: report numbers not found in any table | 0/389 | 25/665 (4%) | 0/566 |
| R3 clean re-run | pass | pass | pass |

- All three reproduce every table and the report on a clean re-run.
- **R1 is a PASS for Sonnet but not a clean one.** The checker tolerates some unmatched numbers. Sonnet's report quotes 25 numbers (examples: 72.7, 266.1, 241.1, 218.8) that the checker could not find in any table. Sonnet keeps a `report_numbers.csv`, so these may be rounded or derived values. CLAUDE.md says every number must come from an output file, so check that these are traceable.
- **The checker did not catch Haiku's weight error** (§3). W3 passes because the check cannot tell whether a point is truly "within 30 m or in the top block".

---

## 3. Correctness

Compared point by point over the 266 survey points:

| Pair | T values agree | Blocks agree | Temperature weight agrees | Pressure weight agrees |
|---|---|---|---|---|
| Opus vs Sonnet | 266 | 266 | **266** | 266 |
| Opus vs Haiku | 266 | 266 | 263 | 266 |
| Sonnet vs Haiku | 266 | 266 | 263 | 266 |

**Haiku's weight error.** Three points are labelled "within 30 m of wellhead or top block of column" although they are deep: WELL-8 at −1300 mRL (layer 63), and WELL-9 at −650 mRL (layer 58) and −1500 mRL (layer 64). Each gets 0.225 where Sonnet and Opus give 0.375 (WELL-8) or 0.6 (WELL-9). The cause is in `Scripts/natural_state_analysis.py` around line 542: the rule tests `r._top_col or r._near_wh`, and one of those flags fires for these points. I did not trace which. The effect is small (total temperature weight 124.2 against 125.0), but it is a real misapplication of §5.5 that the acceptance check does not detect.

Other correctness points:

- Total pressure weight is identical in all three runs (20.65, at 32 feed-zone points). Feed-zone weights match (`weight_P` total 11.4, `weight_T` total 13.75). The tags found (HIGH, MEDIUM, LOW and a non-grade "Synthetic") were mapped the same way.
- All three correctly used IF97 and passed the 0.1 K Tsat check (worst error: Haiku 3.9e-7 K, Sonnet 0.004 K, Opus 0.004 K).
- All three gave WELL-1 "possibly disturbed, low" (feed zone 3.3 bar below trend, outside the 3 bar tolerance), WELL-5 "probably natural state at best, low" (injector, formation hotter than the injectate), and WELL-2 and WELL-3 "natural state, high".
- **Sonnet caught and fixed a method error of its own.** Its first injector rule used a 15 °C window around the injectate temperature, which flagged WELL-8 as injection-affected although its isothermal section is 10–14 °C colder than the injectate. It tightened the window to 5 °C.
- **Opus** worked through a parsing error, a `KeyError`, a missed dot-file in the inventory and a re-run failure (R3), and fixed all of them before finishing. It also changed its data-envelope logic to follow §8.2.

---

## 4. Do the models agree with each other?

| Quantity | Iter 3: wells where all models agree | Iter 4 |
|---|---|---|
| Tmax, T at TD, inversion, pivot pressure | 10/10 (inversion 7/7) | 10/10 |
| 100 °C / 200 °C isotherm | 8/8, 5/5 | 10/10 (spread ≤ 0.05 m) |
| Minimum boiling margin | 9/10 (16.6 °C gap) | 10/10 |
| Shallow gradient | 5/10 (up to 180 °C/km apart) | 10/10 (≤ 0.03) |
| Deep gradient | 1/6 (up to 13 °C/km apart) | 10/10 (≤ 0.04) |
| Water level | 0/10 (up to 18.4 m apart) | 10/10 (≤ 0.05 m) |
| Regime | 9/10 | 10/10 |
| NS confidence | 5/10 | 10/10 |

The full agreement means that the differences between the three runs are no longer in the numbers. They are in the hydrostatic bands (§5), the interpretation (§6) and the figures (§7).

---

## 5. Where the spec is still ambiguous: hydrostatic bands

The window algorithm merges a remnant under 300 m into the window above. That produces windows whose top is above the casing shoe and whose bottom is below it. §4.3 defines only "cased" (bottom at or above the shoe) and "open" (top below the shoe), so these windows fit neither band.

| | Haiku4 | Sonnet4 | Opus4 |
|---|---|---|---|
| Windows | 25 | 26 | 26 |
| Cased / open | 7 / 16 | 9 / 17 | 2 / 15 |
| Straddling windows | 2, labelled `straddles shoe` (triggers WARN S2) | 0 (labelled by the larger share of length, flagged in a `straddles_shoe` column) | 9, left with a blank `band` |
| Cold-well cased ratio (median of cased windows) | 1.176 | 1.174 | 1.473 |
| Open ratio (median) | 1.020 | 1.020 | 1.020 |

- **The open-band result is the same everywhere** (1.020, roughly a 2% field-wide offset).
- **The cased result is not comparable.** Opus's 1.47 rests on only 2 windows (560 m in total). Its process log states that no cold-well window lies wholly behind casing, so its cased-band cold-well component is undefined. Haiku's and Sonnet's values rest on different window sets.
- **Each handling is a departure.** All three recorded it in the process log. The cleanest one for the checker is Sonnet's (no foreign label, so no WARN), but labelling by "larger share" is Sonnet's own rule.
- **Action for the next CLAUDE.md:** state how to band a window that straddles the shoe (for example, split it at the shoe, or add a third band).

---

## 6. Interpretation and report

| | Haiku4 | Sonnet4 | Opus4 |
|---|---|---|---|
| Report words | 3,881 | 5,653 | 5,733 |
| Key conclusions / `**Action:**` lines | 9 / 19 | 8 / 8 | 9 / 9 |
| Claim tags `[obs]` / `[interp]` / `[assumption]` / `[reference model]` | 51 / 22 / 3 / 1 | 61 / 29 / 5 / 1 | 69 / 49 / 4 / 9 |
| Process-log words | 2,044 | 1,352 | 1,751 |
| Script lines | 2,430 | 2,769 | 3,149 |
| Files in `Additional_Output/` (CSV) | 13 | 15 | 15 |
| Figures in `Calibration_Input/` / `Additional_Output/` | 19 / 4 | 21 / 6 | 23 / 10 |

**What the three agree on.** The upflow at WELL-10 and WELL-2, the lateral outflow layer at about −450 to −750 mRL, the cold margins (WELL-4 at 29.0 °C and WELL-8 at 93.1 °C), WELL-10 as near-boiling but not two-phase (minimum margin 2.65 °C), and a field pressure trend P = 32.76 − 0.0877 z through the gauge-verified feed zones of WELL-2 and WELL-3. All three treat a gauge reading on the day production starts as pre-production, and record that as an assumption. This drives the "high" verdicts for WELL-2 and WELL-3.

**What differs.**

- **The cap.** Sonnet and Opus identify a conductive cap from shallow gradients of about 336–692 °C/km, with a labelled top and base in the section. **Haiku concludes the cap is "not resolved by the survey data"** and draws it as absent. The data support a cap (steep shallow gradients above the 200 °C level), so Haiku's conclusion is the weakest of the three.
- **Beneath the outflow.** Only Opus identifies cooler, nearly isothermal water beneath the outflow, as a separate element with its own arrow, and draws it in its conceptual model.
- **Pressure.** Haiku reads the steep cased-section pressure gradients as a conflict with the temperature profiles. Opus finds shallow gradients steeper than any water column (WELL-3, ratio about 1.5) and treats them as something to check. Sonnet focuses on deep pressures agreeing across the field. These are different emphases, and the cased-band ambiguity in §5 feeds into them.
- **Extra analysis.** Opus drew three section lines (A through the two hottest wells, plus B and C chosen to capture the most wells within 500 m) and made two extra conceptual sections. Sonnet drew two (A and B). Haiku drew one (A).
- **Haiku also noted** that wellheads sit 13–105 m below the model surface in 8 wells. It made this a numbered conclusion.

**Report quality.** Haiku's report is about two-thirds the length and has the fewest interpretation tags. Sonnet's and Opus's are close in length. Opus has by far the most `[interp]` tags (49) and `[reference model]` mentions (9), which reflects its fuller conceptual model.

---

## 7. Figures

I looked at F4 (`section_A.png`) and F6 (`conceptual_model.png`) for each run. I did not review the other figures, so this is a partial check against the §8 checklists.

| | Haiku4 | Sonnet4 | Opus4 |
|---|---|---|---|
| F4 colour fill | Yes, but only 0–350 mRL for WELL-8, WELL-10 and WELL-2, since WELL-10 is shallow. The deep parts of WELL-8 and WELL-2 are hatched | Yes, wells W8, W6, W3, W10, W2 (offsets up to 810 m, within 1000 m). Filled down to −2350 mRL where neighbours have deep data | Yes, W8 and W2, filled to −2350 mRL, with a locator inset |
| Offsets printed, fault and river ticks, no-data hatching | Yes | Yes | Yes |
| Vertical exaggeration | ×0.4 | ×1.0 | ×1.5 |
| F4 layout issues | Large empty area below the data; feed-zone and shoe labels small and clipped | Clean, but river labels overlap at the top left | Clean, with the end-distance labels and locator |
| F6 section | Small and cropped: the colour fill only reaches −350 mRL; arrows are tiny; the cap is shown as "absent" | Large and readable; red, orange and blue arrows with numbers; cap top and base labelled; pressure connections labelled | Large and readable; eight numbered elements, including the cooler-water element; cap tops and bases labelled; ΔP connections drawn on both panels |
| F6 plan view | Titles overlap; arrow numbers readable | Clear; cold-region circles and ΔP labels | Clear; three section lines, ΔP network, the hot-zone footprint shaded |
| F6 numbering vs table | 7 elements; legend says cap absent | 6 elements with a key beneath the panel | 8 elements with a legend |

- **Haiku's figures have the least information.** In both F4 and F6 the section is mostly no-data, partly because it chose a section line that includes the shallow WELL-10 as the only neighbour of WELL-8 above −350 mRL. It does satisfy the checker, which cannot test figure content.
- **Sonnet's and Opus's figures follow the §8 checklists closely.** Opus is the most complete, with more section lines and the richest F6.
- **Interpolated companions** are present in all three (Haiku 4 files, Sonnet 6, Opus 10 including two extra conceptual sections), as §8.4 requires. I did not check that every title starts with "[interpolated, not data]".

---

## 8. Effort

| | Iter 3 (H / S / O) | Iter 4 (H / S / O) |
|---|---|---|
| Script lines | 1,362 / 2,456 / 2,301 | 2,430 / 2,769 / 3,149 |
| Report words | 4,247 / 7,792 / 7,390 | 3,881 / 5,653 / 5,733 |
| Process-log words | 2,288 / 2,249 / 2,100 | 2,044 / 1,352 / 1,751 |
| Figures (all folders) | 22 / 23 / 17 | 23 / 27 / 33 |

- **Scripts grew** (Haiku by 78%, Sonnet by 13%, Opus by 37%). Haiku's grew most as it had to build the figure standard and the interpolated companions.
- **Reports got shorter** for all three, by 9% (Haiku), 27% (Sonnet) and 22% (Opus). The calibrator now needs only the single `Calibration_Input/` folder, so detail moved to `Additional_Output/`. This is the intended effect.
- **Figures went up for Sonnet and Opus** and stayed level for Haiku.

---

## 9. Caveats

1. **The CLAUDE.md edit is a confound for Opus.** The commit order is Opus4 (14:42), then "Minor edits to CLAUDE.md" (14:43), then Haiku4 and Sonnet4. That edit added the colour-filled envelope rule and the required interpolated companion set (§8.4). Opus's final commit is a minute before the edit, so it may have run under the earlier v4 text. Its outputs nonetheless include the interpolated companions and a filled envelope, so I cannot tell which version it followed. I treat the runs as comparable, but the difference in figure scope may partly reflect the wording each model saw.
2. **One run per model per iteration.** With no repeats, differences such as Haiku's weight error might not recur.
3. **Haiku's weight error is traced only as far as the rule that applied it.** I did not find which flag was wrong.
4. **Figures were reviewed only for F4 and F6** (§7).
5. **Iteration-3 numbers come from the earlier comparison and used the v3 checker.** Iteration 4 used the stricter v4 checker, so the pass counts are not like for like. The agreement and spread numbers are like for like.

---

## 10. What to take from this

1. **v4 worked.** Pinning the definitions removed the number-level disagreement completely. Pass rates stayed at "no FAIL" under a checker that is nearly 50% bigger.
2. **Next fix for CLAUDE.md:** define how to band a hydrostatic window that straddles the casing shoe. It is now the only place where the models' numbers diverge, and it directly affects the cold-well component.
3. **Next fix for the checker:**
   - test the "within 30 m of wellhead or top block" rule against the track, so Haiku's three mis-weighted points would be flagged;
   - tighten R1, or require Sonnet's `report_numbers.csv`-style traceability for every number;
   - test the content of figures where possible (envelope and hatch checks cannot be read from a PNG).
4. **For calibration,** use the Opus or Sonnet targets. They agree at all 266 points, and Opus adds the cap, cool-water and pressure-network detail. Haiku's targets are within 0.8 weight units of them in total, but its conceptual model omits the cap.
5. **Choice of model is now about interpretation and presentation, not compliance.** Opus ≥ Sonnet > Haiku on the conceptual model and figures. All three are acceptable on the tables.

---

## Files used

| Source | What it gave |
|---|---|
| `origin/Haiku4`, `origin/Sonnet4`, `origin/Opus4` (`Output/`, `Scripts/`) | All iteration-4 numbers and figures |
| `Logs/V4/CHANGELOG.md`, `git diff fa43d29 fd93f74 -- CLAUDE.md` | What changed in v4 and in the later edit |
| `Logs/V3/model_comparison_3/model_comparison_across_iterations.md`, `.../model_agreement_by_iteration.csv`, `.../run_overview.csv` | Iteration-3 figures |
