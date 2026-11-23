# Comparing the Opus1, Sonnet1 and Haiku1 natural-state runs

**What was compared:** three natural-state analyses of `Baseline/Teaching_Model_3D`, each by a different model.

| Branch | Commit | Model |
|---|---|---|
| `origin/Opus1` | `3352c8a` | Opus 5.5 |
| `origin/Sonnet1` | `4c18c22` | Sonnet 5.5 |
| `origin/Haiku1` | `fe15a96` | Haiku 5.5 |

**Purpose:** to judge how useful each run's output is to an engineer or agent calibrating the natural-state model, and to record what the differences mean for future runs.

**Supporting files:**
- `Scripts/model_comparison_check.py`: an independent check of the key numbers against the raw data;
- `Output/model_comparison_check.txt`: its output.

---

## 1. Summary

| | Opus1 | Sonnet1 | Haiku1 |
|---|---|---|---|
| **Overall, for calibration use** | **Best** | Close second | Clearly weakest |
| Strongest point | A complete, usable package: every survey point mapped to a model block with a weight; a clear conceptual model; the best figures | The deepest physical reasoning about pressure, heat flux and enthalpy; the most candid caveats | Data checks: cross-checking the dat files against the JSON, track and liner checks, saturation-function validation, a detailed process log |
| Main weakness | Over-interprets the shallow pressure gradients as "overpressure beneath the cap" | No weights in the targets table; dense, hard-to-scan report; treats the reference model's boundary conditions as "data" | No mapping to the model grid at all; misses the main temperature inversions; would discard 6 of the 10 wells |

- **Facts:** all three runs read the raw data correctly. Every headline number I re-checked agrees with the Baseline data (§3).
- **Where they differ:** interpretation, and what they deliver to the calibrator.
- **Ranking:** Opus produced the output a calibrating agent can use straight away. Sonnet produced the best reasoning but packaged it less usefully. Haiku did the data checks carefully, but the interpretation and the targets it produced are not usable for calibration as delivered.

**Caveat:** this is one run per model. Differences between models are confounded with ordinary run-to-run variation (§8).

## 2. Was the comparison fair?

| Condition | Check | Result |
|---|---|---|
| Same instructions | md5 of `CLAUDE.md` at each run's commit | Identical in all three |
| Same input data | `git diff` of `Baseline/` between the three commits | No differences |
| Same prompt | The process logs | "Analyse the data" in every run |
| Independence | The process logs | Sonnet and Haiku saw that earlier outputs had been deleted (in `git status`). Both say they did not restore or use them |
| Environment | The process logs, and re-running the scripts | Python with numpy and matplotlib only, in every run |

One caveat: Opus ran while the "Fixed well geometry" change (`7b0992b`) was still uncommitted in the working tree. The data it read were the same.

## 3. Accuracy: what they agree on, checked independently

`Scripts/model_comparison_check.py` recomputes the key quantities directly from `Baseline/`. All three runs match it wherever they report the same quantity:

| Quantity | Independent value | Opus1 | Sonnet1 | Haiku1 |
|---|---|---|---|---|
| Tmax, WELL-10 / WELL-2 / WELL-3 (°C) | 268.7 / 266.9 / 261.4 | ✓ | ✓ | ✓ (as T at fixed elevations) |
| Largest inversions: WELL-6 / WELL-3 / WELL-9 (°C) | 98.1 / 86.9 / 57.7 | ✓ | ✓ | ✗ not reported |
| Minimum boiling margin, WELL-10 (°C @ mRL) | 2.6 @ −150 | ✓ ("within 3 °C") | ✓ | ✓ |
| P at −750 mRL, WELL-1 vs WELL-4 (bar) | 89.0 vs 100.5 | ✓ | ✓ | ✓ |
| Survey vs first transient pressure, WELL-2 / WELL-3 / WELL-9 (bar) | 98.2 vs 98.3; 80.7 vs 80.8; 91.25 vs 89.2 | ✓ | ✓ | ✗ not done |
| WELL-10 transient (2020) vs survey (bar) | 48.7 vs 56.7 | Listed only | ✓ ("depleted ~8 bar") | ✗ |
| Wellhead minus model surface, WELL-1 / WELL-8 / WELL-7 (m) | −105 / −84 / −50 | ✓ | ✓ | ✗ not done |
| Fault 1 / Fault 3 intersection | ≈ (2819778, 6335556), about 260 m from WELL-2 | ✓ (qualitative) | ✓ (2819780, 6335480; ~250 m) | ✗ |

No run made an arithmetic or data-reading error that I could find. Every gap in the table is something a run **did not do**, not something it got wrong.

## 4. Interpretation: where the runs differed

| Question | Opus1 | Sonnet1 | Haiku1 | Assessment |
|---|---|---|---|---|
| **Are the 2015 surveys natural state?** | Yes. Tested against the first transient pressures | Yes. Same test, and also noticed that WELL-10's transient series is already depleted | "Possibly disturbed" for 6 of the 10 wells. Recommended excluding the producers and injectors **without running the available test** | Opus and Sonnet are right. Haiku's caution is reasonable in principle but untested, and acting on it would discard most of the hot-zone data |
| **Temperature inversions** | All found, with sizes. Interpreted as lateral outflow over cooler fluid | All found, with sizes. Same interpretation, and called them a strong constraint on the outflow geometry | Reported only gradients of a few °C/100 m *inside feed zones*. **Missed the 46–98 °C inversions** in WELL-3, 6, 7 and 9 | Haiku's main interpretive failure. The inversions are the clearest outflow signal in the dataset |
| **WELL-9 regime** | Outflow (58 °C reversal) | Outflow / cooling | "Hot upflow" | Haiku is wrong: WELL-9 has a 58 °C inversion below its maximum |
| **WELL-4** (12–29 °C to −750 mRL) | Strong cold downflow / recharge | Background conductive (~1.7 °C/100 m) | Cold (injector) | Opus and Sonnet are both defensible. A low gradient next to a hot field fits recharge, but neither run proved it |
| **Gradient excess over hydrostatic** | WELL-2 ≈1.10, taken as upflow. Top 300 m at 1.2–1.8, taken as "overpressure beneath the cap" | Found a **field-wide 2–3 % excess even in cold wells**, and took only the excess above that as meaningful. Concluded that deep pressure is set by a laterally connected source, not by each well's own column | Up to +15 % (WELL-2) and +38 % (WELL-10), from Kell density used above its valid range. Called WELL-10 "not physically possible" and blamed the data | See §4.1 |
| **Enthalpy** | WELL-10 steam enthalpy taken as confirming boiling | WELL-2 first enthalpy (1302 kJ/kg) is 180–260 kJ/kg above the surveyed liquid: needs an explanation | Not used | Sonnet's point is the most useful one for calibration, although strictly it relies on production-period data |
| **Heat flux** | Lower background flux on the margins | ~30–75 mW/m² from the deep gradients, with conductivity **explicitly stated as an assumption**. The reference model's uniform 0.1 W/m² looks high | Not addressed | Sonnet handled the assumption best |

### 4.1 The shallow pressure gradient: an issue none of the runs fully identified

Between the top of the survey and −250 mRL, the measured pressure implies an average fluid density of **1034–1320 kg/m³ in every well**, including the cold wells WELL-4 and WELL-8 (`model_comparison_check.txt`). A static column of water inside a well cannot be denser than cold water.

This points to the profiles being **formation pressures**, perhaps sampled from a model's block pressures, rather than a real static wellbore column.

How each run handled it:
- **Opus** read it as physical overpressure beneath the cap.
- **Sonnet** read it as pressure imposed by a deep, laterally connected source. This is closest, and Sonnet also flagged the cold-well excess as possibly a convention issue.
- **Haiku** flagged WELL-10 as physically impossible and suspected the data. This is the right instinct, but Haiku based it on a density formula used outside its valid range and did not extend it to the other wells.

**For calibration:** do not build cap overpressure into the model on the strength of the shallow pressure gradients alone. `CLAUDE.md` §3.2 now tells runs to treat measured pressure as formation pressure only at feed zones.

## 5. Deliverables for a calibrating agent

| Deliverable | Opus1 | Sonnet1 | Haiku1 |
|---|---|---|---|
| Survey points mapped to model column, layer and block | ✓ 266 points, following deviated tracks | ✓ 250 block intervals, via `gTE14515_wells.json` | ✗ none |
| Suggested weights | ✓ (1.0 / 0.8 / 0.5 / 0.3 / 0) | ✗ described in the text only | ✗ |
| Tsat and boiling margin in the targets table | ✓ | In `hydrostatic_check.csv` | In `downhole_profiles_long.csv` |
| T and P at layer centres | ✓ 34 layers | ✓ inside the targets table | ✗ (fixed elevations only) |
| Feed-zone targets with blocks | ✓ | ✓, plus the liquid-enthalpy equivalent | Without blocks |
| Wellhead vs model surface | ✓ | ✓, plus distances to faults, rivers and lake | ✗ |
| Conceptual-model sketch | ✓ ASCII diagram | Text only | ✗ |
| Cross-sections | ✓ 2 | ✓ 2 | ✗ |
| Recommended permeability, heat and boundary set-up | ✓ explicit | ✓ (spread across the report) | ✗ (targets list only) |
| Data cross-checks (dat vs JSON, track, liner) | Partial | Partial | ✓ the most thorough |

**Figures.** I reviewed the T and P profile figure from each run:
- **Opus:** the most informative. It shows the Tsat curves and marks feed zones well by well.
- **Sonnet:** adds a panel of pressure excess, but it shares the 0–300 bar x-axis with the pressure panel, so the data are squashed against the left edge and unreadable.
- **Haiku:** shades the feed zones as field-wide bands, so you cannot tell which well each one belongs to.

## 6. Reproducibility and code

I re-ran each run's script on a clean copy of its branch, in a scratch folder.

| | Opus1 | Sonnet1 | Haiku1 |
|---|---|---|---|
| Runs without error | ✓ | ✓ (one RuntimeWarning) | ✓ |
| Tables reproduced exactly | ✓ all 6 CSVs | ✓ all 7 CSVs | ✓ all tables. The console log differs only in the output path |
| Script length | 527 lines | 479 lines | 572 lines |
| Run time | 6.9 s | 3.0 s | 2.2 s |
| Saturation (Tsat) | IF97 Region 4 | IF97 Region 4 | IF97 Region 4, validated in the run |
| Liquid density | Saturated-liquid table | Saturated-liquid table plus a compressibility term | Kell polynomial, **used beyond its 150 °C limit** |

All three are fully reproducible. Note that the interpretation reports are written by the agent, not the script, so they cannot be regenerated.

## 7. Compliance with the instructions (original `CLAUDE.md`)

| Rule | Opus1 | Sonnet1 | Haiku1 |
|---|---|---|---|
| `Baseline/` unedited; scripts in `Scripts/`; outputs in `Output/` | ✓ | ✓ | ✓ |
| README ignored | ✓ | ✓ | ✓ |
| Prompt and steps logged | ✓ (short) | ✓ (includes limitations) | ✓ (most detailed, including the bugs it fixed) |
| Reference model not used as data | ✓ | ~ Listed the boundary conditions and heat flux under "what data exist", though it compared them critically | ✓, and arguably too strict: it excluded the grid geometry, which the rules allowed |
| Natural state only | ~ Used WELL-10's production enthalpy as support | ~ Used initial enthalpy and flow rates | ✓ strict |

The original `CLAUDE.md` did not say whether "geometry" meant the model grid. Haiku read it narrowly, as well geometry only, and so never mapped anything to model blocks. This ambiguity in the instructions explains much of the gap in Haiku's deliverables.

## 8. Limitations of this comparison

- **One run per model**, so run-to-run variation cannot be separated from differences between models.
- **No cost data.** Token use and elapsed session time were not recorded, so cost-effectiveness is not assessed.
- **The runs were sequential**, not simultaneous. Later runs could see that earlier outputs existed (they report not using them).
- **My own density check is approximate.** It uses saturated-liquid density with no compressibility. That is enough to test claims at the ±5 % level, not finer.
- **Judgement calls.** "Usefulness for calibration" is my judgement against the purpose stated in `CLAUDE.md`. It has not been checked by actually calibrating the model with each run's outputs.

## 9. Lessons and how `CLAUDE.md` now handles them

| Problem seen in the runs | Run(s) | Rule now in `CLAUDE.md` |
|---|---|---|
| Natural-state status assumed rather than tested | Haiku | §3.1: compare survey pressures with the first transient value, and give a per-well verdict with confidence |
| Unclear whether model geometry may be used | Haiku | §2 table: the `gTE14515.dat` grid sections and `_wells.json` are allowed as geometry |
| Reference model's boundary conditions treated as data | Sonnet | §2: listed as "reference only, never evidence" |
| Kell density used above 150 °C; liquid-density test applied near boiling | Haiku | §3.3: IF97 Region 1 required, and the test is limited to single-phase liquid |
| Shallow wellbore pressure read as reservoir overpressure | Opus (and Haiku's WELL-10 reading) | §3.2: pressure represents the formation only at feed zones and pivot points |
| File names and columns differ between runs | All | §4: fixed file names and a fixed column list for the targets table |
| Targets without weights | Sonnet, Haiku | §4: weights and natural-state confidence required |

### Recommendation on which model to use

- **For production runs of this workflow, use Opus.** Its output was the most complete and the most directly usable.
- **Sonnet is an acceptable alternative** now that `CLAUDE.md` §4 fixes the deliverables: its weakness was packaging, not reasoning.
- **Haiku is suitable only for data-QC steps** (format cross-checks, validation), not for interpretation.

Re-running all three under the revised `CLAUDE.md`, ideally more than once each, would show whether the gaps close once the instructions are tighter.
