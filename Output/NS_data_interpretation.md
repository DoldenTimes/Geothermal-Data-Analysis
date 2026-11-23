# Natural-state data interpretation: TE14515 / Teaching_Model_3D

Scope: natural state only. Inputs: `Baseline/Teaching_Model_3D` (read-only). All numbers below come from the files in `Output/`, written by `Scripts/natural_state_analysis.py`. Units: mRL, bar, °C, kJ/kg, kg/s. Each claim is labelled **[Obs]** (observation in the data), **[Calc]** (derived by the script with a stated method) or **[Interp]** (interpretation; the observation is cited).

---

## 1. Key conclusions for calibration

1. **The data are labelled as synthetic, and this must be resolved first.** Every well's design record reads `Source: "Synthetic data"` (`Data_json/*.json`), the area is named "Synthetic" (`Groups_data.json`), and WELL-10's feed zone quality is "Synthetic". [Obs] If these surveys were generated from a model, calibrating to them tests consistency only. Reference-model results were not used as evidence (CLAUDE §2).

2. **No well can be confirmed as natural state.** The verdicts (`well_summary.csv`) are: possibly disturbed for WELL-2, 3, 4, 5, 8 and 9; probably natural for WELL-1, 6, 7 and 10. Every verdict is low confidence, so every target is weighted at 0.3× its base weight. [Obs] Producers WELL-2, 3 and 9 have production records from 01/01/2014. Injectors WELL-4, 5 and 8 have injection records from 01/01/2014. WELL-1, 6 and 7 have no flow or transient record. WELL-10's production starts 01/12/2019, after its survey. [Calc] `natural_state_verdict()` in the script.

3. **Surveys disagree with the transient gauges by up to about 15 bar.** [Obs] `pressure_consistency_check.csv`:
   - WELL-2 at −750 mRL: survey 98.2 bar; first transient 98.3 bar (01/01/2014), a 0.1 bar match. The transient at the survey date is 88.3 bar, so the survey is **9.9 bar higher**.
   - WELL-3 at −550 mRL: survey 80.7 bar; first transient 80.8 bar. At the survey date the transient is 71.6 bar, so the survey is **9.1 bar higher**.
   - WELL-9 at −675 mRL: survey 91.3 bar; first transient 89.2 bar (+2.1 bar). At the survey date the transient is 76.6 bar, so the survey is **14.7 bar higher**.
   - WELL-10 at −250 mRL: survey 56.7 bar; first transient 48.7 bar (01/01/2020), **8.0 bar higher**. The survey (2015) predates the transient record.
   Either the survey pressures or the transient series is offset. The first readings match the surveys, but they are also the first readings of the production record. [Interp]

4. **Pressure gradients in the upper hot sections are physically impossible for liquid water.** [Calc] `hydrostatic_check.csv`, intervals ≥200 m, all single-phase:
   - WELL-2, 210 → 10 mRL: implied density 1,101 kg/m³; 170 → −75: 1,211; 90 → −150: 1,304; 50 → −150: 1,346 kg/m³.
   - WELL-3, 164 → −75: 1,288 kg/m³; 50 → −150: 1,535 kg/m³.
   The IF97 liquid density for the same intervals is 865–964 kg/m³ (`rho_IF97_kg_m3` in the same file). Wellbore boiling would lower the density, not raise it, and boiling margins there are 30–105 K. [Interp] These pressures cannot be used as formation targets.

5. **A field-wide pressure-gradient excess is present even in cold wells.** [Calc] Median excess over single-phase intervals is +8.0% for the whole field (234 intervals, overlapping). Cold, uniform columns show it clearly: WELL-4 at 170 → −75 mRL (T 14.7 °C) gives 1,074 kg/m³ against about 1,000 for liquid, and 130–190 mRL intervals give +7 to +9%. WELL-8 at 170 → −75 mRL (T 23.3 °C) gives +7.6%, rising to +9.2% at 130 → −75 mRL. WELL-1 at 15 → −250 mRL (T 42.6 °C) gives +7.8%. [Interp] This is consistent with a datum, unit or gauge convention issue rather than overpressure. It should be checked against the logging records before any P target is used.

6. **The hottest points form a cluster, with the near-boiling column at WELL-10.** [Obs] Maximum temperatures: WELL-10 268.7 °C at −250 mRL; WELL-2 266.9 at −550; WELL-3 261.4 at −450; WELL-6 220.5 at −550; WELL-9 217.7 at −750; WELL-7 176.7 at −650; WELL-5 158.0 at −650. At 0 mRL the three hottest wells are WELL-10 219.1 °C, WELL-3 192.8 and WELL-2 171.8 (`T_P_at_common_elevations.csv`). WELL-10 has four survey points with boiling margin below 5 K (−75: 3.6; −150: 2.6; −250: 3.2; −262.5: 4.3), and its feed zone margin is 4.35 K (`well_summary.csv`, `feedzones.csv`).

7. **Each hot profile cools below a peak at about −450 to −750 mRL.** [Obs] `well_summary.csv` gives inversion drops from the interior maximum to the coolest point below: WELL-6 98.1 K (−550 → −1900); WELL-3 86.9 K (−450 → −1700); WELL-9 57.7 K (−750 → −2250); WELL-7 46.0 K (−650 → −1700); WELL-2 30.9 K (−550 → −1500); WELL-5 13.3 K (−650 → −1700). [Interp] This pattern is consistent with lateral outflow above a deeper cooled zone, or with cold downflow. The data alone do not say which. This is the main structural signal for the model.

8. **Two wells are cold and both are injectors.** [Obs] WELL-4 (T ≤ 29.0 °C, all survey points) and WELL-8 (T ≤ 93.1 °C; 12.8 °C at 176 mRL) both have injection records at about 100–110 kg/s with injectate enthalpy 334.9 kJ/kg (≈ 80 °C). WELL-4's survey temperatures are 12–29 °C, colder than its own injectate. A steady injector profile cannot be colder than the injectate, so WELL-4 cannot be read as natural state [Interp]. Cold ground at WELL-4 and WELL-8 therefore cannot be separated from injection cooling with these data.

9. **Deep pressures are connected across the field.** [Obs] At −1500 mRL, survey pressures are 162.0 bar (WELL-1), 163.9 (WELL-3, WELL-9), 164.0 (WELL-5), 164.5 (WELL-2), 164.8 (WELL-7), 165.1 (WELL-6) and 169.2 (WELL-8). Across the field this is a spread of 7 bar, or about ±3.5 bar about the mean (`T_P_at_common_elevations.csv`). At −1000 mRL the spread is 114.4 to 120.0 bar. [Interp] Within the caveat of conclusion 5, the deep pressure field is laterally connected.

10. **Several feed zones fall in layer gaps or on a column disagreement, and two wells have geometry datum problems.** [Obs] `data_integrity_check.csv`: 14 survey points fall in LAYERS gaps (for example 250–260 mRL and 970–980 mRL are in no layer). Three WELL-2 feed zones (Feed 1–3) assign a different column in `gTE14515_wells.json` (ajf) than the max-weight block at the midpoint (aji). WELL-1's track starts 105 m below the model surface, yet its top pressure of 10.8 bar implies a water level of 115.3 mRL (`well_geometry_check.csv`, `well_summary.csv`).

---

## 2. The data

**What exists.** Ten wells, each with a downhole temperature and pressure survey (`Data_json/WELL-n_data.json`). All surveys carry the same date, 29/07/2015 [Obs]. Feed zones: 18 in total (`feedzones.csv`). Transient pressure series: WELL-2 (−750 mRL, 12 values from 2014), WELL-3 (−550, 2014), WELL-9 (−675, 2014) and WELL-10 (−250, 2020). Flow records: production from 01/01/2014 for WELL-2, 3 and 9; production from 01/12/2019 for WELL-10; injection from 01/01/2014 for WELL-4, 5 and 8. WELL-1, 6 and 7 have none [Obs].

**Integrity checks** (`data_integrity_check.csv`) [Calc]:
- The `.dat` T and P files match the JSON at every point for all ten wells (0 difference). WELL-10's files are named `Wel10_`, which the script handles.
- The four transient series match their `.dat` files.
- Track end points match the `WELLS` section of `gTE14515.dat` to 0 m for all wells.
- All 18 feed zones lie within the well track and within a slotted (open) liner interval.
- Every survey point lies inside its track z-range.
- Layer check: 14 survey points fall in LAYERS gaps or outside layers (WELL-1 3, WELL-2 1, WELL-3 2, WELL-5 2, WELL-7 1, WELL-8 2, WELL-9 3).

**Natural-state timing verdicts** (`well_summary.csv`, `ns_verdict` and `ns_confidence`):

| Well | Verdict | Confidence | Basis |
|---|---|---|---|
| WELL-1 | probably natural state | low | no flow or transient record; survey date identical for all wells |
| WELL-2 | possibly disturbed | low | production from 01/01/2014; survey 9.9 bar above transient at survey date |
| WELL-3 | possibly disturbed | low | production from 01/01/2014; survey 9.1 bar above transient at survey date |
| WELL-4 | possibly disturbed | low | injection from 01/01/2014; survey colder than injectate |
| WELL-5 | possibly disturbed | low | injection from 01/01/2014 |
| WELL-6 | probably natural state | low | no record |
| WELL-7 | probably natural state | low | no record |
| WELL-8 | possibly disturbed | low | injection from 01/01/2014 |
| WELL-9 | possibly disturbed | low | production from 01/01/2014; survey 14.7 bar above transient at survey date |
| WELL-10 | probably natural state | low | survey predates production (01/12/2019); unresolved 8.0 bar mismatch with 2020 transient |

Confidence factors used in weights: high 1.0, medium 0.6, low 0.3. All wells are low, so every weight is multiplied by 0.3. If the engineer accepts the surveys as natural state, multiply the weights by 1/0.3 (see §6).

---

## 3. Temperature structure

Table from `well_summary.csv`. Depths are elevations in mRL. "Shallow" gradient = fitted over survey points at or above 0 mRL. "Deep" = fitted over the bottom 500 m of each survey.

| Well | Max T (°C) at mRL | T at TD (°C) | Inversion (mRL top → base, drop K) | T100 crossing (mRL) | T200 crossing (mRL) | Min boiling margin (K) at mRL | Water level (mRL) | Regime |
|---|---|---|---|---|---|---|---|---|
| WELL-1 | 118.6 at −1672 | 118.6 | none | −946 | none | 150.4 at 15 | 115.3 | conductive |
| WELL-2 | 266.9 at −550 | 253.7 (−2348) | −550 → −1500, 30.9 | 128 | −81 | 23.8 at −450 | 224.9 | upflow with lateral outflow |
| WELL-3 | 261.4 at −450 | 188.2 (−2887) | −450 → −1700, 86.9 | 104 | −14; −992 | 12.2 at −150 | 168.0 | upflow with lateral outflow |
| WELL-4 | 29.0 at −750 | 29.0 | none | none | none | 113.8 at 230 | 244.2 | cold |
| WELL-5 | 158.0 at −650 | 154.7 (−2595) | −650 → −1700, 13.3 | −160 | none | 123.6 at 93 | 132.8 | interior max, weak cooling |
| WELL-6 | 220.5 at −550 | 126.1 (−2359) | −550 → −1900, 98.1 | 13 | −329; −709 | 70.0 at −450 | 157.1 | upflow with lateral outflow |
| WELL-7 | 176.7 at −650 | 144.7 (−2762) | −650 → −1700, 46.0 | −114 | none | 80.5 at 170 | 189.9 | upflow with lateral outflow |
| WELL-8 | 93.1 at −2716 | 93.1 | none | none | none | 104.1 at 176 | 184.0 | cold |
| WELL-9 | 217.7 at −750 | 161.6 (−2427) | −750 → −2250, 57.7 | −84 | −496; −940 | 87.8 at −550 | 178.0 | upflow with lateral outflow |
| WELL-10 | 268.7 at −250 | 268.5 (−350) | none (near-boiling) | 165 | 46 | **2.65 at −150** (4 points < 5 K) | 209.9 | two-phase / near-boiling |

Shallow and deep gradients (K/km): WELL-1 140 / 16.5; WELL-2 771 / 31.4; WELL-3 1,057 / 19.7; WELL-4 13.7 / 19.6; WELL-5 328 / 21.2; WELL-6 591 / 7.1; WELL-7 260 / 25.0; WELL-8 80.6 / 30.0; WELL-9 247 / 9.0; WELL-10 949 / 285 (the bottom 500 m of WELL-10 is only 500 m long, so its "deep" value is the full profile).

Notes:
- Water level is estimated from the top survey point assuming a liquid column to atmospheric pressure. In WELL-1 it lies 100 m above the wellhead (see §7).
- Feed-zone boiling margins, from `feedzones.csv`: WELL-10 Feed 1 4.35 K; WELL-3 Feed 1 19.5 K; WELL-2 Feed 1 31.3 K; WELL-2 Feed 2 51.4 K; WELL-6 Feed 1 71.2 K; WELL-3 Feed 2 85.2 K. Every other feed zone exceeds 85 K.
- Isothermal sections (|ΔT| ≤ 1 K over ≥ 50 m, `well_summary.csv`): WELL-2 −1300:−1700; WELL-3 −450:−550, −1100:−1300, −1700:−1900; WELL-5 −550:−850, −1500:−1900; WELL-7 −1300:−1900; WELL-8 −750:−950 and −1100:−1500. WELL-4 is also near-isothermal over 230 → −25 mRL, at 12–17 °C.

**Wellbore reading (CLAUDE §3.2).** The upper section of WELL-2 (above the liner top at −50 mRL) and WELL-3 (above −200) sit in casing. The temperature there is closer to formation temperature than the open-hole interval. The rising shallow gradients in these two wells (771 and 1,057 K/km) are above the casing. [Interp] Two readings are possible: the upper column is strongly warm-upflow above the feed zones, or it is wellbore heating from below. The data cannot separate them.

---

## 4. Pressure structure

Interpolated only; blank = outside the survey range (`T_P_at_common_elevations.csv`, bar and °C).

| Well | T 0 mRL | P 0 mRL | T −250 | P −250 | T −500 | P −500 | T −750 | P −750 | T −1000 | P −1000 | T −1500 | P −1500 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| WELL-1 | 35.0 | 12.3 | 52.4 | 38.6 | 72.7 | 63.5 | 89.8 | 89.0 | 102.3 | 114.4 | 116.2 | 162.0 |
| WELL-2 | 171.8 | 24.1 | 238.2 | 55.0 | 266.1 | 77.5 | 257.7 | 98.2 | 241.1 | 120.0 | 236.0 | 164.5 |
| WELL-3 | 192.8 | 20.2 | 254.7 | 55.0 | 261.0 | 76.6 | 243.7 | 97.4 | 198.4 | 119.1 | 175.6 | 163.9 |
| WELL-4 | 15.3 | 24.5 | 19.4 | 51.1 | 24.1 | 76.3 | 29.0 | 100.5 | – | – | – | – |
| WELL-5 | 57.9 | 15.4 | 118.3 | 44.9 | 152.9 | 71.7 | 157.1 | 95.4 | 154.5 | 118.6 | 144.8 | 164.0 |
| WELL-6 | 106.3 | 18.6 | 186.1 | 49.9 | 218.8 | 75.7 | 193.6 | 97.0 | 153.6 | 119.1 | 126.7 | 165.1 |
| WELL-7 | 68.3 | 18.9 | 129.1 | 46.9 | 163.5 | 73.3 | 173.7 | 96.6 | 153.9 | 118.9 | 130.9 | 164.8 |
| WELL-8 | 26.7 | 19.8 | 48.0 | 46.1 | 61.9 | 71.3 | 67.8 | 95.9 | 66.0 | 119.8 | 64.5 | 169.2 |
| WELL-9 | 63.6 | 18.3 | 153.0 | 49.3 | 200.5 | 75.6 | 217.7 | 97.7 | 195.2 | 119.3 | 182.5 | 163.9 |
| WELL-10 | 219.1 | 26.3 | 268.7 | 56.7 | – | – | – | – | – | – | – | – |

At 200 mRL only WELL-2 (27.4 °C, 2.1 bar), WELL-4 (12.5 °C, 4.4 bar) and WELL-10 (34.9 °C, 1.4 bar) have data.

- **Observation [Obs]:** at 0 mRL pressures run from 12.3 bar (WELL-1) to 26.3 bar (WELL-10); the two cold injectors are at 19.8 bar (WELL-8) and 24.5 bar (WELL-4). At −750 mRL the spread is 95.4 to 100.5 bar. At −1500 mRL it is 162.0 to 169.2 bar.
- **Interpretation [Interp]:** in the hot, upper sections the pressures are not hydrostatic, as in conclusion 4. Below about −750 mRL the pressures agree with each other to a few bar, which makes them more usable as deep pressure targets, subject to the same datum question.

---

## 5. Conceptual model

Built from the data alone. Each feature is linked to its supporting observation.

Schematic, not to scale (positions from `well_geometry_check.csv` tracks; temperatures from the common-elevation table):

| Zone | Wells | Observation |
|---|---|---|
| Upflow core, −250 to −550 mRL | WELL-10, WELL-2, WELL-3 | peaks 268.7, 266.9, 261.4 °C; 0 mRL: 219, 172, 193 °C |
| Lateral outflow / cooling, −550 to −750 mRL and below | WELL-6, WELL-7, WELL-9, WELL-5 | peaks 220.5, 176.7, 217.7, 158.0 °C; cooling 13–98 K below peak |
| Cold ground (injectors) | WELL-4 (NE), WELL-8 (SW) | ≤ 29 °C and ≤ 93 °C; WELL-8 is 235 m from the Tarawera River |
| Deep connected pressure | all wells | 162–169 bar at −1500 mRL |

- **Upflow [Interp]:** the hottest column is centred on WELL-10 and WELL-2, near −250 to −550 mRL. [Obs] WELL-10 has 268.7 °C at −250, WELL-2 266.9 °C at −550, and WELL-3 261.4 °C at −450. At 0 mRL the same three wells are 219, 172 and 193 °C.
- **Near-boiling / two-phase margin [Obs]:** WELL-10 has four points below 5 K margin between −75 and −262.5 mRL, and its feed zone is at 4.35 K. [Interp] This is boiling-point-for-depth behaviour at the top of the upflow. Liquid pressures there remain above saturation.
- **Lateral outflow / cooling below the peaks [Obs + Interp]:** the peaks in WELL-2, 3, 6, 7 and 9 are followed by cooling of 31 to 98 K below. [Interp] The heat and fluid that leave the upflow must go laterally at roughly −550 to −750 mRL, toward WELL-6, 7 and 9 (all three show peaks at −550 to −750). The data do not show where this flow exits.
- **Cap [Obs]:** the 100 °C isotherm is at 104 mRL in WELL-3, 128 mRL in WELL-2 and 165 mRL in WELL-10. The 200 °C isotherm is at −14 mRL in WELL-3, −81 mRL in WELL-2 and 46 mRL in WELL-10. [Interp] The first 200 °C crossings in the three upflow wells lie between 46 and −81 mRL, so the hot column is close to the surface there. The data do not define a cap thickness, and cap thickness in the cold wells (WELL-4, WELL-8) is undefined because they never reach 100 °C.
- **Recharge / cold ground [Obs + Interp]:** WELL-4 (NE, 50 m from the nearest mapped fault) and WELL-8 (SW, 235 m from the Tarawera River) are cold, and both are injectors. [Interp] Cold ground may be injection-cooled, or it may be natural recharge along the fault or the river. The data cannot choose.
- **Faults and water [Obs]:** distance to the nearest of Faults 1–3 (`well_geometry_check.csv`): WELL-2 125 m, WELL-6 113 m, WELL-10 226 m, WELL-3 378 m, WELL-9 269 m, WELL-4 50 m. WELL-8 is 235 m from the Tarawera River and WELL-1 96 m from the unnamed river. No well is inside Lake Tarawera. [Interp] The hot wells lie within about 400 m of a mapped fault, so faults are a plausible control on the upflow. The data do not prove that link, and the fault that is nearest to each well is not recorded in the output.

Strength of claims: the observations (temperatures, boiling margins, inversion depths, distances) are firm. Their interpretation (upflow centre, outflow level, cold ground origin) is not, and should be tested against the survey dates and injection history before it is used as a model constraint.

---

## 6. Recommended calibration set-up

**Targets and weights** (`NS_calibration_targets.csv`, `feedzones.csv`). Base weights follow the CLAUDE §4 table, multiplied by the verdict factor (low = 0.3). The result is:

| Point type | Base | Current factor 0.3 | If surveys accepted (factor 1.0) |
|---|---|---|---|
| Feed zone | 1.0 | 0.30 | 1.00 |
| Open hole (slotted liner) | 0.8 | 0.24 | 0.80 |
| Behind casing | 0.5 | 0.15 | 0.50 |
| Within 50 m of top of track | 0.3 | 0.09 | 0.30 |
| Above model surface | 0 | 0 | 0 |

No point is above the model surface (`above_model_surface` is False for all 266 rows). Counts: 32 feed-zone points, 91 open-hole, 116 behind casing, 27 near top.

Recommendations:
1. **Temperature targets** are the most usable data. Give T at feed zones and open-hole points the highest weight. Use the hot wells first: WELL-10, 2, 3, 6, 7, 9.
2. **Pressure targets:** use P only at feed zones, and only after the datum question (§7.1) is answered. Do not use P above the feed zones in WELL-2 and WELL-3, or the top 50 m of any well.
3. **Injectors (WELL-4, 5, 8):** do not use their temperatures as natural-state targets. Either exclude them or model them as injection-cooled with their injection history. WELL-4 cannot be colder than its injectate in steady state.
4. **Sensitivity runs:** (a) weights as above, (b) survey factor 1.0 for the producers and 0.6 for "probably natural", (c) injectors excluded.
5. **Heat and mass in:** an upflow beneath the WELL-10/2/3 cluster, at about −250 to −550 mRL (§5). [Interp] Needed input: enough heat at the base of the cluster to hold a 200 °C column to about −1000 mRL in WELL-3.
6. **Heat and mass out:** lateral outflow at about −550 to −750 mRL toward WELL-6, 7 and 9, and a cooled zone below −1000 mRL. [Interp] The data give no constraint on where the outflow leaves the field. A boundary or sink is needed, but its location is not determined by these data.
7. **Permeability structure implied [Interp]:** the sharp cooling below the peaks (up to 98 K in WELL-6) suggests permeable lateral pathways at −550 to −750 mRL and low vertical conductivity below. This is an interpretation from temperature shape only. Rock properties are not evidence (CLAUDE §2).
8. **Boundary requirements [Interp]:** deep pressures at −1000 to −1500 mRL agree to about ±3.5 bar across the field (§4). This supports a laterally connected deep system, subject to the pressure caveats.

---

## 7. Data quality issues and open questions

Ranked by how much they affect calibration.

1. **Pressure gradients are physically impossible in the hot upper sections** (WELL-2 and WELL-3, implied 1,100 to 1,535 kg/m³; `hydrostatic_check.csv`). **Question for the data owner:** what gauge, datum and density convention produced these pressures? This affects every P target in the upper parts of the wells.
2. **Field-wide excess of about +8% in cold wells** (WELL-4 +7.4%, WELL-8 +4.0%, WELL-1 +5.5% median; `hydrostatic_check.csv`). Probably a datum or convention issue. **Question:** is the pressure referenced to a different elevation or density?
3. **Survey and transient series disagree by 8 to 15 bar** (`pressure_consistency_check.csv`). **Question:** what is the true survey date? All surveys carry 29/07/2015, which is implausible for ten wells. Were the transients logged on the same datum?
4. **"Synthetic" labelling** (`Well_Design.Source`, `Groups_data.json`, WELL-10 feed-zone quality). **Question:** are these real measurements, or generated from a model? The answer determines how the calibration should be read (§1.1).
5. **Injection-affected cold profiles** (WELL-4 colder than its injectate; WELL-8 at 13 °C near the top). **Question:** when was each injector's survey taken relative to injection, and what is the injection temperature history?
6. **WELL-1 wellhead datum** (track top 105 m below model surface; water level estimated at 115.3 mRL, 100 m above the wellhead). **Question:** is the wellhead elevation or the surface model right? A flowing or artesian condition would also explain it, but there is no discharge record.
7. **Layer assignment** (14 survey points in LAYERS gaps; WELL-2 feed-zone column disagreement). **Question for the model owner:** should the 10 m gaps between layers be closed, and which column is correct for the WELL-2 feed zones?
8. **Timing of WELL-10** (survey 2015, production from December 2019, transient 8 bar lower in 2020). **Question:** was WELL-10 shut in before 2019? If so, what explains the 8 bar difference?
9. **WELL-10 boiling** (four points and the feed zone below 5 K margin). **Question:** could boiling in the wellbore (rather than in the formation) explain the T profile here?
10. **Transient gauges are only at four elevations.** No transient in WELL-1, 4, 5, 6, 7, 8. Coverage is too sparse for a field-wide pressure check.

---

## 8. Output file index

| File | Contents |
|---|---|
| `Scripts/natural_state_analysis.py` | One script that regenerates every output: `python Scripts/natural_state_analysis.py` from the repo root |
| `Output/NS_data_interpretation.md` | This report |
| `Output/NS_calibration_targets.csv` | 266 survey points, all columns in the CLAUDE §4 specification |
| `Output/feedzones.csv` | 18 feed zones with survey T, P, margin, block and weight |
| `Output/T_P_at_layer_centres.csv` | T and P at each model layer centre inside the surveyed range (no extrapolation) |
| `Output/T_P_at_common_elevations.csv` | T and P at 200, 0, −250, −500, −750, −1000 and −1500 mRL, per well |
| `Output/well_summary.csv` | Per-well characteristics (CLAUDE §3.4) |
| `Output/pressure_consistency_check.csv` | Survey pressure vs first and time-matched transient |
| `Output/hydrostatic_check.csv` | Column-density tests, ≥ 200 m intervals, with field-wide median |
| `Output/well_geometry_check.csv` | Wellhead vs model surface, deviation, distances to faults, rivers and lake |
| `Output/data_integrity_check.csv` | .dat vs JSON cross-checks, transients, track vs WELLS section, feed-zone checks, layer checks |
| `Output/figures/well_map_grid_faults_water.png` | Well map with model vertices, faults and water bodies |
| `Output/figures/profile_WELL-n.png` (10 files) | T with Tsat(P), P, casing, liner and feed zones |
| `Output/figures/maps_temperature_by_elevation.png`, `maps_pressure_by_elevation.png` | Interpolated maps at seven elevations, inside the well hull only |
| `Output/figures/section_WELL-10_to_WELL-3.png` | Section through the two hottest wells (see note below) |
| `Output/figures/PT_state_vs_saturation.png` | P–T state against the IF97 saturation curve |
| `Output/run_console_log.txt` | Console output from the final run, including the IF97 validation |
| `Output/process_log.md` | Prompt, rules, steps and errors |

Note on the section figure: the only survey points inside the 400 m corridor are WELL-10 and WELL-3. The colour fill between them is a two-profile linear interpolation, not evidence of lateral structure. Points outside the corridor are drawn as outlines.
