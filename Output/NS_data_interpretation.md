# Natural-state data interpretation: calibration inputs for the natural-state model

Written by `Scripts/natural_state_analysis.py` from the tables in `Output/`. For the reservoir engineer or agent who will calibrate the natural-state model. Claim strength: **[obs]** = observed in the data; **[inf]** = interpretation; **[ref]** = reference-model property (never evidence); **[assump]** = stated assumption. Wells are named W1 to W10 (WELL-1 to WELL-10 in the data).

## 1. Key conclusions for calibration

1. **[obs] Natural-state surveys (gauge-verified).** The survey pressure at the gauge depth matches the earliest gauge reading for W2 (-0.10 bar, gauge at -750 mRL inside feed zone FZ2) and for W3 (-0.10 bar, gauge at -550 mRL between feed zones). Tolerance ±0.2 bar (two reading precisions). **Calibration:** W2 is the highest-confidence target set (confidence high); W3 medium.
2. **[obs] Survey date conflicts with the gauge series.** Every survey is dated 29/07/2015. The W2 gauge drops from 98.3 bar (01/01/2014) to 89.2 bar in January 2015, so the survey pressures match the pre-drawdown reading, not a 2015 reading. **Data question:** confirm the true survey date. The verdict relies on the gauge test, not on the date label.
3. **[obs] W9 fails the gauge test.** The survey is 2.05 bar above the earliest gauge reading (gauge at -675 mRL, between feed zones). Verdict: possibly disturbed, low confidence. Its feed-zone targets are kept but carry the low factor.
4. **[obs] W10 has no pre-production gauge test.** Its gauge begins 01/01/2020, 31 days after its first flow record (01/12/2019). The survey (29/07/2015) predates that record. Survey-gauge difference is 8.0 bar, which is post-production drawdown, not a test. Verdict: probably natural state, low confidence (feed-zone tag 'Synthetic').
5. **[obs] Wells without a gauge agree with the natural-state trend.** Feed-zone pressures of W1, W4, W6 and W7 lie within 3 bar of the mean of W2 and W3 at the same elevations; W1 is -2.9 bar (low confidence). W5 and W8 (injectors) agree with the trend too, with low confidence because their formations are hotter than the injectate.
6. **[obs] Hot upflow cluster.** Tmax is 268.7 °C in WELL-10 at -250 mRL. The hot wells W2, W3, W6, W9 and W10 all peak between -250 and -750 mRL.
7. **[obs] Strong inversions below the peak.** Inversions of 30.9 °C (W2), 86.9 °C (W3), 98.1 °C (W6) and 57.7 °C (W9) lie below the peak. **[inf]** Outflow of hot fluid at the peak level and/or cold downflow below it; the data alone do not separate these (report §5).
8. **[obs] Near-boiling at W10.** Minimum boiling margin is 2.6 °C at -150 mRL; its feed zone FZ1 (mid -262 mRL) has a margin of 4.3 °C. No point reaches zero margin, so two-phase is not demonstrated. **[inf]** W10 sits on a boiling-point-for-depth profile, which indicates upflow.
9. **[obs] Shallow gradients are not conductive.** Shallow gradients (top 300 m of liquid survey) are 490 to 860 °C/km in W10, W2, W3, W6, against 13.9 to 73.2 °C/km in W4 and W1. **[inf]** Wellbore effects (flow in the open hole or the shallow boiling column) are a more likely cause than a conductive reservoir gradient; treat these gradients as diagnostics, not targets.
10. **[obs] Injectors.** Injectate enthalpy is 335.0 kJ/kg, i.e. 79.8 °C (IF97). W4 is colder than the injectate at every depth with no isothermal section near it, which is evidence against an injection-affected survey. W5 and W8 are hotter than the injectate at depth, so injection cooling could be hidden; their confidence is lowered to low.
11. **[obs] Hydrostatic consistency is mixed.** Deep open windows (top below -900 mRL) give ratios of 0.998 to 1.089 (measured ΔP / IF97 ΔP); W1, W5, W6, W7 lie within 1% of 1, and W2, W8 lie above. Cased windows of W1, W3, W5, W6, W7, W9 reach ratios of up to 1.55, implying liquid densities of up to 1341 kg/m³ (see §4). **[inf]** Not explained by datum, gravity or unit convention (these act at every depth), by two-phase flow (which lowers density) or by rounding (±0.1 bar). Treat the upper pressures as unreliable until checked.

## 2. The data

**Inventory.** `Output/input_inventory.csv` assigns every file under `Baseline/` to a role: well_data_dat 38, reference_other 12, well_data_json 10, out_of_scope 7, reference_results 2, geometry_dat 1, mapping 1, reference_bc 1, reference_rock 1, well_blocks_json 1, well_groups 1.

**Source checks.** `Output/dat_json_crosscheck.csv`: 38 of 38 plain-column series match their JSON counterparts. Production-flow `dat` files store flow with the opposite sign (outflow negative), so magnitudes are compared. Transient and flow time axes are years since 01/01/2014 and match the JSON dates. Irregular names: `Wel10_*` for well 10 (not `Well10_*`).

**Provenance labels** (`Output/provenance_labels.csv`). Survey Source and Comment are blank; StableFlag is true for every survey. The well design Source reads 'Synthetic data' for all ten wells. Feed-zone quality tags are HIGH, MEDIUM, LOW and Synthetic (W10 only). Tag mapping for weight_P: HIGH (highest in the dataset) 1.0; MEDIUM (intermediate) 0.8; LOW (lowest grade) 0.5; Synthetic (not a grade) 0.5. Feed-zone top and bottom are labelled 'm' but are mRL.

**Natural-state timing per well** (`Output/pressure_consistency_check.csv`):

| well | gauge mRL | first gauge | first flow record | gauge minus first flow (days) | verdict | confidence |
|---|---|---|---|---|---|---|
| W10 | -250.0 | 01/01/2020 | 01/12/2019 | 31.0 | probably natural state | low |
| W1 | – | – | none | – | probably natural state | low |
| W2 | -750.0 | 01/01/2014 | 01/01/2014 | 0.0 | natural state | high |
| W3 | -550.0 | 01/01/2014 | 01/01/2014 | 0.0 | natural state | medium |
| W4 | – | – | 01/01/2014 | – | probably natural state | medium |
| W5 | – | – | 01/01/2014 | – | probably natural state | low |
| W6 | – | – | none | – | probably natural state | medium |
| W7 | – | – | none | – | probably natural state | medium |
| W8 | – | – | 01/01/2014 | – | probably natural state | low |
| W9 | -675.0 | 01/01/2014 | 01/01/2014 | 0.0 | possibly disturbed | low |

All surveys are dated 29/07/2015, after the start of the flow records (01/01/2014) for producers and injectors. Production before 2014 is not in the dataset, so the pre-production state is available only through the gauge comparison. W1, W6 and W7 have no flow records; W4, W5 and W8 are injectors; W2, W3 and W9 are producers; W10 started producing on 01/12/2019.

## 3. Temperature structure

| well | Tmax_degC | Tmax_elev_mRL | T_at_TD_degC | TD_mRL | inversion_degC | inversion_top_mRL | inversion_base_mRL | shallow_gradient_degC_per_km | deep_gradient_degC_per_km | z_100C_mRL | z_200C_mRL | min_boiling_margin_degC | regime |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| W10 | 268.7 | -250.0 | 268.5 | -350.0 | – | – | – | 760.4 | – | 164.7 | 45.6 | 2.6 | upflow |
| W1 | 118.6 | -1672.0 | 118.6 | -1672.0 | – | – | – | 73.2 | – | -945.9 | – | 150.4 | conductive |
| W2 | 266.9 | -550.0 | 253.7 | -2348.0 | 30.9 | -550.0 | -1500.0 | 664.4 | 21.4 | 128.2 | -81.0 | 23.8 | outflow |
| W3 | 261.4 | -450.0 | 188.2 | -2887.0 | 86.9 | -450.0 | -1700.0 | 860.4 | 14.3 | 104.5 | -13.7 | 12.2 | outflow |
| W4 | 29.0 | -750.0 | 29.0 | -750.0 | – | – | – | 13.9 | – | – | – | 113.8 | recharge/downflow |
| W5 | 158.0 | -650.0 | 154.7 | -2595.0 | 13.3 | -650.0 | -1700.0 | 295.2 | 13.5 | -159.9 | – | 123.6 | outflow |
| W6 | 220.5 | -550.0 | 126.1 | -2359.0 | 98.1 | -550.0 | -1900.0 | 489.8 | 20.2 | 13.0 | -328.5 | 70.0 | outflow |
| W7 | 176.7 | -650.0 | 144.7 | -2762.0 | 46.0 | -650.0 | -1700.0 | 280.2 | 16.4 | -113.7 | – | 80.5 | outflow |
| W8 | 93.1 | -2716.1 | 93.1 | -2716.1 | 4.5 | -850.0 | -1300.0 | 84.6 | 23.9 | – | – | 104.1 | recharge/downflow |
| W9 | 217.7 | -750.0 | 161.6 | -2427.1 | 57.7 | -750.0 | -2250.0 | 318.7 | – | -84.1 | -496.1 | 87.8 | outflow |

Units: T in °C, elevations in mRL, gradients in °C/km (shallow: top 300 m of liquid survey; deep: below both the inversion base and the deepest feed zone, blank if under 300 m). Inversion: largest drop below a running maximum, if at least 3 °C.

**[obs] Isothermal open-hole sections** (≥300 m, T range ≤2 K; possible internal flow, weight 0.5): W2 -1700 to -1100 mRL; W3 -1900 to -1100 mRL; W5 -1900 to -1500 mRL; W6 -2250 to -1700 mRL; W7 -1900 to -1300 mRL; W8 -1500 to -1100 mRL.

**Regime rules** (heuristic, from survey data only): two-phase if any margin ≤0; upflow if minimum margin <5 °C (near boiling-point-for-depth); recharge/downflow if Tmax <100 °C; outflow if an inversion of ≥3 °C exists; upflow if deep gradient ≤10 °C/km with Tmax ≥200 °C; otherwise conductive.

**Isotherm position:** the 100 °C isotherm lies at W10 165, W1 -946, W2 128, W3 104, W5 -160, W6 13, W7 -114, W9 -84 mRL.

**[obs] Boiling margin below 5 °C (every survey point):** W10 at -75.0 mRL, margin 3.60 °C; W10 at -150.0 mRL, margin 2.65 °C; W10 at -250.0 mRL, margin 3.22 °C (feed zone FZ1); W10 at -262.5 mRL, margin 4.35 °C (feed zone FZ1). No point is at or below zero margin.

## 4. Pressure structure

| well | P at 200 mRL (bar) | P at 0 mRL (bar) | P at -500 mRL (bar) | P at -1000 mRL (bar) | P at -1500 mRL (bar) | P at -2000 mRL (bar) |
|---|---|---|---|---|---|---|
| W1 | – | 12.3 | 63.5 | 114.4 | 162.0 | – |
| W10 | 1.4 | 26.3 | – | – | – | – |
| W2 | 2.1 | 24.1 | 77.5 | 120.0 | 164.5 | 208.6 |
| W3 | – | 20.2 | 76.6 | 119.1 | 163.9 | 209.4 |
| W4 | 4.4 | 24.5 | 76.3 | – | – | – |
| W5 | – | 15.4 | 71.7 | 118.6 | 164.0 | 209.7 |
| W6 | – | 18.6 | 75.7 | 119.1 | 165.1 | 211.4 |
| W7 | – | 18.9 | 73.3 | 118.9 | 164.8 | 211.1 |
| W8 | – | 19.8 | 71.3 | 119.8 | 169.2 | 218.4 |
| W9 | – | 18.3 | 75.6 | 119.3 | 163.9 | 209.4 |

Survey pressures at common elevations, interpolated within each survey only (– = outside the survey range; `Output/pressure_at_common_elevations.csv`).

**Feed-zone pressure trend** (natural-state reference: W2 and W3). Feed-zone pressure minus the reference mean at the feed-zone mid-depth (`feedzones.csv`, column `trend_dev_P_bar`):

| well | feedzone | mid_mRL | quality | P_bar | trend_P_bar | trend_dev_P_bar |
|---|---|---|---|---|---|---|
| W10 | FZ1 | -262.5 | Synthetic | 57.7 | 56.14 | 1.56 |
| W1 | FZ1 | -1337.5 | MEDIUM | 146.7 | 149.62 | -2.92 |
| W1 | FZ2 | -1465.0 | MEDIUM | 158.7 | 161.06 | -2.36 |
| W2 | FZ1 | -260.0 | LOW | 55.9 | 55.91 | -0.01 |
| W2 | FZ2 | -745.0 | HIGH | 97.8 | 97.39 | 0.41 |
| W2 | FZ3 | -1910.0 | MEDIUM | 200.7 | 201.0 | -0.3 |
| W3 | FZ1 | -335.0 | HIGH | 62.8 | 62.93 | -0.13 |
| W3 | FZ2 | -845.0 | MEDIUM | 105.6 | 105.98 | -0.39 |
| W4 | FZ1 | -655.0 | MEDIUM | 91.7 | 89.77 | 1.93 |
| W5 | FZ1 | -1325.0 | MEDIUM | 148.1 | 148.49 | -0.39 |
| W5 | FZ2 | -1525.0 | MEDIUM | 166.3 | 166.44 | -0.14 |
| W6 | FZ1 | -475.0 | HIGH | 73.5 | 74.99 | -1.49 |
| W6 | FZ2 | -1125.0 | HIGH | 130.5 | 130.68 | -0.18 |
| W7 | FZ1 | -630.0 | HIGH | 86.0 | 87.7 | -1.7 |
| W7 | FZ2 | -1325.0 | HIGH | 148.7 | 148.49 | 0.21 |
| W8 | FZ1 | -1125.0 | HIGH | 132.0 | 130.68 | 1.32 |
| W9 | FZ1 | -535.0 | LOW | 78.9 | 79.88 | -0.98 |
| W9 | FZ2 | -837.5 | MEDIUM | 105.2 | 105.34 | -0.14 |

**Hydrostatic checks** (`Output/hydrostatic_check.csv`, per-well medians in `Output/hydrostatic_summary.csv`). Windows are ≥300 m, non-overlapping, split at feed zones and the production casing shoe, and single-phase (P ≥1.5 bar and boiling margin ≥10 °C; no boiling test above the critical pressure). Ratio = measured ΔP / ΔP integrated through IF97 Region 1 density over the measured T profile, g = 9.80665 m/s².

**Field-wide component.** Cold wells (Tmax <100 °C: W4, W8) give median ratios of 1.0402 (cased) and 1.0070 (open). A ratio of 1.04 in a cold column implies density near 1040 kg/m³ (water at these temperatures is about 1000 kg/m³). This is treated as a possible datum or convention issue, not overpressure, and is not confirmed. Rounding uncertainty of a window ratio is at most about 0.003 (`rounding_unc_ratio`).

| well | band | n_windows | median_ratio | cold_well_component_ratio | excess_above_component_pct |
|---|---|---|---|---|---|
| W1 | cased | 1 | 1.0633 | 1.0402 | 2.23 |
| W1 | open | 1 | 1.0077 | 1.007 | 0.07 |
| W2 | open | 3 | 1.0916 | 1.007 | 8.4 |
| W3 | cased | 1 | 1.5507 | 1.0402 | 49.08 |
| W3 | open | 2 | 1.0547 | 1.007 | 4.73 |
| W4 | cased | 1 | 1.032 | 1.0402 | -0.78 |
| W5 | cased | 1 | 1.1315 | 1.0402 | 8.78 |
| W5 | open | 2 | 1.0003 | 1.007 | -0.67 |
| W6 | cased | 1 | 1.2962 | 1.0402 | 24.62 |
| W6 | open | 2 | 1.0038 | 1.007 | -0.32 |
| W7 | cased | 1 | 1.1745 | 1.0402 | 12.92 |
| W7 | open | 2 | 1.0139 | 1.007 | 0.69 |
| W8 | cased | 1 | 1.0483 | 1.0402 | 0.78 |
| W8 | open | 2 | 1.007 | 1.007 | 0.0 |
| W9 | cased | 1 | 1.2064 | 1.0402 | 15.98 |
| W9 | open | 2 | 1.0394 | 1.007 | 3.22 |

**Windows (all, with implied densities):**

| well | band | z_top_mRL | z_bottom_mRL | length_m | dP_measured_bar | dP_IF97_bar | ratio_measured_to_IF97 | rho_implied_kgm3 |
|---|---|---|---|---|---|---|---|---|
| W1 | cased | 15.0 | -1000.0 | 1015.0 | 103.6 | 97.434 | 1.0633 | 1040.8 |
| W1 | open | -1100.0 | -1672.0 | 572.0 | 54.0 | 53.589 | 1.0077 | 962.7 |
| W2 | open | -150.0 | -720.0 | 570.0 | 51.768 | 44.474 | 1.164 | 926.1 |
| W2 | open | -720.0 | -1900.0 | 1180.0 | 104.132 | 95.396 | 1.0916 | 899.9 |
| W2 | open | -1900.0 | -2348.0 | 448.0 | 39.3 | 36.099 | 1.0887 | 894.5 |
| W3 | cased | 150.0 | -250.0 | 400.0 | 52.6 | 33.92 | 1.5507 | 1340.9 |
| W3 | open | -250.0 | -820.0 | 570.0 | 48.442 | 44.472 | 1.0893 | 866.6 |
| W3 | open | -820.0 | -2887.0 | 2067.0 | 185.358 | 181.715 | 1.02 | 914.4 |
| W4 | cased | 230.0 | -500.0 | 730.0 | 73.9 | 71.611 | 1.032 | 1032.3 |
| W5 | cased | 93.0 | -1000.0 | 1093.0 | 113.667 | 100.454 | 1.1315 | 1060.5 |
| W5 | open | -1100.0 | -1500.0 | 400.0 | 36.3 | 36.331 | 0.9991 | 925.4 |
| W5 | open | -1500.0 | -2595.0 | 1095.0 | 100.1 | 99.953 | 1.0015 | 932.2 |
| W6 | cased | 147.0 | -440.0 | 587.0 | 68.31 | 52.7 | 1.2962 | 1186.7 |
| W6 | open | -450.0 | -1100.0 | 650.0 | 56.9 | 56.335 | 1.01 | 892.6 |
| W6 | open | -1100.0 | -2359.0 | 1259.0 | 116.6 | 116.893 | 0.9975 | 944.4 |
| W7 | cased | 150.0 | -440.0 | 590.0 | 64.95 | 55.3 | 1.1745 | 1122.6 |
| W7 | open | -450.0 | -1300.0 | 850.0 | 78.3 | 76.204 | 1.0275 | 939.3 |
| W7 | open | -1300.0 | -2762.0 | 1462.0 | 135.1 | 135.043 | 1.0004 | 942.3 |
| W8 | cased | 176.0 | -700.0 | 876.0 | 89.25 | 85.134 | 1.0483 | 1038.9 |
| W8 | open | -750.0 | -1100.0 | 350.0 | 33.6 | 33.789 | 0.9944 | 978.9 |
| W8 | open | -1100.0 | -2716.1 | 1616.1 | 159.1 | 156.037 | 1.0196 | 1003.9 |
| W9 | cased | 170.0 | -400.0 | 570.0 | 64.05 | 53.092 | 1.2064 | 1145.8 |
| W9 | open | -450.0 | -830.0 | 380.0 | 33.757 | 31.967 | 1.056 | 905.9 |
| W9 | open | -830.0 | -2427.1 | 1597.1 | 143.943 | 140.722 | 1.0229 | 919.0 |

Excluded windows (`hydrostatic_check.csv` covers the used windows): WELL-10 cased: single-phase part 0 m, no window >= 300 m; WELL-10 open: single-phase part 0 m, no window >= 300 m; WELL-1 open: no survey points between -1100 and -1000 mRL; window top set to -1100 mRL; WELL-2 cased: single-phase part 290 m, no window >= 300 m; WELL-2 open: no survey points between -150 and -100 mRL; window top set to -150 mRL; WELL-3 cased: survey points above 164 mRL not single-phase; window top set to 150 mRL; WELL-4 open: single-phase part 200 m, no window >= 300 m; WELL-5 open: no survey points between -1100 and -1000 mRL; window top set to -1100 mRL; WELL-6 open: no survey points between -450 and -440 mRL; window top set to -450 mRL; WELL-7 cased: survey points above 170 mRL not single-phase; window top set to 150 mRL; WELL-7 open: no survey points between -450 and -440 mRL; window top set to -450 mRL; WELL-8 open: no survey points between -750 and -700 mRL; window top set to -750 mRL; WELL-9 cased: survey points above 174 mRL not single-phase; window top set to 170 mRL; WELL-9 open: no survey points between -450 and -400 mRL; window top set to -450 mRL.

**[obs]** Deep open windows (top below -900 mRL) give ratios of 0.998 to 1.089; W1, W5, W6, W7 agree with IF97 to within 1%, and W2, W8 are above it. Cased windows of W1, W3, W5, W6, W7, W9 reach ratios of up to 1.55. **[inf]** The upper survey pressures are steeper than any water column with the measured temperatures, so they are inconsistent with the temperature profile. Checks: a datum offset changes pressure, not gradient; a gravity or unit convention would act at every depth (the deep windows show it does not); two-phase would reduce density (the upper windows are single-phase); rounding is ±0.1 bar. The likely causes are a wellbore flow effect on the shallow pressure record or a measurement issue. This needs a check of the raw gauge data (report §7).

## 5. Conceptual model

Sketch: `Output/figures/conceptual_model.png` (plan view, regime per well). Fence section through the four hottest wells: `Output/figures/section_hottest_wells.png` (survey points projected onto one line; not an interpolated section).

- **[obs] Hot cluster.** W2, W3, W6, W9 and W10 lie within a few kilometres of each other, with peak temperatures of 218 to 269 °C. Each lies 113 to 378 m from a mapped fault (`well_geometry_check.csv`). **[inf]** The hot cluster lies on or near the fault system, which is a candidate upflow pathway. The data do not prove this.
- **[inf] Upflow:** the near-boiling profile at W10 (margin down to 2.6 °C) is the strongest upflow indicator. Its shallow feed zone sits at -262 mRL.
- **[inf] Outflow:** the hot wells peak at -450 to -750 mRL, then cool by 31 to 98 °C below that level (W2, W3, W6, W7, W9). This is consistent with lateral outflow at the hot level. Cold downflow would give the same profile; a temperature transect across wells at that depth would separate them.
- **[obs] Recharge:** W4 (29 °C at maximum) and W8 (93 °C) are cold. W8 lies 235 m from the nearest water body. **[inf]** Recharge or downflow is likely at these locations. Their pressures agree with the trend, so the cold does not require a pressure anomaly.
- **[obs] Cap.** The 100 °C isotherm lies at about -84 to 165 mRL in the hot wells and at -946 mRL in the conductive well W1. **[inf]** The hot wells have a 100 °C isotherm close to surface, so they have little cool cap; the conductive well W1 stays below 100 °C down to -946 mRL. Cap base cannot be resolved from the data.
- **[obs] Pressure connection.** Feed-zone pressures agree with the natural-state trend to within 2.9 bar (largest deviation at W1; table in §4). **[inf]** Deep pressures are hydraulically connected across the mapped area, so a single natural-state pressure field is a reasonable starting point.

## 6. Recommended calibration set-up

- **Targets.** Use `NS_calibration_targets.csv`. Temperature weight = location weight × confidence factor (high 1.0, medium 0.75, low 0.5). Pressure targets only at feed zones (`weight_P`), with quality factor HIGH 1.0, MEDIUM 0.8, LOW and Synthetic 0.5, times confidence.
- **Confidence.** Weight W2 targets fully. W3 and the trend-check wells at medium. W5, W8, W9, W10 and W1 at low. Confidence factors vary across wells, so they change relative weights.
- **[inf] Heat in.** Anchor the upflow at the hot cluster. W2 FZ2 (mid -745 mRL, 97.8 bar, quality HIGH, confidence high) is the best-constrained hot feed zone.
- **Heat and mass out.** A lateral outflow at the hot level (about -450 to -750 mRL) and cooling below it; a cold recharge or downflow sink near W4 and W8.
- **[inf] Permeability.** The data do not constrain permeability directly. The hot feed zones sit in the slotted liner and their pressures follow the natural trend, so they are hydraulically connected to the field. Calibrate against the feed-zone pressures in `feedzones.csv`.
- **Boundaries.** Pressure set from the trend-check reference (W2 and W3); reference-model boundary conditions are not evidence and were not used.

## 7. Data quality issues and open questions (ranked by effect on calibration)

0. **All ten well designs are labelled 'Synthetic data'** (`provenance_labels.csv`), and W10's feed zone is tagged 'Synthetic'. The dataset may be a synthetic teaching model rather than field data. The results below are conclusions about the dataset as supplied; confirm the data provenance before using them for a field model.
1. **Upper pressures inconsistent with temperatures** (§4): cased windows of W1, W3, W5, W6, W7, W9 imply densities up to 1341 kg/m³. Ask for the raw gauge and survey records for the top of these wells and whether they were flowing at survey time.
2. **Survey date label** (29/07/2015) conflicts with the gauge drawdown in January 2015. Confirm the date before using it to set a time for the natural-state state.
3. **W9 gauge disagreement** (2.05 bar), with the gauge between feed zones. Confirm the gauge depth and whether the gauge reading is a feed-zone pressure.
4. **W10 'Synthetic' feed-zone tag**, no pre-production gauge and near-boiling margin. Confidence low; feed-zone weights reduced.
5. **Injectors W5 and W8** are hotter than the injectate at depth: injection cooling could be hidden. Their confidence is low.
6. **Shallow gradients** of 490 to 860 °C/km in the hot wells: not conductive. Check for wellbore flow or a wrong shallow depth reference.
7. **Wellhead position.** The wellhead lies below the top of its shallowest model block by up to 16 m, in 8 of 10 wells (`well_geometry_check.csv`, `wellhead_minus_surface_m`). Confirm the model top or the well head elevation.
8. **Flow records start at 01/01/2014** for all wells, so the true first production date is unknown.
9. **Reference files** contain hard-coded paths from another machine (`gTE14515.auxjson`, `currentGeomSetup.json`); they were not used as evidence.
10. **Naming and sign conventions:** `Wel10_*` file names; production-flow `dat` files are sign-reversed relative to the JSON.

## 8. Output file index

- `Output/input_inventory.csv`: role and path of every input file
- `Output/NS_calibration_targets.csv`: one row per survey point with weights
- `Output/feedzones.csv`: feed-zone values, quality factors, trend check
- `Output/T_P_at_layer_centres.csv`: T and P at model layer centres (no extrapolation)
- `Output/well_summary.csv`: per-well characteristics, regime, verdict
- `Output/pressure_consistency_check.csv`: gauge test, verdict and confidence per well
- `Output/hydrostatic_check.csv, hydrostatic_summary.csv`: hydrostatic windows and per-well medians
- `Output/well_geometry_check.csv`: wellhead, model top, faults and water bodies
- `Output/pressure_at_common_elevations.csv`: T and P at common elevations
- `Output/dat_json_crosscheck.csv`: JSON versus plain-column data
- `Output/provenance_labels.csv`: source, comment and quality labels
- `Output/analysis_parameters.csv`: all thresholds used
- `Output/figures/`: well map, profiles, maps, fence section, P–T plot, conceptual sketch

