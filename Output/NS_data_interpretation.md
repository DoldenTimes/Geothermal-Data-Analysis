# Natural-state data interpretation — TE14515 teaching model (10 synthetic wells)

Purpose: give an engineer / calibration agent a compact, data-only description of the natural state, with explicit calibration targets and flagged problems.
Source: `Baseline/Teaching_Model_3D/Data` (surveys, feed zones, tracks, GIS) and `Models/TE14515/gTE14515*` (geometry, well–block mapping, boundary conditions). The supplied natural-state model (`TE14515_NS_000.*`) was **not** used as data. Production history was ignored except for the *first* (initial-condition) values.
Reproduce: `python Scripts/natural_state_analysis.py` (outputs listed at the end).

## 1. What data exist

| Item | Content |
|---|---|
| Wells | 10 (W1–W10), all one downhole **T and P survey dated 29/07/2015** (same elevations for T and P), 17–33 points per well, 5–150 m spacing. No repeat surveys, no shut-in/heat-up series. |
| Geometry | 1 km-ish columns (393), 70 layers: 20 m thick from 1000 to 0 mRL, then 50–100 m down to −1000, 200 m to −2000, 500 m to −3500. Well→block mapping and feed-zone→block mapping provided (250 well/block intersections). |
| Feed zones | 18 zones, quality tags LOW/MEDIUM/HIGH (W10 "Synthetic"). |
| Roles | Producers W2, W3, W9, W10 (initial pressure, enthalpy and rate available); injectors W4, W5, W8 (cold, 335 kJ/kg ≈ 80 °C); W1, W6, W7 monitor-only. W8 and W9 are deviated. |
| Boundary conditions | Atmosphere 1.013 bar, surface T 11–12.7 °C by column; rain 1000 mm/yr × 0.1 infiltration (uniform); **uniform base heat flux 0.1 W/m²** (no per-column entries); 5 lake columns (12 °C, 5–25 m depth) in the far SW/NW corner. |

Data quality notes: pressures are rounded to 0.1 bar, so 20 m-step local gradients are noisy (±5 %); only multi-hundred-metre gradients should be used. Survey tops are at the collar except W4 (230 mRL vs collar 340) and W10 (210 vs 260).

## 2. Is the 2015 survey a natural-state dataset?

Production starts 01/01/2014, but the 2015 survey pressure agrees with the *first* (2014) transient pressure at the gauge depth: W2 98.2 vs 98.3 bar, W3 80.7 vs 80.8, W9 91.2 vs 89.2 (−2.0). So the surveys behave as **undisturbed (pre-production) conditions** and are a legitimate natural-state target. Exception: W10's transient series starts in 2020 (48.7 bar vs 56.7 bar survey): it is already depleted by ~8 bar, so W10's *survey* is the natural-state value and its transient data are not (`pressure_consistency_check.csv`).

## 3. Temperature structure

Maximum survey temperatures and the shape of each profile (`well_summary.csv`, fig 02, 03, 07):

| Well | Tmax °C @ mRL | T at bottom | Shape |
|---|---|---|---|
| W10 | 268.7 @ −250 | 268.5 @ −350 | Hottest, shallowest; follows boiling curve at the base |
| W2 | 266.9 @ −550 | 253.7 @ −2348 | Hot core, inversion to 236 @ −1500, rising again below |
| W3 | 261.4 @ −450 | 188.2 @ −2887 | Hot core, strong inversion (−87 °C) to ~175 |
| W6 | 220.5 @ −550 | 126.1 | Inversion −98 °C, near-isothermal ~125 below −1900 |
| W9 | 217.7 @ −750 | 161.6 | Inversion −58 °C |
| W7 | 176.7 @ −650 | 144.7 | Inversion −46 °C |
| W5 | 158.0 @ −650 | 154.7 | Inversion small (−13 °C) |
| W1 | 118.6 @ −1672 (bottom) | 118.6 | Near-linear, no inversion |
| W8 | 93.1 @ −2716 (bottom) | 93.1 | Max ~68 °C @ −850, slight reversal then rise; cold system |
| W4 | 29.0 @ −750 (bottom) | 29.0 | Background conductive (~1.7 °C/100 m); no thermal anomaly |

Interpretation:
* **Hot upflow/source zone in the W10–W2–W3 corner** (NW of centre). All three reach 260–269 °C between −250 and −550 mRL. This sits within ~250–380 m of Fault 1 (N–S line at x≈2.8196–2.8198 M) and ~125–230 m of Fault 3 (NE–SW line); the two faults intersect at ≈(2819780 E, 6335480 N), ~250 m from W2 (`well_geometry_check.csv`). A fault-intersection-controlled upflow is a reasonable working hypothesis; W9 (270 m from Fault 2, 684 m from Fault 3) is the next hottest producer.
* **Lateral outflow / cooling away from the source.** Tmax drops W2/W3 (~265) → W6/W9 (~220) → W7 (177) → W5 (158) → W1 (<120) → W8 (<70) → W4 (cold) with distance toward S, SE, E, W and N. Outflow is **shallow and layered**: the Tmax depth is at −450…−750 mRL in every hot well, and every well with a hot core shows a **temperature inversion** below it (to ~125–240 °C by −1500…−1900). The inversion is the signature of lateral outflow over cooler water (or a cold downflow beneath), not a vertical-upflow profile; the model needs permeable layers at roughly −400…−800 mRL carrying hot fluid outward.
* **Clay cap / conductive lid** (above ~−100 to 0 mRL near the source): gradients of 138 °C/100 m (W2 @ ~160 mRL), 148 (W3 @ 120), 241 (W10 @ 180) and 62 (W6), 47 (W9), 36 (W7), 34 (W5). Isotherm elevations (mRL): T=100 °C at W10 +165, W3 +104, W2 +128, W6 +13, W9 −84, W7 −114, W5 −160, W1 −946; T=200 °C at W10 +46, W3 −14, W2 −81, W6 −329, W9 −496; none reach 200 °C in W1/W4/W5/W7/W8. Thickness of the cap decreases toward the source — a near-source resistivity-cap equivalent is needed that is thin (~100–250 m) over W10/W2/W3 and thickens outward. Surface-layer temperatures of 25–35 °C in W10/W3/W1 (vs. ~12 °C ambient in W4) suggest weak shallow thermal manifestation there.
* **Near-boiling only at W10**: the W10 profile reaches 2.6 °C of the boiling-point-for-depth curve at −150 mRL and tracks it to the bottom; W3 min margin 12 °C @ −150, W2 24 °C @ −450. All other wells are 70–150 °C sub-boiling. So two-phase conditions are essentially limited to the shallow-hot zone at W10, and the rest of the model should remain liquid.
* **Deep gradients**: deeper than the inversion, bottom temperatures rise again in all wells: ≈1.7 °C/100 m (W4 over its full 980 m), ≈2–3 °C/100 m (W8 −1500…−2716), ≈1 °C/100 m (W3 below −1700), ≈2 °C/100 m (W2 −1900…−2348). These background gradients (≈1–3 °C/100 m) correspond to a conductive heat flux of about 30–75 mW/m² *if* rock conductivity ≈2.5 W/m/K (assumption; that value is only from the reference model, not data). The provided boundary heat flux of 0.1 W/m² (100 mW/m²) therefore looks high for the cold periphery and — being uniform — cannot create a 265 °C anomaly by itself; the hot source must be represented by a localised deep hot-fluid source (mass + enthalpy) or a hot basal area.

## 4. Pressure structure

* **Water table** (extrapolated to ≈1 bar, `well_summary.csv`): W1 115, W5 133, W6 157, W3 168, W9 178, W8 184, W7 190, W10 210, W2 225, W4 244 mRL. Water tables are essentially at the surface in the low, east/south-east part of the field (W1, W5, W6: within 3–7 m of the model surface elevation) and **deeper below the model surface in the high ground** (W8 −76 m, W4 −96 m, W7/W10 −50 m, W9 −42 m, W2 −35 m). Head decreases from north (W4, 244) to east/south-east (W1, 115): about 130 m over ~3.9 km, consistent with groundwater draining toward the Tarawera River / east. W4 and W10 values are extrapolated (survey starts below the collar).
* **Deep pressure is laterally almost uniform**: at −1500 mRL: W1 162.0, W2 164.5, W3 163.9, W5 164.0, W6 165.1, W7 164.8, W9 163.9 bar, and W8 169.2 bar (cold column, highest). Hot and cold columns differ by only ≤5 bar there, although a 20 °C versus 200–260 °C column over ~1.7 km would differ by ~15–25 bar if each well were in hydrostatic equilibrium with its own column.
* **Consequently the surveys are not self-hydrostatic**: integrating the measured temperature with liquid density (`hydrostatic_check.csv`, columns `P_resid`) leaves P exceeding the local-temperature hydrostatic prediction at the base by +30 bar (W2), +27.5 (W3), +17.8 (W9), +17.4 (W10), +16.9 (W6), +14.6 (W5), +10.3 (W7), +7.9 (W8), +7.3 (W1) and +2.1 (W4). Part of this is a *systematic ≈2–3 % steeper-than-fresh-water gradient* that even the cold wells W4/W8 show (implied density ≈1020 vs ≈1000 kg/m³), so only the *excess over that baseline* is meaningful: it is largest for the hot wells (W2, W3, W10, W9, W6) — hot columns are over-pressured relative to their own buoyancy, i.e. **pressure at depth is imposed by a common deep pressure source/lateral connection, not by each column's weight**. For calibration this means: (a) the model's deep pressure must be matched with a deep, laterally connected, pressure-supported region (open or pressure-boundary base/sides, or high-permeability deep layers); (b) the sealed lateral boundaries + low-permeability base typical for starting models would give hotter-is-lower pressure and will not reproduce the data; (c) the ~2–3 % gradient excess should be checked against the model's density/gravity convention before treating it as physical overpressure.
* **Pressure vs saturation curve**: at depth P is 200–290 bar at −1900…−2900 mRL while T ≤ 255 °C, far from boiling (cf. `07_PT_state.png`). Only W10 is near saturation.

## 5. Feed zones, enthalpy and productivity clues (initial values only)

(`feedzones.csv`, `pressure_consistency_check.csv`)

| Well | First h (kJ/kg) | liquid-T equivalent | Survey T at feed zones | Comment |
|---|---|---|---|---|
| W2 | 1301.6 | 292 °C | 240–258 °C (h_f 1036–1124) | **Excess enthalpy ≈ +180–260 kJ/kg** over the surveyed liquid; implies steam/two-phase contribution or hotter inflow than the survey shows. The model must supply an enthalpy-rich feed at W2's −745 mRL feed (HIGH quality). |
| W3 | 1087.9 | 250 °C | 259 / 230 °C | Consistent with a mix of FZ1 (−335) and FZ2 (−845). |
| W9 | 929.6 | 217 °C | 205 / 212 °C | Consistent (slightly higher). |
| W10 | 2681–2713 | steam-dominated (≈93 % steam at ~270 °C) | 268.7 °C @ 57.7 bar (Tsat 273) | Near-dry steam feed at the hot, shallow, near-boiling zone; start 12/2019 at ~0 kg/s. |

Initial rates (kg/s): W2 85, W3 130, W9 200, W10 ≈0.1. Injection initial rates: W4 110, W5 100, W8 100 kg/s at 335 kJ/kg. Injectors are on the periphery in cold/marginal zones: W4 (T=27 °C at feed), W8 (64 °C at feed), W5 (145 °C at feeds) — wells with feeds at −655, −1125 and −1325/−1525 mRL, respectively.

All surveyed feed-zone states are well sub-boiling except W10 (T within ~4 °C of Tsat). Deep feeds of the monitor/injection wells (W1, W5, W6-FZ2, W7-FZ2, W8) sit at 130–166 bar and 64–147 °C at −1100…−1550 mRL.

## 6. Geometry/consistency issues the modeller should know

* Collar elevations are lower than the model surface in every column except W4 and W10 (differences −13 m W6 to −105 m W1; `well_geometry_check.csv`) because of 1 km column smoothing. W1 (collar 15 mRL, water level +115) is *artesian* by ~100 m relative to its collar — the surface of W1's column (120 mRL) is well above the real ground.
* Feed-zone blocks given in the well file: 1–2 blocks per well; the survey-mean T and P over each block interval (and at the layer centre) are tabulated for 243 of 250 blocks (the remaining 7 are above/below surveyed range) in `NS_calibration_targets.csv`.
* The lake (far SW, ≥3.2 km from any well) and rivers matter only at the top boundary: wells W1, W8, W5 and W6 lie 96–424 m from the Tarawera river branches (surface head control), the other wells are >800 m away.

## 7. Recommended calibration targets and weighting

1. **Temperature** at every surveyed well-block interval (use `T_mean_over_well_interval`): weight highest in the hot corridor (W2, W3, W10, W6, W9) and on the cap-base gradient and Tmax elevation; the *inversion depth/amplitude* (W3 −87, W6 −98, W9 −58, W7 −46, W2 −31) is a strong constraint on outflow geometry.
2. **Pressure** at the same intervals (`P_mean_over_well_interval_bar`) with ±1–2 bar tolerance at depth; match the near-uniform ~164 bar at −1500 mRL and the water-table spread 115–244 mRL.
3. **Enthalpy/state**: W10 near-boiling shallow (T within ~3 °C of Tsat at −150…−350 mRL), W2 needs enthalpy ≳1300 kJ/kg at its main feed (or an explanation); others liquid.
4. **Qualitative**: no boiling elsewhere; upflow near the Fault 1/Fault 3 intersection; cold periphery N (W4) and W (W8).
5. **Treat with caution**: W10 transient pressures (depleted), W9 pressure offset (−2 bar vs first transient), 0.1 bar rounding, 2–3 % gradient excess.

## 8. Output files

* `well_summary.csv` – per-well metrics (Tmax, inversion, gradients, boiling margin, water-level, hydrostatic residual)
* `feedzones.csv` – feed-zone T/P (survey), block names
* `NS_calibration_targets.csv` – 250 well-block rows: layer top/centre/bottom and T,P targets
* `T_P_at_common_elevations.csv` – T,P at 200…−2500 mRL for inter-well comparison
* `hydrostatic_check.csv` – T, P, Tsat(P), boiling margin, T-integrated hydrostatic P and residual at every survey point
* `pressure_consistency_check.csv` – survey vs first transient pressure/enthalpy/rate
* `well_geometry_check.csv` – collar vs model surface, distances to faults/rivers/lake
* `figures/01…07` – map, T/P profiles, per-well panels, T and P maps by elevation, two cross-sections, P–T state vs saturation
* `process_log.md` – prompt and all steps
