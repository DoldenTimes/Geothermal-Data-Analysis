# Natural-state interpretation of downhole P/T and well geometry

Dataset: `Baseline/Teaching_Model_3D/Data` (Data_json and Data_dat; the two were cross-checked and are identical).
Script: `Scripts/natural_state_analysis.py`. All numbers below come from the CSV/JSON files in this folder.
Model files under `Baseline/.../Models` were not used as data (per CLAUDE.md).

## 1. What the data are

- 10 wells (WELL-1 to WELL-10). Each has **one** downhole temperature and one downhole pressure survey, dated **29/07/2015**, both flagged stable. These are the only P/T observations in the dataset, so they are the natural-state reference used here.
- Geometry: wellhead elevations 15–340 masl; total depths −350 to −2887 masl. WELL-8 and WELL-9 are deviated (20 and 25 track points); the others are vertical. All measurements and all feed zones lie inside the well tracks, and all feed zones lie inside the slotted liner intervals.
- Feed zones: 18 in total across the wells (LOW / MEDIUM / HIGH quality tags, plus one "Synthetic" in WELL-10). See `feedzones_summary.csv`.

**Timeline caveat (important).** Production or injection records exist before the 2015 survey for WELL-2, 3, 9 (production from 2014), and WELL-4, 5, 8 (injection from 2014). WELL-10 production starts 2019, so its 2015 survey predates production. WELL-1, 6, 7 have no flow records. The 2015 surveys are therefore natural-state only in the strict sense for WELL-1, 6, 7 and (probably) WELL-10. For the other wells they may already reflect early production or injection. The cold profiles in the injection wells (WELL-4, WELL-8) are the clearest sign of this. The engineer should confirm the survey timing against the field history before using these as hard natural-state targets.

## 2. Data checks

| Check | Result | File |
|---|---|---|
| Data_dat vs Data_json (T and P, all wells) | 100% match, max difference 0 | `dat_vs_json_check.csv` |
| Measured points inside well track | All yes | `well_geometry_check.csv` |
| Feed zones inside track and liner | All yes | `well_geometry_check.csv` |
| Saturation-temperature function | Matches IAPWS-IF97 reference values to 0.01 °C (1, 10, 100 bar; 1.01325 bar → 99.97 °C) | run log |
| Pressure gradient vs hydrostatic column (Kell density) | Measured columns steeper than hot-water columns in every well (see §4) | `hydrostatic_check.csv` |
| Pressure resolution | 0.1 bar. Short-segment densities (10–20 m spacing) are noisy (±50 kg/m³); use column-scale numbers | — |

## 3. Thermal regimes

Temperature (°C) at selected elevations (from `T_P_at_common_elevations.csv`):

| Well | 0 masl | −500 | −750 | −1000 | −1500 | Min. Tsat−T (°C) | Regime |
|---|---|---|---|---|---|---|---|
| WELL-2 | 172 | 266 | 258 | 241 | 236 | 24 | Hot upflow |
| WELL-3 | 193 | 261 | 244 | 198 | 176 | 12 | Hot upflow |
| WELL-9 | 64 | 201 | 218 | 195 | 182 | 88 | Hot upflow |
| WELL-10 | 219 | — | — | — | — | **2.6** | Hot, near boiling |
| WELL-6 | 106 | 219 | 194 | 154 | 127 | 70 | Intermediate |
| WELL-5 | 58 | 153 | 157 | 155 | 145 | 124 | Intermediate (injector) |
| WELL-7 | 68 | 164 | 174 | 154 | 131 | 81 | Intermediate |
| WELL-1 | 35 | 73 | 90 | 102 | 116 | 150 | Cold |
| WELL-4 | 15 | 24 | 29 | — | — | 114 | Cold (injector) |
| WELL-8 | 27 | 62 | 68 | 66 | 64 | 104 | Cold (injector) |

Key points for calibration:

1. **Lateral temperature contrast is large.** At −1000 masl temperatures range from 66 °C (WELL-8) to 241 °C (WELL-2). At −1500 masl the range is 64 to 236 °C. The hot wells (WELL-2, 3, 9, 10) sit together on the western side of the field, while the cold wells (WELL-1 to the east, WELL-4 to the north-east, WELL-8 to the south-west) sit on its edges; see `05_T_P_maps_by_elevation.png`.
2. **WELL-10 is within about 3–6 °C of boiling** from −25 to −262 masl (Tsat−T = 2.6–5.8 °C; the minimum of 2.6 °C is at −150 masl). Its top-third gradient is about 132 °C/100 m. Any natural-state model should either hold temperature at or below Tsat(P) here, or allow a boiling/two-phase zone. This is the most constraining single observation.
3. **WELL-3 is also close to boiling** (minimum Tsat−T of 12 °C in the well; the location is in `downhole_profiles_long.csv`).
4. **Temperature inversions** (temperature decreasing with depth) occur in several wells:
   - WELL-3, MEDIUM feed zone −820 to −870 m: −16.8 °C/100 m
   - WELL-2, HIGH feed zone −720 to −770 m: −7.1 °C/100 m
   - WELL-9, −830 to −845 m: −7.1 °C/100 m
   - WELL-6, −1100 to −1150 m: −4.6 °C/100 m
   These suggest cooler inflow or downflow below the feed zones, and they are a useful constraint on the feed-zone locations and the vertical flow direction.
5. **Isothermal sections** (convective or well-mixed) occur at WELL-7 −1300 to −1350 m (131.3 → 131.2 °C) and WELL-5 −1300 to −1550 m (≈147 → 145 °C).

## 4. Pressure: consistency with hydrostatic and across wells

**Hydrostatic test (per well).** Starting from the shallowest measurement, pressure was integrated down the well using the measured temperature profile and the Kell water density. The residual (measured minus predicted) grows with depth in every well:

| Well | Deepest residual (bar) | Column-mean density implied (kg/m³) | Kell density (kg/m³) | Excess |
|---|---|---|---|---|
| WELL-1 | +7.7 | 1012 | 965 | +4.8% |
| WELL-2 | +31.4 | 944 | 819 | **+15.2%** |
| WELL-3 | +29.0 | 961 | 864 | +11.2% |
| WELL-4 | +2.2 | 1021 | 998 | +2.3% |
| WELL-5 | +15.4 | 983 | 925 | +6.3% |
| WELL-6 | +17.7 | 988 | 916 | +7.9% |
| WELL-7 | +11.3 | 969 | 930 | +4.2% |
| WELL-8 | +8.8 | 1011 | 980 | +3.2% |
| WELL-9 | +18.7 | 969 | 895 | +8.2% |
| WELL-10 | +17.6 | 1154 | 834 | **+38.4%** |

Reading this:
- The measured pressure columns are steeper than a column of hot water at the measured temperatures. A 15% excess would need very saline brine (order 20 wt% NaCl). Typical geothermal salinities change density by only a few percent, so salinity alone is unlikely to explain the difference.
- Likely causes, to check in the source data: (a) the temperatures are too high (a cooler column would be denser); (b) the pressure readings have a depth or gauge-reference bias; (c) the datum or depth used in the pressure conversion is inconsistent with the temperature survey.
- WELL-10 (+38%) is the strongest outlier. Its pressure gradient of about 0.12 bar/m is not physically possible for water at its temperatures. Check its pressure record first.
- Kell's polynomial is fitted to 0–150 °C; the wells here reach 270 °C. The errors above that are probably a few percent, smaller than the 10–15% discrepancies in WELL-2 and WELL-3, so the conclusion holds. An IAPWS-IF97 density would confirm this, but the library was not available in the sandbox.

**Cross-well comparison at common elevations.** Comparing pressure at the same elevation avoids the column reduction and is more reliable than the datum-reduced values (`pressure_datum_reduction.csv`; the datum-reduced spread of 17 bar mostly reflects the column mismatch above, not real lateral differences):

| Elevation (masl) | Field range (bar) | Notes |
|---|---|---|
| −250 | 38.6 (WELL-1) to 56.7 (WELL-10) | 18 bar spread. WELL-1 is ~10 bar low; WELL-10 high |
| −500 | 63.5 (WELL-1) to 77.5 (WELL-2) | WELL-1 low by ~8–14 bar |
| −750 | 89.0 (WELL-1) to 100.5 (WELL-4) | Excluding WELL-1, 95–100 bar (±3 bar) |
| −1000 | 114.4 (WELL-1) to 120.0 (WELL-2) | Excluding WELL-1, 118.6–120.0 bar |
| −1500 | 162.0 (WELL-1) to 169.2 (WELL-8) | Field ≈ 164 bar; WELL-8 +5 bar, WELL-1 −2 bar |

So below about −750 masl the pressure field is close to uniform (within about ±3 bar, with WELL-1 and WELL-8 slightly off). Shallower than −500 masl, WELL-1 is consistently low and needs a separate explanation (cold, low-pressure, possibly a different shallow regime or a gauge issue). The field is therefore one connected pressure system at depth, and the shallow contrast is the main pressure-side calibration tension.

**Shallow pressure and free surface.** The shallowest pressures are close to atmospheric (1.0 bar at 210–225 masl in WELL-2 and WELL-10). This places the free surface at about the wellhead elevation there. WELL-1 is different: 10.8 bar at 15 masl implies a free surface near 118 masl, roughly 100 m above its wellhead. This is a useful boundary constraint for the shallow model.

## 5. Feed zones

See `feedzones_summary.csv` for values at each zone top and bottom.

- Feed zones that show a warm-up with depth (T increasing into the zone) are the most likely inflow points. WELL-2 LOW (−250 to −270 m; 238 → 241 °C, 14.9 °C/100 m) and WELL-3 HIGH (−310 to −360 m; 258 → 260 °C) are examples.
- Zones with near-zero gradient and similar temperature at top and bottom (WELL-7 −1300 to −1350 m; WELL-5 −1500 to −1550 m) are weak evidence for flow across the zone.
- The "Synthetic" quality tag on WELL-10 indicates the zone is a placeholder. Its temperature, pressure and boiling margin should be treated accordingly.

## 6. Summary of calibration targets (natural state)

1. Temperatures at the common elevations in `T_P_at_common_elevations.csv` (use the no-extrapolation values only).
2. Temperature and pressure at each feed-zone top and bottom in `feedzones_summary.csv`.
3. Tsat margin constraint: T ≤ Tsat(P) everywhere, with WELL-10 at 3 °C margin (possible two-phase zone) and WELL-3 at 12 °C.
4. Pressure field: uniform to ±3 bar below −750 masl; shallow contrast and WELL-1 offset to be explained.
5. Temperature inversions (WELL-2, 3, 6, 9) as constraints on downflow below feed zones.
6. Free-surface level of about 210–225 masl at WELL-2 and WELL-10; about 118 masl at WELL-1.

Targets to exclude or treat with care: WELL-4, 5, 8 (injection since 2014; cold and possibly injection-affected); WELL-2, 3, 9 (production since 2014); WELL-10 pressures (inconsistent with hydrostatic, see §4).

## 7. Data issues to resolve before calibrating

1. Confirm the survey date against production/injection history for all wells (§1).
2. Confirm the pressure reference (absolute or gauge, and depth reference) and check WELL-10's pressure record.
3. Check whether the temperature profiles in WELL-2 and WELL-3 are too hot for the measured pressures (§4).
4. Replace the "Synthetic" WELL-10 feed zone with a measured one if available.

Not used here: the transient pressure, production and injection flow and enthalpy records. They describe production/injection periods and are outside the natural-state scope (CLAUDE.md Timeline). They will be needed for the production calibration.
