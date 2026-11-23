# Natural-State Data Interpretation: Teaching_Model_3D (TE14515)

**Purpose:** Turn the downwell data and geometry into the targets and conceptual model an engineer or agent needs to calibrate the natural-state model.
**Sources:** `Baseline/Teaching_Model_3D/Data/Data_json/WELL-*_data.json` (stable T/P surveys dated 29/07/2015, feed zones, tracks, casing), `Models/TE14515/gTE14515.dat` (grid columns, layers and surface only), and `Data/Mapping/Mapping.json` (lake, river, faults).
**Not used as data:** the reference model's rock types, boundary conditions, upflows and NS results.
**Script:** `Scripts/natural_state_analysis.py`. Run it from the repository root and it regenerates everything in `Output/`.

---

## 1. Key conclusions (for calibration)

1. **The main upflow is in the NW, at the intersection of Fault 1 and Fault 3, beneath WELL-2 (column `ajf`) and WELL-10 (column `arz`).**
   - WELL-2 is the only well whose temperature rises again at depth: 236 °C at −1500 mRL, 254 °C at TD (−2348 mRL).
   - WELL-2 also shows an upward pressure gradient about 10 % above hydrostatic all the way from −450 to −2348 mRL. This points to deep hot upflow at about 255–270 °C or hotter.
2. **There is a shallow boiling / two-phase zone at WELL-10.** Below about −150 mRL the temperature is within 3 °C of Tsat(P), reaching 268.7 °C at −250 to −350 mRL. WELL-10's produced enthalpy of about 2680–2710 kJ/kg (steam) confirms this. The model must be able to form a small two-phase zone or steam cap at roughly −150 to −350 mRL here.
3. **A hot lateral outflow aquifer lies at about −250 to −900 mRL (model layers 54–60).** It spreads SE (WELL-6, WELL-7), E (WELL-3, then WELL-5) and NE (WELL-9).
   - Every one of these wells shows a temperature reversal, with cooler fluid underneath (see table 2).
   - Peak temperature falls with distance from WELL-2/10: 267 → 261 (W3) → 220 (W6) → 218 (W9) → 177 (W7) → 158 (W5).
4. **Cold recharge comes in from the north and SW.**
   - WELL-4 (north, high ground, wellhead 340 mRL) is near-isothermal at 12–29 °C to −750 mRL. This means strong cold downflow.
   - WELL-8 (SW) is near-isothermal at about 64–68 °C from −650 to −1500 mRL. This points to lateral cool inflow at depth.
5. **The eastern margin (WELL-1) is conductive.** Temperature rises steadily to 118.6 °C at −1672 mRL, the gradient is about 70 °C/km in the top 1 km and about 16 °C/km deeper, and there is no reversal. WELL-1 also has the lowest water level and the lowest pressure at −750 mRL (89.0 bar), next to the river. The natural pressure gradient drives outflow E/SE toward the river.
6. **A low-permeability clay cap is needed over the hot area.** Hot wells show very steep conductive gradients in the top 200–400 m. For example, WELL-3 goes from 35 to 206 °C between 164 and −25 mRL, and WELL-10 from 25 to 216 °C between 210 and 10 mRL. The cap base is about −150 to −250 mRL in the centre and thins at WELL-10.
7. **The profiles can be treated as natural state.** The surveys are dated 2015, after production started in 2014, but the downhole pressures match the first (2014, pre-production) transient pressures:

   | Well | Elevation (mRL) | Pressure, profile vs 2014 transient |
   |---|---|---|
   | WELL-2 | −750 | 98.2 vs 98.3 bar |
   | WELL-3 | −550 | 80.7 vs 80.8 bar |
   | WELL-9 | −675 | 91.2 vs 89.2 bar (2 bar higher; treat with lower confidence) |
   | WELL-10 | −250 | 56.7 bar |

   WELL-10's profile was measured before its own production started in 2019. All profiles are flagged `StableFlag: True`.

## 2. Well-by-well summary

| Well | Role | Column | Wellhead (mRL) | Tmax (°C) at (mRL) | T at TD (°C) | Water level est. (mRL) | Regime |
|---|---|---|---|---|---|---|---|
| WELL-10 | producer | arz | 260 | 268.7 @ −250 | 268.5 (−350) | 186 | Boiling below −150 mRL; two-phase / steam (upflow) |
| WELL-2 | producer | ajf | 225 | 266.9 @ −550 | 253.7 (−2348) | 204 | Upflow: reversal of 31 °C, then reheating at depth |
| WELL-3 | producer | asd | 164 | 261.4 @ −450 | 188.2 (−2887) | 148 | Strong outflow over cool fluid (87 °C reversal) |
| WELL-6 | monitor | aph | 147 | 220.5 @ −550 | 126.1 (−2359) | 144 | Strong outflow (98 °C reversal) |
| WELL-9 | producer | ajj | 174 | 217.7 @ −750 | 161.6 (−2427) | 160 | Outflow (58 °C reversal), deviated to SE |
| WELL-7 | monitor | apk | 190 | 176.7 @ −650 | 144.7 (−2762) | 153 | Outflow (46 °C reversal) |
| WELL-5 | injector | arc | 93 | 158.0 @ −650 | 154.7 (−2595) | 127 | Distal outflow (13 °C reversal), conductive at depth |
| WELL-1 | monitor | arg | 15 | 118.6 @ TD | 118.6 (−1672) | 108 | Conductive E margin |
| WELL-8 | injector | ana | 176 | 93.1 @ TD | 93.1 (−2716) | 176 | Cool; isothermal 64–68 °C (−650 to −1500), lateral cool inflow |
| WELL-4 | injector | aop | 340 | 29.0 @ TD | 29.0 (−750) | 228 | Cold downflow / recharge |

The water levels are extrapolated from the upper pressure profile to 1 bar. In hot wells they are biased low because the hot water column is less dense.

## 3. Pressure

**Pressure at common elevations (bar):**

| Elevation (mRL) | W1 | W2 | W3 | W4 | W5 | W6 | W7 | W8 | W9 | W10 |
|---|---|---|---|---|---|---|---|---|---|---|
| −250 | 38.6 | 55.0 | 55.0 | 51.1 | 44.9 | 49.9 | 46.9 | 46.1 | 49.3 | 56.7 |
| −750 | 89.0 | 98.2 | 97.4 | 100.5 | 95.4 | 97.0 | 96.6 | 95.9 | 97.7 | – |
| −1500 | 162.0 | 164.5 | 163.9 | – | 164.0 | 165.1 | 164.8 | 169.2 | 163.9 | – |
| −1900 | – | 199.8 | 200.4 | – | 200.5 | 202.1 | 201.8 | 208.6 | 200.3 | – |

**What the pressures show:**
- **At −250 mRL** pressure is highest in the upflow (W10, W2, W3, about 55–57 bar) and lowest at W1 (38.6 bar). Shallow fluid flows outward from the upflow.
- **At −750 mRL** the highest pressure is at the northern recharge well W4 (100.5 bar) and the lowest at the eastern margin W1 (89.0 bar). This is the regional N/NW → E/SE gradient.
- **At depth** W8 has the highest pressure, because its column is cold and dense. W2 has the lowest, because its column is hot.
- **Ratio of measured dP/dz to ρ_liquid(T)·g** (`hydrostatic_check.csv`):
  - It is about 1.00–1.03 in most wells below −600 mRL, which is near-hydrostatic.
  - It is about 1.10 at WELL-2 throughout, which indicates upflow.
  - It is 1.2–1.8 in the top 300 m of the hot wells, which indicates overpressure beneath the cap.
  - It is 0.95 at the bottom of WELL-4, which indicates downflow.

**Calibration target:** match the feed-zone pressures to within about 1–2 bar, and reproduce the E–W pressure contrast at −750 mRL (W1 about 9 bar below the centre).

## 4. Conceptual model for the natural state

```
 N (W4: cold recharge, high water table)
        \  downflow
         v
  W9 ─── outflow 218°C @ -750 ──►  NE
   ^
 [F1 x F3]  W2/W10 UPFLOW (≥255–270°C at depth; boiling/steam at W10 −150 to −350 mRL)
   │  clay cap to ~ −150/−250 mRL
   ├──► W3 (261) ──► W5 (158) ──► W1 (conductive, low P, river)    E
   └──► W6 (220) ──► W7 (177)                                       SE
 cool fluid (120–180 °C) beneath outflow below about −1100 mRL
 W8 (SW): cool lateral inflow, 64–68 °C
```

## 5. Recommended calibration set-up

**Targets:**
- `NS_calibration_targets.csv` maps every T/P point to a model block (column + layer) and gives a suggested weight:
  - 1.0 at feed zones
  - 0.8 in open hole
  - 0.5 behind casing
  - 0.3 near the wellhead
  - 0 above the model surface
- `feedzones.csv` lists the 18 feed zones with their T, P and block. Use these as the primary targets.
- `T_P_at_layer_centres.csv` gives T and P interpolated at every layer centre for direct block-by-block comparison.

**Heat/mass input:**
- Concentrate deep upflow in columns around `ajf` and `arz` (F1/F3 intersection).
- Use a source enthalpy consistent with about 270–280 °C liquid (roughly 1180–1240 kJ/kg). W10's near-boiling 268.7 °C sets the lower bound.
- Use a lower background conductive heat flux elsewhere. WELL-1's deep gradient of about 16 °C/km fits a low-flux margin.

**Permeability:**
- High horizontal permeability in layers 54–60 (about −250 to −900 mRL), extending SE, E and NE from the upflow.
- Low permeability in the cap (top layers down to about −150/−250 mRL) over the hot area.
- Moderate/low permeability below −1100 mRL away from W2, so the cool under-flow survives beneath the outflow.
- Permeable downflow paths at W4 (N) and lateral cool inflow at depth near W8 (SW).

**Boundaries:**
- The top boundary must allow the observed water levels (about 110–230 mRL).
- Allow outflow/discharge toward the E/SE and the river near W1.

## 6. Data quality and issues to resolve

- **The model surface is above the wellheads in several columns.** Misfits are: WELL-1 −105 m (surface 120 vs wellhead 15 mRL), WELL-8 −84 m, WELL-7 −50 m, WELL-5 −47 m, WELL-9 −46 m, WELL-2 −35 m. Shallow model temperatures in these columns will be biased. Either correct the SURFA elevations or down-weight the near-surface data.
- **WELL-4's casing record is inconsistent.** It lists a liner to −3251 mRL, but the track and survey end at −750 mRL. The survey also starts at 230 mRL, while the wellhead is at 340 mRL. Treat −750 mRL as TD.
- **WELL-10's feed-zone quality is "Synthetic",** and the well is short (TD −350 mRL). Its boiling-point profile makes it highly diagnostic all the same.
- **WELL-9 is deviated.** Its second feed zone (−837.5 mRL) falls in column `ajq`, not the wellhead column `ajj`. The script maps each point along the track. The same applies to WELL-8.
- **WELL-9's pressure profile is about 2 bar above** its 2014 transient value at −675 mRL. The other wells agree within 0.1 bar.
- **The geometry used includes the WELLS section in `gTE14515.dat`.** It was committed as "Fixed well geometry" (7b0992b) during this analysis. This analysis only read the Baseline folder.

## 7. Output files

| File | Contents |
|---|---|
| `well_summary.csv` | Per-well summary (location, column, surface misfit, Tmax, gradients, water level, boiling margin, regime) |
| `NS_calibration_targets.csv` | Every T/P datum with Tsat, boiling margin, model column/layer/block and suggested weight |
| `feedzones.csv` | Feed zones with T, P and model block |
| `T_P_at_layer_centres.csv` | T and P per well interpolated to the model layer centres |
| `pressure_consistency_check.csv` | Profile P vs first transient P (natural-state check) |
| `hydrostatic_check.csv` | Interval dP/dz vs liquid hydrostatic (vertical-flow indicator) |
| `figures/01_well_map.png` | Wells, grid, lake/river, faults |
| `figures/02_T_P_profiles.png` | All T and P profiles, with Tsat |
| `figures/03_well_panels.png` | Per-well T vs boiling curve, feed zones, casing, layer boundaries |
| `figures/04_T_maps_by_elevation.png` | Temperature maps at 0 to −2000 mRL |
| `figures/05_P_maps_by_elevation.png` | Pressure at −500, −1000, −1500 mRL |
| `figures/06_section_*.png` | Interpolated isotherm cross-sections (indicative only), plus the section-line map |
