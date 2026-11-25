# Standard figures for natural-state data analysis: selection from O'Sullivan et al. (2025)

**Source:** O'Sullivan, J., Gravatt, M., Riffault, J., Renaud, T., O'Sullivan, M., Castillo Ruiz, N., Ayling, B., Mannington, W. (2025). *An updated model of Ohaaki geothermal field, New Zealand.* Geothermics 130, 103339 (CC BY 4.0).

**Purpose.** This note picks the figure types from the paper that are worth standardising as outputs of the data-analysis step, explains why, and states how each one is adapted for data-only use. The resulting specification is §8 of the proposed [CLAUDE_v4_proposed.md](CLAUDE_v4_proposed.md).

## 1. How the paper's figures were judged

Each of the paper's 27 figures was scored against five criteria that follow from CLAUDE.md and from the nine-run comparison:

1. **Data-only.** It can be drawn from the data roles alone, without the reference model.
2. **In scope.** It describes the natural state, or the §4.1 natural-state test.
3. **Calibration value.** An engineer reads a target, a structure or a data problem from it.
4. **Specifiable.** Its layout can be written as rules an agent will follow the same way every time. This is what makes runs comparable.
5. **Honest with sparse data.** It does not invent values between wells. This was a weakness in iteration 3: Opus contoured below its data, and Haiku's sketch carried no structure.

## 2. Catalogue

| Fig. | Content | Data-only | In scope | Decision |
|---|---|---|---|---|
| 1 | Location map: system boundary, river, wellheads, tracks | yes | yes | **Selected → F1** |
| 2 | Annual production and generation | yes | no | Out of scope |
| 3 | Model grid on topography, with section line A–A′ | geometry | yes | **Selected → F1** (grid and section line) |
| 4 | Geological cross-section with feed zones and upflow arrows | needs geology | yes | Symbology used in F6 |
| 5 | Faults in the conceptual model, model extent | yes | yes | **Selected → F1** (faults and extent) |
| 6 | Conceptual model on a vertical slice (two-phase zone, liquid and vapour upflow) | interpretation | yes | **Selected → F6** |
| 7 | 3D alteration and resistivity view | no (MT/DC data) | yes | Not reproducible from these data |
| 8 | Rock types in the bottom layer; 3D faults | reference | — | Reference model only |
| 9–12 | Permeability and porosity on A–A′ | reference | — | Reference only; layout reused in §8.3 |
| 13 | Mass upflow by block, plan view | reference | — | Reference only |
| 14 | Natural-state temperature on A–A′, with flow arrows and surface features labelled on top | model result | yes | **Selected → F4 layout** (data version) |
| 15–16 | Profile small multiples: model line against data dots, inset locator map, shared axes, survey-date legend | yes (data half) | yes | **Selected → F2** |
| 17 | Bubble map of temperature error in an elevation band, in red/blue size classes | yes (adapted) | yes | **Selected → F3** |
| 18 | Natural-state results on geology, with recharge and upflow arrows | model result | yes | Arrow symbology used in F6 |
| 19, 21 | Production flow and enthalpy against time | — | no | Out of scope |
| 20, 23 | Feed-zone T and P against time, with locator map and production shading | yes | partly | **Layout selected → F7** (gauge check) |
| 22 | Pressure decline against time at one elevation, with locator map | yes | partly | Layout used in F7 |
| 24 | A–A′ at five dates | model result | no | Out of scope |
| 25 | Conceptual model 2017 (well types, inflow, cooling, two-phase, isotherms) | interpretation | yes | Legend vocabulary used in F6 |
| 26–27 | Surface CO₂, heat flux and temperature against surface features | needs surface data | yes | Only if the dataset holds surface measurements |

## 3. The selected figures and why

### F1: Base map (from Figs 1, 3 and 5)
The paper uses one base map in every plan-view figure (Figs 1, 5, 13, 15–17, 20–23, 26, 27): the same system boundary, river and well tracks. A reader orients once and then reads everything else against it. Fig. 3 adds the section line A–A′, which ties every later section to the map.
**Why standardise it:** in iteration 3, Haiku's map left out the grid and the three runs drew different layers. One fixed base layer that carries the section lines removes that variation. The checker can confirm it exists, and the reviewer can check labels quickly.

### F2: Profile panels with a locator map (from Figs 15 and 16)
These are the paper's main natural-state calibration figures. Small multiples share one temperature axis (0–350 °C) and one elevation axis. A locator map highlights each well in the colour of its panel, and a legend carries the survey date. Wells are grouped by area, West Bank and East Bank.
**Why standardise it:**
- Shared axes make wells comparable at a glance: hotter, shallower and inverted profiles stand out without reading numbers.
- Grouping by area maps directly onto the `well_groups` input role.
- It is the same figure the engineer will produce after calibration with model lines added, so the data-analysis figure becomes the "data" half of the calibration plot.
- The per-well profile (already required) gains the locator inset and fixed axes. The group panel is new.

### F3: Elevation-band symbol map (from Fig. 17)
Fig. 17 shows, for one elevation band, each well's temperature error as a red or blue circle sized by class, with a cross for "within 10 °C".
**Why standardise it:**
- **It does not interpolate.** With about ten wells, a contoured plan map is mostly invented. A symbol at each well shows only what was measured, which honours "interpolation only, never extrapolate" more strictly than contouring.
- Fixed classes (±10/20/30/40 °C) give the same picture whichever model draws it.
- For data analysis, plotting the deviation from the field median at that depth shows lateral hot and cold zones.
- **The same template becomes the residual map in calibration.** Replace "minus field median" with "model minus data" and it is Fig. 17. That makes the hand-off to calibration seamless.
- It replaces the contoured `maps_temperature*` figures. The file names still match the checker's glob.

### F4: Data section through the hottest wells (from Fig. 14, with Figs 9–12 as layout)
The paper's sections all share one line, A–A′. Surface features (river, springs, fumaroles) are labelled along the top where they cross the line, well tracks are labelled at their bottoms, and the clay-cap outline is overlaid.
**Why standardise it:**
- A fixed line, the same axes and features labelled at their crossings make a section quick to read and comparable across runs.
- The paper's figure is a model result filled everywhere. The data version must contour only inside the data envelope, the rule iteration 3 broke. The spec makes the envelope explicit.
- Stippling near-boiling intervals carries over Fig. 6's two-phase symbol.
- The same line is reused in F6 and in the optional reference-model companion (§8.3, modelled on Figs 9–12). The engineer can lay data, interpretation and starting model over one another.

### F6: Conceptual model (from Figs 6, 14, 18 and 25)
The paper's conceptual figures share one vocabulary: red arrows for deep upflow, blue for meteoric recharge, stipple for the two-phase zone, an outlined cap, and labelled isotherms. Fig. 25's legend sets this out.
**Why standardise it:**
- A fixed legend makes the elements comparable between runs. In iteration 3 one sketch was regime colours only, and the others used different symbols.
- Drawing it on the F4 line, over the data isotherms, ties each interpretation to the data behind it.
- **Solid against dashed arrows and numbered elements** keyed to a table in the report enforce the `[obs]`/`[interp]` distinction CLAUDE.md asks for. In the paper this is implicit; for agents it needs to be explicit.
- A plan-view panel shows upflow location and outflow direction against faults and water bodies, which CLAUDE.md requires.

### F7: Gauge check (layout from Figs 20 and 22)
These panels pair a locator map with time series. They overlay data dots and grey shading for production periods, and print elevation and well name in the panel.
**Why standardise it:** the gauge test decides each well's natural-state verdict and therefore its weights, yet in iteration 3 it existed only as a table. A figure that shows the first reading, the survey value, the precision band and the production start makes the verdict checkable by eye. It also shows straight away when a gauge started after production, the case that split Haiku's confidence grading from Sonnet's and Opus's. Pressure history is in scope here only as this test (§4.1).

### Kept from v3: P–T state (F5)
The paper has no P–T plot, but it remains the clearest view of how close each well is to boiling, so it is kept unchanged.

## 4. Shared conventions taken from the paper

| Convention in the paper | Where it appears | Standard (§8.1) |
|---|---|---|
| One base map under every plan view | Figs 1, 5, 13, 15–17, 20–23 | One base layer |
| Locator inset beside every per-well plot | Figs 15–16, 20–23 | F2 and F7 |
| Shared axes across panels | Figs 15–16 | Fixed elevation and temperature axes |
| Data as dots, model as lines | Figs 15, 16, 19–23 | Data are points; derived curves are lines |
| Red hotter, blue colder, fixed size classes | Fig. 17 | Diverging scale and classes in F3 |
| One section line reused everywhere | Figs 3, 9–12, 14, 18 | F1, F4 and F6 share the line |
| Surface features labelled on top of sections | Figs 9–14 | F4 crossing ticks |
| Fixed flow-arrow colours | Figs 14, 18, 25 | F6 symbology |
| Production periods shaded | Figs 20, 23 | F7 shading |

**Where the standard deliberately differs from the paper:**
- **Colour scale.** The paper uses a rainbow temperature scale. The standard asks for one perceptually uniform sequential scale, so equal colour steps read as equal temperature steps.
- **Data envelope.** The paper's sections are filled everywhere because they show model output. Data sections must stay inside the data.
- **Observed against inferred.** The paper does not mark which arrows are observed and which inferred. The standard does, because the agent's audience must know which claims are data.

## 5. Matching checker changes (for you to make; agents must not edit the checker)

| Check | Change |
|---|---|
| F1 | Add `gauge_check*.png` (1) and `profile_panels_*.png` (1). Count `profile_*.png` excluding `profile_panels_*` for the per-well rule, or the group panels could hide a missing well. |
| F2 (new) | `figures_index.csv` exists, every standard ID F1–F7 appears, and every listed file exists. |
| G4 (new) | `model_surface_mRL` matches the surface section. |
| D1 (new) | Recompute the §4.4 definitions (shallow and deep gradients, water level, minimum boiling margin, isotherm crossings, inversion, regime) from `NS_calibration_targets.csv` within tolerance. |
| W3 (new) | Weight equals location factor × confidence factor; `ns_confidence` agrees across `well_summary.csv`, `pressure_consistency_check.csv`, the targets and `feedzones.csv`. |
| R4 (new) | No `{identifier}` or `nan` text in the report. |
| S3 (new) | LF line endings, integer `feedzone`, ISO dates. |
| — | Add `T_P_at_elevation_bands.csv` and `figures_index.csv` to REQUIRED_FILES and REQUIRED_COLUMNS. |

Figure content (envelope, axes, symbology) cannot be checked automatically from the PNG. `figures_index.csv` and the per-figure checklists in §8.2 give the agent, and a reviewer, a fixed list to check against.
