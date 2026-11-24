"""Natural-state data interpretation for geothermal model calibration (CLAUDE.md v3).

Run from the repository root:

    python Scripts/natural_state_analysis.py

Regenerates every output in Output/ except Output/process_log.md (hand-written).
Thermodynamics: IAPWS-IF97 via the `iapws` package (Region 4 for Tsat, Region 1 for liquid density).
All numbers in NS_data_interpretation.md are written from the output tables produced here.
"""
import ast
import glob
import json
import os
import re
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from iapws import IAPWS97  # noqa: E402
from scipy.interpolate import griddata  # noqa: E402
from scipy.optimize import brentq  # noqa: E402

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE_REL = "Baseline/Teaching_Model_3D"
TM = os.path.join(ROOT, BASE_REL)
OUT = os.path.join(ROOT, "Output")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

# ----------------------------------------------------------------------------------------
# Parameters (all written to Output/analysis_parameters.csv)
# ----------------------------------------------------------------------------------------
PAR = dict(
    tsat_tol_K=0.1, bm_flag_C=5.0, bm_liquid_min_C=10.0, p_liquid_min_bar=1.5, p_critical_bar=220.64,
    window_min_m=300.0, iso_span_m=300.0, iso_range_K=2.0, inv_min_C=3.0,
    shallow_span_m=300.0, deep_min_m=300.0, wellhead_top_m=30.0, cold_Tmax_C=100.0,
    gauge_read_prec_bar=0.1, gauge_tol_bar=0.2, trend_conf_bar=2.0, trend_tol_bar=3.0,
    atm_bar=1.01325, g=9.80665, rounding_unc_dP_bar=0.1,
)
CONF_FACTOR = {"high": 1.0, "medium": 0.75, "low": 0.5}
LOC_WEIGHT = dict(feedzone=1.0, open=0.8, iso_open=0.5, casing=0.5, top=0.3, above=0.0)
QUALITY_RANK = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}
ELEV_COMMON = [200, 0, -500, -1000, -1500, -2000]
MAP_LEVELS = [0, -500, -1000, -1500]
SHORT = lambda n: n.replace("WELL-", "W")  # noqa: E731


# ----------------------------------------------------------------------------------------
# Thermodynamics (IAPWS-IF97)
# ----------------------------------------------------------------------------------------
def tsat(p_bar):
    """Saturation temperature (degC), IF97 Region 4. NaN above the critical pressure (no Tsat)."""
    if p_bar is None or not np.isfinite(p_bar) or p_bar <= 0.0061 or p_bar >= PAR["p_critical_bar"]:
        return np.nan
    return IAPWS97(P=p_bar / 10.0, x=0).T - 273.15


def rho_liq(p_bar, t_c):
    """Liquid density (kg/m3), IF97 Region 1."""
    return IAPWS97(P=p_bar / 10.0, T=t_c + 273.15).rho


def validate_if97():
    rows = []
    for p_bar, t_ref_k in [(1.0, 372.755919), (10.0, 453.035632), (100.0, 584.149488)]:
        t_k = IAPWS97(P=p_bar / 10.0, x=0).T
        rows.append((p_bar, t_k - 273.15, t_k - t_ref_k))
        assert abs(t_k - t_ref_k) < PAR["tsat_tol_K"], f"IF97 Tsat error at {p_bar} bar"
    assert abs(tsat(1.01325) - 99.97) < 0.1 and abs(tsat(100.0) - 311.0) < 0.1
    return rows


def injectate_temperature(h_kjkg):
    """Temperature of compressed liquid with the given enthalpy at 10 bar (IF97 Region 1)."""
    return brentq(lambda t: IAPWS97(P=1.0, T=t + 273.15).h - h_kjkg, 5.0, 150.0)


# ----------------------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------------------
def lit(s):
    return ast.literal_eval(s) if isinstance(s, str) else s


def interp(z_asc, v_asc, z):
    """Linear interpolation; NaN outside the surveyed range (never extrapolates)."""
    if z is None or not np.isfinite(z):
        return np.nan
    if z < z_asc[0] - 1e-9 or z > z_asc[-1] + 1e-9:
        return np.nan
    return float(np.interp(z, z_asc, v_asc))


def dist_to_polyline(px, py, segs):
    best = np.inf
    p = np.array([px, py], float)
    for s in segs:
        if len(s) == 1:
            best = min(best, float(np.hypot(*(s[0] - p))))
            continue
        for a, b in zip(s[:-1], s[1:]):
            ab = b - a
            t = np.clip(np.dot(p - a, ab) / max(np.dot(ab, ab), 1e-12), 0.0, 1.0)
            best = min(best, float(np.hypot(*(a + t * ab - p))))
    return best


def first_crossing(z_desc, t_desc, thr):
    """Elevation of the first crossing of an isotherm, scanning downward from the top survey point."""
    for i in range(len(z_desc) - 1):
        a = t_desc[i] - thr
        b = t_desc[i + 1] - thr
        if a == 0:
            return float(z_desc[i])
        if a * b < 0:
            return float(z_desc[i] + (z_desc[i + 1] - z_desc[i]) * (thr - t_desc[i]) / (t_desc[i + 1] - t_desc[i]))
    return np.nan


def rnd(x, d):
    if x is None:
        return None
    try:
        if not np.isfinite(x):
            return None
    except TypeError:
        return x
    return round(float(x), d)


def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(str(c) for c in cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        vals = []
        for c in cols:
            v = r[c]
            if v is None or (isinstance(v, float) and not np.isfinite(v)):
                vals.append("–")
            else:
                vals.append(str(v))
        lines.append("| " + " | ".join(vals) + " |")
    return "\n".join(lines)


def fmt(x, d=1):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "–"
    return f"{x:.{d}f}"


def is_tf(v):
    return str(v).strip().lower() in ("true", "1", "yes")


# ----------------------------------------------------------------------------------------
# 1. Input inventory (CLAUDE.md 2.1)
# ----------------------------------------------------------------------------------------
def role_for(rel):
    """rel is relative to Baseline/Teaching_Model_3D, forward slashes."""
    if rel.startswith("Data/Data_json/"):
        return "well_data_json", "Per-well JSON: track, design, feed zones, T/P survey, transients, flows (factual)"
    if rel.startswith("Data/Data_dat/"):
        return "well_data_dat", "Plain-column copy of well data; cross-checked against JSON, not a primary source"
    if rel == "Data/Mapping/Mapping.json":
        return "mapping", "Faults, lake and river polylines (factual); river entries are ragged arrays"
    if rel == "Data/Area_data.json":
        return "well_groups", "Area membership: list of wells (factual)"
    if rel == "Data/Groups_data.json":
        return "out_of_scope", "Field-total steam-flow record (production); not used in any calculation"
    if rel.startswith("Data/Plotting/"):
        return "reference_other", "Plotting configuration; not evidence"
    if rel == "Models/TE14515/gTE14515.dat":
        return "geometry_dat", "MULgraph grid; only LAYERS used for layer centres (geometry only)"
    if rel == "Models/TE14515/gTE14515_wells.json":
        return "well_blocks_json", "Well-to-block intersections, feed-zone blocks and interpolation weights; block authority"
    if rel == "Models/TE14515/gTE14515.auxjson":
        return "reference_other", "Scene/plot features copied from Mapping.json; hard-coded path from another machine"
    if rel == "Models/TE14515/gTE14515_boundary_conditions.json":
        return "reference_bc", "Reference-model boundary conditions; reference only, not read"
    if rel == "Models/TE14515/gTE14515_dict_rocktype.json":
        return "reference_rock", "Reference-model rock types; reference only, not read"
    if rel == "Models/TE14515/currentGeomSetup.json":
        return "reference_other", "Model set-up (hard-coded paths from another machine); not read"
    if rel in ("Models/TE14515/NaturalState/TE14515_NS_000.incon", "Models/TE14515/NaturalState/TE14515_NS_000.dat"):
        return "reference_results", "Reference natural-state run output / initial conditions; role from file name, not opened"
    if rel == "Models/TE14515/NaturalState/currentNSSetup.json":
        return "reference_other", "Reference natural-state set-up; not read"
    if rel.startswith("Models/TE14515/Production/") or rel.startswith("Models/TE14515/FutureScenarios/"):
        return "out_of_scope", "Production or future-scenario model set-up; not read"
    if rel.startswith("Parameters/"):
        return "out_of_scope", "Production-history, reinjection and forecast spreadsheets; not read (records used only via well files)"
    if rel.endswith(".DS_Store"):
        return "out_of_scope", "macOS metadata"
    return None, "UNASSIGNED"


def write_inventory():
    rows, unassigned = [], []
    for dirpath, _, files in os.walk(TM):
        for f in sorted(files):
            full = os.path.join(dirpath, f)
            rel = os.path.relpath(full, TM).replace("\\", "/")
            role, note = role_for(rel)
            if role is None:
                unassigned.append(rel)
            rows.append(dict(role=role or "UNASSIGNED", path=f"{BASE_REL}/{rel}", note=note))
    inv = pd.DataFrame(rows).sort_values(["role", "path"]).reset_index(drop=True)
    inv.to_csv(os.path.join(OUT, "input_inventory.csv"), index=False)
    return inv, unassigned


# ----------------------------------------------------------------------------------------
# 2. Load data: JSON primary, block file, grid, mapping
# ----------------------------------------------------------------------------------------
def read_layers(path):
    lines = open(path, encoding="utf-8", errors="ignore").read().splitlines()
    i0 = next(k for k, ln in enumerate(lines) if ln.strip() == "LAYERS")
    L = {}
    for ln in lines[i0 + 1:]:
        p = ln.split()
        if len(p) != 3:
            break
        try:
            L[int(p[0])] = (float(p[1]), float(p[2]))
        except ValueError:
            break
    layers = {i: dict(top=L[i - 1][0], bottom=L[i][0], centre=L[i][1]) for i in sorted(L) if i - 1 in L}
    gaps = [i for i in layers if i + 1 in layers and abs(layers[i]["bottom"] - layers[i + 1]["top"]) > 1e-6]
    assert not gaps, f"layer gaps at {gaps[:5]}"
    return layers, len(L)


def parse_block_name(b):
    m = re.match(r"^([A-Za-z]+)(\d+)$", b)
    return m.group(1), int(m.group(2))


def load_wells():
    wb = json.load(open(os.path.join(TM, "Models/TE14515/gTE14515_wells.json"), encoding="utf-8"))
    wells = {}
    for fp in sorted(glob.glob(os.path.join(TM, "Data", "Data_json", "*.json"))):
        d = json.load(open(fp, encoding="utf-8"))
        name = d["Well_Name"]
        w = dict(name=name, n=int(re.search(r"(\d+)$", name).group(1)), aut2=d.get("AUT2_Name"))
        tr = np.array([lit(s) for s in d["Well_Track"]["Data"]], float)
        w["track"] = tr
        w["wellhead_xy"] = tr[0, :2]
        w["wh"] = float(tr[0, 2])
        o = np.argsort(tr[:, 2])
        w["track_z_asc"], w["track_x_asc"], w["track_y_asc"] = tr[o, 2], tr[o, 0], tr[o, 1]
        des = [lit(s) for s in d["Well_Design"]["Data"]]
        w["design"] = [dict(top=float(a), bottom=float(b), diam=str(c), slotted=bool(s), spec=str(sp))
                       for a, b, c, s, sp in des]
        w["design_source"] = d["Well_Design"].get("Source", "")
        w["shoe"] = min(r["bottom"] for r in w["design"] if not r["slotted"])
        fz = [lit(s) for s in d["Feed_Zones"]["Data"]]
        w["fz"] = [dict(idx=i + 1, top=max(a, b), bottom=min(a, b), q=str(q)) for i, (a, b, q) in enumerate(fz)]
        w["fz_source"] = d["Feed_Zones"].get("Source", "")
        dates = list(d["Downhole_Temperature"].keys())
        assert len(dates) == 1, f"{name}: expected one survey date"
        dt = dates[0]
        Tn, Pn = d["Downhole_Temperature"][dt], d["Downhole_Pressure"][dt]
        assert set(Tn["Data"]) == set(Pn["Data"]), f"{name}: T and P survey depths differ"
        zs = np.array([float(k) for k in Tn["Data"]])
        Ts = np.array([Tn["Data"][k] for k in Tn["Data"]], float)
        Ps = np.array([Pn["Data"][k] for k in Tn["Data"]], float)
        o = np.argsort(zs)
        w["survey_date"] = pd.to_datetime(dt, format="%d/%m/%Y")
        w["zA"], w["TA"], w["PA"] = zs[o], Ts[o], Ps[o]
        w["survey_meta"] = dict(T_source=Tn.get("Source", ""), T_comment=Tn.get("Comment", ""),
                                P_source=Pn.get("Source", ""), P_comment=Pn.get("Comment", ""),
                                stable=str(Tn.get("StableFlag", "")))
        w["gauges"] = {}
        for k, v in d.get("Transient_Pressure", {}).items():
            idx = [pd.to_datetime(x, format="%d/%m/%Y") for x in v["Data"]]
            w["gauges"][float(k)] = dict(series=pd.Series(list(v["Data"].values()), index=idx).sort_index())

        def series(key):
            if key not in d:
                return None
            dd = d[key]["Data"]
            idx = [pd.to_datetime(x, format="%d/%m/%Y") for x in dd]
            return pd.Series(list(dd.values()), index=idx).sort_index()

        w["prod"] = series("Production_Mass_Flow")
        w["inj"] = series("Injection_Mass_Flow")
        if "Transient_Enthalpy" in d:
            w["enth"] = series("Transient_Enthalpy")
        else:
            w["enth"] = series("Injection_Enthalpy")
        w["flow_keys"] = [k for k in ("Production_Mass_Flow", "Injection_Mass_Flow", "Transient_Enthalpy", "Injection_Enthalpy") if k in d]
        b = wb[name]
        w["blocks"] = [dict(block=x["block"], top=max(x["interval"]), bottom=min(x["interval"]),
                            layer=parse_block_name(x["block"])[1], column=parse_block_name(x["block"])[0])
                       for x in b["Blocks"]]
        w["fz_blocks"] = [list(v.keys()) for v in b["Feed_Zones"].values()]
        iw = b.get("Interp_Weights", {})
        iw = json.loads(iw) if isinstance(iw, str) else iw
        w["iw"] = {round(float(k), 3): v for k, v in iw.items()}
        wells[name] = w
    return wells, wb


def load_mapping():
    mp = json.load(open(os.path.join(TM, "Data/Mapping/Mapping.json"), encoding="utf-8"))
    feats = {}
    for k, v in mp.items():
        arr = v["Data"]
        if len(arr) and isinstance(arr[0][0], list):
            segs = [np.array(s, float) for s in arr]
        else:
            segs = [np.array(arr, float)]
        feats[k] = dict(segs=segs, kind="fault" if k.startswith("Fault") else "water")
    return feats


def first_flow(w):
    dates = [s.index.min() for s in (w["prod"], w["inj"], w["enth"]) if s is not None and len(s)]
    return min(dates) if dates else pd.NaT


# ----------------------------------------------------------------------------------------
# 3. Cross-checks: JSON vs plain-column dat, provenance
# ----------------------------------------------------------------------------------------
def crosscheck(wells):
    rows = []
    dat_dir = os.path.join(TM, "Data", "Data_dat")
    ref = pd.Timestamp("2014-01-01")  # transient and flow time axis: years since this date (matched below)
    for w in wells.values():
        n = w["n"]
        ids = ("Well%d" % n, "Wel%d" % n)
        for kind, col in (("downhole_temperature", "T"), ("downhole_pressure", "P")):
            fp = [f for f in glob.glob(os.path.join(dat_dir, f"*_{kind}.dat"))
                  if os.path.basename(f).split("_")[0] in ids]
            assert len(fp) == 1, (w["name"], kind, fp)
            a = np.loadtxt(fp[0])
            jz, jv = w["zA"], (w["TA"] if col == "T" else w["PA"])
            diffs = []
            for val, z in a:
                k = int(np.argmin(np.abs(jz - z)))
                diffs.append(abs(jv[k] - val) if abs(jz[k] - z) < 1e-6 else np.nan)
            diffs = np.array(diffs)
            rows.append(dict(file=os.path.basename(fp[0]), series=kind, n_points=len(a),
                             max_abs_diff=rnd(np.nanmax(diffs), 6), n_unmatched=int(np.isnan(diffs).sum()),
                             status="PASS" if (not np.isnan(diffs).any() and np.nanmax(diffs) < 1e-6) else "CHECK",
                             note="depth matched to 1e-6 m; value matched to 1e-6"))
        for fp in sorted(glob.glob(os.path.join(dat_dir, "*.dat"))):
            base = os.path.basename(fp)
            if base.split("_")[0] not in ids or "downhole" in base:
                continue
            a = np.loadtxt(fp)
            if "transient_pressure" in base:
                ser = next(iter(w["gauges"].values()))["series"] if w["gauges"] else None
                kind = "transient_pressure"
            elif "mass_flow" in base:
                ser = w["prod"] if "production" in base else w["inj"]
                kind = "mass_flow"
            else:
                ser = w["enth"]
                kind = "enthalpy"
            if ser is None:
                rows.append(dict(file=base, series=kind, n_points=len(a), max_abs_diff=None,
                                 n_unmatched=len(a), status="NO JSON SERIES", note=""))
                continue
            dates = [ref + pd.Timedelta(days=float(t) * 365.25) for t in a[:, 0]]
            diffs, miss, signs = [], 0, []
            for dte, val in zip(dates, a[:, 1]):
                k = int(np.argmin(np.abs((ser.index - dte).days)))
                if abs((ser.index[k] - dte).days) <= 3:
                    jv = float(ser.iloc[k])
                    if kind == "mass_flow":
                        signs.append(np.sign(val) * np.sign(jv))
                        diffs.append(abs(abs(jv) - abs(val)))  # compare magnitudes; sign convention reported below
                    else:
                        diffs.append(abs(jv - val))
                else:
                    miss += 1
            md = max(diffs) if diffs else np.nan
            note = ""
            if kind == "mass_flow" and signs:
                note = "dat stores flow with the opposite sign to JSON (outflow negative); magnitudes compared" \
                    if np.mean(signs) < 0 else "same sign as JSON"
            if kind in ("transient_pressure", "enthalpy"):
                note = "time axis: years since 01/01/2014, matched to JSON dates within 3 days"
            rows.append(dict(file=base, series=kind, n_points=len(a), max_abs_diff=rnd(md, 6),
                             n_unmatched=miss, status="PASS" if (miss == 0 and md < 1e-6) else "CHECK", note=note))
    return pd.DataFrame(rows)


def provenance(wells):
    rows = []
    for w in wells.values():
        fz_q = sorted({z["q"] for z in w["fz"]})
        rows.append(dict(well=w["name"], survey_date=w["survey_date"].strftime("%d/%m/%Y"),
                         T_source=w["survey_meta"]["T_source"] or "(blank)",
                         T_comment=w["survey_meta"]["T_comment"] or "(blank)",
                         P_source=w["survey_meta"]["P_source"] or "(blank)",
                         P_comment=w["survey_meta"]["P_comment"] or "(blank)",
                         StableFlag=w["survey_meta"]["stable"],
                         design_source=w["design_source"] or "(blank)",
                         feedzone_source=w["fz_source"] or "(blank)",
                         feedzone_quality_tags=";".join(fz_q),
                         gauges_m=";".join(f"{k:.1f}" for k in w["gauges"]) or "none",
                         flow_records=";".join(w["flow_keys"]) or "none"))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------
# 4. Per-well analysis
# ----------------------------------------------------------------------------------------
def analyse(wells, layers):
    for name, w in wells.items():
        z, T, P = w["zA"], w["TA"], w["PA"]
        Ts = np.array([tsat(p) for p in P])
        bm = Ts - T
        w.update(zd=z[::-1], Td=T[::-1], Pd=P[::-1], Tsd=Ts[::-1], bmd=bm[::-1])
        top_blk = max(w["blocks"], key=lambda b: b["top"])
        w["top_block"] = top_blk["block"]
        w["surface"] = layers[top_blk["layer"]]["top"]
        w["ztop"] = float(z[-1])
        w["TD"] = float(z[0])
        for i, f in enumerate(w["fz"]):
            f["mid"] = 0.5 * (f["top"] + f["bottom"])
            keys = w["fz_blocks"][i]
            if len(keys) == 1:
                f["block"] = keys[0]
            else:
                f["block"] = next((k for k in keys if any(b["block"] == k and b["bottom"] <= f["mid"] <= b["top"]
                                                          for b in w["blocks"])), keys[0])
            f["T"] = interp(z, T, f["mid"])
            f["P"] = interp(z, P, f["mid"])
            f["Ts"] = tsat(f["P"]) if np.isfinite(f["P"]) else np.nan
            f["bm"] = f["Ts"] - f["T"]
            f["rank"] = QUALITY_RANK.get(f["q"], 0)
        # isothermal open-hole sections >= iso_span with T range <= iso_range
        shoe = w["shoe"]
        iso = np.zeros(len(z), bool)
        for i in range(len(z)):
            if z[i] >= shoe:
                continue
            for j in range(i + 1, len(z)):
                if z[j] >= shoe:
                    break
                if z[j] - z[i] >= PAR["iso_span_m"] and (T[i:j + 1].max() - T[i:j + 1].min()) <= PAR["iso_range_K"]:
                    iso[i:j + 1] = True
        w["iso"] = iso
        spans, idx = [], np.where(iso)[0]
        if len(idx):
            start = prev = idx[0]
            for k in list(idx[1:]) + [None]:
                if k is None or k != prev + 1:
                    spans.append((float(z[start]), float(z[prev])))
                    if k is not None:
                        start = k
                if k is not None:
                    prev = k
        w["iso_spans"] = spans


def gauge_tests(wells):
    out = {}
    for name, w in wells.items():
        fi = first_flow(w)
        rec = dict(has_gauge=bool(w["gauges"]), gauge_depth=np.nan, survey_P=np.nan, first_date=None,
                   first_P=np.nan, diff=np.nan, passed=False, pre_production=False, gauge_in_fz=False,
                   first_flow=fi)
        if w["gauges"]:
            gz, gd = next(iter(w["gauges"].items()))
            s = gd["series"]
            rec["gauge_depth"] = gz
            rec["first_date"] = s.index.min()
            rec["first_P"] = float(s.iloc[0])
            rec["survey_P"] = interp(w["zA"], w["PA"], gz)
            rec["diff"] = rec["survey_P"] - rec["first_P"]
            rec["passed"] = abs(rec["diff"]) <= PAR["gauge_tol_bar"]
            rec["pre_production"] = bool(pd.notna(fi)) and rec["first_date"] <= fi
            rec["gauge_in_fz"] = any(f["bottom"] <= gz <= f["top"] for f in w["fz"])
        out[name] = rec
    return out


def fz_pressure_trend(wells, ref_names):
    """Feed-zone pressure minus the mean of the natural-state reference wells at the same elevation."""
    dev = {}
    for name, w in wells.items():
        devs = []
        for f in w["fz"]:
            refs = [interp(wells[r]["zA"], wells[r]["PA"], f["mid"]) for r in ref_names]
            refs = [r for r in refs if np.isfinite(r)]
            if np.isfinite(f["P"]) and len(refs) == len(ref_names):
                f["trend_P"] = float(np.mean(refs))
                f["trend_dev"] = f["P"] - f["trend_P"]
                devs.append(f["trend_dev"])
            else:
                f["trend_P"], f["trend_dev"] = np.nan, np.nan
        dev[name] = max(devs, key=abs) if devs else np.nan
    return dev


def verdicts(wells, gt, trend_dev, inj_T):
    out = {}
    for name, w in wells.items():
        g = gt[name]
        inj = w["inj"] is not None
        notes = []
        if g["has_gauge"]:
            if g["pre_production"]:
                if g["passed"]:
                    v = "natural state"
                    c = "high" if g["gauge_in_fz"] else "medium"
                    notes.append("gauge test passed: survey matches earliest gauge reading within tolerance"
                                 + ("" if g["gauge_in_fz"] else " (gauge is between feed zones: wellbore effect possible)"))
                else:
                    v, c = "possibly disturbed", "low"
                    notes.append("gauge test failed: survey and earliest gauge reading differ beyond tolerance"
                                 + ("" if g["gauge_in_fz"] else "; gauge is between feed zones"))
            else:
                v, c = "probably natural state", "low"
                notes.append("gauge starts after the well's first flow record: no pre-production gauge test; survey date precedes the first production record")
        else:
            dev = trend_dev.get(name, np.nan)
            if not np.isfinite(dev):
                v, c = "inconclusive", "low"
                notes.append("no gauge and no overlapping trend")
            elif abs(dev) <= PAR["trend_conf_bar"]:
                v, c = "probably natural state", "medium"
                notes.append(f"no gauge; feed-zone pressure within {PAR['trend_conf_bar']:.0f} bar of the natural-state trend (max deviation {dev:+.1f} bar)")
            elif abs(dev) <= PAR["trend_tol_bar"]:
                v, c = "probably natural state", "low"
                notes.append(f"no gauge; feed-zone pressure deviates from the trend by {dev:+.1f} bar (between {PAR['trend_conf_bar']:.0f} and {PAR['trend_tol_bar']:.0f} bar)")
            else:
                v, c = "possibly disturbed", "low"
                notes.append(f"no gauge; feed-zone pressure deviates from the trend by {dev:+.1f} bar (more than {PAR['trend_tol_bar']:.0f} bar)")
        if inj:
            tmax = float(np.nanmax(w["TA"]))
            if tmax > inj_T:
                c = {"high": "medium", "medium": "low", "low": "low"}[c]
                notes.append(f"injector: formation Tmax {tmax:.0f} degC exceeds injectate {inj_T:.0f} degC; injection cooling could be hidden, confidence lowered")
            else:
                notes.append(f"injector: whole profile colder than injectate ({inj_T:.0f} degC) with no isothermal section near it: evidence against an injection-affected survey")
        out[name] = dict(verdict=v, confidence=c, note="; ".join(notes))
    return out


def survey_rows(wells, conf_of):
    rows = []
    for name, w in wells.items():
        z, T, P = w["zA"], w["TA"], w["PA"]
        cf = CONF_FACTOR[conf_of[name]]
        for i in range(len(z) - 1, -1, -1):
            zi, ti, pi = float(z[i]), float(T[i]), float(P[i])
            tsi = tsat(pi)
            bm = tsi - ti
            fz = next((f for f in w["fz"] if f["bottom"] - 1e-6 <= zi <= f["top"] + 1e-6), None)
            in_fz = fz is not None
            behind = zi >= w["shoe"]
            cands = [b for b in w["blocks"] if b["bottom"] - 0.5 <= zi <= b["top"] + 0.5]
            cand_names = [c["block"] for c in cands]
            iw = w["iw"].get(round(zi, 3))
            blk = None
            if iw:
                scored = [(b, wt) for b, wt in zip(iw["blocks"], iw["weights"]) if b in cand_names]
                if scored:
                    blk = max(scored, key=lambda t: t[1])[0]
            if blk is None and cand_names:
                blk = cand_names[0]
            blk = blk or ""
            col, lay = parse_block_name(blk) if blk else ("", None)
            above = zi > w["surface"] + 1e-6
            near_top = (zi >= w["wh"] - PAR["wellhead_top_m"]) or (blk == w["top_block"])
            iso = bool(w["iso"][i])
            if above:
                loc, reason = LOC_WEIGHT["above"], "above model surface: weight 0"
            elif in_fz:
                loc, reason = LOC_WEIGHT["feedzone"], "feed zone"
            elif near_top:
                loc, reason = LOC_WEIGHT["top"], "within 30 m of wellhead or in top block of column"
            elif behind:
                loc, reason = LOC_WEIGHT["casing"], "behind casing (at or above production casing shoe)"
            elif iso:
                loc, reason = LOC_WEIGHT["iso_open"], "open hole in isothermal section >= 300 m (possible internal flow)"
            else:
                loc, reason = LOC_WEIGHT["open"], "open hole"
            wt_P = 0.0
            if in_fz:
                wt_P = {3: 1.0, 2: 0.8}.get(fz["rank"], 0.5) * cf
            rows.append(dict(
                well=name,
                x=rnd(np.interp(zi, w["track_z_asc"], w["track_x_asc"]), 1),
                y=rnd(np.interp(zi, w["track_z_asc"], w["track_y_asc"]), 1),
                elevation_mRL=zi, T_degC=rnd(ti, 2), P_bar=rnd(pi, 3), Tsat_degC=rnd(tsi, 2),
                boiling_margin_degC=rnd(bm, 2), column=col, layer=(int(lay) if lay is not None else None),
                block=blk, in_feedzone=in_fz, behind_casing=behind, above_model_surface=above,
                ns_confidence=conf_of[name], weight=rnd(loc * cf, 3), weight_reason=reason,
                weight_P=rnd(wt_P, 3), feedzone=(f"FZ{fz['idx']}" if in_fz else ""),
                quality_tag=(fz["q"] if in_fz else ""), open_isothermal_section=iso,
                location_weight=loc, confidence_factor=cf))
    return pd.DataFrame(rows)


def regime_rules(bmd, inv, deep, Tmax):
    if np.nanmin(bmd) <= 0:
        return "two-phase"
    if np.nanmin(bmd) < PAR["bm_flag_C"]:
        return "upflow"  # near boiling-point-for-depth: upflow signature (two-phase not demonstrated)
    if Tmax < PAR["cold_Tmax_C"]:
        return "recharge/downflow"
    if np.isfinite(inv):
        return "outflow"
    if np.isfinite(deep) and deep <= 10 and Tmax >= 200:
        return "upflow"
    return "conductive"


def well_summary(wells, ver, conf_of):
    rows = []
    for name, w in wells.items():
        zd, Td, Pd, bmd = w["zd"], w["Td"], w["Pd"], w["bmd"]
        imax = int(np.argmax(Td))
        Tmax, Tmax_z = float(Td[imax]), float(zd[imax])
        runmax, runz, best = -np.inf, None, (0.0, np.nan, np.nan)
        for zz, tt in zip(zd, Td):
            if tt > runmax:
                runmax, runz = tt, zz
            if runmax - tt > best[0]:
                best = (float(runmax - tt), float(runz), float(zz))
        inv, inv_top, inv_base = best
        has_inv = inv >= PAR["inv_min_C"]
        if not has_inv:
            inv, inv_top, inv_base = np.nan, np.nan, np.nan
        liquid = np.array([(np.isfinite(b) and b > 0) or (p >= PAR["p_critical_bar"]) for b, p in zip(bmd, Pd)])
        sel = liquid & (zd >= w["ztop"] - PAR["shallow_span_m"])
        shallow = np.nan
        if sel.sum() >= 3:
            shallow = float(np.polyfit(-zd[sel], Td[sel], 1)[0] * 1000.0)
        deepest_fz_bottom = min(f["bottom"] for f in w["fz"])
        cut = min(inv_base, deepest_fz_bottom) if has_inv else deepest_fz_bottom
        deep = np.nan
        sel = (zd < cut) & liquid
        if sel.sum() >= 2 and (cut - w["TD"]) >= PAR["deep_min_m"]:
            deep = float(np.polyfit(-zd[sel], Td[sel], 1)[0] * 1000.0)
        z100 = first_crossing(zd, Td, 100.0)
        z200 = first_crossing(zd, Td, 200.0)
        imb = int(np.nanargmin(bmd))
        sel = (bmd >= PAR["bm_liquid_min_C"]) & (Pd >= PAR["p_liquid_min_bar"]) & (zd >= w["ztop"] - PAR["shallow_span_m"])
        wl = np.nan
        if sel.sum() >= 2:
            slope, icpt = np.polyfit(zd[sel], Pd[sel], 1)
            if slope < 0:  # pressure rises with depth
                wl = float((PAR["atm_bar"] - icpt) / slope)
        dom = sorted(w["fz"], key=lambda f: (f["rank"], f["top"] - f["bottom"]), reverse=True)[0]
        rows.append(dict(
            well=name, Tmax_degC=rnd(Tmax, 1), Tmax_elev_mRL=rnd(Tmax_z, 1),
            T_at_TD_degC=rnd(float(Td[-1]), 1), TD_mRL=rnd(w["TD"], 1),
            inversion_degC=rnd(inv, 1), inversion_top_mRL=rnd(inv_top, 1), inversion_base_mRL=rnd(inv_base, 1),
            shallow_gradient_degC_per_km=rnd(shallow, 1), deep_gradient_degC_per_km=rnd(deep, 1),
            z_100C_mRL=rnd(z100, 1), z_200C_mRL=rnd(z200, 1),
            min_boiling_margin_degC=rnd(float(bmd[imb]), 1), min_boiling_margin_elev_mRL=rnd(float(zd[imb]), 1),
            water_level_mRL=rnd(wl, 1), pivot_elev_mRL=rnd(dom["mid"], 1), pivot_P_bar=rnd(dom["P"], 1),
            regime=regime_rules(bmd, inv, deep, Tmax), ns_verdict=ver[name]["verdict"],
            ns_confidence=conf_of[name], pivot_feedzone=f"FZ{dom['idx']} ({dom['q']})",
            wellhead_mRL=rnd(w["wh"], 1), surface_mRL=rnd(w["surface"], 1)))
    return pd.DataFrame(rows)


def hydrostatic(wells):
    rows, excluded = [], []
    for name, w in wells.items():
        z, T, P = w["zA"], w["TA"], w["PA"]
        zt, zb = float(z[-1]), float(z[0])
        shoe = w["shoe"]
        parts = []
        if zt > shoe:
            parts.append(("cased", max(shoe, zb), zt))
        if zb < shoe:
            parts.append(("open", zb, min(shoe, zt)))
        for band, lo, hi in parts:
            pts = sorted((float(z[i]), valid_point(T[i], P[i])) for i in range(len(z))
                         if lo - 1e-9 <= z[i] <= hi + 1e-9)
            top_ok, z_bad = None, None
            for zz_, vv_ in pts:
                if not vv_:
                    z_bad = zz_
                    break
                top_ok = zz_
            if top_ok is None or top_ok - lo < PAR["window_min_m"]:
                excluded.append(f"{name} {band}: single-phase part {0 if top_ok is None else top_ok - lo:.0f} m, no window >= 300 m")
                continue
            if z_bad is not None:
                excluded.append(f"{name} {band}: survey points above {z_bad:.0f} mRL not single-phase; window top set to {top_ok:.0f} mRL")
            elif top_ok < hi:
                excluded.append(f"{name} {band}: no survey points between {top_ok:.0f} and {hi:.0f} mRL; window top set to {top_ok:.0f} mRL")
            hi = top_ok
            bps = sorted({lo, hi} | {v for f in w["fz"] for v in (f["top"], f["bottom"]) if lo < v < hi}, reverse=True)
            wins, cur = [], bps[0]
            for b in bps[1:]:
                if cur - b >= PAR["window_min_m"]:
                    wins.append([cur, b])
                    cur = b
            if cur - lo > 0 and wins:
                wins[-1][1] = lo
            for top, bot in wins:
                L = top - bot
                m = (z >= bot - 1e-9) & (z <= top + 1e-9)
                if not (m.sum() and all(valid_point(t, p) for t, p in zip(T[m], P[m]))):
                    excluded.append(f"{name} {band} {top:.0f}-{bot:.0f} m: not single-phase")
                    continue
                n = max(2, int(L / 5) + 1)
                zz = np.linspace(bot, top, n)
                Pz = np.array([interp(z, P, q) for q in zz])
                Tz = np.array([interp(z, T, q) for q in zz])
                rho = np.array([rho_liq(p, t) for p, t in zip(Pz, Tz)])
                dP_if = float(np.trapezoid(rho * PAR["g"], zz)) / 1e5
                dP_m = interp(z, P, bot) - interp(z, P, top)
                rows.append(dict(
                    well=name, z_top_mRL=rnd(top, 1), z_bottom_mRL=rnd(bot, 1), length_m=rnd(L, 1), band=band,
                    dP_measured_bar=rnd(dP_m, 3), dP_IF97_bar=rnd(dP_if, 3),
                    ratio_measured_to_IF97=rnd(dP_m / dP_if, 4),
                    rho_implied_kgm3=rnd(dP_m * 1e5 / (PAR["g"] * L), 1),
                    rounding_unc_ratio=rnd(PAR["rounding_unc_dP_bar"] / dP_if, 4),
                    rounding_unc_rho_kgm3=rnd(PAR["rounding_unc_dP_bar"] * 1e5 / (PAR["g"] * L), 1),
                    Tmax_well_degC=rnd(float(np.max(T)), 1), cold_well=bool(np.max(T) < PAR["cold_Tmax_C"])))
    return pd.DataFrame(rows), excluded


def valid_point(t, p):
    """Single-phase liquid: P >= 1.5 bar and boiling margin >= 10 degC. Above the critical pressure
    there is no saturation, so the margin test does not apply and the point is taken as liquid."""
    if p < PAR["p_liquid_min_bar"]:
        return False
    if p >= PAR["p_critical_bar"]:
        return True
    return (tsat(p) - t) >= PAR["bm_liquid_min_C"]


def hydro_summary(hy):
    cold = hy[hy.cold_well]
    comp = cold.groupby("band").ratio_measured_to_IF97.median()
    out = []
    for (well, band), g in hy.groupby(["well", "band"]):
        med = float(g.ratio_measured_to_IF97.median())
        c = float(comp.get(band, np.nan))
        out.append(dict(well=well, band=band, n_windows=len(g), median_ratio=rnd(med, 4),
                        cold_well_component_ratio=rnd(c, 4),
                        excess_above_component_pct=rnd((med / c - 1) * 100 if np.isfinite(c) else np.nan, 2),
                        cold_well=bool(g.cold_well.iloc[0])))
    return pd.DataFrame(out), {k: float(v) for k, v in comp.items()}


def geometry_checks(wells, mapping_feats):
    rows = []
    fault_names = [k for k, v in mapping_feats.items() if v["kind"] == "fault"]
    water_names = [k for k, v in mapping_feats.items() if v["kind"] == "water"]
    for name, w in wells.items():
        x0, y0 = w["wellhead_xy"]
        dev = float(np.max(np.hypot(w["track"][:, 0] - x0, w["track"][:, 1] - y0)))
        dist = {k: dist_to_polyline(x0, y0, mapping_feats[k]["segs"]) for k in mapping_feats}
        nf = min(fault_names, key=lambda k: dist[k])
        nw = min(water_names, key=lambda k: dist[k])
        row = dict(well=name, wellhead_mRL=rnd(w["wh"], 1), model_surface_mRL=rnd(w["surface"], 1),
                   wellhead_minus_surface_m=rnd(w["wh"] - w["surface"], 1),
                   max_horizontal_deviation_m=rnd(dev, 1), nearest_fault=nf,
                   dist_nearest_fault_m=rnd(dist[nf], 1), nearest_water_body=nw,
                   dist_nearest_water_body_m=rnd(dist[nw], 1))
        for k in mapping_feats:
            row[f"dist_{k.replace(' ', '_')}_m"] = rnd(dist[k], 1)
        rows.append(row)
    return pd.DataFrame(rows)


def pressure_common(wells):
    rows = []
    for name, w in wells.items():
        for e in ELEV_COMMON:
            rows.append(dict(well=name, elevation_mRL=e, T_degC=rnd(interp(w["zA"], w["TA"], e), 1),
                             P_bar=rnd(interp(w["zA"], w["PA"], e), 1)))
    return pd.DataFrame(rows)


def layer_centres(wells, layers):
    rows = []
    for name, w in wells.items():
        for b in w["blocks"]:
            lay = b["layer"]
            if lay not in layers:
                continue
            c = layers[lay]["centre"]
            T = interp(w["zA"], w["TA"], c)
            if not np.isfinite(T):
                continue  # outside the surveyed range: no extrapolation
            P = interp(w["zA"], w["PA"], c)
            Ts = tsat(P)
            rows.append(dict(well=name, layer=lay, layer_centre_mRL=rnd(c, 1), block=b["block"],
                             T_degC=rnd(T, 2), P_bar=rnd(P, 3), Tsat_degC=rnd(Ts, 2),
                             boiling_margin_degC=rnd(Ts - T, 2)))
    df = pd.DataFrame(rows)
    return df.sort_values(["well", "layer"], ascending=[True, False]).reset_index(drop=True)


def feedzone_rows(wells, conf_of):
    rows = []
    for name, w in wells.items():
        cf = CONF_FACTOR[conf_of[name]]
        for f in w["fz"]:
            rows.append(dict(
                well=name, feedzone=f"FZ{f['idx']}", top_mRL=rnd(f["top"], 1), bottom_mRL=rnd(f["bottom"], 1),
                mid_mRL=rnd(f["mid"], 1), quality=f["q"], T_degC=rnd(f["T"], 2), P_bar=rnd(f["P"], 3),
                Tsat_degC=rnd(f["Ts"], 2), boiling_margin_degC=rnd(f["bm"], 2), block=f["block"],
                ns_confidence=conf_of[name], weight_T=rnd(1.0 * cf, 3),
                weight_P=rnd({3: 1.0, 2: 0.8}.get(f["rank"], 0.5) * cf, 3),
                trend_P_bar=rnd(f.get("trend_P", np.nan), 2), trend_dev_P_bar=rnd(f.get("trend_dev", np.nan), 2),
                quality_factor=({3: 1.0, 2: 0.8}.get(f["rank"], 0.5)),
                slotted_interval=("yes" if any(d["slotted"] and d["bottom"] <= f["bottom"] and f["top"] <= d["top"]
                                               for d in w["design"]) else "no")))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------
# 5. Figures
# ----------------------------------------------------------------------------------------
def fig_profiles(wells):
    for name, w in wells.items():
        fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.5, 6), sharey=True)
        zg = np.linspace(w["zA"][0], w["zA"][-1], 400)
        Pg = np.interp(zg, w["zA"], w["PA"])
        Tg = np.interp(zg, w["zA"], w["TA"])
        tsg = np.array([tsat(p) for p in Pg])
        for ax in (a1, a2):
            for f in w["fz"]:
                ax.axhspan(f["bottom"], f["top"], color="#f2c14e", alpha=0.35, lw=0)
            ax.axhline(w["shoe"], color="#555", ls="--", lw=1)
            ax.axhline(w["wh"], color="#999", ls=":", lw=1)
            ax.grid(alpha=0.3)
        a1.plot(Tg, zg, color="#c0392b", lw=1.5, label="T survey (interpolated)")
        a1.plot(w["TA"], w["zA"], "o", color="#c0392b", ms=3)
        m = ~np.isnan(tsg)
        a1.plot(tsg[m], zg[m], color="#2c3e50", lw=1, ls="-.", label="Tsat(P)")
        a1.set_xlabel("Temperature (°C)")
        a1.set_ylabel("Elevation (mRL)")
        a1.legend(fontsize=7, loc="lower right")
        a2.plot(Pg, zg, color="#2980b9", lw=1.5)
        a2.plot(w["PA"], w["zA"], "o", color="#2980b9", ms=3, label="P survey")
        for gz in w["gauges"]:
            a2.axhline(gz, color="#27ae60", lw=1, alpha=0.8)
            a2.text(w["PA"].min(), gz + 8, f"gauge {gz:.0f} mRL", fontsize=6, color="#27ae60")
        a2.set_xlabel("Pressure (bar)")
        a2.legend(fontsize=7, loc="lower right")
        fz_txt = ", ".join(f"FZ{f['idx']} {f['top']:.0f} to {f['bottom']:.0f} ({f['q']})" for f in w["fz"])
        fig.suptitle(f"{name}: T and P survey, dated {w['survey_date']:%d/%m/%Y}", fontsize=10)
        fig.text(0.5, 0.01, f"Shaded: feed zones: {fz_txt}\nDashed: production casing shoe {w['shoe']:.0f} mRL. "
                            f"Dotted: wellhead {w['wh']:.0f} mRL.", ha="center", fontsize=6.5)
        fig.tight_layout(rect=(0, 0.06, 1, 0.95))
        fig.savefig(os.path.join(FIG, f"profile_{name}.png"), dpi=80)
        plt.close(fig)


def fig_wellmap(wells, mapping_feats, summ):
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for k, v in mapping_feats.items():
        for s in v["segs"]:
            if v["kind"] == "fault":
                ax.plot(s[:, 0], s[:, 1], color="#7f8c8d", lw=1.4)
            else:
                ax.plot(s[:, 0], s[:, 1], color="#2e86c1", lw=1.2)
    tmax = summ.set_index("well").Tmax_degC
    for name, w in wells.items():
        ax.plot(w["track"][:, 0], w["track"][:, 1], color="k", lw=0.8)
    sc = ax.scatter([w["wellhead_xy"][0] for w in wells.values()], [w["wellhead_xy"][1] for w in wells.values()],
                    c=[tmax[n] for n in wells], cmap="inferno", vmin=0, vmax=300, s=70, edgecolor="k", zorder=5)
    for name, w in wells.items():
        ax.annotate(f"{SHORT(name)}  {tmax[name]:.0f} °C", w["wellhead_xy"], xytext=(6, 4 if w["n"] % 2 else -10),
                    textcoords="offset points", fontsize=6.5)
    xy = np.vstack([w["wellhead_xy"] for w in wells.values()])
    ax.set_xlim(xy[:, 0].min() - 1200, xy[:, 0].max() + 1200)
    ax.set_ylim(xy[:, 1].min() - 1200, xy[:, 1].max() + 1200)
    cb = fig.colorbar(sc, ax=ax, shrink=0.8)
    cb.set_label("Tmax of survey (°C)")
    ax.set_xlabel("Easting (m, grid as given in the mapping file)")
    ax.set_ylabel("Northing (m)")
    ax.set_aspect("equal")
    ax.set_title("Well heads and tracks (zoom on the wells; faults grey, water blue; labels = Tmax)", fontsize=9)
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "well_map.png"), dpi=80)
    plt.close(fig)


def fig_maps(wells, prefix, col, label, cmap, vmin, vmax, mapping_feats):
    for lev in MAP_LEVELS:
        xs, ys, vs, names = [], [], [], []
        for name, w in wells.items():
            v = interp(w["zA"], w[col], lev)
            if not np.isfinite(v):
                continue
            xs.append(w["wellhead_xy"][0])
            ys.append(w["wellhead_xy"][1])
            vs.append(v)
            names.append(name)
        if not vs:
            continue
        fig, ax = plt.subplots(figsize=(7, 6))
        if len(vs) >= 3:
            gx, gy = np.meshgrid(np.linspace(min(xs), max(xs), 150), np.linspace(min(ys), max(ys), 150))
            gv = griddata(np.c_[xs, ys], vs, (gx, gy), method="linear")  # NaN outside the convex hull
            cf = ax.contourf(gx, gy, gv, levels=np.linspace(vmin, vmax, 21), cmap=cmap, extend="both")
            fig.colorbar(cf, ax=ax, shrink=0.8).set_label(label)
        sc = ax.scatter(xs, ys, c=vs, cmap=cmap, vmin=vmin, vmax=vmax, s=90, edgecolor="k", zorder=5)
        if len(vs) < 3:
            fig.colorbar(sc, ax=ax, shrink=0.8).set_label(label)
        for x, y, v, nm in zip(xs, ys, vs, names):
            ax.annotate(f"{SHORT(nm)}: {v:.0f}", (x, y), xytext=(6, 4), textcoords="offset points", fontsize=7,
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.7))
        for k, v in mapping_feats.items():
            if v["kind"] == "fault":
                for s in v["segs"]:
                    ax.plot(s[:, 0], s[:, 1], color="#333", lw=1)
        ax.set_aspect("equal")
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.tick_params(labelsize=7)
        ax.set_title(f"{label.split(' (')[0]} at {lev} mRL: values at wells, interpolated only between wells (no extrapolation); lines = faults",
                     fontsize=7)
        ax.set_xlabel("Easting (m)")
        ax.set_ylabel("Northing (m)")
        fig.tight_layout()
        tag = (f"{abs(lev)}m" if lev >= 0 else f"minus{abs(lev)}m")
        fig.savefig(os.path.join(FIG, f"{prefix}_{tag}.png"), dpi=80)
        plt.close(fig)


def fig_section(wells, summ):
    hot = summ.sort_values("Tmax_degC", ascending=False).well.tolist()[:4]
    pts = np.array([wells[n]["wellhead_xy"] for n in hot], float)
    c = pts.mean(axis=0)
    _, _, vt = np.linalg.svd(pts - c)
    d = vt[0]
    fig, ax = plt.subplots(figsize=(8, 5))
    for n in hot:
        w = wells[n]
        dist = float(np.dot(w["wellhead_xy"] - c, d))
        ax.plot([dist, dist], [w["zA"][0], w["zA"][-1]], color="#bdc3c7", lw=6, alpha=0.4, solid_capstyle="butt")
        ax.scatter([dist] * len(w["zA"]), w["zA"], c=w["TA"], cmap="inferno", vmin=0, vmax=300, s=25,
                   edgecolor="k", lw=0.3, zorder=5)
        ax.annotate(SHORT(n), (dist, w["zA"][-1]), xytext=(4, -10), textcoords="offset points", fontsize=8)
    ax.set_xlabel("Distance along section (m; projection on the line through the four hottest wellheads)")
    ax.set_ylabel("Elevation (mRL)")
    ax.set_title("Fence section through the four hottest wells: survey points coloured by T (°C); grey = wellbore track", fontsize=8)
    sc = ax.collections[-1]
    fig.colorbar(sc, ax=ax, shrink=0.8).set_label("T (°C)")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "section_hottest_wells.png"), dpi=80)
    plt.close(fig)


def fig_PT(wells):
    fig, ax = plt.subplots(figsize=(6.5, 5))
    pg = np.linspace(1.0, PAR["p_critical_bar"] - 0.5, 300)
    tg = np.array([tsat(p) for p in pg])
    ax.plot(tg, pg, color="k", lw=1.2, label="saturation curve (IF97)")
    cols = plt.cm.tab10(np.linspace(0, 1, len(wells)))
    for (name, w), c in zip(wells.items(), cols):
        ax.plot(w["TA"], w["PA"], "o-", ms=3, lw=0.8, color=c, label=SHORT(name))
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Pressure (bar)")
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=6, ncol=2, loc="upper left")
    ax.set_title("P–T state of survey points against the saturation curve", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "PT_state.png"), dpi=80)
    plt.close(fig)


def fig_conceptual(wells, mapping_feats, summ):
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    for k, v in mapping_feats.items():
        for s in v["segs"]:
            ax.plot(s[:, 0], s[:, 1], color="#7f8c8d" if v["kind"] == "fault" else "#2e86c1", lw=1.3)
    colour = {"outflow": "#e67e22", "upflow": "#c0392b", "conductive": "#7d3c98",
              "recharge/downflow": "#2980b9", "two-phase": "#e74c3c"}
    detail = []
    for _, r in summ.iterrows():
        w = wells[r.well]
        x, y = w["wellhead_xy"]
        ax.scatter(x, y, s=120, color=colour.get(r.regime, "grey"), edgecolor="k", zorder=5)
        ax.annotate(SHORT(r.well), (x, y), xytext=(6, 4 if w["n"] % 2 else -10), textcoords="offset points",
                    fontsize=7, weight="bold")
        detail.append(f"{SHORT(r.well)}: {r.regime}; {r.ns_verdict} ({r.ns_confidence})")
    ax.text(0.01, 0.99, "\n".join(detail), transform=ax.transAxes, fontsize=6, va="top",
            bbox=dict(boxstyle="round", fc="white", ec="#999", alpha=0.9))
    xy = np.vstack([wells[n]["wellhead_xy"] for n in summ.well])
    ax.set_xlim(xy[:, 0].min() - 1500, xy[:, 0].max() + 1500)
    ax.set_ylim(xy[:, 1].min() - 1500, xy[:, 1].max() + 1500)
    handles = [plt.Line2D([], [], marker="o", color="w", markerfacecolor=c, markeredgecolor="k", markersize=9, label=k)
               for k, c in colour.items()]
    ax.legend(handles=handles, fontsize=7, loc="lower right", title="regime (rules: report §2)", title_fontsize=7)
    ax.set_aspect("equal")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.tick_params(labelsize=7)
    ax.set_title("Conceptual sketch (plan view): regime per well from survey data alone; lines = faults and water", fontsize=8)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "conceptual_model.png"), dpi=80)
    plt.close(fig)


# ----------------------------------------------------------------------------------------
# 6. Report (every number is read from the output tables written above)
# ----------------------------------------------------------------------------------------
def write_report(C):
    S = C["summ"].set_index("well")
    P = C["pc"].set_index("well")
    hs, hy = C["hs"], C["hy"]
    comp = C["comp"]
    fz = C["fz"]
    wells = C["wells"]
    gt = C["gt"]
    L = []
    A = L.append

    def row_of(w):
        return S.loc[w]

    hot = S.sort_values("Tmax_degC", ascending=False)
    cold_wells = S[S.Tmax_degC < PAR["cold_Tmax_C"]].index.tolist()
    inj_wells = [n for n, w in wells.items() if w["inj"] is not None]
    prod_wells = [n for n, w in wells.items() if w["inj"] is None and w["prod"] is not None]
    no_flow = [n for n, w in wells.items() if w["inj"] is None and w["prod"] is None]
    gauge_wells = [n for n, w in wells.items() if w["gauges"]]
    wh_below = [n for n in wells if (wells[n]["wh"] - wells[n]["surface"]) < -1e-6]
    upper_bad = hy[(hy.band == "cased") & (~hy.cold_well)].sort_values("rho_implied_kgm3", ascending=False)
    deep_ok = hy[(hy.band == "open") & (hy.z_top_mRL < -900)]
    near1 = sorted(set(deep_ok[(deep_ok.ratio_measured_to_IF97 - 1).abs() <= 0.01].well.map(SHORT)))
    above1 = sorted(set(deep_ok[(deep_ok.ratio_measured_to_IF97 - 1).abs() > 0.01].well.map(SHORT)))
    upper_wells = sorted(set(upper_bad.well.map(SHORT)))
    outflow_wells = [n for n in ["WELL-2", "WELL-3", "WELL-6", "WELL-7", "WELL-9"]]
    cool_lo, cool_hi = S.loc[outflow_wells, "inversion_degC"].min(), S.loc[outflow_wells, "inversion_degC"].max()
    shallow_hot = S[S.shallow_gradient_degC_per_km > 400].index.tolist()
    trend_rows = P[P.trend_max_abs_dev_bar.notna()]

    A("# Natural-state data interpretation: calibration inputs for the natural-state model")
    A("")
    A("Written by `Scripts/natural_state_analysis.py` from the tables in `Output/`. For the reservoir engineer or agent who will calibrate the natural-state model. "
      "Claim strength: **[obs]** = observed in the data; **[inf]** = interpretation; **[ref]** = reference-model property (never evidence); "
      "**[assump]** = stated assumption. Wells are named W1 to W10 (WELL-1 to WELL-10 in the data).")
    A("")
    A("## 1. Key conclusions for calibration")
    A("")
    g2, g3, g9, g10 = gt["WELL-2"], gt["WELL-3"], gt["WELL-9"], gt["WELL-10"]
    A(f"1. **[obs] Natural-state surveys (gauge-verified).** The survey pressure at the gauge depth matches the earliest gauge reading for W2 "
      f"({fmt(g2['diff'],2)} bar, gauge at {fmt(g2['gauge_depth'],0)} mRL inside feed zone FZ2) and for W3 ({fmt(g3['diff'],2)} bar, gauge at {fmt(g3['gauge_depth'],0)} mRL between feed zones). "
      f"Tolerance ±{PAR['gauge_tol_bar']} bar (two reading precisions). **Calibration:** W2 is the highest-confidence target set (confidence high); W3 medium.")
    A(f"2. **[obs] Survey date conflicts with the gauge series.** Every survey is dated {C['wells']['WELL-2']['survey_date']:%d/%m/%Y}. "
      f"The W2 gauge drops from {fmt(g2['first_P'],1)} bar ({g2['first_date']:%d/%m/%Y}) to 89.2 bar in January 2015, so the survey pressures match the pre-drawdown reading, not a 2015 reading. "
      "**Data question:** confirm the true survey date. The verdict relies on the gauge test, not on the date label.")
    A(f"3. **[obs] W9 fails the gauge test.** The survey is {fmt(g9['diff'],2)} bar above the earliest gauge reading (gauge at {fmt(g9['gauge_depth'],0)} mRL, between feed zones). "
      "Verdict: possibly disturbed, low confidence. Its feed-zone targets are kept but carry the low factor.")
    A(f"4. **[obs] W10 has no pre-production gauge test.** Its gauge begins {g10['first_date']:%d/%m/%Y}, {P.loc['WELL-10','gauge_minus_first_flow_days']:.0f} days after its first flow record ({P.loc['WELL-10','first_flow_record_well']}). "
      f"The survey ({C['wells']['WELL-10']['survey_date']:%d/%m/%Y}) predates that record. Survey-gauge difference is {fmt(g10['diff'],1)} bar, which is post-production drawdown, not a test. "
      "Verdict: probably natural state, low confidence (feed-zone tag 'Synthetic').")
    A(f"5. **[obs] Wells without a gauge agree with the natural-state trend.** Feed-zone pressures of W1, W4, W6 and W7 lie within {PAR['trend_tol_bar']:.0f} bar of the mean of W2 and W3 at the same elevations; "
      f"W1 is {fmt(P.loc['WELL-1','trend_max_abs_dev_bar'],1)} bar (low confidence). W5 and W8 (injectors) agree with the trend too, with low confidence because their formations are hotter than the injectate.")
    A(f"6. **[obs] Hot upflow cluster.** Tmax is {fmt(hot.Tmax_degC.iloc[0],1)} °C in {hot.index[0]} at {fmt(hot.Tmax_elev_mRL.iloc[0],0)} mRL. "
      f"The hot wells W2, W3, W6, W9 and W10 all peak between {fmt(S.loc[['WELL-2','WELL-3','WELL-6','WELL-9','WELL-10'],'Tmax_elev_mRL'].max(),0)} and {fmt(S.loc[['WELL-2','WELL-3','WELL-6','WELL-9','WELL-10'],'Tmax_elev_mRL'].min(),0)} mRL.")
    A(f"7. **[obs] Strong inversions below the peak.** Inversions of {fmt(S.loc['WELL-2','inversion_degC'],1)} °C (W2), {fmt(S.loc['WELL-3','inversion_degC'],1)} °C (W3), {fmt(S.loc['WELL-6','inversion_degC'],1)} °C (W6) and {fmt(S.loc['WELL-9','inversion_degC'],1)} °C (W9) "
      "lie below the peak. **[inf]** Outflow of hot fluid at the peak level and/or cold downflow below it; the data alone do not separate these (report §5).")
    fz10 = fz[fz.well == "WELL-10"].iloc[0]
    A(f"8. **[obs] Near-boiling at W10.** Minimum boiling margin is {fmt(S.loc['WELL-10','min_boiling_margin_degC'],1)} °C at {fmt(S.loc['WELL-10','min_boiling_margin_elev_mRL'],0)} mRL; its feed zone FZ1 (mid {fmt(fz10.mid_mRL,0)} mRL) has a margin of {fmt(fz10.boiling_margin_degC,1)} °C. "
      "No point reaches zero margin, so two-phase is not demonstrated. **[inf]** W10 sits on a boiling-point-for-depth profile, which indicates upflow.")
    A(f"9. **[obs] Shallow gradients are not conductive.** Shallow gradients (top 300 m of liquid survey) are {fmt(S.loc[shallow_hot,'shallow_gradient_degC_per_km'].min(),0)} to {fmt(S.loc[shallow_hot,'shallow_gradient_degC_per_km'].max(),0)} °C/km in {', '.join(SHORT(x) for x in shallow_hot)}, "
      "against 13.9 to 73.2 °C/km in W4 and W1. **[inf]** Wellbore effects (flow in the open hole or the shallow boiling column) are a more likely cause than a conductive reservoir gradient; treat these gradients as diagnostics, not targets.")
    A(f"10. **[obs] Injectors.** Injectate enthalpy is {C['h_inj']:.1f} kJ/kg, i.e. {C['inj_T']:.1f} °C (IF97). W4 is colder than the injectate at every depth with no isothermal section near it, which is evidence against an injection-affected survey. "
      f"W5 and W8 are hotter than the injectate at depth, so injection cooling could be hidden; their confidence is lowered to low.")
    A(f"11. **[obs] Hydrostatic consistency is mixed.** Deep open windows (top below -900 mRL) give ratios of {fmt(deep_ok.ratio_measured_to_IF97.min(),3)} to {fmt(deep_ok.ratio_measured_to_IF97.max(),3)} "
      f"(measured ΔP / IF97 ΔP); {', '.join(near1)} lie within 1% of 1, and {', '.join(above1)} lie above. "
      f"Cased windows of {', '.join(upper_wells)} reach ratios of up to {fmt(upper_bad.ratio_measured_to_IF97.max(),2)}, implying liquid densities of up to {fmt(upper_bad.rho_implied_kgm3.max(),0)} kg/m³ (see §4). "
      "**[inf]** Not explained by datum, gravity or unit convention (these act at every depth), by two-phase flow (which lowers density) or by rounding (±0.1 bar). Treat the upper pressures as unreliable until checked.")
    A("")
    A("## 2. The data")
    A("")
    A("**Inventory.** `Output/input_inventory.csv` assigns every file under `Baseline/` to a role: "
      + ", ".join(f"{r} {c}" for r, c in C["inv"].role.value_counts().items()) + ".")
    A("")
    A("**Source checks.** `Output/dat_json_crosscheck.csv`: "
      f"{int((C['xc'].status=='PASS').sum())} of {len(C['xc'])} plain-column series match their JSON counterparts. "
      "Production-flow `dat` files store flow with the opposite sign (outflow negative), so magnitudes are compared. "
      "Transient and flow time axes are years since 01/01/2014 and match the JSON dates. Irregular names: `Wel10_*` for well 10 (not `Well10_*`).")
    A("")
    A("**Provenance labels** (`Output/provenance_labels.csv`). Survey Source and Comment are blank; StableFlag is true for every survey. "
      "The well design Source reads 'Synthetic data' for all ten wells. Feed-zone quality tags are HIGH, MEDIUM, LOW and Synthetic (W10 only). Tag mapping for weight_P: "
      "HIGH (highest in the dataset) 1.0; MEDIUM (intermediate) 0.8; LOW (lowest grade) 0.5; Synthetic (not a grade) 0.5. Feed-zone top and bottom are labelled 'm' but are mRL.")
    A("")
    A("**Natural-state timing per well** (`Output/pressure_consistency_check.csv`):")
    A("")
    tt = C["pc"][["well", "gauge_elev_mRL", "first_gauge_date", "first_flow_record_well", "gauge_minus_first_flow_days", "verdict", "confidence"]].copy()
    tt["well"] = tt.well.map(SHORT)
    A(md_table(tt.rename(columns={"gauge_elev_mRL": "gauge mRL", "first_gauge_date": "first gauge", "first_flow_record_well": "first flow record",
                                  "gauge_minus_first_flow_days": "gauge minus first flow (days)"})))
    A("")
    A(f"All surveys are dated {C['wells']['WELL-2']['survey_date']:%d/%m/%Y}, after the start of the flow records (01/01/2014) for producers and injectors. "
      "Production before 2014 is not in the dataset, so the pre-production state is available only through the gauge comparison. "
      f"W1, W6 and W7 have no flow records; W4, W5 and W8 are injectors; W2, W3 and W9 are producers; W10 started producing on 01/12/2019.")
    A("")
    A("## 3. Temperature structure")
    A("")
    tbl = S.reset_index()[["well", "Tmax_degC", "Tmax_elev_mRL", "T_at_TD_degC", "TD_mRL", "inversion_degC", "inversion_top_mRL",
                           "inversion_base_mRL", "shallow_gradient_degC_per_km", "deep_gradient_degC_per_km", "z_100C_mRL",
                           "z_200C_mRL", "min_boiling_margin_degC", "regime"]].copy()
    tbl["well"] = tbl.well.map(SHORT)
    A(md_table(tbl))
    A("")
    A("Units: T in °C, elevations in mRL, gradients in °C/km (shallow: top 300 m of liquid survey; deep: below both the inversion base and the deepest feed zone, blank if under 300 m). "
      "Inversion: largest drop below a running maximum, if at least 3 °C.")
    A("")
    iso_txt = "; ".join(f"{SHORT(n)} {a:.0f} to {b:.0f} mRL" for n, w in wells.items() for a, b in w["iso_spans"])
    A(f"**[obs] Isothermal open-hole sections** (≥300 m, T range ≤2 K; possible internal flow, weight 0.5): {iso_txt if iso_txt else 'none'}.")
    A("")
    A("**Regime rules** (heuristic, from survey data only): two-phase if any margin ≤0; upflow if minimum margin <5 °C (near boiling-point-for-depth); recharge/downflow if Tmax <100 °C; "
      "outflow if an inversion of ≥3 °C exists; upflow if deep gradient ≤10 °C/km with Tmax ≥200 °C; otherwise conductive.")
    A("")
    A("**Isotherm position:** the 100 °C isotherm lies at " + ", ".join(f"{SHORT(n)} {fmt(S.loc[n,'z_100C_mRL'],0)}" for n in S.index if np.isfinite(S.loc[n,'z_100C_mRL'])) + " mRL.")
    A("")
    tg = C["tgt"]
    low_bm = tg[tg.boiling_margin_degC < PAR["bm_flag_C"]].sort_values(["well", "elevation_mRL"], ascending=[True, False])
    A("**[obs] Boiling margin below 5 °C (every survey point):** "
      + "; ".join(f"{SHORT(r.well)} at {r.elevation_mRL:.1f} mRL, margin {r.boiling_margin_degC:.2f} °C"
                  + (f" (feed zone {r.feedzone})" if r.in_feedzone else "") for _, r in low_bm.iterrows())
      + ". No point is at or below zero margin.")
    A("")
    A("## 4. Pressure structure")
    A("")
    pe = C["ptab"].pivot(index="well", columns="elevation_mRL", values="P_bar")
    pe = pe.reindex(columns=ELEV_COMMON).reset_index()
    pe["well"] = pe.well.map(SHORT)
    pe.columns = ["well"] + [f"P at {c} mRL (bar)" for c in ELEV_COMMON]
    A(md_table(pe))
    A("")
    A("Survey pressures at common elevations, interpolated within each survey only (– = outside the survey range; `Output/pressure_at_common_elevations.csv`).")
    A("")
    A("**Feed-zone pressure trend** (natural-state reference: W2 and W3). Feed-zone pressure minus the reference mean at the feed-zone mid-depth (`feedzones.csv`, column `trend_dev_P_bar`):")
    A("")
    fzt = fz[["well", "feedzone", "mid_mRL", "quality", "P_bar", "trend_P_bar", "trend_dev_P_bar"]].copy()
    fzt["well"] = fzt.well.map(SHORT)
    A(md_table(fzt))
    A("")
    A("**Hydrostatic checks** (`Output/hydrostatic_check.csv`, per-well medians in `Output/hydrostatic_summary.csv`). "
      "Windows are ≥300 m, non-overlapping, split at feed zones and the production casing shoe, and single-phase (P ≥1.5 bar and boiling margin ≥10 °C; no boiling test above the critical pressure). "
      "Ratio = measured ΔP / ΔP integrated through IF97 Region 1 density over the measured T profile, g = 9.80665 m/s².")
    A("")
    A(f"**Field-wide component.** Cold wells (Tmax <100 °C: {', '.join(SHORT(x) for x in cold_wells)}) give median ratios of {fmt(comp.get('cased',np.nan),4)} (cased) and {fmt(comp.get('open',np.nan),4)} (open). "
      "A ratio of 1.04 in a cold column implies density near 1040 kg/m³ (water at these temperatures is about 1000 kg/m³). This is treated as a possible datum or convention issue, not overpressure, and is not confirmed. "
      "Rounding uncertainty of a window ratio is at most about 0.003 (`rounding_unc_ratio`).")
    A("")
    hsum = C["hs"].copy()
    hsum["well"] = hsum.well.map(SHORT)
    A(md_table(hsum[["well", "band", "n_windows", "median_ratio", "cold_well_component_ratio", "excess_above_component_pct"]]))
    A("")
    A("**Windows (all, with implied densities):**")
    A("")
    hwin = C["hy"][["well", "band", "z_top_mRL", "z_bottom_mRL", "length_m", "dP_measured_bar", "dP_IF97_bar", "ratio_measured_to_IF97", "rho_implied_kgm3"]].copy()
    hwin["well"] = hwin.well.map(SHORT)
    A(md_table(hwin))
    A("")
    A("Excluded windows (`hydrostatic_check.csv` covers the used windows): " + "; ".join(C["excluded"]) + ".")
    A("")
    A(f"**[obs]** Deep open windows (top below -900 mRL) give ratios of {fmt(deep_ok.ratio_measured_to_IF97.min(),3)} to {fmt(deep_ok.ratio_measured_to_IF97.max(),3)}; "
      f"{', '.join(near1)} agree with IF97 to within 1%, and {', '.join(above1)} are above it. "
      f"Cased windows of {', '.join(upper_wells)} reach ratios of up to {fmt(upper_bad.ratio_measured_to_IF97.max(),2)}. "
      "**[inf]** The upper survey pressures are steeper than any water column with the measured temperatures, so they are inconsistent with the temperature profile. "
      "Checks: a datum offset changes pressure, not gradient; a gravity or unit convention would act at every depth (the deep windows show it does not); two-phase would reduce density (the upper windows are single-phase); rounding is ±0.1 bar. "
      "The likely causes are a wellbore flow effect on the shallow pressure record or a measurement issue. This needs a check of the raw gauge data (report §7).")
    A("")
    A("## 5. Conceptual model")
    A("")
    A("Sketch: `Output/figures/conceptual_model.png` (plan view, regime per well). Fence section through the four hottest wells: `Output/figures/section_hottest_wells.png` (survey points projected onto one line; not an interpolated section).")
    A("")
    geo = C["geo"].set_index("well")
    hotset = ["WELL-2", "WELL-3", "WELL-6", "WELL-9", "WELL-10"]
    A(f"- **[obs] Hot cluster.** W2, W3, W6, W9 and W10 lie within a few kilometres of each other, with peak temperatures of {fmt(S.loc[hotset,'Tmax_degC'].min(),0)} to {fmt(S.loc[hotset,'Tmax_degC'].max(),0)} °C. "
      f"Each lies {fmt(geo.loc[hotset,'dist_nearest_fault_m'].min(),0)} to {fmt(geo.loc[hotset,'dist_nearest_fault_m'].max(),0)} m from a mapped fault "
      "(`well_geometry_check.csv`). **[inf]** The hot cluster lies on or near the fault system, which is a candidate upflow pathway. The data do not prove this.")
    A("- **[inf] Upflow:** the near-boiling profile at W10 (margin down to " + fmt(S.loc['WELL-10','min_boiling_margin_degC'],1) + " °C) is the strongest upflow indicator. "
      "Its shallow feed zone sits at " + fmt(fz.loc[(fz.well=='WELL-10'),'mid_mRL'].iloc[0],0) + " mRL.")
    A(f"- **[inf] Outflow:** the hot wells peak at {fmt(S.loc[outflow_wells,'Tmax_elev_mRL'].max(),0)} to {fmt(S.loc[outflow_wells,'Tmax_elev_mRL'].min(),0)} mRL, then cool by {fmt(cool_lo,0)} to {fmt(cool_hi,0)} °C below that level (W2, W3, W6, W7, W9). "
      "This is consistent with lateral outflow at the hot level. Cold downflow would give the same profile; a temperature transect across wells at that depth would separate them.")
    A(f"- **[obs] Recharge:** W4 ({fmt(S.loc['WELL-4','Tmax_degC'],0)} °C at maximum) and W8 ({fmt(S.loc['WELL-8','Tmax_degC'],0)} °C) are cold. W8 lies {fmt(geo.loc['WELL-8','dist_nearest_water_body_m'],0)} m from the nearest water body. "
      "**[inf]** Recharge or downflow is likely at these locations. Their pressures agree with the trend, so the cold does not require a pressure anomaly.")
    A(f"- **[obs] Cap.** The 100 °C isotherm lies at about {fmt(S.loc[hotset,'z_100C_mRL'].min(),0)} to {fmt(S.loc[hotset,'z_100C_mRL'].max(),0)} mRL in the hot wells and at {fmt(S.loc['WELL-1','z_100C_mRL'],0)} mRL in the conductive well W1. "
      f"**[inf]** The hot wells have a 100 °C isotherm close to surface, so they have little cool cap; the conductive well W1 stays below 100 °C down to {fmt(S.loc['WELL-1','z_100C_mRL'],0)} mRL. "
      "Cap base cannot be resolved from the data.")
    A(f"- **[obs] Pressure connection.** Feed-zone pressures agree with the natural-state trend to within {fz.trend_dev_P_bar.abs().max():.1f} bar (largest deviation at W1; table in §4). "
      "**[inf]** Deep pressures are hydraulically connected across the mapped area, so a single natural-state pressure field is a reasonable starting point.")
    A("")
    A("## 6. Recommended calibration set-up")
    A("")
    A("- **Targets.** Use `NS_calibration_targets.csv`. Temperature weight = location weight × confidence factor (high 1.0, medium 0.75, low 0.5). "
      "Pressure targets only at feed zones (`weight_P`), with quality factor HIGH 1.0, MEDIUM 0.8, LOW and Synthetic 0.5, times confidence.")
    A("- **Confidence.** Weight W2 targets fully. W3 and the trend-check wells at medium. W5, W8, W9, W10 and W1 at low. Confidence factors vary across wells, so they change relative weights.")
    w2f2 = fz[(fz.well == "WELL-2") & (fz.feedzone == "FZ2")].iloc[0]
    A(f"- **[inf] Heat in.** Anchor the upflow at the hot cluster. W2 FZ2 (mid {fmt(w2f2.mid_mRL,0)} mRL, {fmt(w2f2.P_bar,1)} bar, quality {w2f2.quality}, confidence {w2f2.ns_confidence}) is the best-constrained hot feed zone.")
    A("- **Heat and mass out.** A lateral outflow at the hot level (about -450 to -750 mRL) and cooling below it; a cold recharge or downflow sink near W4 and W8.")
    A("- **[inf] Permeability.** The data do not constrain permeability directly. The hot feed zones sit in the slotted liner and their pressures follow the natural trend, so they are hydraulically connected to the field. Calibrate against the feed-zone pressures in `feedzones.csv`.")
    A("- **Boundaries.** Pressure set from the trend-check reference (W2 and W3); reference-model boundary conditions are not evidence and were not used.")
    A("")
    A("## 7. Data quality issues and open questions (ranked by effect on calibration)")
    A("")
    A("0. **All ten well designs are labelled 'Synthetic data'** (`provenance_labels.csv`), and W10's feed zone is tagged 'Synthetic'. The dataset may be a synthetic teaching model rather than field data. The results below are conclusions about the dataset as supplied; confirm the data provenance before using them for a field model.")
    A("1. **Upper pressures inconsistent with temperatures** (§4): cased windows of " + ", ".join(upper_wells) + " imply densities up to " + fmt(upper_bad.rho_implied_kgm3.max(),0) + " kg/m³. Ask for the raw gauge and survey records for the top of these wells and whether they were flowing at survey time.")
    A(f"2. **Survey date label** ({C['wells']['WELL-2']['survey_date']:%d/%m/%Y}) conflicts with the gauge drawdown in January 2015. Confirm the date before using it to set a time for the natural-state state.")
    A(f"3. **W9 gauge disagreement** ({fmt(g9['diff'],2)} bar), with the gauge between feed zones. Confirm the gauge depth and whether the gauge reading is a feed-zone pressure.")
    A("4. **W10 'Synthetic' feed-zone tag**, no pre-production gauge and near-boiling margin. Confidence low; feed-zone weights reduced.")
    A("5. **Injectors W5 and W8** are hotter than the injectate at depth: injection cooling could be hidden. Their confidence is low.")
    A(f"6. **Shallow gradients** of {fmt(S.loc[shallow_hot,'shallow_gradient_degC_per_km'].min(),0)} to {fmt(S.loc[shallow_hot,'shallow_gradient_degC_per_km'].max(),0)} °C/km in the hot wells: not conductive. Check for wellbore flow or a wrong shallow depth reference.")
    A(f"7. **Wellhead position.** The wellhead lies below the top of its shallowest model block by up to {fmt(-geo.wellhead_minus_surface_m.min(),0)} m, in {len(wh_below)} of 10 wells "
      "(`well_geometry_check.csv`, `wellhead_minus_surface_m`). Confirm the model top or the well head elevation.")
    A("8. **Flow records start at 01/01/2014** for all wells, so the true first production date is unknown.")
    A("9. **Reference files** contain hard-coded paths from another machine (`gTE14515.auxjson`, `currentGeomSetup.json`); they were not used as evidence.")
    A("10. **Naming and sign conventions:** `Wel10_*` file names; production-flow `dat` files are sign-reversed relative to the JSON.")
    A("")
    A("## 8. Output file index")
    A("")
    idx = [("input_inventory.csv", "role and path of every input file"),
           ("NS_calibration_targets.csv", "one row per survey point with weights"),
           ("feedzones.csv", "feed-zone values, quality factors, trend check"),
           ("T_P_at_layer_centres.csv", "T and P at model layer centres (no extrapolation)"),
           ("well_summary.csv", "per-well characteristics, regime, verdict"),
           ("pressure_consistency_check.csv", "gauge test, verdict and confidence per well"),
           ("hydrostatic_check.csv, hydrostatic_summary.csv", "hydrostatic windows and per-well medians"),
           ("well_geometry_check.csv", "wellhead, model top, faults and water bodies"),
           ("pressure_at_common_elevations.csv", "T and P at common elevations"),
           ("dat_json_crosscheck.csv", "JSON versus plain-column data"),
           ("provenance_labels.csv", "source, comment and quality labels"),
           ("analysis_parameters.csv", "all thresholds used"),
           ("figures/", "well map, profiles, maps, fence section, P–T plot, conceptual sketch")]
    for f_, d_ in idx:
        A(f"- `Output/{f_}`: {d_}")
    A("")
    with open(os.path.join(OUT, "NS_data_interpretation.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")


# ----------------------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------------------
def main():
    print("== Natural-state analysis (CLAUDE.md v3) ==")
    print("Repository root:", ROOT)
    for p_bar, t_c, err in validate_if97():
        print(f"IF97 Region 4 check: Tsat({p_bar:g} bar) = {t_c:.4f} degC; error vs IF97 table {err:+.2e} K (tolerance 0.1 K)")
    print(f"IF97 Tsat(1.01325 bar) = {tsat(1.01325):.2f} degC (reference 99.97); Tsat(100 bar) = {tsat(100.0):.2f} degC (reference 311.0)")

    inv, unassigned = write_inventory()
    print(f"Input inventory: {len(inv)} files; roles: " + ", ".join(f"{k}={v}" for k, v in inv.role.value_counts().items()))
    if unassigned:
        raise SystemExit(f"Unassigned inputs: {unassigned}")

    wells, _ = load_wells()
    feats = load_mapping()
    layers, nlay = read_layers(os.path.join(TM, "Models/TE14515/gTE14515.dat"))
    print(f"Wells: {len(wells)}; grid layers in LAYERS section: {nlay}")

    h_inj = float(np.median([float(w["enth"].iloc[-1]) for w in wells.values() if w["inj"] is not None]))
    inj_T = injectate_temperature(h_inj)
    print(f"Injectate enthalpy {h_inj:.1f} kJ/kg -> injectate temperature {inj_T:.2f} degC (IF97, liquid at 10 bar)")

    xc = crosscheck(wells)
    prov = provenance(wells)
    xc.to_csv(os.path.join(OUT, "dat_json_crosscheck.csv"), index=False)
    prov.to_csv(os.path.join(OUT, "provenance_labels.csv"), index=False)
    print(f"JSON vs dat cross-check: {int((xc.status == 'PASS').sum())} of {len(xc)} series agree")

    analyse(wells, layers)
    gt = gauge_tests(wells)
    ref_names = [n for n in wells if gt[n]["has_gauge"] and gt[n]["passed"] and gt[n]["pre_production"]]
    print("Natural-state reference wells (gauge test passed):", ref_names)
    trend_dev = fz_pressure_trend(wells, ref_names)
    ver = verdicts(wells, gt, trend_dev, inj_T)
    conf_of = {n: ver[n]["confidence"] for n in wells}

    tgt = survey_rows(wells, conf_of)
    fz = feedzone_rows(wells, conf_of)
    summ = well_summary(wells, ver, conf_of)
    pc_rows = []
    for n, w in wells.items():
        g = gt[n]
        flow = g["first_flow"]
        pc_rows.append(dict(
            well=n,
            gauge_elev_mRL=rnd(g["gauge_depth"], 1) if g["has_gauge"] else None,
            survey_P_at_gauge_bar=rnd(g["survey_P"], 2) if g["has_gauge"] else None,
            first_gauge_date=g["first_date"].strftime("%d/%m/%Y") if g["has_gauge"] else None,
            first_gauge_P_bar=rnd(g["first_P"], 2) if g["has_gauge"] else None,
            diff_survey_minus_first_bar=rnd(g["diff"], 2) if g["has_gauge"] else None,
            verdict=ver[n]["verdict"], confidence=conf_of[n],
            first_flow_record_well=(flow.strftime("%d/%m/%Y") if pd.notna(flow) else "none"),
            gauge_minus_first_flow_days=(int((g["first_date"] - flow).days) if g["has_gauge"] and pd.notna(flow) else None),
            gauge_in_feedzone=bool(g["gauge_in_fz"]) if g["has_gauge"] else None,
            trend_max_abs_dev_bar=(rnd(trend_dev.get(n, np.nan), 2) if not g["has_gauge"] else None),
            note=ver[n]["note"]))
    pc = pd.DataFrame(pc_rows)

    hy, excluded = hydrostatic(wells)
    hs, comp = hydro_summary(hy)
    geo = geometry_checks(wells, feats)
    lc = layer_centres(wells, layers)
    ptab = pressure_common(wells)

    tgt.to_csv(os.path.join(OUT, "NS_calibration_targets.csv"), index=False)
    fz.to_csv(os.path.join(OUT, "feedzones.csv"), index=False)
    summ.to_csv(os.path.join(OUT, "well_summary.csv"), index=False)
    pc.to_csv(os.path.join(OUT, "pressure_consistency_check.csv"), index=False)
    hy.to_csv(os.path.join(OUT, "hydrostatic_check.csv"), index=False)
    hs.to_csv(os.path.join(OUT, "hydrostatic_summary.csv"), index=False)
    geo.to_csv(os.path.join(OUT, "well_geometry_check.csv"), index=False)
    lc.to_csv(os.path.join(OUT, "T_P_at_layer_centres.csv"), index=False)
    ptab.to_csv(os.path.join(OUT, "pressure_at_common_elevations.csv"), index=False)

    par_rows = [dict(parameter=k, value=v) for k, v in PAR.items()]
    par_rows += [dict(parameter="injectate_enthalpy_kJkg", value=round(h_inj, 2)),
                 dict(parameter="injectate_T_degC", value=round(inj_T, 2)),
                 dict(parameter="cold_well_component_cased", value=round(comp.get("cased", np.nan), 4)),
                 dict(parameter="cold_well_component_open", value=round(comp.get("open", np.nan), 4))]
    pd.DataFrame(par_rows).to_csv(os.path.join(OUT, "analysis_parameters.csv"), index=False)

    fig_wellmap(wells, feats, summ)
    fig_profiles(wells)
    fig_maps(wells, "maps_temperature", "TA", "Temperature (°C)", "inferno", 0, 300, feats)
    fig_maps(wells, "maps_pressure", "PA", "Pressure (bar)", "viridis", 0, 200, feats)
    fig_section(wells, summ)
    fig_PT(wells)
    fig_conceptual(wells, feats, summ)

    print("\nGauge tests and verdicts:")
    print(pc[["well", "gauge_elev_mRL", "diff_survey_minus_first_bar", "trend_max_abs_dev_bar", "verdict", "confidence"]].to_string(index=False))
    print("\nWell summary:")
    print(summ[["well", "Tmax_degC", "inversion_degC", "shallow_gradient_degC_per_km", "deep_gradient_degC_per_km",
                "water_level_mRL", "regime", "ns_verdict", "ns_confidence"]].to_string(index=False))
    print(f"\nHydrostatic windows used: {len(hy)}; excluded or trimmed notes: {len(excluded)}")
    for e in excluded:
        print("  ", e)
    print("\nHydrostatic summary:")
    print(hs.to_string(index=False))

    C = dict(summ=summ, pc=pc, hs=hs, hy=hy, fz=fz, comp=comp, wells=wells, gt=gt, inv=inv, xc=xc, tgt=tgt,
             ptab=ptab, geo=geo, excluded=excluded, inj_T=inj_T, h_inj=h_inj)
    write_report(C)
    print("\nReport written: Output/NS_data_interpretation.md")


if __name__ == "__main__":
    main()
