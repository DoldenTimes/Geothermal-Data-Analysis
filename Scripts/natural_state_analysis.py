"""Natural-state data interpretation for geothermal model calibration.

Run from the repository root:

    python Scripts/natural_state_analysis.py

Reads only Baseline/ (read-only) and regenerates every output in
Output/Calibration_Input/ and Output/Additional_Output/ (CLAUDE.md sections 2-8), including
the input inventory, all tables, all figures and the report NS_data_interpretation.md.
The process log (Output/Additional_Output/process_log.md) is written by hand, not by this script.
"""
import ast
import datetime as dt
import glob
import io
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402
from iapws.iapws97 import _TSat_P, _PSat_T, _Region1  # noqa: E402

# ------------------------------------------------------------------------------------------
# Paths, console log
# ------------------------------------------------------------------------------------------
ROOT = os.getcwd()
BASE = "Baseline"
CAL = os.path.join("Output", "Calibration_Input")
ADD = os.path.join("Output", "Additional_Output")
CALFIG = os.path.join(CAL, "figures")
ADDFIG = os.path.join(ADD, "figures")
for d in (CAL, ADD, CALFIG, ADDFIG):
    os.makedirs(d, exist_ok=True)

_LOG = io.StringIO()


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    _LOG.write(s + "\n")


def rel(p):
    return os.path.relpath(p, ROOT).replace("\\", "/")


# ------------------------------------------------------------------------------------------
# Parameters (every threshold used; written to parameters_used.csv)
# ------------------------------------------------------------------------------------------
P_ATM = 1.01325          # bar
P_CRIT = 220.64          # bar
G = 9.80665              # m/s2
LIQ_P_MIN = 1.5          # bar, liquid-filled / hydrostatic eligibility
HYD_BM_MIN = 10.0        # degC, hydrostatic eligibility
HYD_WIN = 300.0          # m, hydrostatic window length
ISO_SPAN = 300.0         # m, isothermal section minimum span
ISO_DT = 5.0             # degC, isothermal section maximum range
NEAR_BOIL = 5.0          # degC
INV_MIN = 3.0            # degC
OUTFLOW_INV = 10.0       # degC
COLD_TMAX = 100.0        # degC
UPFLOW_TMAX = 200.0      # degC
SHALLOW_SPAN = 300.0     # m
DEEP_SPAN = 300.0        # m
WELLHEAD_ZONE = 30.0     # m
GAUGE_PRECISION = 0.1    # bar: gauge and survey pressures are reported to 0.1 bar
TOL_MIN = 3.0            # bar
BAND = 500.0             # m, F3 elevation bands
SECTION_HALFWIDTH = 500.0  # m, projection distance for sections
CONF_FACTOR = {"high": 1.0, "medium": 0.75, "low": 0.5}
LOC = {"above_surface": 0.0, "wellhead_or_top_block": 0.3, "feed_zone": 1.0,
       "behind_casing": 0.5, "open_hole_isothermal": 0.5, "open_hole": 0.8}

# ------------------------------------------------------------------------------------------
# IAPWS-IF97 (Region 4 saturation, Region 1 liquid) with verification
# ------------------------------------------------------------------------------------------


def tsat(p_bar):
    """Saturation temperature (degC), IF97 Region 4; NaN above the critical pressure."""
    if p_bar is None or not np.isfinite(p_bar) or p_bar > P_CRIT or p_bar < 0.00611213:
        return float("nan")
    return _TSat_P(p_bar / 10.0) - 273.15


def psat(t_c):
    return _PSat_T(t_c + 273.15) * 10.0


def rho_liq(t_c, p_bar):
    """Liquid density (kg/m3), IF97 Region 1."""
    return 1.0 / _Region1(t_c + 273.15, p_bar / 10.0)["v"]


def h_liq(t_c, p_bar):
    """Liquid enthalpy (kJ/kg), IF97 Region 1."""
    return _Region1(t_c + 273.15, p_bar / 10.0)["h"]


def validate_if97():
    errs = []
    for p_mpa, t_k in [(0.1, 372.755919), (1.0, 453.035632), (10.0, 584.149488)]:
        errs.append(("Tsat", p_mpa, abs(_TSat_P(p_mpa) - t_k)))
    for t_k, p_mpa in [(300.0, 0.353658941e-2), (500.0, 0.263889776e1), (600.0, 0.123443146e2)]:
        e = abs(_PSat_T(t_k) - p_mpa) / p_mpa
        assert e < 1e-8, "IF97 Region 4 Psat verification failed"
    for t_k, p_mpa, v, h in [(300.0, 3.0, 0.100215168e-2, 0.115331273e3),
                             (300.0, 80.0, 0.971180894e-3, 0.184142828e3),
                             (500.0, 3.0, 0.120241800e-2, 0.975542239e3)]:
        r = _Region1(t_k, p_mpa)
        assert abs(r["v"] - v) / v < 1e-8 and abs(r["h"] - h) / h < 1e-8, \
            "IF97 Region 1 verification failed"
    errs.append(("Tsat_ref", 1.01325, abs(tsat(1.01325) - 99.97)))
    errs.append(("Tsat_ref", 100.0, abs(tsat(100.0) - 311.0)))
    worst = max(e for _, _, e in errs)
    if worst > 0.1:
        raise SystemExit(f"IF97 validation failed: max error {worst:.4f} K > 0.1 K")
    log(f"IF97 validation passed: max Tsat error {worst:.4f} K (Tsat(1.01325 bar) = "
        f"{tsat(1.01325):.2f} degC, Tsat(100 bar) = {tsat(100.0):.2f} degC); Region 1 v, h to 1e-8")
    return worst


IF97_MAX_ERR = validate_if97()

# ------------------------------------------------------------------------------------------
# Small helpers
# ------------------------------------------------------------------------------------------


def r(x, d=1):
    """Round for output; NaN stays NaN; -0.0 becomes 0.0."""
    if x is None:
        return np.nan
    try:
        if not np.isfinite(x):
            return np.nan
    except TypeError:
        return x
    v = round(float(x), d)
    return 0.0 if v == 0 else v


def fmt(x, d=1):
    """Format a number exactly as it is written to the CSV tables (for the report)."""
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/a"
    v = r(x, d)
    if d == 0:
        return str(int(round(v)))
    return f"{v:.{d}f}"


def write_csv(df, path):
    df.to_csv(path, index=False, lineterminator="\n", encoding="utf-8")
    log(f"wrote {rel(path)} ({len(df)} rows)")


def parse_date(s):
    return dt.datetime.strptime(s.strip(), "%d/%m/%Y").date()


def iso(d):
    return d.isoformat() if d is not None else ""


def lit(v):
    """Python-literal strings such as "[top, bottom, 'QUALITY']" -> list."""
    if isinstance(v, str):
        return ast.literal_eval(v)
    return v


# ------------------------------------------------------------------------------------------
# 1. Input inventory (CLAUDE.md 2.1)
# ------------------------------------------------------------------------------------------


def classify(path):
    """Role and note for one Baseline file, from its name, location and content."""
    p = path.replace("\\", "/")
    name = os.path.basename(p)
    low = p.lower()
    if "/productionhistory/" in low or "/futurescenarios/" in low or "/production/" in low:
        return "out_of_scope", "production / injection / forecast model input (not read)"
    if name.endswith(".dat"):
        with open(path, errors="ignore") as f:
            head = f.read(200000)
        if "VERTICES" in head and "LAYERS" in head and "SURFA" in head:
            return "geometry_dat", "MULgraph grid file: vertices, grid, layers, surface, well tracks (geometry only)"
        if "/data_dat/" in low:
            return "well_data_dat", "plain-column copies of well data (cross-check only)"
        if "ROCKS" in head:
            return "reference_rock", ("[reference model] TOUGH2 natural-state input: rock types and permeabilities; "
                                      "also holds generators (heat/mass sources) - reference only, not read")
    if name.endswith(".incon"):
        return "reference_results", "[reference model] initial conditions / natural-state results (reference only, not read)"
    if name.endswith(".json"):
        try:
            js = json.load(open(path))
        except Exception:
            js = None
        if isinstance(js, dict):
            if "Well_Name" in js:
                return "well_data_json", "per-well data: track, design, feed zones, surveys, transients, flows"
            if "Areas" in js:
                return "well_groups", "well list and groupings"
            vals = list(js.values())
            if vals and all(isinstance(v, dict) and "Blocks" in v and "Feed_Zones" in v for v in vals):
                return "well_blocks_json", "well-to-block intersections (authority for block assignment)"
            if "/mapping/" in low:
                return "mapping", "polylines of faults, rivers and lake"
            if "rocktype" in name.lower():
                return "reference_rock", "[reference model] rock-type dictionary (reference only)"
            if "boundary" in name.lower():
                return "reference_bc", "[reference model] boundary conditions (reference only, not read)"
            if "Total" in js and isinstance(js["Total"], dict):
                return "out_of_scope", "field-total steam flow and power output (production; not read)"
            return "reference_other", "[reference model] set-up / configuration file (reference only, not read)"
        if isinstance(js, list):
            return "reference_other", "plotting configuration with hard-coded paths from another machine (not read)"
    if name.endswith(".auxjson"):
        return "reference_other", "[reference model] scene/plotting features with hard-coded source paths (not read)"
    return "out_of_scope", "not needed"


def build_inventory():
    files = sorted(os.path.join(dp, f) for dp, _, fs in os.walk(BASE) for f in fs)
    rows = [(classify(f) + (rel(f),)) for f in files]
    df = pd.DataFrame(rows, columns=["role", "note", "path"])[["role", "path", "note"]]
    # collapse directories in which every file has the same role into one directory row
    out, done = [], set()
    for d in sorted({os.path.dirname(p) for p in df.path}, key=lambda s: (s.count("/"), s)):
        sub = df[df.path.str.startswith(d + "/")]
        if len(sub) > 1 and sub.role.nunique() == 1 and d.count("/") >= 3 and \
                not any(d.startswith(x + "/") or d == x for x in done):
            role = sub.role.iloc[0]
            if role in ("well_data_json", "well_data_dat", "reference_other", "out_of_scope"):
                out.append((role, d, f"directory ({len(sub)} files): {sub.note.iloc[0]}"))
                done.add(d)
    covered = lambda p: any(p.startswith(d + "/") for d in done)  # noqa: E731
    for _, rw in df.iterrows():
        if not covered(rw.path):
            out.append((rw.role, rw.path, rw.note))
    inv = pd.DataFrame(out, columns=["role", "path", "note"])
    order = ["well_data_json", "well_data_dat", "well_groups", "mapping", "geometry_dat", "well_blocks_json",
             "reference_rock", "reference_bc", "reference_results", "reference_other", "out_of_scope"]
    inv["_o"] = inv.role.map(order.index)
    inv = inv.sort_values(["_o", "path"]).drop(columns="_o").reset_index(drop=True)
    return inv, df


INV, INV_FILES = build_inventory()
write_csv(INV, os.path.join(ADD, "input_inventory.csv"))


def role_paths(role):
    return list(INV[INV.role == role].path)


def role_files(role):
    return sorted(INV_FILES[INV_FILES.role == role].path)


# ------------------------------------------------------------------------------------------
# 2. Load inputs
# ------------------------------------------------------------------------------------------
GROUPS_JS = json.load(open(role_files("well_groups")[0]))
GROUPS = {}
for area, gs in GROUPS_JS["Areas"].items():
    for gname, wl in gs.items():
        GROUPS[gname] = list(wl)
WELLS = []
for wl in GROUPS.values():
    for w in wl:
        if w not in WELLS:
            WELLS.append(w)
log(f"well_groups: {len(WELLS)} wells in {len(GROUPS)} group(s): {', '.join(WELLS)}")


def load_well(path):
    js = json.load(open(path))
    w = {"name": js["Well_Name"], "aut2": js.get("AUT2_Name", ""), "raw": js, "file": rel(path)}
    w["track"] = np.array([lit(v) for v in js["Well_Track"]["Data"]], float)
    w["track_units"] = js["Well_Track"].get("Units", {})
    w["design"] = [lit(v) for v in js["Well_Design"]["Data"]]
    w["design_source"] = js["Well_Design"].get("Source", "")
    fz = [lit(v) for v in js["Feed_Zones"]["Data"]]
    w["fz"] = [dict(top=float(max(a[0], a[1])), bottom=float(min(a[0], a[1])), quality=str(a[2])) for a in fz]
    w["fz_units"] = js["Feed_Zones"].get("Units", {})
    for key, var in (("Downhole_Temperature", "T"), ("Downhole_Pressure", "P")):
        surveys = js.get(key, {})
        dates = sorted(surveys, key=parse_date)
        w[var + "_dates"] = dates
        if dates:
            s = surveys[dates[0]]
            w[var + "_date"] = parse_date(dates[0])
            w[var + "_meta"] = {k: s.get(k) for k in s if k not in ("Data",)}
            w[var] = {float(k): float(v) for k, v in s["Data"].items()}
    w["gauges"] = {}
    for zkey, s in js.get("Transient_Pressure", {}).items():
        if zkey in ("Units", "Source", "Comment", "Data"):
            continue
        ser = sorted((parse_date(k), float(v)) for k, v in s["Data"].items())
        w["gauges"][float(zkey)] = dict(series=ser, source=s.get("Source", ""), comment=s.get("Comment", ""))
    for key in ("Production_Mass_Flow", "Injection_Mass_Flow", "Transient_Enthalpy", "Injection_Enthalpy"):
        if key in js:
            w[key] = sorted((parse_date(k), float(v)) for k, v in js[key]["Data"].items())
            w[key + "_meta"] = {k: js[key].get(k) for k in ("Source", "Comment")}
    return w


WD = {}
for f in role_files("well_data_json"):
    w = load_well(f)
    WD[w["name"]] = w
missing = [w for w in WELLS if w not in WD]
if missing:
    raise SystemExit(f"wells in groups without data files: {missing}")
WELL_COLOR = {w: plt.get_cmap("tab10")(i % 10) for i, w in enumerate(WELLS)}

# mapping
MAP = json.load(open(role_files("mapping")[0]))
FAULTS, WATER = {}, {}
for name, f in MAP.items():
    data = f["Data"]
    if f.get("type") == "multiplot":
        segs = [np.array(s, float) for s in data]
    else:
        segs = [np.array(data, float)]
    if "fault" in name.lower() or str(f.get("legend", "")).lower() == "fault":
        FAULTS[name] = segs
    else:
        WATER[name] = dict(segs=segs, fill=f.get("type") == "fill", legend=f.get("legend", ""))
log(f"mapping: faults {list(FAULTS)}; water bodies {list(WATER)}")

# geometry (MULgraph)
GEOM_PATH = role_files("geometry_dat")[0]
GEOM_NOTES = []


def read_geometry(path):
    lines = open(path).read().splitlines()
    sec, data = None, {}
    heads = ("VERTICES", "GRID", "CONNECTIONS", "LAYERS", "SURFA", "WELLS")
    for ln in lines:
        key = ln.strip()
        if any(key.startswith(h) and len(key.split()) == 1 for h in heads):
            sec = [h for h in heads if key.startswith(h)][0]
            data[sec] = []
            continue
        if sec:
            data[sec].append(ln)
    # vertices: whitespace split fails (fixed-width, fields run together) -> fixed-width fallback
    verts = {}
    ws_ok = True
    for ln in data["VERTICES"]:
        if not ln.strip():
            break
        p = ln.split()
        if len(p) == 3:
            verts[p[0]] = (float(p[1]), float(p[2]))
        else:
            ws_ok = False
            m = re.match(r"^(.{3})\s*(-?\d+\.\d{2})\s*(-?\d+\.\d{2})\s*$", ln)
            verts[m.group(1).strip()] = (float(m.group(2)), float(m.group(3)))
    if not ws_ok:
        GEOM_NOTES.append("VERTICES section is fixed-width (name 3 chars, x and y 10 chars each, no separating "
                          "whitespace); whitespace parsing fails, so a fixed-width regex was used for this section only")
    cols, cur, n = {}, None, 0
    for ln in data["GRID"]:
        if not ln.strip():
            if cur is None:
                continue
            break
        if cur is None or len(cols[cur]) == n:
            cur = ln[:3].strip()
            n = int(ln.split()[-1])
            cols[cur] = []
        else:
            cols[cur].append(ln.strip())
    layers = {}
    for ln in data["LAYERS"]:
        p = ln.split()
        if len(p) != 3:
            break
        layers[int(p[0])] = (float(p[1]), float(p[2]))
    surf = {}
    for ln in data["SURFA"]:
        p = ln.split()
        if len(p) != 2:
            break
        surf[p[0]] = float(p[1])
    return verts, cols, layers, surf


VERTS, COLS, LAYERS_RAW, SURF = read_geometry(GEOM_PATH)
GRID_TOP = LAYERS_RAW[0][0]
LAYERS = {i: dict(top=LAYERS_RAW[i - 1][0], bottom=LAYERS_RAW[i][0], centre=LAYERS_RAW[i][1])
          for i in sorted(LAYERS_RAW) if i > 0}
COLPOLY = {c: np.array([VERTS[v] for v in vs]) for c, vs in COLS.items()}
COLPATH = {c: MPath(p) for c, p in COLPOLY.items()}
COLCEN = {c: p.mean(0) for c, p in COLPOLY.items()}
COL_WIDTH = len(next(iter(SURF)))
log(f"geometry: {len(VERTS)} vertices, {len(COLS)} columns, {len(LAYERS)} layers "
    f"(top {GRID_TOP} mRL, bottom {LAYERS[max(LAYERS)]['bottom']} mRL), {len(SURF)} surface entries")
for c in COLS:
    if c not in SURF:
        GEOM_NOTES.append(f"column {c} has no surface entry; grid top {GRID_TOP} mRL used")
# layer tiling check
for i in LAYERS:
    if i > 1 and LAYERS[i]["top"] != LAYERS[i - 1]["bottom"]:
        GEOM_NOTES.append(f"layer {i} does not tile")
    if not (LAYERS[i]["bottom"] < LAYERS[i]["centre"] < LAYERS[i]["top"]):
        GEOM_NOTES.append(f"layer {i} centre outside the layer")


def surface_of(col):
    return SURF.get(col, GRID_TOP)


def layer_of(z):
    for i in sorted(LAYERS):
        if LAYERS[i]["bottom"] <= z <= LAYERS[i]["top"]:
            return i
    return None


def surface_layer(col):
    """Top active layer of a column: the layer containing its surface elevation."""
    s = surface_of(col)
    for i in sorted(LAYERS):
        if LAYERS[i]["bottom"] < s <= LAYERS[i]["top"]:
            return i
    return None


def column_at(x, y):
    for c, pth in COLPATH.items():
        if pth.contains_point((x, y)):
            return c
    return None


# model outer boundary: edges used by exactly one column
_edges = {}
for c, vs in COLS.items():
    for a, b in zip(vs, vs[1:] + vs[:1]):
        k = tuple(sorted((a, b)))
        _edges[k] = _edges.get(k, 0) + 1
BOUNDARY_EDGES = [k for k, n in _edges.items() if n == 1]

# well-to-block file
WB_PATH = role_files("well_blocks_json")[0]
WB = json.load(open(WB_PATH))
for wname, d in WB.items():
    for k, v in list(d.items()):
        if isinstance(v, str):
            try:
                d[k] = ast.literal_eval(v)
            except (ValueError, SyntaxError):
                pass


def split_block(b):
    m = re.match(r"^(.*?)(\d+)$", b)
    return m.group(1).strip(), int(m.group(2))


def block_at(well, z):
    """Block from the well-to-block file whose elevation interval contains z (upper block on a boundary)."""
    for b in WB[well]["Blocks"]:
        top, bot = max(b["interval"]), min(b["interval"])
        if bot - 1e-6 <= z <= top + 1e-6:
            return b["block"]
    return None


# ------------------------------------------------------------------------------------------
# 3. Per-well geometry
# ------------------------------------------------------------------------------------------
def track_xy(w, z):
    t = WD[w]["track"]
    zz = t[:, 2][::-1]
    return float(np.interp(z, zz, t[:, 0][::-1])), float(np.interp(z, zz, t[:, 1][::-1]))


for wname, w in WD.items():
    t = w["track"]
    w["wellhead"] = float(t[0, 2])
    w["td"] = float(t[-1, 2])
    casing = [d for d in w["design"] if not bool(d[3])]
    slotted = [d for d in w["design"] if bool(d[3])]
    w["shoe"] = float(min(min(d[0], d[1]) for d in casing)) if casing else w["wellhead"]
    w["slotted"] = [(float(max(d[0], d[1])), float(min(d[0], d[1]))) for d in slotted]
    w["open_top"] = max(s[0] for s in w["slotted"]) if w["slotted"] else w["shoe"]
    w["head_block"] = max(WB[wname]["Blocks"], key=lambda b: max(b["interval"]))["block"]
    w["head_col"] = split_block(w["head_block"])[0]
    w["model_surface"] = surface_of(w["head_col"])
    w["max_dev"] = float(np.max(np.hypot(t[:, 0] - t[0, 0], t[:, 1] - t[0, 1])))
    for f in w["fz"]:
        f["mid"] = 0.5 * (f["top"] + f["bottom"])

# quality tags
QUALITY_ORDER = ["HIGH", "MEDIUM", "LOW"]   # quality grades, highest first
TAGS = sorted({f["quality"] for w in WD.values() for f in w["fz"]})
GRADES = [q for q in QUALITY_ORDER if q in {t.upper() for t in TAGS}]


def q_factor(tag):
    t = tag.upper()
    if t not in GRADES:
        return 0.5, "not a quality grade -> lowest factor"
    if t == GRADES[0]:
        return 1.0, "highest tag in the dataset"
    if t == GRADES[-1]:
        return 0.5, "lowest tag in the dataset"
    return 0.8, "intermediate tag"


TAG_MAP = {t: q_factor(t) for t in TAGS}
log("feed-zone quality tags:", {t: v[0] for t, v in TAG_MAP.items()})

# ------------------------------------------------------------------------------------------
# 4. Survey table
# ------------------------------------------------------------------------------------------
rows = []
for wname in WELLS:
    w = WD[wname]
    zs = sorted(set(w["T"]) | set(w["P"]), reverse=True)
    for z in zs:
        T, P = w["T"].get(z, np.nan), w["P"].get(z, np.nan)
        x, y = track_xy(wname, z)
        rows.append(dict(well=wname, x=x, y=y, elevation_mRL=z, T_degC=T, P_bar=P))
S = pd.DataFrame(rows)
S["Tsat_degC"] = [r(tsat(p), 2) for p in S.P_bar]
S["boiling_margin_degC"] = [r(ts - t, 2) if np.isfinite(ts) else np.nan for ts, t in zip(S.Tsat_degC, S.T_degC)]
S["liquid"] = [(p >= LIQ_P_MIN) and ((bm >= 0) if np.isfinite(bm) else (p > P_CRIT))
               for p, bm in zip(S.P_bar, S.boiling_margin_degC)]


def wsurvey(wname):
    return S[S.well == wname].sort_values("elevation_mRL", ascending=False).reset_index(drop=True)


def interp_w(wname, var, z):
    g = wsurvey(wname)
    zz, vv = g.elevation_mRL.to_numpy()[::-1], g[var].to_numpy()[::-1]
    if z < zz.min() - 1e-9 or z > zz.max() + 1e-9:
        return np.nan
    return float(np.interp(z, zz, vv))


# ------------------------------------------------------------------------------------------
# 5. Integrity checks (CLAUDE.md 4.1 step 1)
# ------------------------------------------------------------------------------------------
INTEG = []   # (well, check, status, detail)


def integ(well, check, ok, detail):
    INTEG.append(dict(well=well, check=check, status="PASS" if ok else "ISSUE", detail=detail))


# 5a. plain-column files: match to wells by content and name
DAT_ROWS, DAT_NAME_NOTES = [], []
REF_DATE = None


def read_dat(path):
    a = np.loadtxt(path, ndmin=2)
    return a


def series_for(w, kind, zkey=None):
    if kind == "downhole_temperature":
        return sorted(w["T"].items())
    if kind == "downhole_pressure":
        return sorted(w["P"].items())
    if kind == "transient_pressure":
        return w["gauges"].get(zkey, {}).get("series")
    key = {"production_mass_flow": "Production_Mass_Flow", "transient_enthalpy": "Transient_Enthalpy",
           "injection_enthalpy": "Injection_Enthalpy", "injection_mass_flow": "Injection_Mass_Flow"}[kind]
    return w.get(key)


dat_files = role_files("well_data_dat")
pattern = re.compile(r"^(?P<prefix>[A-Za-z]+\d+)_(?:(?P<z>-?\d+(?:\.\d+)?)_)?(?P<kind>[a-z_]+)\.dat$")
kinds_time = {"transient_pressure", "production_mass_flow", "transient_enthalpy", "injection_enthalpy",
              "injection_mass_flow"}
# reference date: match a transient series against the JSON (first value of each series)
cands = {}
for f in dat_files:
    m = pattern.match(os.path.basename(f))
    if not m or m.group("kind") not in kinds_time:
        continue
    a = read_dat(f)
    for wname, w in WD.items():
        s = series_for(w, m.group("kind"), float(m.group("z")) if m.group("z") else None)
        if s and len(s) == len(a) and np.allclose(np.abs(a[:, 1]), [v for _, v in s], atol=0.051):
            d0 = s[0][0]
            ref = d0 - dt.timedelta(days=round(a[0, 0] * 365.25))
            cands[ref] = cands.get(ref, 0) + 1
REF_DATE = max(cands, key=lambda k: (cands[k], -k.toordinal())) if cands else None
log(f"plain-column transient files: reference date (t = 0) = {iso(REF_DATE)} "
    f"(agreed by {cands.get(REF_DATE, 0)} of {sum(cands.values())} matched series)")

for f in dat_files:
    base = os.path.basename(f)
    m = pattern.match(base)
    if not m:
        DAT_ROWS.append(dict(file=rel(f), well="", kind="", n_dat=0, n_json=0, max_value_diff=np.nan,
                             max_time_diff_yr=np.nan, status="unmatched name"))
        continue
    kind, prefix = m.group("kind"), m.group("prefix")
    zkey = float(m.group("z")) if m.group("z") else None
    a = read_dat(f)
    best = None
    for wname, w in WD.items():
        s = series_for(w, kind, zkey)
        if not s or len(s) != len(a):
            continue
        if kind.startswith("downhole"):
            dj = dict(s)
            diffs = [abs(dj.get(float(z), np.nan) - v) for v, z in a]
            tdiff = np.nan
        else:
            diffs = [abs(abs(v) - sv) for (t, v), (_, sv) in zip(a, s)]
            tdiff = max(abs(t - (d - REF_DATE).days / 365.25) for (t, _), (d, _) in zip(a, s))
        md = float(np.nanmax(diffs)) if np.isfinite(diffs).all() else np.inf
        name_ok = prefix.upper() == w["aut2"].upper()
        score = (md, not name_ok)
        if best is None or score < best[0]:
            best = (score, wname, md, tdiff, name_ok)
    if best is None:
        DAT_ROWS.append(dict(file=rel(f), well="", kind=kind, n_dat=len(a), n_json=0, max_value_diff=np.nan,
                             max_time_diff_yr=np.nan, status="no JSON series of the same length"))
        continue
    _, wname, md, tdiff, name_ok = best
    sign = ""
    if kind == "production_mass_flow" and (a[:, 1] < 0).all():
        sign = "; production stored as negative in .dat, positive in JSON"
    status = "match" if md <= 0.051 else "MISMATCH"
    DAT_ROWS.append(dict(file=rel(f), well=wname, kind=kind, n_dat=len(a), n_json=len(series_for(WD[wname], kind, zkey)),
                         max_value_diff=r(md, 3), max_time_diff_yr=r(tdiff, 3),
                         status=status + ("" if name_ok else "; file prefix differs from AUT2 name") + sign))
    if not re.match(r"^Well\d+$", prefix):
        DAT_NAME_NOTES.append(f"{base}: prefix '{prefix}' does not follow the 'Well<N>' pattern of the other files "
                              f"(matched to {wname} by content; AUT2 name {WD[wname]['aut2']})")
DATX = pd.DataFrame(DAT_ROWS)
write_csv(DATX, os.path.join(ADD, "dat_json_crosscheck.csv"))
for wname in WELLS:
    have = set(DATX[DATX.well == wname].kind)
    need = {"downhole_temperature", "downhole_pressure"}
    w = WD[wname]
    need |= {"transient_pressure"} if w["gauges"] else set()
    for k, kk in (("Production_Mass_Flow", "production_mass_flow"), ("Transient_Enthalpy", "transient_enthalpy"),
                  ("Injection_Enthalpy", "injection_enthalpy"), ("Injection_Mass_Flow", "injection_mass_flow")):
        if k in w:
            need.add(kk)
    sub = DATX[DATX.well == wname]
    ok = need <= have and (sub.status.str.startswith("match")).all()
    integ(wname, "dat_vs_json", ok, f"{len(sub)} .dat files; kinds matched {sorted(have)}"
          + ("" if need <= have else f"; missing {sorted(need - have)}"))

# 5b. surveys, feed zones and design against the track
for wname in WELLS:
    w = WD[wname]
    g = wsurvey(wname)
    zt, zb = w["wellhead"], w["td"]
    out = g[(g.elevation_mRL > zt + 1e-6) | (g.elevation_mRL < zb - 1e-6)]
    integ(wname, "survey_within_track", out.empty,
          f"survey {g.elevation_mRL.max()} to {g.elevation_mRL.min()} mRL; track {zt} to {zb} mRL"
          + ("" if out.empty else f"; outside: {list(out.elevation_mRL)}"))
    gap = zt - g.elevation_mRL.max()
    integ(wname, "survey_starts_at_wellhead", gap <= 1.0,
          f"first survey point {g.elevation_mRL.max()} mRL is {r(gap, 1)} m below the wellhead {zt} mRL")
    if set(w["T"]) != set(w["P"]):
        integ(wname, "T_P_same_elevations", False, "T and P surveys at different elevations")
    else:
        integ(wname, "T_P_same_elevations", True, f"{len(w['T'])} paired points, survey date "
              f"{iso(w['T_date'])}")
    for i, f in enumerate(w["fz"], 1):
        ok_t = zb - 1e-6 <= f["bottom"] and f["top"] <= zt + 1e-6
        integ(wname, f"feedzone{i}_within_track", ok_t, f"{f['top']} to {f['bottom']} mRL")
        ok_s = any(top + 1e-6 >= f["top"] and bot - 1e-6 <= f["bottom"] for top, bot in w["slotted"])
        integ(wname, f"feedzone{i}_within_slotted", ok_s,
              f"{f['top']} to {f['bottom']} mRL; slotted {w['slotted']}")
        has_mid = any(abs(z - f["mid"]) < 1e-6 for z in w["T"])
        integ(wname, f"feedzone{i}_survey_point_at_mid", has_mid, f"mid {f['mid']} mRL")
    d_top = max(max(d[0], d[1]) for d in w["design"])
    d_bot = min(min(d[0], d[1]) for d in w["design"])
    integ(wname, "design_top_at_wellhead", abs(d_top - zt) < 0.5, f"design top {d_top}, wellhead {zt}")
    integ(wname, "design_bottom_at_TD", abs(d_bot - zb) < 0.5, f"design bottom {d_bot}, track TD {zb}")
    laps = []
    ds = sorted(w["design"], key=lambda d: -max(d[0], d[1]))
    for a, b in zip(ds, ds[1:]):
        if min(b[0], b[1]) < min(a[0], a[1]) and max(b[0], b[1]) >= min(a[0], a[1]):
            laps.append(f"{b[4]} lapped {max(b[0], b[1]) - min(a[0], a[1]):+.0f} m into {a[4]}")
        else:
            laps.append(f"gap/no lap between {a[4]} and {b[4]}")
    integ(wname, "casing_liner_lap", all("lapped" in s for s in laps), "; ".join(laps))
    # well-to-block file vs track
    wb = WB[wname]["Blocks"]
    integ(wname, "blockfile_spans_track",
          abs(max(max(b["interval"]) for b in wb) - zt) < 0.5 and abs(min(min(b["interval"]) for b in wb) - zb) < 0.5,
          f"block file {max(max(b['interval']) for b in wb)} to {min(min(b['interval']) for b in wb)} mRL")
    # geometry file well track vs JSON track (WELLS section of the grid file)
INTEGDF = pd.DataFrame(INTEG)
write_csv(INTEGDF, os.path.join(ADD, "integrity_checks.csv"))
for _, rw in INTEGDF[INTEGDF.status != "PASS"].iterrows():
    log(f"integrity ISSUE {rw.well} {rw.check}: {rw.detail}")

# provenance fields
PROV = []
for wname in WELLS:
    w, js = WD[wname], WD[wname]["raw"]
    for key, v in js.items():
        if isinstance(v, dict):
            if key in ("Downhole_Temperature", "Downhole_Pressure"):
                for d, s in v.items():
                    PROV.append(dict(well=wname, item=f"{key} {iso(parse_date(d))}", field="Source", value=s.get("Source", "")))
                    PROV.append(dict(well=wname, item=f"{key} {iso(parse_date(d))}", field="Comment", value=s.get("Comment", "")))
                    PROV.append(dict(well=wname, item=f"{key} {iso(parse_date(d))}", field="StableFlag", value=str(s.get("StableFlag", ""))))
            elif key == "Transient_Pressure":
                for zk, s in v.items():
                    PROV.append(dict(well=wname, item=f"{key} {zk} mRL", field="Source", value=s.get("Source", "")))
                    PROV.append(dict(well=wname, item=f"{key} {zk} mRL", field="Comment", value=s.get("Comment", "")))
            else:
                for fld in ("Source", "Comment"):
                    if fld in v:
                        PROV.append(dict(well=wname, item=key, field=fld, value=v[fld]))
                if "Units" in v:
                    PROV.append(dict(well=wname, item=key, field="Units", value=json.dumps(v["Units"], sort_keys=True)))
    for i, f in enumerate(w["fz"], 1):
        PROV.append(dict(well=wname, item=f"Feed zone {i}", field="Quality", value=f["quality"]))
PROVDF = pd.DataFrame(PROV)
write_csv(PROVDF, os.path.join(ADD, "provenance_fields.csv"))

# ------------------------------------------------------------------------------------------
# 6. Point attributes: blocks, casing, feed zones, surface, isothermal sections
# ------------------------------------------------------------------------------------------
S["block"] = [block_at(wl, z) for wl, z in zip(S.well, S.elevation_mRL)]
S["column"] = [split_block(b)[0] if b else "" for b in S.block]
S["layer"] = [split_block(b)[1] if b else np.nan for b in S.block]
S["model_surface_mRL"] = [surface_of(c) if c else np.nan for c in S.column]
S["above_model_surface"] = S.elevation_mRL > S.model_surface_mRL + 1e-9
S["behind_casing"] = [z >= WD[wl]["shoe"] - 1e-9 for wl, z in zip(S.well, S.elevation_mRL)]


def fz_index(wl, z):
    for i, f in enumerate(WD[wl]["fz"], 1):
        if f["bottom"] - 1e-9 <= z <= f["top"] + 1e-9:
            return i
    return 0


S["feedzone_no"] = [fz_index(wl, z) for wl, z in zip(S.well, S.elevation_mRL)]
S["in_feedzone"] = S.feedzone_no > 0
S["top_block"] = [bool(c) and int(lay) == surface_layer(c) for c, lay in zip(S.column, S.layer.fillna(-1))]
S["near_wellhead"] = [WD[wl]["wellhead"] - z <= WELLHEAD_ZONE + 1e-9 for wl, z in zip(S.well, S.elevation_mRL)]


def isothermal_flags(g, shoe):
    """Open-hole points inside an isothermal section (>= ISO_SPAN m, range <= ISO_DT degC)."""
    z, t = g.elevation_mRL.to_numpy(), g.T_degC.to_numpy()
    flag = np.zeros(len(z), bool)
    sections = []
    idx = [i for i in range(len(z)) if z[i] < shoe - 1e-9]
    for a in range(len(idx)):
        best = None
        for b in range(a, len(idx)):
            seg = t[idx[a]:idx[b] + 1]
            if seg.max() - seg.min() > ISO_DT + 1e-9:
                break
            if z[idx[a]] - z[idx[b]] >= ISO_SPAN - 1e-9:
                best = b
        if best is not None:
            flag[idx[a]:idx[best] + 1] = True
    # merge flagged runs into sections
    i = 0
    while i < len(z):
        if flag[i]:
            j = i
            while j + 1 < len(z) and flag[j + 1]:
                j += 1
            sections.append((z[i], z[j], float(t[i:j + 1].min()), float(t[i:j + 1].max())))
            i = j + 1
        else:
            i += 1
    return flag, sections


S["isothermal"] = False
ISO_SECTIONS = {}
for wname in WELLS:
    m = S.well == wname
    g = S[m].sort_values("elevation_mRL", ascending=False)
    fl, secs = isothermal_flags(g.reset_index(drop=True), WD[wname]["shoe"])
    S.loc[g.index, "isothermal"] = fl
    ISO_SECTIONS[wname] = secs
log("isothermal open-hole sections:", {k: [(a, b) for a, b, _, _ in v] for k, v in ISO_SECTIONS.items() if v})

# ------------------------------------------------------------------------------------------
# 7. Natural-state timing test (CLAUDE.md 4.1 step 2)
# ------------------------------------------------------------------------------------------
def first_flow(series):
    for d, v in series or []:
        if abs(v) > 0:
            return d
    return None


for wname in WELLS:
    w = WD[wname]
    w["first_prod"] = first_flow(w.get("Production_Mass_Flow"))
    w["first_inj"] = first_flow(w.get("Injection_Mass_Flow"))
    w["first_flow"] = min([d for d in (w["first_prod"], w["first_inj"]) if d], default=None)
    w["kind"] = "injector" if w["first_inj"] else ("producer" if w["first_prod"] else "no flow record")
FIELD_FIRST_PROD = min(d for d in (WD[w]["first_prod"] for w in WELLS) if d)
FIELD_FIRST_INJ = min((d for d in (WD[w]["first_inj"] for w in WELLS) if d), default=None)
SURVEY_DATES = sorted({WD[w]["T_date"] for w in WELLS} | {WD[w]["P_date"] for w in WELLS})
log(f"first production anywhere in the field: {iso(FIELD_FIRST_PROD)}; first injection: {iso(FIELD_FIRST_INJ)}; "
    f"survey dates: {[iso(d) for d in SURVEY_DATES]}")

# injectate temperature from injectate enthalpy (IF97 Region 1, at an assumed 10 bar)
INJ_P_ASSUMED = 10.0


def t_from_h(h, p_bar=INJ_P_ASSUMED):
    lo, hi = 1.0, 300.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if h_liq(mid, p_bar) < h:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


for wname in WELLS:
    w = WD[wname]
    if w.get("Injection_Enthalpy"):
        hs = np.array([v for _, v in w["Injection_Enthalpy"]])
        w["inj_h"] = float(np.median(hs))
        w["inj_h_range"] = (float(hs.min()), float(hs.max()))
        w["inj_T"] = t_from_h(w["inj_h"])

GAUGE_ROWS = []
for wname in WELLS:
    w = WD[wname]
    if not w["gauges"]:
        continue
    for zg, gdat in sorted(w["gauges"].items(), reverse=True):
        ser = gdat["series"]
        d0, p0 = ser[0]
        ps = interp_w(wname, "P_bar", zg)
        diff = ps - p0
        pre = d0 <= FIELD_FIRST_PROD and (w["first_flow"] is None or d0 <= w["first_flow"])
        # survey elevation at which the survey pressure equals the first gauge reading
        g = wsurvey(wname)
        zz, pp = g.elevation_mRL.to_numpy(), g.P_bar.to_numpy()
        z_match = float(np.interp(p0, pp, zz)) if pp.min() <= p0 <= pp.max() else np.nan
        # along-hole hypothesis: gauge "elevation" = wellhead - measured depth
        t = w["track"]
        md = np.concatenate([[0], np.cumsum(np.sqrt(np.sum(np.diff(t, axis=0) ** 2, axis=1)))])
        md_g = w["wellhead"] - zg
        z_md = float(np.interp(md_g, md, t[:, 2])) if md_g <= md[-1] else np.nan
        p_md = interp_w(wname, "P_bar", z_md) if np.isfinite(z_md) else np.nan
        GAUGE_ROWS.append(dict(well=wname, gauge_elev_mRL=zg, survey_P_at_gauge_bar=r(ps, 2),
                               first_gauge_date=d0, first_gauge_P_bar=p0, diff=diff, pre=pre,
                               z_match=z_match, z_md=z_md, p_md=p_md, series=ser))

# trend of gauge-verified feed zones
VERIFIED = [g["well"] for g in GAUGE_ROWS if g["pre"] and abs(g["diff"]) <= GAUGE_PRECISION + 1e-9]
VERIFIED = [w for w in WELLS if w in VERIFIED]
fzp = []
for wname in VERIFIED:
    for f in WD[wname]["fz"]:
        fzp.append((f["mid"], interp_w(wname, "P_bar", f["mid"]), wname))
_z = np.array([a for a, _, _ in fzp])
_p = np.array([b for _, b, _ in fzp])
TREND_B, TREND_A = np.polyfit(_z, _p, 1)   # P = A + B z
TREND_RES = _p - (TREND_A + TREND_B * _z)
TREND_RMS = float(np.sqrt(np.mean(TREND_RES ** 2)))
TOL = max(TOL_MIN, 2 * TREND_RMS)
TREND_ZRANGE = (float(_z.min()), float(_z.max()))
TREND_RHO = -TREND_B * 1e5 / G
log(f"pressure trend of gauge-verified feed zones ({', '.join(VERIFIED)}; n = {len(fzp)}): "
    f"P = {TREND_A:.3f} + {TREND_B:.5f} z bar (z in mRL); RMS {TREND_RMS:.3f} bar; tolerance {TOL:.2f} bar")


def trend(z):
    return TREND_A + TREND_B * z


NS = {}
for wname in WELLS:
    w = WD[wname]
    gs = [g for g in GAUGE_ROWS if g["well"] == wname]
    fzres = [interp_w(wname, "P_bar", f["mid"]) - trend(f["mid"]) for f in w["fz"]]
    maxres = max((abs(x) for x in fzres), default=np.nan)
    on_trend = bool(fzres) and maxres <= TOL + 1e-9
    match = [g for g in gs if g["pre"] and abs(g["diff"]) <= GAUGE_PRECISION + 1e-9]
    miss = [g for g in gs if g["pre"] and abs(g["diff"]) > GAUGE_PRECISION + 1e-9]
    fz_T = [interp_w(wname, "T_degC", f["mid"]) for f in w["fz"]]
    inj_hot = w["kind"] == "injector" and "inj_T" in w and max(fz_T) > w["inj_T"]
    # injection signature: open-hole isothermal section within 10 degC of the injectate temperature
    inj_sig = w["kind"] == "injector" and "inj_T" in w and any(
        abs(0.5 * (lo + hi) - w["inj_T"]) <= 10 for _, _, lo, hi in ISO_SECTIONS[wname])
    if match:
        v, c = "natural state", "high"
        why = (f"survey P at the gauge ({fmt(match[0]['gauge_elev_mRL'], 1)} mRL) "
               f"{fmt(match[0]['survey_P_at_gauge_bar'], 2)} bar vs first gauge reading "
               f"{fmt(match[0]['first_gauge_P_bar'], 1)} bar on {iso(match[0]['first_gauge_date'])} "
               f"(difference {fmt(match[0]['diff'], 2)} bar, within the {GAUGE_PRECISION} bar reading precision)")
    elif inj_hot:
        v, c = "probably natural state at best", "low"
        why = (f"injector; feed-zone formation temperature {fmt(max(fz_T), 1)} degC is hotter than the injectate "
               f"({fmt(w['inj_T'], 1)} degC), so injection cooling could be hidden")
    elif miss and on_trend:
        g0 = miss[0]
        v, c = "probably natural state", "medium"
        why = (f"gauge test misses: survey {fmt(g0['survey_P_at_gauge_bar'], 2)} bar vs first reading "
               f"{fmt(g0['first_gauge_P_bar'], 1)} bar (difference {fmt(g0['diff'], 2)} bar); feed-zone pressures "
               f"within {fmt(maxres, 2)} bar of the verified trend (tolerance {fmt(TOL, 1)} bar)")
    elif not match and not miss and on_trend and not inj_sig:
        v, c = "probably natural state", "medium"
        reason = ("gauge series starts after field production began" if gs else "no pressure gauge")
        why = (f"{reason}; feed-zone pressures within {fmt(maxres, 2)} bar of the verified trend "
               f"(tolerance {fmt(TOL, 1)} bar)")
    else:
        v, c = "possibly disturbed", "low"
        why = (f"feed-zone pressures off the verified trend by up to {fmt(maxres, 2)} bar "
               f"(tolerance {fmt(TOL, 1)} bar)" if not on_trend else "profile signature of injection")
    if w["kind"] == "injector" and "inj_T" in w and not inj_hot:
        why += (f"; injector with feed zones colder ({fmt(max(fz_T), 1)} degC) than the injectate "
                f"({fmt(w['inj_T'], 1)} degC) and no isothermal section near the injectate temperature: "
                "evidence against an injection-affected survey")
    NS[wname] = dict(verdict=v, confidence=c, reason=why, fz_res=fzres, max_fz_res=maxres, inj_sig=inj_sig,
                     inj_hot=inj_hot)
    log(f"natural state {wname}: {v} ({c}) - {why}")

# pressure_consistency_check.csv
pc_rows = []
for wname in WELLS:
    w = WD[wname]
    gs = [g for g in GAUGE_ROWS if g["well"] == wname]
    base = dict(verdict=NS[wname]["verdict"], confidence=NS[wname]["confidence"],
                well_kind=w["kind"], survey_date=iso(w["P_date"]),
                own_first_production_or_injection_date=iso(w["first_flow"]),
                field_first_production_date=iso(FIELD_FIRST_PROD),
                max_feedzone_residual_from_trend_bar=r(NS[wname]["max_fz_res"], 2), tolerance_bar=r(TOL, 2),
                reason=NS[wname]["reason"])
    if gs:
        g = gs[0]
        row = dict(well=wname, gauge_elev_mRL=g["gauge_elev_mRL"], survey_P_at_gauge_bar=g["survey_P_at_gauge_bar"],
                   first_gauge_date=iso(g["first_gauge_date"]), first_gauge_P_bar=g["first_gauge_P_bar"],
                   diff_survey_minus_first_bar=r(g["diff"], 2))
        row.update(base)
        row.update(gauge_start_vs_own_first_flow_days=(g["first_gauge_date"] - w["first_flow"]).days if w["first_flow"] else np.nan,
                   gauge_start_vs_field_first_production_days=(g["first_gauge_date"] - FIELD_FIRST_PROD).days,
                   pre_production_reading=g["pre"],
                   survey_elev_matching_first_reading_mRL=r(g["z_match"], 1),
                   match_elev_minus_gauge_elev_m=r(g["z_match"] - g["gauge_elev_mRL"], 1),
                   elev_if_gauge_depth_is_along_hole_mRL=r(g["z_md"], 1),
                   survey_P_at_along_hole_elev_bar=r(g["p_md"], 2),
                   second_gauge_date=iso(g["series"][1][0]) if len(g["series"]) > 1 else "",
                   second_gauge_P_bar=g["series"][1][1] if len(g["series"]) > 1 else np.nan,
                   drop_first_to_second_bar=r(g["series"][0][1] - g["series"][1][1], 1) if len(g["series"]) > 1
                   else np.nan,
                   last_gauge_date=iso(g["series"][-1][0]), last_gauge_P_bar=g["series"][-1][1])
    else:
        row = dict(well=wname, gauge_elev_mRL=np.nan, survey_P_at_gauge_bar=np.nan, first_gauge_date="",
                   first_gauge_P_bar=np.nan, diff_survey_minus_first_bar=np.nan)
        row.update(base)
    if w.get("inj_T") is not None:
        row.update(injectate_enthalpy_kJkg=r(w["inj_h"], 1), injectate_T_degC=r(w["inj_T"], 1))
    pc_rows.append(row)
PC = pd.DataFrame(pc_rows)
lead = ["well", "gauge_elev_mRL", "survey_P_at_gauge_bar", "first_gauge_date", "first_gauge_P_bar",
        "diff_survey_minus_first_bar", "verdict", "confidence"]
PC = PC[lead + [c for c in PC.columns if c not in lead]]
for c in ("gauge_start_vs_own_first_flow_days", "gauge_start_vs_field_first_production_days"):
    if c in PC:
        PC[c] = PC[c].astype("Int64")
write_csv(PC, os.path.join(ADD, "pressure_consistency_check.csv"))

# ------------------------------------------------------------------------------------------
# 8. Weights (CLAUDE.md 5.5)
# ------------------------------------------------------------------------------------------


def loc_weight(row):
    if row["above_model_surface"]:
        return LOC["above_surface"], "above model surface"
    if row["near_wellhead"] or row["top_block"]:
        return LOC["wellhead_or_top_block"], ("within 30 m of wellhead" if row["near_wellhead"]
                                               else "top block of its column")
    if row["in_feedzone"]:
        return LOC["feed_zone"], f"feed zone {row['feedzone_no']}"
    if row["behind_casing"]:
        return LOC["behind_casing"], "behind casing"
    if row["isothermal"]:
        return LOC["open_hole_isothermal"], "open hole, isothermal section"
    return LOC["open_hole"], "open hole"


S["ns_confidence"] = [NS[w]["confidence"] for w in S.well]
S["ns_verdict"] = [NS[w]["verdict"] for w in S.well]
lw = [loc_weight(rw) for _, rw in S.iterrows()]
S["loc_weight"] = [a for a, _ in lw]
S["weight"] = [r(a * CONF_FACTOR[c], 4) for (a, _), c in zip(lw, S.ns_confidence)]
S["weight_reason"] = [f"{b} ({a}) x {c} confidence ({CONF_FACTOR[c]})" for (a, b), c in zip(lw, S.ns_confidence)]


def wp(row):
    if not row["in_feedzone"] or row["above_model_surface"]:
        return 0.0
    q = WD[row["well"]]["fz"][row["feedzone_no"] - 1]["quality"]
    return r(TAG_MAP[q][0] * CONF_FACTOR[row["ns_confidence"]], 4)


S["weight_P"] = [wp(rw) for _, rw in S.iterrows()]

TGT_COLS = ["well", "x", "y", "elevation_mRL", "T_degC", "P_bar", "Tsat_degC", "boiling_margin_degC", "column",
            "layer", "block", "in_feedzone", "behind_casing", "above_model_surface", "ns_confidence", "weight",
            "weight_reason", "weight_P"]
TGT = S.copy()
TGT["x"] = TGT.x.round(1)
TGT["y"] = TGT.y.round(1)
TGT["layer"] = TGT.layer.astype("Int64")
TGT["feedzone"] = TGT.feedzone_no.replace(0, np.nan).astype("Int64")
TGT["liquid_filled"] = TGT.liquid
TGT["isothermal_section"] = TGT.isothermal
TGT["near_wellhead_or_top_block"] = TGT.near_wellhead | TGT.top_block
TGT["P_trend_bar"] = [r(trend(z), 2) for z in TGT.elevation_mRL]
TGT["P_minus_trend_bar"] = [r(p - trend(z), 2) for p, z in zip(TGT.P_bar, TGT.elevation_mRL)]
TGT = TGT[TGT_COLS + ["feedzone", "model_surface_mRL", "liquid_filled", "isothermal_section",
                      "near_wellhead_or_top_block", "loc_weight", "P_trend_bar", "P_minus_trend_bar", "ns_verdict"]]
write_csv(TGT, os.path.join(CAL, "NS_calibration_targets.csv"))

# ------------------------------------------------------------------------------------------
# 9. Feed zones
# ------------------------------------------------------------------------------------------
FZ_ROWS = []
for wname in WELLS:
    w = WD[wname]
    fzwb = list(WB[wname]["Feed_Zones"].values())
    for i, f in enumerate(w["fz"], 1):
        T = interp_w(wname, "T_degC", f["mid"])
        P = interp_w(wname, "P_bar", f["mid"])
        ts = r(tsat(P), 2)
        cands_b = list(fzwb[i - 1].keys()) if i <= len(fzwb) else []
        b_mid = block_at(wname, f["mid"])
        blk = b_mid if b_mid in cands_b else (max(cands_b, key=lambda b: fzwb[i - 1][b].get("Fraction", 0))
                                              if cands_b else b_mid)
        col = split_block(blk)[0]
        pt = dict(above_model_surface=f["mid"] > surface_of(col),
                  near_wellhead=w["wellhead"] - f["mid"] <= WELLHEAD_ZONE,
                  top_block=split_block(blk)[1] == surface_layer(col), in_feedzone=True, feedzone_no=i,
                  behind_casing=f["mid"] >= w["shoe"], isothermal=False)
        lwt, _ = loc_weight(pt)
        cf = CONF_FACTOR[NS[wname]["confidence"]]
        qf = TAG_MAP[f["quality"]][0]
        xy = track_xy(wname, f["mid"])
        FZ_ROWS.append(dict(well=wname, feedzone=i, top_mRL=f["top"], bottom_mRL=f["bottom"], mid_mRL=f["mid"],
                            quality=f["quality"], T_degC=r(T, 2), P_bar=r(P, 2), Tsat_degC=ts,
                            boiling_margin_degC=r(ts - T, 2) if np.isfinite(ts) else np.nan, block=blk,
                            ns_confidence=NS[wname]["confidence"], weight_T=r(lwt * cf, 4), weight_P=r(qf * cf, 4),
                            thickness_m=f["top"] - f["bottom"], quality_factor=qf,
                            blockfile_blocks="/".join(cands_b), P_trend_bar=r(trend(f["mid"]), 2),
                            P_minus_trend_bar=r(P - trend(f["mid"]), 2),
                            x=r(xy[0], 1), y=r(xy[1], 1), ns_verdict=NS[wname]["verdict"]))
FZ = pd.DataFrame(FZ_ROWS)
write_csv(FZ, os.path.join(CAL, "feedzones.csv"))

# ------------------------------------------------------------------------------------------
# 10. Per-well characteristics (CLAUDE.md 4.4), from the targets table
# ------------------------------------------------------------------------------------------


def crossing(z, t, x):
    for i in range(len(t) - 1):
        if t[i] < x <= t[i + 1]:
            return z[i] + (x - t[i]) * (z[i + 1] - z[i]) / (t[i + 1] - t[i])
    return np.nan


def summarise(wname):
    g = TGT[TGT.well == wname].sort_values("elevation_mRL", ascending=False)
    z, t = g.elevation_mRL.to_numpy(float), g.T_degC.to_numpy(float)
    p, bm = g.P_bar.to_numpy(float), g.boiling_margin_degC.to_numpy(float)
    liq = g.liquid_filled.to_numpy(bool)
    e = dict(well=wname)
    i = int(np.argmax(t))
    e["Tmax_degC"], e["Tmax_elev_mRL"] = t[i], z[i]
    e["T_at_TD_degC"], e["TD_mRL"] = t[-1], z[-1]
    run_max, run_z, best = -np.inf, np.nan, (0.0, np.nan, np.nan)
    for zi, ti in zip(z, t):
        if ti > run_max:
            run_max, run_z = ti, zi
        if run_max - ti > best[0]:
            best = (run_max - ti, run_z, zi)
    inv = best if best[0] >= INV_MIN else (np.nan, np.nan, np.nan)
    e["inversion_degC"], e["inversion_top_mRL"], e["inversion_base_mRL"] = r(inv[0], 1), inv[1], inv[2]
    zl, tl, pl, bml = z[liq], t[liq], p[liq], bm[liq]
    e["shallow_gradient_degC_per_km"] = np.nan
    if len(zl):
        sel = zl >= zl[0] - SHALLOW_SPAN
        if sel.sum() >= 3:
            e["shallow_gradient_degC_per_km"] = r(-np.polyfit(zl[sel], tl[sel], 1)[0] * 1000, 1)
    fzb = min((f["bottom"] for f in WD[wname]["fz"]), default=np.nan)
    bounds = [b for b in (inv[2], fzb) if np.isfinite(b)]
    sel = z <= min(bounds) + 1e-9 if bounds else np.ones(len(z), bool)
    e["deep_gradient_degC_per_km"] = np.nan
    n_deep, span_deep = int(sel.sum()), (float(z[sel].max() - z[sel].min()) if sel.any() else 0.0)
    if n_deep >= 3 and span_deep >= DEEP_SPAN:
        e["deep_gradient_degC_per_km"] = r(-np.polyfit(z[sel], t[sel], 1)[0] * 1000, 1)
    e["z_100C_mRL"], e["z_200C_mRL"] = r(crossing(z, t, 100.0), 1), r(crossing(z, t, 200.0), 1)
    ok = np.isfinite(bml)
    if ok.any():
        j = int(np.argmin(np.where(ok, bml, np.inf)))
        e["min_boiling_margin_degC"], e["min_boiling_margin_elev_mRL"] = bml[j], zl[j]
    else:
        e["min_boiling_margin_degC"] = e["min_boiling_margin_elev_mRL"] = np.nan
    e["water_level_mRL"], wl_extra = np.nan, np.nan
    if len(zl) >= 3:
        b, a = np.polyfit(zl[:3], pl[:3], 1)
        if b != 0:
            e["water_level_mRL"] = r((P_ATM - a) / b, 1)
            wl_extra = r(e["water_level_mRL"] - z[0], 1)
    fzs = FZ[FZ.well == wname]
    if len(fzs):
        cand = fzs.sort_values(["quality_factor", "thickness_m", "mid_mRL"], ascending=[False, False, True])
        pv = cand.iloc[0]
        e["pivot_elev_mRL"], e["pivot_P_bar"] = pv.mid_mRL, pv.P_bar
        pivot_fz = int(pv.feedzone)
    else:
        e["pivot_elev_mRL"] = e["pivot_P_bar"] = np.nan
        pivot_fz = np.nan
    if ok.any() and np.nanmin(bml) <= 1:
        reg = "two-phase"
    elif e["Tmax_degC"] < COLD_TMAX:
        reg = "recharge/downflow"
    elif np.isfinite(inv[0]) and inv[0] >= OUTFLOW_INV:
        reg = "outflow"
    elif e["Tmax_degC"] >= UPFLOW_TMAX and not np.isfinite(inv[0]):
        reg = "upflow"
    else:
        reg = "conductive"
    e["regime"] = reg
    e["ns_verdict"], e["ns_confidence"] = NS[wname]["verdict"], NS[wname]["confidence"]
    e["near_boiling"] = bool(ok.any() and np.nanmin(bml) < NEAR_BOIL)
    e["deep_fit_n"], e["deep_fit_span_m"] = n_deep, r(span_deep, 1)
    e["water_level_extrapolation_m"] = wl_extra
    e["pivot_feedzone"] = pivot_fz
    e["shoe_mRL"] = WD[wname]["shoe"]
    e["wellhead_mRL"] = WD[wname]["wellhead"]
    e["well_kind"] = WD[wname]["kind"]
    return e


WS = pd.DataFrame([summarise(w) for w in WELLS])
WS["pivot_feedzone"] = WS.pivot_feedzone.astype("Int64")
WS_COLS = ["well", "Tmax_degC", "Tmax_elev_mRL", "T_at_TD_degC", "TD_mRL", "inversion_degC", "inversion_top_mRL",
           "inversion_base_mRL", "shallow_gradient_degC_per_km", "deep_gradient_degC_per_km", "z_100C_mRL",
           "z_200C_mRL", "min_boiling_margin_degC", "min_boiling_margin_elev_mRL", "water_level_mRL",
           "pivot_elev_mRL", "pivot_P_bar", "regime", "ns_verdict", "ns_confidence"]


def regime_note(wname, e):
    notes = []
    if e["regime"] == "conductive" and e["Tmax_elev_mRL"] == e["TD_mRL"]:
        notes.append("temperature still rising at TD; conductive above the deepest point")
    if e["regime"] == "upflow" and e["TD_mRL"] > -500:
        notes.append(f"shallow well (TD {fmt(e['TD_mRL'], 1)} mRL): upflow inferred from a hot, near-boiling, "
                     "isothermal base; the well does not reach the deep upflow")
    if e["regime"] == "recharge/downflow" and WD[wname]["kind"] == "injector":
        notes.append("injector, but feed zones colder than the injectate: cold profile is natural")
    if e["regime"] == "outflow" and np.isfinite(e["deep_gradient_degC_per_km"]):
        notes.append(f"temperature rises again slowly below the inversion base "
                     f"({fmt(e['deep_gradient_degC_per_km'], 1)} degC/km)")
    if e["regime"] == "outflow" and e["T_at_TD_degC"] >= 0.9 * float(WS.Tmax_degC.max()):
        notes.append(f"T at TD {fmt(e['T_at_TD_degC'], 1)} degC is within 10 % of the field Tmax: the rule gives "
                     "outflow because of the shallower inversion, but the deep reheating suggests the well lies "
                     "above the upflow; treated as an upflow well in the conceptual model")
    if e["near_boiling"]:
        notes.append(f"near boiling (min margin {fmt(e['min_boiling_margin_degC'], 2)} degC at "
                     f"{fmt(e['min_boiling_margin_elev_mRL'], 1)} mRL)")
    return "; ".join(notes)


WS["regime_note"] = [regime_note(w, WS.iloc[i]) for i, w in enumerate(WS.well)]
WS["ns_reason"] = [NS[w]["reason"] for w in WS.well]
WS = WS[WS_COLS + [c for c in WS.columns if c not in WS_COLS]]
write_csv(WS, os.path.join(CAL, "well_summary.csv"))
WSI = WS.set_index("well")

# ------------------------------------------------------------------------------------------
# 11. T and P at model layer centres (interpolation only)
# ------------------------------------------------------------------------------------------
LC_ROWS = []
for wname in WELLS:
    g = wsurvey(wname)
    zmax, zmin = g.elevation_mRL.max(), g.elevation_mRL.min()
    for i in sorted(LAYERS):
        zc = LAYERS[i]["centre"]
        if zmin - 1e-9 <= zc <= zmax + 1e-9:
            T, P = interp_w(wname, "T_degC", zc), interp_w(wname, "P_bar", zc)
            ts = r(tsat(P), 2)
            blk = block_at(wname, zc)
            LC_ROWS.append(dict(well=wname, layer=i, layer_centre_mRL=zc, block=blk, T_degC=r(T, 2), P_bar=r(P, 2),
                                Tsat_degC=ts, boiling_margin_degC=r(ts - T, 2) if np.isfinite(ts) else np.nan,
                                P_minus_trend_bar=r(P - trend(zc), 2)))
LC = pd.DataFrame(LC_ROWS)
write_csv(LC, os.path.join(CAL, "T_P_at_layer_centres.csv"))

# common elevations (every 250 m), interpolation only
COMMON_Z = list(np.arange(0.0, -3000.0, -250.0))
CE_ROWS = []
for z in COMMON_Z:
    for wname in WELLS:
        T, P = interp_w(wname, "T_degC", z), interp_w(wname, "P_bar", z)
        if np.isfinite(T):
            ts = tsat(P)
            CE_ROWS.append(dict(elevation_mRL=z, well=wname, T_degC=r(T, 1), P_bar=r(P, 1),
                                P_trend_bar=r(trend(z), 1), P_minus_trend_bar=r(P - trend(z), 1),
                                boiling_margin_degC=r(ts - T, 1) if np.isfinite(ts) else np.nan))
CE = pd.DataFrame(CE_ROWS)
write_csv(CE, os.path.join(ADD, "T_P_at_common_elevations.csv"))

# ------------------------------------------------------------------------------------------
# 12. Elevation bands (F3 data)
# ------------------------------------------------------------------------------------------
zall = S.elevation_mRL
b_top = math.ceil(zall.max() / BAND) * BAND
b_bot = math.floor(zall.min() / BAND) * BAND
cand_bands = [(t, t - BAND) for t in np.arange(b_top, b_bot, -BAND)]


def in_band(z, band):
    return band[1] <= z < band[0]


nwell_band = {bd: sum(any(in_band(z, bd) for z in wsurvey(w).elevation_mRL) for w in WELLS) for bd in cand_bands}
ok_b = [bd for bd in cand_bands if nwell_band[bd] >= 3]
BANDS = [bd for bd in cand_bands if ok_b[0][0] >= bd[0] >= ok_b[-1][0]]
EB_ROWS = []
for bd in BANDS:
    rows_b = []
    for wname in WELLS:
        g = wsurvey(wname)
        g = g[[in_band(z, bd) for z in g.elevation_mRL]]
        if g.empty:
            continue
        xm, ym = track_xy(wname, 0.5 * (bd[0] + bd[1]))
        rows_b.append(dict(well=wname, band_top_mRL=bd[0], band_bottom_mRL=bd[1], n_points=len(g),
                           T_mean_degC=g.T_degC.mean(), P_mean_bar=g.P_bar.mean(),
                           P_minus_trend_bar=float(np.mean(g.P_bar - trend(g.elevation_mRL))),
                           z_mean_mRL=g.elevation_mRL.mean(), x_mid=xm, y_mid=ym))
    med = float(np.median([x["T_mean_degC"] for x in rows_b]))
    for x in rows_b:
        x["T_field_median_degC"] = med
        x["T_minus_field_median_degC"] = x["T_mean_degC"] - med
    EB_ROWS += rows_b
EB = pd.DataFrame(EB_ROWS)
for c, d in (("T_mean_degC", 1), ("T_field_median_degC", 1), ("T_minus_field_median_degC", 1), ("P_mean_bar", 1),
             ("P_minus_trend_bar", 1), ("z_mean_mRL", 1), ("x_mid", 1), ("y_mid", 1)):
    EB[c] = [r(v, d) for v in EB[c]]
EB = EB[["well", "band_top_mRL", "band_bottom_mRL", "n_points", "T_mean_degC", "T_field_median_degC",
         "T_minus_field_median_degC", "P_mean_bar", "P_minus_trend_bar", "z_mean_mRL", "x_mid", "y_mid"]]
EB["band_top_mRL"] = EB.band_top_mRL.astype(int)
EB["band_bottom_mRL"] = EB.band_bottom_mRL.astype(int)
write_csv(EB, os.path.join(ADD, "T_P_at_elevation_bands.csv"))

# ------------------------------------------------------------------------------------------
# 13. Hydrostatic check (CLAUDE.md 4.3)
# ------------------------------------------------------------------------------------------


def dp_if97(zs, ts, ps, z0, z1, step=1.0):
    """Integrate dP/dz = -rho(T(z), P) g downward from z0 (P measured) to z1, T linear in z."""
    zz, tt = zs[::-1], ts[::-1]
    p = float(np.interp(z0, zz, ps[::-1]))
    p_start = p
    n = max(1, int(math.ceil((z0 - z1) / step)))
    h = (z0 - z1) / n
    z = z0
    for _ in range(n):
        t_a = float(np.interp(z, zz, tt))
        t_b = float(np.interp(z - h, zz, tt))
        k1 = rho_liq(t_a, p) * G * h / 1e5
        k2 = rho_liq(t_b, p + k1) * G * h / 1e5
        p += 0.5 * (k1 + k2)
        z -= h
    return p - p_start


HY_ROWS = []
for wname in WELLS:
    w = WD[wname]
    g = wsurvey(wname)
    z, t, p, bm = (g.elevation_mRL.to_numpy(), g.T_degC.to_numpy(), g.P_bar.to_numpy(),
                   g.boiling_margin_degC.to_numpy())
    elig = (p >= LIQ_P_MIN) & np.where(np.isfinite(bm), bm >= HYD_BM_MIN, p > P_CRIT)
    feats = sorted([w["shoe"]] + [f["mid"] for f in w["fz"]], reverse=True)
    # eligible intervals = runs of consecutive eligible points
    runs, i = [], 0
    while i < len(z):
        if elig[i]:
            j = i
            while j + 1 < len(z) and elig[j + 1]:
                j += 1
            if j > i:
                runs.append((i, j))
            i = j + 1
        else:
            i += 1
    for a, b in runs:
        wins, start = [], a
        for k in range(a + 1, b + 1):
            length = z[start] - z[k]
            reached = any((z[k] <= f < z[k - 1]) or abs(f - z[k]) < 1e-9 for f in feats)
            if length >= HYD_WIN - 1e-9 and reached:
                wins.append([start, k])
                start = k
        if start < b:
            if z[start] - z[b] < HYD_WIN and wins:
                wins[-1][1] = b
            else:
                wins.append([start, b])
        for s0, s1 in wins:
            z0, z1 = z[s0], z[s1]
            dpm = p[s1] - p[s0]
            dpi = dp_if97(z, t, p, z0, z1)
            L = z0 - z1
            band = "cased" if z1 >= w["shoe"] - 1e-9 else ("open" if z0 < w["shoe"] - 1e-9 else "")
            HY_ROWS.append(dict(well=wname, z_top_mRL=z0, z_bottom_mRL=z1, length_m=r(L, 1), band=band,
                                dP_measured_bar=r(dpm, 2), dP_IF97_bar=r(dpi, 2),
                                ratio_measured_to_IF97=r(dpm / dpi, 4), rho_implied_kgm3=r(dpm * 1e5 / (G * L), 1),
                                rho_IF97_mean_kgm3=r(dpi * 1e5 / (G * L), 1),
                                ratio_rounding_uncertainty=r(2 * 0.05 / dpi, 4),
                                rho_rounding_uncertainty_kgm3=r(0.1e5 / (G * L), 1),
                                T_mean_degC=r(float(np.mean(t[s0:s1 + 1])), 1), n_points=s1 - s0 + 1,
                                band_note="" if band else "window spans the production casing shoe (neither band)"))
HY = pd.DataFrame(HY_ROWS)
write_csv(HY, os.path.join(ADD, "hydrostatic_check.csv"))

COLD_WELLS = [w for w in WELLS if WSI.loc[w, "Tmax_degC"] < COLD_TMAX]
HS_ROWS = []
cold_med = {}
for band in ("cased", "open", "all"):
    sub = HY[HY.well.isin(COLD_WELLS)] if band == "all" else HY[HY.well.isin(COLD_WELLS) & (HY.band == band)]
    cold_med[band] = float(sub.ratio_measured_to_IF97.median()) if len(sub) else np.nan
for wname in WELLS + ["FIELD"]:
    for band in ("cased", "open", "all"):
        sub = HY if wname == "FIELD" else HY[HY.well == wname]
        if band != "all":
            sub = sub[sub.band == band]
        if sub.empty:
            continue
        med = float(sub.ratio_measured_to_IF97.median())
        cm = cold_med[band]
        HS_ROWS.append(dict(well=wname, band=band, median_ratio=r(med, 4),
                            min_ratio=r(sub.ratio_measured_to_IF97.min(), 4),
                            max_ratio=r(sub.ratio_measured_to_IF97.max(), 4),
                            cold_well_median_ratio=r(cm, 4),
                            excess_over_cold_component=r(med - cm, 4) if np.isfinite(cm) else np.nan,
                            cold_well=wname in COLD_WELLS, n_windows=len(sub)))
HS = pd.DataFrame(HS_ROWS)
write_csv(HS, os.path.join(ADD, "hydrostatic_summary.csv"))
log(f"hydrostatic: {len(HY)} windows; cold wells {COLD_WELLS}; cold-well median ratio {cold_med}")

# ------------------------------------------------------------------------------------------
# 14. Well geometry check
# ------------------------------------------------------------------------------------------


def dist_to_polyline(pt, segs):
    best = np.inf
    px, py = pt
    for s in segs:
        if len(s) == 1:
            best = min(best, math.hypot(px - s[0, 0], py - s[0, 1]))
            continue
        a, b = s[:-1], s[1:]
        d = b - a
        L2 = (d ** 2).sum(1)
        L2[L2 == 0] = 1e-12
        tt = np.clip(((px - a[:, 0]) * d[:, 0] + (py - a[:, 1]) * d[:, 1]) / L2, 0, 1)
        cx, cy = a[:, 0] + tt * d[:, 0], a[:, 1] + tt * d[:, 1]
        best = min(best, float(np.min(np.hypot(px - cx, py - cy))))
    return best


def dist_feature(pt, name):
    if name in FAULTS:
        return dist_to_polyline(pt, FAULTS[name])
    wb = WATER[name]
    if wb["fill"] and MPath(wb["segs"][0]).contains_point(pt):
        return 0.0
    return dist_to_polyline(pt, wb["segs"])


def safe(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


GEO_ROWS = []
for wname in WELLS:
    w = WD[wname]
    pt = tuple(w["track"][0, :2])
    dfa = {n: dist_feature(pt, n) for n in FAULTS}
    dwa = {n: dist_feature(pt, n) for n in WATER}
    nf = min(dfa, key=dfa.get)
    nw = min(dwa, key=dwa.get)
    grid_col = column_at(*pt)
    gw = [l for l in open(GEOM_PATH).read().split("WELLS", 1)[1].splitlines() if l.startswith(w["aut2"])]
    gtop = float(gw[0].split()[3]) if gw else np.nan
    row = dict(well=wname, wellhead_mRL=w["wellhead"], model_surface_mRL=w["model_surface"],
               wellhead_minus_surface_m=r(w["wellhead"] - w["model_surface"], 1),
               max_horizontal_deviation_m=r(w["max_dev"], 1), nearest_fault=nf, dist_nearest_fault_m=r(dfa[nf], 0),
               nearest_water_body=nw, dist_nearest_water_body_m=r(dwa[nw], 0))
    for n in list(FAULTS) + list(WATER):
        row[f"dist_{safe(n)}_m"] = r(dfa.get(n, dwa.get(n)), 0)
    row.update(wellhead_column=w["head_col"], wellhead_block=w["head_block"],
               point_in_polygon_column=grid_col or "", column_agrees=(grid_col == w["head_col"]),
               first_track_point_mRL=float(w["track"][0, 2]), geometry_file_wellhead_mRL=gtop,
               shoe_mRL=w["shoe"], TD_mRL=w["td"], first_survey_mRL=wsurvey(wname).elevation_mRL.max())
    GEO_ROWS.append(row)
GEO = pd.DataFrame(GEO_ROWS)
write_csv(GEO, os.path.join(ADD, "well_geometry_check.csv"))
GEOI = GEO.set_index("well")

# ------------------------------------------------------------------------------------------
# 15. Section lines (F4)
# ------------------------------------------------------------------------------------------
HOT2 = list(WS.sort_values(["Tmax_degC", "well"], ascending=[False, True]).well[:2])
_allv = np.array(list(VERTS.values()))
MODEL_BBOX = (_allv[:, 0].min(), _allv[:, 0].max(), _allv[:, 1].min(), _allv[:, 1].max())


def in_model(x, y):
    if not (MODEL_BBOX[0] <= x <= MODEL_BBOX[1] and MODEL_BBOX[2] <= y <= MODEL_BBOX[3]):
        return False
    return column_at(x, y) is not None


def extend_line(p, u, step=20.0):
    """Ends of the line through p with direction u, clipped to the model extent."""
    ends = []
    for sgn in (-1, 1):
        t = 0.0
        while in_model(p[0] + sgn * (t + step) * u[0], p[1] + sgn * (t + step) * u[1]):
            t += step
        ends.append(np.array([p[0] + sgn * t * u[0], p[1] + sgn * t * u[1]]))
    return ends


WH = {w: WD[w]["track"][0, :2] for w in WELLS}


def wells_near(a, u, L):
    out = {}
    for w in WELLS:
        d = WH[w] - a
        s, off = float(d @ u), float(d[0] * u[1] - d[1] * u[0])
        if abs(off) <= SECTION_HALFWIDTH and 0 <= s <= L:
            out[w] = (s, off)
    return out


def make_line(name, p, q):
    u = (q - p) / np.linalg.norm(q - p)
    if u[1] < 0 or (u[1] == 0 and u[0] < 0):   # A at the southern end
        u = -u
    a, b = extend_line(p, u)
    L = float(np.linalg.norm(b - a))
    return dict(name=name, a=a, b=b, u=u, L=L, wells=wells_near(a, u, L))


LINES = [make_line("A", WH[HOT2[0]], WH[HOT2[1]])]
LINES[0]["basis"] = f"through the two hottest wells ({HOT2[0]}, {HOT2[1]})"
# additional lines: lines through pairs of wells that capture the most wells within the half-width
cands = []
for i, w1 in enumerate(WELLS):
    for w2 in WELLS[i + 1:]:
        ln = make_line("?", WH[w1], WH[w2])
        ang = math.degrees(math.atan2(ln["u"][1], ln["u"][0])) % 180
        hot = sum(WSI.loc[w, "Tmax_degC"] for w in ln["wells"])
        cands.append((len(ln["wells"]), hot, w1, w2, ang, ln))
cands.sort(key=lambda c: (-c[0], -c[1], c[2], c[3]))
used_ang = [math.degrees(math.atan2(LINES[0]["u"][1], LINES[0]["u"][0])) % 180]
for n, hot, w1, w2, ang, ln in cands:
    if len(LINES) >= 3:
        break
    if all(min(abs(ang - a), 180 - abs(ang - a)) > 25 for a in used_ang) and n >= 3:
        ln["name"] = "ABC"[len(LINES)]
        ln["basis"] = f"through {w1} and {w2}; captures {n} wells within {SECTION_HALFWIDTH:.0f} m"
        LINES.append(ln)
        used_ang.append(ang)
SL_ROWS = []
for ln in LINES:
    SL_ROWS.append(dict(line=ln["name"], start_label=ln["name"], end_label=ln["name"] + "'",
                        x_start=r(ln["a"][0], 0), y_start=r(ln["a"][1], 0), x_end=r(ln["b"][0], 0),
                        y_end=r(ln["b"][1], 0), length_m=r(ln["L"], 0), basis=ln["basis"],
                        wells="; ".join(f"{w} (s={r(s, 0):.0f} m, offset {r(o, 0):+.0f} m)"
                                        for w, (s, o) in sorted(ln["wells"].items(), key=lambda kv: kv[1][0]))))
SLDF = pd.DataFrame(SL_ROWS)
write_csv(SLDF, os.path.join(ADD, "section_lines.csv"))

# ------------------------------------------------------------------------------------------
# 16. Figure helpers (house style, CLAUDE.md 8.1)
# ------------------------------------------------------------------------------------------
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import BoundaryNorm  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Polygon as MPoly, FancyArrowPatch  # noqa: E402
import matplotlib.tri as mtri  # noqa: E402

plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "savefig.dpi": 150, "figure.dpi": 150})
FIELD_TMAX = float(WS.Tmax_degC.max())
T_TOP = math.ceil(FIELD_TMAX / 50.0) * 50.0
TCMAP = plt.get_cmap("plasma")
T_LEVELS = np.arange(0, T_TOP + 10, 10)
TNORM = BoundaryNorm(T_LEVELS, TCMAP.N)
DCMAP = plt.get_cmap("RdBu_r")
Z_TOP = math.ceil(max(WD[w]["wellhead"] for w in WELLS) / 100.0) * 100.0
Z_BOT = math.floor(S.elevation_mRL.min() / 100.0) * 100.0
ISO_LEVELS = list(np.arange(50, T_TOP + 1, 50))
_tr = np.vstack([WD[w]["track"][:, :2] for w in WELLS])
WELL_EXTENT = (_tr[:, 0].min() - 1300, _tr[:, 0].max() + 1300, _tr[:, 1].min() - 1300, _tr[:, 1].max() + 1300)
FULL_EXTENT = (MODEL_BBOX[0] - 300, MODEL_BBOX[1] + 300, MODEL_BBOX[2] - 300, MODEL_BBOX[3] + 300)
FIGS_CAL, FIGS_ADD = [], []
DPI = 150
P_DEV_MAX = math.ceil(np.nanmax(np.abs(TGT.P_minus_trend_bar)))
P_LEVELS = np.arange(-P_DEV_MAX, P_DEV_MAX + 1, 1.0)
PNORM = BoundaryNorm(P_LEVELS, DCMAP.N)


def savefig(fig, folder, name):
    path = os.path.join(folder, name)
    fig.savefig(path, dpi=DPI, metadata={"Software": None})
    plt.close(fig)
    (FIGS_CAL if folder == CALFIG else FIGS_ADD).append(rel(path))
    log(f"figure {rel(path)}")
    return rel(path)


_EDGES_ALL = sorted(_edges)
_SEG_ALL = np.array([[VERTS[a], VERTS[b]] for a, b in _EDGES_ALL])
_SEG_BND = np.array([[VERTS[a], VERTS[b]] for a, b in sorted(BOUNDARY_EDGES)])


def label_point(segs, extent):
    pts = np.vstack(segs)
    inside = pts[(pts[:, 0] > extent[0]) & (pts[:, 0] < extent[1]) & (pts[:, 1] > extent[2]) & (pts[:, 1] < extent[3])]
    if not len(inside):
        return None
    c = np.array([(extent[0] + extent[1]) / 2, (extent[2] + extent[3]) / 2])
    # a point inside the view, 70 % of the way from the view centre to the far end of the feature
    d = np.hypot(*(inside - c).T)
    return inside[np.argsort(d)[int(0.7 * (len(d) - 1))]]


def fault_label_point(segs, extent):
    pts = []
    for s_ in segs:
        for a_, b_ in zip(s_[:-1], s_[1:]):
            for t_ in np.linspace(0, 1, 60):
                pts.append(a_ + t_ * (b_ - a_))
    pts = np.array(pts)
    mx, my = 0.06 * (extent[1] - extent[0]), 0.06 * (extent[3] - extent[2])
    ok = (pts[:, 0] > extent[0] + mx) & (pts[:, 0] < extent[1] - mx) & (pts[:, 1] > extent[2] + my) & \
         (pts[:, 1] < extent[3] - my)
    pts = pts[ok]
    if not len(pts):
        return None
    wh = np.array([WH[w] for w in WELLS])
    d = np.min(np.hypot(pts[:, None, 0] - wh[None, :, 0], pts[:, None, 1] - wh[None, :, 1]), axis=1)
    return pts[int(np.argmax(np.minimum(d, 0.25 * (extent[1] - extent[0]))))] if len(pts) else None


def well_label_offsets():
    off = {}
    for w in WELLS:
        others = [WH[o] for o in WELLS if o != w]
        d = np.array([np.linalg.norm(WH[w] - o) for o in others])
        v = WH[w] - others[int(np.argmin(d))]
        v = v / (np.linalg.norm(v) or 1)
        off[w] = v
    return off


LABEL_OFF = well_label_offsets()


_PLACED = {}
_DIRS = [(1, 0), (1, 1), (0, 1), (-1, 1), (-1, 0), (-1, -1), (0, -1), (1, -1)]


def _choose_dir(ax, w, xy, fs, dist):
    """Greedy label placement: the direction whose label box overlaps least with other wellheads and labels."""
    ax.apply_aspect()
    tr = ax.transData
    px, py = tr.transform(xy)
    k = fig_dpi = ax.figure.dpi / 72.0
    wd, ht = len(w) * fs * 0.62 * k, fs * 1.25 * k
    pts = [tr.transform(WH[o]) for o in WELLS if o != w]
    boxes = _PLACED.setdefault(id(ax), [])
    best = None
    for n, (dx, dy) in enumerate(_DIRS):
        ox, oy = px + dx * dist * fig_dpi, py + dy * dist * fig_dpi
        x0 = ox if dx > 0 else (ox - wd if dx < 0 else ox - wd / 2)
        y0 = oy if dy > 0 else (oy - ht if dy < 0 else oy - ht / 2)
        bx = (x0 - 2, y0 - 2, x0 + wd + 2, y0 + ht + 2)
        sc = 0.0
        for qx, qy in pts:
            if bx[0] - 4 <= qx <= bx[2] + 4 and bx[1] - 4 <= qy <= bx[3] + 4:
                sc += 1000
        for b in boxes:
            ox_ = max(0, min(bx[2], b[2]) - max(bx[0], b[0]))
            oy_ = max(0, min(bx[3], b[3]) - max(bx[1], b[1]))
            sc += ox_ * oy_
        sc += n * 0.01
        if best is None or sc < best[0]:
            best = (sc, dx, dy, bx)
    boxes.append(best[3])
    return best[1], best[2]


def label_well(ax, w, xy=None, color="k", fs=7, dist=7):
    x, y = WH[w] if xy is None else xy
    dx, dy = _choose_dir(ax, w, (x, y), fs, dist)
    v = np.array([dx, dy], float)
    ha = "left" if v[0] > 0 else ("right" if v[0] < 0 else "center")
    va = "bottom" if v[1] > 0 else ("top" if v[1] < 0 else "center")
    ax.annotate(w, (x, y), xytext=(v[0] * dist, v[1] * dist), textcoords="offset points", ha=ha, va=va,
                fontsize=fs, color=color, fontweight="bold", zorder=12,
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7))


def scale_bar(ax, extent, length=None):
    w = extent[1] - extent[0]
    if length is None:
        length = 1000.0 if w < 12000 else 2000.0
    x0 = extent[0] + 0.05 * w
    y0 = extent[2] + 0.05 * (extent[3] - extent[2])
    ax.plot([x0, x0 + length], [y0, y0], color="k", lw=3, solid_capstyle="butt", zorder=13)
    ax.text(x0 + length / 2, y0 + 0.015 * (extent[3] - extent[2]), f"{length / 1000:.0f} km", ha="center",
            va="bottom", fontsize=7, zorder=13)


def base_layer(ax, extent, label_wells=True, highlight=None, fs=7, fault_labels=True, small=False):
    ax.add_collection(LineCollection(_SEG_ALL, colors="0.85", linewidths=0.3, zorder=1))
    ax.add_collection(LineCollection(_SEG_BND, colors="0.35", linewidths=0.9, linestyles="--", zorder=2))
    for name, wb in WATER.items():
        for s in wb["segs"]:
            if wb["fill"]:
                ax.add_patch(MPoly(s, closed=True, fc="#cfe8fb", ec="#2b7bba", lw=0.6, zorder=2))
            else:
                ax.plot(s[:, 0], s[:, 1], color="#2b7bba", lw=0.9, zorder=3)
        p = label_point(wb["segs"], extent)
        if p is not None and not small:
            ax.text(p[0], p[1], name, color="#1f5f99", fontsize=fs - 1, style="italic", zorder=4,
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.6))
    for name, segs in FAULTS.items():
        for s in segs:
            ax.plot(s[:, 0], s[:, 1], color="k", lw=1.4, zorder=5)
        p = fault_label_point(segs, extent)
        if p is not None and fault_labels:
            ax.text(p[0], p[1], name, color="k", fontsize=fs - 1, fontweight="bold", zorder=6, ha="center",
                    bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.7))
    for w in WELLS:
        t = WD[w]["track"]
        col = WELL_COLOR[w] if highlight == w else "0.5"
        ax.plot(t[:, 0], t[:, 1], color=col, lw=1.2 if highlight == w else 0.8, zorder=7)
        ms = 5 if highlight == w else 3
        ax.plot(t[0, 0], t[0, 1], "o", color=WELL_COLOR[w] if highlight == w else "k", ms=ms, zorder=8)
        if label_wells:
            label_well(ax, w, fs=fs, color=WELL_COLOR[w] if highlight == w else "k")
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_aspect("equal")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
    ax.yaxis.set_major_locator(matplotlib.ticker.MaxNLocator(5))
    if not small:
        ax.set_xlabel("Easting (m)")
        ax.set_ylabel("Northing (m)")
        scale_bar(ax, extent)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
        scale_bar(ax, extent)


BASE_HANDLES = [Line2D([], [], color="0.85", lw=0.8, label="model grid columns"),
                Line2D([], [], color="0.35", lw=0.9, ls="--", label="model extent"),
                Line2D([], [], color="k", lw=1.4, label="fault"),
                Patch(fc="#cfe8fb", ec="#2b7bba", label="lake"),
                Line2D([], [], color="#2b7bba", lw=0.9, label="river"),
                Line2D([], [], marker="o", color="k", ls="", ms=3, label="wellhead"),
                Line2D([], [], color="0.5", lw=0.8, label="well track")]


def clip_to_view(a, b, ext, inset=0.04):
    """End points of segment a-b clipped to the (slightly inset) view box."""
    x0, x1, y0, y1 = ext
    mx, my = inset * (x1 - x0), inset * (y1 - y0)
    x0, x1, y0, y1 = x0 + mx, x1 - mx, y0 + my, y1 - my
    ts = np.linspace(0, 1, 4001)
    pts = a[None, :] + ts[:, None] * (b - a)[None, :]
    ok = (pts[:, 0] >= x0) & (pts[:, 0] <= x1) & (pts[:, 1] >= y0) & (pts[:, 1] <= y1)
    if not ok.any():
        return None
    return pts[np.argmax(ok)], pts[len(ok) - 1 - np.argmax(ok[::-1])]


def draw_lines(ax, lines=None, fs=8):
    ext = (*ax.get_xlim(), *ax.get_ylim())
    for ln in (lines or LINES):
        ax.plot([ln["a"][0], ln["b"][0]], [ln["a"][1], ln["b"][1]], color="#7a1fa2", lw=1.3, zorder=9)
        cl = clip_to_view(ln["a"], ln["b"], ext)
        if cl is None:
            continue
        for p, lab in ((cl[0], ln["name"]), (cl[1], ln["name"] + "'")):
            ax.text(p[0], p[1], lab, fontsize=fs + 1, fontweight="bold", color="#7a1fa2", zorder=14, ha="center",
                    va="center", bbox=dict(boxstyle="circle,pad=0.15", fc="white", ec="#7a1fa2", lw=0.8))


# ------------------------------------------------------------------------------------------
# 17. Section data and drawing (F4, F6, companions)
# ------------------------------------------------------------------------------------------
def proj(ln, x, y):
    d = np.array([x, y]) - ln["a"]
    return float(d @ ln["u"]), float(d[0] * ln["u"][1] - d[1] * ln["u"][0])


def section_wells(ln):
    out = []
    for w in sorted(ln["wells"], key=lambda k: ln["wells"][k][0]):
        g = wsurvey(w)
        s = np.array([proj(ln, x, y)[0] for x, y in zip(g.x, g.y)])
        ts = np.array([proj(ln, x, y)[0] for x, y in WD[w]["track"][:, :2]])
        out.append(dict(well=w, s=s, z=g.elevation_mRL.to_numpy(), T=g.T_degC.to_numpy(),
                        dP=(g.P_bar - trend(g.elevation_mRL)).to_numpy(), bm=g.boiling_margin_degC.to_numpy(),
                        off=ln["wells"][w][1], s_head=ln["wells"][w][0], track_s=ts, track_z=WD[w]["track"][:, 2]))
    return out


def surface_profile(ln, step=50.0):
    ss = np.arange(0, ln["L"] + step, step)
    zs = []
    for s in ss:
        p = ln["a"] + s * ln["u"]
        c = column_at(p[0], p[1])
        zs.append(surface_of(c) if c else np.nan)
    return ss, np.array(zs)


def envelope_meshes(sw, var, nu=30, dz=5.0):
    """Quad meshes between neighbouring wells over their common elevation range (no extrapolation)."""
    meshes = []
    brk = sorted({float(v) for d in sw for v in (d["z"].max(), d["z"].min())}, reverse=True)
    pairs = []
    for zt_, zb_ in zip(brk, brk[1:]):
        zm = 0.5 * (zt_ + zb_)
        have = [d for d in sw if d["z"].min() <= zm <= d["z"].max()]
        have.sort(key=lambda d: float(np.interp(zm, d["z"][::-1], d["s"][::-1])))
        for a, b in zip(have, have[1:]):
            pairs.append((a, b, zt_, zb_))
    for a, b, top, bot in pairs:
        if top - bot < 1e-6:
            continue
        zz = np.linspace(top, bot, max(2, int(math.ceil((top - bot) / dz)) + 1))
        sa = np.interp(zz, a["z"][::-1], a["s"][::-1])
        sb = np.interp(zz, b["z"][::-1], b["s"][::-1])
        va = np.interp(zz, a["z"][::-1], a[var][::-1])
        vb = np.interp(zz, b["z"][::-1], b[var][::-1])
        u = np.linspace(0, 1, nu)
        SS = sa[:, None] + u[None, :] * (sb - sa)[:, None]
        ZZ = np.repeat(zz[:, None], nu, axis=1)
        VV = va[:, None] + u[None, :] * (vb - va)[:, None]
        outline = np.vstack([np.c_[sa, zz], np.c_[sb[::-1], zz[::-1]], [[sa[0], zz[0]]]])
        meshes.append(dict(S=SS, Z=ZZ, V=VV, outline=outline, wells=(a["well"], b["well"])))
    return meshes


def hatch_nodata(ax, ln, ss, zs):
    ok = np.isfinite(zs)
    if not ok.any():
        return
    poly = np.vstack([np.c_[ss[ok], zs[ok]], [[ss[ok][-1], Z_BOT - 500], [ss[ok][0], Z_BOT - 500]]])
    ax.add_patch(MPoly(poly, closed=True, fc="0.93", ec="none", hatch="///", zorder=0))
    ax.plot(ss[ok], zs[ok], color="saddlebrown", lw=1.2, zorder=6)


def surface_ticks(ax, ln):
    feats = []
    a, b = ln["a"], ln["b"]
    for name, segs in list(FAULTS.items()) + [(n, wb["segs"]) for n, wb in WATER.items()]:
        for s in segs:
            for p, q in zip(s[:-1], s[1:]):
                d1, d2 = b - a, q - p
                den = d1[0] * d2[1] - d1[1] * d2[0]
                if abs(den) < 1e-9:
                    continue
                t = ((p[0] - a[0]) * d2[1] - (p[1] - a[1]) * d2[0]) / den
                v = ((p[0] - a[0]) * d1[1] - (p[1] - a[1]) * d1[0]) / den
                if 0 <= t <= 1 and 0 <= v <= 1:
                    feats.append((t * ln["L"], name))
    feats.sort()
    kept = []
    for s, n in feats:   # one tick per feature per 150 m
        if not any(n == kn and abs(s - ks) < 150 for ks, kn in kept):
            kept.append((s, n))
    x0_, x1_ = section_xlim(ln)
    kept = [(s_, n_) for s_, n_ in kept if x0_ <= s_ <= x1_]
    for s, n in kept:
        col = "k" if n in FAULTS else "#2b7bba"
        ax.annotate("", xy=(s, Z_TOP), xytext=(s, Z_TOP + 0.03 * (Z_TOP - Z_BOT)),
                    arrowprops=dict(arrowstyle="-", color=col, lw=2), annotation_clip=False)
        ax.text(s, Z_TOP + 0.035 * (Z_TOP - Z_BOT), n, rotation=45, fontsize=6.5, color=col, ha="left",
                va="bottom", clip_on=False)
    return kept


def draw_wells_on_section(ax, ln, sw, var, show_points=True, norm=None, cmap=None):
    for d in sw:
        w = d["well"]
        ax.plot(d["track_s"], d["track_z"], color="0.25", lw=1.0, zorder=7)
        sh = WD[w]["shoe"]
        s_sh = np.interp(sh, d["track_z"][::-1], d["track_s"][::-1])
        ax.plot([s_sh - 70, s_sh + 70], [sh, sh], color="k", lw=2, zorder=8)
        for f in WD[w]["fz"]:
            s_f = np.interp(f["mid"], d["track_z"][::-1], d["track_s"][::-1])
            ax.add_patch(plt.Rectangle((s_f - 45, f["bottom"]), 90, f["top"] - f["bottom"], fc="none", ec="lime",
                                       lw=1.6, zorder=8))
        # near-boiling stipple
        nb = np.isfinite(d["bm"]) & (d["bm"] < NEAR_BOIL)
        i = 0
        while i < len(nb):
            if nb[i]:
                j = i
                while j + 1 < len(nb) and nb[j + 1]:
                    j += 1
                z0 = d["z"][i] + (10 if i == j else 0)
                z1 = d["z"][j] - (10 if i == j else 0)
                ax.add_patch(plt.Rectangle((d["s"][i] - 110, z1), 220, z0 - z1, fc="none", ec="k", hatch="...",
                                           lw=0.4, zorder=9))
                i = j + 1
            else:
                i += 1
        if show_points:
            ax.scatter(d["s"], d["z"], c=d[var], cmap=cmap, norm=norm, s=16, edgecolors="k", linewidths=0.4,
                       zorder=10)
        ax.text(d["track_s"][-1] + 60, d["track_z"][-1] + 30, f"{w}\noffset {d['off']:+.0f} m", ha="left",
                va="bottom",
                fontsize=6.5, zorder=11, color=WELL_COLOR[w], fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.1", fc="white", ec="none", alpha=0.8))


def vertical_exaggeration(ax, fig):
    bb = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    xr = np.diff(ax.get_xlim())[0] / bb.width
    yr = np.diff(ax.get_ylim())[0] / bb.height
    return xr / yr


SECTION_MARGIN = 2000.0


def section_xlim(ln):
    ss_ = [v for d in section_wells(ln) for v in d["track_s"]]
    return max(0.0, min(ss_) - SECTION_MARGIN), min(ln["L"], max(ss_) + SECTION_MARGIN)


def section_axes_setup(ax, ln):
    x0, x1 = section_xlim(ln)
    ax.set_xlim(x0, x1)
    ax.set_ylim(Z_BOT, Z_TOP)
    ax.set_xlabel(f"Distance along line {ln['name']}-{ln['name']}' from {ln['name']} (m); "
                  f"view cropped to the wells +/- {SECTION_MARGIN:.0f} m (full line {ln['L']:.0f} m)")
    ax.set_ylabel("Elevation (mRL)")
    for p, lab, ha in ((x0, ln["name"], "left"), (x1, ln["name"] + "'", "right")):
        dist_end = p if lab == ln["name"] else ln["L"] - p
        txt = (f"{lab}" if dist_end < 1 else (f"<- {lab} ({dist_end / 1000:.1f} km)" if ha == "left"
                                               else f"{lab} ({dist_end / 1000:.1f} km) ->"))
        ax.text(p, Z_BOT + 0.02 * (Z_TOP - Z_BOT), txt, fontsize=10, fontweight="bold", color="#7a1fa2",
                ha=ha, va="bottom", zorder=15,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#7a1fa2", lw=0.8))


def draw_section(ax, fig, ln, var="T", for_concept=False):
    sw = section_wells(ln)
    ss, zs = surface_profile(ln)
    hatch_nodata(ax, ln, ss, zs)
    meshes = envelope_meshes(sw, "T" if var == "T" else "dP")
    if var == "T":
        cmap, norm, levels = TCMAP, TNORM, T_LEVELS
    else:
        cmap, norm, levels = DCMAP, PNORM, P_LEVELS
    for m in meshes:
        ax.contourf(m["S"], m["Z"], m["V"], levels=levels, cmap=cmap, norm=norm, zorder=3, extend="both")
        if var == "T":
            cs = ax.contour(m["S"], m["Z"], m["V"], levels=ISO_LEVELS, colors="k",
                            linewidths=[1.6 if lv in (100, 200) else 0.6 for lv in ISO_LEVELS], zorder=4)
            ax.clabel(cs, fmt="%d", fontsize=6, inline=True)
        else:
            ax.contour(m["S"], m["Z"], m["V"], levels=[l for l in np.arange(-P_DEV_MAX, P_DEV_MAX + 1, 5) if l != 0],
                       colors="k", linewidths=0.5, zorder=4)
            ax.contour(m["S"], m["Z"], m["V"], levels=[0], colors="k", linewidths=1.4, zorder=4)
        ax.plot(m["outline"][:, 0], m["outline"][:, 1], color="k", ls=":", lw=0.8, zorder=5)
    draw_wells_on_section(ax, ln, sw, "T" if var == "T" else "dP", norm=norm, cmap=cmap)
    section_axes_setup(ax, ln)
    kept = surface_ticks(ax, ln)
    return sw, meshes, kept


def section_legend_handles(var="T"):
    h = [Patch(fc="0.93", ec="0.6", hatch="///", label="no data"),
         Line2D([], [], color="saddlebrown", lw=1.2, label="model surface"),
         Line2D([], [], color="k", ls=":", lw=0.8, label="data envelope"),
         Line2D([], [], color="0.25", lw=1.0, label="well track (projected)"),
         Line2D([], [], color="k", lw=2, label="production casing shoe"),
         Patch(fc="none", ec="lime", lw=1.6, label="feed zone"),
         Patch(fc="none", ec="k", hatch="...", lw=0.4, label="boiling margin < 5 degC"),
         Line2D([], [], marker="o", color="w", mec="k", ms=5, ls="", label="survey point")]
    if var == "T":
        h += [Line2D([], [], color="k", lw=1.6, label="100 / 200 degC isotherm"),
              Line2D([], [], color="k", lw=0.6, label="50 degC isotherms")]
    else:
        h += [Line2D([], [], color="k", lw=1.4, label="0 bar (on trend)"),
              Line2D([], [], color="k", lw=0.5, label="5 bar contours")]
    return h


# ------------------------------------------------------------------------------------------
# 18. Standard figures
# ------------------------------------------------------------------------------------------
FIG_INDEX = []


def index_fig(fid, path, shows, wells, line, sources, claim):
    FIG_INDEX.append(dict(figure_id=fid, file=path, shows=shows, wells=wells, elevation_or_line=line,
                          data_sources=sources, claim_type=claim))


SRC_WELL = "well_data_json"
SRC_BASE = "mapping; geometry_dat; well_data_json (tracks)"

# F1 base map
fig, axs = plt.subplots(1, 2, figsize=(15, 7.2), gridspec_kw=dict(width_ratios=[1.05, 1]))
base_layer(axs[0], FULL_EXTENT, fs=6)
draw_lines(axs[0])
axs[0].plot([WELL_EXTENT[0], WELL_EXTENT[1], WELL_EXTENT[1], WELL_EXTENT[0], WELL_EXTENT[0]],
            [WELL_EXTENT[2], WELL_EXTENT[2], WELL_EXTENT[3], WELL_EXTENT[3], WELL_EXTENT[2]], color="r", lw=0.8)
axs[0].set_title("Model extent with faults, water bodies, wells and section lines (plan view, m)")
base_layer(axs[1], WELL_EXTENT, fs=7.5)
draw_lines(axs[1])
axs[1].set_title("Well field (red box at left): wellheads, tracks and section lines (plan view, m)")
fig.legend(handles=BASE_HANDLES + [Line2D([], [], color="#7a1fa2", lw=1.3, label="section line (F4, F6)")],
           loc="lower center", ncol=8, fontsize=7, frameon=False)
fig.suptitle("F1 Base map: model grid, mapped faults and water bodies, wells and section lines", fontsize=10)
fig.tight_layout(rect=(0, 0.05, 1, 0.96))
p = savefig(fig, CALFIG, "well_map.png")
index_fig("F1", p, "base map: grid, extent, faults, water bodies, wellheads, tracks, section lines",
          "all", "plan view; lines " + ", ".join(l["name"] for l in LINES), SRC_BASE, "obs")


# F2 well profiles
def locator(ax, highlight=None, group=None):
    base_layer(ax, WELL_EXTENT, label_wells=False, highlight=highlight, small=True, fault_labels=True, fs=6)
    for w in WELLS:
        if group and w in group:
            ax.plot(*WH[w], "o", color=WELL_COLOR[w], ms=6, zorder=10)
            label_well(ax, w, color=WELL_COLOR[w], fs=6, dist=5)
        elif highlight == w:
            label_well(ax, w, color=WELL_COLOR[w], fs=7, dist=5)


def profile_T(ax, w, small=False):
    g = wsurvey(w)
    wd = WD[w]
    col = WELL_COLOR[w]
    for f in wd["fz"]:
        ax.axhspan(f["bottom"], f["top"], color="lime", alpha=0.35, zorder=1)
        if not small:
            ax.text(T_TOP * 0.99, f["mid"], f"feed zone ({f['quality']})", fontsize=6.5, ha="right", va="center",
                    zorder=6)
    nb = (g.boiling_margin_degC < NEAR_BOIL).to_numpy()
    for i in range(len(g)):
        if nb[i]:
            z0 = g.elevation_mRL[i - 1] if i > 0 and nb[i - 1] else g.elevation_mRL[i] + 5
            z1 = g.elevation_mRL[i + 1] if i + 1 < len(g) and nb[i + 1] else g.elevation_mRL[i] - 5
            ax.axhspan(min(z0, z1), max(z0, z1), color="orange", alpha=0.35, zorder=1, lw=0)
    ax.axhline(wd["shoe"], color="k", lw=1.0, ls="-", zorder=2)
    ax.axhline(wd["model_surface"], color="saddlebrown", lw=1.0, ls="--", zorder=2)
    ax.axhline(wd["wellhead"], color="0.3", lw=1.0, ls="-.", zorder=2)
    if not small:
        ax.text(2, wd["shoe"], "shoe", fontsize=7, va="bottom")
        ax.text(2, wd["model_surface"], f"model surface {fmt(wd['model_surface'], 0)} mRL", fontsize=6.5,
                va="bottom", color="saddlebrown")
        ax.text(T_TOP * 0.98, wd["wellhead"], f"wellhead {fmt(wd['wellhead'], 0)} mRL", fontsize=6.5, va="top",
                ha="right", color="0.3")
    ok = g.P_bar <= P_CRIT
    ax.plot(g.Tsat_degC[ok], g.elevation_mRL[ok], color="tab:red", lw=1.0, ls="--", zorder=3)
    ax.plot(g.T_degC, g.elevation_mRL, "o", color=col, ms=3.5 if small else 4.5, mec="k", mew=0.3, zorder=5)
    for iso_t, z_iso in ((100, WSI.loc[w, "z_100C_mRL"]), (200, WSI.loc[w, "z_200C_mRL"])):
        if np.isfinite(z_iso):
            ax.plot([iso_t], [z_iso], marker="D", color="k", ms=4, mfc="w", zorder=6)
            if not small:
                ax.text(iso_t + 4, z_iso, f"{iso_t} degC at {fmt(z_iso, 1)} mRL", fontsize=6.5, va="center")
    ax.set_xlim(0, T_TOP)
    ax.set_ylim(Z_BOT, Z_TOP)
    ax.grid(alpha=0.3, lw=0.4)


T_PROFILE_HANDLES = [Line2D([], [], marker="o", color="0.5", mec="k", ls="", label="survey temperature"),
                     Line2D([], [], color="tab:red", ls="--", label="Tsat(P), IF97 (<= 220.64 bar)"),
                     Line2D([], [], color="k", lw=1, label="production casing shoe"),
                     Line2D([], [], color="saddlebrown", ls="--", label="model surface"),
                     Line2D([], [], color="0.3", ls="-.", label="wellhead"),
                     Patch(color="lime", alpha=0.35, label="feed zone"),
                     Patch(color="orange", alpha=0.35, label="boiling margin < 5 degC"),
                     Line2D([], [], marker="D", color="k", mfc="w", ls="", label="100 / 200 degC crossing")]

for w in WELLS:
    wd = WD[w]
    g = wsurvey(w)
    fig = plt.figure(figsize=(14, 8))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.0, 1.25, 1.25], height_ratios=[1, 1.15])
    axm = fig.add_subplot(gs[0, 0])
    locator(axm, highlight=w)
    axm.set_title("Location (plan view)", fontsize=8)
    axl = fig.add_subplot(gs[1, 0])
    axl.axis("off")
    axt = fig.add_subplot(gs[:, 1])
    axp = fig.add_subplot(gs[:, 2], sharey=axt)
    profile_T(axt, w)
    axt.set_xlabel("Temperature (degC)")
    axt.set_ylabel("Elevation (mRL)")
    axt.set_title("Temperature")
    # pressure panel
    for f in wd["fz"]:
        axp.axhspan(f["bottom"], f["top"], color="lime", alpha=0.35, zorder=1)
    axp.axhline(wd["shoe"], color="k", lw=1.0)
    axp.axhline(wd["model_surface"], color="saddlebrown", lw=1.0, ls="--")
    axp.axhline(wd["wellhead"], color="0.3", lw=1.0, ls="-.")
    zt = np.array([Z_BOT, Z_TOP])
    axp.plot(trend(zt), zt, color="0.4", lw=0.8, ls=":", zorder=2)
    axp.plot(g.P_bar, g.elevation_mRL, "o", color=WELL_COLOR[w], ms=4.5, mec="k", mew=0.3, zorder=5)
    fzs = FZ[FZ.well == w]
    axp.plot(fzs.P_bar, fzs.mid_mRL, "s", color=WELL_COLOR[w], ms=9, mec="k", mew=0.8, zorder=6)
    pv_z, pv_p = WSI.loc[w, "pivot_elev_mRL"], WSI.loc[w, "pivot_P_bar"]
    axp.plot([pv_p], [pv_z], "*", color="gold", ms=16, mec="k", zorder=7)
    axp.text(pv_p + 6, pv_z, f"pivot {fmt(pv_p, 1)} bar at {fmt(pv_z, 1)} mRL", fontsize=7, va="center")
    wlv = WSI.loc[w, "water_level_mRL"]
    if np.isfinite(wlv):
        axp.axhline(wlv, color="tab:blue", lw=1.0, ls="--")
        axp.text(2, wlv, f"water level {fmt(wlv, 1)} mRL", fontsize=7, va="bottom", color="tab:blue")
    axp.set_xlim(0, math.ceil(S.P_bar.max() / 50) * 50)
    axp.set_xlabel("Pressure (bar)")
    axp.set_title("Pressure")
    axp.grid(alpha=0.3, lw=0.4)
    plt.setp(axp.get_yticklabels(), visible=False)
    hp = [Line2D([], [], marker="o", color=WELL_COLOR[w], mec="k", ls="", label="survey pressure"),
          Line2D([], [], marker="s", color=WELL_COLOR[w], mec="k", ms=8, ls="", label="feed-zone pressure"),
          Line2D([], [], marker="*", color="gold", mec="k", ms=12, ls="", label="pivot (dominant feed zone)"),
          Line2D([], [], color="tab:blue", ls="--", label="water level (fit, 3 shallowest liquid points)"),
          Line2D([], [], color="0.4", ls=":", label="trend of gauge-verified feed zones")]
    info = [Line2D([], [], ls="", label=f"survey date (label): {iso(wd['T_date'])}"),
            Line2D([], [], ls="", label=f"natural state: {NS[w]['verdict']} ({NS[w]['confidence']} confidence)"),
            Line2D([], [], ls="", label=f"regime (rule-based): {WSI.loc[w, 'regime']}")]
    axl.legend(handles=info + T_PROFILE_HANDLES + hp, loc="upper left", fontsize=7, frameon=False)
    fig.suptitle(f"F2 {w}: downhole temperature and pressure against elevation (mRL); degC, bar", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = savefig(fig, CALFIG, f"profile_{w}.png")
    index_fig("F2", p, f"{w} T and P profiles with Tsat, shoe, feed zones, surface, pivot, water level", w,
              f"{fmt(Z_TOP, 0)} to {fmt(Z_BOT, 0)} mRL", SRC_WELL + "; geometry_dat (surface)", "obs")

for gname, gw in GROUPS.items():
    n = len(gw)
    ncol = min(5, n)
    nrow = int(math.ceil(n / ncol))
    fig = plt.figure(figsize=(3.0 * ncol + 4, 4.2 * nrow + 0.6))
    gs = fig.add_gridspec(nrow, ncol + 1, width_ratios=[1.6] + [1] * ncol)
    axm = fig.add_subplot(gs[:, 0])
    locator(axm, group=gw)
    axm.set_title(f"Group '{gname}' (plan view)", fontsize=8)
    first = None
    for i, w in enumerate(gw):
        ax = fig.add_subplot(gs[i // ncol, 1 + i % ncol], sharex=first, sharey=first)
        first = first or ax
        profile_T(ax, w, small=True)
        ax.set_title(f"{w} ({WSI.loc[w, 'regime']})", fontsize=8, color=WELL_COLOR[w], fontweight="bold")
        if i % ncol == 0:
            ax.set_ylabel("Elevation (mRL)")
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        if i // ncol == nrow - 1:
            ax.set_xlabel("Temperature (degC)")
    fig.legend(handles=T_PROFILE_HANDLES, loc="lower center", ncol=8, fontsize=7, frameon=False)
    fig.suptitle(f"F2 Temperature profiles, group '{gname}': temperature (degC) against elevation (mRL), "
                 "identical axes", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 0.96))
    p = savefig(fig, CALFIG, f"profile_panels_{safe(gname)}.png")
    index_fig("F2", p, f"temperature small multiples for group {gname}", ", ".join(gw),
              f"{fmt(Z_TOP, 0)} to {fmt(Z_BOT, 0)} mRL", SRC_WELL, "obs")


# F3 elevation-band symbol maps
def band_symbols(ax, bd, var):
    sub = EB[(EB.band_top_mRL == bd[0])]
    zmid = 0.5 * (bd[0] + bd[1])
    step = 10.0 if var == "T" else 2.0
    col = "T_minus_field_median_degC" if var == "T" else "P_minus_trend_bar"
    for w in WELLS:
        x, y = track_xy(w, zmid)
        rw = sub[sub.well == w]
        if rw.empty:
            ax.plot(x, y, "o", mfc="none", mec="0.5", ms=5, zorder=10)
            continue
        v = float(rw[col].iloc[0])
        k = int(min(4, abs(v) // step))
        if k == 0:
            ax.plot(x, y, "x", color="k", ms=7, mew=1.5, zorder=11)
        else:
            ax.plot(x, y, "o", color="tab:red" if v > 0 else "tab:blue", ms=4 + 3 * k, mec="k", mew=0.5,
                    alpha=0.8, zorder=11)
        ax.annotate(f"{v:+.1f}", (x, y), xytext=(0, -7 - 1.5 * k), textcoords="offset points", ha="center",
                    va="top", fontsize=6, zorder=12)


def band_legend(var):
    u = "degC" if var == "T" else "bar"
    s = 10 if var == "T" else 2
    h = [Line2D([], [], marker="x", color="k", ls="", ms=7, mew=1.5, label=f"within +/-{s} {u}")]
    for k in range(1, 5):
        lab = f"{s * k}-{s * (k + 1)} {u}" if k < 4 else f"> {s * 4} {u}"
        h.append(Line2D([], [], marker="o", color="tab:red", mec="k", ls="", ms=4 + 3 * k,
                        label=(f"+{lab}" if k < 4 else lab) + " (higher)"))
    for k in range(1, 5):
        lab = f"{s * k}-{s * (k + 1)} {u}" if k < 4 else f"> {s * 4} {u}"
        h.append(Line2D([], [], marker="o", color="tab:blue", mec="k", ls="", ms=4 + 3 * k,
                        label=(f"-{lab}" if k < 4 else lab) + " (lower)"))
    h.append(Line2D([], [], marker="o", mfc="none", mec="0.5", ls="", ms=5, label="no data in band"))
    return h


def band_maps(var, folder, fname, interpolated=False):
    nb = len(BANDS)
    ncol = 4
    nrow = int(math.ceil(nb / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.3 * ncol, 4.3 * nrow + 1.2))
    axs = np.atleast_1d(axs).ravel()
    col = "T_minus_field_median_degC" if var == "T" else "P_minus_trend_bar"
    vmax = 50.0 if var == "T" else 10.0
    for i, bd in enumerate(BANDS):
        ax = axs[i]
        base_layer(ax, WELL_EXTENT, label_wells=False, small=True, fs=6)
        sub = EB[EB.band_top_mRL == bd[0]]
        if interpolated and len(sub) >= 3:
            tri = mtri.Triangulation(sub.x_mid.to_numpy(), sub.y_mid.to_numpy())
            xi = np.linspace(WELL_EXTENT[0], WELL_EXTENT[1], 200)
            yi = np.linspace(WELL_EXTENT[2], WELL_EXTENT[3], 200)
            XI, YI = np.meshgrid(xi, yi)
            ZI = mtri.LinearTriInterpolator(tri, sub[col].to_numpy())(XI, YI)
            levels = np.arange(-vmax, vmax + 0.01, 10.0 if var == "T" else 1.0)
            ax.add_patch(plt.Rectangle((WELL_EXTENT[0], WELL_EXTENT[2]), WELL_EXTENT[1] - WELL_EXTENT[0],
                                       WELL_EXTENT[3] - WELL_EXTENT[2], fc="none", ec="0.75", hatch="///", lw=0,
                                       zorder=0.5))
            cf = ax.contourf(XI, YI, ZI, levels=levels, cmap=DCMAP, extend="both", zorder=1.5, alpha=0.9)
            hull = sub[["x_mid", "y_mid"]].to_numpy()
            from scipy.spatial import ConvexHull
            hh = ConvexHull(hull)
            hp_ = hull[hh.vertices]
            ax.plot(np.r_[hp_[:, 0], hp_[0, 0]], np.r_[hp_[:, 1], hp_[0, 1]], "k--", lw=0.8, zorder=9)
        band_symbols(ax, bd, var)
        if var == "T":
            med = float(sub.T_field_median_degC.iloc[0])
            ax.set_title(f"{bd[0]:.0f} to {bd[1]:.0f} mRL: field median {med:.1f} degC", fontsize=8)
        else:
            ax.set_title(f"{bd[0]:.0f} to {bd[1]:.0f} mRL: P minus trend (bar)", fontsize=8)
    for ax in axs[nb:]:
        ax.axis("off")
    hl = band_legend(var) + BASE_HANDLES
    if interpolated:
        hl = [Patch(fc="none", ec="0.6", hatch="///", label="no data (outside convex hull)"),
              Line2D([], [], color="k", ls="--", label="convex hull of wells (data envelope)")] + hl
        fig.colorbar(cf, ax=list(axs[:nb]), shrink=0.6, label=("T minus field median (degC)" if var == "T"
                                                                  else "P minus trend (bar)"))
    fig.legend(handles=hl, loc="lower center", ncol=6, fontsize=7, frameon=False)
    if var == "T":
        ttl = ("F3 Temperature by 500 m elevation band: band-mean T minus the field median of band means "
               "(degC), symbols at each track's mid-band position")
    else:
        ttl = ("F3 Pressure by 500 m elevation band: band-mean of (P minus trend of gauge-verified feed zones) "
               "(bar), symbols at each track's mid-band position")
    if interpolated:
        ttl = ("[interpolated, not data] " + ttl.replace("F3 ", "") +
               "\nlinear Delaunay interpolation between wells; colours are not constrained by neighbouring wells "
               "except at the symbols")
    fig.suptitle(ttl, fontsize=9.5)
    fig.subplots_adjust(left=0.03, right=0.9 if interpolated else 0.98, top=0.9, bottom=0.14, wspace=0.08,
                        hspace=0.18)
    return savefig(fig, folder, fname)


p = band_maps("T", CALFIG, "maps_temperature_bands.png")
index_fig("F3", p, "band-mean temperature minus field median per 500 m band (symbol classes of 10 degC)", "all",
          f"bands {BANDS[0][0]:.0f} to {BANDS[-1][1]:.0f} mRL", SRC_WELL + "; T_P_at_elevation_bands.csv", "obs")
p = band_maps("P", CALFIG, "maps_pressure_bands.png")
index_fig("F3", p, "band-mean pressure minus trend of gauge-verified feed zones (symbol classes of 2 bar)", "all",
          f"bands {BANDS[0][0]:.0f} to {BANDS[-1][1]:.0f} mRL", SRC_WELL + "; T_P_at_elevation_bands.csv", "mixed")


# F4 data sections
def section_figure(ln, var, folder, fname, interpolated=False):
    fig = plt.figure(figsize=(15, 8.2))
    ax = fig.add_axes([0.06, 0.1, 0.62, 0.74])
    if not interpolated:
        sw, meshes, kept = draw_section(ax, fig, ln, var)
    else:
        sw = section_wells(ln)
        ss, zs = surface_profile(ln)
        hatch_nodata(ax, ln, ss, zs)
        allS = np.concatenate([d["s"] for d in sw])
        allZ = np.concatenate([d["z"] for d in sw])
        allV = np.concatenate([d["T" if var == "T" else "dP"] for d in sw])
        tri = mtri.Triangulation(allS, allZ)
        SI, ZI = np.meshgrid(np.linspace(*section_xlim(ln), 300), np.linspace(Z_BOT, Z_TOP, 300))
        VI = mtri.LinearTriInterpolator(tri, allV)(SI, ZI)
        # never colour above the model surface
        surf = np.interp(SI, ss[np.isfinite(zs)], zs[np.isfinite(zs)])
        VI = np.ma.masked_where(ZI > surf, VI)
        if var == "T":
            ax.contourf(SI, ZI, VI, levels=T_LEVELS, cmap=TCMAP, norm=TNORM, zorder=3, extend="both")
            cs = ax.contour(SI, ZI, VI, levels=ISO_LEVELS, colors="k",
                            linewidths=[1.6 if lv in (100, 200) else 0.6 for lv in ISO_LEVELS], zorder=4)
            ax.clabel(cs, fmt="%d", fontsize=6)
        else:
            ax.contourf(SI, ZI, VI, levels=P_LEVELS, cmap=DCMAP, norm=PNORM, zorder=3, extend="both")
            ax.contour(SI, ZI, VI, levels=[0], colors="k", linewidths=1.4, zorder=4)
        for m in envelope_meshes(sw, "T"):
            ax.plot(m["outline"][:, 0], m["outline"][:, 1], color="k", ls="--", lw=1.0, zorder=5)
        draw_wells_on_section(ax, ln, sw, "T" if var == "T" else "dP", norm=TNORM if var == "T" else PNORM,
                              cmap=TCMAP if var == "T" else DCMAP)
        section_axes_setup(ax, ln)
        surface_ticks(ax, ln)
    ve = vertical_exaggeration(ax, fig)
    cax = fig.add_axes([0.70, 0.1, 0.012, 0.74])
    if var == "T":
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=TNORM, cmap=TCMAP), cax=cax)
        cb.set_label("Temperature (degC), 10 degC steps")
    else:
        cb = fig.colorbar(plt.cm.ScalarMappable(norm=PNORM, cmap=DCMAP), cax=cax)
        cb.set_label("P minus trend of gauge-verified feed zones (bar), 1 bar steps; red = higher")
    axm = fig.add_axes([0.76, 0.52, 0.22, 0.36])
    base_layer(axm, WELL_EXTENT, label_wells=False, small=True, fault_labels=False)
    draw_lines(axm, [ln], fs=6)
    for d in sw:
        axm.plot(*WH[d["well"]], "o", color=WELL_COLOR[d["well"]], ms=5, zorder=12)
        label_well(axm, d["well"], color=WELL_COLOR[d["well"]], fs=6, dist=4)
    hl = section_legend_handles(var)
    if interpolated:
        hl = [Line2D([], [], color="k", ls="--", lw=1.0,
                     label="data envelope (F4); colours outside it are\nnot constrained by neighbouring wells"),
              Patch(fc="0.93", ec="0.6", hatch="///", label="no data (outside convex hull)")] + hl[1:]
    fig.legend(handles=hl, loc="lower left", bbox_to_anchor=(0.745, 0.02), fontsize=6.5, frameon=False)
    what = "Temperature (degC)" if var == "T" else "Pressure minus hydrostatic trend of gauge-verified feed zones (bar)"
    ttl = (f"{what} on section {ln['name']}-{ln['name']}' ({fmt(ln['a'][0], 0)}, {fmt(ln['a'][1], 0)}) to "
           f"({fmt(ln['b'][0], 0)}, {fmt(ln['b'][1], 0)}); wells within {SECTION_HALFWIDTH:.0f} m projected; "
           f"vertical exaggeration {ve:.1f}x")
    if interpolated:
        ttl = "[interpolated, not data] " + ttl + "\nlinear Delaunay interpolation of all survey points on the line"
    else:
        ttl = "F4 " + ttl + "\ncolour only inside the data envelope (linear between neighbouring wells)"
    fig.suptitle(ttl, fontsize=9)
    ln.setdefault("ve", {})[(var, interpolated)] = ve
    return savefig(fig, folder, fname)


for ln in LINES:
    wl = ", ".join(sorted(ln["wells"], key=lambda k: ln["wells"][k][0]))
    p = section_figure(ln, "T", CALFIG, f"section_{ln['name']}.png")
    index_fig("F4", p, f"temperature section {ln['name']}-{ln['name']}' inside the data envelope", wl,
              f"line {ln['name']}-{ln['name']}'", SRC_WELL + "; geometry_dat; mapping", "mixed")
    p = section_figure(ln, "P", CALFIG, f"section_{ln['name']}_P.png")
    index_fig("F4", p, f"pressure-minus-trend section {ln['name']}-{ln['name']}' inside the data envelope", wl,
              f"line {ln['name']}-{ln['name']}'", SRC_WELL + "; geometry_dat; mapping", "mixed")
    section_figure(ln, "T", ADDFIG, f"section_{ln['name']}_interpolated.png", interpolated=True)
    section_figure(ln, "P", ADDFIG, f"section_{ln['name']}_P_interpolated.png", interpolated=True)
band_maps("T", ADDFIG, "maps_temperature_interpolated.png", interpolated=True)
band_maps("P", ADDFIG, "maps_pressure_interpolated.png", interpolated=True)

# F5 P-T state
fig, ax = plt.subplots(figsize=(10, 7.5))
tt = np.linspace(0.01, 373.946, 400)
ax.plot(tt, [psat(t) for t in tt], color="k", lw=1.5, label="IF97 saturation curve (Region 4)")
ax.plot([373.946], [P_CRIT], "k*", ms=12, label="critical point (373.95 degC, 220.64 bar)")
for w in WELLS:
    g = wsurvey(w)
    ax.plot(g.T_degC, g.P_bar, "o", color=WELL_COLOR[w], ms=4, mec="k", mew=0.2, label=w)
    f = FZ[FZ.well == w]
    ax.plot(f.T_degC, f.P_bar, "s", color=WELL_COLOR[w], ms=10, mec="k", mew=1.0)
ax.plot([], [], "s", color="w", mec="k", ms=10, label="feed zone (larger square)")
ax.set_xlabel("Temperature (degC)")
ax.set_ylabel("Pressure (bar)")
ax.set_xlim(0, 400)
ax.set_ylim(0, math.ceil(S.P_bar.max() / 50) * 50)
ax.text(300, 40, "two-phase / steam\n(left of curve: liquid)", fontsize=8, ha="center")
ax.grid(alpha=0.3)
ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=7, frameon=False)
ax.set_title("F5 P-T state of all survey points (bar, degC) against the IF97 saturation curve")
fig.tight_layout()
p = savefig(fig, CALFIG, "PT_state.png")
index_fig("F5", p, "survey P-T points and feed zones against the IF97 saturation curve", "all", "all elevations",
          SRC_WELL, "obs")

# ------------------------------------------------------------------------------------------
# 19. Conceptual model elements (data-derived), pressure connections
# ------------------------------------------------------------------------------------------
UPFLOW_WELLS = [w for w in WELLS if WSI.loc[w, "regime"] == "upflow" or
                WSI.loc[w, "T_at_TD_degC"] >= 0.9 * FIELD_TMAX]
OUTFLOW_WELLS = [w for w in WELLS if WSI.loc[w, "regime"] == "outflow" and w not in UPFLOW_WELLS]
COLD_RECHARGE = [w for w in WELLS if WSI.loc[w, "regime"] == "recharge/downflow"]
CAP_WELLS = [w for w in WELLS if WSI.loc[w, "Tmax_degC"] >= UPFLOW_TMAX and
             np.isfinite(WSI.loc[w, "shallow_gradient_degC_per_km"])]
UP_CENTRE = np.mean([WH[w] for w in UPFLOW_WELLS], axis=0)


def cap_interval(w):
    """Cap = from the shallowest liquid-filled point to where T first reaches 90 % of the well's Tmax."""
    g = TGT[TGT.well == w].sort_values("elevation_mRL", ascending=False)
    z, t = g.elevation_mRL.to_numpy(), g.T_degC.to_numpy()
    top = float(g[g.liquid_filled].elevation_mRL.max())
    base = crossing(z, t, 0.9 * WSI.loc[w, "Tmax_degC"])
    return top, base


CAP = {w: cap_interval(w) for w in CAP_WELLS}
NB_INT = {}
for w in WELLS:
    g = wsurvey(w)
    nb = g[(g.boiling_margin_degC < NEAR_BOIL)]
    if len(nb):
        NB_INT[w] = (float(nb.elevation_mRL.max()), float(nb.elevation_mRL.min()), float(nb.boiling_margin_degC.min()))

# pressure connections between neighbouring wells (Delaunay of wellheads) at a common elevation
PCON = []
_tri = mtri.Triangulation(np.array([WH[w][0] for w in WELLS]), np.array([WH[w][1] for w in WELLS]))
_pairs = sorted({tuple(sorted((int(a), int(b)))) for a, b in _tri.edges})
for i, j in _pairs:
    w1, w2 = WELLS[i], WELLS[j]
    g1, g2 = wsurvey(w1), wsurvey(w2)
    top = min(g1.elevation_mRL.max(), g2.elevation_mRL.max())
    bot = max(g1.elevation_mRL.min(), g2.elevation_mRL.min())
    cands_z = [z for z in np.arange(0.0, -3001.0, -250.0) if bot <= z <= top]
    if not cands_z:
        continue
    zc = min(cands_z, key=lambda z: (abs(z + 1000.0), z))
    p1, p2 = interp_w(w1, "P_bar", zc), interp_w(w2, "P_bar", zc)
    PCON.append(dict(well_1=w1, well_2=w2, elevation_mRL=zc, P1_bar=r(p1, 1), P2_bar=r(p2, 1),
                     dP_2_minus_1_bar=r(p2 - p1, 1),
                     distance_m=r(float(np.linalg.norm(WH[w1] - WH[w2])), 0),
                     P1_minus_trend_bar=r(p1 - trend(zc), 1), P2_minus_trend_bar=r(p2 - trend(zc), 1)))
PCONDF = pd.DataFrame(PCON)
write_csv(PCONDF, os.path.join(ADD, "pressure_connections.csv"))

CM_ROWS = []


def cm(no, element, where, observation, claim, symbol):
    CM_ROWS.append(dict(no=no, element=element, location=where, supporting_observation=observation,
                        claim_type=claim, symbol=symbol))


cm(1, "Deep hot upflow", f"beneath {' and '.join(UPFLOW_WELLS)}",
   "; ".join(f"{w}: Tmax {fmt(WSI.loc[w, 'Tmax_degC'], 1)} degC at {fmt(WSI.loc[w, 'Tmax_elev_mRL'], 1)} mRL, "
             f"T at TD {fmt(WSI.loc[w, 'T_at_TD_degC'], 1)} degC at {fmt(WSI.loc[w, 'TD_mRL'], 1)} mRL, regime "
             f"{WSI.loc[w, 'regime']}" for w in UPFLOW_WELLS), "interp",
   "red upward arrow (solid where the hot base is measured, dashed below the deepest data)")
cm(2, "Lateral outflow", "away from the upflow towards " + ", ".join(OUTFLOW_WELLS),
   "; ".join(f"{w}: inversion {fmt(WSI.loc[w, 'inversion_degC'], 1)} degC below "
             f"{fmt(WSI.loc[w, 'inversion_top_mRL'], 1)} mRL" for w in OUTFLOW_WELLS), "interp",
   "orange horizontal arrow at the inversion-top elevation")
cm(3, "Cooler water beneath the outflow (lateral inflow / downflow at depth)",
   "below the inversion bases of " + ", ".join(OUTFLOW_WELLS),
   "; ".join(f"{w}: {fmt(WSI.loc[w, 'inversion_degC'] * 0 + wsurvey(w).set_index('elevation_mRL').T_degC[WSI.loc[w, 'inversion_base_mRL']], 1)} degC at "
             f"{fmt(WSI.loc[w, 'inversion_base_mRL'], 1)} mRL" for w in OUTFLOW_WELLS), "interp",
   "blue arrow (dashed)")
cm(4, "Cold recharge / downflow", ", ".join(COLD_RECHARGE),
   "; ".join(f"{w}: Tmax {fmt(WSI.loc[w, 'Tmax_degC'], 1)} degC, shallow gradient "
             f"{fmt(WSI.loc[w, 'shallow_gradient_degC_per_km'], 1)} degC/km" for w in COLD_RECHARGE), "interp",
   "blue downward arrow (dashed)")
cm(5, "Cap (low permeability, conductive gradient)", ", ".join(CAP_WELLS),
   "; ".join(f"{w}: shallow gradient {fmt(WSI.loc[w, 'shallow_gradient_degC_per_km'], 1)} degC/km, cap "
             f"{fmt(CAP[w][0], 1)} to {fmt(CAP[w][1], 1)} mRL" for w in CAP_WELLS), "interp",
   "grey band, top and base labelled")
cm(6, "Near-boiling zone", ", ".join(NB_INT) or "none",
   "; ".join(f"{w}: boiling margin < 5 degC from {fmt(a, 1)} to {fmt(b, 1)} mRL (min {fmt(m, 2)} degC)"
             for w, (a, b, m) in NB_INT.items()) or "no point with boiling margin < 5 degC", "obs", "stipple")
cm(7, "Isotherms 100 and 200 degC", "all wells",
   "z_100C_mRL and z_200C_mRL in well_summary.csv", "obs", "heavy contour lines from F4")
pc_ok = PCONDF[(PCONDF.P1_minus_trend_bar.abs() <= TOL) & (PCONDF.P2_minus_trend_bar.abs() <= TOL)]
cm(8, "Pressure connection between wells", "field-wide at depth",
   f"{len(pc_ok)} of {len(PCONDF)} neighbouring-well pairs have both pressures within {fmt(TOL, 1)} bar of the "
   f"gauge-verified trend at the comparison elevation; feed-zone residuals in feedzones.csv", "interp",
   "thin dashed line labelled with the pressure difference")
CMDF = pd.DataFrame(CM_ROWS)
write_csv(CMDF, os.path.join(ADD, "conceptual_model_elements.csv"))
CAPDF = pd.DataFrame([dict(well=w, cap_top_mRL=r(a, 1), cap_base_mRL=r(b, 1), cap_thickness_m=r(a - b, 1),
                           shallow_gradient_degC_per_km=WSI.loc[w, "shallow_gradient_degC_per_km"])
                      for w, (a, b) in CAP.items()])
write_csv(CAPDF, os.path.join(ADD, "cap_intervals.csv"))


def num_badge(ax, x, y, n, color="k"):
    ax.text(x, y, str(n), fontsize=8, fontweight="bold", color="white", ha="center", va="center", zorder=20,
            bbox=dict(boxstyle="circle,pad=0.2", fc=color, ec="k", lw=0.5))


def arrow(ax, p0, p1, color, solid=True, lw=2.5, z=18):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=16, color=color, lw=lw,
                                 linestyle="-" if solid else (0, (4, 2)), zorder=z))


def concept_section(ax, fig, ln):
    sw, meshes, kept = draw_section(ax, fig, ln, "T")
    pos = {d["well"]: d for d in sw}
    s_up = [pos[w]["s_head"] for w in UPFLOW_WELLS if w in pos]
    s_up_c = np.mean(s_up) if s_up else float(np.dot(UP_CENTRE - ln["a"], ln["u"]))
    # 1 upflow
    for w in UPFLOW_WELLS:
        if w not in pos:
            continue
        d = pos[w]
        zb = float(d["z"].min())
        arrow(ax, (d["s_head"] - 120, zb), (d["s_head"] - 120, min(zb + 500, Z_TOP - 300)), "red", solid=True)
        arrow(ax, (d["s_head"] - 120, Z_BOT + 150), (d["s_head"] - 120, zb - 30), "red", solid=False)
        num_badge(ax, d["s_head"] - 260, zb - 200, 1, "red")
    # 2 outflow, 3 cooler water beneath
    for w in OUTFLOW_WELLS:
        if w not in pos:
            continue
        d = pos[w]
        sgn = 1 if d["s_head"] >= s_up_c else -1
        zt_ = WSI.loc[w, "inversion_top_mRL"]
        zb_ = WSI.loc[w, "inversion_base_mRL"]
        arrow(ax, (d["s_head"] - sgn * 450, zt_), (d["s_head"] + sgn * 450, zt_), "darkorange", solid=True)
        num_badge(ax, d["s_head"] + sgn * 550, zt_, 2, "darkorange")
        arrow(ax, (d["s_head"] + sgn * 450, zb_), (d["s_head"] - sgn * 350, zb_), "royalblue", solid=False)
        num_badge(ax, d["s_head"] + sgn * 550, zb_, 3, "royalblue")
    for w in COLD_RECHARGE:
        if w in pos:
            d = pos[w]
            arrow(ax, (d["s_head"] + 150, Z_TOP - 250), (d["s_head"] + 150, Z_TOP - 1000), "royalblue", solid=False)
            num_badge(ax, d["s_head"] + 300, Z_TOP - 300, 4, "royalblue")
    # 5 cap band between cap wells on the line
    cw = sorted([w for w in CAP_WELLS if w in pos], key=lambda w: pos[w]["s_head"])
    if cw:
        ss_ = [pos[w]["s_head"] for w in cw]
        tops = [CAP[w][0] for w in cw]
        bases = [CAP[w][1] for w in cw]
        if len(cw) == 1:
            ss_ = [ss_[0] - 150, ss_[0] + 150]
            tops, bases = tops * 2, bases * 2
        ax.fill_between(ss_, bases, tops, color="0.4", alpha=0.35, zorder=8, lw=0)
        for k_, w in enumerate(cw):
            sx, ha_ = (pos[w]["s_head"] + 60, "left") if k_ % 2 == 0 else (pos[w]["s_head"] - 60, "right")
            ax.text(sx, CAP[w][1], f"{w} cap base {fmt(CAP[w][1], 0)}", fontsize=6, va="top", ha=ha_,
                    zorder=19, bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.5))
            ax.text(sx, CAP[w][0] - (0 if k_ % 2 == 0 else 60), f"cap top {fmt(CAP[w][0], 0)}", fontsize=6,
                    va="bottom", ha=ha_, zorder=19, bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.5))
        num_badge(ax, ss_[0] - 250, np.mean(tops[:1] + bases[:1]), 5, "0.4")
    for w, (a, b, m) in NB_INT.items():
        if w in pos:
            num_badge(ax, pos[w]["s_head"] + 260, 0.5 * (a + b), 6, "0.2")
    # 7 isotherm badge at the 100 degC crossing of the first well with one
    for d in sw:
        z1 = WSI.loc[d["well"], "z_100C_mRL"]
        if np.isfinite(z1):
            num_badge(ax, d["s_head"] - 200, z1, 7, "k")
            break
    # 8 pressure connections between neighbouring wells on the line
    for k_, (a_, b_) in enumerate(zip(sw, sw[1:])):
        top = min(a_["z"].max(), b_["z"].max())
        bot = max(a_["z"].min(), b_["z"].min())
        cz = [z for z in np.arange(0.0, -3001.0, -250.0) if bot <= z <= top]
        if not cz:
            continue
        zc = min(cz, key=lambda z: (abs(z + 1000.0), z))
        dp = interp_w(b_["well"], "P_bar", zc) - interp_w(a_["well"], "P_bar", zc)
        sa = float(np.interp(zc, a_["z"][::-1], a_["s"][::-1]))
        sb = float(np.interp(zc, b_["z"][::-1], b_["s"][::-1]))
        ax.plot([sa, sb], [zc, zc], color="k", ls="--", lw=0.9, zorder=17)
        ax.text(0.5 * (sa + sb), zc + 25 + 90 * (k_ % 2), f"dP {dp:+.1f} bar @ {zc:.0f} mRL", fontsize=6,
                ha="center",
                va="bottom", zorder=19, bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.5))
        num_badge(ax, 0.5 * (sa + sb), zc - 70, 8, "k")
    return sw


def concept_plan(ax):
    base_layer(ax, WELL_EXTENT, fs=7)
    draw_lines(ax)
    for w in UPFLOW_WELLS:
        ax.plot(*WH[w], marker="^", color="red", ms=15, mec="k", zorder=15)
    num_badge(ax, UP_CENTRE[0] - 350, UP_CENTRE[1] + 150, 1, "red")
    for w in OUTFLOW_WELLS:
        v = WH[w] - UP_CENTRE
        L = np.linalg.norm(v)
        p0 = UP_CENTRE + v / L * 250
        p1 = WH[w] - v / L * 120
        arrow(ax, tuple(p0), tuple(p1), "darkorange", solid=False, lw=2)
        mid = 0.5 * (p0 + p1)
        ax.text(mid[0], mid[1], f"{fmt(WSI.loc[w, 'inversion_top_mRL'], 0)} mRL", fontsize=6, color="darkorange",
                zorder=19, bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.4))
    if OUTFLOW_WELLS:
        w0 = OUTFLOW_WELLS[0]
        mid = UP_CENTRE + 0.65 * (WH[w0] - UP_CENTRE)
        num_badge(ax, mid[0] + 150, mid[1] - 150, 2, "darkorange")
        for w in OUTFLOW_WELLS:
            ax.plot(*WH[w], marker="o", ms=11, mfc="none", mec="royalblue", mew=1.6, ls="", zorder=14)
        num_badge(ax, WH[w0][0] + 250, WH[w0][1] - 250, 3, "royalblue")
    for w in COLD_RECHARGE:
        ax.plot(*WH[w], marker="v", color="royalblue", ms=14, mec="k", zorder=15)
        num_badge(ax, WH[w][0] + 300, WH[w][1] + 250, 4, "royalblue")
    if len(CAP_WELLS) >= 3:
        from scipy.spatial import ConvexHull
        pts = np.vstack([WH[w] + 350 * np.array([np.cos(a), np.sin(a)]) for w in CAP_WELLS
                         for a in np.linspace(0, 2 * np.pi, 16, endpoint=False)])
        h = ConvexHull(pts)
        hp_ = pts[h.vertices]
        ax.add_patch(MPoly(hp_, closed=True, fc="0.5", alpha=0.18, ec="0.3", ls="--", lw=1.2, zorder=6))
        num_badge(ax, hp_[:, 0].min() + 100, hp_[:, 1].min() + 100, 5, "0.4")
    for w in NB_INT:
        ax.add_patch(plt.Circle(tuple(WH[w]), 220, fc="none", ec="k", hatch="...", lw=0.6, zorder=13))
        num_badge(ax, WH[w][0] - 330, WH[w][1] - 250, 6, "0.2")
    for _, rw in PCONDF.iterrows():
        a_, b_ = WH[rw.well_1], WH[rw.well_2]
        ax.plot([a_[0], b_[0]], [a_[1], b_[1]], color="k", ls="--", lw=0.5, alpha=0.6, zorder=10)
        m = 0.5 * (a_ + b_)
        ax.text(m[0], m[1], f"{rw.dP_2_minus_1_bar:+.1f}", fontsize=5.5, ha="center", va="center", zorder=11,
                bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.3))
    num_badge(ax, WELL_EXTENT[1] - 500, WELL_EXTENT[3] - 500, 8, "k")


CONCEPT_HANDLES = [
    FancyArrowPatch((0, 0), (1, 0), color="red"), ]
CONCEPT_LEGEND = [Line2D([], [], color="red", lw=2.5, marker="^", ms=8, label="1 deep hot upflow (solid: measured; dashed: inferred)"),
                  Line2D([], [], color="darkorange", lw=2.5, marker=">", ms=7, label="2 lateral outflow at its depth"),
                  Line2D([], [], color="royalblue", lw=2.5, ls="--", marker="<", ms=7, label="3 cooler water beneath outflow"),
                  Line2D([], [], color="royalblue", lw=0, marker="v", ms=9, mec="k", label="4 cold recharge / downflow"),
                  Patch(fc="0.4", alpha=0.35, label="5 cap (top and base labelled)"),
                  Patch(fc="none", ec="k", hatch="...", label="6 two-phase / near-boiling zone"),
                  Line2D([], [], color="k", lw=1.6, label="7 isotherms 100 / 200 degC (from F4)"),
                  Line2D([], [], color="k", lw=0.9, ls="--", label="8 pressure connection, dP at common elevation (bar)")]

F6_LINE = sorted(LINES, key=lambda l: (-len(l["wells"]), -sum(WSI.loc[w, "Tmax_degC"] for w in l["wells"]),
                                        l["name"]))[0]


def concept_figure(ln, folder, fname, plan=True):
    fig = plt.figure(figsize=(17, 8.5))
    if plan:
        axp = fig.add_axes([0.06, 0.12, 0.28, 0.76])
        concept_plan(axp)
        axp.set_title("Plan view (numbers key to the conceptual-model table)", fontsize=8)
        ax = fig.add_axes([0.38, 0.12, 0.47, 0.72])
    else:
        ax = fig.add_axes([0.06, 0.12, 0.75, 0.72])
    concept_section(ax, fig, ln)
    ve = vertical_exaggeration(ax, fig)
    ax.set_title(f"Section {ln['name']}-{ln['name']}' (same line as F4), vertical exaggeration {ve:.1f}x", fontsize=8,
                 pad=62)
    fig.legend(handles=CONCEPT_LEGEND + section_legend_handles("T")[:3], loc="center right",
               bbox_to_anchor=(0.995, 0.5), fontsize=6.5, frameon=False)
    fig.suptitle(f"F6 Conceptual model from the data alone: plan view and section {ln['name']}-{ln['name']}' "
                 f"(temperature background in degC, elevation in mRL); solid = measured feature, dashed = inferred",
                 fontsize=9.5)
    return savefig(fig, folder, fname)


p = concept_figure(F6_LINE, CALFIG, "conceptual_model.png")
index_fig("F6", p, "conceptual model: upflow, outflow, recharge, cap, near-boiling zone, isotherms, pressure "
          "connections (plan and section)", "all", f"plan view; line {F6_LINE['name']}-{F6_LINE['name']}'",
          SRC_WELL + "; mapping; derived tables", "interp")
for ln in LINES:
    if ln is not F6_LINE:
        concept_figure(ln, ADDFIG, f"conceptual_model_section_{ln['name']}.png", plan=False)

# ------------------------------------------------------------------------------------------
# F7 gauge check
# ------------------------------------------------------------------------------------------
gauged = [g for g in GAUGE_ROWS]
if gauged:
    n = len(gauged)
    ncol = 2
    nrow = int(math.ceil(n / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(14, 4.2 * nrow + 0.8), squeeze=False)
    PERIOD = (min(min(g["first_gauge_date"] for g in gauged), FIELD_FIRST_PROD) - dt.timedelta(days=180),
              max(g["series"][-1][0] for g in gauged) + dt.timedelta(days=180))
    for ax, g in zip(axs.ravel(), gauged):
        w = g["well"]
        wd = WD[w]
        flow = wd.get("Production_Mass_Flow") or wd.get("Injection_Mass_Flow") or []
        for (d0, v0), (d1, _) in zip(flow, flow[1:] + [(flow[-1][0] + dt.timedelta(days=30), 0)] if flow else []):
            if abs(v0) > 0:
                ax.axvspan(d0, d1, color="0.88", lw=0, zorder=0)
        ds = [d for d, _ in g["series"]]
        ps = [v for _, v in g["series"]]
        ax.plot(ds, ps, "o-", color=WELL_COLOR[w], ms=5, mec="k", lw=1, zorder=5, label="gauge pressure")
        ax.axhspan(g["first_gauge_P_bar"] - GAUGE_PRECISION, g["first_gauge_P_bar"] + GAUGE_PRECISION,
                   color="tab:green", alpha=0.3, zorder=1, label=f"first reading +/- {GAUGE_PRECISION} bar")
        ax.axhline(g["survey_P_at_gauge_bar"], color="k", ls="--", lw=1.2, zorder=4,
                   label=f"survey P at gauge elevation ({fmt(g['survey_P_at_gauge_bar'], 2)} bar), survey dated "
                         f"{iso(wd['P_date'])}")
        ax.axvline(wd["first_flow"], color="tab:purple", lw=1.5, zorder=3,
                   label=f"{w} first {'production' if wd['first_prod'] else 'injection'} {iso(wd['first_flow'])}")
        ax.axvline(FIELD_FIRST_PROD, color="tab:red", lw=1.2, ls=":", zorder=3,
                   label=f"first field production {iso(FIELD_FIRST_PROD)}")
        ax.axvline(wd["P_date"], color="0.3", lw=0.8, ls="-.", zorder=3, label="survey date label")
        ax.set_xlim(*PERIOD)
        ax.set_ylabel("Pressure (bar)")
        ax.set_title(f"{w} gauge at {fmt(g['gauge_elev_mRL'], 1)} mRL: {NS[w]['verdict']} "
                     f"({NS[w]['confidence']} confidence)", fontsize=8.5)
        ax.legend(fontsize=6.5, loc="center right", frameon=True)
        ax.grid(alpha=0.3)
    for ax in axs.ravel()[n:]:
        ax.axis("off")
    fig.suptitle(f"F7 Natural-state gauge check: transient pressure (bar) against date, "
                 f"{iso(PERIOD[0])} to {iso(PERIOD[1])}; grey = well flowing", fontsize=10)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    p = savefig(fig, CALFIG, "gauge_check.png")
    index_fig("F7", p, "gauge pressure history vs survey pressure at gauge elevation, with production timing",
              ", ".join(g["well"] for g in gauged), "gauge elevations " + ", ".join(fmt(g["gauge_elev_mRL"], 1)
                                                                             for g in gauged) + " mRL",
              SRC_WELL + " (transient pressure, flows)", "obs")

FIDX = pd.DataFrame(FIG_INDEX)
write_csv(FIDX, os.path.join(CAL, "figures_index.csv"))

# ------------------------------------------------------------------------------------------
# 20. Supporting tables written after the figures
# ------------------------------------------------------------------------------------------
for ln in LINES:
    i = SLDF.index[SLDF.line == ln["name"]][0]
    ve = ln.get("ve", {})
    SLDF.loc[i, "vertical_exaggeration_T"] = r(ve.get(("T", False), np.nan), 1)
    SLDF.loc[i, "vertical_exaggeration_P"] = r(ve.get(("P", False), np.nan), 1)
    SLDF.loc[i, "used_in_F6"] = ln is F6_LINE
write_csv(SLDF, os.path.join(ADD, "section_lines.csv"))

PCI = PC.set_index("well")
hs = HS.set_index(["well", "band"])
n_iso_pts = int(TGT.isothermal_section.sum())
PARAMS = [
    ("P_atm", P_ATM, "bar", "CLAUDE.md 4.4 water level"), ("P_critical", P_CRIT, "bar", "IF97"),
    ("g", G, "m/s2", "CLAUDE.md 4.3"), ("liquid_P_min", LIQ_P_MIN, "bar", "CLAUDE.md 4.3/4.4"),
    ("hydrostatic_boiling_margin_min", HYD_BM_MIN, "degC", "CLAUDE.md 4.3"),
    ("hydrostatic_window", HYD_WIN, "m", "CLAUDE.md 4.3"), ("isothermal_span", ISO_SPAN, "m", "CLAUDE.md 4.2"),
    ("isothermal_range", ISO_DT, "degC", "CLAUDE.md 4.2"), ("near_boiling_margin", NEAR_BOIL, "degC", "CLAUDE.md 4.3"),
    ("inversion_min", INV_MIN, "degC", "CLAUDE.md 4.4"), ("outflow_inversion", OUTFLOW_INV, "degC", "CLAUDE.md 4.4"),
    ("cold_Tmax", COLD_TMAX, "degC", "CLAUDE.md 4.3/4.4"), ("upflow_Tmax", UPFLOW_TMAX, "degC", "CLAUDE.md 4.4"),
    ("shallow_fit_span", SHALLOW_SPAN, "m", "CLAUDE.md 4.4"), ("deep_fit_span", DEEP_SPAN, "m", "CLAUDE.md 4.4"),
    ("wellhead_zone", WELLHEAD_ZONE, "m", "CLAUDE.md 5.5"), ("gauge_precision", GAUGE_PRECISION, "bar", "reading precision of the data"),
    ("tolerance_min", TOL_MIN, "bar", "CLAUDE.md 4.1"), ("tolerance_used", r(TOL, 2), "bar", "max(3 bar, 2 x RMS)"),
    ("trend_A", r(TREND_A, 2), "bar", "P = A + B z, gauge-verified feed zones"),
    ("trend_B", r(TREND_B, 5), "bar/m", "P = A + B z"), ("trend_rms", r(TREND_RMS, 2), "bar", "RMS residual"),
    ("trend_2xrms", r(2 * TREND_RMS, 2), "bar", "2 x RMS"), ("trend_implied_density", r(TREND_RHO, 1), "kg/m3", "-B/g"),
    ("trend_n_feedzones", len(fzp), "-", "gauge-verified feed zones"),
    ("trend_z_top", TREND_ZRANGE[1], "mRL", "shallowest verified feed zone"),
    ("trend_z_bottom", TREND_ZRANGE[0], "mRL", "deepest verified feed zone"),
    ("band_thickness", BAND, "m", "CLAUDE.md 8.2 F3"), ("section_half_width", SECTION_HALFWIDTH, "m", "CLAUDE.md 8.2 F4"),
    ("section_view_margin", SECTION_MARGIN, "m", "view crop around the wells"),
    ("injectate_pressure_assumed", INJ_P_ASSUMED, "bar", "[assumption] for T(h) of the injectate"),
    ("injectate_near_threshold", 10.0, "degC", "isothermal section 'near' the injectate temperature"),
    ("cap_base_fraction_of_Tmax", 0.9, "-", "cap base = first crossing of 0.9 x Tmax"),
    ("upflow_TD_fraction_of_field_Tmax", 0.9, "-", "upflow well if T at TD >= 0.9 x field Tmax"),
    ("field_Tmax", FIELD_TMAX, "degC", "well_summary"), ("T_scale_top", T_TOP, "degC", "figure colour scale"),
    ("elevation_axis_top", Z_TOP, "mRL", "figures"), ("elevation_axis_bottom", Z_BOT, "mRL", "figures"),
    ("P_deviation_scale", P_DEV_MAX, "bar", "figure colour scale"),
    ("IF97_max_validation_error", r(IF97_MAX_ERR, 4), "K", "validation"),
    ("Tsat_1.01325bar", r(tsat(1.01325), 2), "degC", "validation"), ("Tsat_100bar", r(tsat(100.0), 2), "degC", "validation"),
    ("Tsat_ref_1", 99.97, "degC", "CLAUDE.md 4.3"), ("Tsat_ref_2", 311.0, "degC", "CLAUDE.md 4.3"),
    ("conf_factor_high", 1.0, "-", "CLAUDE.md 5.5"), ("conf_factor_medium", 0.75, "-", "CLAUDE.md 5.5"),
    ("conf_factor_low", 0.5, "-", "CLAUDE.md 5.5"), ("n_survey_points", len(TGT), "-", "targets"),
    ("n_isothermal_points", n_iso_pts, "-", "targets"), ("n_feed_zones", len(FZ), "-", "feedzones"),
    ("n_hydrostatic_windows", len(HY), "-", "hydrostatic_check"), ("n_layers", len(LAYERS), "-", "geometry"),
    ("n_columns", len(COLS), "-", "geometry"), ("grid_top", GRID_TOP, "mRL", "geometry"),
    ("grid_bottom", LAYERS[max(LAYERS)]["bottom"], "mRL", "geometry"),
    ("cold_median_ratio_open", r(cold_med["open"], 4), "-", "hydrostatic"),
    ("cold_median_ratio_all", r(cold_med["all"], 4), "-", "hydrostatic"),
]
PARDF = pd.DataFrame(PARAMS, columns=["parameter", "value", "unit", "source"])
write_csv(PARDF, os.path.join(ADD, "parameters_used.csv"))

# ------------------------------------------------------------------------------------------
# 21. Report (CLAUDE.md 6)
# ------------------------------------------------------------------------------------------
L = []
A = L.append


def table(df, cols, heads=None, fmts=None):
    heads = heads or cols
    fmts = fmts or {}
    out = ["| " + " | ".join(heads) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, rw in df.iterrows():
        cells = []
        for c in cols:
            v = rw[c]
            if c in fmts:
                cells.append(fmt(v, fmts[c]) if pd.notna(v) else "")
            else:
                cells.append("" if (isinstance(v, float) and not np.isfinite(v)) else str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


def fig_md(path, caption):
    A(f"![{caption}](figures/{os.path.basename(path)})")
    A("")
    A(f"*{caption}*")
    A("")


def add_fig_md(name, caption):
    A(f"[{caption}](../Additional_Output/figures/{name})")
    A("")


def W(w, c):
    return WSI.loc[w, c]


ws_sorted = WS.sort_values("Tmax_degC", ascending=False)
g2 = [g for g in GAUGE_ROWS if g["well"] in VERIFIED]
n_conf = WS.ns_confidence.value_counts()
up = UPFLOW_WELLS
w_up = [w for w in WELLS if WSI.loc[w, "regime"] == "upflow"]
fz_res = FZ.assign(absres=FZ.P_minus_trend_bar.abs())
worst = fz_res.sort_values("absres", ascending=False).iloc[0]
others_max = fz_res[fz_res.well != worst.well].absres.max()
cold_desc = ", ".join(f"{w} (Tmax {fmt(W(w, 'Tmax_degC'), 1)} degC)" for w in COLD_RECHARGE)
hot_open = HY[(HY.band == "open")]
shallow_hot = HY[(HY.band != "open") & (~HY.well.isin(COLD_WELLS))].sort_values("ratio_measured_to_IF97",
                                                                                  ascending=False)
sh0 = shallow_hot.iloc[0]
inj = [w for w in WELLS if WD[w]["kind"] == "injector"]
geo_low = GEO.sort_values("wellhead_minus_surface_m").iloc[0]

A("# Natural-state data interpretation for model calibration")
A("")
A(f"Generated by `Scripts/natural_state_analysis.py` on {dt.date.today().isoformat()} from `Baseline/` only.")
A("")
A("Claim tags: [obs] observation, [interp] interpretation, [assumption] stated assumption, "
  "[reference model] description of the starting model. Units: mRL, bar, degC, kJ/kg, kg/s.")
A("")
# ---------------------------------------------------------------------------------- 1
A("## 1. Key conclusions for calibration")
A("")
A(f"1. **The surveys represent the pre-production natural state, although they are labelled "
  f"{iso(SURVEY_DATES[0])}.** " +
  " ".join(f"{g['well']}: survey {fmt(g['survey_P_at_gauge_bar'], 2)} bar at the gauge ({fmt(g['gauge_elev_mRL'], 1)} mRL) "
           f"against the first reading {fmt(g['first_gauge_P_bar'], 1)} bar on {iso(g['first_gauge_date'])}, a "
           f"difference of {fmt(g['diff'], 2)} bar; by {PCI.loc[g['well'], 'second_gauge_date']} the gauge had already "
           f"fallen by {fmt(PCI.loc[g['well'], 'drop_first_to_second_bar'], 1)} bar." for g in g2) +
  " A survey run on the labelled date would therefore have read several bar lower [obs]. The label is a data "
  "question, not a reason to downgrade [interp].")
A("   **Action:** use every survey as a natural-state target with the weights in `NS_calibration_targets.csv`; "
  "ask the data owner for the true survey dates.")
A(f"2. **The upflow is beneath {' and '.join(up)} in the west of the field, near Fault 1 and Fault 3.** "
  + "; ".join(f"{w} reaches {fmt(W(w, 'Tmax_degC'), 1)} degC at {fmt(W(w, 'Tmax_elev_mRL'), 1)} mRL and "
              f"{fmt(W(w, 'T_at_TD_degC'), 1)} degC at TD ({fmt(W(w, 'TD_mRL'), 1)} mRL)" for w in up) +
  f" [obs]. {', '.join(w_up)} has no inversion and a boiling margin as low as "
  f"{fmt(W(w_up[0], 'min_boiling_margin_degC'), 2)} degC at {fmt(W(w_up[0], 'min_boiling_margin_elev_mRL'), 1)} mRL "
  f"[obs]. " + "; ".join(f"{w} is {fmt(GEOI.loc[w, 'dist_Fault_1_m'], 0)} m from Fault 1 and "
                         f"{fmt(GEOI.loc[w, 'dist_Fault_3_m'], 0)} m from Fault 3" for w in up) +
  " [obs]. That the faults (or their intersection) carry the upflow is an interpretation [interp].")
A("   **Action:** put the deep heat and mass source under columns " +
  ", ".join(f"`{WD[w]['head_col']}` ({w})" for w in up) +
  " and test whether a source along the Fault 1 / Fault 3 intersection reproduces the near-boiling top of "
  f"{w_up[0]}.")
A(f"3. **Hot water flows out laterally at about {fmt(min(W(w, 'inversion_top_mRL') for w in OUTFLOW_WELLS), 0)} to "
  f"{fmt(max(W(w, 'inversion_top_mRL') for w in OUTFLOW_WELLS), 0)} mRL, away from the upflow.** Temperature "
  "inversions start at " + ", ".join(f"{fmt(W(w, 'inversion_top_mRL'), 1)} mRL in {w} ({fmt(W(w, 'inversion_degC'), 1)} degC)"
                                     for w in OUTFLOW_WELLS) +
  " [obs]; Tmax falls with distance from the upflow [obs]. This is a lateral outflow tongue [interp].")
A("   **Action:** give the model a laterally permeable layer at the inversion-top elevations, connecting the upflow to "
  + ", ".join(OUTFLOW_WELLS) + ", with lower vertical permeability below it so that the inversions survive.")
A("4. **Cooler water lies beneath the outflow.** Below the inversion bases the profiles are cool and nearly isothermal, "
  "with low gradients (" + ", ".join(f"{w} {fmt(W(w, 'deep_gradient_degC_per_km'), 1)} degC/km" for w in OUTFLOW_WELLS
                  if np.isfinite(W(w, 'deep_gradient_degC_per_km'))) + f"); {n_iso_pts} open-hole survey points lie in "
  "isothermal sections [obs]. Cool lateral inflow at depth, or internal flow in the open hole, would explain this "
  "[interp]; the data cannot separate the two (a repeat log at a different date, or a static log after a long "
  "shut-in, would).")
A("   **Action:** allow cool lateral inflow at depth below the outflow (side boundary or deep permeable unit); "
  "the isothermal open-hole points already carry the reduced weight 0.5.")
A(f"5. **A conductive cap overlies the hot area.** Shallow gradients are " +
  ", ".join(f"{fmt(W(w, 'shallow_gradient_degC_per_km'), 1)} degC/km in {w}" for w in CAP_WELLS) +
  " [obs]; the cap (shallowest liquid point to 90 % of Tmax) is " +
  ", ".join(f"{fmt(CAPDF.set_index('well').loc[w, 'cap_thickness_m'], 1)} m thick in {w}" for w in CAP_WELLS) +
  " [interp].")
A("   **Action:** use low-permeability cap layers from the surface down to the cap bases in "
  "`Additional_Output/cap_intervals.csv`; the shallow temperatures there constrain the cap conductivity.")
A(f"6. **The margins are cold recharge areas.** {cold_desc} [obs]. The low gradients point to downflow of cold "
  "water [interp].")
A("   **Action:** keep the north-east and south-west margins open to cold recharge (atmospheric top boundary, "
  "permeable near-surface layers) and use these wells' temperatures to tune the recharge.")
A(f"7. **Deep pressures form one connected system.** The feed-zone pressures of {', '.join(VERIFIED)} define "
  f"P = {fmt(TREND_A, 2)} {TREND_B:+.5f} z bar (z in mRL), an implied density of {fmt(TREND_RHO, 1)} kg/m3, RMS "
  f"{fmt(TREND_RMS, 2)} bar [obs]. All other feed zones lie within {fmt(others_max, 2)} bar of it except "
  f"{worst.well} ({fmt(worst.P_minus_trend_bar, 2)} bar at {fmt(worst.mid_mRL, 1)} mRL) [obs].")
A(f"   **Action:** calibrate feed-zone pressures (`weight_P`) as one pressure system; treat the {worst.well} "
  "deficit as either an eastern lateral gradient or early interference, and do not force it.")
A(f"8. **Shallow pressure gradients in the hot wells are steeper than any water column.** The steepest window "
  f"({sh0.well}, {fmt(sh0.z_top_mRL, 1)} to {fmt(sh0.z_bottom_mRL, 1)} mRL) gives a measured/IF97 ratio of "
  f"{fmt(sh0.ratio_measured_to_IF97, 4)} (implied density {fmt(sh0.rho_implied_kgm3, 1)} kg/m3), while open-hole "
  f"windows give a field median of {fmt(hs.loc[('FIELD', 'open'), 'median_ratio'], 4)} [obs]. A datum offset cannot "
  "change a gradient, a convention would act at every depth, two-phase fluid would be lighter and rounding is "
  "below 0.01 [obs]. The excess is a data question [interp].")
A("   **Action:** do not use pressures between feed zones as targets (the targets table already sets "
  "`weight_P` to 0 there); ask whether the shallow pressures were measured or derived.")
A(f"9. **Natural-state confidence varies between wells, so it changes relative weights.** high: "
  f"{', '.join(WS[WS.ns_confidence == 'high'].well)}; medium: {', '.join(WS[WS.ns_confidence == 'medium'].well)}; "
  f"low: {', '.join(WS[WS.ns_confidence == 'low'].well)} [obs].")
A("   **Action:** use `weight` and `weight_P` as written; they already carry the confidence factors.")
A("")
# ---------------------------------------------------------------------------------- 2
A("## 2. The data")
A("")
A("### Inventory")
A("")
A("Every file under `Baseline/` was assigned a role (`../Additional_Output/input_inventory.csv`). "
  "Data read: per-well JSON files, the plain-column `.dat` copies (cross-check only), the well list, the mapping "
  "file, the grid file (geometry sections only) and the well-to-block file. Reference-model files were classified "
  "but not used as evidence; production history workbooks, forecast and production set-ups are out of scope [obs].")
A("")
A(table(INV.assign(path=INV.path.map(lambda x: "`" + x + "`")), ["role", "path", "note"]))
A("")
A("### Checks")
A("")
n_iss = int((INTEGDF.status != "PASS").sum())
A(f"- {len(INTEGDF)} integrity checks, {n_iss} issues (`../Additional_Output/integrity_checks.csv`) [obs]:")
for _, rw in INTEGDF[INTEGDF.status != "PASS"].iterrows():
    A(f"  - {rw.well}, {rw.check}: {rw.detail} [obs]")
A(f"- All {len(DATX)} plain-column files match the JSON values exactly (largest difference "
  f"{fmt(DATX.max_value_diff.max(), 3)}) [obs]. The transient files count time in years from "
  f"{iso(REF_DATE)}, in nominal monthly steps of 1/12 year (largest offset from exact dates "
  f"{fmt(DATX.max_time_diff_yr.max(), 3)} year) [obs]; production rates are negative in the `.dat` files and positive "
  "in the JSON [obs].")
for n_ in DAT_NAME_NOTES:
    A(f"- Irregular file name: {n_} [obs].")
_n_col = int(GEO.column_agrees.sum())
_n_wh = int((GEO.wellhead_mRL == GEO.geometry_file_wellhead_mRL).sum())
A(f"- Survey-point and feed-zone blocks are taken from the well-to-block file; the wellhead column found by "
  f"point-in-polygon agrees with it for {_n_col} of {len(GEO)} wells, and the wellhead (first track point) equals "
  f"the grid file's well track for {_n_wh} of {len(GEO)} wells (`../Additional_Output/well_geometry_check.csv`) [obs].")
A(f"- Every wellhead lies at or below the model surface of its column; {geo_low.well} is "
  f"{fmt(-geo_low.wellhead_minus_surface_m, 1)} m below it [obs]. No survey point lies above the model surface, "
  "and no survey point falls in the top block of its column [obs].")
A(f"- IAPWS-IF97 was verified against the IF97 tables (largest Tsat error {fmt(IF97_MAX_ERR, 4)} K; "
  f"Tsat(1.01325 bar) = {fmt(tsat(1.01325), 2)} degC, Tsat(100 bar) = {fmt(tsat(100.0), 2)} degC) [obs].")
A("")
A("### Provenance labels")
A("")
src = sorted(set(PROVDF[PROVDF.field == "Source"].value) - {""})
A(f"- Non-empty `Source` fields: {', '.join(repr(s_) for s_ in src) or 'none'} (casing design of every well) [obs]. "
  "All survey, transient and flow `Source` and `Comment` fields are empty, so the data carry no provenance beyond "
  "this label [obs].")
A(f"- Every survey carries `StableFlag` = {', '.join(sorted(set(PROVDF[PROVDF.field == 'StableFlag'].value)))} [obs].")
A(f"- Feed-zone quality tags found: {', '.join(TAGS)}. Mapping for `weight_P`: " +
  "; ".join(f"{t} -> {v[0]} ({v[1]})" for t, v in TAG_MAP.items()) + " [obs].")
A("- Track units are labelled `masl` and feed-zone units `m`; both were read as elevations in mRL [obs].")
A("")
A("### Natural-state timing verdict")
A("")
A(f"First production anywhere in the field: {iso(FIELD_FIRST_PROD)}; first injection: {iso(FIELD_FIRST_INJ)}; "
  f"all surveys are labelled {iso(SURVEY_DATES[0])} [obs]. A gauge reading dated on the day production started is "
  "treated as pre-production [assumption]. Reading precision: "
  f"{GAUGE_PRECISION} bar (values are reported to one decimal) [obs]. Pressure trend of gauge-verified feed zones: "
  f"P = {fmt(TREND_A, 2)} {TREND_B:+.5f} z bar, RMS {fmt(TREND_RMS, 2)} bar; tolerance = "
  f"max(3 bar, 2 x RMS = {fmt(2 * TREND_RMS, 2)} bar) = {fmt(TOL, 2)} bar [obs]. Injectate temperature: median "
  f"injectate enthalpy {fmt(WD[inj[0]]['inj_h'], 1)} kJ/kg converted with IF97 Region 1 at "
  f"{fmt(INJ_P_ASSUMED, 1)} bar [assumption] gives {fmt(WD[inj[0]]['inj_T'], 1)} degC; an isothermal section within "
  "10 degC of this would mark injection [assumption].")
A("")
pc_tab = PC[["well", "well_kind", "gauge_elev_mRL", "diff_survey_minus_first_bar",
             "max_feedzone_residual_from_trend_bar", "verdict", "confidence", "reason"]]
A(table(pc_tab, list(pc_tab.columns), ["well", "kind", "gauge (mRL)", "survey - first reading (bar)",
                                       "max feed-zone residual (bar)", "verdict", "confidence", "reason"],
        {"gauge_elev_mRL": 1, "diff_survey_minus_first_bar": 2, "max_feedzone_residual_from_trend_bar": 2}))
A("")
for g in GAUGE_ROWS:
    if g["pre"] and abs(g["diff"]) > GAUGE_PRECISION:
        A(f"- {g['well']}: the survey matches the first gauge reading at {fmt(g['z_match'], 1)} mRL, not at the stated "
          f"gauge elevation {fmt(g['gauge_elev_mRL'], 1)} mRL. If the gauge depth were along-hole (wellhead minus "
          f"measured depth), its elevation would be {fmt(g['z_md'], 1)} mRL, where the survey reads "
          f"{fmt(g['p_md'], 2)} bar against {fmt(g['first_gauge_P_bar'], 1)} bar [obs]. A gauge depth recorded "
          "along hole in this deviated well would explain the miss [interp]; the verdict stays at the rubric value.")
    elif not g["pre"]:
        A(f"- {g['well']}: the gauge series starts {iso(g['first_gauge_date'])}, after field production started "
          f"({iso(FIELD_FIRST_PROD)}), so it cannot test the survey [obs].")
A(f"- {worst.well}: the deeper feed-zone pressures are below the trend by up to {fmt(worst.absres, 2)} bar "
  f"(tolerance {fmt(TOL, 2)} bar), so the rubric gives '{NS[worst.well]['verdict']}' [obs]. The well has no flow "
  "record; interference from field production or a natural east-west pressure gradient could each explain it "
  "[interp]. A gauge reading in this well before production would decide.")
A("")
fig_md(os.path.join(CALFIG, "gauge_check.png"),
       "F7 Gauge check: transient pressure against date for each gauged well, with the survey pressure at the gauge "
       "elevation, the first-reading precision band and the production start dates; claim type obs.")
# ---------------------------------------------------------------------------------- 3
A("## 3. Temperature structure")
A("")
A("Per-well characteristics (definitions in CLAUDE.md 4.4; `well_summary.csv`) [obs]; the regime is rule-based "
  "[interp]:")
A("")
tcols = ["well", "Tmax_degC", "Tmax_elev_mRL", "T_at_TD_degC", "TD_mRL", "inversion_degC", "inversion_top_mRL",
         "inversion_base_mRL", "shallow_gradient_degC_per_km", "deep_gradient_degC_per_km", "z_100C_mRL",
         "z_200C_mRL", "min_boiling_margin_degC", "regime"]
A(table(WS, tcols, ["well", "Tmax (degC)", "at (mRL)", "T at TD (degC)", "TD (mRL)", "inversion (degC)",
                    "inv. top (mRL)", "inv. base (mRL)", "shallow grad. (degC/km)", "deep grad. (degC/km)",
                    "100 degC (mRL)", "200 degC (mRL)", "min boiling margin (degC)", "regime"],
        {c: (2 if c == "min_boiling_margin_degC" else 1) for c in tcols[1:-1]}))
A("")
for _, rw in WS.iterrows():
    if rw.regime_note:
        A(f"- {rw.well} ({rw.regime}): {rw.regime_note} [interp].")
A("")
A("Wellbore effects [interp]: the inversions start below the production casing shoe in every outflow well, "
  "within the open hole; internal flow between feed zones could flatten the open-hole profile (hence the reduced "
  "weight on isothermal sections), but it cannot create an inversion of tens of degC sustained over hundreds of "
  "metres below the deepest feed zone, so the inversions are taken as formation features.")
A("")
fig_md(os.path.join(CALFIG, f"profile_panels_{safe(list(GROUPS)[0])}.png"),
       f"F2 Temperature profiles of all wells on identical axes ({fmt(Z_TOP, 0)} to {fmt(Z_BOT, 0)} mRL, 0 to "
       f"{fmt(T_TOP, 0)} degC), with Tsat(P), shoe, feed zones, surface and near-boiling intervals; claim type obs.")
A("Individual well profiles (F2, temperature and pressure): " +
  ", ".join(f"[{w}](figures/profile_{w}.png)" for w in WELLS) + ".")
A("")
fig_md(os.path.join(CALFIG, "maps_temperature_bands.png"),
       f"F3 Band-mean temperature minus the field median in 500 m bands from {BANDS[0][0]:.0f} to "
       f"{BANDS[-1][1]:.0f} mRL, at each track's mid-band position; claim type obs.")
for ln in LINES:
    fig_md(os.path.join(CALFIG, f"section_{ln['name']}.png"),
           f"F4 Temperature section {ln['name']}-{ln['name']}' from ({fmt(ln['a'][0], 0)}, {fmt(ln['a'][1], 0)}) to "
           f"({fmt(ln['b'][0], 0)}, {fmt(ln['b'][1], 0)}), {ln['basis']}; colour only inside the data envelope; "
           f"vertical exaggeration {fmt(ln['ve'][('T', False)], 1)}; claim type mixed.")
fig_md(os.path.join(CALFIG, "PT_state.png"),
       "F5 P-T state of every survey point against the IF97 saturation curve; claim type obs.")
A("Interpolated companions (whole-field picture, not data; no value is quoted from them) [interp]:")
A("")
for ln in LINES:
    add_fig_md(f"section_{ln['name']}_interpolated.png",
               f"[interp] Temperature on line {ln['name']}-{ln['name']}', linear Delaunay interpolation")
add_fig_md("maps_temperature_interpolated.png", "[interp] Band temperature deviation, linear Delaunay interpolation")
# ---------------------------------------------------------------------------------- 4
A("## 4. Pressure structure")
A("")
A("Pressure at common elevations (interpolation within each survey only; "
  "`../Additional_Output/T_P_at_common_elevations.csv`) [obs]:")
A("")
ce_sel = CE[CE.elevation_mRL.isin([-250.0, -500.0, -1000.0, -1500.0, -2000.0])]
piv = ce_sel.pivot(index="well", columns="elevation_mRL", values="P_bar").reindex(WELLS)
dev = ce_sel.pivot(index="well", columns="elevation_mRL", values="P_minus_trend_bar").reindex(WELLS)
zs_ = sorted(piv.columns, reverse=True)
rows_ = []
for w in WELLS:
    rw = {"well": w}
    for z in zs_:
        v, d = piv.loc[w, z], dev.loc[w, z]
        rw[z] = "" if not np.isfinite(v) else f"{fmt(v, 1)} ({d:+.1f})"
    rows_.append(rw)
A(table(pd.DataFrame(rows_), ["well"] + zs_, ["well"] + [f"P at {z:.0f} mRL, bar (minus trend)" for z in zs_]))
A("")
_d1000 = dev[-1000.0].abs()
_d2000 = dev[-2000.0]
A(f"At -1000 mRL every well with data is within {fmt(_d1000.max(), 1)} bar of the gauge-verified trend "
  f"({_d1000.idxmax()} the largest) [obs]. The deviations grow above and below that level: the trend is a straight "
  f"line through hot-well feed zones between {fmt(TREND_ZRANGE[0], 1)} and {fmt(TREND_ZRANGE[1], 1)} mRL, while each "
  f"well column follows its own temperature, so the cold {_d2000.idxmax()} is {fmt(_d2000.max(), 1)} bar above the "
  "trend at -2000 mRL and the shallow wellbore pressures lie far below its upward extrapolation [obs]. Only "
  "feed-zone pressures are used as targets [interp].")
A("")
A("Water levels from the three shallowest liquid-filled points [obs]: " +
  ", ".join(f"{w} {fmt(W(w, 'water_level_mRL'), 1)} mRL" for w in WELLS) +
  f". {', '.join(w for w in WELLS if WSI.loc[w, 'water_level_extrapolation_m'] > 50)} needs an extrapolation of more "
  "than 50 m above its top survey point (`water_level_extrapolation_m`), so treat its water level as uncertain "
  "[interp].")
A("")
A("Hydrostatic check (CLAUDE.md 4.3; `../Additional_Output/hydrostatic_check.csv`, "
  "`hydrostatic_summary.csv`) [obs]:")
A("")
hs_tab = HS[HS.band != "all"]
A(table(hs_tab, ["well", "band", "n_windows", "median_ratio", "min_ratio", "max_ratio", "excess_over_cold_component"],
        ["well", "band", "windows", "median ratio measured/IF97", "min", "max", "excess over cold-well median"],
        {"median_ratio": 4, "min_ratio": 4, "max_ratio": 4, "excess_over_cold_component": 4}))
A("")
n_sp = int((HY.band.fillna("") == "").sum())
A(f"- Cold wells (Tmax < 100 degC): {', '.join(COLD_WELLS)}. Their median ratio is "
  f"{fmt(cold_med['open'], 4)} in open hole and {fmt(cold_med['all'], 4)} over all windows; no cold-well window "
  "lies wholly behind casing, so the cased-band component is undefined [obs]. This field-wide component is "
  "attributed to a datum, gravity or density convention, not to physical overpressure [interp].")
A(f"- {n_sp} windows span a casing shoe (the window rule closes only at a feature once 300 m long, and a short "
  "remnant joins the window above); they have no band and count only in the overall statistics [obs].")
A(f"- Open-hole windows range from {fmt(hs.loc[('FIELD', 'open'), 'min_ratio'], 4)} to "
  f"{fmt(hs.loc[('FIELD', 'open'), 'max_ratio'], 4)} (field median {fmt(hs.loc[('FIELD', 'open'), 'median_ratio'], 4)}); "
  f"{hs.xs('open', level='band').drop('FIELD').median_ratio.idxmax()} has the largest open-hole excess over the "
  "cold-well component. The windows that start in the shallow part of the hot wells reach "
  f"{fmt(shallow_hot.ratio_measured_to_IF97.max(), 4)} [obs]. Checks: a datum offset cannot change a gradient; a "
  "convention would act at every depth and in the cold wells too; two-phase fluid would lower the ratio; the "
  f"rounding uncertainty of the ratio is at most {fmt(HY.ratio_rounding_uncertainty.max(), 4)} [obs]. No wellbore "
  "explanation is left; dynamic overpressure from upflow through the cap would raise formation gradients but "
  "cannot be measured inside a static well column, so the shallow pressures are a data question [interp].")
A("")
fig_md(os.path.join(CALFIG, "maps_pressure_bands.png"),
       f"F3 Band-mean pressure minus the trend of gauge-verified feed zones in 500 m bands from {BANDS[0][0]:.0f} to "
       f"{BANDS[-1][1]:.0f} mRL; claim type mixed (the trend is derived).")
for ln in LINES:
    fig_md(os.path.join(CALFIG, f"section_{ln['name']}_P.png"),
           f"F4 Pressure-minus-trend section {ln['name']}-{ln['name']}' (same line and ends as the temperature "
           f"section), vertical exaggeration {fmt(ln['ve'][('P', False)], 1)}; claim type mixed.")
A("Interpolated companions [interp]:")
A("")
for ln in LINES:
    add_fig_md(f"section_{ln['name']}_P_interpolated.png",
               f"[interp] Pressure minus trend on line {ln['name']}-{ln['name']}', linear Delaunay interpolation")
add_fig_md("maps_pressure_interpolated.png", "[interp] Band pressure deviation, linear Delaunay interpolation")
# ---------------------------------------------------------------------------------- 5
A("## 5. Conceptual model")
A("")
fig_md(os.path.join(CALFIG, "conceptual_model.png"),
       f"F6 Conceptual model from the data alone: plan view and section {F6_LINE['name']}-{F6_LINE['name']}' "
       "(the F4 line with the most wells); numbers key to the table below; solid = measured feature, dashed = "
       "inferred; claim type interp.")
A(f"The same elements on the other section lines: " +
  ", ".join(f"[line {l['name']}](../Additional_Output/figures/conceptual_model_section_{l['name']}.png)"
            for l in LINES if l is not F6_LINE) + " [interp].")
A("")
A(table(CMDF.assign(claim_type=CMDF.claim_type.map(lambda c: f"[{c}]")),
        ["no", "element", "location", "supporting_observation", "claim_type"],
        ["no", "element", "location", "supporting observation", "claim type"]))
A("")
A("Relation to mapped features [interp]: the upflow wells sit close to Fault 1 and Fault 3 (distances in "
  "`../Additional_Output/well_geometry_check.csv`), the outflow runs south-east and east across Fault 2 towards "
  f"the river corridor, and the conductive eastern well {', '.join(w for w in WELLS if WSI.loc[w, 'regime'] == 'conductive')} "
  "lies next to the river. Whether Fault 2 is a barrier cannot be decided from these data: the outflow wells on "
  "either side of it share one pressure trend, which argues against a strong barrier; a pressure-interference test "
  "across Fault 2 would decide.")
A("")
# ---------------------------------------------------------------------------------- 6
A("## 6. Recommended calibration set-up")
A("")
A("**Targets and weights.**")
A("")
A(f"- Temperature: all {len(TGT)} survey points in `NS_calibration_targets.csv`, weight = location factor x "
  "natural-state confidence factor (defaults of CLAUDE.md 5.5, no departure) [obs].")
A(f"- Pressure: the {len(FZ)} feed zones in `feedzones.csv` (`weight_P` = quality factor x confidence factor); "
  "`weight_P` is zero away from feed zones [obs].")
A("- Layer-centre values for direct block comparison: `T_P_at_layer_centres.csv` (interpolation only) [obs].")
A("")
wsum = TGT.groupby("weight_reason").size().reset_index(name="points")
A(table(wsum, ["weight_reason", "points"], ["weight rule", "points"]))
A("")
A("**Where heat and mass enter and leave** [interp]:")
A("")
A("- Enter: deep upflow beneath " + ", ".join(f"{w} (column `{WD[w]['head_col']}`)" for w in up) +
  "; cold recharge from the surface near " + ", ".join(f"{w} (column `{WD[w]['head_col']}`)" for w in COLD_RECHARGE) +
  "; cool lateral inflow at depth beneath the outflow.")
A("- Leave: lateral outflow towards " + ", ".join(OUTFLOW_WELLS) + " and beyond; conductive loss through the cap.")
A("")
A("**Permeability structure the data imply** [interp]: a permeable vertical path at the upflow; a laterally "
  "permeable outflow layer near the inversion tops; a low-permeability cap above it; lower vertical permeability "
  "between the outflow layer and the cooler water below; permeable feed zones at the elevations in `feedzones.csv`.")
A("")
A("**Boundary requirements** [interp]: an atmospheric top boundary that admits cold recharge at the margins; "
  "a basal heat and mass source confined to the upflow area; side boundaries that can supply cool water at depth "
  "and accept the lateral outflow; a pressure datum consistent with the feed-zone trend (implied density "
  f"{fmt(TREND_RHO, 1)} kg/m3).")
A("")
A("[reference model] The starting model's rock types, boundary conditions and results were catalogued in the "
  "inventory but not read; no property of the reference model is used above.")
A("")
# ---------------------------------------------------------------------------------- 7
A("## 7. Data quality issues and open questions")
A("")
A("Ranked by effect on calibration:")
A("")
A(f"1. Survey dates: every survey is labelled {iso(SURVEY_DATES[0])}, after production started, but the gauge test "
  "shows pre-production pressures [obs]. Confirm the true dates.")
A(f"2. {worst.well} feed-zone pressures sit {fmt(worst.absres, 2)} bar below the verified trend [obs]; interference or a "
  "natural gradient [interp]. This lowers its confidence and weight.")
A("3. Shallow pressure gradients in the hot wells exceed any water column (section 4) [obs]; confirm whether "
  "shallow pressures were measured or derived.")
for g in GAUGE_ROWS:
    if g["pre"] and abs(g["diff"]) > GAUGE_PRECISION:
        A(f"4. {g['well']} gauge elevation: the first reading matches the survey {fmt(g['z_match'] - g['gauge_elev_mRL'], 1)} m "
          "above the stated gauge elevation, consistent with an along-hole depth [interp]. Confirm the gauge elevation.")
A(f"5. Field-wide hydrostatic component: cold-well ratio {fmt(cold_med['open'], 4)} (open hole) to "
  f"{fmt(cold_med['all'], 4)} (all windows) [obs]; check the pressure datum and units with the data owner.")
for k_, (_, rw) in enumerate(INTEGDF[INTEGDF.status != "PASS"].iterrows()):
    A(f"{6 + k_}. {rw.well}: {rw.check} - {rw.detail} [obs].")
A(f"1. Wellheads lie below the model surface by up to {fmt(-geo_low.wellhead_minus_surface_m, 1)} m ({geo_low.well}) "
  "[obs]; the model topography is smoother than the wellhead elevations [interp].")
inj_hot = [w for w in WELLS if NS[w]["inj_hot"]]
if inj_hot:
    A(f"1. {', '.join(inj_hot)}: injector whose feed zones are hotter than the injectate, so injection cooling could be "
      "hidden; low confidence by rule [obs].")
A(f"1. Feed-zone quality tag 'Synthetic' ({', '.join(FZ[FZ.quality == 'Synthetic'].well)}) is not a quality grade and "
  "gets the lowest pressure factor [obs].")
A("1. All provenance fields except the casing-design source are empty; the casing-design source reads "
  "'Synthetic data' [obs].")
A("")
# ---------------------------------------------------------------------------------- 8
A("## 8. Output file index")
A("")
A("`Output/Calibration_Input/`:")
A("")
for f in sorted(set(os.listdir(CAL)) | {"NS_data_interpretation.md"}):
    if f != "figures":
        A(f"- `{f}`")
A("- `figures/`: " + ", ".join(f"`{os.path.basename(p_)}`" for p_ in FIGS_CAL))
A("")
A("`Output/Additional_Output/`:")
A("")
add_files = sorted(set(os.listdir(ADD)) | {"run_console_log.txt", "acceptance_check.txt", "process_log.md"})
for f in add_files:
    if f != "figures":
        A(f"- `{f}`")
A("- `figures/`: " + ", ".join(f"`{os.path.basename(p_)}`" for p_ in FIGS_ADD))
A("")
REPORT = "\n".join(L) + "\n"
with open(os.path.join(CAL, "NS_data_interpretation.md"), "w", encoding="utf-8", newline="\n") as f:
    f.write(REPORT)
log(f"wrote {rel(os.path.join(CAL, 'NS_data_interpretation.md'))}")
with open(os.path.join(ADD, "run_console_log.txt"), "w", encoding="utf-8", newline="\n") as f:
    f.write(_LOG.getvalue())
