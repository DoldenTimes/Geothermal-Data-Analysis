"""Natural-state data interpretation for geothermal model calibration (CLAUDE.md v4).

Run from the repository root:

    python Scripts/natural_state_analysis.py

Reads only Baseline/ and regenerates everything in Output/Calibration_Input/ and
Output/Additional_Output/ (except process_log.md and acceptance_check.txt): the input
inventory, every table, every figure and the report NS_data_interpretation.md.
"""
import ast
import datetime as dt
import fnmatch
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
from matplotlib.path import Path as MplPath  # noqa: E402
from matplotlib.collections import LineCollection, PolyCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# ------------------------------------------------------------------------------------------
# Paths and console log
# ------------------------------------------------------------------------------------------
ROOT = os.getcwd()
BASE = os.path.join(ROOT, "Baseline")
CAL = os.path.join(ROOT, "Output", "Calibration_Input")
ADD = os.path.join(ROOT, "Output", "Additional_Output")
CALFIG = os.path.join(CAL, "figures")
ADDFIG = os.path.join(ADD, "figures")
for d in (CAL, ADD, CALFIG, ADDFIG):
    os.makedirs(d, exist_ok=True)
for d in (CALFIG,):   # standard figure folder holds only this run's standard figures
    for f in os.listdir(d):
        if f.endswith(".png"):
            os.remove(os.path.join(d, f))

_CONSOLE = io.StringIO()


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    _CONSOLE.write(s + "\n")


def rel(p):
    return os.path.relpath(p, ROOT).replace("\\", "/")


def write_csv(df, path):
    df.to_csv(path, index=False, lineterminator="\n", encoding="utf-8")
    log(f"wrote {rel(path)} ({len(df)} rows)")


def write_text(text, path):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    log(f"wrote {rel(path)}")


G = 9.80665
P_ATM = 1.01325
P_CRIT = 220.64
T_CRIT = 373.946

# ------------------------------------------------------------------------------------------
# IAPWS-IF97: Region 4 (saturation) and Region 1 (compressed liquid)
# ------------------------------------------------------------------------------------------
_N4 = [None, 0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5,
       -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6,
       -0.23855557567849, 0.65017534844798e3]


def tsat(p_bar):
    """Saturation temperature (degC) from pressure (bar), IF97 eq. 31. NaN above Pcrit."""
    p = p_bar / 10.0
    if not (611.213e-6 <= p <= 22.064):
        return float("nan")
    n = _N4
    b = p ** 0.25
    e = b * b + n[3] * b + n[6]
    f = n[1] * b * b + n[4] * b + n[7]
    g = n[2] * b * b + n[5] * b + n[8]
    d = 2 * g / (-f - math.sqrt(f * f - 4 * e * g))
    return (n[10] + d - math.sqrt((n[10] + d) ** 2 - 4 * (n[9] + n[10] * d))) / 2 - 273.15


def psat(t_c):
    """Saturation pressure (bar) from temperature (degC), IF97 eq. 30."""
    n = _N4
    T = t_c + 273.15
    th = T + n[9] / (T - n[10])
    A = th * th + n[1] * th + n[2]
    B = n[3] * th * th + n[4] * th + n[5]
    C = n[6] * th * th + n[7] * th + n[8]
    return (2 * C / (-B + math.sqrt(B * B - 4 * A * C))) ** 4 * 10.0


_R1 = [(0, -2, 0.14632971213167), (0, -1, -0.84548187169114), (0, 0, -0.37563603672040e1),
       (0, 1, 0.33855169168385e1), (0, 2, -0.95791963387872), (0, 3, 0.15772038513228),
       (0, 4, -0.16616417199501e-1), (0, 5, 0.81214629983568e-3), (1, -9, 0.28319080123804e-3),
       (1, -7, -0.60706301565874e-3), (1, -1, -0.18990068218419e-1), (1, 0, -0.32529748770505e-1),
       (1, 1, -0.21841717175414e-1), (1, 3, -0.52838357969930e-4), (2, -3, -0.47184321073267e-3),
       (2, 0, -0.30001780793026e-3), (2, 1, 0.47661393906987e-4), (2, 3, -0.44141845330846e-5),
       (2, 17, -0.72694996297594e-15), (3, -4, -0.31679644845054e-4), (3, 0, -0.28270797985312e-5),
       (3, 6, -0.85205128120103e-9), (4, -5, -0.22425281908000e-5), (4, -2, -0.65171222895601e-6),
       (4, 10, -0.14341729937924e-12), (5, -8, -0.40516996860117e-6), (8, -11, -0.12734301741641e-8),
       (8, -6, -0.17424871230634e-9), (21, -29, -0.68762131295531e-18), (23, -31, 0.14478307828521e-19),
       (29, -38, 0.26335781662795e-22), (30, -39, -0.11947622640071e-22), (31, -40, 0.18228094581404e-23),
       (32, -41, -0.93537087292458e-25)]
_R = 0.461526  # kJ/(kg K)


def _r1_derivs(t_c, p_bar):
    T = t_c + 273.15
    p = p_bar / 10.0
    pi, tau = p / 16.53, 1386.0 / T
    gp = gt = 0.0
    for I, J, n in _R1:
        gp += -n * I * (7.1 - pi) ** (I - 1) * (tau - 1.222) ** J
        gt += n * (7.1 - pi) ** I * J * (tau - 1.222) ** (J - 1)
    return T, p, pi, tau, gp, gt


def rho_liq(t_c, p_bar):
    """Liquid density (kg/m3), IF97 Region 1."""
    T, p, pi, tau, gp, _ = _r1_derivs(t_c, p_bar)
    v = _R * T / (p * 1000.0) * pi * gp   # m3/kg (R in kJ/kgK, p in MPa -> kPa)
    return 1.0 / v


def h_liq(t_c, p_bar):
    """Liquid specific enthalpy (kJ/kg), IF97 Region 1."""
    T, p, pi, tau, gp, gt = _r1_derivs(t_c, p_bar)
    return _R * T * tau * gt


def validate_if97():
    errs = []
    # Region 4 verification (IF97 Tables 35 and 36)
    for p_mpa, tk in [(0.1, 372.755919), (1.0, 453.035632), (10.0, 584.149488)]:
        errs.append(("Tsat", p_mpa, abs(tsat(p_mpa * 10) + 273.15 - tk)))
    for tk, p_mpa in [(300, 0.353658941e-2), (500, 0.263889776e1), (600, 0.123443146e2)]:
        errs.append(("Psat", tk, abs(psat(tk - 273.15) / 10 - p_mpa) / p_mpa))
    # Region 1 verification (IF97 Table 5): v and h
    for tk, p_mpa, v, h in [(300, 3, 0.100215168e-2, 0.115331273e3), (300, 80, 0.971180894e-3, 0.184142828e3),
                            (500, 3, 0.120241800e-2, 0.975542239e3)]:
        errs.append(("v1", tk, abs(1 / rho_liq(tk - 273.15, p_mpa * 10) - v) / v))
        errs.append(("h1", tk, abs(h_liq(tk - 273.15, p_mpa * 10) - h) / h))
    ref = [("Tsat(1.01325 bar)", tsat(1.01325), 99.97), ("Tsat(100 bar)", tsat(100.0), 311.0)]
    bad = [e for e in errs if (e[0] == "Tsat" and e[2] > 1e-4) or (e[0] != "Tsat" and e[2] > 1e-7)]
    bad += [r for r in ref if abs(r[1] - r[2]) > 0.1]
    for r in ref:
        log(f"IF97 check: {r[0]} = {r[1]:.3f} degC (reference {r[2]})")
    if bad:
        raise SystemExit(f"IF97 implementation failed verification: {bad}")
    log("IF97 Region 1 and Region 4 pass the verification tables (Tsat error < 1e-4 K)")


validate_if97()

# ------------------------------------------------------------------------------------------
# 1. Input inventory (CLAUDE.md 2.1)
# ------------------------------------------------------------------------------------------
INV_RULES = [  # (glob relative to Baseline/, role, note); first match wins; directories end with /
    ("*/Data/Data_json/", "well_data_json",
     "Per-well JSON: track, casing/liner design, feed zones, downhole T and P surveys, transient pressure, "
     "production/injection flow and enthalpy. Primary data source"),
    ("*/Data/Area_data.json", "well_data_json",
     "Field-total steam flow and power output (production record). Read only for the date of first field "
     "production (4.1 timing test)"),
    ("*/Data/Data_dat/", "well_data_dat",
     "Plain-column copies of the well data; cross-check only. Irregular file prefix 'Wel10' for WELL-10"),
    ("*/Data/Groups_data.json", "well_groups", "Well list and grouping (one group 'All'); defines well order"),
    ("*/Data/Mapping/", "mapping", "Faults (polylines), lake (filled polygon) and rivers (ragged multi-segment)"),
    ("*/Data/Plotting/", "reference_other",
     "Plotting configuration with hard-coded paths from another machine; not read"),
    ("*/Data/Make_up_json/", "out_of_scope", "Empty directory"),
    ("*/Models/*/gTE14515.dat", "geometry_dat",
     "MULgraph grid: VERTICES, GRID, LAYERS, SURFA, WELLS sections used (geometry only)"),
    ("*/Models/*/gTE14515_wells.json", "well_blocks_json",
     "Well-to-block intersections, feed-zone blocks and interpolation weights; authority for block assignment"),
    ("*/Models/*/gTE14515_boundary_conditions.json", "reference_bc",
     "Reference model boundary conditions (atmosphere, rain, basal heat flux, lake columns); reference only"),
    ("*/Models/*/gTE14515_dict_rocktype.json", "reference_rock",
     "Reference model rock-type legend (lithology, faults, zones); reference only"),
    ("*/Models/*/gTE14515.auxjson", "reference_other",
     "Viewer scene file with hard-coded source paths; not read"),
    ("*/Models/*/currentGeomSetup.json", "reference_other", "Model set-up and plotting configuration; not read"),
    ("*/Models/*/NaturalState/TE14515_NS_000.dat", "reference_other",
     "Reference natural-state simulator input (rocks, generators, parameters); reference only, not evidence"),
    ("*/Models/*/NaturalState/TE14515_NS_000.incon", "reference_results",
     "Reference model initial conditions; reference only, not read"),
    ("*/Models/*/NaturalState/currentNSSetup.json", "reference_other", "Natural-state run configuration; not read"),
    ("*/Models/*/NaturalState/Reports/", "out_of_scope", "Empty directory"),
    ("*/Models/*/Production/", "out_of_scope", "Production model set-up; out of scope"),
    ("*/Models/*/FutureScenarios/", "out_of_scope", "Future-scenario set-up; out of scope"),
    ("*/Parameters/ProductionHistory/", "out_of_scope",
     "Production, reinjection and monitoring spreadsheets for the production model; not read"),
]


def build_inventory():
    rows, covered, unclassified = [], set(), []
    allfiles = []
    for dp, dns, fns in os.walk(BASE):
        dns.sort()
        for f in sorted(fns):
            allfiles.append(os.path.relpath(os.path.join(dp, f), BASE).replace("\\", "/"))
    alldirs = []
    for dp, dns, _ in os.walk(BASE):
        for d in sorted(dns):
            alldirs.append(os.path.relpath(os.path.join(dp, d), BASE).replace("\\", "/") + "/")
    for pat, role, note in INV_RULES:
        targets = alldirs if pat.endswith("/") else allfiles
        for t in sorted(targets):
            if fnmatch.fnmatch(t, pat):
                rows.append((role, "Baseline/" + t.rstrip("/"), note))
                covered.add(t)
    for f in allfiles:
        if not any(f == c or (c.endswith("/") and f.startswith(c)) for c in covered):
            unclassified.append(f)
            rows.append(("out_of_scope", "Baseline/" + f, "Unclassified by the inventory rules; not read"))
    inv = pd.DataFrame(rows, columns=["role", "path", "note"])
    if unclassified:
        log("WARNING: files not matched by any inventory rule:", unclassified)
    empty = [d for d in alldirs if not any(f.startswith(d) for f in allfiles)]
    return inv, allfiles, empty


INV, ALLFILES, EMPTY_DIRS = build_inventory()
write_csv(INV, os.path.join(ADD, "input_inventory.csv"))
log(f"inventory: {len(ALLFILES)} files under Baseline, {len(INV)} inventory rows; empty directories: {EMPTY_DIRS}")


def inv_paths(role):
    return [os.path.join(ROOT, p) for p in INV[INV.role == role].path]


GEOM_DAT = inv_paths("geometry_dat")[0]
WBLOCKS = inv_paths("well_blocks_json")[0]
JSON_DIR = [p for p in inv_paths("well_data_json") if os.path.isdir(p)][0]
AREA_JSON = [p for p in inv_paths("well_data_json") if os.path.isfile(p)]
DAT_DIR = inv_paths("well_data_dat")[0]
GROUPS = inv_paths("well_groups")[0]
MAPPING = [os.path.join(p, f) for p in inv_paths("mapping") for f in sorted(os.listdir(p)) if f.endswith(".json")][0]


# ------------------------------------------------------------------------------------------
# 2. Load well data
# ------------------------------------------------------------------------------------------
def lit(v):
    return ast.literal_eval(v) if isinstance(v, str) else v


def ddmmyyyy(s):
    return dt.datetime.strptime(s, "%d/%m/%Y").date()


groups = json.load(open(GROUPS))
GROUP_WELLS = {}
for area, gs in groups["Areas"].items():
    for gname, wl in gs.items():
        GROUP_WELLS[gname] = wl
WELLS = [w for gname in GROUP_WELLS for w in GROUP_WELLS[gname]]
WELLS = list(dict.fromkeys(WELLS))
log("wells (well_groups order):", WELLS)

W = {}
PROVENANCE = []
for w in WELLS:
    d = json.load(open(os.path.join(JSON_DIR, f"{w}_data.json")))
    rec = {"name": w, "aut2": d.get("AUT2_Name", "")}
    rec["track"] = np.array([lit(r) for r in d["Well_Track"]["Data"]], float)
    rec["design"] = [lit(r) for r in d["Well_Design"]["Data"]]
    rec["fz"] = [lit(r) for r in d["Feed_Zones"]["Data"]]
    tsv = d["Downhole_Temperature"]
    psv = d["Downhole_Pressure"]
    if len(tsv) != 1 or len(psv) != 1:
        log(f"NOTE {w}: {len(tsv)} T and {len(psv)} P survey dates; the first of each is used")
    tdate, tdat = next(iter(tsv.items()))
    pdate, pdat = next(iter(psv.items()))
    rec["survey_date"] = ddmmyyyy(tdate)
    rec["survey_date_P"] = ddmmyyyy(pdate)
    T = {float(k): float(v) for k, v in tdat["Data"].items()}
    P = {float(k): float(v) for k, v in pdat["Data"].items()}
    zs = sorted(set(T) | set(P), reverse=True)
    rec["survey"] = pd.DataFrame({"z": zs, "T": [T.get(z, np.nan) for z in zs], "P": [P.get(z, np.nan) for z in zs]})
    for key, blk in (("Downhole_Temperature", tdat), ("Downhole_Pressure", pdat)):
        PROVENANCE.append((w, key, blk.get("Source", ""), blk.get("Comment", ""), blk.get("StableFlag", "")))
    PROVENANCE.append((w, "Well_Design", d["Well_Design"].get("Source", ""), "", ""))
    rec["gauges"] = {}
    for k, v in d.get("Transient_Pressure", {}).items():
        if isinstance(v, dict) and "Data" in v:
            rec["gauges"][float(k)] = {ddmmyyyy(a): float(b) for a, b in v["Data"].items()}
            PROVENANCE.append((w, f"Transient_Pressure {k}", v.get("Source", ""), v.get("Comment", ""), ""))
    rec["flows"] = {}
    for key in ("Production_Mass_Flow", "Injection_Mass_Flow", "Transient_Enthalpy", "Injection_Enthalpy"):
        if key in d:
            rec["flows"][key] = {ddmmyyyy(a): float(b) for a, b in d[key]["Data"].items()}
            PROVENANCE.append((w, key, d[key].get("Source", ""), d[key].get("Comment", ""), ""))
    rec["raw_json"] = d
    W[w] = rec

# derived well geometry
for w, r in W.items():
    tr = r["track"]
    r["wellhead"] = tr[0, 2]
    casings = [e for e in r["design"] if not e[3]]
    liners = [e for e in r["design"] if e[3]]
    r["shoe"] = min(e[1] for e in casings) if casings else r["wellhead"]
    r["slotted"] = [(max(e[0], e[1]), min(e[0], e[1])) for e in liners]
    r["casing_spec"] = "; ".join(f"{e[4]} {e[0]} to {e[1]} mRL" for e in r["design"])
    r["td_track"] = tr[-1, 2]
    r["fz"] = [(float(max(a, b)), float(min(a, b)), str(q)) for a, b, q in r["fz"]]
    r["injector"] = "Injection_Mass_Flow" in r["flows"]
    r["producer"] = "Production_Mass_Flow" in r["flows"]


def track_xy(w, z):
    tr = W[w]["track"]
    zz = tr[::-1, 2]
    return float(np.interp(z, zz, tr[::-1, 0])), float(np.interp(z, zz, tr[::-1, 1]))


def interp_profile(zs, vs, z):
    zs, vs = np.asarray(zs, float), np.asarray(vs, float)
    o = np.argsort(zs)
    return float(np.interp(z, zs[o], vs[o]))


# ------------------------------------------------------------------------------------------
# 3. Geometry: MULgraph grid and the well-to-block file
# ------------------------------------------------------------------------------------------
GEO_LINES = open(GEOM_DAT).read().split("\n")
GEO_HEADERS = {"VERTICES", "GRID", "CONNECTIONS", "LAYERS", "SURFA", "WELLS"}


def section_lines(header):
    out, on = [], False
    for line in GEO_LINES:
        if not on:
            if line.startswith(header):
                on = True
            continue
        if line.strip() == "" or line.strip() in GEO_HEADERS:
            break
        out.append(line)
    return out


LAYERS = {}
for line in section_lines("LAYERS"):
    p = line.split()
    LAYERS[int(p[0])] = (float(p[1]), float(p[2]))
LAYER = {i: {"top": LAYERS[i - 1][0], "bottom": LAYERS[i][0], "centre": LAYERS[i][1]} for i in sorted(LAYERS) if i > 0}
GRID_TOP = LAYERS[0][0]
SURF = {}
for line in section_lines("SURFA"):
    p = line.split()
    SURF[p[0]] = float(p[1])
VERT = {}
for line in section_lines("VERTICES"):
    # vertex lines carry no whitespace between name and coordinates (e.g. 'aab2812872.706327567.88'),
    # so whitespace parsing is impossible here: the name is the leading letters and each coordinate
    # has exactly two decimals (checked: every vertex line matches this pattern)
    m = re.fullmatch(r"([a-z]+)\s*(-?\d+\.\d{2})\s*(-?\d+\.\d{2})", line.strip())
    VERT[m.group(1)] = (float(m.group(2)), float(m.group(3)))
COLS = {}
cur = None
for line in section_lines("GRID"):
    p = line.split()
    if len(p) >= 2:
        cur = re.match(r"^[a-z]+", p[0]).group(0)
        COLS[cur] = []
    elif len(p) == 1:
        COLS[cur].append(VERT[p[0]])
COL_POLY = {c: np.array(v) for c, v in COLS.items()}
COL_PATH = {c: MplPath(v) for c, v in COL_POLY.items()}
log(f"grid: {len(VERT)} vertices, {len(COLS)} columns, {len(LAYER)} layers, {len(SURF)} surface entries; grid top {GRID_TOP} mRL")
missing_surf = sorted(set(COLS) - set(SURF))
log(f"columns without a surface entry (grid top used): {len(missing_surf)}")


def col_surface(c):
    return SURF.get(c, GRID_TOP)


def col_top_layer(c):
    s = col_surface(c)
    return min(i for i in LAYER if LAYER[i]["bottom"] < s)


def column_at(x, y):
    for c, pth in COL_PATH.items():
        if pth.contains_point((x, y)):
            return c
    return None


# model boundary: edges used by exactly one column
_edge_count = {}
for c, v in COL_POLY.items():
    n = len(v)
    for i in range(n):
        a, b = tuple(v[i]), tuple(v[(i + 1) % n])
        k = (a, b) if a < b else (b, a)
        _edge_count[k] = _edge_count.get(k, 0) + 1
BOUNDARY_EDGES = [k for k, n in _edge_count.items() if n == 1]


def inside_model(x, y):
    return column_at(x, y) is not None


WB = json.load(open(WBLOCKS))
for w, d in WB.items():
    for k, v in list(d.items()):
        if isinstance(v, str):
            try:
                d[k] = ast.literal_eval(v)
            except (ValueError, SyntaxError):
                pass


def split_block(b):
    m = re.match(r"^([a-z]+)(\d+)$", b)
    return m.group(1), int(m.group(2))


def block_at(w, z):
    """Block from the well-to-block file whose elevation interval contains z."""
    bl = WB[w]["Blocks"]
    for b in bl:
        top, bot = max(b["interval"]), min(b["interval"])
        if bot < z <= top:
            return b["block"]
    for b in bl:
        top, bot = max(b["interval"]), min(b["interval"])
        if bot - 0.5 <= z <= top + 0.5:
            return b["block"]
    return ""


def fz_blocks(w):
    out = []
    for _, blocks in WB[w]["Feed_Zones"].items():
        best = max(blocks.items(), key=lambda kv: kv[1].get("Fraction", 0))
        out.append(best[0])
    return out


# ------------------------------------------------------------------------------------------
# 4. Mapping
# ------------------------------------------------------------------------------------------
MAP = json.load(open(MAPPING))
FAULTS, WATER = {}, {}
for name, v in MAP.items():
    d = v["Data"]
    if v.get("type") == "multiplot":
        segs = [np.array(s, float) for s in d]
    else:
        segs = [np.array(d, float)]
    if "fault" in (v.get("legend", "") + name).lower():
        FAULTS[name] = segs
    else:
        WATER[name] = {"segs": segs, "fill": v.get("type") == "fill", "legend": v.get("legend", "")}
log("mapping: faults", list(FAULTS), "water bodies", list(WATER))


def dist_to_segs(x, y, segs, closed=False):
    best = np.inf
    for s in segs:
        pts = np.vstack([s, s[:1]]) if closed else s
        a, b = pts[:-1], pts[1:]
        ab = b - a
        L2 = (ab ** 2).sum(1)
        L2[L2 == 0] = 1e-12
        t = np.clip(((x - a[:, 0]) * ab[:, 0] + (y - a[:, 1]) * ab[:, 1]) / L2, 0, 1)
        px, py = a[:, 0] + t * ab[:, 0], a[:, 1] + t * ab[:, 1]
        best = min(best, float(np.min(np.hypot(px - x, py - y))))
    if closed:
        for s in segs:
            if MplPath(s).contains_point((x, y)):
                return 0.0
    return best


# ------------------------------------------------------------------------------------------
# 5. Cross-check JSON against the plain-column files (4.1.1)
# ------------------------------------------------------------------------------------------
REF_DATES = []
xrows = []
datfiles = sorted(f for f in os.listdir(DAT_DIR) if f.endswith(".dat"))
aut2 = {W[w]["aut2"].lower(): w for w in W}
irregular = []
pending = []
for f in datfiles:
    m = re.match(r"^([A-Za-z]+\d+)_(.*)\.dat$", f)
    prefix, qty = m.group(1), m.group(2)
    w = aut2.get(prefix.lower())
    if w is None:
        xrows.append(("", qty, f, 0, 0, 0, np.nan, "no matching well"))
        continue
    if prefix != W[w]["aut2"] and prefix.lower() != w.replace("-", "").lower():
        irregular.append(f"{f} -> {w} (AUT2 name {W[w]['aut2']})")
    a = np.loadtxt(os.path.join(DAT_DIR, f), ndmin=2)
    r = W[w]
    if qty in ("downhole_temperature", "downhole_pressure"):
        col = "T" if qty.endswith("temperature") else "P"
        js = dict(zip(r["survey"].z, r["survey"][col]))
        dd = {float(z): float(v) for v, z in a}
        common = set(js) & set(dd)
        mx = max(abs(js[z] - dd[z]) for z in common) if common else np.nan
        xrows.append((w, qty, f, len(js), len(dd), len(common), mx, "value elevation_mRL columns"))
        continue
    # time series: years since a reference date
    if qty.endswith("transient_pressure"):
        elev = float(qty.split("_")[0])
        js = r["gauges"].get(elev, {})
    else:
        key = {"production_mass_flow": "Production_Mass_Flow", "transient_enthalpy": "Transient_Enthalpy",
               "injection_mass_flow": "Injection_Mass_Flow", "injection_enthalpy": "Injection_Enthalpy"}[qty]
        js = r["flows"].get(key, {})
    jd = sorted(js)
    pending.append((w, qty, f, js, jd, a))
from collections import Counter
cands = Counter()
for w, qty, f, js, jd, a in pending:
    for k in (math.floor, math.ceil):
        cands[jd[0] - dt.timedelta(days=k(a[0, 0] * 365.25))] += 1
REF_DATE = max(cands, key=lambda d: (cands[d], -d.toordinal()))
for w, qty, f, js, jd, a in pending:
    sign, mx, nmatch = 1.0, 0.0, 0
    for t, v in a:
        day = REF_DATE + dt.timedelta(days=t * 365.25)
        cand = [j for j in jd if abs((dt.datetime.combine(j, dt.time()) - dt.datetime.combine(REF_DATE, dt.time())
                                      - dt.timedelta(days=t * 365.25)).total_seconds()) <= 3.01 * 86400]
        if cand:
            nmatch += 1
            jv = js[cand[0]]
            if abs(v + jv) < abs(v - jv) and jv != 0:
                sign = -1.0
            mx = max(mx, abs(abs(v) - abs(jv)))
    note = f"time = years (365.25 d) since {REF_DATE.isoformat()}"
    if sign < 0:
        note += "; dat values are the negative of the JSON (generation sign convention)"
    xrows.append((w, qty, f, len(js), len(a), nmatch, mx, note))
XCHK = pd.DataFrame(xrows, columns=["well", "quantity", "dat_file", "n_json", "n_dat", "n_matched",
                                     "max_abs_diff", "note"])
XCHK["max_abs_diff"] = XCHK.max_abs_diff.round(4)
write_csv(XCHK, os.path.join(ADD, "dat_json_crosscheck.csv"))
log("dat reference date:", REF_DATE.isoformat(), "; irregular names:", irregular)

# ------------------------------------------------------------------------------------------
# 6. Integrity checks (4.1.1)
# ------------------------------------------------------------------------------------------
irows = []
for w in WELLS:
    r = W[w]
    s = r["survey"]
    wh, td = r["wellhead"], r["td_track"]
    def add(check, ok, detail):
        irows.append((w, check, "OK" if ok else "ISSUE", detail))
    add("survey within track", s.z.max() <= wh + 0.01 and s.z.min() >= td - 0.01,
        f"survey {s.z.max()} to {s.z.min()} mRL; track {wh} to {td} mRL")
    add("survey T and P at same elevations", s["T"].notna().all() and s["P"].notna().all(),
        f"{len(s)} points")
    add("T and P survey dates equal", r["survey_date"] == r["survey_date_P"], r["survey_date"].isoformat())
    for i, (top, bot, q) in enumerate(r["fz"], 1):
        add(f"feed zone {i} within track", top <= wh and bot >= td - 0.01, f"{top} to {bot} mRL")
        ins = any(top <= st + 0.01 and bot >= sb - 0.01 for st, sb in r["slotted"])
        add(f"feed zone {i} within slotted interval", ins,
            f"slotted {[(a, b) for a, b in r['slotted']]}")
    des = r["design"]
    add("casing top at wellhead", abs(max(e[0] for e in des) - wh) < 0.5,
        f"top of design {max(e[0] for e in des)} mRL; wellhead {wh} mRL")
    deepest = min(min(e[0], e[1]) for e in des)
    add("deepest liner bottom at track TD", abs(deepest - td) < 0.5,
        f"design bottom {deepest} mRL; track TD {td} mRL")
    laps = all(any(lt >= r["shoe"] and lb < r["shoe"] for lt, lb in r["slotted"]) for _ in [0])
    add("liner overlaps production casing shoe", laps, f"shoe {r['shoe']} mRL; liners {r['slotted']}")
    # wellhead vs geometry file WELLS section and block file
    blk = WB[w]["Blocks"]
    add("wellhead = first block-file interval top", abs(max(blk[0]["interval"]) - wh) < 0.5,
        f"block file top {max(blk[0]['interval'])} mRL")
INTEG = pd.DataFrame(irows, columns=["well", "check", "status", "detail"])
write_csv(INTEG, os.path.join(ADD, "integrity_checks.csv"))
for _, rr in INTEG[INTEG.status != "OK"].iterrows():
    log(f"INTEGRITY ISSUE {rr.well}: {rr.check}: {rr.detail}")

# ------------------------------------------------------------------------------------------
# 7. Flow timing and injectate temperature
# ------------------------------------------------------------------------------------------
def first_flow(r):
    out = None
    for key in ("Production_Mass_Flow", "Injection_Mass_Flow"):
        for d_, v in sorted(r["flows"].get(key, {}).items()):
            if abs(v) > 0:
                out = d_ if out is None or d_ < out else out
                break
    return out


FIRST_FLOW = {w: first_flow(W[w]) for w in WELLS}
prod_starts = [first_flow({"flows": {"Production_Mass_Flow": W[w]["flows"]["Production_Mass_Flow"]}})
               for w in WELLS if W[w]["producer"]]
area_start = None
for p in AREA_JSON:
    a = json.load(open(p))
    for area, qs in a.items():
        for q, v in qs.items():
            ds = sorted(ddmmyyyy(k) for k, x in v["Data"].items() if x > 0)
            if ds:
                area_start = ds[0] if area_start is None else min(area_start, ds[0])
FIELD_START = min([d for d in prod_starts if d] + ([area_start] if area_start else []))
log(f"first production anywhere in the field: {FIELD_START} (wells: {min(d for d in prod_starts if d)}, field totals: {area_start})")

inj_h = [v for w in WELLS for v in W[w]["flows"].get("Injection_Enthalpy", {}).values()]
H_INJ = float(np.median(inj_h)) if inj_h else np.nan


def t_from_h(h, p_bar=P_ATM):
    lo, hi = 0.01, 99.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if h_liq(mid, p_bar) < h:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


T_INJ = t_from_h(H_INJ) if inj_h else np.nan
log(f"injectate: median enthalpy {H_INJ:.1f} kJ/kg -> T = {T_INJ:.1f} degC (IF97 Region 1 inverse at 1.01325 bar)")

# ------------------------------------------------------------------------------------------
# 8. Survey points: geometry flags (targets table without confidence yet)
# ------------------------------------------------------------------------------------------
QUAL_ORDER = ["HIGH", "MEDIUM", "LOW"]
TAGS = sorted({q for w in WELLS for _, _, q in W[w]["fz"]})
grades = [q for q in QUAL_ORDER if q in [t.upper() for t in TAGS]]


def qfactor(tag):
    t = str(tag).upper()
    if grades and t == grades[0]:
        return 1.0
    if t in grades and t != grades[-1]:
        return 0.8
    return 0.5


TAG_MAP = {t: qfactor(t) for t in TAGS}
log("feed-zone quality tags and weight_P factors:", TAG_MAP)


def isothermal_flags(z, t, openhole):
    """Points in runs of consecutive open-hole points spanning >= 300 m with T range <= 5 degC."""
    n = len(z)
    flag = np.zeros(n, bool)
    for i in range(n):
        if not openhole[i]:
            continue
        lo = hi = t[i]
        for j in range(i + 1, n):
            if not openhole[j]:
                break
            lo, hi = min(lo, t[j]), max(hi, t[j])
            if hi - lo > 5.0:
                break
            if z[i] - z[j] >= 300:
                flag[i:j + 1] = True
    return flag


rows = []
for w in WELLS:
    r = W[w]
    s = r["survey"]
    z, t, p = s.z.to_numpy(), s["T"].to_numpy(), s.P.to_numpy()
    openhole = z < r["shoe"]
    iso = isothermal_flags(z, t, openhole)
    r["iso_flags"] = iso
    for i in range(len(z)):
        x, y = track_xy(w, z[i])
        blk = block_at(w, z[i])
        col, lay = split_block(blk) if blk else ("", np.nan)
        ts = tsat(p[i])
        fzi = [k for k, (top, bot, q) in enumerate(r["fz"], 1) if bot <= z[i] <= top]
        above = bool(z[i] > col_surface(col)) if col else False
        rows.append(dict(well=w, x=round(x, 1), y=round(y, 1), elevation_mRL=z[i], T_degC=t[i], P_bar=p[i],
                         Tsat_degC=round(ts, 2) if not math.isnan(ts) else np.nan,
                         column=col, layer=lay, block=blk, in_feedzone=bool(fzi),
                         behind_casing=bool(z[i] >= r["shoe"]), above_model_surface=above,
                         _fz=fzi[0] if fzi else 0, _iso=bool(iso[i]),
                         _top=bool(col and lay == col_top_layer(col)),
                         _near_wh=bool(r["wellhead"] - z[i] <= 30)))
TG = pd.DataFrame(rows)
TG["boiling_margin_degC"] = (TG.Tsat_degC - TG.T_degC).round(2)
TG["layer"] = TG.layer.astype("Int64")

# ------------------------------------------------------------------------------------------
# 9. Feed zones (survey values at mid elevation)
# ------------------------------------------------------------------------------------------
frows = []
fzb = {w: fz_blocks(w) for w in WELLS}
for w in WELLS:
    r = W[w]
    s = r["survey"]
    if len(fzb[w]) != len(r["fz"]):
        log(f"WARNING {w}: {len(r['fz'])} feed zones in the data, {len(fzb[w])} in the block file")
    for i, (top, bot, q) in enumerate(r["fz"], 1):
        mid = 0.5 * (top + bot)
        T_ = round(interp_profile(s.z, s["T"], mid), 1)
        P_ = round(interp_profile(s.z, s.P, mid), 1)
        ts = tsat(P_)
        col, lay = split_block(fzb[w][i - 1])
        frows.append(dict(well=w, feedzone=i, top_mRL=top, bottom_mRL=bot, mid_mRL=mid, quality=q,
                          T_degC=T_, P_bar=P_, Tsat_degC=round(ts, 2),
                          boiling_margin_degC=round(round(ts, 2) - T_, 2), block=fzb[w][i - 1],
                          _qf=qfactor(q), _thick=top - bot, _above=mid > col_surface(col),
                          _near_wh=r["wellhead"] - mid <= 30, _top=lay == col_top_layer(col),
                          survey_point_at_mid=bool(np.isclose(s.z, mid).any())))
FZ = pd.DataFrame(frows)

# ------------------------------------------------------------------------------------------
# 10. Natural-state test (4.1.2)
# ------------------------------------------------------------------------------------------
def decimals(vals):
    d = 0
    for v in vals:
        s = repr(float(v))
        if "." in s:
            d = max(d, len(s.split(".")[1].rstrip("0")))
    return d


gauge_rows = {}
for w in WELLS:
    r = W[w]
    if not r["gauges"]:
        continue
    ge = sorted(r["gauges"])[0]
    series = r["gauges"][ge]
    d0 = min(series)
    p0 = series[d0]
    prec = 10.0 ** (-decimals(list(series.values())))
    sp = round(interp_profile(r["survey"].z, r["survey"].P, ge), 2)
    gauge_rows[w] = dict(gauge_elev_mRL=ge, survey_P_at_gauge_bar=sp, first_gauge_date=d0, first_gauge_P_bar=p0,
                         diff=round(sp - p0, 2), precision=prec, preprod=d0 <= FIELD_START,
                         own_start=FIRST_FLOW[w], series=series)

verdict, conf, reason = {}, {}, {}
# step 1: gauge-verified wells
for w, g in gauge_rows.items():
    if g["preprod"] and abs(g["diff"]) <= g["precision"] + 1e-9:
        verdict[w], conf[w] = "natural state", "high"
        reason[w] = (f"survey P at the gauge ({g['gauge_elev_mRL']:.0f} mRL) {g['survey_P_at_gauge_bar']:.2f} bar vs "
                     f"first gauge reading {g['first_gauge_P_bar']:.1f} bar on {g['first_gauge_date'].isoformat()} "
                     f"(difference {g['diff']:+.2f} bar, within the {g['precision']:.1f} bar reading precision); "
                     f"the reading is not later than the first field production ({FIELD_START.isoformat()})")
VERIFIED = [w for w in WELLS if conf.get(w) == "high"]
vf = FZ[FZ.well.isin(VERIFIED)]
if len(vf) >= 2:
    TREND_B, TREND_A = np.polyfit(vf.mid_mRL, vf.P_bar, 1)
else:
    raise SystemExit("fewer than two gauge-verified feed zones; cannot build the pressure trend")
res_v = vf.P_bar - (TREND_A + TREND_B * vf.mid_mRL)
RMS_V = float(np.sqrt(np.mean(res_v ** 2)))
TOL = max(3.0, 2 * RMS_V)
log(f"verified-trend: P = {TREND_A:.3f} + {TREND_B:.5f} z (bar, mRL) from {len(vf)} feed zones of {VERIFIED}; "
    f"RMS {RMS_V:.3f} bar; tolerance {TOL:.2f} bar")


def trend(z):
    return TREND_A + TREND_B * np.asarray(z, float)


FZ["P_minus_trend_bar"] = (FZ.P_bar - trend(FZ.mid_mRL)).round(2)
fz_maxres = FZ.groupby("well").P_minus_trend_bar.apply(lambda s: float(s.abs().max())).to_dict()

for w in WELLS:
    if w in verdict:
        continue
    r = W[w]
    g = gauge_rows.get(w)
    on_trend = fz_maxres.get(w, np.inf) <= TOL
    resid_txt = (f"feed-zone pressures lie within {fz_maxres.get(w, np.nan):.2f} bar of the verified trend "
                 f"(tolerance {TOL:.2f} bar)")
    fzT = FZ[FZ.well == w].T_degC.max()
    iso_T = [t for t, f in zip(r["survey"]["T"], r["iso_flags"]) if f]
    inj_sig = r["injector"] and any(abs(t - T_INJ) <= 10 for t in iso_T)
    if r["injector"] and fzT > T_INJ:
        verdict[w], conf[w] = "probably natural state at best", "low"
        reason[w] = (f"injector since {FIRST_FLOW[w].isoformat()}; feed-zone temperatures (up to {fzT:.1f} degC) are "
                     f"hotter than the injectate ({T_INJ:.1f} degC), so injection cooling could be hidden; "
                     + resid_txt)
    elif inj_sig:
        verdict[w], conf[w] = "possibly disturbed", "low"
        reason[w] = f"injector with an isothermal section near the injectate temperature ({T_INJ:.1f} degC)"
    elif g is not None and g["preprod"] and on_trend:
        verdict[w], conf[w] = "probably natural state", "medium"
        reason[w] = (f"gauge test misses: survey P at {g['gauge_elev_mRL']:.0f} mRL is {g['survey_P_at_gauge_bar']:.2f} bar "
                     f"vs first reading {g['first_gauge_P_bar']:.1f} bar (difference {g['diff']:+.2f} bar, more than the "
                     f"{g['precision']:.1f} bar precision; recorded as a data question); " + resid_txt)
    elif on_trend:
        verdict[w], conf[w] = "probably natural state", "medium"
        why = ("no gauge" if g is None else
               f"gauge series starts {g['first_gauge_date'].isoformat()}, after first field production "
               f"({FIELD_START.isoformat()})")
        reason[w] = f"no pre-production gauge reading ({why}); " + resid_txt
        if r["injector"]:
            reason[w] += (f"; injector whose feed zones (up to {fzT:.1f} degC) are colder than the injectate "
                          f"({T_INJ:.1f} degC) with no isothermal section near it, which is evidence against an "
                          f"injection-affected survey")
    else:
        verdict[w], conf[w] = "possibly disturbed", "low"
        reason[w] = (f"feed-zone pressures off the verified trend by up to {fz_maxres.get(w, np.nan):.2f} bar "
                     f"(tolerance {TOL:.2f} bar)")
CONF_F = {"high": 1.0, "medium": 0.75, "low": 0.5}
for w in WELLS:
    log(f"NS verdict {w}: {verdict[w]} / {conf[w]} -- {reason[w]}")

pcrows = []
for w in WELLS:
    g = gauge_rows.get(w)
    pcrows.append(dict(
        well=w,
        gauge_elev_mRL=g["gauge_elev_mRL"] if g else np.nan,
        survey_P_at_gauge_bar=g["survey_P_at_gauge_bar"] if g else np.nan,
        first_gauge_date=g["first_gauge_date"].isoformat() if g else "",
        first_gauge_P_bar=g["first_gauge_P_bar"] if g else np.nan,
        diff_survey_minus_first_bar=g["diff"] if g else np.nan,
        verdict=verdict[w], confidence=conf[w],
        survey_date=W[w]["survey_date"].isoformat(),
        gauge_precision_bar=g["precision"] if g else np.nan,
        gauge_before_or_at_first_field_production=g["preprod"] if g else "",
        first_own_flow_date=FIRST_FLOW[w].isoformat() if FIRST_FLOW[w] else "",
        well_type="injector" if W[w]["injector"] else ("producer" if W[w]["producer"] else "monitoring/none"),
        first_field_production_date=FIELD_START.isoformat(),
        fz_max_abs_P_minus_trend_bar=round(fz_maxres.get(w, np.nan), 2),
        trend_tolerance_bar=round(TOL, 2),
        reason=reason[w]))
PC = pd.DataFrame(pcrows)
write_csv(PC, os.path.join(ADD, "pressure_consistency_check.csv"))

# ------------------------------------------------------------------------------------------
# 11. Weights (5.5) and the targets table
# ------------------------------------------------------------------------------------------
def loc_weight(above, near_wh, top, in_fz, behind, iso):
    if above:
        return 0.0, "above model surface"
    if near_wh or top:
        return 0.3, "within 30 m of wellhead" if near_wh else "top block of column"
    if in_fz:
        return 1.0, "feed zone"
    if behind:
        return 0.5, "behind casing"
    if iso:
        return 0.5, "open hole, isothermal section"
    return 0.8, "open hole"


qf_by = {(a, b): c for a, b, c in zip(FZ.well, FZ.feedzone, FZ._qf)}
TG["ns_confidence"] = TG.well.map(conf)
wts, why, wps = [], [], []
for rr in TG.rename(columns=lambda c: c.lstrip("_")).itertuples():
    lw, rs = loc_weight(rr.above_model_surface, rr.near_wh, rr.top, rr.in_feedzone, rr.behind_casing, rr.iso)
    cf = CONF_F[rr.ns_confidence]
    wts.append(round(lw * cf, 4))
    why.append(f"{rs} ({lw}) x {rr.ns_confidence} confidence ({cf})")
    wps.append(round(qf_by[(rr.well, rr.fz)] * cf, 4) if rr.in_feedzone and not rr.above_model_surface else 0.0)
TG["weight"], TG["weight_reason"], TG["weight_P"] = wts, why, wps
TG["isothermal_section"] = TG._iso
TG["feedzone"] = TG._fz.replace(0, np.nan).astype("Int64")
TG["liquid_filled"] = (TG.P_bar >= 1.5) & ((TG.boiling_margin_degC >= 0) | TG.boiling_margin_degC.isna())
TG["near_boiling_lt5C"] = TG.boiling_margin_degC < 5
TG["P_minus_trend_bar"] = (TG.P_bar - trend(TG.elevation_mRL)).round(2)
TCOLS = ["well", "x", "y", "elevation_mRL", "T_degC", "P_bar", "Tsat_degC", "boiling_margin_degC", "column", "layer",
         "block", "in_feedzone", "behind_casing", "above_model_surface", "ns_confidence", "weight", "weight_reason",
         "weight_P", "feedzone", "isothermal_section", "liquid_filled", "near_boiling_lt5C", "P_minus_trend_bar"]
TG = TG[TCOLS]
write_csv(TG, os.path.join(CAL, "NS_calibration_targets.csv"))

FZ["ns_confidence"] = FZ.well.map(conf)
FZ["weight_T"] = [round(loc_weight(a, n, t, True, False, False)[0] * CONF_F[c], 4)
                  for a, n, t, c in zip(FZ._above, FZ._near_wh, FZ._top, FZ.ns_confidence)]
FZ["weight_P"] = [round(q * CONF_F[c], 4) for q, c in zip(FZ._qf, FZ.ns_confidence)]
FZ["quality_factor"] = FZ._qf
FZ["column"] = [split_block(b)[0] for b in FZ.block]
FZ["layer"] = [split_block(b)[1] for b in FZ.block]
FZCOLS = ["well", "feedzone", "top_mRL", "bottom_mRL", "mid_mRL", "quality", "T_degC", "P_bar", "Tsat_degC",
          "boiling_margin_degC", "block", "ns_confidence", "weight_T", "weight_P", "quality_factor",
          "P_minus_trend_bar", "column", "layer", "survey_point_at_mid"]
FZOUT = FZ[FZCOLS]
write_csv(FZOUT, os.path.join(CAL, "feedzones.csv"))

# ------------------------------------------------------------------------------------------
# 12. Per-well characteristics (4.4), computed from the targets table
# ------------------------------------------------------------------------------------------
def fit_slope(z, v):
    return np.polyfit(np.asarray(z, float), np.asarray(v, float), 1)


def crossing(z, t, x):
    for i in range(len(t) - 1):
        if t[i] < x <= t[i + 1]:
            return z[i] + (x - t[i]) * (z[i + 1] - z[i]) / (t[i + 1] - t[i])
    return np.nan


srows = []
for w in WELLS:
    g = TG[TG.well == w].sort_values("elevation_mRL", ascending=False)
    z, t = g.elevation_mRL.to_numpy(float), g.T_degC.to_numpy(float)
    p, bm = g.P_bar.to_numpy(float), g.boiling_margin_degC.to_numpy(float)
    e = {"well": w}
    i = int(np.nanargmax(t))
    e["Tmax_degC"], e["Tmax_elev_mRL"] = t[i], z[i]
    e["T_at_TD_degC"], e["TD_mRL"] = t[-1], z[-1]
    run_max, run_z, best = -np.inf, np.nan, (0.0, np.nan, np.nan)
    for zi, ti in zip(z, t):
        if ti > run_max:
            run_max, run_z = ti, zi
        if run_max - ti > best[0]:
            best = (run_max - ti, run_z, zi)
    inv = best if best[0] >= 3 else (np.nan, np.nan, np.nan)
    e["inversion_degC"] = round(inv[0], 2) if not np.isnan(inv[0]) else np.nan
    e["inversion_top_mRL"], e["inversion_base_mRL"] = inv[1], inv[2]
    liq = (p >= 1.5) & ((bm >= 0) | (np.isnan(bm) & (p > P_CRIT)))
    zl, tl, pl, bml = z[liq], t[liq], p[liq], bm[liq]
    e["shallow_gradient_degC_per_km"] = np.nan
    if len(zl):
        m = zl >= zl[0] - 300
        if m.sum() >= 3:
            e["shallow_gradient_degC_per_km"] = round(-fit_slope(zl[m], tl[m])[0] * 1000, 1)
    fzw = FZ[FZ.well == w]
    fzbot = fzw.bottom_mRL.min() if len(fzw) else np.nan
    bounds = [b for b in (inv[2], fzbot) if not np.isnan(b)]
    sel = z <= min(bounds) if bounds else np.ones_like(z, bool)
    e["deep_gradient_degC_per_km"] = np.nan
    span = z[sel].max() - z[sel].min() if sel.sum() else 0.0
    if sel.sum() >= 3 and span >= 300:
        e["deep_gradient_degC_per_km"] = round(-fit_slope(z[sel], t[sel])[0] * 1000, 1)
    e["z_100C_mRL"] = round(crossing(z, t, 100.0), 1)
    e["z_200C_mRL"] = round(crossing(z, t, 200.0), 1)
    ok = ~np.isnan(bml)
    if ok.any():
        j = int(np.argmin(np.where(ok, bml, np.inf)))
        e["min_boiling_margin_degC"], e["min_boiling_margin_elev_mRL"] = bml[j], zl[j]
    else:
        e["min_boiling_margin_degC"] = e["min_boiling_margin_elev_mRL"] = np.nan
    e["water_level_mRL"] = np.nan
    wl_ex = np.nan
    if len(zl) >= 3:
        b, a = fit_slope(zl[:3], pl[:3])
        if b != 0:
            e["water_level_mRL"] = round((P_ATM - a) / b, 1)
            wl_ex = round(e["water_level_mRL"] - z[0], 1)
    if len(fzw):
        cand = fzw.assign(th=fzw.top_mRL - fzw.bottom_mRL).sort_values(
            ["_qf", "th", "mid_mRL"], ascending=[False, False, True])
        pv = cand.iloc[0]
        e["pivot_elev_mRL"], e["pivot_P_bar"] = pv.mid_mRL, pv.P_bar
        pivot_fz = int(pv.feedzone)
    else:
        e["pivot_elev_mRL"] = e["pivot_P_bar"] = np.nan
        pivot_fz = np.nan
    near_b = ok.any() and np.nanmin(bml) < 5
    if ok.any() and np.nanmin(bml) <= 1:
        reg = "two-phase"
    elif e["Tmax_degC"] < 100:
        reg = "recharge/downflow"
    elif not np.isnan(inv[0]) and inv[0] >= 10:
        reg = "outflow"
    elif e["Tmax_degC"] >= 200 and np.isnan(inv[0]):
        reg = "upflow"
    else:
        reg = "conductive"
    e["regime"] = reg
    e["ns_verdict"], e["ns_confidence"] = verdict[w], conf[w]
    e["deep_fit_n"] = int(sel.sum())
    e["deep_fit_span_m"] = round(span, 1)
    e["water_level_extrapolation_m"] = wl_ex
    e["near_boiling"] = bool(near_b)
    e["pivot_feedzone"] = pivot_fz
    # cap base: first top-down pair of points below which the local gradient drops below 50 degC/km
    cap = np.nan
    if e["Tmax_degC"] >= 200:
        for k in range(len(z) - 1):
            gk = (t[k + 1] - t[k]) / (z[k] - z[k + 1]) * 1000
            if z[k] < z[0] - 100 and gk < 50:
                cap = z[k]
                break
    e["cap_base_mRL"] = cap
    e["wellhead_mRL"] = W[w]["wellhead"]
    e["casing_shoe_mRL"] = W[w]["shoe"]
    srows.append(e)
WS = pd.DataFrame(srows)
REGIME_NOTE = {}
for rr in WS.itertuples():
    note = ""
    if rr.regime == "upflow" and rr.TD_mRL > -500:
        note = (f"shallow well (TD {rr.TD_mRL:.0f} mRL): no inversion is seen within its depth, so 'upflow' is the "
                f"rule result; the well may sit on the top of the upflow or on the near-source part of the outflow")
    if rr.regime == "recharge/downflow" and not np.isnan(rr.inversion_degC):
        note = (f"cold well with a {rr.inversion_degC:.1f} degC inversion between {rr.inversion_top_mRL:.0f} and "
                f"{rr.inversion_base_mRL:.0f} mRL: cool lateral inflow or downflow at that depth")
    if rr.regime == "conductive" and rr.Tmax_degC >= 100:
        note = ("monotonic profile with no inversion; warm conductive margin, consistent with distal heating by the "
                "outflow above or a conductive flank")
    if rr.regime == "outflow" and rr.well in VERIFIED and not np.isnan(rr.deep_gradient_degC_per_km) and rr.deep_gradient_degC_per_km > 0:
        note = "outflow inversion above a deep temperature increase towards TD: heat source also at depth"
    REGIME_NOTE[rr.well] = note
WS["regime_note"] = WS.well.map(REGIME_NOTE)
WS["ns_reason"] = WS.well.map(reason)
WCOLS = ["well", "Tmax_degC", "Tmax_elev_mRL", "T_at_TD_degC", "TD_mRL", "inversion_degC", "inversion_top_mRL",
         "inversion_base_mRL", "shallow_gradient_degC_per_km", "deep_gradient_degC_per_km", "z_100C_mRL", "z_200C_mRL",
         "min_boiling_margin_degC", "min_boiling_margin_elev_mRL", "water_level_mRL", "pivot_elev_mRL", "pivot_P_bar",
         "regime", "ns_verdict", "ns_confidence", "deep_fit_n", "deep_fit_span_m", "water_level_extrapolation_m",
         "near_boiling", "regime_note", "pivot_feedzone", "cap_base_mRL", "wellhead_mRL", "casing_shoe_mRL", "ns_reason"]
WS = WS[WCOLS]
write_csv(WS, os.path.join(CAL, "well_summary.csv"))
WSI = WS.set_index("well")

# ------------------------------------------------------------------------------------------
# 13. T and P at layer centres (no extrapolation)
# ------------------------------------------------------------------------------------------
lrows = []
for w in WELLS:
    s = W[w]["survey"]
    zmax, zmin = s.z.max(), s.z.min()
    for lay in sorted(LAYER):
        c = LAYER[lay]["centre"]
        if not (zmin <= c <= zmax):
            continue
        T_ = round(interp_profile(s.z, s["T"], c), 2)
        P_ = round(interp_profile(s.z, s.P, c), 2)
        ts = tsat(P_)
        lrows.append(dict(well=w, layer=lay, layer_centre_mRL=c, block=block_at(w, c), T_degC=T_, P_bar=P_,
                          Tsat_degC=round(ts, 2) if not math.isnan(ts) else np.nan,
                          boiling_margin_degC=round(round(ts, 2) - T_, 2) if not math.isnan(ts) else np.nan))
LC = pd.DataFrame(lrows)
write_csv(LC, os.path.join(CAL, "T_P_at_layer_centres.csv"))

# ------------------------------------------------------------------------------------------
# 14. Hydrostatic check (4.3)
# ------------------------------------------------------------------------------------------
def dp_if97(zt, zb, zs, ts, ps):
    n = max(2, int(math.ceil((zt - zb) / 5.0)) + 1)
    zz = np.linspace(zt, zb, n)
    rho = np.array([rho_liq(interp_profile(zs, ts, q), interp_profile(zs, ps, q)) for q in zz])
    return float(np.trapezoid(rho, -zz) * G / 1e5)


hrows = []
for w in WELLS:
    g = TG[TG.well == w].sort_values("elevation_mRL", ascending=False).reset_index(drop=True)
    shoe = W[w]["shoe"]
    z, t, p, bm = g.elevation_mRL.to_numpy(), g.T_degC.to_numpy(), g.P_bar.to_numpy(), g.boiling_margin_degC.to_numpy()
    infz = g.in_feedzone.to_numpy()
    elig = (p >= 1.5) & ((bm >= 10) | np.isnan(bm))
    runs, cur = [], []
    for k in range(len(z)):
        if elig[k]:
            cur.append(k)
        else:
            if cur:
                runs.append(cur)
            cur = []
    if cur:
        runs.append(cur)
    parts = []
    for run in runs:
        cased = [k for k in run if z[k] >= shoe]
        opn = [k for k in run if z[k] < shoe]
        for part in (cased, opn):
            if len(part) >= 2:
                parts.append(part)
    for part in parts:
        wins, start = [], part[0]
        for k in part[1:]:
            if z[start] - z[k] >= 300 and (infz[k] or k == part[-1]):
                wins.append([start, k])
                start = k
        if start != part[-1]:
            if wins and z[start] - z[part[-1]] < 300:
                wins[-1][1] = part[-1]
            else:
                wins.append([start, part[-1]])
        for a, b in wins:
            dpm = p[b] - p[a]
            dpi = dp_if97(z[a], z[b], z, t, p)
            L = z[a] - z[b]
            band = "cased" if z[b] >= shoe else "open"
            hrows.append(dict(well=w, z_top_mRL=z[a], z_bottom_mRL=z[b], length_m=round(L, 1), band=band,
                              dP_measured_bar=round(dpm, 2), dP_IF97_bar=round(dpi, 3),
                              ratio_measured_to_IF97=round(dpm / dpi, 4),
                              rho_implied_kgm3=round(dpm * 1e5 / (G * L), 1),
                              ratio_rounding_uncertainty=round(0.1 / dpi, 4),
                              T_top_degC=t[a], T_bottom_degC=t[b], n_points=b - a + 1,
                              cold_well=bool(WSI.loc[w, "Tmax_degC"] < 100)))
HY = pd.DataFrame(hrows)
cold_comp = HY[HY.cold_well].groupby("band").ratio_measured_to_IF97.median().to_dict()
HY["cold_well_component"] = HY.band.map(cold_comp).round(4)
HY["excess_over_cold_component"] = (HY.ratio_measured_to_IF97 - HY.cold_well_component).round(4)
write_csv(HY, os.path.join(ADD, "hydrostatic_check.csv"))
hsum = []
for w in WELLS + ["ALL"]:
    for band in ("cased", "open"):
        h = HY[(HY.band == band) & ((HY.well == w) | (w == "ALL"))]
        if len(h):
            hsum.append(dict(well=w, band=band, n_windows=len(h), median_ratio=round(h.ratio_measured_to_IF97.median(), 4),
                             cold_well_component=round(cold_comp.get(band, np.nan), 4),
                             median_excess=round(h.ratio_measured_to_IF97.median() - cold_comp.get(band, np.nan), 4),
                             median_rho_implied_kgm3=round(h.rho_implied_kgm3.median(), 1)))
HYS = pd.DataFrame(hsum)
HYS = HYS[[c for c in HYS.columns if c != "n_windows"] + ["n_windows"]]
write_csv(HYS, os.path.join(ADD, "hydrostatic_summary.csv"))
log("cold-well component (median ratio) by band:", {k: round(v, 4) for k, v in cold_comp.items()})

# ------------------------------------------------------------------------------------------
# 15. Well geometry check
# ------------------------------------------------------------------------------------------
grows = []
for w in WELLS:
    r = W[w]
    tr = r["track"]
    x0, y0 = tr[0, 0], tr[0, 1]
    top_block = max(WB[w]["Blocks"], key=lambda b: max(b["interval"]))["block"]
    col = split_block(top_block)[0]
    ms = col_surface(col)
    e = dict(well=w, wellhead_mRL=r["wellhead"], model_surface_mRL=ms,
             wellhead_minus_surface_m=round(r["wellhead"] - ms, 1),
             max_horizontal_deviation_m=round(float(np.max(np.hypot(tr[:, 0] - x0, tr[:, 1] - y0))), 1))
    fd = {n: dist_to_segs(x0, y0, s) for n, s in FAULTS.items()}
    wd = {n: dist_to_segs(x0, y0, v["segs"], closed=v["fill"]) for n, v in WATER.items()}
    nf = min(fd, key=fd.get)
    nw = min(wd, key=wd.get)
    e.update(nearest_fault=nf, dist_nearest_fault_m=round(fd[nf], 1), nearest_water_body=nw,
             dist_nearest_water_body_m=round(wd[nw], 1))
    for n in list(FAULTS) + list(WATER):
        e[f"dist_{n.replace(' ', '_')}_m"] = round(fd.get(n, wd.get(n)), 1)
    e["wellhead_column"] = col
    e["top_block_in_block_file"] = top_block
    e["column_top_layer"] = col_top_layer(col)
    e["column_at_wellhead_xy"] = column_at(x0, y0)
    wells_sec = [l.split() for l in section_lines("WELLS") if l.split() and l.split()[0] == r["aut2"]]
    e["wellhead_in_grid_WELLS_section_mRL"] = float(wells_sec[0][3]) if wells_sec else np.nan
    grows.append(e)
GEO = pd.DataFrame(grows)
write_csv(GEO, os.path.join(ADD, "well_geometry_check.csv"))

# ------------------------------------------------------------------------------------------
# 16. Elevation bands (F3 data) and common elevations
# ------------------------------------------------------------------------------------------
zmax_all, zmin_all = TG.elevation_mRL.max(), TG.elevation_mRL.min()
band_tops = []
bt = math.ceil(zmax_all / 500) * 500
while bt > zmin_all:
    band_tops.append(bt)
    bt -= 500
brows = []
for bt in band_tops:
    bb = bt - 500
    sub = TG[(TG.elevation_mRL <= bt) & (TG.elevation_mRL > bb)]
    if sub.empty:
        continue
    for w in WELLS:
        s = sub[sub.well == w]
        if len(s):
            brows.append(dict(well=w, band_top_mRL=bt, band_bottom_mRL=bb, n_points=len(s),
                              T_mean_degC=round(s.T_degC.mean(), 2), P_mean_bar=round(s.P_bar.mean(), 2),
                              P_minus_trend_bar=round(float((s.P_bar - trend(s.elevation_mRL)).mean()), 2),
                              z_mean_mRL=round(s.elevation_mRL.mean(), 1)))
BA = pd.DataFrame(brows)
nwb = BA.groupby("band_top_mRL").well.nunique()
deepest_ok = min(b for b, n in nwb.items() if n >= 3)
BA = BA[BA.band_top_mRL >= deepest_ok].copy()
med = BA.groupby("band_top_mRL").T_mean_degC.median()
BA["T_field_median_degC"] = BA.band_top_mRL.map(med).round(2)
BA["T_minus_field_median_degC"] = (BA.T_mean_degC - BA.T_field_median_degC).round(2)
BA = BA[["well", "band_top_mRL", "band_bottom_mRL", "n_points", "T_mean_degC", "T_field_median_degC",
         "T_minus_field_median_degC", "P_mean_bar", "P_minus_trend_bar", "z_mean_mRL"]]
BA["_o"] = BA.well.map({w: i for i, w in enumerate(WELLS)})
BA = BA.sort_values(["band_top_mRL", "_o"], ascending=[False, True]).drop(columns="_o")
write_csv(BA, os.path.join(ADD, "T_P_at_elevation_bands.csv"))
BANDS = sorted(BA.band_top_mRL.unique(), reverse=True)

COMMON_Z = [0, -250, -500, -750, -1000, -1500, -2000, -2500]
crow = []
for w in WELLS:
    s = W[w]["survey"]
    for zc in COMMON_Z:
        if s.z.min() <= zc <= s.z.max():
            T_ = round(interp_profile(s.z, s["T"], zc), 1)
            P_ = round(interp_profile(s.z, s.P, zc), 1)
            crow.append(dict(well=w, elevation_mRL=zc, T_degC=T_, P_bar=P_,
                             P_minus_trend_bar=round(P_ - float(trend(zc)), 1),
                             open_hole=bool(zc < W[w]["shoe"])))
CE = pd.DataFrame(crow)
write_csv(CE, os.path.join(ADD, "T_P_at_common_elevations.csv"))

# provenance table
PROV = pd.DataFrame(PROVENANCE, columns=["well", "item", "Source", "Comment", "StableFlag"])
write_csv(PROV, os.path.join(ADD, "provenance_fields.csv"))

log("analysis tables done")

# ==========================================================================================
# FIGURES (CLAUDE.md 8)
# ==========================================================================================
plt.rcParams.update({"font.size": 9, "figure.dpi": 100})
DPI = 160
WCOL = {w: plt.cm.tab10(i % 10) for i, w in enumerate(WELLS)}
T_LIM = math.ceil(TG.T_degC.max() / 50.0) * 50
CMAP_T = plt.cm.viridis
NORM_T = matplotlib.colors.Normalize(0, T_LIM)
Z_TOP = math.ceil(max(W[w]["wellhead"] for w in WELLS) / 100.0) * 100
Z_BOT = math.floor(TG.elevation_mRL.min() / 100.0) * 100
ISO_LEVELS = list(range(50, int(T_LIM) + 1, 50))
P_LIM = math.ceil(TG.P_bar.max() / 50.0) * 50
FIG_INDEX = []
_allx = np.concatenate([W[w]["track"][:, 0] for w in WELLS])
_ally = np.concatenate([W[w]["track"][:, 1] for w in WELLS])
FIELD_VIEW = (math.floor((_allx.min() - 1500) / 500) * 500, math.ceil((_allx.max() + 1500) / 500) * 500,
              math.floor((_ally.min() - 1500) / 500) * 500, math.ceil((_ally.max() + 1500) / 500) * 500)
_vx = np.array(list(VERT.values()))
MODEL_VIEW = (_vx[:, 0].min() - 300, _vx[:, 0].max() + 300, _vx[:, 1].min() - 300, _vx[:, 1].max() + 300)
LABEL_OFF = {w: (90, 90) for w in WELLS}


def save(fig, path, fid, shows, wells, elev, sources, claim):
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    if fid:
        FIG_INDEX.append(dict(figure_id=fid, file=rel(path), shows=shows, wells=wells, elevation_or_line=elev,
                              data_sources=sources, claim_type=claim))
    log(f"figure {rel(path)}")


def scale_bar(ax, length=None):
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    if length is None:
        span = x1 - x0
        length = 1000 if span < 12000 else 2000 if span < 25000 else 5000
    bx, by = x0 + 0.05 * (x1 - x0), y0 + 0.05 * (y1 - y0)
    ax.plot([bx, bx + length], [by, by], color="k", lw=3, solid_capstyle="butt", zorder=20)
    ax.text(bx + length / 2, by + 0.015 * (y1 - y0), f"{length / 1000:g} km", ha="center", va="bottom",
            fontsize=7, zorder=20)


def _in_view(pt, view):
    return view[0] <= pt[0] <= view[1] and view[2] <= pt[1] <= view[3]


def _shrink(view, f):
    dx, dy = view[1] - view[0], view[3] - view[2]
    return (view[0] + f * dx, view[1] - f * dx, view[2] + f * dy, view[3] - f * dy)


def base_map(ax, view=None, highlight=None, labels=True, fs=7, wells_on=True, feature_labels=True):
    view = view or FIELD_VIEW
    ax.add_collection(PolyCollection(list(COL_POLY.values()), facecolors="none", edgecolors="#c8c8c8",
                                     linewidths=0.35, zorder=1))
    ax.add_collection(LineCollection([list(e) for e in BOUNDARY_EDGES], colors="#555555", linewidths=1.0,
                                     linestyles="--", zorder=2))
    for name, v in WATER.items():
        for s in v["segs"]:
            if v["fill"]:
                ax.fill(s[:, 0], s[:, 1], color="#a6d8f5", ec="#1f78b4", lw=0.8, zorder=2)
            else:
                ax.plot(s[:, 0], s[:, 1], color="#1f78b4", lw=1.2, zorder=3)
        if feature_labels:
            pts = np.vstack(v["segs"])
            inside = [p for p in pts if _in_view(p, _shrink(view, 0.1))]
            if inside:
                p = inside[len(inside) // 3]
                ax.text(p[0], p[1], name, color="#1f78b4", fontsize=fs, style="italic", zorder=6,
                        bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.5))
    for name, segs in FAULTS.items():
        for s in segs:
            ax.plot(s[:, 0], s[:, 1], color="k", lw=1.6, zorder=4)
        if feature_labels:
            s = segs[0]
            ts = np.linspace(0, 1, 401)
            pts = s[0] + np.outer(ts, s[-1] - s[0])
            vin = [p for p in pts if _in_view(p, _shrink(view, 0.04))]
            if vin:
                p = vin[int(len(vin) * 0.08)]
                ax.text(p[0], p[1], " " + name, fontsize=fs, fontweight="bold", zorder=6,
                        bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.4))
    if wells_on:
        for w in WELLS:
            tr = W[w]["track"]
            c = WCOL[w] if highlight == w else "grey"
            ax.plot(tr[:, 0], tr[:, 1], color=c, lw=2.2 if highlight == w else 1.0, zorder=7)
            ax.plot(tr[0, 0], tr[0, 1], "o", color="k", ms=3.5, zorder=8)
            if highlight == w:
                ax.plot(tr[0, 0], tr[0, 1], "o", color=WCOL[w], ms=8, mec="k", zorder=9)
            if labels:
                dx, dy = LABEL_OFF[w]
                ax.annotate(w, (tr[0, 0], tr[0, 1]), xytext=(tr[0, 0] + dx, tr[0, 1] + dy), fontsize=fs,
                            zorder=10, color="k", bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.3))
    ax.set_xlim(view[0], view[1])
    ax.set_ylim(view[2], view[3])
    ax.set_aspect("equal")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.set_xlabel("Easting (m)", fontsize=7)
    ax.set_ylabel("Northing (m)", fontsize=7)
    ax.tick_params(labelsize=6)
    for lab in ax.get_xticklabels():
        lab.set_rotation(30)
    scale_bar(ax)


# ------------------------------------------------------------------------------------------
# Section lines
# ------------------------------------------------------------------------------------------
def model_mask(xs, ys):
    pts = np.column_stack([xs, ys])
    inside = np.zeros(len(pts), bool)
    for pth in COL_PATH.values():
        inside |= pth.contains_points(pts)
    return inside


def make_line(name, w1, w2, why, maxdist=500.0):
    p1, p2 = W[w1]["track"][0, :2], W[w2]["track"][0, :2]
    u = (p2 - p1) / np.linalg.norm(p2 - p1)
    if abs(u[0]) >= abs(u[1]):
        u = u if u[0] > 0 else -u
    else:
        u = u if u[1] > 0 else -u
    ts = np.arange(-30000, 30000, 10.0)
    ins = model_mask(p1[0] + ts * u[0], p1[1] + ts * u[1])
    t_in = ts[ins]
    ta, tb = float(t_in.min()), float(t_in.max())
    A, B = p1 + ta * u, p1 + tb * u
    nvec = np.array([-u[1], u[0]])
    wells = []
    for w in WELLS:
        off = float(np.dot(W[w]["track"][0, :2] - A, nvec))
        if abs(off) <= maxdist:
            wells.append((w, off))
    wells.sort(key=lambda x: float(np.dot(W[x[0]]["track"][0, :2] - A, u)))
    return dict(name=name, A=A, B=B, u=u, n=nvec, length=float(tb - ta), wells=wells, w1=w1, w2=w2,
                maxdist=maxdist, why=why)


def proj(line, x, y):
    return float(np.dot(np.array([x, y]) - line["A"], line["u"]))


HOT = WS.sort_values("Tmax_degC", ascending=False).well.tolist()
LINES = [make_line("A", HOT[0], HOT[1], "through the two hottest wells by Tmax (CLAUDE.md default)")]
_cond = WS[(WS.regime == "conductive")].sort_values("Tmax_degC").well.tolist()
if _cond:
    LINES.append(make_line("B", HOT[0], _cond[0], "from the hottest well towards the conductive margin "
                           "(eastward decline of temperature)"))
_out = WS[WS.regime == "outflow"].copy()
if len(_out):
    _out["y"] = [W[w]["track"][0, 1] for w in _out.well]
    LINES.append(make_line("C", HOT[0], _out.sort_values("y").well.iloc[0],
                           "from the hottest well to the southernmost outflow well"))
LINE_ROWS = []
for L in LINES:
    LINE_ROWS.append(dict(line=L["name"], through_wells=f"{L['w1']}; {L['w2']}", reason=L["why"],
                          start_label=L["name"], start_x=round(L["A"][0], 1), start_y=round(L["A"][1], 1),
                          end_label=L["name"] + "'", end_x=round(L["B"][0], 1), end_y=round(L["B"][1], 1),
                          length_m=round(L["length"], 1), max_offset_m=L["maxdist"],
                          wells_projected="; ".join(f"{w} ({o:+.0f} m)" for w, o in L["wells"])))
    log(f"section {L['name']}: {LINE_ROWS[-1]}")
LINES_DF = pd.DataFrame(LINE_ROWS)
LINE = {L["name"]: L for L in LINES}


def surface_along(line, n=600):
    ts = np.linspace(0, line["length"], n)
    xs, ys = line["A"][0] + ts * line["u"][0], line["A"][1] + ts * line["u"][1]
    out = np.full(n, np.nan)
    for c, pth in COL_PATH.items():
        out[pth.contains_points(np.column_stack([xs, ys]))] = col_surface(c)
    return ts, out


def seg_intersections(P, Q, segs):
    out = []
    d = Q - P
    for s in segs:
        a, b = s[:-1], s[1:]
        e = b - a
        den = d[0] * e[:, 1] - d[1] * e[:, 0]
        ok = np.abs(den) > 1e-9
        dd = np.where(ok, den, 1.0)
        t = ((a[:, 0] - P[0]) * e[:, 1] - (a[:, 1] - P[1]) * e[:, 0]) / dd
        v = ((a[:, 0] - P[0]) * d[1] - (a[:, 1] - P[1]) * d[0]) / dd
        hit = ok & (t >= 0) & (t <= 1) & (v >= 0) & (v <= 1)
        out += list(t[hit] * np.linalg.norm(d))
    return sorted(out)


def line_crossings(line):
    res = []
    for name, segs in FAULTS.items():
        for s in seg_intersections(line["A"], line["B"], segs):
            res.append((s, name, "fault"))
    for name, v in WATER.items():
        segs = [np.vstack([s, s[:1]]) for s in v["segs"]] if v["fill"] else v["segs"]
        last = -1e9
        for s in seg_intersections(line["A"], line["B"], segs):
            if s - last > 300:
                res.append((s, name, "water"))
            last = s
    return sorted(res)


def well_profile_on_line(line, w):
    g = TG[TG.well == w].sort_values("elevation_mRL", ascending=False)
    s = np.array([proj(line, x, y) for x, y in zip(g.x, g.y)])
    return g, s


def envelope_grid(line, wa, wb, field, nz=150, nu=30):
    """Linear interpolation between neighbouring wells, only where both have data (no extrapolation)."""
    ga, sa = well_profile_on_line(line, wa)
    gb, sb = well_profile_on_line(line, wb)
    ztop = min(ga.elevation_mRL.max(), gb.elevation_mRL.max())
    zbot = max(ga.elevation_mRL.min(), gb.elevation_mRL.min())
    if ztop - zbot < 20:
        return None
    zz = np.linspace(ztop, zbot, nz)
    fa = np.array([interp_profile(ga.elevation_mRL, ga[field], z) for z in zz])
    fb = np.array([interp_profile(gb.elevation_mRL, gb[field], z) for z in zz])
    s_a = np.array([interp_profile(ga.elevation_mRL, sa, z) for z in zz])
    s_b = np.array([interp_profile(gb.elevation_mRL, sb, z) for z in zz])
    u = np.linspace(0, 1, nu)
    S = s_a[:, None] + (s_b - s_a)[:, None] * u[None, :]
    Z = np.repeat(zz[:, None], nu, axis=1)
    F = fa[:, None] + (fb - fa)[:, None] * u[None, :]
    return S, Z, F


PDEV_LIM = max(4.0, math.ceil(np.nanmax(np.abs(TG.P_minus_trend_bar)) / 2.0) * 2.0)


def draw_section(ax, line, mode="T", fig=None, labels=True):
    """mode 'T': temperature; 'P': deviation from the verified hydrostatic trend."""
    ts, surf = surface_along(line)
    ax.plot(ts, surf, color="saddlebrown", lw=1.5, zorder=3)
    if mode == "T":
        cmap, norm, field, levels = CMAP_T, NORM_T, "T_degC", ISO_LEVELS
    else:
        cmap, norm, field = plt.cm.RdBu_r, matplotlib.colors.TwoSlopeNorm(0, -PDEV_LIM, PDEV_LIM), "P_minus_trend_bar"
        levels = [float(v) for v in np.arange(-PDEV_LIM, PDEV_LIM + 0.1, 2.0)]
    wl = [w for w, _ in line["wells"]]
    for i in range(len(wl) - 1):
        grid = envelope_grid(line, wl[i], wl[i + 1], field)
        if grid is None:
            continue
        S, Z, F = grid
        if mode == "T":
            lw = [1.8 if lv in (100, 200) else 0.7 for lv in levels]
            cs = ax.contour(S, Z, F, levels=levels, colors="k", linewidths=lw, zorder=4)
        else:
            cs = ax.contour(S, Z, F, levels=levels, colors="0.3", linewidths=0.8, zorder=4)
        ax.clabel(cs, fmt="%g", fontsize=6)
        poly = np.vstack([np.column_stack([S[:, 0], Z[:, 0]]), np.column_stack([S[::-1, -1], Z[::-1, -1]])])
        ax.plot(np.append(poly[:, 0], poly[0, 0]), np.append(poly[:, 1], poly[0, 1]), color="#999999", lw=0.6,
                ls=":", zorder=2)
    sc = None
    for w, off in line["wells"]:
        g, s = well_profile_on_line(line, w)
        tr = W[w]["track"]
        st = np.array([proj(line, x, y) for x, y in tr[:, :2]])
        ax.plot(st, tr[:, 2], color="0.35", lw=1.0, zorder=5)
        sc = ax.scatter(s, g.elevation_mRL, c=g[field], cmap=cmap, norm=norm, s=24, edgecolors="k",
                        linewidths=0.3, zorder=7)
        shoe = W[w]["shoe"]
        ss = interp_profile(tr[:, 2], st, shoe)
        ax.plot([ss - 0.012 * line["length"], ss + 0.012 * line["length"]], [shoe, shoe], color="k", lw=2, zorder=6)
        if labels:
            ax.text(ss + 0.014 * line["length"], shoe, "shoe", fontsize=6, va="center", zorder=6)
        for top, bot, q in W[w]["fz"]:
            s1, s2 = interp_profile(tr[:, 2], st, top), interp_profile(tr[:, 2], st, bot)
            ax.plot([s1, s2], [top, bot], color="darkorange", lw=8, alpha=0.75, solid_capstyle="butt", zorder=6)
        mask = g.near_boiling_lt5C.to_numpy()
        if mask.any():
            nb = g[mask]
            z1, z2 = nb.elevation_mRL.max() + 15, nb.elevation_mRL.min() - 15
            ax.fill_betweenx([z2, z1], s[mask].min() - 0.015 * line["length"], s[mask].max() + 0.015 * line["length"],
                             facecolor="none", hatch="....", edgecolor="red", lw=0.0, zorder=5)
        ax.text(st[-1] + 0.012 * line["length"], tr[-1, 2] + 30, f"{w}\noffset {off:+.0f} m", ha="left",
                va="bottom", fontsize=7,
                color=WCOL[w], fontweight="bold", zorder=8, bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.3))
    for s, name, kind in line_crossings(line):
        colr = "k" if kind == "fault" else "#1f78b4"
        ax.plot([s, s], [Z_TOP - 0.03 * (Z_TOP - Z_BOT), Z_TOP], color=colr, lw=2.5, zorder=9)
        ax.text(s, Z_TOP + 0.01 * (Z_TOP - Z_BOT), name, rotation=90, fontsize=7, ha="center", va="bottom",
                color=colr)
    ax.set_xlim(0, line["length"])
    ax.set_ylim(Z_BOT, Z_TOP)
    ax.set_xlabel(f"Distance along section {line['name']}-{line['name']}' (m)")
    ax.set_ylabel("Elevation (mRL)")
    ax.text(0.005, 0.01, line["name"], transform=ax.transAxes, fontsize=13, fontweight="bold")
    ax.text(0.995, 0.01, line["name"] + "'", transform=ax.transAxes, fontsize=13, fontweight="bold", ha="right")
    ax.text(0.5, 0.01, "blank = no data: contours only between neighbouring wells within their common depth range",
            transform=ax.transAxes, fontsize=7, ha="center", color="0.35")
    if fig is not None and sc is not None:
        cb = fig.colorbar(sc, ax=ax, pad=0.01, fraction=0.03)
        cb.set_label("Temperature (°C)" if mode == "T" else "P minus verified hydrostatic trend (bar)")
    return sc


def vexag(ax, fig):
    fig.canvas.draw()
    bb = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    return ((x1 - x0) / bb.width) / ((y1 - y0) / bb.height)


def draw_lines_on_map(ax):
    view = ax.get_xlim() + ax.get_ylim()
    for L in LINES:
        ax.plot([L["A"][0], L["B"][0]], [L["A"][1], L["B"][1]], color="crimson", lw=1.4, ls="-.", zorder=11)
        for t_, lab in ((0.0, L["name"]), (1.0, L["name"] + "'")):
            tsq = np.linspace(0, 1, 2001) if t_ == 0 else np.linspace(1, 0, 2001)
            p = L["A"]
            for tt in tsq:
                p = L["A"] + tt * (L["B"] - L["A"])
                if _in_view(p, _shrink(view, 0.03)):
                    break
            ax.text(p[0], p[1], lab, color="crimson", fontsize=10, fontweight="bold", zorder=12,
                    ha="center", va="center", bbox=dict(fc="white", ec="crimson", pad=0.2, alpha=0.9))


# ------------------------------------------------------------------------------------------
# F1 base map
# ------------------------------------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(17, 8.5), gridspec_kw=dict(width_ratios=[1.1, 1]))
base_map(axs[0], view=MODEL_VIEW, labels=False, fs=7)
fv = FIELD_VIEW
axs[0].plot([fv[0], fv[1], fv[1], fv[0], fv[0]], [fv[2], fv[2], fv[3], fv[3], fv[2]], color="red", lw=1)
axs[0].text(fv[0], fv[3] + 150, "extent of panel (b)", color="red", fontsize=7)
axs[0].set_title("(a) Whole model: grid columns, model extent, faults, water bodies, wells", fontsize=9)
base_map(axs[1], view=FIELD_VIEW, labels=True, fs=8)
for ax in axs:
    draw_lines_on_map(ax)
axs[1].set_title("(b) Well field: wellheads, tracks and section lines", fontsize=9)
handles = [Line2D([], [], color="#c8c8c8", lw=1, label="model grid columns"),
           Line2D([], [], color="#555555", ls="--", label="model extent"),
           Line2D([], [], color="k", lw=1.6, label="fault"),
           Patch(fc="#a6d8f5", ec="#1f78b4", label="lake"), Line2D([], [], color="#1f78b4", label="river"),
           Line2D([], [], color="grey", label="well track"),
           Line2D([], [], marker="o", color="k", ls="", ms=4, label="wellhead"),
           Line2D([], [], color="crimson", ls="-.", label="section line")]
fig.legend(handles=handles, loc="lower center", ncol=8, frameon=False, bbox_to_anchor=(0.5, -0.03))
fig.suptitle("F1 Base map: model grid, mapped features, wells and section lines (plan view; easting and northing in m)",
             fontsize=11)
save(fig, os.path.join(CALFIG, "well_map.png"), "F1", "Base map: grid, model extent, faults, water bodies, "
     "wellheads, tracks, section lines", "all", "plan view; lines " + ", ".join(L["name"] for L in LINES),
     "geometry_dat; mapping; well_data_json (tracks)", "obs")


# ------------------------------------------------------------------------------------------
# F2 well profiles
# ------------------------------------------------------------------------------------------
def nb_intervals(g):
    out, cur = [], None
    for z, f in zip(g.elevation_mRL, g.near_boiling_lt5C):
        if f:
            cur = [z, z] if cur is None else [cur[0], z]
        elif cur is not None:
            out.append(cur)
            cur = None
    if cur is not None:
        out.append(cur)
    return out


GEOI = GEO.set_index("well")


def t_panel(ax, w, small=False):
    g = TG[TG.well == w].sort_values("elevation_mRL", ascending=False)
    r = W[w]
    for lv in ISO_LEVELS:
        heavy = lv in (100, 200)
        ax.axvline(lv, color="0.35" if heavy else "0.75", lw=1.1 if heavy else 0.4, zorder=0)
    ms = GEOI.loc[w, "model_surface_mRL"]
    ax.axhline(r["wellhead"], color="k", lw=0.9)
    ax.axhline(ms, color="saddlebrown", lw=0.9, ls="--")
    ax.axhline(r["shoe"], color="k", lw=1.4, ls="-.")
    for top, bot, q in r["fz"]:
        ax.axhspan(bot, top, color="darkorange", alpha=0.35, lw=0)
    for a, b in nb_intervals(g):
        ax.axhspan(b - 10, a + 10, color="red", alpha=0.18, lw=0)
    tsl = np.where(g.P_bar <= P_CRIT, g.Tsat_degC, np.nan)
    ax.plot(tsl, g.elevation_mRL, color="crimson", lw=1.0, ls="--", label="Tsat(P), IF97 (line)")
    ax.plot(g.T_degC, g.elevation_mRL, "o", color=WCOL[w], mec="k", mew=0.3, ms=3.5 if small else 5,
            label=f"T survey {r['survey_date'].isoformat()} (points)")
    for col, lv in (("z_100C_mRL", 100), ("z_200C_mRL", 200)):
        zc = WSI.loc[w, col]
        if not pd.isna(zc):
            ax.plot(lv, zc, marker="D", color="k", ms=4 if small else 6, zorder=6)
    ax.set_xlim(0, T_LIM)
    ax.set_ylim(Z_BOT, Z_TOP)
    if not small:
        xt = T_LIM * 0.99
        ax.text(xt, r["wellhead"], "wellhead", ha="right", va="top", fontsize=7)
        ax.text(T_LIM * 0.45, ms, "model surface", ha="right", va="bottom", fontsize=7, color="saddlebrown")
        ax.text(xt, r["shoe"], "shoe", ha="right", va="bottom", fontsize=7)
        for top, bot, q in r["fz"]:
            ax.text(3, 0.5 * (top + bot), f"feed zone ({q})", fontsize=7, va="center", color="saddlebrown")
        for col, lv in (("z_100C_mRL", 100), ("z_200C_mRL", 200)):
            zc = WSI.loc[w, col]
            if not pd.isna(zc):
                ax.text(lv + 4, zc + 20, f"{lv} °C at {zc:.0f} mRL", fontsize=7, va="bottom")
    ax.set_xlabel("Temperature (°C)")
    ax.grid(axis="y", lw=0.3, color="0.85")


for w in WELLS:
    r = W[w]
    g = TG[TG.well == w].sort_values("elevation_mRL", ascending=False)
    fig = plt.figure(figsize=(15.5, 8.5))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.3, 1.1], wspace=0.22)
    axm = fig.add_subplot(gs[0])
    base_map(axm, view=FIELD_VIEW, highlight=w, labels=False, fs=6, feature_labels=True)
    tr = r["track"]
    axm.annotate(w, (tr[0, 0], tr[0, 1]), xytext=(tr[0, 0] + 150, tr[0, 1] + 150), fontsize=8, color=WCOL[w],
                 fontweight="bold", zorder=12, bbox=dict(fc="white", ec=WCOL[w], pad=0.3))
    axm.set_title("Locator: this well in colour, others grey", fontsize=8)
    axt = fig.add_subplot(gs[1])
    t_panel(axt, w)
    axt.set_ylabel("Elevation (mRL)")
    axt.set_title("Temperature (°C)", fontsize=9)
    axp = fig.add_subplot(gs[2], sharey=axt)
    axp.axhline(r["wellhead"], color="k", lw=0.9)
    axp.axhline(r["shoe"], color="k", lw=1.4, ls="-.")
    for top, bot, q in r["fz"]:
        axp.axhspan(bot, top, color="darkorange", alpha=0.35, lw=0)
    zz = np.linspace(Z_BOT, Z_TOP, 50)
    axp.plot(trend(zz), zz, color="0.5", ls=":", lw=1, label="verified feed-zone P trend (line)")
    axp.plot(g.P_bar, g.elevation_mRL, "o", color=WCOL[w], mec="k", mew=0.3, ms=5, label="P survey (points)")
    f = FZ[FZ.well == w]
    axp.plot(f.P_bar, f.mid_mRL, "s", color=WCOL[w], mec="k", ms=10, label="feed-zone P")
    pv = WSI.loc[w]
    if not pd.isna(pv.pivot_elev_mRL):
        axp.plot(pv.pivot_P_bar, pv.pivot_elev_mRL, "*", color="gold", mec="k", ms=17, zorder=8,
                 label=f"pivot: {pv.pivot_P_bar:.1f} bar at {pv.pivot_elev_mRL:.0f} mRL")
    if not pd.isna(pv.water_level_mRL):
        axp.axhline(pv.water_level_mRL, color="#1f78b4", lw=1.2, ls="--",
                    label=f"water level {pv.water_level_mRL:.0f} mRL")
    axp.set_xlim(0, P_LIM)
    axp.set_xlabel("Pressure (bar)")
    axp.set_title("Pressure (bar)", fontsize=9)
    axp.grid(lw=0.3, color="0.85")
    plt.setp(axp.get_yticklabels(), visible=False)
    h1, _ = axt.get_legend_handles_labels()
    h2, _ = axp.get_legend_handles_labels()
    extra = [Patch(color="darkorange", alpha=0.35, label="feed zone (quality tag)"),
             Line2D([], [], color="k", ls="-.", lw=1.4, label=f"production casing shoe {r['shoe']:.0f} mRL"),
             Line2D([], [], color="k", lw=0.9, label=f"wellhead {r['wellhead']:.0f} mRL"),
             Line2D([], [], color="saddlebrown", ls="--", lw=0.9,
                    label=f"model surface {GEOI.loc[w, 'model_surface_mRL']:.0f} mRL"),
             Line2D([], [], marker="D", color="k", ls="", label="100 / 200 °C crossing"),
             Patch(color="red", alpha=0.18, label="boiling margin < 5 °C"),
             Line2D([], [], color="w", label=f"survey {r['survey_date'].isoformat()}; NS verdict: {verdict[w]} "
                                             f"({conf[w]} confidence)")]
    fig.legend(handles=h1 + h2 + extra, loc="lower center", ncol=4, frameon=False, fontsize=8,
               bbox_to_anchor=(0.55, -0.09))
    fig.suptitle(f"F2 {w}: temperature and pressure vs elevation ({Z_TOP} to {Z_BOT} mRL; °C, bar). "
                 f"Regime (rule): {WSI.loc[w, 'regime']}", fontsize=11)
    save(fig, os.path.join(CALFIG, f"profile_{w}.png"), "F2", f"{w} T and P profiles with Tsat, shoe, feed zones, "
         "surface, isotherm crossings, pivot and water level", w, f"{Z_TOP} to {Z_BOT} mRL",
         "well_data_json; IF97; geometry_dat (surface)", "obs")

for gname, gw in GROUP_WELLS.items():
    n = len(gw)
    ncol = min(5, n)
    nrow = math.ceil(n / ncol)
    fig = plt.figure(figsize=(3.0 * ncol + 5, 5.2 * nrow))
    gs = fig.add_gridspec(nrow, ncol + 1, width_ratios=[1.8] + [1] * ncol, wspace=0.12, hspace=0.28)
    axm = fig.add_subplot(gs[:, 0])
    base_map(axm, view=FIELD_VIEW, labels=False, fs=6)
    for w in gw:
        tr = W[w]["track"]
        axm.plot(tr[:, 0], tr[:, 1], color=WCOL[w], lw=2.2, zorder=9)
        axm.plot(tr[0, 0], tr[0, 1], "o", color=WCOL[w], mec="k", ms=7, zorder=10)
        dx, dy = LABEL_OFF[w]
        axm.annotate(w, (tr[0, 0], tr[0, 1]), xytext=(tr[0, 0] + dx, tr[0, 1] + dy), fontsize=7,
                     color=WCOL[w], fontweight="bold", zorder=11, bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.2))
    axm.set_title(f"Group '{gname}': locator (wells in their colours)", fontsize=8)
    first = None
    for i, w in enumerate(gw):
        ax = fig.add_subplot(gs[i // ncol, 1 + i % ncol], sharey=first)
        first = first or ax
        t_panel(ax, w, small=True)
        ax.set_title(f"{w} ({conf[w]} confidence)", fontsize=8, color=WCOL[w], fontweight="bold")
        if i % ncol:
            plt.setp(ax.get_yticklabels(), visible=False)
        else:
            ax.set_ylabel("Elevation (mRL)")
    extra = [Line2D([], [], marker="o", color="0.5", ls="", label="T survey (points, well colour)"),
             Line2D([], [], color="crimson", ls="--", label="Tsat(P)"),
             Patch(color="darkorange", alpha=0.35, label="feed zone"),
             Line2D([], [], color="k", ls="-.", label="casing shoe"),
             Line2D([], [], color="k", lw=0.9, label="wellhead"),
             Line2D([], [], color="saddlebrown", ls="--", label="model surface"),
             Line2D([], [], marker="D", color="k", ls="", label="100/200 °C crossing"),
             Patch(color="red", alpha=0.18, label="boiling margin < 5 °C")]
    fig.legend(handles=extra, loc="lower center", ncol=8, frameon=False, bbox_to_anchor=(0.55, -0.03))
    fig.suptitle(f"F2 group '{gname}': temperature profiles on identical axes (0 to {T_LIM:.0f} °C, {Z_TOP} to {Z_BOT} mRL)",
                 fontsize=11)
    save(fig, os.path.join(CALFIG, f"profile_panels_{gname}.png"), "F2", f"Small-multiple temperature profiles, group {gname}",
         "; ".join(gw), f"{Z_TOP} to {Z_BOT} mRL", "well_data_json; IF97", "obs")


# ------------------------------------------------------------------------------------------
# F3 elevation-band symbol maps
# ------------------------------------------------------------------------------------------
def band_maps(kind):
    nb = len(BANDS)
    ncol = 4
    nrow = math.ceil((nb + 1) / ncol)
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.7 * ncol, 5.0 * nrow))
    axs = axs.ravel()
    col, step, unit = ("T_minus_field_median_degC", 10.0, "°C") if kind == "T" else ("P_minus_trend_bar", 2.0, "bar")
    sizes = [30, 70, 130, 210]
    for k, bt in enumerate(BANDS):
        ax = axs[k]
        base_map(ax, view=FIELD_VIEW, labels=True, fs=5, feature_labels=(k == 0))
        sub = BA[BA.band_top_mRL == bt]
        zmid = bt - 250
        for w in WELLS:
            tr = W[w]["track"]
            zq = min(max(zmid, tr[:, 2].min()), tr[:, 2].max())
            x, y = track_xy(w, zq)
            row = sub[sub.well == w]
            if row.empty:
                ax.plot(x, y, "o", mfc="none", mec="0.45", ms=6, zorder=12)
                continue
            d = float(row[col].iloc[0])
            if abs(d) < step:
                ax.plot(x, y, "x", color="k", ms=8, mew=1.8, zorder=12)
            else:
                cls = min(int(abs(d) // step), 4)
                ax.scatter(x, y, s=sizes[cls - 1], color="#d7301f" if d > 0 else "#2b6cb0", edgecolors="k",
                           linewidths=0.5, zorder=12)
        if kind == "T":
            ax.set_title(f"{bt} to {bt - 500} mRL: field median {sub.T_field_median_degC.iloc[0]:.1f} °C", fontsize=8)
        else:
            ax.set_title(f"{bt} to {bt - 500} mRL: P minus verified trend", fontsize=8)
    for ax in axs[nb:]:
        ax.axis("off")
    lab = [f"{step * i:g} to {step * (i + 1):g} {unit}" for i in range(1, 4)] + [f"> {step * 4:g} {unit}"]
    hs = [Line2D([], [], marker="x", color="k", ls="", ms=8, mew=1.8, label=f"within ±{step:g} {unit}")]
    hot_lab, cold_lab = ("hotter", "colder") if kind == "T" else ("higher", "lower")
    for i in range(4):
        hs.append(Line2D([], [], marker="o", color="#d7301f", mec="k", ls="", ms=math.sqrt(sizes[i]) * 0.9,
                         label=f"{hot_lab}: {lab[i]}"))
    for i in range(4):
        hs.append(Line2D([], [], marker="o", color="#2b6cb0", mec="k", ls="", ms=math.sqrt(sizes[i]) * 0.9,
                         label=f"{cold_lab}: {lab[i]}"))
    hs.append(Line2D([], [], marker="o", mfc="none", mec="0.45", ls="", ms=6, label="no data in band"))
    axs[nb].legend(handles=hs, loc="center", fontsize=8, frameon=False,
                   title=("Well mean T minus band field median" if kind == "T" else
                          "Well mean P minus verified hydrostatic trend"))
    ttl = ("F3 Temperature: well mean minus band field median in 500 m elevation bands (plan view; °C)" if kind == "T"
           else "F3 Pressure: well mean minus hydrostatic trend of gauge-verified feed zones, 500 m bands (plan view; bar)")
    fig.suptitle(ttl, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


save(band_maps("T"), os.path.join(CALFIG, "maps_temperature_bands.png"), "F3",
     "Well mean T minus field median per 500 m band", "all", f"bands {BANDS[0]} to {BANDS[-1] - 500} mRL",
     "T_P_at_elevation_bands.csv", "obs")
save(band_maps("P"), os.path.join(CALFIG, "maps_pressure_bands.png"), "F3",
     "Well mean P minus verified hydrostatic trend per 500 m band", "all", f"bands {BANDS[0]} to {BANDS[-1] - 500} mRL",
     "T_P_at_elevation_bands.csv", "mixed")

# ------------------------------------------------------------------------------------------
# F4 sections
# ------------------------------------------------------------------------------------------
VE = {}
for L in LINES:
    for mode in ("T", "P"):
        fig, ax = plt.subplots(figsize=(15, 7.5))
        draw_section(ax, L, mode=mode, fig=fig)
        ve = vexag(ax, fig)
        VE[(L["name"], mode)] = round(ve, 1)
        hs = [Line2D([], [], color="saddlebrown", lw=1.5, label="model surface"),
              Line2D([], [], color="0.35", lw=1, label="well track (projected)"),
              Line2D([], [], color="k", lw=2, label="production casing shoe"),
              Line2D([], [], color="darkorange", lw=6, alpha=0.75, label="feed zone"),
              Patch(fc="none", ec="red", hatch="....", label="boiling margin < 5 °C"),
              Line2D([], [], color="#999999", ls=":", label="data envelope"),
              Line2D([], [], color="k", lw=2.5, label="fault crossing (tick, top)"),
              Line2D([], [], color="#1f78b4", lw=2.5, label="water-body crossing (tick, top)")]
        if mode == "T":
            hs[:0] = [Line2D([], [], color="k", lw=1.8, label="100 / 200 °C isotherm"),
                      Line2D([], [], color="k", lw=0.7, label="50 °C isotherms")]
        else:
            hs[:0] = [Line2D([], [], color="0.3", lw=0.8, label="contours of P minus trend, 2 bar")]
        ax.legend(handles=hs, loc="upper left", bbox_to_anchor=(1.08, 1.0), fontsize=8, frameon=False)
        what = "temperature (°C)" if mode == "T" else "pressure minus verified hydrostatic trend (bar)"
        ax.set_title(f"F4 Section {L['name']}-{L['name']}': {what}; wells within {L['maxdist']:.0f} m projected; "
                     f"vertical exaggeration {ve:.1f}", fontsize=10, pad=75)
        fn = f"section_{L['name']}.png" if mode == "T" else f"section_{L['name']}_P.png"
        save(fig, os.path.join(CALFIG, fn), "F4", f"Section {L['name']}-{L['name']}': {what}",
             "; ".join(w for w, _ in L["wells"]), f"line {L['name']}-{L['name']}' ({L['A'][0]:.0f}, {L['A'][1]:.0f}) "
             f"to ({L['B'][0]:.0f}, {L['B'][1]:.0f})", "NS_calibration_targets.csv; geometry_dat; mapping",
             "obs" if mode == "T" else "mixed")
LINES_DF["vertical_exaggeration_T"] = [VE[(L["name"], "T")] for L in LINES]
write_csv(LINES_DF, os.path.join(ADD, "section_lines.csv"))

# ------------------------------------------------------------------------------------------
# F5 P-T state
# ------------------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(10, 7.5))
tt = np.linspace(0.01, T_CRIT - 1e-6, 400)
ax.plot(tt, [psat(t) for t in tt], color="crimson", lw=1.5, label="IF97 saturation curve (to critical point)")
ax.plot(T_CRIT, P_CRIT, "*", color="crimson", ms=14, mec="k", label=f"critical point ({T_CRIT} °C, {P_CRIT} bar)")
for w in WELLS:
    g = TG[TG.well == w]
    ax.plot(g.T_degC, g.P_bar, "o", color=WCOL[w], ms=4, mec="none", label=w)
    f = FZ[FZ.well == w]
    ax.plot(f.T_degC, f.P_bar, "s", color=WCOL[w], ms=11, mec="k", mew=1.2)
ax.plot([], [], "s", color="w", mec="k", ms=11, label="feed zone (large square)")
ax.set_xlim(0, 400)
ax.set_ylim(0, P_LIM)
ax.set_xlabel("Temperature (°C)")
ax.set_ylabel("Pressure (bar)")
ax.grid(lw=0.3, color="0.85")
ax.text(300, 30, "steam side of the curve", fontsize=8, color="crimson")
ax.text(40, 200, "compressed liquid", fontsize=8, color="0.3")
ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8, frameon=False)
ax.set_title("F5 P-T state of all survey points against the IF97 saturation curve (all elevations; °C, bar)")
save(fig, os.path.join(CALFIG, "PT_state.png"), "F5", "All survey points and feed zones on P-T axes with the "
     "saturation curve and critical point", "all", "all elevations", "NS_calibration_targets.csv; feedzones.csv; IF97", "obs")

# ------------------------------------------------------------------------------------------
# F7 gauge check
# ------------------------------------------------------------------------------------------
GW = [w for w in WELLS if w in gauge_rows]
GAUGE_PERIOD = None
if GW:
    t0 = dt.date(min(min(g["series"]) for g in gauge_rows.values()).year - 1, 7, 1)
    t1 = dt.date(max(max(g["series"]) for g in gauge_rows.values()).year, 7, 1)
    GAUGE_PERIOD = (t0, t1)
    fig, axs = plt.subplots(1, len(GW), figsize=(4.9 * len(GW), 5.4))
    axs = np.atleast_1d(axs)
    for ax, w in zip(axs, GW):
        g = gauge_rows[w]
        ds = sorted(g["series"])
        key = "Injection_Mass_Flow" if W[w]["injector"] else "Production_Mass_Flow"
        fl = sorted(W[w]["flows"].get(key, {}).items())
        for (d1, v1), (d2, _) in zip(fl[:-1], fl[1:]):
            if abs(v1) > 0:
                ax.axvspan(d1, d2, color="0.87", lw=0, zorder=0)
        ax.axhspan(g["first_gauge_P_bar"] - g["precision"], g["first_gauge_P_bar"] + g["precision"],
                   color="green", alpha=0.3, label=f"first reading ± {g['precision']:.1f} bar")
        ax.plot(ds, [g["series"][d] for d in ds], "o", color=WCOL[w], mec="k", ms=6, label="gauge reading (points)",
                zorder=5)
        ax.axhline(g["survey_P_at_gauge_bar"], color="k", lw=1.3, ls="--",
                   label=f"survey P at gauge, survey {W[w]['survey_date'].isoformat()}: {g['survey_P_at_gauge_bar']:.2f} bar")
        if FIRST_FLOW[w]:
            ax.axvline(FIRST_FLOW[w], color="purple", lw=1.6, label=f"own first flow {FIRST_FLOW[w].isoformat()}")
        ax.axvline(FIELD_START, color="red", lw=1.2, ls=":", label=f"first field production {FIELD_START.isoformat()}")
        ax.set_xlim(t0, t1)
        ax.set_xlabel("Date")
        ax.set_ylabel(f"Pressure at {g['gauge_elev_mRL']:.0f} mRL (bar)")
        ax.set_title(f"{w}: {verdict[w]} ({conf[w]})", fontsize=9)
        ax.legend(fontsize=6.5, loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False)
        for lab in ax.get_xticklabels():
            lab.set_rotation(30)
    fig.suptitle(f"F7 Gauge check at the gauge elevation, {t0.isoformat()} to {t1.isoformat()} "
                 "(grey shading = well flowing; bar)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, os.path.join(CALFIG, "gauge_check.png"), "F7", "Gauge transient pressure vs survey P at the gauge "
         "elevation, precision band and flow periods", "; ".join(GW), "gauge elevations " +
         ", ".join(f"{w} {gauge_rows[w]['gauge_elev_mRL']:.0f} mRL" for w in GW),
         "well_data_json (transient pressure, flows); pressure_consistency_check.csv", "obs")

log("standard figures F1-F5 and F7 done")

# ==========================================================================================
# F6 CONCEPTUAL MODEL (data only)
# ==========================================================================================
def fmt(v, d=1, unit=""):
    if v is None or (isinstance(v, float) and math.isnan(v)) or (not isinstance(v, str) and pd.isna(v)):
        return "blank"
    return f"{v:.{d}f}{unit}"


# fault intersections (upflow location candidates)
FAULT_X = []
fn = list(FAULTS)
for i in range(len(fn)):
    for j in range(i + 1, len(fn)):
        a = FAULTS[fn[i]][0]
        for k in range(len(a) - 1):
            P, Q = a[k], a[k + 1]
            for s_ in seg_intersections(P, Q, FAULTS[fn[j]]):
                pt = P + (Q - P) / np.linalg.norm(Q - P) * s_
                FAULT_X.append((fn[i], fn[j], round(float(pt[0]), 1), round(float(pt[1]), 1)))
log("fault intersections:", FAULT_X)

UPW = WS[WS.regime == "upflow"].well.tolist() or [HOT[0]]
UP_WELLS = list(dict.fromkeys(UPW + [HOT[0], HOT[1]]))
UPC = np.mean([W[w]["track"][0, :2] for w in UP_WELLS], axis=0)
near_x = None
if FAULT_X:
    near_x = min(FAULT_X, key=lambda f: math.hypot(f[2] - UPC[0], f[3] - UPC[1]))
    near_x = near_x + (round(math.hypot(near_x[2] - UPC[0], near_x[3] - UPC[1]), 1),)

# pressure connection at a common deep elevation (open hole only)
Z_CONN = -1300.0
pconn = {}
for w in WELLS:
    s = W[w]["survey"]
    if s.z.min() <= Z_CONN <= s.z.max() and Z_CONN < W[w]["shoe"]:
        pconn[w] = round(interp_profile(s.z, s.P, Z_CONN), 1)
cw = list(pconn)
edges = []
if len(cw) >= 2:
    inn = {cw[0]}
    while len(inn) < len(cw):
        best = None
        for a in inn:
            for b in cw:
                if b in inn:
                    continue
                d = math.dist(W[a]["track"][0, :2], W[b]["track"][0, :2])
                if best is None or d < best[0]:
                    best = (d, a, b)
        edges.append(best)
        inn.add(best[2])
PCONN = pd.DataFrame([dict(well_a=a, well_b=b, elevation_mRL=Z_CONN, P_a_bar=pconn[a], P_b_bar=pconn[b],
                           dP_b_minus_a_bar=round(pconn[b] - pconn[a], 1), distance_m=round(d, 1))
                      for d, a, b in edges])
write_csv(PCONN, os.path.join(ADD, "pressure_connections.csv"))

LB = LINE.get("B", LINES[0])
on_B = [w for w, _ in LB["wells"]]
EL = []


def add_el(element, symbol, location, observation, claim, style):
    EL.append(dict(no=len(EL) + 1, element=element, symbol=symbol, location=location, observation=observation,
                   claim_type=claim, line_style=style))
    return len(EL)


up_obs = "; ".join(f"{w} Tmax {WSI.loc[w, 'Tmax_degC']:.1f} °C at {WSI.loc[w, 'Tmax_elev_mRL']:.0f} mRL, "
                   f"200 °C crossing {fmt(WSI.loc[w, 'z_200C_mRL'])} mRL, regime {WSI.loc[w, 'regime']}" for w in UP_WELLS)
loc_up = f"beneath {', '.join(UP_WELLS)}"
if near_x:
    loc_up += f"; {near_x[0]} / {near_x[1]} intersection ({near_x[2]:.0f}, {near_x[3]:.0f}) is {near_x[4]:.0f} m from their centre"
N_UP = add_el("Deep hot upflow", "red upward arrow", loc_up, up_obs, "interp", "dashed (inferred)")
OUT_WELLS = WS[WS.regime == "outflow"].sort_values("Tmax_degC", ascending=False).well.tolist()
N_OUT = {}
for w in OUT_WELLS:
    r_ = WSI.loc[w]
    N_OUT[w] = add_el(f"Lateral outflow towards {w}", "orange horizontal arrow",
                      f"from the upflow to {w} at {r_.inversion_top_mRL:.0f} mRL",
                      f"{w}: inversion {r_.inversion_degC:.1f} °C, peak {r_.Tmax_degC:.1f} °C at {r_.inversion_top_mRL:.0f} mRL, "
                      f"base {r_.inversion_base_mRL:.0f} mRL", "interp",
                      "solid in the section where the well is on the line; dashed path in plan")
REC_WELLS = WS[WS.regime == "recharge/downflow"].well.tolist()
N_REC = {}
for w in REC_WELLS:
    r_ = WSI.loc[w]
    N_REC[w] = add_el(f"Cold recharge or downflow at {w}", "blue arrow", f"{w} (field margin)",
                      f"{w}: Tmax {r_.Tmax_degC:.1f} °C at {r_.Tmax_elev_mRL:.0f} mRL, shallow gradient "
                      f"{fmt(r_.shallow_gradient_degC_per_km)} °C/km" +
                      (f", inversion {r_.inversion_degC:.1f} °C {r_.inversion_top_mRL:.0f} to {r_.inversion_base_mRL:.0f} mRL"
                       if not pd.isna(r_.inversion_degC) else ""), "interp", "dashed (inferred)")
deep_cool = [w for w in OUT_WELLS if not pd.isna(WSI.loc[w, "inversion_base_mRL"]) and WSI.loc[w, "inversion_degC"] >= 30]
N_DEEP = add_el("Cool water beneath the outflow (deep lateral inflow)", "blue arrow",
                "below the inversion bases, " + ", ".join(deep_cool),
                "; ".join(f"{w} falls to {WSI.loc[w, 'Tmax_degC'] - WSI.loc[w, 'inversion_degC']:.1f} °C at "
                          f"{WSI.loc[w, 'inversion_base_mRL']:.0f} mRL" for w in deep_cool), "interp", "dashed (inferred)")
CAPW = WS[WS.cap_base_mRL.notna()]
N_CAP = add_el("Cap (low permeability, conductive gradient)", "grey band",
               f"{', '.join(CAPW.well)}: base {CAPW.cap_base_mRL.max():.0f} to {CAPW.cap_base_mRL.min():.0f} mRL",
               "shallow gradients " + "; ".join(f"{w} {WSI.loc[w, 'shallow_gradient_degC_per_km']:.1f} °C/km"
                                                for w in CAPW.well), "interp", "outlined")
NBW = WS[WS.near_boiling].well.tolist()
N_NB = add_el("Near-boiling zone", "stipple", ", ".join(NBW) or "none observed",
              "; ".join(f"{w} minimum boiling margin {WSI.loc[w, 'min_boiling_margin_degC']:.2f} °C at "
                        f"{WSI.loc[w, 'min_boiling_margin_elev_mRL']:.0f} mRL" for w in NBW) or "no margin below 5 °C",
              "obs", "stipple")
N_ISO = add_el("Isotherms 100 °C and 200 °C", "black lines (from F4)", f"section {LB['name']}-{LB['name']}'",
               "interpolated between neighbouring wells inside the data envelope", "interp", "solid")
N_PC = add_el("Pressure connection between wells", "thin dashed line, labelled dP",
              f"open-hole pressures at {Z_CONN:.0f} mRL",
              "; ".join(f"{r.well_a}-{r.well_b} {r.dP_b_minus_a_bar:+.1f} bar" for r in PCONN.itertuples()),
              "mixed", "thin dashed")
ELDF = pd.DataFrame(EL)
write_csv(ELDF, os.path.join(ADD, "conceptual_model_elements.csv"))


def num_label(ax, x, y, n, color="k"):
    ax.text(x, y, str(n), fontsize=8, fontweight="bold", color="white", ha="center", va="center", zorder=30,
            bbox=dict(boxstyle="circle,pad=0.25", fc=color, ec="k", lw=0.6))


fig = plt.figure(figsize=(20, 9.5))
gsx = fig.add_gridspec(1, 2, width_ratios=[1, 1.45], wspace=0.12)
axp = fig.add_subplot(gsx[0])
base_map(axp, view=FIELD_VIEW, labels=True, fs=7)
draw_lines_on_map(axp)
# cap outline (convex hull of cap wells)
if len(CAPW) >= 3:
    pts = np.array([W[w]["track"][0, :2] for w in CAPW.well])
    c0 = pts.mean(0)
    ang = np.arctan2(pts[:, 1] - c0[1], pts[:, 0] - c0[0])
    hull = pts[np.argsort(ang)]
    hull = c0 + (hull - c0) * 1.25
    axp.fill(hull[:, 0], hull[:, 1], fc="0.6", alpha=0.18, ec="0.3", ls="--", lw=1.0, zorder=5)
    num_label(axp, hull[:, 0].min() + 150, c0[1], N_CAP, "0.45")
# near-boiling
for w in NBW:
    x0, y0 = W[w]["track"][0, :2]
    axp.add_patch(matplotlib.patches.Circle((x0, y0), 260, fc="none", ec="red", hatch="....", lw=0.5, zorder=13))
    num_label(axp, x0 - 330, y0 - 200, N_NB, "red")
# upflow
axp.plot(UPC[0], UPC[1], marker="^", ms=26, mfc="red", mec="k", mew=1.0, ls="", zorder=14, alpha=0.85)
axp.add_patch(matplotlib.patches.Circle(UPC, 420, fc="none", ec="red", ls="--", lw=1.5, zorder=14))
num_label(axp, UPC[0] - 520, UPC[1] + 330, N_UP, "red")
# outflow
for w in OUT_WELLS:
    xy = W[w]["track"][0, :2]
    if math.dist(xy, UPC) < 300:
        num_label(axp, xy[0] + 250, xy[1] - 250, N_OUT[w], "darkorange")
        continue
    axp.annotate("", xy=xy, xytext=UPC, arrowprops=dict(arrowstyle="-|>", color="darkorange", lw=2.2, ls="--",
                                                       shrinkA=18, shrinkB=8, mutation_scale=18), zorder=15)
    m = UPC + 0.6 * (xy - UPC)
    num_label(axp, m[0], m[1], N_OUT[w], "darkorange")
    axp.text(m[0] + 120, m[1] - 160, f"{WSI.loc[w, 'inversion_top_mRL']:.0f} mRL", fontsize=6.5, color="darkorange",
             zorder=16, bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.2))
# recharge
for w in REC_WELLS:
    xy = W[w]["track"][0, :2]
    v = (xy - UPC) / np.linalg.norm(xy - UPC)
    start = xy + v * 900
    axp.annotate("", xy=xy, xytext=start, arrowprops=dict(arrowstyle="-|>", color="#2b6cb0", lw=2.4, ls="--",
                                                         shrinkB=6, mutation_scale=18), zorder=15)
    num_label(axp, start[0], start[1], N_REC[w], "#2b6cb0")
# pressure connections
for r_ in PCONN.itertuples():
    a, b = W[r_.well_a]["track"][0, :2], W[r_.well_b]["track"][0, :2]
    axp.plot([a[0], b[0]], [a[1], b[1]], color="purple", lw=0.8, ls=(0, (4, 3)), zorder=12)
    m = 0.5 * (a + b)
    axp.text(m[0], m[1], f"{r_.dP_b_minus_a_bar:+.1f}", fontsize=6, color="purple", zorder=16,
             bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.1))
if len(PCONN):
    a = W[PCONN.well_a.iloc[0]]["track"][0, :2]
    num_label(axp, a[0] - 150, a[1] - 350, N_PC, "purple")
axp.set_title(f"(a) Plan view: numbered elements (dP labels: bar at {Z_CONN:.0f} mRL, second well minus first)",
              fontsize=9)

axs_ = fig.add_subplot(gsx[1])
draw_section(axs_, LB, mode="T", fig=fig)
sx = {w: proj(LB, *W[w]["track"][0, :2]) for w, _ in LB["wells"]}
# upflow arrow at the upflow wells on the line
for w in [u for u in UP_WELLS if u in sx]:
    zt = WSI.loc[w, "TD_mRL"]
    axs_.annotate("", xy=(sx[w] - 250, WSI.loc[w, "Tmax_elev_mRL"] + 150), xytext=(sx[w] - 250, zt - 400),
                  arrowprops=dict(arrowstyle="-|>", color="red", lw=3, ls="--", mutation_scale=22), zorder=20)
    num_label(axs_, sx[w] - 450, zt - 300, N_UP, "red")
# outflow arrows (solid: inversion measured in a well on the line)
prev = [u for u in UP_WELLS if u in sx]
src = sx[prev[0]] if prev else 0.0
for w in sorted([w for w in OUT_WELLS if w in sx], key=lambda q: sx[q]):
    zi = WSI.loc[w, "inversion_top_mRL"]
    axs_.annotate("", xy=(sx[w] + 300, zi), xytext=(src, zi),
                  arrowprops=dict(arrowstyle="-|>", color="darkorange", lw=3, mutation_scale=20), zorder=20)
    num_label(axs_, 0.5 * (src + sx[w]), zi + 110, N_OUT[w], "darkorange")
    src = sx[w]
# deep cool inflow
dc = [w for w in deep_cool if w in sx]
if dc:
    w = dc[0]
    zb = WSI.loc[w, "inversion_base_mRL"]
    farther = [q for q, _ in LB["wells"] if sx[q] > sx[w]]
    xs0 = sx[farther[-1]] if farther else sx[w] + 2000
    axs_.annotate("", xy=(sx[w] + 150, zb - 150), xytext=(xs0, zb - 150),
                  arrowprops=dict(arrowstyle="-|>", color="#2b6cb0", lw=3, ls="--", mutation_scale=20), zorder=20)
    num_label(axs_, 0.5 * (xs0 + sx[w]), zb - 30, N_DEEP, "#2b6cb0")
# cap band between cap wells on the line
cw_line = [w for w in CAPW.well if w in sx]
if cw_line:
    xs_ = sorted(sx[w] for w in cw_line)
    xa, xb = xs_[0] - 400, xs_[-1] + 400
    tops = [WSI.loc[w, "wellhead_mRL"] for w in cw_line]
    bases = [WSI.loc[w, "cap_base_mRL"] for w in cw_line]
    xq = [sx[w] for w in cw_line]
    o = np.argsort(xq)
    xq, tops, bases = np.array(xq)[o], np.array(tops)[o], np.array(bases)[o]
    xx = np.linspace(xa, xb, 50)
    axs_.fill_between(xx, np.interp(xx, xq, bases), np.interp(xx, xq, tops), color="0.5", alpha=0.25, zorder=1)
    for w in cw_line:
        right = sx[w] == max(sx[q] for q in cw_line)
        xt_, ha_ = (sx[w] + 450, "left") if right else (sx[w] - 450, "right")
        axs_.annotate(f"cap base {WSI.loc[w, 'cap_base_mRL']:.0f} mRL", (sx[w], WSI.loc[w, "cap_base_mRL"]),
                      xytext=(xt_, WSI.loc[w, "cap_base_mRL"] - 60), fontsize=6.5, ha=ha_, color="0.2", zorder=21,
                      arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6),
                      bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.2))
        axs_.annotate(f"cap top {WSI.loc[w, 'wellhead_mRL']:.0f} mRL", (sx[w], WSI.loc[w, "wellhead_mRL"]),
                      xytext=(xt_, WSI.loc[w, "wellhead_mRL"] + 60), fontsize=6.5, ha=ha_, color="0.2", zorder=21,
                      arrowprops=dict(arrowstyle="-", color="0.4", lw=0.6),
                      bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.2))
    num_label(axs_, xa - 250, 0.5 * (tops[0] + bases[0]), N_CAP, "0.45")
for w in NBW:
    if w in sx:
        num_label(axs_, sx[w] + 350, WSI.loc[w, "min_boiling_margin_elev_mRL"], N_NB, "red")
# pressure connection along the line
pl = [w for w, _ in LB["wells"] if w in pconn]
for a, b in zip(pl[:-1], pl[1:]):
    axs_.plot([sx[a], sx[b]], [Z_CONN, Z_CONN], color="purple", lw=0.9, ls=(0, (4, 3)), zorder=19)
    axs_.text(0.5 * (sx[a] + sx[b]), Z_CONN + 30, f"dP {pconn[b] - pconn[a]:+.1f} bar", fontsize=7, color="purple",
              ha="center", zorder=21, bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.2))
if len(pl) >= 2:
    num_label(axs_, sx[pl[0]] + 250, Z_CONN - 90, N_PC, "purple")
num_label(axs_, 0.97 * LB["length"], 0.6 * Z_BOT, N_ISO, "k")
axs_.set_title(f"(b) Section {LB['name']}-{LB['name']}' (as F4) with F4 isotherms; vertical exaggeration "
               f"{VE[(LB['name'], 'T')]:.1f}", fontsize=9, pad=75)
leg = [Line2D([], [], marker="^", color="red", ls="", ms=12, label="deep hot upflow (red upward arrow)"),
       Line2D([], [], color="darkorange", lw=2.5, label="lateral outflow at its depth (orange)"),
       Line2D([], [], color="#2b6cb0", lw=2.5, ls="--", label="cold recharge / downflow (blue)"),
       Patch(fc="0.5", alpha=0.3, ec="0.3", ls="--", label="cap (grey band, top/base labelled)"),
       Patch(fc="none", ec="red", hatch="....", label="two-phase / near-boiling (stipple)"),
       Line2D([], [], color="k", lw=1.8, label="isotherms 100 / 200 °C (from F4)"),
       Line2D([], [], color="purple", lw=0.9, ls=(0, (4, 3)), label="pressure connection, dP at common elevation"),
       Line2D([], [], color="0.3", lw=2, label="solid = measured feature at that well and depth"),
       Line2D([], [], color="0.3", lw=2, ls="--", label="dashed = inferred"),
       Line2D([], [], marker="o", color="0.3", ls="", ms=9, label="numbered: see report section 5 table")]
fig.legend(handles=leg, loc="lower center", ncol=5, frameon=False, fontsize=8, bbox_to_anchor=(0.5, -0.06))
fig.suptitle("F6 Conceptual model from the data alone (plan view and section "
             f"{LB['name']}-{LB['name']}'; elevations in mRL, temperatures in °C)", fontsize=11)
save(fig, os.path.join(CALFIG, "conceptual_model.png"), "F6", "Conceptual model: upflow, outflow, recharge, cap, "
     "near-boiling zone, isotherms, pressure connections (numbered)", "all",
     f"plan view and line {LB['name']}-{LB['name']}'", "well_summary.csv; conceptual_model_elements.csv; "
     "pressure_connections.csv; mapping", "interp")

FIDX = pd.DataFrame(FIG_INDEX)
FIDX["_o"] = FIDX.figure_id.str[1:].astype(int)
FIDX = FIDX.sort_values(["_o"], kind="stable").drop(columns="_o")
write_csv(FIDX, os.path.join(CAL, "figures_index.csv"))

# ==========================================================================================
# PARAMETERS USED (so that every number in the report is traceable to an output table)
# ==========================================================================================
H_UP = h_liq(WSI.loc[HOT[0], "Tmax_degC"], max(float(TG[TG.well == HOT[0]].P_bar.max()),
                                                psat(WSI.loc[HOT[0], "Tmax_degC"]) + 1))
REFBC = {}
for p in inv_paths("reference_bc"):
    b = json.load(open(p))
    REFBC["heat_flux_default_Wm2"] = b.get("Heat", {}).get("Default", {}).get("Heat_Flux", np.nan)
    REFBC["rain_default_mm_per_yr"] = b.get("Rain", {}).get("Default", {}).get("Rainfall", np.nan)
    REFBC["infiltration_fraction_default"] = b.get("Rain", {}).get("Default", {}).get("Infiltration_Fraction", np.nan)
    REFBC["atmosphere_T_default_degC"] = b.get("Atmosphere", {}).get("Default", {}).get("Temperature", np.nan)
    REFBC["n_rain_columns"] = len(b.get("Rain", {}).get("Columns", {}))
    REFBC["n_lake_columns"] = sum(len(z.get("Columns", {})) for z in b.get("Wet_Columns", {}).get("Zones", {}).values())
params = [
    ("IF97 Tsat at 1.01325 bar (degC)", round(tsat(P_ATM), 2), "validation"),
    ("IF97 Tsat at 100 bar (degC)", round(tsat(100.0), 2), "validation"),
    ("reference Tsat 1.01325 bar (degC)", 99.97, "CLAUDE.md"), ("reference Tsat 100 bar (degC)", 311.0, "CLAUDE.md"),
    ("critical pressure (bar)", P_CRIT, "IF97"), ("critical temperature (degC)", T_CRIT, "IF97"),
    ("atmospheric pressure (bar)", P_ATM, "definition"), ("g (m/s2)", G, "CLAUDE.md"),
    ("hydrostatic window minimum length (m)", 300, "CLAUDE.md 4.3"),
    ("isothermal section minimum span (m)", 300, "CLAUDE.md 4.2"), ("isothermal section T range (degC)", 5, "CLAUDE.md 4.2"),
    ("near-wellhead distance (m)", 30, "CLAUDE.md 5.5"), ("elevation band thickness (m)", 500, "CLAUDE.md 8.2"),
    ("section maximum offset (m)", 500, "CLAUDE.md 8.2"), ("verified-trend tolerance floor (bar)", 3.0, "CLAUDE.md 4.1"),
    ("verified trend intercept (bar)", round(TREND_A, 3), "feed zones of gauge-verified wells"),
    ("verified trend slope (bar/m, P increasing downward)", round(-TREND_B, 5), "feed zones of gauge-verified wells"),
    ("verified trend implied density (kg/m3)", round(-TREND_B * 1e5 / G, 1), "slope / g"),
    ("verified trend RMS residual (bar)", round(RMS_V, 3), "feed zones of gauge-verified wells"),
    ("tolerance used (bar)", round(TOL, 2), "max(3 bar, 2 x RMS)"),
    ("gauge reading precision (bar)", 0.1, "decimals of gauge readings"),
    ("injectate enthalpy, median (kJ/kg)", round(H_INJ, 1), "Injection_Enthalpy records"),
    ("injectate temperature (degC)", round(T_INJ, 1), "IF97 Region 1 inverse at 1.01325 bar"),
    ("cold-well component, cased (ratio)", round(cold_comp.get("cased", np.nan), 4), "median of cold-well windows"),
    ("cold-well component, open (ratio)", round(cold_comp.get("open", np.nan), 4), "median of cold-well windows"),
    ("cap-base gradient threshold (degC/km)", 50, "this analysis [assumption]"),
    ("pressure-connection elevation (mRL)", Z_CONN, "this analysis"),
    ("liquid enthalpy at hottest-well Tmax (kJ/kg)", round(H_UP, 1), f"IF97 Region 1 at {HOT[0]} Tmax"),
    ("figure elevation axis top (mRL)", Z_TOP, "CLAUDE.md 8.1"), ("figure elevation axis bottom (mRL)", Z_BOT, "CLAUDE.md 8.1"),
    ("temperature colour-scale maximum (degC)", T_LIM, "CLAUDE.md 8.1"),
    ("first field production", FIELD_START.isoformat(), "production records"),
    ("dat-file time reference date", REF_DATE.isoformat(), "matched against JSON"),
] + [(f"[reference model] {k}", v, "reference_bc") for k, v in REFBC.items()]
if near_x:
    params += [(f"{near_x[0]}/{near_x[1]} intersection easting (m)", near_x[2], "mapping"),
               (f"{near_x[0]}/{near_x[1]} intersection northing (m)", near_x[3], "mapping"),
               ("distance from intersection to upflow-well centre (m)", near_x[4], "mapping; wells")]
PAR = pd.DataFrame(params, columns=["parameter", "value", "source"])
write_csv(PAR, os.path.join(ADD, "parameters_used.csv"))

# ==========================================================================================
# REPORT
# ==========================================================================================
PCI = PC.set_index("well")
HYI = HYS.set_index(["well", "band"])
R = []
A = R.append


def wl(ws):
    return ", ".join(ws)


hi_w = [w for w in WELLS if conf[w] == "high"]
DROP1 = {}
for w_, g_ in gauge_rows.items():
    if g_["preprod"]:
        ds_ = sorted(g_["series"])
        DROP1[w_] = round(g_["series"][ds_[0]] - g_["series"][ds_[1]], 1)
GZ = {}
for w_, g_ in gauge_rows.items():
    s_ = W[w_]["survey"]
    GZ[w_] = round(interp_profile(s_.P, s_.z, g_["first_gauge_P_bar"]), 1)
PAR_EXTRA = [(f"{w_} gauge drop first to second reading (bar)", v_, "transient pressure") for w_, v_ in DROP1.items()]
PAR_EXTRA += [(f"{w_} survey elevation at first gauge pressure (mRL)", v_, "survey P interpolated") for w_, v_ in GZ.items()]
PAR = pd.concat([PAR, pd.DataFrame(PAR_EXTRA, columns=PAR.columns)], ignore_index=True)
write_csv(PAR, os.path.join(ADD, "parameters_used.csv"))
med_w = [w for w in WELLS if conf[w] == "medium"]
low_w = [w for w in WELLS if conf[w] == "low"]
gw9 = [w for w in GW if w in verdict and conf[w] == "medium" and gauge_rows[w]["preprod"]]
late_g = [w for w in GW if not gauge_rows[w]["preprod"]]
hot_out = [w for w in OUT_WELLS]
cased_max = HY[HY.band == "cased"].sort_values("ratio_measured_to_IF97").iloc[-1]
open_hot = HY[(HY.band == "open") & (~HY.cold_well)]
wh_below = GEO[GEO.wellhead_minus_surface_m < 0].sort_values("wellhead_minus_surface_m")
issues_int = INTEG[INTEG.status != "OK"]
UPW1 = UP_WELLS[0]

A("# Natural-state data interpretation for model calibration")
A("")
A(f"Generated by `Scripts/natural_state_analysis.py` on {dt.date.today().isoformat()}.")
A("")
A("Claim tags: [obs] observation; [interp] interpretation; [assumption] assumed value; [reference model] description "
  "of the starting model, never evidence. All elevations in mRL, pressures in bar (absolute), temperatures in °C.")
A("")
# ---------------- 1
A("## 1. Key conclusions for calibration")
A("")
k = 1
A(f"{k}. **The upflow lies beneath {wl(UP_WELLS)}.** {HOT[0]} reaches {WSI.loc[HOT[0], 'Tmax_degC']:.1f} °C at "
  f"{WSI.loc[HOT[0], 'Tmax_elev_mRL']:.0f} mRL with no inversion and the shallowest 200 °C crossing "
  f"({WSI.loc[HOT[0], 'z_200C_mRL']:.1f} mRL); {HOT[1]} reaches {WSI.loc[HOT[1], 'Tmax_degC']:.1f} °C at "
  f"{WSI.loc[HOT[1], 'Tmax_elev_mRL']:.0f} mRL [obs]. "
  + (f"The {near_x[0]} / {near_x[1]} intersection lies {near_x[4]:.0f} m from the centre of these wells [obs], so the "
     f"faults are a plausible conduit [interp]. " if near_x else "")
  + f"Liquid at {WSI.loc[HOT[0], 'Tmax_degC']:.1f} °C carries about {H_UP:.1f} kJ/kg (IF97 Region 1) [obs].")
A(f"   **Action:** put the deep mass and heat input below the columns of {wl(UP_WELLS)} (columns "
  f"{', '.join(sorted({split_block(block_at(w, WSI.loc[w, 'TD_mRL']))[0] for w in UP_WELLS}))}), with an enthalpy of "
  f"at least {H_UP:.1f} kJ/kg; calibrate its rate against the temperatures of the outflow wells.")
k += 1
tops = WS[WS.regime == "outflow"].inversion_top_mRL
A(f"{k}. **A lateral outflow spreads from the upflow at {tops.max():.0f} to {tops.min():.0f} mRL.** "
  + "; ".join(f"{w} has a {WSI.loc[w, 'inversion_degC']:.1f} °C inversion below its {WSI.loc[w, 'Tmax_degC']:.1f} °C "
              f"peak at {WSI.loc[w, 'inversion_top_mRL']:.0f} mRL" for w in OUT_WELLS) + " [obs]. "
  "Peak temperatures fall, and peak depths deepen, away from the upflow wells, which is the signature of a hot "
  "tongue moving laterally above cooler water [interp].")
A(f"   **Action:** give the model a laterally permeable layer between about {tops.max():.0f} and {tops.min():.0f} mRL "
  "connecting the upflow to these wells, with lower vertical permeability below it; weight the peak and the "
  "inversion temperatures fully (they are the most diagnostic targets).")
k += 1
A(f"{k}. **Cool water underlies the outflow and enters from the margins.** "
  + "; ".join(f"{w} falls to {WSI.loc[w, 'Tmax_degC'] - WSI.loc[w, 'inversion_degC']:.1f} °C at "
              f"{WSI.loc[w, 'inversion_base_mRL']:.0f} mRL" for w in deep_cool) + " [obs]. "
  + "; ".join(f"{w} reaches only {WSI.loc[w, 'Tmax_degC']:.1f} °C (at {WSI.loc[w, 'Tmax_elev_mRL']:.0f} mRL)"
              for w in REC_WELLS) + " [obs]. "
  "Cold recharge or downflow occurs at the margins and a cool deep inflow runs beneath the outflow [interp].")
A("   **Action:** allow cold water to enter at depth from the lateral boundaries (or from the north and southwest "
  "margins) rather than sealing them; the deep isothermal sections are targets for that inflow, not for conduction.")
k += 1
A(f"{k}. **A conductive cap overlies the hot wells.** Shallow gradients in the hot wells are "
  + ", ".join(f"{WSI.loc[w, 'shallow_gradient_degC_per_km']:.1f} ({w})" for w in CAPW.well)
  + f" °C/km, and the gradient drops below 50 °C/km at {CAPW.cap_base_mRL.max():.0f} to {CAPW.cap_base_mRL.min():.0f} mRL "
  "[obs]. The cap base deepens away from the upflow [interp].")
A(f"   **Action:** use low permeability from the surface to the cap base in these columns; calibrate it against the "
  "shallow gradients (open-hole weights) rather than against shallow pressures.")
k += 1
A(f"{k}. **Deep pressures are connected across the field.** Feed-zone pressures of all wells lie within "
  f"{FZ.P_minus_trend_bar.abs().max():.2f} bar of one hydrostatic trend fitted to the gauge-verified feed zones "
  f"(RMS {RMS_V:.3f} bar, implied density {-TREND_B * 1e5 / G:.1f} kg/m3) [obs]. Only "
  + wl([w for w in WELLS if fz_maxres.get(w, 0) > TOL]) + f" exceeds the {TOL:.2f} bar tolerance "
  f"({fz_maxres.get('WELL-1', np.nan):.2f} bar below the trend) [obs].")
A("   **Action:** use feed-zone pressures as the pressure targets (weight_P); a single connected reservoir is "
  "adequate; treat the WELL-1 offset as a test of the permeability between WELL-1 and the main field (see section 7).")
k += 1
A(f"{k}. **The surveys represent the natural state, despite their date label.** All surveys are dated "
  f"{W[WELLS[0]]['survey_date'].isoformat()}, after first field production on {FIELD_START.isoformat()} [obs]. "
  f"Yet {wl(hi_w)} match their first gauge readings (taken on {FIELD_START.isoformat()}) to within 0.1 bar, whereas "
  f"the same gauges fell by {min(DROP1.values()):.1f} to {max(DROP1.values()):.1f} bar by their second reading one year "
  f"later [obs]. "
  + (f"{late_g[0]}'s survey pressure at its gauge ({PCI.loc[late_g[0], 'survey_P_at_gauge_bar']:.2f} bar) is "
     f"{PCI.loc[late_g[0], 'survey_P_at_gauge_bar'] - PCI.loc[late_g[0], 'first_gauge_P_bar']:.1f} bar above its "
     f"first reading in {PCI.loc[late_g[0], 'first_gauge_date'][:4]} [obs]. " if late_g else "")
  + "The surveys therefore pre-date production and the date is a label problem [interp].")
A(f"   **Action:** use every survey as a natural-state target with the confidence weights below "
  f"(high: {wl(hi_w)}; medium: {wl(med_w)}; low: {wl(low_w)}); ask the data owner to confirm the survey dates.")
k += 1
if NBW:
    w = NBW[0]
    A(f"{k}. **The top of the upflow is close to boiling.** {w} has a minimum boiling margin of "
      f"{WSI.loc[w, 'min_boiling_margin_degC']:.2f} °C at {WSI.loc[w, 'min_boiling_margin_elev_mRL']:.0f} mRL [obs]; "
      "no survey point is two-phase (margin ≤ 1 °C) [obs].")
    A(f"   **Action:** the model must keep {w} liquid but within a few degrees of saturation; an upflow that is "
      "too hot or too shallow will boil here, so use this margin as a constraint on the upflow rate and depth.")
    k += 1
A(f"{k}. **Pressure gradients behind casing exceed a static water column.** Measured to IF97 ratios in cased "
  f"windows reach {cased_max.ratio_measured_to_IF97:.4f} ({cased_max.well}), against a cold-well component of "
  f"{cold_comp.get('cased', np.nan):.4f}; in open hole the cold-well component is {cold_comp.get('open', np.nan):.4f} "
  f"and hot-well ratios lie between {open_hot.ratio_measured_to_IF97.min():.4f} and "
  f"{open_hot.ratio_measured_to_IF97.max():.4f} [obs]. Datum, gravity or unit conventions and two-phase fluid cannot "
  "explain an excess that is confined to the cased, hot, shallow part of the profiles (section 4) [interp].")
A("   **Action:** do not use pressures above the feed zones as targets (their weight_P is zero); if the model's cap "
  "pressures matter, treat the excess as a sign of upward flow through the cap and resolve its origin with the "
  "data owner.")
A("")
# ---------------- 2
A("## 2. The data")
A("")
A("### 2.1 What exists")
A("")
A("The input inventory (`../Additional_Output/input_inventory.csv`) assigns every file under `Baseline/` a role [obs]:")
A("")
A("| Role | Path | Note |")
A("|---|---|---|")
for r_ in INV.itertuples():
    A(f"| {r_.role} | `{r_.path}` | {r_.note} |")
A("")
A(f"There are {len(WELLS)} wells in one group ('{list(GROUP_WELLS)[0]}'). Each has one temperature and one pressure "
  "survey on the same date and at the same elevations, a track, a casing/liner design and feed zones [obs]. "
  f"{wl(GW)} carry a transient pressure gauge; {wl([w for w in WELLS if W[w]['producer']])} have production records; "
  f"{wl([w for w in WELLS if W[w]['injector']])} have injection records [obs].")
A("")
A("### 2.2 How it was checked")
A("")
A(f"- **JSON against plain-column files:** {len(XCHK)} files compared; every value matches exactly "
  f"(largest difference {XCHK.max_abs_diff.max():.1f}) [obs]. Transient files use years since "
  f"{REF_DATE.isoformat()} (established by matching), production mass flow is stored with a negative sign, and "
  f"the WELL-10 files use the irregular prefix 'Wel10' [obs]. Details: `../Additional_Output/dat_json_crosscheck.csv`.")
A(f"- **Survey, feed zones, casing and track:** {len(INTEG)} checks in `../Additional_Output/integrity_checks.csv`; "
  f"{len(issues_int)} issue(s): " + "; ".join(f"{r_.well} {r_.check} ({r_.detail})" for r_ in issues_int.itertuples())
  + " [obs]. All survey points and feed zones lie within the tracks and every feed zone lies within a slotted liner [obs].")
A(f"- **Geometry:** survey and feed-zone blocks are taken from the well-to-block file; the wellhead column found by "
  "point-in-polygon on the grid agrees with the block file for every well; model surface values come from the grid "
  "surface section [obs]. "
  + f"{len(wh_below)} wellheads lie below the model surface, by up to {-wh_below.wellhead_minus_surface_m.min():.1f} m "
  f"({wh_below.well.iloc[0]}) [obs]; no survey point lies above the model surface [obs].")
A(f"- **Thermodynamics:** IAPWS-IF97 Region 4 and Region 1 pass the verification tables; Tsat(1.01325 bar) = "
  f"{tsat(P_ATM):.2f} °C and Tsat(100 bar) = {tsat(100.0):.2f} °C [obs].")
A("")
A("### 2.3 Provenance labels")
A("")
src_vals = sorted(set(PROV.Source.astype(str)) - {""})
A(f"Every `Source` and `Comment` field of the surveys, gauges and flow records is empty; the only source label is "
  f"'{', '.join(src_vals)}' on the well designs; every survey carries `StableFlag` = True [obs]. "
  f"Feed-zone quality tags found: {', '.join(TAGS)}; mapped to weight_P factors "
  + ", ".join(f"{t} = {TAG_MAP[t]}" for t in TAGS) + " ('Synthetic' is not a quality grade, so it takes the lowest "
  "factor) [obs]. The data are therefore labelled synthetic; nothing in the files records who measured them or how "
  "[obs].")
A("")
A("### 2.4 Natural-state timing verdict for each well")
A("")
A(f"First production anywhere in the field: {FIELD_START.isoformat()} [obs]. Injectate enthalpy {H_INJ:.1f} kJ/kg "
  f"gives {T_INJ:.1f} °C (IF97 Region 1 inverse at 1.01325 bar) [obs]. Reading precision of the gauges: 0.1 bar. "
  f"Trend tolerance: {TOL:.2f} bar (3 bar exceeds twice the {RMS_V:.3f} bar RMS of the gauge-verified feed zones) [obs].")
A("")
A("| Well | Gauge (mRL) | Survey P at gauge (bar) | First reading (date, bar) | Difference (bar) | Max feed-zone offset from trend (bar) | Verdict | Confidence |")
A("|---|---|---|---|---|---|---|---|")
for r_ in PC.itertuples():
    A(f"| {r_.well} | {fmt(r_.gauge_elev_mRL, 0)} | {fmt(r_.survey_P_at_gauge_bar, 2)} | "
      f"{(r_.first_gauge_date + ', ' + fmt(r_.first_gauge_P_bar, 1)) if r_.first_gauge_date else 'no gauge'} | "
      f"{fmt(r_.diff_survey_minus_first_bar, 2)} | {fmt(r_.fz_max_abs_P_minus_trend_bar, 2)} | {r_.verdict} | {r_.confidence} |")
A("")
A("Reason for each verdict [obs unless stated]:")
A("")
for w in WELLS:
    A(f"- **{w}:** {reason[w]}.")
A("")
if gw9:
    w = gw9[0]
    A(f"{w}'s survey pressure at the gauge is {PCI.loc[w, 'diff_survey_minus_first_bar']:+.2f} bar above the first "
      f"reading. A disturbed survey would read lower, not higher; the survey value {PCI.loc[w, 'first_gauge_P_bar']:.1f} "
      "bar is reached at "
      f"{GZ[w]:.1f} mRL, {GZ[w] - PCI.loc[w, 'gauge_elev_mRL']:.1f} m above the stated gauge elevation, so a gauge-depth or datum error in a "
      "deviated well is the likelier cause [interp]. This is recorded as a data question (section 7).")
    A("")
A("Confidence differs between wells, so the factor does change relative weights [obs].")
A("")
A("![F7 gauge check](figures/gauge_check.png)")
A("")
A(f"*F7. Gauge pressure against time at each gauge elevation ({', '.join(f'{w} {gauge_rows[w]['gauge_elev_mRL']:.0f} mRL' for w in GW)}), "
  "with the survey pressure at the gauge, the ±0.1 bar band around the first reading, the first own flow and first "
  "field production; grey = flowing. Claim type: obs.*")
A("")
# ---------------- 3
A("## 3. Temperature structure")
A("")
A("Per-well characteristics from `well_summary.csv` (definitions of CLAUDE.md §4.4; blank = not defined) [obs]:")
A("")
A("| Well | Tmax °C (mRL) | T at TD °C (mRL) | Inversion °C (top, base mRL) | Shallow / deep gradient °C/km | z100 / z200 mRL | Min boiling margin °C (mRL) | Water level mRL | Regime | NS confidence |")
A("|---|---|---|---|---|---|---|---|---|---|")
for w in WELLS:
    r_ = WSI.loc[w]
    A(f"| {w} | {r_.Tmax_degC:.1f} ({r_.Tmax_elev_mRL:.0f}) | {r_.T_at_TD_degC:.1f} ({r_.TD_mRL:.1f}) | "
      f"{fmt(r_.inversion_degC)} ({fmt(r_.inversion_top_mRL, 0)}, {fmt(r_.inversion_base_mRL, 0)}) | "
      f"{fmt(r_.shallow_gradient_degC_per_km)} / {fmt(r_.deep_gradient_degC_per_km)} | "
      f"{fmt(r_.z_100C_mRL)} / {fmt(r_.z_200C_mRL)} | {fmt(r_.min_boiling_margin_degC, 2)} ({fmt(r_.min_boiling_margin_elev_mRL, 0)}) | "
      f"{fmt(r_.water_level_mRL)} | {r_.regime} | {r_.ns_confidence} |")
A("")
A("Regime notes [interp]: " + " ".join(f"{w}: {REGIME_NOTE[w]}." for w in WELLS if REGIME_NOTE[w]))
A("")
A("Wellbore effects considered before interpreting features [interp]:")
A("")
iso_w = sorted(set(TG[TG.isothermal_section].well), key=WELLS.index)
A(f"- Open-hole isothermal sections (≥ 300 m within 5 °C) occur in {wl(iso_w)} [obs]. They lie between feed zones "
  "or below the deepest one, where internal flow in the well could homogenise temperature, so these points carry "
  "weight 0.5 rather than 0.8 [interp]. The deep cool sections are nevertheless continuous with the inversions above "
  "them and recur in several wells at similar depths, so they are not purely a wellbore artefact [interp].")
inv_cased = [w for w in OUT_WELLS if WSI.loc[w, "inversion_top_mRL"] >= W[w]["shoe"]]
A(f"- Inversion peaks lie below the production casing shoe in {wl([w for w in OUT_WELLS if w not in inv_cased])}"
  + (f" and at or above it in {wl(inv_cased)}" if inv_cased else "") + " [obs]. A peak near a feed zone could be a "
  "wellbore inflow; in the outflow wells the peak is broad (hundreds of metres) and repeats from well to well at "
  "similar depths, which a single-well effect would not produce [interp].")
A(f"- The steep shallow gradients are measured behind casing in the hot wells, where temperatures are closer to the "
  "conductive formation temperature [interp].")
A("")
A("![F2 group profiles](figures/profile_panels_All.png)")
A("")
A(f"*F2. Temperature profiles of all wells on identical axes ({Z_TOP} to {Z_BOT} mRL, 0 to {T_LIM:.0f} °C), with "
  "Tsat(P), casing shoe, feed zones, wellhead and model surface. Individual well figures with pressure panels are "
  "`figures/profile_<well>.png`. Claim type: obs.*")
A("")
for w in WELLS:
    A(f"![F2 {w}](figures/profile_{w}.png)")
    A("")
    A(f"*F2. {w}: T and P against elevation ({Z_TOP} to {Z_BOT} mRL) with locator map, survey {W[w]['survey_date'].isoformat()}, "
      f"verdict {verdict[w]} ({conf[w]}). Claim type: obs.*")
    A("")
med_rows = BA.groupby("band_top_mRL").T_field_median_degC.first()
A("Lateral pattern in 500 m bands (`../Additional_Output/T_P_at_elevation_bands.csv`) [obs]: field medians "
  + ", ".join(f"{med_rows[b]:.1f} °C ({b} to {b - 500} mRL)" for b in BANDS) + ". "
  + f"The upflow and outflow wells are the hot outliers above -1000 mRL, and {wl(REC_WELLS)} the cold ones at every depth [obs].")
A("")
A("![F3 temperature](figures/maps_temperature_bands.png)")
A("")
A(f"*F3. Well mean temperature minus the band field median, 500 m bands from {BANDS[0]} to {BANDS[-1] - 500} mRL, plotted "
  "at each well's track position at mid-band; crosses within ±10 °C, open circles = no data. Claim type: obs.*")
A("")
for L in LINES:
    A(f"![F4 section {L['name']}](figures/section_{L['name']}.png)")
    A("")
    A(f"*F4. Section {L['name']}-{L['name']}' from ({L['A'][0]:.1f}, {L['A'][1]:.1f}) to ({L['B'][0]:.1f}, {L['B'][1]:.1f}), "
      f"{L['why']}; wells within {L['maxdist']:.0f} m: "
      + ", ".join(f"{w} ({o:+.0f} m)" for w, o in L["wells"])
      + f". Temperature points and isotherms within the data envelope; vertical exaggeration {VE[(L['name'], 'T')]:.1f}. "
        "Claim type: obs (points), interp (isotherms).*")
    A("")
A("![F5 P-T state](figures/PT_state.png)")
A("")
A("*F5. All survey points and feed zones on P-T axes with the IF97 saturation curve, ending at the critical point. "
  "Claim type: obs.*")
A("")
# ---------------- 4
A("## 4. Pressure structure")
A("")
A(f"Survey pressures at common elevations (interpolated, never extrapolated; "
  "`../Additional_Output/T_P_at_common_elevations.csv`), and the deviation from the verified trend "
  f"P = {TREND_A:.3f} + {-TREND_B:.5f} × depth below 0 mRL [obs]:")
A("")
zs_ = [z for z in COMMON_Z if z in set(CE.elevation_mRL)]
A("| Well | " + " | ".join(f"{z} mRL" for z in zs_) + " |")
A("|---|" + "---|" * len(zs_))
for w in WELLS:
    cells = []
    for z in zs_:
        q = CE[(CE.well == w) & (CE.elevation_mRL == z)]
        cells.append(f"{q.P_bar.iloc[0]:.1f} ({q.P_minus_trend_bar.iloc[0]:+.1f}{'' if q.open_hole.iloc[0] else ', cased'})"
                     if len(q) else "")
    A(f"| {w} | " + " | ".join(cells) + " |")
A("")
A("Values in brackets are P minus the verified trend; 'cased' marks wellbore-column values that are not formation "
  "pressure (§4.2) [obs]. Below the feed zones the deviation grows with depth in cold wells (for example WELL-8) "
  "because their colder water column is denser than the hot trend fluid; this is a density effect, not overpressure "
  "[interp].")
A("")
A("Feed-zone pressures (the formation pressure targets) [obs]:")
A("")
A("| Well | Feed zone | Mid (mRL) | Quality | P (bar) | P minus trend (bar) | T (°C) | Block |")
A("|---|---|---|---|---|---|---|---|")
for r_ in FZOUT.itertuples():
    A(f"| {r_.well} | {r_.feedzone} | {r_.mid_mRL:.1f} | {r_.quality} | {r_.P_bar:.1f} | {r_.P_minus_trend_bar:+.2f} | "
      f"{r_.T_degC:.1f} | {r_.block} |")
A("")
A(f"Pivot pressures (dominant feed zone, CLAUDE.md §4.4): " + "; ".join(
    f"{w} {WSI.loc[w, 'pivot_P_bar']:.1f} bar at {WSI.loc[w, 'pivot_elev_mRL']:.1f} mRL" for w in WELLS) + " [obs]. "
  f"Water levels from the three shallowest liquid points range from {WS.water_level_mRL.min():.1f} to "
  f"{WS.water_level_mRL.max():.1f} mRL; WELL-1's ({WSI.loc['WELL-1', 'water_level_mRL']:.1f} mRL) is "
  f"{WSI.loc['WELL-1', 'water_level_extrapolation_m']:.1f} m above its top survey point, i.e. the well is artesian "
  "[obs].")
A("")
A("**Hydrostatic check by band** (`../Additional_Output/hydrostatic_check.csv`, `hydrostatic_summary.csv`; windows "
  "≥ 300 m in single-phase liquid intervals; ratio = measured ΔP / ΔP from IF97 Region 1 density over the measured "
  "T profile, g = 9.80665 m/s2; rounding uncertainty of the ratio ≤ "
  f"{HY.ratio_rounding_uncertainty.max():.4f}) [obs]:")
A("")
A("| Well | Cased: n, median ratio, excess over cold component | Open: n, median ratio, excess over cold component |")
A("|---|---|---|")
for w in WELLS + ["ALL"]:
    cells = []
    for band in ("cased", "open"):
        if (w, band) in HYI.index:
            h = HYI.loc[(w, band)]
            cells.append(f"{int(h.n_windows)}, {h.median_ratio:.4f}, {h.median_excess:+.4f}")
        else:
            cells.append("no window")
    A(f"| {w} | {cells[0]} | {cells[1]} |")
A("")
A(f"The cold wells ({wl([w for w in WELLS if WSI.loc[w, 'Tmax_degC'] < 100])}) give a field-wide component of "
  f"{cold_comp.get('cased', np.nan):.4f} (cased) and {cold_comp.get('open', np.nan):.4f} (open) [obs]. "
  "This is presented as a possible datum or convention issue, not physical overpressure: a few percent excess in "
  "cold columns could come from a pressure-gauge or depth convention [interp]. Checks on the larger excess in the "
  "hot cased windows [interp]: a datum offset cannot change a gradient; a gravity or unit convention would act at "
  f"every depth, but the open-hole windows of the same wells are close to the cold component; two-phase fluid would "
  f"make the column lighter, not heavier; rounding contributes at most {HY.ratio_rounding_uncertainty.max():.4f}. "
  "None of these explains it. The excess is real in the data and confined to the shallow hot section; it is "
  "consistent with these 'pressures' being formation pressures with upward flow through a low-permeability cap "
  "rather than a static wellbore column, which is a question for the data owner (section 7).")
A("")
A("![F3 pressure](figures/maps_pressure_bands.png)")
A("")
A("*F3 (pressure). Well mean pressure minus the verified hydrostatic trend in 500 m bands, classes of 2 bar. "
  "Shallow bands include cased wellbore pressures. Claim type: mixed (trend is fitted).*")
A("")
for L in LINES:
    A(f"![F4 section {L['name']} pressure](figures/section_{L['name']}_P.png)")
    A("")
    A(f"*F4 (pressure). Section {L['name']}-{L['name']}': survey pressure minus the verified hydrostatic trend, "
      "contoured at 2 bar within the data envelope. Claim type: mixed.*")
    A("")
# ---------------- 5
A("## 5. Conceptual model")
A("")
A(f"![F6 conceptual model](figures/conceptual_model.png)")
A("")
A(f"*F6. Conceptual model from the data alone: plan view and section {LB['name']}-{LB['name']}' (the F4 line from "
  f"the hottest well {HOT[0]} towards the conductive margin, chosen because it crosses the upflow, two outflow wells "
  "and the margin; line A crosses only the upflow wells and WELL-8). Numbers refer to the table below. "
  "Claim type: interp.*")
A("")
A("| No. | Element | Symbol | Location | Supporting observation | Claim type |")
A("|---|---|---|---|---|---|")
for r_ in ELDF.itertuples():
    A(f"| {r_.no} | {r_.element} | {r_.symbol} ({r_.line_style}) | {r_.location} | {r_.observation} | {r_.claim_type} |")
A("")
A("Relation to mapped features [obs unless stated]:")
A("")
for w in WELLS:
    g_ = GEOI.loc[w]
    A(f"- {w}: nearest fault {g_.nearest_fault} at {g_.dist_nearest_fault_m:.1f} m; nearest water body "
      f"{g_.nearest_water_body} at {g_.dist_nearest_water_body_m:.1f} m.")
A("")
A("What the data cannot decide [interp]:")
A("")
A(f"- Whether the upflow rises along {near_x[0] if near_x else 'a fault'} or {near_x[1] if near_x else ''} or is "
  f"distributed: no well penetrates below {WSI.loc[UPW1, 'TD_mRL']:.0f} mRL at {UPW1}. A deep well, or a deep "
  "temperature log, near the fault intersection would decide it.")
A(f"- Whether the cold water beneath the outflow comes from the south, east or north: deep temperatures in "
  f"{wl(REC_WELLS)} are colder than any inflow water seen elsewhere, but no well lies between them and the outflow wells "
  "at depth. Deep wells or deep temperature logs on the margins would decide it.")
A("- Whether WELL-1's lower feed-zone pressures reflect a flow barrier (it lies east of Fault 2 and near the river) or "
  "a pressure-gauge problem: an interference test or a repeat survey would decide it.")
A(f"- Whether {wl(REC_WELLS)} show downflow from the surface (low gradients near the top) or lateral inflow (inversion "
  "at depth in WELL-8): both are consistent with the profiles.")
A("")
# ---------------- 6
A("## 6. Recommended calibration set-up")
A("")
A("**Targets and weights** [obs, rules of CLAUDE.md §5.5]:")
A("")
nt = len(TG)
A(f"- Temperature targets: all {nt} survey points in `NS_calibration_targets.csv` (`weight`), feed-zone temperatures "
  f"in `feedzones.csv` (`weight_T`), and layer-centre values in `T_P_at_layer_centres.csv` for block-by-block comparison.")
A(f"- Pressure targets: only the {int((TG.weight_P > 0).sum())} survey points inside feed zones (`weight_P`) and the "
  f"{len(FZOUT)} feed-zone pressures; all other pressures have weight_P = 0.")
A("- Location weights: above model surface 0; within 30 m of the wellhead or in the column's top block 0.3; feed zone "
  "1.0; behind casing 0.5; open-hole isothermal section 0.5; other open hole 0.8. No point lies above the model "
  "surface or in a top block, so the 0 and top-block cases are unused [obs].")
counts = TG.groupby(TG.weight_reason.str.split(" \\(").str[0]).size()
A("- Points by location class: " + ", ".join(f"{k_} {v}" for k_, v in counts.items()) + " [obs].")
A(f"- Confidence factors: high 1.0 ({wl(hi_w)}), medium 0.75 ({wl(med_w)}), low 0.5 ({wl(low_w)}). Pressure quality "
  "factors: HIGH 1.0, MEDIUM 0.8, LOW and 'Synthetic' 0.5. No departure from the defaults [obs].")
A("")
A("**Heat and mass in and out** [interp]:")
A("")
A(f"- In: a deep hot upflow below {wl(UP_WELLS)} (≥ {H_UP:.1f} kJ/kg), plus cold water at depth from the lateral "
  f"margins, and possibly downflow at {wl(REC_WELLS)}.")
A(f"- Out: lateral outflow at {tops.max():.0f} to {tops.min():.0f} mRL towards {wl(OUT_WELLS)}; ultimately to the "
  "surface or the rivers (the data do not show the discharge point: the files hold no surface thermal features or spring flows).")
A("")
A("**Permeability structure implied by the data** [interp]:")
A("")
A(f"- Low permeability cap from the surface to {CAPW.cap_base_mRL.max():.0f} to {CAPW.cap_base_mRL.min():.0f} mRL over "
  "the hot area.")
A(f"- A laterally permeable outflow layer at about {tops.max():.0f} to {tops.min():.0f} mRL, with lower vertical "
  "permeability beneath it, so that cooler water can persist below.")
A(f"- A vertical conduit beneath {wl(UP_WELLS)}"
  + (f", possibly along {near_x[0]} and {near_x[1]}" if near_x else "") + ".")
A("- Good lateral pressure connection at feed-zone depths across the field (all feed zones within the tolerance of "
  "one trend except WELL-1).")
A("")
A("**Boundary requirements** [interp]:")
A("")
A("- Lateral boundaries that admit cold water at depth (the deep cool inflow and the cold margin wells require it).")
A("- A top boundary that allows the cap to remain conductive and liquid; WELL-1's artesian water level "
  f"({WSI.loc['WELL-1', 'water_level_mRL']:.1f} mRL) suggests confined pressure in the east.")
A("")
A("**Starting model for comparison** [reference model]: the reference boundary conditions use a default basal heat "
  f"flux of {REFBC.get('heat_flux_default_Wm2', np.nan)} W/m2, a default rainfall of "
  f"{REFBC.get('rain_default_mm_per_yr', np.nan):.1f} mm/yr with infiltration fraction "
  f"{REFBC.get('infiltration_fraction_default', np.nan)}, an atmosphere at {REFBC.get('atmosphere_T_default_degC', np.nan)} °C "
  f"and {REFBC.get('n_lake_columns', 0)} lake columns. These are properties of the starting model, not evidence.")
A("")
# ---------------- 7
A("## 7. Data quality issues and open questions")
A("")
A("Ranked by their effect on calibration:")
A("")
q = 1
A(f"{q}. **Survey date label** [obs]: every survey is dated {W[WELLS[0]]['survey_date'].isoformat()}, after production "
  f"and injection began ({FIELD_START.isoformat()}), but the gauge test shows pre-production pressures. "
  "Effect: decides whether the surveys can be natural-state targets at all. Resolve: confirm the survey dates.")
q += 1
A(f"{q}. **Super-hydrostatic shallow pressure gradients** [obs]: cased-window ratios up to "
  f"{cased_max.ratio_measured_to_IF97:.4f} ({cased_max.well}). Effect: pressures above the reservoir cannot be used, "
  "and if they are formation pressures they constrain the cap permeability. Resolve: ask how the pressures were "
  "obtained (wellbore log or formation estimate).")
q += 1
if gw9:
    w = gw9[0]
    A(f"{q}. **{w} gauge mismatch** [obs]: survey {PCI.loc[w, 'survey_P_at_gauge_bar']:.2f} bar vs first reading "
      f"{PCI.loc[w, 'first_gauge_P_bar']:.1f} bar at {PCI.loc[w, 'gauge_elev_mRL']:.0f} mRL in a deviated well "
      f"(horizontal deviation {GEOI.loc[w, 'max_horizontal_deviation_m']:.1f} m). Effect: {w} is medium rather than "
      "high confidence. Resolve: confirm whether the gauge elevation is true vertical or along-hole.")
    q += 1
A(f"{q}. **WELL-1 off the pressure trend** [obs]: feed zones {FZ[FZ.well == 'WELL-1'].P_minus_trend_bar.min():.2f} and "
  f"{FZ[FZ.well == 'WELL-1'].P_minus_trend_bar.max():.2f} bar from the trend (tolerance {TOL:.2f} bar), so WELL-1 is "
  "'possibly disturbed' (low) by the rubric, although it has no flow records [obs]. A natural lateral pressure "
  "difference across Fault 2 is an alternative explanation [interp]. Effect: halves WELL-1 weights. Resolve: check "
  "whether WELL-1 was ever flowed or whether it is separated from the main field.")
q += 1
A(f"{q}. **WELL-5 is an injector hotter than the injectate** [obs]: injection cooling could be hidden, so its weights "
  "are halved. Resolve: confirm the survey pre-dates injection.")
q += 1
A(f"{q}. **Wellheads below the model surface** [obs]: by up to {-wh_below.wellhead_minus_surface_m.min():.1f} m "
  f"({wh_below.well.iloc[0]}); the model surface is smoother than the topography at the wells. Effect: the shallowest "
  "points sit deeper in the model's column than in reality. Resolve: accept, or refine the surface.")
q += 1
for r_ in issues_int.itertuples():
    A(f"{q}. **{r_.well} {r_.check}** [obs]: {r_.detail}. Effect: none on the natural-state targets (the survey stops "
      "at the track TD). Resolve: correct the well design or the track.")
    q += 1
A(f"{q}. **Quality tag 'Synthetic'** on WELL-10's feed zone is not a quality grade [obs]; it takes factor 0.5. Resolve: "
  "supply a grade.")
q += 1
A(f"{q}. **Empty provenance fields** [obs]: no Source or Comment on surveys and gauges; well designs say "
  f"'{', '.join(src_vals)}'. Effect: the trust in every value rests on the labels alone.")
q += 1
A(f"{q}. **File naming** [obs]: the WELL-10 plain-column files use 'Wel10'; production flows in those files have the "
  "opposite sign to the JSON. No effect on the analysis (JSON is primary).")
A("")
# ---------------- 8
A("## 8. Output file index")
A("")
A("`Output/Calibration_Input/` (what the calibrator needs):")
A("")
for f_ in ["NS_data_interpretation.md", "NS_calibration_targets.csv", "feedzones.csv", "T_P_at_layer_centres.csv",
           "well_summary.csv", "figures_index.csv"]:
    A(f"- `{f_}`")
A(f"- `figures/` ({len(FIDX)} standard figures F1 to F7; see `figures_index.csv`)")
A("")
A("`Output/Additional_Output/` (supporting detail):")
A("")
ADD_FILES = ["input_inventory.csv", "pressure_consistency_check.csv", "hydrostatic_check.csv", "hydrostatic_summary.csv",
             "well_geometry_check.csv", "T_P_at_elevation_bands.csv", "T_P_at_common_elevations.csv",
             "dat_json_crosscheck.csv", "integrity_checks.csv", "provenance_fields.csv", "section_lines.csv",
             "pressure_connections.csv", "conceptual_model_elements.csv", "parameters_used.csv",
             "run_console_log.txt", "acceptance_check.txt", "process_log.md"]
for f_ in ADD_FILES:
    A(f"- `{f_}`")
A("")
REPORT = "\n".join(R) + "\n"
write_text(REPORT, os.path.join(CAL, "NS_data_interpretation.md"))

# console log last
log("run complete")
with open(os.path.join(ADD, "run_console_log.txt"), "w", encoding="utf-8", newline="\n") as f:
    f.write(_CONSOLE.getvalue())
