"""Natural-state data interpretation for geothermal model calibration (CLAUDE.md v4).

Run from the repository root:   python Scripts/natural_state_analysis.py
It reads only Baseline/ and regenerates every output under Output/ (except process_log.md,
which is written by hand), including the report, the input inventory and every figure.
"""
import ast
import bisect
import datetime as dt
import glob
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patches
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.path import Path as MPath
from matplotlib.tri import Triangulation, LinearTriInterpolator
from scipy.spatial import Delaunay

ROOT = os.getcwd()
BASE = os.path.join(ROOT, "Baseline")
CAL = os.path.join(ROOT, "Output", "Calibration_Input")
ADD = os.path.join(ROOT, "Output", "Additional_Output")
CALFIG = os.path.join(CAL, "figures")
ADDFIG = os.path.join(ADD, "figures")
for _d in (CAL, ADD, CALFIG, ADDFIG):
    os.makedirs(_d, exist_ok=True)

G = 9.80665
P_ATM = 1.01325
P_CRIT = 220.64
DPI = 150
PREC_P = 0.1          # reading precision of survey and gauge pressures [bar] (values quoted to 0.1 bar)
SECTION_MAX_OFFSET = 600.0   # m: wells projected onto a section line must lie within this distance
CONF_FACTOR = {"high": 1.0, "medium": 0.75, "low": 0.5}
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.append(s)


def write_csv(df, path, float_format=None):
    df.to_csv(path, index=False, lineterminator="\n", encoding="utf-8", float_format=float_format)


def write_text(path, text):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)


def rel(p):
    return os.path.relpath(p, ROOT).replace("\\", "/")


# =========================================================================================
# 1. IAPWS-IF97 (Region 4 saturation, Region 1 liquid) with validation
# =========================================================================================
_N4 = [None, 0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5,
       -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6,
       -0.23855557567849, 0.65017534844798e3]


def tsat_c_scalar(p_bar):
    """Saturation temperature [degC] at p [bar]; NaN above the critical pressure or below the triple point."""
    p = p_bar / 10.0
    if not np.isfinite(p) or p < 611.213e-6 or p > 22.064:
        return float("nan")
    n = _N4
    b = p ** 0.25
    e = b * b + n[3] * b + n[6]
    f = n[1] * b * b + n[4] * b + n[7]
    g = n[2] * b * b + n[5] * b + n[8]
    d = 2 * g / (-f - math.sqrt(f * f - 4 * e * g))
    t = (n[10] + d - math.sqrt((n[10] + d) ** 2 - 4 * (n[9] + n[10] * d))) / 2
    return t - 273.15


def psat_mpa(t_c):
    """Saturation pressure [MPa] at T [degC] (IF97 eq. 30)."""
    n = _N4
    th = (t_c + 273.15) + n[9] / ((t_c + 273.15) - n[10])
    a = th * th + n[1] * th + n[2]
    b = n[3] * th * th + n[4] * th + n[5]
    c = n[6] * th * th + n[7] * th + n[8]
    return (2 * c / (-b + math.sqrt(b * b - 4 * a * c))) ** 4


tsat_c = np.vectorize(tsat_c_scalar, otypes=[float])

_R1 = [(0, -2, 0.14632971213167), (0, -1, -0.84548187169114), (0, 0, -0.37563603672040e1),
       (0, 1, 0.33855169168385e1), (0, 2, -0.95791963387872), (0, 3, 0.15772038513228),
       (0, 4, -0.16616417199501e-1), (0, 5, 0.81214629983568e-3), (1, -9, 0.28319080123804e-3),
       (1, -7, -0.60706301565874e-3), (1, -1, -0.18990068218419e-1), (1, 0, -0.32529748770505e-1),
       (1, 1, -0.21841717175414e-1), (1, 3, -0.52838357969930e-4), (2, -3, -0.47184321073267e-3),
       (2, 0, -0.30001780793026e-3), (2, 1, 0.47661393906987e-4), (2, 3, -0.44141845330846e-5),
       (2, 17, -0.72694996297594e-15), (3, -4, -0.31679644845054e-4), (3, 0, -0.28270797985312e-5),
       (3, 6, -0.85205128120103e-9), (4, -5, -0.22425281908000e-5), (4, -2, -0.65171222895601e-6),
       (4, 10, -0.14341729937924e-12), (5, -8, -0.40516996860117e-6), (8, -11, -0.12734301741641e-8),
       (8, -6, -0.17424871230634e-9), (21, -29, -0.68762131295531e-18),
       (23, -31, 0.14478307828521e-19), (29, -38, 0.26335781662795e-22),
       (30, -39, -0.11947622640071e-22), (31, -40, 0.18228094581404e-23),
       (32, -41, -0.93537087292458e-25)]
_R = 0.461526  # kJ/(kg K)


def region1(t_c, p_bar):
    """IF97 Region 1: returns (specific volume m3/kg, specific enthalpy kJ/kg)."""
    T = t_c + 273.15
    pi = (p_bar / 10.0) / 16.53
    tau = 1386.0 / T
    gp = gt = 0.0
    for i, j, n in _R1:
        gp += -n * i * (7.1 - pi) ** (i - 1) * (tau - 1.222) ** j
        gt += n * (7.1 - pi) ** i * j * (tau - 1.222) ** (j - 1)
    v = _R * T / (p_bar / 10.0) * pi * gp * 1e-3
    h = _R * T * tau * gt
    return v, h


def rho_liquid(t_c, p_bar):
    """Liquid density [kg/m3] from IF97 Region 1 (valid to 350 degC / 1000 bar; clipped to 350 degC)."""
    t = np.minimum(np.asarray(t_c, float), 350.0)
    p = np.maximum(np.asarray(p_bar, float), 0.5)
    out = np.empty(np.broadcast(t, p).shape)
    it = np.nditer([t, p], flags=["multi_index"])
    for a, b in it:
        out[it.multi_index] = 1.0 / region1(float(a), float(b))[0]
    return out


def validate_if97():
    """Make the script fail if any IF97 implementation error exceeds 0.1 K (or the equivalent)."""
    errs = []
    # Region 4 verification table (IF97 Table 35): T [K] -> p [MPa]
    for tk, pm in [(300.0, 0.353658941e-2), (500.0, 0.263889776e1), (600.0, 0.123443146e2)]:
        errs.append(abs(tsat_c_scalar(pm * 10.0) + 273.15 - tk))
    # Region 4 verification table (IF97 Table 36): p [MPa] -> T [K]
    for pm, tk in [(0.1, 0.372755919e3), (1.0, 0.453035632e3), (10.0, 0.584149488e3)]:
        errs.append(abs(tsat_c_scalar(pm * 10.0) + 273.15 - tk))
    errs.append(abs(tsat_c_scalar(1.01325) - 99.97))
    errs.append(abs(tsat_c_scalar(100.0) - 311.0))
    # psat(T) round trip
    for tc in (50.0, 150.0, 250.0, 300.0):
        errs.append(abs(tsat_c_scalar(psat_mpa(tc) * 10.0) - tc))
    # Region 1 verification (IF97 Table 5): (T K, p MPa) -> v m3/kg, h kJ/kg
    r1 = [(300.0, 3.0, 0.100215168e-2, 0.115331273e3), (300.0, 80.0, 0.971180894e-3, 0.184142828e3),
          (500.0, 3.0, 0.120241800e-2, 0.975542239e3)]
    r1err = []
    for tk, pm, v, h in r1:
        vv, hh = region1(tk - 273.15, pm * 10.0)
        r1err.append(max(abs(vv - v) / v, abs(hh - h) / h))
    log("IF97 validation: max Region 4 error = %.2e K; max Region 1 relative error = %.2e" % (max(errs), max(r1err)))
    if max(errs) > 0.1 or max(r1err) > 1e-6:
        raise SystemExit("IF97 validation FAILED (error above 0.1 K or Region 1 table mismatch)")


# =========================================================================================
# 2. Input inventory (CLAUDE.md 2.1)
# =========================================================================================
def _sniff_json(path):
    """Role of a .json file from its content; None if not recognised."""
    try:
        if os.path.getsize(path) > 30e6:
            return None
        d = json.load(open(path, encoding="utf-8"))
    except Exception:
        return None
    if isinstance(d, dict):
        if "Well_Track" in d and "Downhole_Temperature" in d:
            return "well_data_json"
        if "Areas" in d:
            return "well_groups"
        vals = [v for v in d.values() if isinstance(v, dict)]
        if vals and all("Blocks" in v and "Interp_Weights" in v for v in vals):
            return "well_blocks_json"
        if vals and all("Data" in v and "type" in v for v in vals):
            return "mapping"
    return None


def classify(path):
    """(role, note) for one file under Baseline/ by location and content."""
    rp = rel(path)
    low = rp.lower()
    name = os.path.basename(path).lower()
    ext = os.path.splitext(name)[1]
    if "/parameters/" in low or "/production/" in low or "/futurescenarios/" in low or ext == ".xlsx" \
            or name == ".ds_store":
        return "out_of_scope", "production / injection / forecast set-up or history; not needed for the natural state"
    if "/naturalstate/" in low and "currentnssetup" not in name:
        return "reference_results", "reference natural-state model input/initial conditions; never evidence"
    if "/naturalstate/" in low:
        return "reference_other", "reference natural-state set-up; never evidence"
    if "/data/plotting/" in low:
        return "reference_other", "plotting configuration; hard-coded paths from another machine; not read"
    if ext == ".json":
        role = _sniff_json(path)
        if role == "well_data_json":
            return role, "per-well data (track, design, feed zones, surveys, transients, flows)"
        if role == "well_groups":
            return role, "well list and grouping"
        if role == "well_blocks_json":
            return role, "well-to-block intersections; authority for block assignment"
        if role == "mapping":
            return role, "faults, rivers and lake polylines"
        if "rocktype" in name:
            return "reference_rock", "reference rock types and permeabilities; never evidence"
        if "boundary" in name:
            return "reference_bc", "reference boundary conditions; never evidence"
        if "area_data" in name:
            return "out_of_scope", "field production totals; not needed (production used only via well files)"
        return "reference_other", "model set-up / configuration; not evidence"
    if ext == ".dat":
        with open(path, "r", errors="ignore") as f:
            head = [next(f, "").strip() for _ in range(5)]
        if "VERTICES" in head:
            return "geometry_dat", "MULgraph grid file: vertices, grid, layers, surface and well tracks"
        if "/data/data_dat/" in low:
            return "well_data_dat", "plain-column copy of well data; cross-check only"
        return "reference_other", "unrecognised .dat; not read"
    if ext == ".auxjson":
        return "reference_other", "reference model auxiliary geometry; not read"
    if ext == ".incon":
        return "reference_results", "reference initial conditions; never evidence"
    return "out_of_scope", "unclassified file; not read"


def build_inventory():
    """Walk Baseline/, assign roles, collapse directories whose files share one role, write the CSV."""
    files = []
    for dp, _, fs in os.walk(BASE):
        for f in sorted(fs):
            files.append(os.path.join(dp, f))
    files.sort()
    roles = {f: classify(f) for f in files}
    rows, used = [], set()

    def collapse(d):
        """Topmost directory below which every file has one role; returns role or None."""
        sub = [f for f in files if f.startswith(d + os.sep)]
        r = {roles[f][0] for f in sub}
        return (r.pop(), sub) if sub and len(r) == 1 else (None, sub)

    def walk(d):
        r, sub = collapse(d)
        if r is not None and d != BASE and r not in ("geometry_dat", "well_blocks_json", "well_groups"):
            rows.append((r, rel(d), roles[sub[0]][1] + f" ({len(sub)} files)"))
            used.update(sub)
            return
        for e in sorted(os.listdir(d)):
            p = os.path.join(d, e)
            if os.path.isdir(p):
                walk(p)
            elif p not in used:
                rows.append((roles[p][0], rel(p), roles[p][1]))
                used.add(p)

    walk(BASE)
    inv = pd.DataFrame(rows, columns=["role", "path", "note"])
    write_csv(inv, os.path.join(ADD, "input_inventory.csv"))
    return inv, roles


def inv_paths(inv, role):
    return [os.path.join(ROOT, p) for p in inv[inv.role == role].path]


# =========================================================================================
# 3. Readers: grid file, well data, well-to-block file, mapping
# =========================================================================================
def lit(x):
    return ast.literal_eval(x) if isinstance(x, str) else x


def read_mulgraph(path):
    """Vertices, grid columns (polygons), layers and the surface section, by section header."""
    L = open(path, errors="ignore").read().split("\n")
    heads = {l.strip(): i for i, l in enumerate(L) if l.strip() in ("VERTICES", "GRID", "CONNECTIONS", "LAYERS",
                                                                     "SURFA", "WELLS")}
    order = sorted(heads.items(), key=lambda kv: kv[1])
    sec = {}
    for k, (name, i) in enumerate(order):
        j = order[k + 1][1] if k + 1 < len(order) else len(L)
        sec[name] = [x for x in L[i + 1:j] if x.strip()]
    verts = {}
    for l in sec["VERTICES"]:
        verts[l[:3]] = (float(l[3:13]), float(l[13:23]))
    # GRID: 'name <flag> <nvert>' header followed by nvert vertex names
    cols, i, lines = {}, 0, sec["GRID"]
    while i < len(lines):
        head = lines[i].split()
        name = lines[i][:3]
        nv = int(head[-1])
        cols[name] = [lines[i + 1 + k].strip()[:3] for k in range(nv)]
        i += 1 + nv
    layers = {}
    for l in sec["LAYERS"]:
        p = l.split()
        layers[int(p[0])] = (float(p[1]), float(p[2]))
    surface = {}
    for l in sec["SURFA"]:
        p = l.split()
        surface[p[0]] = float(p[1])
    return verts, cols, layers, surface


def layer_table(layers):
    """Layer i (i>0): top = bottom of layer i-1; layer 0 gives the top surface."""
    return {i: dict(top=layers[i - 1][0], bottom=layers[i][0], centre=layers[i][1]) for i in sorted(layers) if i > 0}


def date_of(s):
    return dt.datetime.strptime(s, "%d/%m/%Y").date()


def read_well_json(path):
    d = json.load(open(path, encoding="utf-8"))
    w = dict(name=d["Well_Name"], aut2=d.get("AUT2_Name"), raw=d)
    w["fz"] = [list(lit(x)) for x in d["Feed_Zones"]["Data"]]
    w["fz_meta"] = {k: v for k, v in d["Feed_Zones"].items() if k not in ("Data",)}
    w["track"] = np.array([lit(x) for x in d["Well_Track"]["Data"]], float)
    w["design"] = [list(lit(x)) for x in d["Well_Design"]["Data"]]
    w["design_source"] = d["Well_Design"].get("Source", "")

    def surveys(key, unit):
        out = {}
        for ds, v in d.get(key, {}).items():
            z = np.array([float(k) for k in v["Data"].keys()])
            val = np.array([float(x) for x in v["Data"].values()])
            o = np.argsort(-z)
            out[date_of(ds)] = dict(z=z[o], v=val[o], source=v.get("Source", ""), comment=v.get("Comment", ""),
                                    stable=v.get("StableFlag", None))
        return out

    w["T_surveys"] = surveys("Downhole_Temperature", "T")
    w["P_surveys"] = surveys("Downhole_Pressure", "P")

    def series(key):
        v = d.get(key)
        if not v:
            return None
        return dict(dates=[date_of(k) for k in v["Data"].keys()], values=np.array(list(v["Data"].values()), float),
                    source=v.get("Source", ""), comment=v.get("Comment", ""))

    w["prod"] = series("Production_Mass_Flow")
    w["inj"] = series("Injection_Mass_Flow")
    w["inj_h"] = series("Injection_Enthalpy")
    w["prod_h"] = series("Transient_Enthalpy")
    w["transient"] = {}
    for zk, v in d.get("Transient_Pressure", {}).items():
        w["transient"][float(zk)] = dict(dates=[date_of(k) for k in v["Data"].keys()],
                                         values=np.array(list(v["Data"].values()), float),
                                         source=v.get("Source", ""), comment=v.get("Comment", ""))
    return w


def interp_down(z_query, z, v):
    """Linear interpolation on a profile stored top-down (z descending); NaN outside the survey range."""
    z_query = np.asarray(z_query, float)
    out = np.interp(z_query, z[::-1], v[::-1])
    out = np.where((z_query > z.max() + 1e-9) | (z_query < z.min() - 1e-9), np.nan, out)
    return out


# ----- geometry helpers -----
def point_seg_dist(px, py, ax, ay, bx, by):
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / L2))
    return math.hypot(px - (ax + t * dx), py - (ay + t * dy))


def dist_to_polyline(px, py, pts):
    pts = np.asarray(pts, float)
    if len(pts) == 1:
        return math.hypot(px - pts[0, 0], py - pts[0, 1])
    return min(point_seg_dist(px, py, *pts[i], *pts[i + 1]) for i in range(len(pts) - 1))


def mapping_features(path):
    """Features as {name: dict(kind, segments=[Nx2 arrays], closed)} using the names in the file."""
    files = sorted(glob.glob(os.path.join(path, "*.json"))) if os.path.isdir(path) else [path]
    m = {}
    for f in files:
        m.update(json.load(open(f, encoding="utf-8")))
    feats = {}
    for name, v in m.items():
        data = v["Data"]
        segs = []
        if data and isinstance(data[0][0], (list, tuple)):      # ragged list of segments
            segs = [np.array(s, float) for s in data]
        else:
            segs = [np.array(data, float)]
        kind = "fault" if name.lower().startswith("fault") else "water"
        feats[name] = dict(kind=kind, segs=segs, style=v.get("type", ""), colour=v.get("c", ""))
    return feats


def dist_to_feature(px, py, feat):
    return min(dist_to_polyline(px, py, s) for s in feat["segs"])


# =========================================================================================
# 4. Analysis
# =========================================================================================
class Ctx:
    pass


def quality_map(tags):
    """Map feed-zone quality tags to the weight_P factors of 5.5 (highest grade 1.0, intermediate 0.8,
    lowest grade or any tag that is not a quality grade 0.5)."""
    rank = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}
    found = sorted(set(tags))
    grades = sorted([t for t in found if t.upper() in rank], key=lambda t: -rank[t.upper()])
    out = {}
    for k, t in enumerate(grades):
        if k == 0:
            out[t] = 1.0
        elif k == len(grades) - 1:
            out[t] = 0.5
        else:
            out[t] = 0.8
    for t in found:
        out.setdefault(t, 0.5)
    return out, found


def run_isothermal(z, t, span=300.0, tol=5.0):
    """Mask of points inside a run of consecutive points spanning >= span m with all T within tol of each other."""
    n = len(z)
    mask = np.zeros(n, bool)
    for i in range(n):
        j = i
        while j + 1 < n and (max(t[i:j + 2]) - min(t[i:j + 2])) <= tol:
            j += 1
        if j > i and abs(z[i] - z[j]) >= span:
            mask[i:j + 1] = True
    return mask


def injectate_temperature(h_kjkg, p_bar):
    """Temperature [degC] at which IF97 Region 1 liquid at p has enthalpy h (bisection)."""
    lo, hi = 1.0, 200.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if region1(mid, p_bar)[1] < h_kjkg:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def analyse():
    C = Ctx()
    validate_if97()
    inv, roles = build_inventory()
    C.inv = inv
    geo_path = inv_paths(inv, "geometry_dat")[0]
    blk_path = inv_paths(inv, "well_blocks_json")[0]
    map_path = inv_paths(inv, "mapping")[0]
    grp_path = inv_paths(inv, "well_groups")[0]
    json_dir = inv_paths(inv, "well_data_json")[0]
    dat_dir = inv_paths(inv, "well_data_dat")[0]
    C.paths = dict(geo=geo_path, blk=blk_path, map=map_path, grp=grp_path, json=json_dir, dat=dat_dir)

    # ---- groups ----
    groups = json.load(open(grp_path, encoding="utf-8"))["Areas"]
    C.groups = {}
    for area, gd in groups.items():
        for gname, wl in gd.items():
            C.groups[f"{area}_{gname}"] = list(wl)
    C.wells_order = []
    for wl in C.groups.values():
        for w in wl:
            if w not in C.wells_order:
                C.wells_order.append(w)

    # ---- geometry ----
    verts, cols, layers, surface = read_mulgraph(geo_path)
    C.verts, C.cols, C.layers_raw, C.surface = verts, cols, layers, surface
    C.layers = layer_table(layers)
    C.grid_top = layers[0][0]
    C.col_poly = {c: np.array([verts[v] for v in vs]) for c, vs in cols.items()}
    C.col_centre = {c: p.mean(axis=0) for c, p in C.col_poly.items()}
    wb = json.load(open(blk_path, encoding="utf-8"))
    for w, d in wb.items():
        for k, v in list(d.items()):
            if isinstance(v, str):
                try:
                    d[k] = ast.literal_eval(v)
                except (ValueError, SyntaxError):
                    pass
    C.wb = wb
    C.features = mapping_features(map_path)

    # ---- wells ----
    wells = {}
    for f in sorted(glob.glob(os.path.join(json_dir, "*.json"))):
        w = read_well_json(f)
        wells[w["name"]] = w
    missing = [w for w in C.wells_order if w not in wells]
    if missing:
        log("WARNING: wells in the group list without data files:", missing)
    C.wells_order = [w for w in C.wells_order if w in wells]
    C.wells = wells
    cmap = plt.get_cmap("tab10")
    C.colour = {w: cmap(i % 10) for i, w in enumerate(C.wells_order)}

    quality_tags = [t for w in wells.values() for *_, t in w["fz"]]
    C.qmap, C.qfound = quality_map(quality_tags)

    # ---- per-well geometry, survey table, flags ----
    def blocks_for(well, z, tol=0.5):
        hit = []
        for b in wb[well]["Blocks"]:
            top, bot = max(b["interval"]), min(b["interval"])
            if bot - tol <= z <= top + tol:
                hit.append((b["block"], top, bot))
        if not hit:
            return None
        inside = [h for h in hit if h[2] < z <= h[1]]
        return (inside or hit)[0][0]

    C.blocks_for = blocks_for
    rows = []
    C.info = {}
    for w in C.wells_order:
        d = wells[w]
        tr = d["track"]
        wh = tr[0]
        top_block = wb[w]["Blocks"][0]["block"]
        col = top_block[:3]
        surf = surface.get(col, C.grid_top)
        surf_note = "" if col in surface else "column has no entry in the surface section; grid top used"
        dates = sorted(set(d["T_surveys"]) & set(d["P_surveys"]))
        sd = dates[0]
        T = d["T_surveys"][sd]
        P = d["P_surveys"][sd]
        same_z = len(T["z"]) == len(P["z"]) and np.allclose(T["z"], P["z"])
        z = T["z"]
        t = T["v"]
        p = P["v"] if same_z else interp_down(z, P["z"], P["v"])
        ts = tsat_c(p)
        bm = ts - t
        slotted = [(max(a, b), min(a, b)) for a, b, dia, sl, spec in d["design"] if sl]
        cased = [(max(a, b), min(a, b)) for a, b, dia, sl, spec in d["design"] if not sl]
        shoe = min(b for a, b in cased) if cased else np.nan
        fzs = [(max(a, b), min(a, b), q) for a, b, q in d["fz"]]
        info = dict(well=w, x=wh[0], y=wh[1], wellhead=wh[2], column=col, top_block=top_block, surface=surf,
                    surf_note=surf_note, survey_date=sd, n_dates=len(dates), shoe=shoe, cased=cased,
                    slotted=slotted, fz=fzs, track=tr, same_z=same_z)
        top_layer = next((i for i in sorted(C.layers) if C.layers[i]["bottom"] < surf), min(C.layers))
        info["top_layer"] = top_layer
        C.info[w] = info
        iso_open = np.zeros(len(z), bool)
        in_cased = np.array([any(b - 1e-6 <= zz <= a + 1e-6 for a, b in cased) for zz in z])
        open_idx = np.where(~in_cased)[0]
        if len(open_idx):
            m = run_isothermal(z[open_idx], t[open_idx])
            iso_open[open_idx[m]] = True
        for i in range(len(z)):
            zz = z[i]
            if len(tr) >= 2:
                o = np.argsort(tr[:, 2])
                xx = float(np.interp(zz, tr[o, 2], tr[o, 0]))
                yy = float(np.interp(zz, tr[o, 2], tr[o, 1]))
            else:
                xx, yy = wh[0], wh[1]
            blk = blocks_for(w, zz)
            in_fz = any(b - 1e-6 <= zz <= a + 1e-6 for a, b, q in fzs)
            fzi = next((k for k, (a, b, q) in enumerate(fzs) if b - 1e-6 <= zz <= a + 1e-6), None)
            rows.append(dict(
                well=w, x=round(xx, 1), y=round(yy, 1), elevation_mRL=zz, T_degC=t[i], P_bar=p[i],
                Tsat_degC=ts[i], boiling_margin_degC=bm[i], column=(blk[:3] if blk else ""),
                layer=(int(blk[3:]) if blk else np.nan), block=blk or "", in_feedzone=bool(in_fz),
                behind_casing=bool(in_cased[i]), above_model_surface=bool(zz > surf + 1e-6),
                isothermal_open=bool(iso_open[i]), fz_index=(fzi + 1 if fzi is not None else 0),
                near_wellhead=bool(zz >= wh[2] - 30.0),
                in_top_block=bool(blk is not None and int(blk[3:]) == top_layer)))
    pts = pd.DataFrame(rows)
    C.pts = pts

    # ---- injectate temperatures ----
    C.inj = {}
    for w in C.wells_order:
        d = wells[w]
        if d["inj"] is not None and d["inj_h"] is not None:
            h = float(np.median(d["inj_h"]["values"]))
            g = pts[pts.well == w]
            fzP = [float(interp_down(0.5 * (a + b), g.elevation_mRL.values, g.P_bar.values)) for a, b, q in C.info[w]["fz"]]
            pr = float(np.nanmedian(fzP)) if len(fzP) else 20.0
            C.inj[w] = dict(h=h, T=injectate_temperature(h, pr), p=pr)

    # ---- feed zones ----
    fzrows = []
    for w in C.wells_order:
        g = pts[pts.well == w]
        for k, (a, b, q) in enumerate(C.info[w]["fz"], 1):
            mid = 0.5 * (a + b)
            tt = float(interp_down(mid, g.elevation_mRL.values, g.T_degC.values))
            pp = float(interp_down(mid, g.elevation_mRL.values, g.P_bar.values))
            ts_ = float(tsat_c(pp))
            fzd = list(wb[w]["Feed_Zones"].values())[k - 1] if k - 1 < len(wb[w]["Feed_Zones"]) else {}
            blk = list(fzd.keys())[0] if fzd else (blocks_for(w, mid) or "")
            fzrows.append(dict(well=w, feedzone=k, top_mRL=a, bottom_mRL=b, mid_mRL=mid, quality=q, T_degC=tt,
                               P_bar=pp, Tsat_degC=ts_, boiling_margin_degC=ts_ - tt, block=blk,
                               thickness_m=a - b, q_factor=C.qmap[q],
                               all_blocks="/".join(fzd.keys()) if fzd else ""))
    C.fz = pd.DataFrame(fzrows)

    # ---- gauges and natural-state verdicts ----
    def first_flow(d):
        out = {}
        for key, nm in (("prod", "production"), ("inj", "injection")):
            s = d[key]
            if s is not None:
                pos = [dd for dd, v in zip(s["dates"], s["values"]) if v > 0]
                if pos:
                    out[nm] = min(pos)
        return out

    prod_firsts = [first_flow(wells[w]).get("production") for w in C.wells_order]
    prod_firsts = [x for x in prod_firsts if x]
    first_prod_field = min(prod_firsts)
    C.first_prod_field = first_prod_field
    inj_firsts = [first_flow(wells[w]).get("injection") for w in C.wells_order]
    inj_firsts = [x for x in inj_firsts if x]
    C.first_inj_field = min(inj_firsts) if inj_firsts else None
    C.first_flow = {w: first_flow(wells[w]) for w in C.wells_order}

    gauge, gz = [], {}
    for w in C.wells_order:
        d = wells[w]
        g = pts[pts.well == w]
        ff = C.first_flow[w]
        own_first = min(ff.values()) if ff else None
        rowd = dict(well=w, gauge_elev_mRL=np.nan, survey_P_at_gauge_bar=np.nan, first_gauge_date="",
                    first_gauge_P_bar=np.nan, diff_survey_minus_first_bar=np.nan,
                    own_first_flow_date=(own_first.isoformat() if own_first else ""),
                    own_flow_type=("/".join(ff.keys()) if ff else ""),
                    field_first_production_date=first_prod_field.isoformat(),
                    pre_production_gauge="", survey_date=C.info[w]["survey_date"].isoformat())
        if d["transient"]:
            ze, ser = sorted(d["transient"].items())[0]
            sp = float(interp_down(ze, g.elevation_mRL.values, g.P_bar.values))
            rowd.update(gauge_elev_mRL=ze, survey_P_at_gauge_bar=sp, first_gauge_date=ser["dates"][0].isoformat(),
                        first_gauge_P_bar=float(ser["values"][0]),
                        diff_survey_minus_first_bar=sp - float(ser["values"][0]))
            pre = own_first is None or ser["dates"][0] <= own_first
            rowd["pre_production_gauge"] = bool(pre and ser["dates"][0] <= first_prod_field)
            gz[w] = ser
        gauge.append(rowd)
    C.gauge = pd.DataFrame(gauge)
    C.gauge_series = gz

    verified = []
    for r in gauge:
        if r["pre_production_gauge"] is True and abs(r["diff_survey_minus_first_bar"]) <= PREC_P + 1e-9:
            verified.append(r["well"])
    C.verified = verified
    fzv = C.fz[C.fz.well.isin(verified)]
    if len(fzv) >= 2:
        b1, a1 = np.polyfit(fzv.mid_mRL.values, fzv.P_bar.values, 1)
    else:
        b1, a1 = 0.0, P_ATM
    C.trend = (a1, b1)
    resid = fzv.P_bar.values - (a1 + b1 * fzv.mid_mRL.values)
    C.trend_rms = float(np.sqrt(np.mean(resid ** 2))) if len(resid) else 0.0
    C.tol = max(3.0, 2.0 * C.trend_rms)
    C.fz["P_minus_trend_bar"] = C.fz.P_bar - (a1 + b1 * C.fz.mid_mRL)
    C.trend_fn = lambda zz: a1 + b1 * np.asarray(zz, float)

    verdict = {}
    for r in gauge:
        w = r["well"]
        fzw = C.fz[C.fz.well == w]
        maxdev = float(np.abs(fzw.P_minus_trend_bar).max()) if len(fzw) else np.nan
        on_trend = bool(len(fzw) and maxdev <= C.tol)
        g = pts[pts.well == w]
        has_gauge = not np.isnan(r["gauge_elev_mRL"])
        v, c, why, sig = None, None, "", False
        if w in C.inj:
            ti = C.inj[w]["T"]
            op = g[~g.behind_casing & ~g.above_model_surface]
            iso = op[op.isothermal_open]
            sig = bool(len(iso) and abs(iso.T_degC.mean() - ti) <= 5.0)
            fzT = float(fzw.T_degC.max()) if len(fzw) else float(g.T_degC.max())
            hotter = fzT > ti + 5.0
            if sig:
                v, c = "possibly disturbed", "low"
                why = (f"injector with an isothermal open-hole section within 5 degC of the injectate temperature "
                       f"({ti:.1f} degC from h = {C.inj[w]['h']:.1f} kJ/kg)")
            elif hotter:
                v, c = "probably natural state at best", "low"
                why = (f"injector; formation at the feed zones ({fzT:.1f} degC) is hotter than the injectate "
                       f"({ti:.1f} degC from h = {C.inj[w]['h']:.1f} kJ/kg), so injection cooling could be hidden; "
                       f"no isothermal section near the injectate temperature")
            else:
                why = (f"injector, but no isothermal section within 5 degC of the injectate temperature ({ti:.1f} degC) and the "
                       f"feed-zone temperature ({fzT:.1f} degC) is not hotter than the injectate: evidence against an "
                       f"injection-affected survey; ")
        if v is None:
            if r["pre_production_gauge"] is True and has_gauge and abs(r["diff_survey_minus_first_bar"]) <= PREC_P + 1e-9:
                v, c = "natural state", "high"
                why = (f"survey P at the gauge elevation ({r['gauge_elev_mRL']:.0f} mRL) matches the first gauge reading "
                       f"({r['first_gauge_date']}) to {abs(r['diff_survey_minus_first_bar']):.2f} bar (precision {PREC_P} bar), "
                       f"and the gauge starts no later than the first flow of the well ({r['own_first_flow_date']})")
            elif has_gauge and r["pre_production_gauge"] is True and on_trend:
                v, c = "probably natural state", "medium"
                why = (f"gauge test misses by {abs(r['diff_survey_minus_first_bar']):.2f} bar (> {PREC_P} bar; recorded as a data "
                       f"question) but feed-zone pressures lie on the verified trend (max deviation {maxdev:.1f} bar <= {C.tol:.1f} bar)")
            elif on_trend:
                v, c = "probably natural state", "medium"
                if has_gauge:
                    why += (f"no pre-production gauge reading (gauge starts {r['first_gauge_date']}, after the first flow of the "
                            f"well {r['own_first_flow_date']}); ")
                else:
                    why += "no gauge; "
                why += f"feed-zone pressures lie on the verified trend (max deviation {maxdev:.1f} bar <= {C.tol:.1f} bar)"
            else:
                v, c = "possibly disturbed", "low"
                why += f"feed-zone pressures are off the verified trend by up to {maxdev:.1f} bar (> {C.tol:.1f} bar)"
        verdict[w] = dict(verdict=v, confidence=c, reason=why, maxdev=maxdev, on_trend=on_trend, injector=w in C.inj,
                          inj_signature=sig)
    C.verdict = verdict
    pc = C.gauge.copy()
    pc["verdict"] = [verdict[w]["verdict"] for w in pc.well]
    pc["confidence"] = [verdict[w]["confidence"] for w in pc.well]
    pc["reason"] = [verdict[w]["reason"] for w in pc.well]
    pc["fz_max_dev_from_trend_bar"] = [verdict[w]["maxdev"] for w in pc.well]
    C.pc = pc
    return C


# =========================================================================================
# 5. Tables
# =========================================================================================
def loc_weight(r):
    """Location weight of 5.5 (first matching row) and its reason."""
    if r.above_model_surface:
        return 0.0, "above the model surface"
    if r.near_wellhead:
        return 0.3, "within 30 m of the wellhead"
    if r.in_top_block:
        return 0.3, "in the top block of its column"
    if r.in_feedzone:
        return 1.0, "feed zone"
    if r.behind_casing:
        return 0.5, "behind casing"
    if r.isothermal_open:
        return 0.5, "open hole within an isothermal section (possible internal flow)"
    return 0.8, "open hole"


def build_targets(C):
    pts = C.pts.copy()
    pts["Tsat_degC"] = pts.Tsat_degC.round(2)
    pts["boiling_margin_degC"] = pts.boiling_margin_degC.round(2)
    # the margin is Tsat - T; recompute from the rounded values so that the two columns are consistent
    pts["boiling_margin_degC"] = (pts.Tsat_degC - pts.T_degC).round(2)
    conf = {w: C.verdict[w]["confidence"] for w in C.wells_order}
    ver = {w: C.verdict[w]["verdict"] for w in C.wells_order}
    pts["ns_confidence"] = pts.well.map(conf)
    pts["ns_verdict"] = pts.well.map(ver)
    ws, wr, wp = [], [], []
    for _, r in pts.iterrows():
        lw, why = loc_weight(r)
        cf = CONF_FACTOR[r.ns_confidence]
        ws.append(round(lw * cf, 4))
        wr.append(f"{why} (location {lw:g} x confidence {cf:g})")
        if r.in_feedzone and not r.above_model_surface:
            q = C.fz[(C.fz.well == r.well) & (C.fz.feedzone == r.fz_index)].q_factor.iloc[0]
            wp.append(round(q * cf, 4))
        else:
            wp.append(0.0)
    pts["weight"], pts["weight_reason"], pts["weight_P"] = ws, wr, wp
    lead = ["well", "x", "y", "elevation_mRL", "T_degC", "P_bar", "Tsat_degC", "boiling_margin_degC", "column", "layer",
            "block", "in_feedzone", "behind_casing", "above_model_surface", "ns_confidence", "weight", "weight_reason",
            "weight_P"]
    extra = ["ns_verdict", "fz_index", "isothermal_open", "near_wellhead", "in_top_block"]
    pts["layer"] = pts.layer.astype("Int64")
    tg = pts[lead + extra].copy()
    C.targets = tg
    return tg


def fit_line(z, v):
    b, a = np.polyfit(np.asarray(z, float), np.asarray(v, float), 1)
    return b, a


def crossing(z, t, x):
    for i in range(len(t) - 1):
        if t[i] < x <= t[i + 1]:
            return z[i] + (x - t[i]) * (z[i + 1] - z[i]) / (t[i + 1] - t[i])
    return float("nan")


def build_fz(C):
    """feedzones.csv with weights."""
    fz = C.fz.copy()
    fz["Tsat_degC"] = fz.Tsat_degC.round(2)
    fz["T_degC"] = fz.T_degC.round(2)
    fz["P_bar"] = fz.P_bar.round(2)
    fz["boiling_margin_degC"] = (fz.Tsat_degC - fz.T_degC).round(2)
    fz["ns_confidence"] = fz.well.map({w: C.verdict[w]["confidence"] for w in C.wells_order})
    cf = fz.ns_confidence.map(CONF_FACTOR)
    fz["weight_T"] = (1.0 * cf).round(4)
    fz["weight_P"] = (fz.q_factor * cf).round(4)
    lead = ["well", "feedzone", "top_mRL", "bottom_mRL", "mid_mRL", "quality", "T_degC", "P_bar", "Tsat_degC",
            "boiling_margin_degC", "block", "ns_confidence", "weight_T", "weight_P"]
    fz["P_minus_trend_bar"] = fz.P_minus_trend_bar.round(2)
    fz = fz[lead + ["thickness_m", "q_factor", "P_minus_trend_bar", "all_blocks"]]
    C.fz_out = fz
    return fz


def expected_summary(g, fz_bottom):
    """Per-well characteristics from one well's target rows (CLAUDE.md 4.4)."""
    g = g.sort_values("elevation_mRL", ascending=False)
    g = g.groupby("elevation_mRL", sort=False).mean(numeric_only=True).reset_index()
    z, t = g.elevation_mRL.to_numpy(float), g.T_degC.to_numpy(float)
    p, bm = g.P_bar.to_numpy(float), g.boiling_margin_degC.to_numpy(float)
    e = {}
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
    e["inversion_degC"], e["inversion_top_mRL"], e["inversion_base_mRL"] = inv
    liq = (p >= 1.5) & ((bm >= 0) | (np.isnan(bm) & (p > P_CRIT)))
    zl, tl, pl, bml = z[liq], t[liq], p[liq], bm[liq]
    e["shallow_gradient_degC_per_km"] = np.nan
    e["shallow_fit_n"] = 0
    if len(zl):
        w = zl >= zl[0] - 300
        e["shallow_fit_n"] = int(w.sum())
        if w.sum() >= 3:
            e["shallow_gradient_degC_per_km"] = -fit_line(zl[w], tl[w])[0] * 1000
    bounds = [b for b in (inv[2], fz_bottom) if not np.isnan(b)]
    sel = z <= min(bounds) if bounds else np.ones_like(z, bool)
    e["deep_gradient_degC_per_km"] = np.nan
    e["deep_fit_n"] = int(sel.sum())
    e["deep_fit_span_m"] = float(z[sel].max() - z[sel].min()) if sel.sum() else 0.0
    if sel.sum() >= 3 and e["deep_fit_span_m"] >= 300:
        e["deep_gradient_degC_per_km"] = -fit_line(z[sel], t[sel])[0] * 1000
    e["z_100C_mRL"], e["z_200C_mRL"] = crossing(z, t, 100.0), crossing(z, t, 200.0)
    ok = ~np.isnan(bml)
    if ok.any():
        j = int(np.argmin(np.where(ok, bml, np.inf)))
        e["min_boiling_margin_degC"], e["min_boiling_margin_elev_mRL"] = bml[j], zl[j]
    else:
        e["min_boiling_margin_degC"] = e["min_boiling_margin_elev_mRL"] = np.nan
    e["water_level_mRL"] = np.nan
    e["water_level_extrapolation_m"] = np.nan
    if len(zl) >= 3:
        b, a = fit_line(zl[:3], pl[:3])
        if b != 0:
            e["water_level_mRL"] = (P_ATM - a) / b
            e["water_level_extrapolation_m"] = e["water_level_mRL"] - z[0]
    e["near_boiling"] = bool(ok.any() and np.nanmin(bml) < 5)
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
    return e


def build_summary(C):
    tg, fz = C.targets, C.fz_out
    rows = []
    for w in C.wells_order:
        g = tg[tg.well == w]
        fzw = fz[fz.well == w]
        e = expected_summary(g, fzw.bottom_mRL.min() if len(fzw) else np.nan)
        # pivot: highest quality factor, then thickest, then deepest
        piv = fzw.assign(q=fzw.q_factor, th=(fzw.top_mRL - fzw.bottom_mRL).abs()).sort_values(
            ["q", "th", "mid_mRL"], ascending=[False, False, True]).iloc[0]
        e["pivot_elev_mRL"], e["pivot_P_bar"] = piv.mid_mRL, piv.P_bar
        e["pivot_feedzone"] = int(piv.feedzone)
        e["ns_verdict"] = C.verdict[w]["verdict"]
        e["ns_confidence"] = C.verdict[w]["confidence"]
        e["well"] = w
        regnote = ""
        rows.append(e)
    S = pd.DataFrame(rows)
    # regime notes
    notes = []
    for _, r in S.iterrows():
        n = []
        if r.regime == "upflow" and r.near_boiling:
            n.append(f"near-boiling: minimum margin {r.min_boiling_margin_degC:.1f} degC at {r.min_boiling_margin_elev_mRL:.0f} mRL, above the 1 degC two-phase threshold")
        if r.regime == "recharge/downflow":
            n.append("Tmax < 100 degC: cold column; rule 2")
        if r.regime == "outflow":
            n.append(f"inversion of {r.inversion_degC:.1f} degC between {r.inversion_top_mRL:.0f} and {r.inversion_base_mRL:.0f} mRL")
        notes.append("; ".join(n))
    S["regime_note"] = notes
    lead = ["well", "Tmax_degC", "Tmax_elev_mRL", "T_at_TD_degC", "TD_mRL", "inversion_degC", "inversion_top_mRL",
            "inversion_base_mRL", "shallow_gradient_degC_per_km", "deep_gradient_degC_per_km", "z_100C_mRL",
            "z_200C_mRL", "min_boiling_margin_degC", "min_boiling_margin_elev_mRL", "water_level_mRL",
            "pivot_elev_mRL", "pivot_P_bar", "regime", "ns_verdict", "ns_confidence"]
    extra = ["near_boiling", "regime_note", "deep_fit_n", "deep_fit_span_m", "shallow_fit_n",
             "water_level_extrapolation_m", "pivot_feedzone"]
    S = S[lead + extra]
    num = S.select_dtypes("number").columns
    S[num] = S[num].round(2)
    C.summary = S
    return S


def build_layer_centres(C):
    rows = []
    for w in C.wells_order:
        g = C.targets[C.targets.well == w]
        z, t, p = g.elevation_mRL.values, g.T_degC.values, g.P_bar.values
        for lay, d in C.layers.items():
            c = d["centre"]
            if z.min() - 1e-9 <= c <= z.max() + 1e-9:
                tt = float(interp_down(c, z, t))
                pp = float(interp_down(c, z, p))
                ts = float(tsat_c(pp))
                blk = C.blocks_for(w, c) or ""
                rows.append(dict(well=w, layer=lay, layer_centre_mRL=c, block=blk, T_degC=round(tt, 2), P_bar=round(pp, 2),
                                 Tsat_degC=round(ts, 2), boiling_margin_degC=round(ts - tt, 2) if ts == ts else np.nan))
    C.layer_centres = pd.DataFrame(rows)
    return C.layer_centres


def build_hydrostatic(C):
    """Window algorithm of 4.3. Returns the window table (band cased/open)."""
    rows = []
    for w in C.wells_order:
        g = C.targets[C.targets.well == w].sort_values("elevation_mRL", ascending=False)
        g = g[~g.above_model_surface]
        z, t, p = g.elevation_mRL.values, g.T_degC.values, g.P_bar.values
        bm = g.boiling_margin_degC.values
        elig = (p >= 1.5) & ((bm >= 10) | (np.isnan(bm) & (p > P_CRIT)))
        shoe = C.info[w]["shoe"]
        fz_tops = sorted([a for a, b, q in C.info[w]["fz"]], reverse=True)
        markers = sorted(set(fz_tops + ([shoe] if shoe == shoe else [])), reverse=True)
        i = 0
        while i < len(z):
            if not elig[i]:
                i += 1
                continue
            j = i
            while j + 1 < len(z) and elig[j + 1]:
                j += 1
            ztop, zend = z[i], z[j]
            wins, start = [], ztop
            while True:
                nxt = [m for m in markers if m < start - 1e-9 and m > zend + 1e-9 and start - m >= 300]
                if nxt:
                    wins.append([start, nxt[0]])
                    start = nxt[0]
                    continue
                if start - zend >= 300 or (not wins and False):
                    wins.append([start, zend])
                elif wins:
                    wins[-1][1] = zend
                break
            for a, b in wins:
                zz = np.arange(b, a + 0.5, 1.0)
                zz = np.unique(np.append(zz, [a, b]))
                tz = np.interp(zz, z[i:j + 1][::-1], t[i:j + 1][::-1])
                pz = np.interp(zz, z[i:j + 1][::-1], p[i:j + 1][::-1])
                rho = rho_liquid(tz, pz)
                mid_rho = 0.5 * (rho[1:] + rho[:-1])
                dp_if97 = float(np.sum(mid_rho * G * np.diff(zz)) / 1e5)
                dp_meas = float(np.interp(b, z[i:j + 1][::-1], p[i:j + 1][::-1]) - np.interp(a, z[i:j + 1][::-1], p[i:j + 1][::-1]))
                length = a - b
                rho_imp = dp_meas * 1e5 / (G * length)
                ratio = dp_meas / dp_if97
                unc = (math.sqrt(2) * PREC_P / 2) / dp_if97
                cased_len = max(0.0, min(a, shoe) - b) if False else None
                if shoe == shoe:
                    over_cased = max(0.0, a - max(b, shoe))
                    band = "cased" if over_cased >= 0.5 * length else "open"
                    straddles = bool(a > shoe + 1e-9 and b < shoe - 1e-9)
                else:
                    band, straddles = "open", False
                rows.append(dict(well=w, z_top_mRL=round(a, 1), z_bottom_mRL=round(b, 1), length_m=round(length, 1), band=band,
                                 dP_measured_bar=round(dp_meas, 2), dP_IF97_bar=round(dp_if97, 2),
                                 ratio_measured_to_IF97=round(ratio, 4), rho_implied_kgm3=round(rho_imp, 1),
                                 ratio_uncertainty=round(unc, 4), T_mean_degC=round(float(np.mean(tz)), 1),
                                 straddles_shoe=straddles, shoe_mRL=shoe,
                                 Tmax_well_degC=round(float(g.T_degC.max()), 1)))
            i = j + 1
    H = pd.DataFrame(rows)
    C.hydro = H
    # well x band, cold component
    cold = {w for w in C.wells_order if C.targets[C.targets.well == w].T_degC.max() < 100}
    C.cold_wells = sorted(cold)
    comp = {}
    for band in ("cased", "open"):
        c = H[(H.band == band) & H.well.isin(cold)]
        comp[band] = (float(c.ratio_measured_to_IF97.median()) if len(c) else np.nan, len(c))
    C.cold_comp = comp
    rows = []
    for (w, band), g in H.groupby(["well", "band"], sort=False):
        med = float(g.ratio_measured_to_IF97.median())
        rows.append(dict(well=w, band=band, n_windows=len(g), length_m=round(g.length_m.sum(), 1),
                         ratio_median=round(med, 4), cold_component=round(comp[band][0], 4) if comp[band][0] == comp[band][0] else np.nan,
                         excess_over_cold_component=round(med - comp[band][0], 4) if comp[band][0] == comp[band][0] else np.nan,
                         cold_well=w in cold))
    C.hydro_wb = pd.DataFrame(rows)
    allb = []
    for band in ("cased", "open", "all"):
        g = H if band == "all" else H[H.band == band]
        allb.append(dict(group="all wells", band=band, n_windows=len(g),
                         ratio_median=round(float(g.ratio_measured_to_IF97.median()), 4) if len(g) else np.nan,
                         ratio_min=round(float(g.ratio_measured_to_IF97.min()), 4) if len(g) else np.nan,
                         ratio_max=round(float(g.ratio_measured_to_IF97.max()), 4) if len(g) else np.nan))
        if band != "all":
            c = g[g.well.isin(cold)]
            allb.append(dict(group="cold wells (Tmax < 100 degC)", band=band, n_windows=len(c),
                             ratio_median=round(float(c.ratio_measured_to_IF97.median()), 4) if len(c) else np.nan,
                             ratio_min=round(float(c.ratio_measured_to_IF97.min()), 4) if len(c) else np.nan,
                             ratio_max=round(float(c.ratio_measured_to_IF97.max()), 4) if len(c) else np.nan))
    C.hydro_summary = pd.DataFrame(allb)
    return H


def build_geometry_check(C):
    rows = []
    for w in C.wells_order:
        i = C.info[w]
        tr = i["track"]
        dev = float(np.max(np.hypot(tr[:, 0] - tr[0, 0], tr[:, 1] - tr[0, 1])))
        d = {k: dist_to_feature(i["x"], i["y"], f) for k, f in C.features.items()}
        faults = {k: v for k, v in d.items() if C.features[k]["kind"] == "fault"}
        waters = {k: v for k, v in d.items() if C.features[k]["kind"] == "water"}
        nf = min(faults, key=faults.get)
        nw = min(waters, key=waters.get)
        # own point-in-polygon column search vs the well-to-block file at every survey point
        g = C.targets[C.targets.well == w]
        bad = 0
        for _, r in g.iterrows():
            hit = [c for c, poly in C.col_poly.items() if MPath(poly).contains_point((r.x, r.y))]
            if r.column not in hit and hit:
                bad += 1
        wh_hit = [c for c, poly in C.col_poly.items() if MPath(poly).contains_point((i["x"], i["y"]))]
        row = dict(well=w, wellhead_mRL=i["wellhead"], model_surface_mRL=i["surface"],
                   wellhead_minus_surface_m=round(i["wellhead"] - i["surface"], 1), max_horizontal_deviation_m=round(dev, 1),
                   nearest_fault=nf, dist_nearest_fault_m=round(faults[nf], 1), nearest_water_body=nw,
                   dist_nearest_water_body_m=round(waters[nw], 1))
        for k, v in d.items():
            row["dist_" + re.sub(r"\W+", "_", k) + "_m"] = round(v, 1)
        row.update(x=i["x"], y=i["y"], top_block=i["top_block"], column=i["column"],
                   column_from_point_in_polygon="/".join(wh_hit), pip_column_agrees=bool(i["column"] in wh_hit),
                   survey_points_pip_disagree=bad, survey_points=len(g),
                   first_track_point_mRL=float(tr[0, 2]), wellhead_equals_first_track_point=bool(abs(tr[0, 2] - i["wellhead"]) < 1e-6),
                   survey_top_mRL=float(g.elevation_mRL.max()), survey_bottom_mRL=float(g.elevation_mRL.min()),
                   track_bottom_mRL=float(tr[:, 2].min()))
        rows.append(row)
    C.geo = pd.DataFrame(rows)
    return C.geo


def find_dat_files(C):
    """Group the plain-column files by well number from their names (irregular names are recorded)."""
    out, irregular = {}, []
    for f in sorted(glob.glob(os.path.join(C.paths["dat"], "*.dat"))):
        b = os.path.basename(f)
        m = re.match(r"^(Well|Wel|Wl)(\d+)_(.*)\.dat$", b)
        if not m:
            irregular.append((b, "name does not match Well<N>_<kind>.dat"))
            continue
        if m.group(1) != "Well":
            irregular.append((b, f"prefix '{m.group(1)}' instead of 'Well'"))
        out.setdefault(int(m.group(2)), {})[m.group(3)] = f
    return out, irregular


def build_crosscheck(C):
    """JSON vs plain-column files, track/design/feed-zone/survey consistency."""
    rows = []

    def add(w, item, status, **kw):
        d = dict(well=w, item=item, status=status)
        d.update(kw)
        rows.append(d)

    dat, irregular = find_dat_files(C)
    C.dat_irregular = irregular
    ref = None
    for w in C.wells_order:
        d = C.wells[w]
        n = int(re.search(r"(\d+)$", w).group(1))
        files = dat.get(n, {})
        # surveys
        for kind, key, sv in (("downhole_temperature", "T", d["T_surveys"]), ("downhole_pressure", "P", d["P_surveys"])):
            f = files.get(kind)
            sd = C.info[w]["survey_date"]
            if f is None:
                add(w, f"survey {key}: JSON vs .dat", "no .dat file")
                continue
            a = np.loadtxt(f)
            a = a.reshape(-1, 2)
            o = np.argsort(-a[:, 1])
            a = a[o]
            js = sv[sd]
            same = len(a) == len(js["z"])
            md = float(np.max(np.abs(a[:, 0] - js["v"]))) if same else np.nan
            mz = float(np.max(np.abs(a[:, 1] - js["z"]))) if same else np.nan
            add(w, f"survey {key}: JSON vs .dat", "OK" if same and md < 1e-6 and mz < 1e-6 else "MISMATCH",
                n_json=len(js["z"]), n_dat=len(a), max_abs_diff=md, max_elev_diff=mz, file=os.path.basename(f),
                note=f"survey date {sd.isoformat()}; {len(sv)} survey date(s) in the JSON")
        # transient pressure
        for ze, ser in d["transient"].items():
            f = None
            for k, v in files.items():
                if "transient_pressure" in k and str(float(ze)) in k:
                    f = v
            if f is None:
                add(w, f"transient P at {ze:.0f} mRL: JSON vs .dat", "no .dat file")
                continue
            a = np.loadtxt(f).reshape(-1, 2)
            same = len(a) == len(ser["values"])
            md = float(np.max(np.abs(a[:, 1] - ser["values"]))) if same else np.nan
            r0 = ser["dates"][0] - dt.timedelta(days=round(a[0, 0] * 365.25))
            ref = r0 if ref is None else ref
            td = float(np.max(np.abs(a[:, 0] - np.array([(x - r0).days / 365.25 for x in ser["dates"]])))) if same else np.nan
            add(w, f"transient P at {ze:.0f} mRL: JSON vs .dat", "OK" if same and md < 1e-6 else "MISMATCH", n_json=len(ser["values"]),
                n_dat=len(a), max_abs_diff=md, max_time_diff_years=td, file=os.path.basename(f),
                note=f"reference date of the .dat time axis (years since): {r0.isoformat()}; first reading {ser['dates'][0].isoformat()}")
        # flows / enthalpy
        for kind, key in (("production_mass_flow", "prod"), ("transient_enthalpy", "prod_h"),
                          ("injection_mass_flow", "inj"), ("injection_enthalpy", "inj_h")):
            s = d[key]
            f = files.get(kind)
            if s is None and f is None:
                continue
            if s is None or f is None:
                add(w, f"{kind}: JSON vs .dat", "MISMATCH", note="present in one source only", file=os.path.basename(f) if f else "")
                continue
            a = np.loadtxt(f).reshape(-1, 2)
            sign_note = ""
            if key == "prod" and np.all(a[:, 1] <= 0):
                a[:, 1] = -a[:, 1]
                sign_note = "; the .dat stores production as negative (generation-rate sign convention), sign reversed for the comparison"
            same = len(a) == len(s["values"])
            md = float(np.max(np.abs(a[:, 1] - s["values"]))) if same else np.nan
            r0 = ref or dt.date(2014, 1, 1)
            td = float(np.max(np.abs(a[:, 0] - np.array([(x - r0).days / 365.25 for x in s["dates"]])))) if same else np.nan
            add(w, f"{kind}: JSON vs .dat", "OK" if same and md < 1e-6 and td < 0.05 else "CHECK", n_json=len(s["values"]),
                n_dat=len(a), max_abs_diff=md, max_time_diff_years=td, file=os.path.basename(f),
                note="time axis is years since 2014-01-01 (monthly steps in the .dat; dates in the JSON)" + sign_note)
        # geometry: survey and feed zones within track, within slotted interval, design vs track
        i = C.info[w]
        tr = i["track"]
        zt, zb = tr[:, 2].max(), tr[:, 2].min()
        g = C.targets[C.targets.well == w]
        out = g[(g.elevation_mRL > zt + 1e-6) | (g.elevation_mRL < zb - 1e-6)]
        add(w, "survey points within the well track", "OK" if not len(out) else "OUTSIDE", n_json=len(g),
            note=f"track {zt:.0f} to {zb:.0f} mRL; survey {g.elevation_mRL.max():.0f} to {g.elevation_mRL.min():.0f} mRL")
        fzo = [k + 1 for k, (a, b, q) in enumerate(i["fz"]) if a > zt + 1e-6 or b < zb - 1e-6]
        add(w, "feed zones within the well track", "OK" if not fzo else "OUTSIDE", note=f"feed zones outside: {fzo}" if fzo else "")
        fzs = [k + 1 for k, (a, b, q) in enumerate(i["fz"]) if not any(sb - 1e-6 <= b and a <= st + 1e-6 for st, sb in i["slotted"])]
        add(w, "feed zones within the slotted/open interval", "OK" if not fzs else "OUTSIDE",
            note=f"slotted {i['slotted']}; feed zones outside: {fzs}" if fzs else f"slotted {i['slotted']}")
        dtop = max(a for a, b, *_ in d["design"])
        dbot = min(b for a, b, *_ in d["design"])
        ok = abs(dtop - tr[0, 2]) < 1e-6 and abs(dbot - zb) < 1.0
        add(w, "casing/liner design consistent with the track", "OK" if ok else "INCONSISTENT",
            note=f"design {dtop:.0f} to {dbot:.0f} mRL; track {tr[0, 2]:.0f} to {zb:.0f} mRL; source: {C.wells[w]['design_source']}")
        cs = sorted([(max(a, b), min(a, b)) for a, b, dia, sl, spec in d["design"] if not sl])
        ls = sorted([(max(a, b), min(a, b)) for a, b, dia, sl, spec in d["design"] if sl])
        if cs and ls:
            ov = min(b for a, b in cs) - max(a for a, b in ls)
            add(w, "casing shoe vs liner top", "OK", note=f"shoe {min(b for a, b in cs):.0f} mRL, liner top {max(a for a, b in ls):.0f} mRL; overlap {ov:.0f} m")
        if i["wellhead"] != g.elevation_mRL.max():
            add(w, "survey starts at the wellhead", "INFO", note=f"wellhead {i['wellhead']:.0f} mRL; shallowest survey point {g.elevation_mRL.max():.0f} mRL")
        for src, nm in ((d["T_surveys"], "T"), (d["P_surveys"], "P")):
            for sd, sv in src.items():
                if sv["source"] or sv["comment"] or sv["stable"] is not True:
                    add(w, f"survey {nm} provenance", "INFO", note=f"source '{sv['source']}', comment '{sv['comment']}', StableFlag {sv['stable']}")
        add(w, "survey date label", "INFO", note=f"all surveys dated {i['survey_date'].isoformat()}; first flow {C.first_flow[w]}")
    for b, why in irregular:
        add("(field)", f"irregular .dat file name: {b}", "INFO", note=why, file=b)
    for w in C.wells_order:
        n = int(re.search(r"(\d+)$", w).group(1))
        for kind, f in dat.get(n, {}).items():
            if os.path.basename(f).startswith("Wel") and not os.path.basename(f).startswith("Well"):
                pass
    C.dat_ref = ref
    X = pd.DataFrame(rows)
    C.crosscheck = X
    return X


def build_bands(C):
    tg = C.targets[~C.targets.above_model_surface]
    tg = tg.assign(Pres=tg.P_bar - C.trend_fn(tg.elevation_mRL))
    tg = tg.assign(band_top=np.floor(tg.elevation_mRL / 500.0) * 500.0 + 500.0)
    # band = [top-500, top): bottom inclusive
    rows = []
    for w in C.wells_order:
        g = tg[tg.well == w]
        for bt, gg in g.groupby("band_top"):
            rows.append(dict(well=w, band_top_mRL=bt, band_bottom_mRL=bt - 500.0, n_points=len(gg), T_mean_degC=gg.T_degC.mean(),
                             P_mean_bar=gg.P_bar.mean(), P_minus_trend_mean=gg.Pres.mean()))
    B = pd.DataFrame(rows)
    B["order"] = B.well.map({w: i for i, w in enumerate(C.wells_order)})
    B = B.sort_values(["order", "band_top_mRL"], ascending=[True, False]).drop(columns="order").reset_index(drop=True)
    med = B.groupby("band_top_mRL").T_mean_degC.median()
    B["T_field_median_degC"] = B.band_top_mRL.map(med)
    B["T_minus_field_median_degC"] = B.T_mean_degC - B.T_field_median_degC
    B["P_minus_trend_bar"] = B.P_minus_trend_mean
    B = B[["well", "band_top_mRL", "band_bottom_mRL", "n_points", "T_mean_degC", "T_field_median_degC", "T_minus_field_median_degC",
           "P_mean_bar", "P_minus_trend_bar"]]
    B = B.round(2)
    C.bands = B
    nw = B.groupby("band_top_mRL").well.nunique()
    good = nw[nw >= 3].index
    deepest = min(good) if len(good) else min(B.band_top_mRL)
    C.band_tops = [b for b in sorted(B.band_top_mRL.unique(), reverse=True) if b >= deepest]
    return B


def common_elevations(C):
    els = [0, -250, -500, -750, -1000, -1250, -1500, -2000, -2500]
    rows = []
    for w in C.wells_order:
        g = C.targets[C.targets.well == w]
        for e in els:
            t = float(interp_down(e, g.elevation_mRL.values, g.T_degC.values))
            p = float(interp_down(e, g.elevation_mRL.values, g.P_bar.values))
            rows.append(dict(well=w, elevation_mRL=e, T_degC=round(t, 1) if t == t else np.nan,
                             P_bar=round(p, 1) if p == p else np.nan,
                             P_minus_trend_bar=round(p - float(C.trend_fn(e)), 1) if p == p else np.nan))
    C.common = pd.DataFrame(rows)
    return C.common


# =========================================================================================
# 6. Figures: house style, base layer, F1, F2, F3, F5, F7
# =========================================================================================
import matplotlib.patheffects as pe
from matplotlib.collections import LineCollection

FIGS = []          # records for figures_index.csv (standard figures only)


def style_setup(C):
    C.tmax_lim = float(math.ceil(C.targets.T_degC.max() / 50.0) * 50.0)
    C.zlim = (math.floor(C.targets.elevation_mRL.min() / 100.0) * 100.0,
              math.ceil(max(C.info[w]["wellhead"] for w in C.wells_order) / 100.0) * 100.0)
    xs = np.array([v[0] for v in C.verts.values()])
    ys = np.array([v[1] for v in C.verts.values()])
    C.extent = (xs.min(), xs.max(), ys.min(), ys.max())
    C.tcmap = plt.get_cmap("inferno")
    C.tlevels = np.arange(0, C.tmax_lim + 0.1, 10.0)
    C.tnorm = BoundaryNorm(C.tlevels, C.tcmap.N)
    C.tcmap_d = ListedColormap(C.tcmap(np.linspace(0.0, 0.95, len(C.tlevels) - 1)))
    C.tnorm_d = BoundaryNorm(C.tlevels, len(C.tlevels) - 1)
    C.dcmap = plt.get_cmap("RdBu_r")
    C.plevels = np.arange(-20, 21, 1.0)
    C.pnorm = BoundaryNorm(C.plevels, C.dcmap.N, extend="both")
    C.isotherms = list(np.arange(50, C.tmax_lim, 50))
    # plan-view limits: model extent with a 3 % margin
    dx, dy = C.extent[1] - C.extent[0], C.extent[3] - C.extent[2]
    C.xlim = (C.extent[0] - 0.03 * dx, C.extent[1] + 0.03 * dx)
    C.ylim = (C.extent[2] - 0.03 * dy, C.extent[3] + 0.03 * dy)
    # grid edges as segments
    segs = []
    for c, poly in C.col_poly.items():
        pp = np.vstack([poly, poly[:1]])
        segs += [pp[i:i + 2] for i in range(len(pp) - 1)]
    C.grid_segs = segs


def fmt_axes_en(ax, small=False):
    ax.set_aspect("equal")
    ax.ticklabel_format(style="plain", useOffset=False)
    ax.tick_params(labelsize=6 if small else 8)
    if not small:
        ax.set_xlabel("Easting (m)", fontsize=8)
        ax.set_ylabel("Northing (m)", fontsize=8)
    else:
        for lab in ax.get_xticklabels():
            lab.set_rotation(90)


def feature_label_pos(C, name, f):
    """In-view point on the feature farthest from every wellhead (so the label does not sit on a well)."""
    wx = np.array([[C.info[w]["x"], C.info[w]["y"]] for w in C.wells_order])
    best = None
    for sg in f["segs"]:
        for k in range(len(sg) - 1):
            for t in np.linspace(0, 1, 6, endpoint=False):
                p = sg[k] + t * (sg[k + 1] - sg[k])
                if not (C.xlim[0] + 300 < p[0] < C.xlim[1] - 1500 and C.ylim[0] + 600 < p[1] < C.ylim[1] - 300):
                    continue
                d = float(np.min(np.hypot(wx[:, 0] - p[0], wx[:, 1] - p[1])))
                if best is None or d > best[0]:
                    best = (d, p)
    return None if best is None else best[1]


def base_layer(C, ax, labels=True, small=False, highlight=None, grey_others=False, track_lw=0.8):
    """The one base layer of 8.1: grid, extent (dashed), faults (black, named), water (blue), wellheads, tracks."""
    ax.add_collection(LineCollection(C.grid_segs, colors="#d4d4d4", linewidths=0.3, zorder=1))
    e = C.extent
    ax.plot([e[0], e[1], e[1], e[0], e[0]], [e[2], e[2], e[3], e[3], e[2]], ls="--", color="0.45", lw=0.8, zorder=2)
    for name, f in C.features.items():
        if f["kind"] == "water":
            if f["style"] == "fill":
                poly = f["segs"][0]
                ax.fill(poly[:, 0], poly[:, 1], color="#7fb2e5", alpha=0.7, zorder=2, lw=0)
            else:
                for s in f["segs"]:
                    ax.plot(s[:, 0], s[:, 1], color="#2b7bc0", lw=1.3, zorder=2)
        else:
            for s in f["segs"]:
                ax.plot(s[:, 0], s[:, 1], color="k", lw=1.4, zorder=3)
        if labels:
            p = feature_label_pos(C, name, f)
            if p is not None:
                ax.text(p[0], p[1], name, fontsize=5.5 if small else 7, color="k" if f["kind"] == "fault" else "#1f5f99",
                        rotation=0, ha="left", va="bottom", zorder=6, clip_on=True,
                        path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    for w in C.wells_order:
        tr = C.info[w]["track"]
        col = "0.55"
        lw = track_lw
        if highlight is not None and w in (highlight if isinstance(highlight, (list, tuple, set)) else [highlight]):
            col, lw = C.colour[w], track_lw * 2.5
        ax.plot(tr[:, 0], tr[:, 1], color=col, lw=lw, zorder=4)
        ax.plot(tr[0, 0], tr[0, 1], "o", color="k", ms=2.5 if small else 4, zorder=5)
    ax.set_xlim(*C.xlim)
    ax.set_ylim(*C.ylim)
    fmt_axes_en(ax, small)
    L = 2000.0
    x0 = C.xlim[0] + 0.04 * (C.xlim[1] - C.xlim[0])
    y0 = C.ylim[0] + 0.04 * (C.ylim[1] - C.ylim[0])
    ax.plot([x0, x0 + L], [y0, y0], "k-", lw=2.5, zorder=7)
    ax.text(x0 + L / 2, y0 + 0.012 * (C.ylim[1] - C.ylim[0]), "2 km", ha="center", va="bottom", fontsize=5.5 if small else 7, zorder=7)


def label_wells(C, ax, fs=7, only=None):
    off = {w: (6, 4) for w in C.wells_order}
    for w in C.wells_order:
        if only is not None and w not in only:
            continue
        i = C.info[w]
        ax.annotate(w.replace("WELL-", "W"), (i["x"], i["y"]), xytext=off[w], textcoords="offset points", fontsize=fs,
                    fontweight="bold", color=C.colour[w], zorder=8, ha="left", va="bottom",
                    path_effects=[pe.withStroke(linewidth=2, foreground="white")])


def short(w):
    return w.replace("WELL-", "W")


def save_fig(fig, path):
    fig.savefig(path, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def add_index(fid, path, shows, wells, line, sources, claim):
    FIGS.append(dict(figure_id=fid, file=rel(path), shows=shows, wells=wells, elevation_or_line=line,
                     data_sources=sources, claim_type=claim))


# ----------------------------------------------------------------- F1
def fig_well_map(C, lines):
    fig, ax = plt.subplots(figsize=(9.5, 8.5))
    base_layer(C, ax)
    label_wells(C, ax, fs=9)
    for ln in lines:
        p0, p1 = ln["p_ext0"], ln["p_ext1"]
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#7b3294", lw=1.6, ls="-.", zorder=6)
        for p, lab in ((p0, ln["id"]), (p1, ln["id"] + "′")):
            ax.text(p[0], p[1], lab, fontsize=11, fontweight="bold", color="#7b3294", ha="center", va="center", zorder=9,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#7b3294", lw=0.8))
    handles = [Line2D([], [], color="#d4d4d4", lw=1, label="model grid columns"),
               Line2D([], [], color="0.45", ls="--", label="model extent"),
               Line2D([], [], color="k", lw=1.4, label="faults (named)"),
               Line2D([], [], color="#2b7bc0", lw=1.3, label="rivers / lake (blue)"),
               Line2D([], [], marker="o", color="k", ls="", ms=4, label="wellheads"),
               Line2D([], [], color="0.55", lw=1, label="well tracks (plan)"),
               Line2D([], [], color="#7b3294", ls="-.", lw=1.6, label="section lines (F4, F6)")]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8, frameon=False)
    ax.set_title("F1  Well map: wells, faults, water bodies and section lines (plan view, NZGM easting/northing in m)", fontsize=10)
    path = os.path.join(CALFIG, "well_map.png")
    save_fig(fig, path)
    add_index("F1", path, "Base map: model grid, extent, faults, water bodies, wellheads, tracks and section lines",
              ";".join(C.wells_order), "plan view; lines " + ",".join(l["id"] for l in lines),
              "gTE14515.dat (grid, surface); Mapping.json; well_data_json (tracks)", "obs")


# ----------------------------------------------------------------- F2
def tsat_curve_pts(g):
    """Tsat(P) at the survey points where P is below the critical pressure (no clamped curve above it)."""
    m = g.Tsat_degC.notna() & (g.P_bar < P_CRIT)
    return g[m]


def draw_T_panel(C, ax, w, detail=True, xmax=None):
    g = C.targets[C.targets.well == w]
    i = C.info[w]
    S = C.summary.set_index("well").loc[w]
    col = C.colour[w]
    ax.plot(g.T_degC, g.elevation_mRL, "-", color=col, lw=0.8, alpha=0.6, zorder=2)
    ax.plot(g.T_degC, g.elevation_mRL, "o", color=col, ms=3.5 if detail else 2.5, zorder=3, label="survey T (points)")
    ts = tsat_curve_pts(g)
    ax.plot(ts.Tsat_degC, ts.elevation_mRL, "-", color="0.25", lw=1.2, zorder=2, label="Tsat(P), IF97")
    # near-boiling intervals (< 5 degC)
    z = g.elevation_mRL.values
    nb = (g.boiling_margin_degC.values < 5) & (g.P_bar.values >= 1.5)
    k = 0
    first = True
    while k < len(z):
        if nb[k]:
            j = k
            while j + 1 < len(z) and nb[j + 1]:
                j += 1
            top = z[k] if k == 0 else 0.5 * (z[k] + z[k - 1])
            bot = z[j] if j == len(z) - 1 else 0.5 * (z[j] + z[j + 1])
            ax.axhspan(bot, top, color="#f4a6a6", alpha=0.55, zorder=0, label="boiling margin < 5 degC" if first else None)
            first = False
            k = j + 1
        else:
            k += 1
    for k, (a, b, q) in enumerate(i["fz"], 1):
        ax.axhspan(b, a, color="#6fb6e8", alpha=0.45, zorder=1, label="feed zone" if k == 1 else None)
        if detail:
            ax.text(0.99, 0.5 * (a + b), q, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=6.5, color="#0b4a7a")
    ax.axhline(i["shoe"], color="k", lw=1.0, ls="--", zorder=2)
    ax.text(0.01, i["shoe"], "shoe", transform=ax.get_yaxis_transform(), fontsize=7, va="bottom", ha="left")
    ax.axhline(i["surface"], color="#2e7d32", lw=0.9, ls=":", zorder=2)
    ax.text(0.01, i["surface"], f"model surface {i['surface']:.0f}", transform=ax.get_yaxis_transform(), fontsize=6.5,
            va="bottom", color="#2e7d32")
    ax.axhline(i["wellhead"], color="0.3", lw=0.9, ls=":", zorder=2)
    ax.text(0.99, i["wellhead"], f"wellhead {i['wellhead']:.0f}", transform=ax.get_yaxis_transform(), fontsize=6.5,
            va="top", ha="right", color="0.3")
    for lev, zc in ((100, S.z_100C_mRL), (200, S.z_200C_mRL)):
        if zc == zc:
            ax.plot([lev], [zc], marker="D", ms=6, mfc="white", mec="k", zorder=4, ls="")
            ax.annotate(f"{lev} °C", (lev, zc), xytext=(4, -9), textcoords="offset points", fontsize=6.5)
    ax.set_xlim(0, xmax or C.tmax_lim)
    ax.set_ylim(*C.zlim)
    ax.grid(alpha=0.25, lw=0.4)


def fig_profile_well(C, w):
    fig = plt.figure(figsize=(13, 7.8))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1, 1], wspace=0.38)
    axm = fig.add_subplot(gs[0])
    base_layer(C, axm, labels=False, small=True, highlight=w)
    for fn, f in C.features.items():
        pass
    i = C.info[w]
    axm.annotate(short(w), (i["x"], i["y"]), xytext=(6, 4), textcoords="offset points", fontsize=9, fontweight="bold",
                 color=C.colour[w], path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    axm.set_title("Locator (plan)", fontsize=8)
    axT = fig.add_subplot(gs[1])
    axP = fig.add_subplot(gs[2], sharey=axT)
    draw_T_panel(C, axT, w)
    axT.set_xlabel("Temperature (°C)", fontsize=9)
    axT.set_ylabel("Elevation (mRL)", fontsize=9)
    axT.set_title("Temperature", fontsize=9)
    # pressure panel
    g = C.targets[C.targets.well == w]
    S = C.summary.set_index("well").loc[w]
    fzw = C.fz_out[C.fz_out.well == w]
    axP.plot(g.P_bar, g.elevation_mRL, "-", color=C.colour[w], lw=0.8, alpha=0.6)
    axP.plot(g.P_bar, g.elevation_mRL, "o", color=C.colour[w], ms=3.5, label="survey P (points)")
    axP.plot(fzw.P_bar, fzw.mid_mRL, "s", color="#0b4a7a", ms=8, mfc="none", mew=1.8, label="feed-zone pressure")
    zz = np.array(C.zlim)
    axP.plot(C.trend_fn(zz), zz, color="0.4", ls="--", lw=1, label="verified feed-zone trend")
    axP.plot([S.pivot_P_bar], [S.pivot_elev_mRL], marker="*", color="gold", mec="k", ms=14, ls="", label="pivot (dominant feed zone)")
    if S.water_level_mRL == S.water_level_mRL:
        axP.plot([P_ATM], [S.water_level_mRL], marker="v", color="#2b7bc0", ms=8, ls="", mec="k",
                 label=f"water level {S.water_level_mRL:.0f} mRL (extrapolated)")
    axP.axhline(C.info[w]["shoe"], color="k", lw=1.0, ls="--")
    for k, (a, b, q) in enumerate(C.info[w]["fz"], 1):
        axP.axhspan(b, a, color="#6fb6e8", alpha=0.45, zorder=0)
    axP.set_xlim(0, max(300, math.ceil(C.targets.P_bar.max() / 50) * 50))
    axP.set_ylim(*C.zlim)
    axP.grid(alpha=0.25, lw=0.4)
    axP.set_xlabel("Pressure (bar)", fontsize=9)
    axP.set_title("Pressure", fontsize=9)
    plt.setp(axP.get_yticklabels(), visible=False)
    v = C.verdict[w]
    h1, l1 = axT.get_legend_handles_labels()
    h2, l2 = axP.get_legend_handles_labels()
    handles = h1 + h2
    labels = l1 + l2
    handles.append(Line2D([], [], ls="", label=f"survey date {C.info[w]['survey_date'].isoformat()}"))
    labels.append(f"survey date {C.info[w]['survey_date'].isoformat()}")
    handles.append(Line2D([], [], ls="", label=""))
    labels.append(f"natural state: {v['verdict']} ({v['confidence']} confidence)")
    fig.legend(handles, labels, loc="lower center", ncol=4, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.07))
    reg = S.regime
    fig.suptitle(f"F2  {w}: temperature and pressure profiles vs elevation (mRL; °C, bar). Regime: {reg}", fontsize=11, y=0.96)
    path = os.path.join(CALFIG, f"profile_{w}.png")
    save_fig(fig, path)
    add_index("F2", path, f"{w} temperature and pressure profiles with casing shoe, feed zones, surface, isotherm crossings, boiling margin",
              w, "profile (elevation axis " + f"{C.zlim[0]:.0f} to {C.zlim[1]:.0f} mRL)", "well_data_json (survey, design, feed zones); IF97", "obs")


def fig_profile_panels(C, gname, wl):
    n = len(wl)
    ncol = 4
    nrow = int(math.ceil((n + 2) / ncol))
    fig = plt.figure(figsize=(3.4 * ncol, 4.6 * nrow))
    gs = fig.add_gridspec(nrow, ncol, wspace=0.12, hspace=0.28)
    axes = []
    for k, w in enumerate(wl):
        ax = fig.add_subplot(gs[k // ncol, k % ncol], sharey=axes[0] if axes else None)
        draw_T_panel(C, ax, w, detail=False)
        ax.set_title(f"{w}  ({C.verdict[w]['verdict']}, {C.verdict[w]['confidence']})", fontsize=8, color=C.colour[w])
        if k % ncol == 0:
            ax.set_ylabel("Elevation (mRL)", fontsize=8)
        else:
            plt.setp(ax.get_yticklabels(), visible=False)
        ax.set_xlabel("T (°C)", fontsize=8)
        ax.tick_params(labelsize=7)
        axes.append(ax)
    # locator map spans the remaining cells of the last row
    k0 = n
    r, c0 = k0 // ncol, k0 % ncol
    axm = fig.add_subplot(gs[r:, c0:])
    base_layer(C, axm, labels=False, small=True)
    for w in wl:
        i = C.info[w]
        tr = i["track"]
        axm.plot(tr[:, 0], tr[:, 1], color=C.colour[w], lw=2.5, zorder=6)
        axm.plot(tr[0, 0], tr[0, 1], "o", color=C.colour[w], mec="k", ms=6, zorder=7)
        axm.annotate(short(w), (i["x"], i["y"]), xytext=(6, 4), textcoords="offset points", fontsize=8, fontweight="bold",
                     color=C.colour[w], path_effects=[pe.withStroke(linewidth=2, foreground="white")], zorder=9)
    axm.set_title("Locator (plan)", fontsize=9)
    h = [Line2D([], [], marker="o", color="k", ls="", ms=4, label="survey T (points)"), Line2D([], [], color="0.25", label="Tsat(P)"),
         Line2D([], [], color="k", ls="--", label="casing shoe"),
         patches.Patch(color="#6fb6e8", alpha=0.45, label="feed zone"), patches.Patch(color="#f4a6a6", alpha=0.55, label="boiling margin < 5 °C"),
         Line2D([], [], color="#2e7d32", ls=":", label="model surface"), Line2D([], [], color="0.3", ls=":", label="wellhead"),
         Line2D([], [], marker="D", mfc="white", mec="k", ls="", label="100 / 200 °C crossing")]
    fig.legend(handles=h, loc="lower center", ncol=8, fontsize=8, frameon=False, bbox_to_anchor=(0.5, 0.02))
    fig.suptitle(f"F2  Group '{gname}': temperature profiles, all panels on the same axes "
                 f"({C.zlim[0]:.0f} to {C.zlim[1]:.0f} mRL; 0 to {C.tmax_lim:.0f} °C)", fontsize=12, y=0.96)
    path = os.path.join(CALFIG, f"profile_panels_{re.sub(r'[^A-Za-z0-9]+', '_', gname)}.png")
    save_fig(fig, path)
    add_index("F2", path, f"Small multiples of temperature profiles for group {gname}", ";".join(wl),
              f"profiles, common axes {C.zlim[0]:.0f} to {C.zlim[1]:.0f} mRL", "well_data_json; IF97", "obs")


# ----------------------------------------------------------------- F3
T_CLASSES = [(10, 20, 30, 40)]


def class_sym(delta, edges):
    """(marker, size, colour) for a signed deviation; 'x' cross inside +-edges[0]."""
    a = abs(delta)
    cm = plt.get_cmap("RdBu_r")
    if a <= edges[0]:
        return "P", 60, "0.2"
    k = 1 + sum(a > e for e in edges[1:])    # 1..len(edges)
    size = 40 + 55 * k
    col = cm(0.88 if delta > 0 else 0.12)
    return "o", size, col


def band_positions(C, w, zmid):
    tr = C.info[w]["track"]
    if len(tr) >= 2:
        o = np.argsort(tr[:, 2])
        return float(np.interp(zmid, tr[o, 2], tr[o, 0])), float(np.interp(zmid, tr[o, 2], tr[o, 1]))
    return tr[0, 0], tr[0, 1]


def fig_band_maps(C, kind):
    """kind 'T' or 'P'. Symbol maps on the base layer (F3)."""
    B = C.bands
    col = "T_minus_field_median_degC" if kind == "T" else "P_minus_trend_bar"
    edges = (10, 20, 30, 40) if kind == "T" else (2, 4, 6, 8)
    unit = "°C" if kind == "T" else "bar"
    tops = C.band_tops
    n = len(tops)
    ncol = 3
    nrow = int(math.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.6 * ncol, 5.3 * nrow + 1.2))
    axes = np.atleast_1d(axes).ravel()
    for ax, bt in zip(axes, tops):
        base_layer(C, ax, labels=False, small=True, track_lw=0.4)
        zmid = bt - 250.0
        bb = B[B.band_top_mRL == bt]
        have = set(bb.well)
        for w in C.wells_order:
            x, y = band_positions(C, w, zmid)
            if w in have:
                d = float(bb[bb.well == w][col].iloc[0])
                m, s, c = class_sym(d, edges)
                ax.scatter([x], [y], s=s, marker=m, c=[c], edgecolors="k" if m == "o" else None, linewidths=0.6, zorder=8)
                ax.annotate(short(w), (x, y), xytext=(5, 4), textcoords="offset points", fontsize=6.5, color="k", zorder=9,
                            path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
            else:
                ax.scatter([x], [y], s=16, facecolors="none", edgecolors="0.45", linewidths=0.8, zorder=8)
        med = bb.T_field_median_degC.iloc[0] if len(bb) else float("nan")
        if kind == "T":
            ax.set_title(f"{bt:.0f} to {bt - 500:.0f} mRL; field median T = {med:.1f} °C ({len(have)} wells)", fontsize=8)
        else:
            ax.set_title(f"{bt:.0f} to {bt - 500:.0f} mRL; deviation from feed-zone trend ({len(have)} wells)", fontsize=8)
    for ax in axes[n:]:
        ax.axis("off")
    # legend (outside the data area)
    cm = plt.get_cmap("RdBu_r")
    lab = [f"within ±{edges[0]} {unit} (cross)"]
    h = [Line2D([], [], marker="P", color="0.2", ls="", ms=7)]
    names = [f"{edges[0]}–{edges[1]}", f"{edges[1]}–{edges[2]}", f"{edges[2]}–{edges[3]}", f"> {edges[3]}"]
    for k, nm in enumerate(names, 1):
        for sign, cc, word in ((1, 0.88, "higher" if kind == "P" else "hotter"), (-1, 0.12, "lower" if kind == "P" else "colder")):
            h.append(Line2D([], [], marker="o", color=cm(cc), mec="k", ls="", ms=3.2 + 2.2 * k))
            lab.append(f"{word} {nm} {unit}")
    h.append(Line2D([], [], marker="o", mfc="none", mec="0.45", ls="", ms=4))
    lab.append("no data in band (open circle)")
    fig.legend(h, lab, loc="lower center", ncol=5, fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.02))
    if kind == "T":
        ttl = ("F3  Temperature relative to the field median in 500 m elevation bands (plan view at the band mid-elevation, "
               "no interpolation; °C)")
    else:
        ttl = ("F3  Pressure relative to the hydrostatic trend of gauge-verified feed zones in 500 m elevation bands "
               "(plan view at band mid-elevation, no interpolation; bar)")
    fig.suptitle(ttl, fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0.04, 1, 0.97))
    path = os.path.join(CALFIG, "maps_temperature_bands.png" if kind == "T" else "maps_pressure_bands.png")
    save_fig(fig, path)
    add_index("F3", path, ("Temperature deviation from the field median" if kind == "T" else "Pressure deviation from the feed-zone trend")
              + " by elevation band", ";".join(C.wells_order), "bands " + ", ".join(f"{b:.0f} to {b - 500:.0f} mRL" for b in tops),
              "T_P_at_elevation_bands.csv", "obs")


# ----------------------------------------------------------------- F5
def fig_PT(C):
    fig, ax = plt.subplots(figsize=(9.5, 7.5))
    Tc = np.linspace(1, 373.9, 400)
    Ps = np.array([psat_mpa(t) * 10 for t in Tc])
    ax.plot(Tc, Ps, "k-", lw=1.5, label="IF97 saturation curve (Region 4)", zorder=2)
    ax.plot([373.946], [220.64], marker="*", ms=14, color="k", ls="", label="critical point (373.9 °C, 220.64 bar)", zorder=5)
    for w in C.wells_order:
        g = C.targets[C.targets.well == w]
        ax.plot(g.T_degC, g.P_bar, "o", ms=3.5, color=C.colour[w], alpha=0.8, label=short(w), zorder=3)
    for w in C.wells_order:
        f = C.fz_out[C.fz_out.well == w]
        ax.plot(f.T_degC, f.P_bar, "s", ms=9, mfc="none", mec=C.colour[w], mew=2, zorder=4)
    ax.plot([], [], "s", ms=9, mfc="none", mec="k", mew=2, label="feed zones (larger, open squares)")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Pressure (bar)")
    ax.set_xlim(0, 400)
    ax.set_ylim(0, 300)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1.0), fontsize=8, frameon=False)
    ax.set_title("F5  Pressure–temperature state of all survey points against the saturation curve (°C, bar)", fontsize=10)
    path = os.path.join(CALFIG, "PT_state.png")
    save_fig(fig, path)
    add_index("F5", path, "P vs T of every survey point and feed zone with the IF97 saturation curve to the critical point",
              ";".join(C.wells_order), "all elevations", "NS_calibration_targets.csv; IF97", "obs")


# ----------------------------------------------------------------- F7
def fig_gauge(C):
    wl = [w for w in C.wells_order if w in C.gauge_series]
    if not wl:
        return
    n = len(wl)
    ncol = 2
    nrow = int(math.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(7.2 * ncol, 4.6 * nrow + 0.8))
    axes = np.atleast_1d(axes).ravel()
    gz = C.gauge.set_index("well")
    for ax, w in zip(axes, wl):
        ser = C.gauge_series[w]
        r = gz.loc[w]
        d0 = dt.date(2013, 7, 1)
        d1 = dt.date(2025, 12, 31)
        d = C.wells[w]
        for key in ("prod", "inj"):
            s = d[key]
            if s is None:
                continue
            on = [(a, b) for a, b, v in zip(s["dates"][:-1], s["dates"][1:], s["values"][:-1]) if v > 0]
            if on:
                ax.axvspan(on[0][0], on[-1][1], color="0.85", zorder=0, label="well flowing (grey)")
                break
        ax.plot(ser["dates"], ser["values"], "o-", color=C.colour[w], ms=4, lw=1, label="gauge pressure (transient)", zorder=3)
        ax.plot([ser["dates"][0]], [ser["values"][0]], "o", ms=9, mfc="none", mec="k", mew=1.5, label="first reading", zorder=4)
        ax.axhspan(r.first_gauge_P_bar - PREC_P, r.first_gauge_P_bar + PREC_P, color="#f7c948", alpha=0.7, zorder=2,
                   label=f"first reading ±{PREC_P} bar")
        sd = C.info[w]["survey_date"]
        ax.axhline(r.survey_P_at_gauge_bar, color="#c0392b", lw=1.4, ls="--", zorder=3,
                   label=f"survey P at {r.gauge_elev_mRL:.0f} mRL = {r.survey_P_at_gauge_bar:.2f} bar (survey dated {sd.isoformat()})")
        own = r.own_first_flow_date
        if own:
            ax.axvline(dt.date.fromisoformat(own), color="#8e44ad", lw=1.3, ls="-", label=f"first {r.own_flow_type} of this well ({own})")
        fp = C.first_prod_field
        ax.axvline(fp, color="#16a085", lw=1.3, ls=":", label=f"first production in the field ({fp.isoformat()})")
        ax.set_xlim(d0, d1)
        ax.set_ylabel("Pressure (bar)", fontsize=8)
        ax.set_title(f"{w}: {C.verdict[w]['verdict']} ({C.verdict[w]['confidence']} confidence); survey − first reading = "
                     f"{r.diff_survey_minus_first_bar:+.2f} bar", fontsize=8.5)
        ax.grid(alpha=0.3)
        ax.tick_params(labelsize=7)
        ax.legend(loc="upper left", bbox_to_anchor=(0.0, -0.1), fontsize=6.5, frameon=False, ncol=1)
    for ax in axes[n:]:
        ax.axis("off")
    fig.suptitle("F7  Natural-state gauge check: transient pressure at the gauge elevation vs the survey pressure interpolated to it "
                 "(2014–2025; bar)", fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    path = os.path.join(CALFIG, "gauge_check.png")
    save_fig(fig, path)
    add_index("F7", path, "Gauge series vs survey pressure at the gauge elevation, production start lines and precision band",
              ";".join(wl), "gauge elevations " + ", ".join(f"{gz.loc[w].gauge_elev_mRL:.0f}" for w in wl) + " mRL",
              "well_data_json (Transient_Pressure, flows); pressure_consistency_check.csv", "obs")


# =========================================================================================
# 7. Sections (F4), interpolated companions (8.4), conceptual model (F6)
# =========================================================================================
from scipy.spatial import cKDTree


def clip_line_to_extent(p, d, e):
    """Parameter range [t0, t1] of the infinite line p + t d inside the rectangle e = (xmin, xmax, ymin, ymax)."""
    t0, t1 = -1e9, 1e9
    for k, (lo, hi) in enumerate(((e[0], e[1]), (e[2], e[3]))):
        if abs(d[k]) < 1e-12:
            continue
        a, b = (lo - p[k]) / d[k], (hi - p[k]) / d[k]
        t0, t1 = max(t0, min(a, b)), min(t1, max(a, b))
    return t0, t1


def make_lines(C):
    """Line A through the two hottest wells (by Tmax); line B through the hottest well, perpendicular to A.
    Both are extended to the model extent."""
    S = C.summary.set_index("well")
    hot = list(S.Tmax_degC.sort_values(ascending=False).index)
    lines = []

    def build(lid, p, d, max_off, desc):
        d = np.asarray(d, float)
        d /= np.hypot(*d)
        if abs(d[1]) >= abs(d[0]):
            d = d if d[1] > 0 else -d      # run south to north
        else:
            d = d if d[0] > 0 else -d      # run west to east
        t0, t1 = clip_line_to_extent(p, d, C.extent)
        p0, p1 = p + t0 * d, p + t1 * d
        ln = dict(id=lid, p=np.asarray(p, float), d=d, p_ext0=p0, p_ext1=p1, length=t1 - t0, max_offset=max_off, desc=desc)
        lines.append(ln)
        return ln

    a, b = hot[0], hot[1]
    pa = np.array([C.info[a]["x"], C.info[a]["y"]])
    pb = np.array([C.info[b]["x"], C.info[b]["y"]])
    build("A", pa, pb - pa, 1000.0, f"through the two hottest wells ({short(a)} and {short(b)})")
    dA = lines[0]["d"]
    build("B", pa, np.array([dA[1], -dA[0]]), 600.0, f"through the hottest well ({short(a)}), perpendicular to A")
    for ln in lines:
        members = []
        for w in C.wells_order:
            s, off = project(ln, C.info[w]["x"], C.info[w]["y"])
            if abs(off) <= ln["max_offset"]:
                members.append((w, s, off))
        ln["members"] = members
    return lines


def project(ln, x, y):
    v = np.array([x, y]) - ln["p_ext0"]
    s = float(v @ ln["d"])
    off = float(v[0] * (-ln["d"][1]) + v[1] * ln["d"][0])      # positive to the left of the direction of travel
    return s, off


def section_data(C, ln):
    """Projected profiles of the wells within max_offset, sorted along the line."""
    wells = []
    for w, s0, off in ln["members"]:
        g = C.targets[C.targets.well == w].sort_values("elevation_mRL")
        ss = np.array([project(ln, x, y)[0] for x, y in zip(g.x, g.y)])
        wells.append(dict(well=w, s_head=s0, offset=off, z=g.elevation_mRL.values, T=g.T_degC.values, P=g.P_bar.values,
                          s=ss, bm=g.boiling_margin_degC.values, s_mean=float(np.mean(ss))))
    wells.sort(key=lambda d: d["s_mean"])
    tree = cKDTree(np.array([C.col_centre[c] for c in C.cols]))
    names = list(C.cols)
    sg = np.arange(0, ln["length"] + 1, 50.0)
    xy = ln["p_ext0"][None, :] + sg[:, None] * ln["d"][None, :]
    idx = tree.query(xy)[1]
    surf = np.array([C.surface.get(names[i], C.grid_top) for i in idx])
    # surface features crossing the line
    cross = []
    for name, f in C.features.items():
        for seg in f["segs"]:
            for k in range(len(seg) - 1):
                a, b = seg[k], seg[k + 1]
                da, db = project(ln, *a), project(ln, *b)
                if da[1] * db[1] < 0:
                    t = da[1] / (da[1] - db[1])
                    cross.append((name, da[0] + t * (db[0] - da[0]), f["kind"]))
    cross = [c for c in cross if 0 <= c[1] <= ln["length"]]
    return dict(wells=wells, sg=sg, surf=surf, cross=cross, ln=ln)


def envelope_meshes(C, SD, kind):
    """Colour-filled data envelope between neighbouring wells (never below either well's deepest point or above either's
    shallowest). Returns a list of (S, Z, V) arrays masked above the model surface."""
    out = []
    W = SD["wells"]
    for a, b in zip(W[:-1], W[1:]):
        zlo, zhi = max(a["z"].min(), b["z"].min()), min(a["z"].max(), b["z"].max())
        if zhi - zlo < 5:
            continue
        zz = np.arange(zlo, zhi + 0.1, 5.0)
        f = np.linspace(0, 1, 41)
        sa, sb = np.interp(zz, a["z"], a["s"]), np.interp(zz, b["z"], b["s"])
        if kind == "T":
            va, vb = np.interp(zz, a["z"], a["T"]), np.interp(zz, b["z"], b["T"])
        else:
            va = np.interp(zz, a["z"], a["P"]) - C.trend_fn(zz)
            vb = np.interp(zz, b["z"], b["P"]) - C.trend_fn(zz)
        S = sa[:, None] * (1 - f)[None, :] + sb[:, None] * f[None, :]
        V = va[:, None] * (1 - f)[None, :] + vb[:, None] * f[None, :]
        Z = np.repeat(zz[:, None], len(f), axis=1)
        surf = np.interp(S, SD["sg"], SD["surf"])
        V = np.ma.masked_where(Z > surf, V)
        out.append((S, Z, V))
    return out


HALO = [pe.withStroke(linewidth=2.2, foreground="white")]


def draw_section_into(C, SD, ax, kind, size_in, interp=False):
    """Draw the common section background (F4 / F6 / companion) into ax. size_in = (width, height) of ax in inches."""
    ln = SD["ln"]
    W = SD["wells"]
    smin = max(0.0, min(w["s"].min() for w in W) - 1200)
    smax = min(ln["length"], max(w["s"].max() for w in W) + 1200)
    ax.set_xlim(smin, smax)
    ax.set_ylim(*C.zlim)
    # model surface, hatched no-data below it
    sg, surf = SD["sg"], SD["surf"]
    ax.fill_between(sg, C.zlim[0], surf, step="mid", facecolor="#ececec", edgecolor="#9a9a9a", hatch="////", lw=0, zorder=0)
    ax.step(sg, surf, where="mid", color="#2e7d32", lw=1.6, zorder=6)
    cm = C.tcmap_d if kind == "T" else C.dcmap
    norm = C.tnorm_d if kind == "T" else C.pnorm
    levels = C.tlevels if kind == "T" else C.plevels
    meshes = envelope_meshes(C, SD, kind)
    if not interp:
        for S, Z, V in meshes:
            ax.contourf(S, Z, V, levels=levels, cmap=cm, norm=norm, extend="neither" if kind == "T" else "both", zorder=2)
            ax.plot(S[:, 0], Z[:, 0], ":", color="k", lw=0.7, zorder=5)
            ax.plot(S[:, -1], Z[:, -1], ":", color="k", lw=0.7, zorder=5)
            ax.plot(S[0], Z[0], ":", color="k", lw=0.7, zorder=5)
            ax.plot(S[-1], Z[-1], ":", color="k", lw=0.7, zorder=5)
            try:
                if kind == "T":
                    cs1 = ax.contour(S, Z, V, levels=[l for l in C.isotherms if l not in (100, 200)], colors="k", linewidths=0.6, zorder=4)
                    cs2 = ax.contour(S, Z, V, levels=[100, 200], colors="k", linewidths=1.8, zorder=4)
                    for cs in (cs1, cs2):
                        for c in cs.collections if hasattr(cs, "collections") else []:
                            c.set_path_effects(HALO)
                    ax.clabel(cs2, fmt="%d °C", fontsize=7, inline=True)
                else:
                    ax.contour(S, Z, V, levels=[l for l in np.arange(-20, 21, 2) if l != 0], colors="k", linewidths=0.5, zorder=4)
                    ax.contour(S, Z, V, levels=[0], colors="k", linewidths=1.6, linestyles="--", zorder=4)
            except Exception:
                pass
    else:
        s_all = np.concatenate([w["s"] for w in W])
        z_all = np.concatenate([w["z"] for w in W])
        v_all = np.concatenate([w["T"] if kind == "T" else w["P"] - C.trend_fn(w["z"]) for w in W])
        tri = Triangulation(s_all, z_all)
        ax.tricontourf(tri, v_all, levels=levels, cmap=cm, norm=norm, extend="neither" if kind == "T" else "both", zorder=2)
        if kind == "T":
            cs1 = ax.tricontour(tri, v_all, levels=[l for l in C.isotherms if l not in (100, 200)], colors="k", linewidths=0.6, zorder=4)
            cs2 = ax.tricontour(tri, v_all, levels=[100, 200], colors="k", linewidths=1.8, zorder=4)
        # air above the model surface stays white
        ax.fill_between(sg, surf, C.zlim[1] + 100, step="mid", facecolor="white", edgecolor="none", zorder=3)
        ax.step(sg, surf, where="mid", color="#2e7d32", lw=1.6, zorder=6)
        # outline the F4 data envelope with a dashed line
        for S, Z, V in meshes:
            for arr in ((S[:, 0], Z[:, 0]), (S[:, -1], Z[:, -1]), (S[0], Z[0]), (S[-1], Z[-1])):
                ax.plot(arr[0], arr[1], "--", color="k", lw=0.9, zorder=5)
    # wells: tracks, points, shoe, feed zones, boiling margin stipple, labels
    for w in W:
        i = C.info[w["well"]]
        col = C.colour[w["well"]]
        ax.plot(w["s"], w["z"], "-", color="0.15", lw=0.8, zorder=5)
        vals = w["T"] if kind == "T" else w["P"] - C.trend_fn(w["z"])
        ax.scatter(w["s"], w["z"], c=vals, cmap=cm, norm=norm, s=22, edgecolors="k", linewidths=0.5, zorder=7)
        for a, b, q in i["fz"]:
            sa = np.interp([a, b], w["z"], w["s"])
            ax.plot(sa, [a, b], "-", color="#00bfff", lw=5, solid_capstyle="butt", zorder=6, path_effects=[pe.withStroke(linewidth=7, foreground="k")])
        if i["shoe"] == i["shoe"] and w["z"].min() <= i["shoe"] <= w["z"].max():
            ss = float(np.interp(i["shoe"], w["z"], w["s"]))
            ax.plot([ss - 60, ss + 60], [i["shoe"]] * 2, "-", color="k", lw=2.5, zorder=8)
        nb = w["bm"] < 5
        if nb.any():
            idx = np.where(nb)[0]
            lo = max(w["z"].min(), w["z"][idx.min()] - 0.5 * (w["z"][idx.min()] - w["z"][max(idx.min() - 1, 0)]))
            hi = min(w["z"].max(), w["z"][idx.max()] + 0.5 * (w["z"][min(idx.max() + 1, len(w["z"]) - 1)] - w["z"][idx.max()]))
            ss = float(np.interp(0.5 * (lo + hi), w["z"], w["s"]))
            ax.add_patch(patches.Rectangle((ss - 70, lo), 140, hi - lo, facecolor="none", edgecolor="k", hatch="....", lw=0.6, zorder=9))
        zb = w["z"].min()
        sb = float(w["s"][np.argmin(w["z"])])
        ax.annotate(f"{short(w['well'])}\n({w['offset']:+.0f} m)", (sb, zb), xytext=(0, -6), textcoords="offset points", ha="center",
                    va="top", fontsize=8, fontweight="bold", color=col, zorder=10, path_effects=HALO)
    # features crossing the line: labelled ticks on the top edge
    ztop = C.zlim[1]
    for name, s, kd in SD["cross"]:
        if smin <= s <= smax:
            ax.plot([s, s], [ztop - 40, ztop], "-", color="k" if kd == "fault" else "#2b7bc0", lw=2.2, zorder=10, clip_on=False)
            ax.text(s, ztop + 20, name, rotation=90, ha="center", va="bottom", fontsize=6.5, color="k" if kd == "fault" else "#1f5f99", clip_on=False)
    ax.set_xlabel(f"Distance along line {ln['id']} from {ln['id']} (m)", fontsize=9)
    ax.set_ylabel("Elevation (mRL)", fontsize=9)
    ax.grid(alpha=0.2, lw=0.4)
    xr = smax - smin
    yr = C.zlim[1] - C.zlim[0]
    ve = (xr / size_in[0]) / (yr / size_in[1]) ** 1
    ve = 1.0 / ve if False else (xr / size_in[0]) / (yr / size_in[1])
    ve = 1.0 / ve
    meta = dict(smin=smin, smax=smax, ve=ve, cm=cm, norm=norm)
    return meta


def section_axes(C, SD, kind, interp=False):
    fig = plt.figure(figsize=(15, 8.2))
    ax = fig.add_axes([0.06, 0.12, 0.70, 0.73])
    meta = draw_section_into(C, SD, ax, kind, (0.70 * 15, 0.73 * 8.2), interp)
    return fig, ax, meta


def section_legend(C, fig, kind, interp=False):
    cax = fig.add_axes([0.80, 0.55, 0.018, 0.33])
    sm = plt.cm.ScalarMappable(norm=C.tnorm_d if kind == "T" else C.pnorm, cmap=C.tcmap_d if kind == "T" else C.dcmap)
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("Temperature (°C), 10 °C steps" if kind == "T" else "P − feed-zone trend (bar), 1 bar steps; red above trend", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    if kind == "T":
        cb.set_ticks(np.arange(0, C.tmax_lim + 1, 50))
    else:
        cb.set_ticks(np.arange(-20, 21, 5))
    h = [Line2D([], [], marker="o", color="w", mfc="0.5", mec="k", ms=5, label="survey points (coloured on the same scale)"),
         patches.Patch(facecolor="#ececec", edgecolor="#9a9a9a", hatch="////", label="no data (hatched grey)"),
         Line2D([], [], color="#2e7d32", lw=1.6, label="model surface"),
         Line2D([], [], color="#00bfff", lw=5, label="feed zone"),
         Line2D([], [], color="k", lw=2.5, label="casing shoe"),
         patches.Patch(facecolor="none", edgecolor="k", hatch="....", label="boiling margin < 5 °C (stipple)"),
         Line2D([], [], color="k", ls=":", lw=0.9, label="data envelope" + (" (F4 outline)" if interp else ""))]
    if kind == "T":
        h += [Line2D([], [], color="k", lw=0.6, label="isotherms every 50 °C"), Line2D([], [], color="k", lw=1.8, label="100 and 200 °C (heavy)")]
    else:
        h += [Line2D([], [], color="k", lw=0.5, label="contours every 2 bar"), Line2D([], [], color="k", lw=1.6, ls="--", label="zero deviation")]
    h += [Line2D([], [], color="k", lw=2.2, label="fault crossing (tick, top edge)"), Line2D([], [], color="#2b7bc0", lw=2.2, label="river / lake crossing")]
    fig.legend(handles=h, loc="upper left", bbox_to_anchor=(0.78, 0.50), fontsize=7.5, frameon=False)


def fig_section(C, SD, kind, interp=False):
    ln = SD["ln"]
    fig, ax, meta = section_axes(C, SD, kind, interp)
    section_legend(C, fig, kind, interp)
    p0, p1 = ln["p_ext0"], ln["p_ext1"]
    ends = (f"{ln['id']} ({p0[0]:.0f} E, {p0[1]:.0f} N) to {ln['id']}′ ({p1[0]:.0f} E, {p1[1]:.0f} N)")
    what = "Temperature" if kind == "T" else "Pressure deviation from the feed-zone hydrostatic trend"
    unit = "°C" if kind == "T" else "bar"
    pre = "[interpolated, not data] " if interp else "F4  "
    if interp:
        sub = ("Method: linear interpolation on a Delaunay triangulation of all survey points (distance along line, elevation); "
               "no extrapolation: colour only inside the convex hull of the data. Colours outside the dashed envelope are not constrained by "
               "neighbouring wells.")
    else:
        sub = ("Colour only inside the data envelope (between neighbouring wells, within both wells' survey range); everything else is hatched.")
    import textwrap
    ttl = (f"{pre}{what} along section line {ends}; wells within {ln['max_offset']:.0f} m (offset in brackets); {unit}; "
           f"vertical exaggeration ×{meta['ve']:.1f}")
    fig.suptitle("\n".join(textwrap.wrap(ttl, 135)), fontsize=10, y=0.995, va="top")
    fig.text(0.06, 0.885, "\n".join(textwrap.wrap(sub, 190)), fontsize=7.5, color="0.25", va="bottom")
    name = f"section_{ln['id']}" + ("_P" if kind == "P" else "") + ("_interpolated" if interp else "") + ".png"
    path = os.path.join(ADDFIG if interp else CALFIG, name)
    fig.savefig(path, dpi=DPI, facecolor="white")
    plt.close(fig)
    if not interp:
        add_index("F4", path, f"{what} section along line {ln['id']} with data envelope, isotherms/contours, shoes, feed zones, surface features",
                  ";".join(w["well"] for w in SD["wells"]), f"line {ends}; VE x{meta['ve']:.1f}",
                  "NS_calibration_targets.csv; Mapping.json; gTE14515.dat", "obs")
    return meta


# ----------------------------------------------------------------- interpolated plan views (8.4)
def fig_band_maps_interp(C, kind):
    B = C.bands
    col = "T_minus_field_median_degC" if kind == "T" else "P_minus_trend_bar"
    tops = C.band_tops
    n = len(tops)
    ncol = 3
    nrow = int(math.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(5.6 * ncol, 5.3 * nrow + 1.4))
    axes = np.atleast_1d(axes).ravel()
    if kind == "T":
        levels = np.arange(-150, 151, 10.0)
        norm = BoundaryNorm(levels, C.dcmap.N, extend="both")
    else:
        levels = C.plevels
        norm = C.pnorm
    e = C.extent
    for ax, bt in zip(axes, tops):
        base_layer(C, ax, labels=False, small=True, track_lw=0.4)
        ax.fill([e[0], e[1], e[1], e[0]], [e[2], e[2], e[3], e[3]], facecolor="#ececec", edgecolor="none", hatch="////", zorder=0.5)
        bb = B[B.band_top_mRL == bt]
        zmid = bt - 250.0
        xs, ys, vs, names = [], [], [], []
        for w in bb.well:
            x, y = band_positions(C, w, zmid)
            xs.append(x), ys.append(y), vs.append(float(bb[bb.well == w][col].iloc[0])), names.append(w)
        xs, ys, vs = map(np.array, (xs, ys, vs))
        ok = False
        if len(xs) >= 3:
            try:
                tri = Triangulation(xs, ys)
                ax.tricontourf(tri, vs, levels=levels, cmap=C.dcmap, norm=norm, extend="both", zorder=1.5)
                ok = True
            except Exception:
                ok = False
        ax.scatter(xs, ys, c=vs, cmap=C.dcmap, norm=norm, s=45, edgecolors="k", linewidths=0.7, zorder=8)
        for x, y, w in zip(xs, ys, names):
            ax.annotate(short(w), (x, y), xytext=(5, 4), textcoords="offset points", fontsize=6.5, zorder=9,
                        path_effects=[pe.withStroke(linewidth=1.8, foreground="white")])
        med = bb.T_field_median_degC.iloc[0] if len(bb) else float("nan")
        extra = f"; field median T = {med:.1f} °C" if kind == "T" else ""
        ax.set_title(f"[interpolated, not data] {bt:.0f} to {bt - 500:.0f} mRL{extra}" + ("" if ok else " (symbols only, < 3 wells)"), fontsize=7.5)
    for ax in axes[n:]:
        ax.axis("off")
    cax = fig.add_axes([0.25, 0.03, 0.5, 0.015])
    sm = plt.cm.ScalarMappable(norm=norm, cmap=C.dcmap)
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label(("T minus field median (°C); red hotter" if kind == "T" else "P minus feed-zone trend (bar); red above trend") +
                 ". Linear interpolation on a Delaunay triangulation of the wells' band values; colour only inside the convex hull; "
                 "hatched grey = no data (not constrained).", fontsize=7)
    fig.suptitle("[interpolated, not data] " + ("Temperature relative to the field median" if kind == "T" else "Pressure relative to the "
                 "feed-zone trend") + " in 500 m elevation bands (plan views at band mid-elevation; " + ("°C" if kind == "T" else "bar") + ")",
                 fontsize=11, y=0.995)
    fig.tight_layout(rect=(0, 0.06, 1, 0.97))
    path = os.path.join(ADDFIG, "maps_temperature_interpolated.png" if kind == "T" else "maps_pressure_interpolated.png")
    save_fig(fig, path)


# ----------------------------------------------------------------- conceptual model elements
def pressure_connections(C, line_pairs):
    """P difference between neighbouring wells at the deepest common candidate elevation."""
    cands = [-2000, -1500, -1000, -750, -500, -250]
    rows = []
    for a, b in line_pairs:
        ga = C.targets[C.targets.well == a]
        gb = C.targets[C.targets.well == b]
        for z in cands:
            pa = float(interp_down(z, ga.elevation_mRL.values, ga.P_bar.values))
            pb = float(interp_down(z, gb.elevation_mRL.values, gb.P_bar.values))
            if pa == pa and pb == pb:
                rows.append(dict(well_a=a, well_b=b, elevation_mRL=z, P_a_bar=round(pa, 1), P_b_bar=round(pb, 1),
                                 dP_a_minus_b_bar=round(pa - pb, 1),
                                 dev_a_from_trend_bar=round(pa - float(C.trend_fn(z)), 1),
                                 dev_b_from_trend_bar=round(pb - float(C.trend_fn(z)), 1),
                                 T_a_degC=round(float(interp_down(z, ga.elevation_mRL.values, ga.T_degC.values)), 1),
                                 T_b_degC=round(float(interp_down(z, gb.elevation_mRL.values, gb.T_degC.values)), 1)))
                break
    return pd.DataFrame(rows)


def conceptual_elements(C, A):
    """Numbered elements of the conceptual model, derived from the tables; returns (list, context)."""
    S = C.summary.set_index("well")
    W = [w["well"] for w in A["wells"]]
    hot2 = list(S.Tmax_degC.sort_values(ascending=False).index[:2])
    ups = [w for w in C.wells_order if S.loc[w, "regime"] == "upflow"]
    hot_c = np.mean([[C.info[w]["x"], C.info[w]["y"]] for w in (ups or hot2)] + [[C.info[w]["x"], C.info[w]["y"]] for w in hot2], axis=0)
    outf = [w for w in C.wells_order if S.loc[w, "regime"] == "outflow"]
    cold = [w for w in C.wells_order if S.loc[w, "regime"] == "recharge/downflow"]
    nb = [w for w in C.wells_order if S.loc[w, "near_boiling"]]
    els = []
    ctx = dict(hot_c=hot_c, outf=outf, cold=cold, nb=nb, ups=ups, hot2=hot2)
    # 1 upflow
    obs = "; ".join(f"{short(w)} Tmax {S.loc[w, 'Tmax_degC']:.1f} degC at {S.loc[w, 'Tmax_elev_mRL']:.0f} mRL" for w in hot2)
    if ups:
        obs += "; " + "; ".join(f"{short(w)} regime upflow (no inversion, shallow gradient {S.loc[w, 'shallow_gradient_degC_per_km']:.0f} degC/km)" for w in ups)
    els.append(dict(n=1, kind="upflow", name="Deep hot upflow", wells=hot2 + [w for w in ups if w not in hot2], solid=False, claim="interp",
                    obs=obs + ". [obs] The two hottest wells lie in the same corner of the field; the upflow itself is inferred from their temperature"
                    " and the absence of an inversion at " + ", ".join(short(w) for w in ups or hot2) + "."))
    # 2 outflow solid on line, 3 outflow elsewhere
    on_line = [w for w in outf if w in W]
    off_line = [w for w in outf if w not in W]
    if on_line:
        obs = "; ".join(f"{short(w)} inversion {S.loc[w, 'inversion_degC']:.1f} degC from {S.loc[w, 'inversion_top_mRL']:.0f} to {S.loc[w, 'inversion_base_mRL']:.0f} mRL"
                        for w in on_line)
        els.append(dict(n=2, kind="outflow_on", name="Lateral outflow along the section line", wells=on_line, solid=True, claim="mixed",
                        obs=obs + ". [obs] Hottest water at the top of the inversion, cooler water below it."))
    if off_line:
        obs = "; ".join(f"{short(w)} Tmax {S.loc[w, 'Tmax_degC']:.1f} degC at {S.loc[w, 'Tmax_elev_mRL']:.0f} mRL with inversion {S.loc[w, 'inversion_degC']:.1f} degC"
                        for w in off_line)
        els.append(dict(n=3, kind="outflow_off", name="Lateral outflow away from the section line (direction inferred)", wells=off_line, solid=False,
                        claim="interp", obs=obs + ". [obs] Each well off the line shows the same hot-over-cool profile; the flow direction is not measured."))
    # 4 cold recharge / downflow
    if cold:
        obs = "; ".join(f"{short(w)} Tmax {S.loc[w, 'Tmax_degC']:.1f} degC, T at TD {S.loc[w, 'T_at_TD_degC']:.1f} degC at {S.loc[w, 'TD_mRL']:.0f} mRL" for w in cold)
        els.append(dict(n=4, kind="cold", name="Cold recharge or downflow", wells=cold, solid=True, claim="mixed",
                        obs=obs + ". [obs] Whole surveyed column below 100 degC."))
    # 5 cap
    capw = [w for w in W if S.loc[w, "z_200C_mRL"] == S.loc[w, "z_200C_mRL"]]
    if capw:
        obs = "; ".join(f"{short(w)}: surface {C.info[w]['surface']:.0f} mRL, 200 degC at {S.loc[w, 'z_200C_mRL']:.0f} mRL, shallow gradient {S.loc[w, 'shallow_gradient_degC_per_km']:.0f} degC/km"
                        for w in capw)
        els.append(dict(n=5, kind="cap", name="Cap (conductive gradient, low permeability)", wells=capw, solid=True, claim="mixed",
                        obs=obs + ". [obs] Steep, near-linear temperature rise from the surface to the 200 degC level."))
    # 6 near boiling
    if nb:
        obs = "; ".join(f"{short(w)} minimum margin {S.loc[w, 'min_boiling_margin_degC']:.2f} degC at {S.loc[w, 'min_boiling_margin_elev_mRL']:.0f} mRL" for w in nb)
        els.append(dict(n=6, kind="nearboil", name="Near-boiling zone (not two-phase: margin above 1 degC)", wells=nb, solid=True, claim="obs",
                        obs=obs + ". [obs] Temperature follows boiling-point-for-depth to within 5 degC."))
    return els, ctx


def fig_conceptual(C, A, els, ctx, conns):
    S = C.summary.set_index("well")
    ln = A["ln"]
    SD = A
    fig = plt.figure(figsize=(17, 10.5))
    axp = fig.add_axes([0.04, 0.06, 0.34, 0.80])
    base_layer(C, axp, labels=True, small=True)
    label_wells(C, axp, fs=8)
    p0, p1 = ln["p_ext0"], ln["p_ext1"]
    axp.plot([p0[0], p1[0]], [p0[1], p1[1]], color="#7b3294", lw=1.4, ls="-.", zorder=6)
    for p, lab in ((p0, ln["id"]), (p1, ln["id"] + "′")):
        axp.text(p[0], p[1], lab, fontsize=10, fontweight="bold", color="#7b3294", ha="center", va="center", zorder=9,
                 bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="#7b3294", lw=0.8))
    hc = ctx["hot_c"]

    def num(ax, x, y, n, dx=0, dy=0, fs=9):
        ax.annotate(str(n), (x, y), xytext=(dx, dy), textcoords="offset points", fontsize=fs, fontweight="bold", color="white", ha="center",
                    va="center", zorder=12, bbox=dict(boxstyle="circle,pad=0.25", fc="k", ec="white", lw=0.8))

    # plan: pressure connections (thin dashed black)
    for _, r in conns.iterrows():
        a, b = C.info[r.well_a], C.info[r.well_b]
        axp.plot([a["x"], b["x"]], [a["y"], b["y"]], "k--", lw=0.8, zorder=5)
        axp.text(0.5 * (a["x"] + b["x"]), 0.5 * (a["y"] + b["y"]), f"{r.dP_a_minus_b_bar:+.1f} bar", fontsize=6.5, color="k", ha="center",
                 zorder=9, path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    for e in els:
        if e["kind"] == "upflow":
            axp.plot([hc[0]], [hc[1]], marker="^", color="red", ms=20, mec="k", zorder=10, ls="")
            num(axp, hc[0], hc[1], e["n"], 0, -1, 8)
        elif e["kind"] in ("outflow_on", "outflow_off"):
            for w in e["wells"]:
                x, y = C.info[w]["x"], C.info[w]["y"]
                vec = np.array([x, y]) - hc
                if np.hypot(*vec) < 200:
                    continue
                end = hc + 0.9 * vec
                axp.annotate("", xy=end, xytext=hc, arrowprops=dict(arrowstyle="-|>", color="#ff8c00", lw=2.6,
                             ls="-" if e["solid"] else "--", shrinkA=10, shrinkB=0, mutation_scale=16), zorder=8)
            mid = hc + 0.5 * (np.array([C.info[e["wells"][0]]["x"], C.info[e["wells"][0]]["y"]]) - hc)
            num(axp, mid[0], mid[1], e["n"], 0, 9)
        elif e["kind"] == "cold":
            for w in e["wells"]:
                x, y = C.info[w]["x"], C.info[w]["y"]
                vec = hc - np.array([x, y])
                u = vec / np.hypot(*vec)
                st = np.array([x, y]) - 0.0 * u
                axp.annotate("", xy=st + 900 * u, xytext=st - 100 * u, arrowprops=dict(arrowstyle="-|>", color="#1f77b4", lw=2.6, ls="--",
                             shrinkA=0, shrinkB=0, mutation_scale=16), zorder=8)
                num(axp, x, y, e["n"], -12, -10)
        elif e["kind"] == "nearboil":
            for w in e["wells"]:
                x, y = C.info[w]["x"], C.info[w]["y"]
                axp.add_patch(patches.Circle((x, y), 260, facecolor="none", edgecolor="k", hatch="....", lw=0.6, zorder=7))
                num(axp, x, y, e["n"], 14, -14)
    axp.set_title("Plan view (base layer of F1)", fontsize=10)
    # section panel
    axs = fig.add_axes([0.46, 0.38, 0.52, 0.48])
    meta = draw_section_into(C, SD, axs, "T", (0.52 * 17, 0.48 * 10.5))
    W = SD["wells"]
    spos = {w["well"]: float(np.mean(w["s"])) for w in W}
    # pressure connections along the line
    for _, r in conns.iterrows():
        if r.well_a in spos and r.well_b in spos:
            sa, sb = spos[r.well_a], spos[r.well_b]
            axs.plot([sa, sb], [r.elevation_mRL] * 2, "k--", lw=0.9, zorder=9)
            axs.text(0.5 * (sa + sb), r.elevation_mRL + (25 if _ % 2 == 0 else -70), f"ΔP {r.dP_a_minus_b_bar:+.1f} bar at {r.elevation_mRL:.0f} mRL", fontsize=7, ha="center", zorder=11,
                     path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    hs = float(np.mean([spos[w] for w in ctx["hot2"] if w in spos])) if any(w in spos for w in ctx["hot2"]) else None
    for e in els:
        if e["kind"] == "upflow" and hs is not None:
            axs.annotate("", xy=(hs, -150), xytext=(hs, -1100), arrowprops=dict(arrowstyle="-|>", color="red", lw=4, ls="--", mutation_scale=22), zorder=10)
            num(axs, hs, -650, e["n"], 12, 0)
        elif e["kind"] in ("outflow_on",) and hs is not None:
            for w in e["wells"]:
                if w not in spos:
                    continue
                z = float(S.loc[w, "Tmax_elev_mRL"])
                if abs(spos[w] - hs) < 300:
                    continue
                axs.annotate("", xy=(spos[w], z), xytext=(hs, z), arrowprops=dict(arrowstyle="-|>", color="#ff8c00", lw=3.2, mutation_scale=18), zorder=10)
            w0 = [w for w in e["wells"] if w in spos and abs(spos[w] - hs) >= 300]
            if w0:
                zz = float(S.loc[w0[0], "Tmax_elev_mRL"])
                num(axs, 0.5 * (hs + spos[w0[0]]), zz, e["n"], 0, 12)
        elif e["kind"] == "cold":
            for w in e["wells"]:
                if w in spos:
                    zt = float(S.loc[w, "TD_mRL"])
                    axs.annotate("", xy=(spos[w] + 80, zt + 200), xytext=(spos[w] + 80, zt + 1100), arrowprops=dict(arrowstyle="-|>", color="#1f77b4", lw=4, mutation_scale=20), zorder=10)
                    num(axs, spos[w] + 80, zt + 650, e["n"], 12, 0)
        elif e["kind"] == "cap":
            ws = [w for w in e["wells"] if w in spos]
            if len(ws) >= 2:
                xs_ = [spos[w] for w in ws]
                top = [C.info[w]["surface"] for w in ws]
                base = [float(S.loc[w, "z_200C_mRL"]) for w in ws]
                axs.fill_between(xs_, top, base, color="0.35", alpha=0.45, zorder=3, label="cap")
                axs.plot(xs_, base, color="0.15", lw=1.2, ls="--", zorder=8)
                for xx, tt, bb in zip(xs_, top, base):
                    axs.text(xx, tt - 5, f"top {tt:.0f}", fontsize=7, ha="center", va="top", zorder=11, path_effects=HALO)
                    axs.text(xx, bb + 5, f"base {bb:.0f}", fontsize=7, ha="center", va="bottom", zorder=11, path_effects=HALO)
                num(axs, 0.5 * (xs_[0] + xs_[-1]), 0.5 * (np.mean(top) + np.mean(base)), e["n"], 0, 0)
        elif e["kind"] == "nearboil":
            for w in e["wells"]:
                if w in spos:
                    z0 = float(S.loc[w, "min_boiling_margin_elev_mRL"])
                    num(axs, spos[w] + 110, z0, e["n"], 8, 0)
    # solid / dashed key and symbology legend (always shown)
    hnd = [Line2D([], [], color="red", lw=4, marker=">", ms=0, label="Deep hot upflow (red arrow)"),
           Line2D([], [], color="#ff8c00", lw=3, label="Lateral outflow (orange arrow at its depth)"),
           Line2D([], [], color="#1f77b4", lw=3, label="Cold recharge / downflow (blue arrow)"),
           patches.Patch(facecolor="0.35", alpha=0.45, label="Cap (grey band; top and base labelled)"),
           patches.Patch(facecolor="none", edgecolor="k", hatch="....", label="Two-phase / near-boiling zone (stipple)"),
           Line2D([], [], color="k", lw=1.8, label="100 / 200 °C isotherms (as F4)"),
           Line2D([], [], color="k", lw=0.9, ls="--", label="Pressure connection (dashed, ΔP labelled)"),
           Line2D([], [], color="0.2", lw=2.5, ls="-", label="solid arrow: supported by a measured feature"),
           Line2D([], [], color="0.2", lw=2.5, ls="--", label="dashed arrow: inferred"),
           patches.Patch(facecolor="#ececec", edgecolor="#9a9a9a", hatch="////", label="no data")]
    fig.legend(handles=hnd, loc="lower left", bbox_to_anchor=(0.44, 0.04), ncol=2, fontsize=8.5, frameon=False, title="Fixed symbology (F6)")
    # element key
    txt = "Elements (numbers key to the conceptual-model table, report section 5):\n" + "\n".join(
        f"{e['n']}  {e['name']}  [{e['claim']}; {'solid' if e['solid'] else 'dashed'}]" for e in els)
    fig.text(0.46, 0.335, txt, fontsize=8, va="top", family="monospace")
    fig.suptitle(f"F6  Conceptual model (data only): plan view and section along line {ln['id']} "
                 f"({p0[0]:.0f} E, {p0[1]:.0f} N to {p1[0]:.0f} E, {p1[1]:.0f} N); elevation in mRL, temperature in °C; "
                 f"section vertical exaggeration ×{meta['ve']:.1f}", fontsize=11, y=0.985)
    axs.set_title("Section (F4 background: colour-filled data envelope, isotherms, no-data hatching)", fontsize=9, pad=70)
    path = os.path.join(CALFIG, "conceptual_model.png")
    fig.savefig(path, dpi=DPI, facecolor="white")
    plt.close(fig)
    add_index("F6", path, "Conceptual model: plan view and section A with numbered upflow, outflow, recharge, cap, near-boiling zone and pressure connections",
              ";".join(w["well"] for w in SD["wells"]), f"line {ln['id']}", "all tables; F1 and F4 backgrounds", "mixed")


# =========================================================================================
# 8. Report
# =========================================================================================
def fv(x, d=1):
    """Format a number for the report; blank for missing values."""
    try:
        if x is None or (isinstance(x, float) and math.isnan(x)) or pd.isna(x):
            return ""
    except (TypeError, ValueError):
        pass
    return f"{float(x):.{d}f}"


def md_table(df, cols=None, fmt=None, header=None):
    cols = cols or list(df.columns)
    fmt = fmt or {}
    head = header or cols
    out = ["| " + " | ".join(str(h) for h in head) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            d = fmt.get(c)
            if d is not None and not isinstance(v, str):
                cells.append(fv(v, d))
            else:
                cells.append("" if (not isinstance(v, str) and pd.isna(v)) else str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


class Quoted:
    """Collects every derived number quoted in the report text so that it is traceable to an output table."""

    def __init__(self):
        self.rows = []

    def __call__(self, name, value, d=1):
        s = fv(value, d)
        self.rows.append(dict(quantity=name, value=s))
        return s


def build_report(C, lines, SDs, els, ctx, conns, Q):
    S = C.summary.set_index("well")
    PC = C.pc.set_index("well")
    FZ = C.fz_out
    T = C.targets
    G = C.geo.set_index("well")
    W = C.wells_order
    sh = short
    q = Q
    hot = S.Tmax_degC.sort_values(ascending=False)
    h1, h2 = hot.index[0], hot.index[1]
    A, B = SDs[0], SDs[1]
    LA, LB = lines[0], lines[1]
    nconf = pd.Series({w: C.verdict[w]["confidence"] for w in W})
    cnt = nconf.value_counts()
    out = []
    ad = out.append
    ad("# Natural-state data interpretation for model calibration")
    ad("")
    ad(f"Generated by `Scripts/natural_state_analysis.py` from `Baseline/` only. Run date: {dt.date.today().isoformat()}")
    ad("")
    ad("Every statement is tagged with exactly one of [obs] (a measured value or a value computed directly from the data), [interp] "
       "(an interpretation), [assumption] (a stated assumption) or [reference model] (a property of the starting model). "
       "Elevations are in mRL, pressures in bar, temperatures in °C, enthalpy in kJ/kg, flow in kg/s.")
    ad("")
    # ---------------------------------------------------------------- 1
    ad("## 1. Key conclusions for calibration")
    ad("")
    gv = [w for w in W if C.verdict[w]["verdict"] == "natural state"]
    n_disturbed = [w for w in W if C.verdict[w]["verdict"] == "possibly disturbed"]
    ad(f"1. **All {len(W)} surveys carry the date {C.info[W[0]]['survey_date'].isoformat()}, but production and injection began on "
       f"{C.first_prod_field.isoformat()}.** The gauge test (not the date) decides: {', '.join(sh(w) for w in gv)} match their first (pre-production) gauge "
       f"reading to within {q('gauge precision', PREC_P, 1)} bar, so they represent the natural state [obs]. "
       f"The verdicts of the other wells are: {'; '.join(f'{sh(w)} {C.verdict[w][chr(118)+chr(101)+chr(114)+chr(100)+chr(105)+chr(99)+chr(116)]} ({C.verdict[w][chr(99)+chr(111)+chr(110)+chr(102)+chr(105)+chr(100)+chr(101)+chr(110)+chr(99)+chr(101)]})' for w in W if w not in gv)} [interp]. "
       f"The date label conflicts with the gauge test and is a data question, not a reason to downgrade [obs].")
    ad("   **Action:** use the surveys of all wells as natural-state targets, weighted by the confidence factor (section 6); ask the data owner which date applies.")
    ad("")
    ad(f"2. **Hot upflow lies at the {sh(h1)}–{sh(h2)} corner of the field.** {sh(h1)} reaches {q('Tmax h1', S.loc[h1, 'Tmax_degC'])} °C at {q('Tmax elev h1', S.loc[h1, 'Tmax_elev_mRL'], 0)} mRL "
       f"(regime {S.loc[h1, 'regime']}) and {sh(h2)} reaches {q('Tmax h2', S.loc[h2, 'Tmax_degC'])} °C at {q('Tmax elev h2', S.loc[h2, 'Tmax_elev_mRL'], 0)} mRL [obs]. "
       f"Both lie within {q('dist W2 fault3', G.loc[h2, 'dist_nearest_fault_m'], 0)} m ({sh(h2)}) and {q('dist W10 fault3', G.loc[h1, 'dist_nearest_fault_m'], 0)} m ({sh(h1)}) of {G.loc[h2, 'nearest_fault']} [obs]. "
       f"A deep hot upflow here, possibly along that fault, is the simplest explanation [interp].")
    ad("   **Action:** put the deep heat and mass input beneath these two wells; check that the model reaches about 270 °C at -250 to -550 mRL there.")
    ad("")
    outf = ctx["outf"]
    ad(f"3. **A lateral outflow layer lies at about -450 to -750 mRL.** {', '.join(sh(w) for w in outf)} have a temperature inversion of at least 10 °C: the hottest "
       f"water sits at {q('outflow min', min(S.loc[w, 'Tmax_elev_mRL'] for w in outf), 0)} to {q('outflow max', max(S.loc[w, 'Tmax_elev_mRL'] for w in outf), 0)} mRL and cooler water lies below it "
       f"(largest inversion {q('max inv', S.loc[outf, 'inversion_degC'].max())} °C at {sh(S.loc[outf, 'inversion_degC'].idxmax())}) [obs]. "
       f"The outflow spreads away from the upflow to the south, east and north-east [interp].")
    ad("   **Action:** make the layer around -450 to -750 mRL permeable laterally and weight the temperature points in it at their table values.")
    ad("")
    capw = [w for w in W if S.loc[w, "z_200C_mRL"] == S.loc[w, "z_200C_mRL"]]
    g_cap = S.loc[capw, "shallow_gradient_degC_per_km"]
    ad(f"4. **A conductive cap overlies the hot zone.** Shallow gradients are {q('gap min', g_cap.min(), 0)} to {q('gap max', g_cap.max(), 0)} °C/km in the wells that reach 200 °C, "
       f"and the 200 °C level lies at {q('z200 max', S.loc[capw, 'z_200C_mRL'].max(), 0)} to {q('z200 min', S.loc[capw, 'z_200C_mRL'].min(), 0)} mRL [obs]. "
       f"With an assumed conductivity of 2 W/m/K the gradients imply a conductive flux of about {q('flux min', 2 * g_cap.min() / 1000, 1)} to {q('flux max', 2 * g_cap.max() / 1000, 1)} W/m2 [assumption]. "
       f"The cap is a low-permeability layer above about -300 to 0 mRL near the upflow [interp].")
    ad("   **Action:** give the cap a low vertical permeability and compare the modelled shallow gradient with the table values; do not copy the flux as a boundary condition.")
    ad("")
    cold = ctx["cold"]
    ad(f"5. **Cold water occupies the north-east and south-west margins.** {', '.join(sh(w) for w in cold)} stay below 100 °C over their whole survey "
       f"(Tmax {', '.join(q('Tmax ' + w, S.loc[w, 'Tmax_degC']) for w in cold)} °C) [obs]. "
       f"{sh('WELL-4')} lies {q('W4 fault', G.loc['WELL-4', 'dist_Fault_3_m'], 0)} m from Fault 3 and {sh('WELL-8')} {q('W8 river', G.loc['WELL-8', 'dist_Tarawera_River_m'], 0)} m from the Tarawera River [obs]. "
       f"They are possibly cold recharge or downflow zones associated with the fault and the river respectively [interp].")
    ad("   **Action:** allow cold water to enter at these margins, and test whether the cold downflow reaches the full surveyed depth.")
    ad("")
    ad(f"6. **Formation pressure is defined only at feed zones, and deep pressures agree across the field.** The trend through the gauge-verified feed zones is "
       f"P = {q('trend a', C.trend[0], 2)} {q('trend b', C.trend[1], 4)} z; the tolerance is {q('tolerance', C.tol, 1)} bar (the larger of 3 bar and twice the {q('rms', C.trend_rms, 2)} bar RMS residual) [obs]. "
       f"All feed-zone pressures except {', '.join(sh(w) for w in n_disturbed) or 'none'} lie within it [obs]. "
       f"The wellbore column in cased intervals is steeper than a column of water in the hot wells (ratio up to {q('max cased ratio', C.hydro[C.hydro.band == 'cased'].ratio_measured_to_IF97.max(), 2)}), so shallow wellbore pressure is not reservoir pressure [obs].")
    ad("   **Action:** use `weight_P` at feed zones only; ignore shallow wellbore pressure; treat the cold-well hydrostatic ratio as a datum or convention question (section 4).")
    ad("")
    _g10 = T[(T.well == "WELL-10") & (T.boiling_margin_degC < 5)]
    nb_span = float(_g10.elevation_mRL.max() - _g10.elevation_mRL.min())
    ad(f"7. **{sh('WELL-10')} is near boiling at shallow depth.** Its minimum boiling margin is {q('W10 minbm', S.loc['WELL-10', 'min_boiling_margin_degC'], 2)} °C at {q('W10 minbm z', S.loc['WELL-10', 'min_boiling_margin_elev_mRL'], 0)} mRL, "
       f"above the 1 °C two-phase threshold, so its regime is {S.loc['WELL-10', 'regime']} [obs]. The temperature follows boiling-point-for-depth over about {q('nb span', nb_span, 0)} m (margin below 5 °C) [interp]. "
       f"Its gauge starts after its first flow, so it has no pre-production gauge test [obs].")
    ad("   **Action:** the model must reproduce boiling-point-for-depth here; keep the model from boiling more than the data allow.")
    ad("")
    ad(f"8. **Confidence differs between wells, so the factor changes relative weights.** Confidence counts: " +
       ", ".join(f"{k} {int(v)}" for k, v in cnt.items()) + f". Mean temperature weight ranges from {q('wmin', T.groupby('well').weight.mean().min(), 2)} to {q('wmax', T.groupby('well').weight.mean().max(), 2)} [obs].")
    ad("   **Action:** use the `weight` and `weight_P` columns as given; revisit them if the gauge questions in section 7 are resolved.")
    ad("")
    # ---------------------------------------------------------------- 2
    ad("## 2. The data")
    ad("")
    ad("### 2.1 What exists")
    ad("")
    ad("Roles assigned to every file under `Baseline/` (full table: `../Additional_Output/input_inventory.csv`):")
    ad("")
    inv = C.inv.copy()
    inv["path"] = inv.path.str.replace("Baseline/Teaching_Model_3D/", "", regex=False)
    ad(md_table(inv, ["role", "path"], header=["role", "path (below Baseline/Teaching_Model_3D)"]))
    ad("")
    ad("[obs] The reference model files (rock types, boundary conditions, natural-state results, set-up) and production, injection and forecast files were not used as evidence. "
       "[assumption] Production and injection records inside the well files were used only to test whether the surveys represent the natural state.")
    ad("")
    X = C.crosscheck
    nok = int((X.status == "OK").sum())
    bad = X[X.status.isin(["MISMATCH", "OUTSIDE", "INCONSISTENT", "CHECK"])]
    ad(f"### 2.2 How the data were checked")
    ad("")
    ad(f"[obs] {q('checks ok', nok, 0)} of {q('checks all', len(X), 0)} cross-checks passed (`../Additional_Output/data_crosscheck.csv`). The plain-column files reproduce every JSON survey point and transient value exactly; "
       f"the time axis of the transient files is years since {C.dat_ref.isoformat() if C.dat_ref else 'n/a'} [obs]. Surveys lie within the well tracks, and every feed zone lies within a slotted interval [obs]. "
       f"The model block of each survey point agrees with the well-to-block file, and a point-in-polygon search reproduces the same columns for all wells [obs].")
    ad("")
    if len(bad):
        ad("Items that did not pass cleanly:")
        ad("")
        for _, r in bad.iterrows():
            note = str(r.note) if isinstance(r.note, str) else ""
            ad(f"- {r.well}: {r['item']} ({r.status}). {note[:170]} [obs]")
        ad("")
    ad("Irregular file names: " + "; ".join(b for b, _ in C.dat_irregular[:6]) + " (prefix `Wel` instead of `Well`; matched to WELL-10 by content) [obs].")
    ad("")
    ad("### 2.3 Provenance labels")
    ad("")
    ad(f"[obs] The casing and liner design carries the source label '{C.wells[W[0]]['design_source']}' for every well. The survey records have an empty Source and Comment field and StableFlag true. "
       f"Feed-zone Source fields hold the unit label 'Text'. Quality tags found: {', '.join(C.qfound)}.")
    ad("")
    ad("Feed-zone quality mapping to the `weight_P` factor: " + ", ".join(f"{k} = {v:g}" for k, v in sorted(C.qmap.items(), key=lambda kv: -kv[1])) +
       ". 'Synthetic' is not a quality grade and takes the lowest factor [assumption].")
    ad("")
    ad("### 2.4 Is each survey the natural state?")
    ad("")
    ad(f"[obs] First production in the field: {C.first_prod_field.isoformat()}. First injection: {C.first_inj_field.isoformat()}. Tolerance for feed-zone pressure against the verified trend: {q('tol again', C.tol, 1)} bar. "
       f"Injectate temperature from the enthalpy of {q('inj h', list(C.inj.values())[0]['h'], 1)} kJ/kg (IF97 Region 1 at the feed-zone pressure): "
       + ", ".join(f"{sh(w)} {q('Tinj ' + w, v['T'], 1)} °C" for w, v in C.inj.items()) + ".")
    ad("")
    t = C.pc[["well", "gauge_elev_mRL", "survey_P_at_gauge_bar", "first_gauge_date", "first_gauge_P_bar", "diff_survey_minus_first_bar",
              "own_first_flow_date", "verdict", "confidence"]].copy()
    ad(md_table(t, fmt={"gauge_elev_mRL": 0, "survey_P_at_gauge_bar": 2, "first_gauge_P_bar": 1, "diff_survey_minus_first_bar": 2},
                header=["well", "gauge mRL", "survey P at gauge", "first gauge date", "first gauge P", "survey minus first", "own first flow", "verdict", "confidence"]))
    ad("")
    ad("Reason for each verdict:")
    ad("")
    for w in W:
        ad(f"- {w}: {C.verdict[w]['verdict']}, {C.verdict[w]['confidence']}. {C.verdict[w]['reason']}. [obs]")
    ad("")
    ad(f"[interp] {sh('WELL-2')} and {sh('WELL-3')} match the January 2014 gauge reading, not the reading of 2015 (W2 had already fallen by about {q('drawdown W2', PC.loc['WELL-2', 'first_gauge_P_bar'] - C.gauge_series['WELL-2']['values'][1], 0)} bar). "
       f"So the surveys represent the pre-production state although they are labelled July 2015. "
       f"{sh('WELL-1')} is flagged possibly disturbed only because its feed zone lies {q('W1 dev', abs(FZ[FZ.well == 'WELL-1'].P_minus_trend_bar).max(), 1)} bar from the trend, "
       f"just above the {q('tol W1', C.tol, 1)} bar tolerance; it has no gauge, no production and no injection, so the miss is more likely a lateral pressure difference or its datum (section 7).")
    ad("")
    ad("![F7: transient gauge pressure at the gauge elevation against the survey pressure interpolated to it, 2014-2025. Claim type: obs.](figures/gauge_check.png)")
    ad("")
    ad("![F1: base map with faults, water bodies, wells and section lines A and B (plan view). Claim type: obs.](figures/well_map.png)")
    ad("")
    ad(f"[obs] Model surface and wellhead: the wellhead lies below the model surface of its column in {int((C.geo.wellhead_minus_surface_m < 0).sum())} wells, by up to "
       f"{q('max below', -C.geo.wellhead_minus_surface_m.min(), 0)} m ({sh(C.geo.set_index('well').wellhead_minus_surface_m.idxmin())}); "
       f"survey points above the model surface: {int(T.above_model_surface.sum())} (`../Additional_Output/well_geometry_check.csv`).")
    ad("")
    # ---------------------------------------------------------------- 3
    ad("## 3. Temperature structure")
    ad("")
    cols = ["well", "Tmax_degC", "Tmax_elev_mRL", "T_at_TD_degC", "TD_mRL", "inversion_degC", "inversion_top_mRL", "inversion_base_mRL",
            "shallow_gradient_degC_per_km", "deep_gradient_degC_per_km", "z_100C_mRL", "z_200C_mRL", "min_boiling_margin_degC", "regime"]
    ad(md_table(C.summary, cols, fmt={c: 1 for c in cols if c not in ("well", "regime")} | {"Tmax_elev_mRL": 0, "TD_mRL": 0, "inversion_top_mRL": 0,
                "inversion_base_mRL": 0, "shallow_gradient_degC_per_km": 0, "deep_gradient_degC_per_km": 0, "z_100C_mRL": 0, "z_200C_mRL": 0},
                header=["well", "Tmax", "at mRL", "T at TD", "TD mRL", "inversion", "top mRL", "base mRL", "shallow grad °C/km", "deep grad °C/km",
                        "100 °C mRL", "200 °C mRL", "min margin", "regime"]))
    ad("")
    ad("[obs] The table is computed from `NS_calibration_targets.csv` (definitions of section 4.4). Blank deep gradients mean fewer than three points or a span under 300 m below the deeper of the inversion base and the deepest feed zone.")
    ad("")
    ad(f"[obs] Isothermal sections of at least 300 m in open hole occur in: "
       + ", ".join(f"{sh(w)} ({int(T[(T.well == w) & T.isothermal_open].shape[0])} points)" for w in W if T[(T.well == w) & T.isothermal_open].shape[0]) + ".")
    ad(f"[interp] Open-hole isothermal sections in {sh('WELL-5')} and {sh('WELL-8')} lie across or below the feed zones, so internal flow in the wellbore may have flattened them; they carry the reduced weight of 0.5.")
    ad("")
    ad(f"[interp] The hot wells are the ones closest to Fault 3 ({sh('WELL-2')} {q('d2', G.loc['WELL-2', 'dist_Fault_3_m'], 0)} m, {sh('WELL-10')} {q('d10', G.loc['WELL-10', 'dist_Fault_3_m'], 0)} m). "
       f"{sh('WELL-6')} ({q('d6', G.loc['WELL-6', 'dist_Fault_1_m'], 0)} m from Fault 1) is hot but cooler. The data cannot decide whether the faults control the upflow; a deeper well beside {sh('WELL-10')} (it stops at {q('W10 TD', S.loc['WELL-10', 'TD_mRL'], 0)} mRL in a {q('W10 TD T', S.loc['WELL-10', 'T_at_TD_degC'])} °C isothermal section) would decide it.")
    ad("")
    ad("![F2: temperature profiles of all wells on common axes (-2900 to 400 mRL). Claim type: obs.](figures/profile_panels_Synthetic_All.png)")
    ad("")
    ad("Individual profiles with feed zones, shoe, pressure and verdict: " + ", ".join(f"[{sh(w)}](figures/profile_{w}.png)" for w in W) + " (F2, claim type obs).")
    ad("")
    ad("![F3: temperature relative to the field median in 500 m elevation bands, plan view at band mid-elevation. Claim type: obs.](figures/maps_temperature_bands.png)")
    ad("")
    for SD in SDs:
        ln = SD["ln"]
        ends = f"{ln['id']} ({LA['p_ext0'][0] if ln is LA else LB['p_ext0'][0]:.0f} E, {LA['p_ext0'][1] if ln is LA else LB['p_ext0'][1]:.0f} N) to {ln['id']}′ ({ln['p_ext1'][0]:.0f} E, {ln['p_ext1'][1]:.0f} N)"
        offs = ", ".join(f"{sh(w['well'])} ({w['offset']:+.0f} m)" for w in SD["wells"])
        ad(f"[obs] Section line {ends}, {ln['desc']}; included wells within {q('maxoff' + ln['id'], ln['max_offset'], 0)} m with offsets: {offs}.")
        ad("")
        ad(f"![F4: temperature section along line {ln['id']} ({ends}); data envelope colour-filled, isotherms every 50 °C, no-data hatched. Claim type: obs.](figures/section_{ln['id']}.png)")
        ad("")
    ad("![Interpolated companion: temperature on line A, [interp] linear interpolation, not data.](../Additional_Output/figures/section_A_interpolated.png)")
    ad("")
    ad("![Interpolated companion: temperature on line B, [interp] linear interpolation, not data.](../Additional_Output/figures/section_B_interpolated.png)")
    ad("")
    ad("![Interpolated companion: plan views of temperature relative to the field median, [interp] not data.](../Additional_Output/figures/maps_temperature_interpolated.png)")
    ad("")
    ad("[interp] The interpolated figures show colours between wells for orientation only; no value is quoted from them.")
    ad("")
    # ---------------------------------------------------------------- 4
    ad("## 4. Pressure structure")
    ad("")
    ad("Temperature and pressure at common elevations (interpolated within each survey; blank where the elevation lies outside it). "
       "Pressure deviation from the feed-zone trend in brackets.")
    ad("")
    els_c = sorted(C.common.elevation_mRL.unique(), reverse=True)
    rows = []
    for w in W:
        g = C.common[C.common.well == w].set_index("elevation_mRL")
        r = {"well": w}
        for e in els_c:
            if e in (0, -500, -1000, -1500, -2000):
                t_, p_, d_ = g.loc[e, "T_degC"], g.loc[e, "P_bar"], g.loc[e, "P_minus_trend_bar"]
                r[f"{e:.0f} mRL"] = "" if pd.isna(t_) else f"{t_:.1f} °C / {p_:.1f} bar ({d_:+.1f})"
        rows.append(r)
    ad(md_table(pd.DataFrame(rows)))
    ad("")
    ad("[obs] Values are from `../Additional_Output/T_P_at_common_elevations.csv`.")
    ad("")
    hs = C.hydro_summary
    ad("### Hydrostatic checks (single-phase intervals, IF97 Region 1, g = 9.80665 m/s2)")
    ad("")
    ad(md_table(hs, ["group", "band", "n_windows", "ratio_median", "ratio_min", "ratio_max"], fmt={"ratio_median": 3, "ratio_min": 3, "ratio_max": 3},
                header=["group", "band", "windows", "median ratio", "min", "max"]))
    ad("")
    ad(f"[obs] Ratio = measured pressure difference / pressure difference integrated through IF97 liquid density over the measured temperature. The rounding uncertainty of a window is below {q('max ratio unc', C.hydro.ratio_uncertainty.max(), 3)}. "
       f"Cold-well (Tmax < 100 °C: {', '.join(sh(w) for w in C.cold_wells)}) median ratio: {q('cold cased', C.cold_comp['cased'][0], 3)} (cased), {q('cold open', C.cold_comp['open'][0], 3)} (open). "
       f"This is the field-wide component; it is a possible datum or convention issue, not physical overpressure [interp].")
    ad("")
    wb = C.hydro_wb.copy()
    ad(md_table(wb, ["well", "band", "n_windows", "length_m", "ratio_median", "cold_component", "excess_over_cold_component"],
                fmt={"length_m": 0, "ratio_median": 3, "cold_component": 3, "excess_over_cold_component": 3},
                header=["well", "band", "windows", "length m", "median ratio", "cold component", "excess above component"]))
    ad("")
    hc = C.hydro[C.hydro.band == "cased"]
    ad(f"[obs] In the cased band the hot wells have implied densities up to {q('max rho', hc.rho_implied_kgm3.max(), 0)} kg/m3, heavier than any liquid water, with the largest excess at {sh(C.hydro_wb[C.hydro_wb.band == 'cased'].sort_values('excess_over_cold_component').iloc[-1].well)}. "
       f"A datum offset cannot change a gradient. A gravity or unit convention would act at every depth, but the open-hole ratios are near 1 (median {q('open med', hs[(hs.group == 'all wells') & (hs.band == 'open')].ratio_median.iloc[0], 3)}). "
       f"Two-phase fluid would make the column lighter, not heavier, and rounding is below {q('rounding', C.hydro.ratio_uncertainty.max(), 3)}.")
    ad(f"[interp] What is left is that the shallow pressure of the hot wells (steep temperature gradient, cap) is not a static wellbore column in equilibrium with the survey temperature; it cannot be used as formation pressure. It is not used for `weight_P`.")
    ad("")
    ad("![F3: pressure relative to the hydrostatic trend of gauge-verified feed zones in 500 m bands. Claim type: obs.](figures/maps_pressure_bands.png)")
    ad("")
    for SD in SDs:
        ln = SD["ln"]
        ad(f"![F4: pressure deviation from the feed-zone trend along line {ln['id']} (1 bar steps; same wells and offsets as the temperature section). Claim type: obs.](figures/section_{ln['id']}_P.png)")
        ad("")
    ad("![Interpolated companion: pressure deviation on line A, [interp] not data.](../Additional_Output/figures/section_A_P_interpolated.png)")
    ad("")
    ad("![Interpolated companion: pressure deviation on line B, [interp] not data.](../Additional_Output/figures/section_B_P_interpolated.png)")
    ad("")
    ad("![Interpolated companion: pressure deviation in plan views, [interp] not data.](../Additional_Output/figures/maps_pressure_interpolated.png)")
    ad("")
    ad("![F5: pressure against temperature of all survey points with the IF97 saturation curve to the critical point. Claim type: obs.](figures/PT_state.png)")
    ad("")
    ad("Feed-zone pressure targets (survey values at the feed-zone mid elevation):")
    ad("")
    ad(md_table(FZ, ["well", "feedzone", "mid_mRL", "quality", "T_degC", "P_bar", "boiling_margin_degC", "block", "P_minus_trend_bar"],
                fmt={"mid_mRL": 1, "T_degC": 1, "P_bar": 1, "boiling_margin_degC": 1, "P_minus_trend_bar": 1}))
    ad("")
    ad(f"[obs] The deep pressure of the field is connected: at the deepest common elevation the neighbouring wells on line A differ by "
       + ", ".join(f"{sh(r.well_a)}–{sh(r.well_b)} {q('dP' + r.well_a + r.well_b, r.dP_a_minus_b_bar)} bar at {r.elevation_mRL:.0f} mRL" for _, r in conns.iterrows()) + " (`../Additional_Output/pressure_connections.csv`).")
    ad("[interp] Differences of a few bar between wells of unlike temperature are consistent with different water-column densities; a continuous permeable connection at depth is the simplest explanation.")
    ad("")
    # ---------------------------------------------------------------- 5
    ad("## 5. Conceptual model")
    ad("")
    ad(f"![F6: conceptual model, plan view and section along line A, elements numbered as in the table below. Solid arrows are supported by a measured feature, dashed arrows are inferred. Claim type: mixed.](figures/conceptual_model.png)")
    ad("")
    ad("| No. | Element | Supporting observation | Claim type | Arrow |")
    ad("|---|---|---|---|---|")
    for e in els:
        ad(f"| {e['n']} | {e['name']} | {e['obs']} | {e['claim']} | {'solid' if e['solid'] else 'dashed'} |")
    ad("")
    ad(f"| {len(els) + 1} | Pressure connection between wells | " + "; ".join(f"{sh(r.well_a)}–{sh(r.well_b)} {r.dP_a_minus_b_bar:+.1f} bar at {r.elevation_mRL:.0f} mRL" for _, r in conns.iterrows()) + " | obs | dashed line |")
    ad("")
    ad("Relation to faults and surface water, and what would decide between explanations:")
    ad("")
    ad(f"- The hottest wells ({sh(h1)}, {sh(h2)}) are nearest Fault 3 and the coldest ({sh('WELL-4')}) is also close to Fault 3 ({q('W4 f3', G.loc['WELL-4', 'dist_Fault_3_m'], 0)} m) [obs]. "
       f"Fault 3 may carry hot upflow in the south-west of its trace and cold recharge in the north-east, or the fault may have no role [interp].")
    ad(f"- Decisive observation: a temperature survey of a well drilled on Fault 3 between {sh('WELL-2')} and {sh('WELL-4')}, or a deeper {sh('WELL-10')} [interp].")
    ad(f"- {sh('WELL-1')} ({q('W1 river', G.loc['WELL-1', 'dist_river_m'], 0)} m from the river) is conductive with {q('W1 grad', S.loc['WELL-1', 'shallow_gradient_degC_per_km'], 0)} °C/km and a wellhead {q('W1 wh', -G.loc['WELL-1', 'wellhead_minus_surface_m'], 0)} m below the model surface [obs]. "
       f"It shows no outflow, so the outflow does not reach the east of the field at the depth of its survey [interp].")
    ad(f"- Pressure connection and the deep cooler water below the outflow (inversion bases at {q('inv base min', S.inversion_base_mRL.min(), 0)} to {q('inv base max', S.inversion_base_mRL.max(), 0)} mRL) can be explained by deep cold downflow or by a deep lateral cool flow; one deep well in the north-west would decide it [interp].")
    ad("")
    # ---------------------------------------------------------------- 6
    ad("## 6. Recommended calibration set-up")
    ad("")
    ad(f"**Targets.** `NS_calibration_targets.csv` holds {q('n targets', len(T), 0)} survey points in {len(W)} wells; `feedzones.csv` holds {q('n fz', len(FZ), 0)} feed zones; `T_P_at_layer_centres.csv` holds {q('n lc', len(C.layer_centres), 0)} layer-centre values without extrapolation [obs].")
    ad("")
    ws = C.weights_summary
    ad(md_table(ws, ["category", "n_points", "mean_weight", "max_weight"], fmt={"n_points": 0, "mean_weight": 2, "max_weight": 2},
                header=["location category", "points", "mean weight", "max weight"]))
    ad("")
    ad("[assumption] Weights are the defaults of the specification (location factor times confidence factor; temperature weight 0 above the model surface, 0.3 near the wellhead or in the top block, 1.0 at feed zones, 0.5 behind casing or in open-hole isothermal sections, 0.8 elsewhere in open hole). "
       "No departure from the defaults was made. `weight_P` is non-zero only at feed zones.")
    ad("")
    ad("Feed-zone weights:")
    ad("")
    ad(md_table(FZ, ["well", "feedzone", "quality", "ns_confidence", "weight_T", "weight_P"], fmt={"weight_T": 3, "weight_P": 3}))
    ad("")
    ad("**Heat and mass.** Deep hot upflow beneath the " + f"{sh(h1)}–{sh(h2)} corner [interp]; the observed temperature there is about {q('upflow T', S.loc[h1, 'Tmax_degC'], 0)} °C [obs]. Heat is lost by the conductive cap and by lateral outflow, and cold water enters at the margins (elements 1 to 4) [interp]. "
       "Natural discharge locations are not observed in the data; the Tarawera River and the lake are the candidates [interp].")
    ad("")
    ad("**Permeability structure implied.** A low-permeability cap above the 200 °C level; a laterally permeable layer at the outflow depth; permeable vertical paths beneath the hot wells, possibly Fault 3; cold-water access near the north-east and south-west margins [interp].")
    ad("")
    ad(f"**Boundary requirements.** A top boundary at the model surface with the water level near the surface (computed water levels in `well_summary.csv`, column `water_level_mRL`: {q('n extrap', int((S.water_level_extrapolation_m > 0).sum()), 0)} wells extrapolate above the top survey point, the others lie within the survey); a heat/mass source at the base under the hot wells; open or cold-water sides near W4 and W8 [interp].")
    ad("")
    # ---------------------------------------------------------------- 7
    ad("## 7. Data quality issues and open questions")
    ad("")
    iss = []
    iss.append(f"1. **Survey date label.** All surveys are dated {C.info[W[0]]['survey_date'].isoformat()}, after production began ({C.first_prod_field.isoformat()}); the gauge test shows {', '.join(sh(w) for w in gv)} are pre-production [obs]. Resolve which date applies; it controls every confidence factor.")
    iss.append(f"2. **{sh('WELL-9')} gauge miss.** The survey is {q('W9 diff', PC.loc['WELL-9', 'diff_survey_minus_first_bar'], 2)} bar above the first gauge reading; its feed zones lie on the trend, so confidence is medium [obs]. Check the gauge elevation or the first reading.")
    iss.append(f"3. **{sh('WELL-1')} datum.** Wellhead {q('W1 wh2', C.info['WELL-1']['wellhead'], 0)} mRL is {q('W1 below', -G.loc['WELL-1', 'wellhead_minus_surface_m'], 0)} m below the model surface of its column and the survey pressure at the wellhead is {q('W1 Ptop', T[T.well == 'WELL-1'].P_bar.iloc[0], 1)} bar, which puts the water level at {q('W1 wl', S.loc['WELL-1', 'water_level_mRL'], 0)} mRL, above the wellhead [obs]. A wellhead or datum error is possible; its feed-zone pressure is {q('W1 dev2', abs(FZ[FZ.well == 'WELL-1'].P_minus_trend_bar).max(), 1)} bar off the trend, so it is rated possibly disturbed, low.")
    iss.append(f"4. **Shallow hydrostatic ratios above the density of water** in the cased band of the hot wells (up to {q('max cased2', hc.ratio_measured_to_IF97.max(), 2)}), and a cold-well component of {q('cold cased2', C.cold_comp['cased'][0], 3)}; cause not determined [obs].")
    iss.append(f"5. **{sh('WELL-10')} gauge and flow.** The gauge starts {PC.loc['WELL-10', 'first_gauge_date']}, after the well began to flow ({PC.loc['WELL-10', 'own_first_flow_date']}); it has no pre-production test [obs].")
    iss.append(f"6. **{sh('WELL-4')} design.** The liner is listed to {q('W4 liner', min(b for a, b, *_ in C.wells['WELL-4']['design']), 0)} mRL but the track and survey end at {q('W4 track', C.info['WELL-4']['track'][:, 2].min(), 0)} mRL; the survey starts at {q('W4 top', T[T.well == 'WELL-4'].elevation_mRL.max(), 0)} mRL, below the wellhead at {q('W4 wh', C.info['WELL-4']['wellhead'], 0)} mRL [obs].")
    iss.append(f"7. **Water levels** are extrapolations above the shallowest survey point by up to {q('max extrap', S.water_level_extrapolation_m.abs().max(), 0)} m and are uncertain [obs].")
    iss.append("8. **Irregular file names** (`Wel10_*`) and the feed-zone quality tag 'Synthetic' for WELL-10 [obs].")
    for i in iss:
        ad(i)
        ad("")
    # ---------------------------------------------------------------- 8
    ad("## 8. Output file index")
    ad("")
    ad("`Output/Calibration_Input/`:")
    ad("")
    for f_ in ("NS_data_interpretation.md", "NS_calibration_targets.csv", "feedzones.csv", "T_P_at_layer_centres.csv", "well_summary.csv", "figures_index.csv",
               "figures/ (F1 to F7: well_map, profile_<well>, profile_panels_<group>, maps_temperature_bands, maps_pressure_bands, section_A, section_A_P, section_B, section_B_P, PT_state, conceptual_model, gauge_check)"):
        ad(f"- `{f_}`")
    ad("")
    ad("`Output/Additional_Output/`:")
    ad("")
    for f_ in ("input_inventory.csv", "pressure_consistency_check.csv", "hydrostatic_check.csv", "hydrostatic_by_well_band.csv", "hydrostatic_summary.csv",
               "well_geometry_check.csv", "T_P_at_elevation_bands.csv", "T_P_at_common_elevations.csv", "pressure_connections.csv", "data_crosscheck.csv",
               "conceptual_model_key.csv", "section_offsets.csv", "weights_summary.csv", "analysis_parameters.csv", "report_numbers.csv", "run_console_log.txt", "acceptance_check.txt",
               "process_log.md", "figures/ (interpolated companions: section_A_interpolated, section_A_P_interpolated, section_B_interpolated, section_B_P_interpolated, maps_temperature_interpolated, maps_pressure_interpolated)"):
        ad(f"- `{f_}`")
    ad("")
    return "\n".join(out) + "\n"


# =========================================================================================
# 9. Main
# =========================================================================================
def main():
    C = analyse()
    tg = build_targets(C)
    fz = build_fz(C)
    S = build_summary(C)
    lc = build_layer_centres(C)
    H = build_hydrostatic(C)
    geo = build_geometry_check(C)
    X = build_crosscheck(C)
    B = build_bands(C)
    com = common_elevations(C)
    style_setup(C)
    lines = make_lines(C)
    SDs = [section_data(C, l) for l in lines]
    A = SDs[0]
    pairs = [(a["well"], b["well"]) for a, b in zip(A["wells"][:-1], A["wells"][1:])]
    conns = pressure_connections(C, pairs)
    els, ctx = conceptual_elements(C, A)

    # weights summary
    def cat(r):
        return loc_weight(r)[1]
    tg2 = tg.copy()
    tg2["category"] = [cat(r) for _, r in tg2.iterrows()]
    C.weights_summary = tg2.groupby("category").agg(n_points=("weight", "size"), mean_weight=("weight", "mean"),
                                                   max_weight=("weight", "max")).reset_index().round(3)

    # ---- tables ----
    write_csv(tg, os.path.join(CAL, "NS_calibration_targets.csv"))
    write_csv(fz, os.path.join(CAL, "feedzones.csv"))
    write_csv(lc, os.path.join(CAL, "T_P_at_layer_centres.csv"))
    write_csv(S, os.path.join(CAL, "well_summary.csv"))
    pc = C.pc.copy()
    lead = ["well", "gauge_elev_mRL", "survey_P_at_gauge_bar", "first_gauge_date", "first_gauge_P_bar", "diff_survey_minus_first_bar", "verdict", "confidence"]
    pc = pc[lead + [c for c in pc.columns if c not in lead]]
    num = pc.select_dtypes("number").columns
    pc[num] = pc[num].round(3)
    write_csv(pc, os.path.join(ADD, "pressure_consistency_check.csv"))
    write_csv(H, os.path.join(ADD, "hydrostatic_check.csv"))
    write_csv(C.hydro_wb, os.path.join(ADD, "hydrostatic_by_well_band.csv"))
    write_csv(C.hydro_summary, os.path.join(ADD, "hydrostatic_summary.csv"))
    write_csv(geo, os.path.join(ADD, "well_geometry_check.csv"))
    write_csv(B, os.path.join(ADD, "T_P_at_elevation_bands.csv"))
    write_csv(com, os.path.join(ADD, "T_P_at_common_elevations.csv"))
    write_csv(conns, os.path.join(ADD, "pressure_connections.csv"))
    xx = X.copy()
    for c in xx.select_dtypes("number").columns:
        xx[c] = xx[c].round(4)
    write_csv(xx, os.path.join(ADD, "data_crosscheck.csv"))
    write_csv(C.weights_summary, os.path.join(ADD, "weights_summary.csv"))
    write_csv(pd.DataFrame([dict(line=l["id"], x0=round(l["p_ext0"][0], 1), y0=round(l["p_ext0"][1], 1), x1=round(l["p_ext1"][0], 1),
                                 y1=round(l["p_ext1"][1], 1), length_m=round(l["length"], 1), well=w["well"], offset_m=round(w["offset"], 1),
                                 distance_along_line_m=round(w["s_mean"], 1)) for l, SD in zip(lines, SDs) for w in SD["wells"]]),
              os.path.join(ADD, "section_offsets.csv"))
    write_csv(pd.DataFrame([dict(element=e["n"], name=e["name"], wells=";".join(e["wells"]), claim_type=e["claim"],
                                 arrow="solid" if e["solid"] else "dashed", supporting_observation=e["obs"]) for e in els]),
              os.path.join(ADD, "conceptual_model_key.csv"))
    params = [("pressure_reading_precision_bar", PREC_P), ("feed_zone_trend_intercept_bar", round(C.trend[0], 4)),
              ("feed_zone_trend_slope_bar_per_m", round(C.trend[1], 6)), ("trend_rms_bar", round(C.trend_rms, 3)),
              ("tolerance_bar", round(C.tol, 3)), ("verified_wells", ";".join(C.verified)), ("section_offset_A_m", lines[0]["max_offset"]),
              ("section_offset_B_m", lines[1]["max_offset"]), ("T_colour_limit_degC", C.tmax_lim), ("elevation_axis_mRL", f"{C.zlim[0]:.0f} to {C.zlim[1]:.0f}"),
              ("gravity_m_s2", G), ("first_production_field", C.first_prod_field.isoformat()),
              ("first_injection_field", C.first_inj_field.isoformat() if C.first_inj_field else ""),
              ("injectate_T_degC", ";".join(f"{w}:{v['T']:.2f}" for w, v in C.inj.items()))]
    write_csv(pd.DataFrame(params, columns=["parameter", "value"]), os.path.join(ADD, "analysis_parameters.csv"))

    # ---- figures ----
    for f in glob.glob(os.path.join(CALFIG, "*.png")) + glob.glob(os.path.join(ADDFIG, "*.png")):
        os.remove(f)
    FIGS.clear()
    fig_well_map(C, lines)
    for w in C.wells_order:
        fig_profile_well(C, w)
    for gname, wl in C.groups.items():
        fig_profile_panels(C, gname, [w for w in wl if w in C.wells])
    fig_band_maps(C, "T")
    fig_band_maps(C, "P")
    for SD in SDs:
        for kind in ("T", "P"):
            fig_section(C, SD, kind)
            fig_section(C, SD, kind, interp=True)
    fig_band_maps_interp(C, "T")
    fig_band_maps_interp(C, "P")
    fig_PT(C)
    fig_conceptual(C, A, els, ctx, conns)
    fig_gauge(C)
    order = {"F1": 1, "F2": 2, "F3": 3, "F4": 4, "F5": 5, "F6": 6, "F7": 7}
    idx = pd.DataFrame(FIGS)
    idx["o"] = idx.figure_id.map(order)
    idx = idx.sort_values(["o", "file"], kind="stable").drop(columns="o")
    write_csv(idx, os.path.join(CAL, "figures_index.csv"))

    # ---- report ----
    Q = Quoted()
    rep = build_report(C, lines, SDs, els, ctx, conns, Q)
    write_csv(pd.DataFrame(Q.rows), os.path.join(ADD, "report_numbers.csv"))
    write_text(os.path.join(CAL, "NS_data_interpretation.md"), rep)
    log(f"wells: {len(C.wells_order)}; survey points: {len(tg)}; feed zones: {len(fz)}; figures indexed: {len(idx)}")
    log("verdicts: " + "; ".join(f"{w}={C.verdict[w]['verdict']}/{C.verdict[w]['confidence']}" for w in C.wells_order))
    write_text(os.path.join(ADD, "run_console_log.txt"), "\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
