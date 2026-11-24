#!/usr/bin/env python
"""Natural-state data interpretation for geothermal model calibration.

Run from the repository root:   python Scripts/natural_state_analysis.py

Reads only Baseline/. Writes the input inventory, every table, every figure and the main report
(Output/NS_data_interpretation.md). Output/process_log.md is written by hand; the acceptance
check is Scripts/check_deliverables.py.

Sections of this file
  1  constants and I/O helpers
  2  IAPWS-IF97 (Region 4 saturation, Region 1 liquid) with self-validation
  3  input inventory (CLAUDE.md 2.1)
  4  readers: grid, well-to-block file, well JSON, plain-column files, mapping
  5  per-well analysis (profiles, weights, summaries)
  6  field-wide analysis (integrity, natural-state timing, hydrostatics, structure)
  7  tables, figures, report
"""
import ast
import csv
import datetime as dt
import glob
import json
import math
import os
import re
import sys
import warnings

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.path import Path as MplPath          # noqa: E402
import matplotlib.tri as mtri                        # noqa: E402

warnings.filterwarnings("ignore", category=RuntimeWarning)

# ======================================================================================
# 1. constants and helpers
# ======================================================================================
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(ROOT, "Baseline")
OUT = os.path.join(ROOT, "Output")
FIG = os.path.join(OUT, "figures")

G = 9.80665                 # m/s2
P_CRIT_BAR = 220.64
T_CRIT_C = 373.946
P_ATM_BAR = 1.01325

# weights (CLAUDE.md section 5)
W_FZ, W_OPEN, W_ISO, W_CASED, W_TOP, W_ABOVE = 1.0, 0.8, 0.5, 0.5, 0.3, 0.0
CONF_FACTOR = {"high": 1.0, "medium": 0.75, "low": 0.5}
# analysis parameters (all reported in Output/analysis_parameters.csv)
PRM = dict(
    gauge_match_bar=0.15,        # reading precision of the gauge series is 0.1 bar (one decimal)
    survey_precision_bar=0.1,
    iso_range_degC=5.0,          # isothermal section: T range of the section
    iso_min_length_m=300.0,
    inversion_min_degC=3.0,
    shallow_fit_m=300.0,
    deep_min_length_m=300.0,
    hydro_min_margin_degC=10.0,
    hydro_min_P_bar=1.5,
    hydro_min_length_m=300.0,
    bpd_margin_degC=5.0,
    top_block_m=30.0,
    wl_points=3,
    fz_trend_ok_bar=3.0,
    fz_trend_flag_bar=6.0,
    cold_well_Tmax_degC=100.0,
    k_rock_W_mK=2.0,
    hot_core_Tmax_degC=255.0,
    hot_layer_drop_degC=15.0,
)


def log(*a):
    print(*a, flush=True)


def fmt(x, n=1):
    """Format a number for tables/report; blank for NaN."""
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return f"{x:.{n}f}"


def is_nan(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


def interp_z(z, v, zq):
    """Linear interpolation of v(z) at zq; NaN outside the sampled range (never extrapolates).
    z may be in any order."""
    z = np.asarray(z, float)
    v = np.asarray(v, float)
    o = np.argsort(z)
    z, v = z[o], v[o]
    zq_a = np.atleast_1d(np.asarray(zq, float))
    r = np.interp(zq_a, z, v)
    r[(zq_a < z[0] - 1e-9) | (zq_a > z[-1] + 1e-9)] = np.nan
    return r if np.ndim(zq) else float(r[0])


# ======================================================================================
# 2. IAPWS-IF97
# ======================================================================================
_N4 = [None, 0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5,
       -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6,
       -0.23855557567849, 0.65017534844798e3]


def tsat_c(p_bar):
    """Region 4 (eq. 31): saturation temperature in degC for pressure in bar.
    NaN above the critical pressure (no Tsat) or below the triple-point pressure."""
    p = p_bar / 10.0
    if not (611.213e-6 <= p <= 22.064):
        return float("nan")
    n = _N4
    b = p ** 0.25
    e = b * b + n[3] * b + n[6]
    f = n[1] * b * b + n[4] * b + n[7]
    g = n[2] * b * b + n[5] * b + n[8]
    d = 2 * g / (-f - math.sqrt(f * f - 4 * e * g))
    t = (n[10] + d - math.sqrt((n[10] + d) ** 2 - 4 * (n[9] + n[10] * d))) / 2
    return t - 273.15


def psat_bar(t_c):
    """Region 4 (eq. 30): saturation pressure in bar for T in degC (273.15 K .. 647.096 K)."""
    t = t_c + 273.15
    n = _N4
    th = t + n[9] / (t - n[10])
    a = th * th + n[1] * th + n[2]
    b = n[3] * th * th + n[4] * th + n[5]
    c = n[6] * th * th + n[7] * th + n[8]
    return 10.0 * (2 * c / (-b + math.sqrt(b * b - 4 * a * c))) ** 4


# Region 1 coefficients (IF97 Table 2): I, J, n
_R1 = [
    (0, -2, 0.14632971213167), (0, -1, -0.84548187169114), (0, 0, -0.37563603672040e1),
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
_R1I = np.array([r[0] for r in _R1], float)
_R1J = np.array([r[1] for r in _R1], float)
_R1N = np.array([r[2] for r in _R1], float)
_RGAS = 0.461526  # kJ/kg/K


def region1(t_c, p_bar):
    """Region 1 specific volume (m3/kg) and enthalpy (kJ/kg) at T [degC], p [bar]."""
    t = t_c + 273.15
    p = p_bar / 10.0
    pi, tau = p / 16.53, 1386.0 / t
    a, b = 7.1 - pi, tau - 1.222
    gp = np.sum(-_R1N * _R1I * a ** (_R1I - 1) * b ** _R1J)
    gt = np.sum(_R1N * a ** _R1I * _R1J * b ** (_R1J - 1))
    v = _RGAS * t / (p * 1000.0) * pi * gp
    h = _RGAS * t * tau * gt
    return float(v), float(h)


def rho_liq(t_c, p_bar):
    """Liquid density (kg/m3), Region 1."""
    return 1.0 / region1(t_c, p_bar)[0]


def t_from_h_liq(h, p_bar):
    """Temperature (degC) of compressed liquid with enthalpy h (kJ/kg), by bisection (Region 1)."""
    lo, hi = 0.0, 350.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if region1(mid, p_bar)[1] < h:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def validate_if97():
    """Check the implementation against IF97 verification tables; fail if any error exceeds 0.1 K.
    Returns a list of (what, expected, computed, error) rows for the log."""
    rows = []
    # Region 4, saturation temperature (IF97 Table 36): p [MPa] -> T [K]
    for p, tk in [(0.1, 372.755919), (1.0, 453.035632), (10.0, 584.149488)]:
        c = tsat_c(p * 10) + 273.15
        rows.append((f"Tsat({p * 10:g} bar) [K]", tk, c, abs(c - tk)))
    # reference values quoted in the task
    rows.append(("Tsat(1.01325 bar) [degC]", 99.97, tsat_c(1.01325), abs(tsat_c(1.01325) - 99.97)))
    rows.append(("Tsat(100 bar) [degC]", 311.0, tsat_c(100.0), abs(tsat_c(100.0) - 311.0)))
    # Region 4, saturation pressure (IF97 Table 35): T [K] -> p [MPa]
    for tk, p in [(300.0, 0.353658941e-2), (500.0, 0.263889776e1), (600.0, 0.123443146e2)]:
        c = psat_bar(tk - 273.15) / 10
        # express the pressure error as an equivalent temperature error
        te = abs(tsat_c(c * 10) + 273.15 - tk)
        rows.append((f"Psat({tk:g} K) [MPa]", p, c, te))
    # Region 1 (IF97 Table 5): v [m3/kg], h [kJ/kg]
    for tk, p, v, h in [(300.0, 3.0, 0.100215168e-2, 0.115331273e3),
                        (300.0, 80.0, 0.971180894e-3, 0.184142828e3),
                        (500.0, 3.0, 0.120241800e-2, 0.975542239e3)]:
        cv, ch = region1(tk - 273.15, p * 10)
        # express as equivalent temperature error through h (dh/dT ~ cp ~ 4.2 kJ/kg/K)
        rows.append((f"Region1 v({tk:g} K,{p:g} MPa) [m3/kg]", v, cv, abs(cv - v) / v * 1000.0))
        rows.append((f"Region1 h({tk:g} K,{p:g} MPa) [kJ/kg]", h, ch, abs(ch - h) / 4.2))
    worst = max(r[3] for r in rows)
    if worst > 0.1:
        raise SystemExit(f"IF97 implementation validation FAILED: worst equivalent error {worst:.4f} K")
    return rows


# ======================================================================================
# 3. input inventory (CLAUDE.md 2.1)
# ======================================================================================
def rel(p):
    return os.path.relpath(p, ROOT).replace("\\", "/")


def is_grid_file(path):
    """A MULgraph grid file holds VERTICES, GRID, LAYERS and SURFA sections."""
    with open(path, errors="ignore") as fh:
        head = fh.read(3000)
    return "VERTICES" in head and "GENER" in head


def build_inventory():
    """Assign a role to every file and directory under Baseline/. Rules use file kinds, not names of
    wells or dates; anything unmatched becomes out_of_scope and is reported."""
    rows = []
    covered = []

    def add(role, path, note):
        rows.append((role, rel(path), note))
        covered.append(os.path.normpath(path))

    all_files = []
    for d, _, fs in os.walk(BASE):
        for f in fs:
            all_files.append(os.path.normpath(os.path.join(d, f)))
    all_dirs = [os.path.normpath(os.path.join(d, x)) for d, ds, _ in os.walk(BASE) for x in ds]

    def jsons_in(dirname):
        return [p for p in all_dirs if os.path.basename(p) == dirname]

    # directories of per-well data
    for d in all_dirs:
        n = os.path.basename(d).lower()
        files = [f for f in all_files if os.path.dirname(f) == d]
        if n == "data_json" and files:
            add("well_data_json", d, "per-well JSON: track, design, feed zones, T/P surveys, transients, flows")
        elif n == "data_dat" and files:
            add("well_data_dat", d, "plain-column copies of the well data (cross-check only)")
        elif n == "mapping" and files:
            add("mapping", d, "faults, rivers and lake outline polylines")
        elif n == "plotting" and files:
            add("reference_other", d, "plotting set-up with hard-coded paths from other machines; not data")
        elif n == "futurescenarios" or n == "productionhistory" or n == "production":
            add("out_of_scope", d, "production / forecast set-up; not read")
        elif not files and not any(os.path.dirname(x) == d for x in all_dirs):
            add("out_of_scope", d, "empty directory")
    # individual files
    for f in all_files:
        if any(f.startswith(c + os.sep) or f == c for c in covered):
            continue
        b = os.path.basename(f)
        low = b.lower()
        if low == "groups_data.json":
            add("well_groups", f, "list of wells and their groups")
        elif low == "area_data.json":
            add("out_of_scope", f, "field-total steam flow and power history (production); not read")
        elif low.endswith(".dat") and is_grid_file(f):
            add("geometry_dat", f, "MULgraph model grid: vertices, grid, layers, surface, well tracks")
        elif low.endswith("_wells.json"):
            add("well_blocks_json", f, "well-to-block intersections, feed-zone blocks, interpolation weights")
        elif low.endswith("_boundary_conditions.json"):
            add("reference_bc", f, "reference model boundary conditions (atmosphere, rain, heat, wet columns)")
        elif low.endswith("_dict_rocktype.json"):
            add("reference_rock", f, "reference model rock-type dictionary")
        elif low.endswith(".auxjson"):
            add("reference_other", f, "viewer scene (GIS layers, slices, display settings); hard-coded paths")
        elif low.endswith(".incon"):
            add("reference_results", f, "reference natural-state initial conditions / results")
        elif low.endswith(".dat") and "naturalstate" in f.lower().replace("\\", "/"):
            add("reference_rock", f, "reference natural-state model input (rock types and permeabilities, parameters)")
        elif low.startswith("current") and low.endswith(".json"):
            if any(k in f.lower() for k in ("futurescen", "production", "currentprsetup", "currentfssetup")):
                add("out_of_scope", f, "production / forecast set-up; not read")
            else:
                add("reference_other", f, "reference model set-up / configuration")
        else:
            add("out_of_scope", f, "not needed for the natural-state analysis (e.g. production spreadsheets, OS files)")
    # make sure every file is covered
    for f in all_files:
        if not any(f.startswith(c + os.sep) or f == c for c in covered):
            add("out_of_scope", f, "unclassified; not read")
    rows.sort(key=lambda r: (r[0], r[1]))
    return rows


def write_inventory(rows):
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, "input_inventory.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["role", "path", "note"])
        w.writerows(rows)


def role_paths(rows, role):
    return [os.path.join(ROOT, r[1]) for r in rows if r[0] == role]


# ======================================================================================
# 4. readers
# ======================================================================================
def read_grid(path):
    """MULgraph grid file. Sections are found by their keyword line. Vertex lines have the two
    coordinates run together without a separator, so they are read with a fixed pattern
    (3-character name, then two numbers with two decimals); the other sections are split on
    whitespace."""
    lines = open(path, encoding="utf-8", errors="ignore").read().splitlines()
    idx = {}
    for i, ln in enumerate(lines):
        if ln.strip() in ("VERTICES", "GRID", "CONNECTIONS", "LAYERS", "SURFA", "WELLS"):
            idx[ln.strip()] = i
    order = sorted(idx.items(), key=lambda kv: kv[1])

    def section(name):
        i = idx[name]
        nxt = [j for _, j in order if j > i]
        return lines[i + 1:(nxt[0] if nxt else len(lines))]

    verts = {}
    for ln in section("VERTICES"):
        if not ln.strip():
            continue
        m = re.match(r"^(\w{3})(\d{7}\.\d{2})(\d{7}\.\d{2})\s*$", ln)
        if m:
            verts[m.group(1)] = (float(m.group(2)), float(m.group(3)))
        else:
            p = ln.split()
            if len(p) == 3:
                verts[p[0]] = (float(p[1]), float(p[2]))
            else:
                raise ValueError(f"cannot parse vertex line: {ln!r}")
    cols, g = {}, section("GRID")
    i = 0
    while i < len(g):
        ln = g[i]
        if not ln.strip():
            i += 1
            continue
        if len(ln.strip()) > 3:                       # header: name, then vertex count
            name = ln[:3]
            n = int(ln[3:].split()[-1])
            cols[name] = [g[i + 1 + k].strip() for k in range(n)]
            i += 1 + n
        else:
            i += 1
    layers = {}
    for ln in section("LAYERS"):
        p = ln.split()
        if len(p) == 3:
            layers[int(p[0])] = (float(p[1]), float(p[2]))
    lay = {}
    for i in sorted(layers):
        if i > 0:
            lay[i] = dict(top=layers[i - 1][0], bottom=layers[i][0], centre=layers[i][1])
    surf = {}
    for ln in section("SURFA"):
        p = ln.split()
        if len(p) == 2:
            surf[p[0]] = float(p[1])
    wells = {}
    for ln in section("WELLS"):
        p = ln.split()
        if len(p) == 4:
            wells.setdefault(p[0], []).append((float(p[1]), float(p[2]), float(p[3])))
    polys = {c: [verts[v] for v in vs] for c, vs in cols.items()}
    paths = {c: MplPath(np.array(pl + [pl[0]])) for c, pl in polys.items()}
    return dict(verts=verts, cols=cols, polys=polys, paths=paths, layers=lay, surf=surf, wells=wells,
                top_surface=layers[0][0])


def column_at(grid, x, y):
    for c, p in grid["paths"].items():
        if p.contains_point((x, y)):
            return c
    return None


def layer_at(grid, z):
    """Layer index containing z (layer i spans bottom_i .. top_i); the upper layer on a boundary."""
    for i, L in grid["layers"].items():
        if L["bottom"] < z <= L["top"]:
            return i
    return None


def read_well_blocks(path):
    return json.load(open(path, encoding="utf-8"))


def block_candidates(wb, well, z, tol=0.5):
    """Blocks of the well-to-block file whose interval contains z (within tol)."""
    hits = []
    for b in wb[well]["Blocks"]:
        top, bot = max(b["interval"]), min(b["interval"])
        if bot - tol <= z <= top + tol:
            hits.append((b["block"], top, bot))
    return hits


def block_at(wb, well, z):
    """Block for elevation z from the well-to-block file (the authority). On a block boundary the
    upper block is used."""
    hits = block_candidates(wb, well, z)
    if not hits:
        return None
    strict = [h for h in hits if h[2] < z <= h[1] + 1e-9]
    return (strict or hits)[0][0]


def split_block(name):
    m = re.match(r"^([a-z]{3})(\d+)$", name or "")
    return (m.group(1), int(m.group(2))) if m else (None, None)


def pdate(s):
    return dt.datetime.strptime(s, "%d/%m/%Y").date()


def years_since(d, ref):
    return (d - ref).days / 365.25


def lit(s):
    """Literal evaluation of the list-valued strings in the well JSON files."""
    return ast.literal_eval(s) if isinstance(s, str) else s


def read_well_json(path):
    d = json.load(open(path, encoding="utf-8"))
    w = dict(raw=d, name=d["Well_Name"], file=rel(path))
    w["num"] = int(re.search(r"(\d+)$", w["name"]).group(1))
    w["fz"] = []
    for k, s in enumerate(d["Feed_Zones"]["Data"], 1):
        top, bot, q = lit(s)
        w["fz"].append(dict(n=k, top=float(top), bottom=float(bot), mid=0.5 * (top + bot),
                            quality=str(q), length=abs(top - bot)))
    w["track"] = np.array([lit(s) for s in d["Well_Track"]["Data"]], float)
    w["design"] = []
    for s in d["Well_Design"]["Data"]:
        top, bot, dia, slot, spec = lit(s)
        w["design"].append(dict(top=float(top), bottom=float(bot), diameter=str(dia),
                                slotted=bool(slot), spec=str(spec)))
    w["design_source"] = d["Well_Design"].get("Source", "")
    w["wellhead_z"] = float(w["track"][0, 2])
    w["td_z"] = float(w["track"][-1, 2])
    w["wellhead_xy"] = (float(w["track"][0, 0]), float(w["track"][0, 1]))
    cas = [x for x in w["design"] if not x["slotted"]]
    slo = [x for x in w["design"] if x["slotted"]]
    w["shoe_z"] = min(x["bottom"] for x in cas) if cas else float("nan")
    w["casing_top_z"] = max(x["top"] for x in cas) if cas else float("nan")
    w["liner_top_z"] = max(x["top"] for x in slo) if slo else float("nan")
    w["liner_bottom_z"] = min(x["bottom"] for x in slo) if slo else float("nan")
    # surveys: keyed by date, then by elevation (string)
    w["surveys"] = {}
    for var, key in (("T", "Downhole_Temperature"), ("P", "Downhole_Pressure")):
        for date, s in d.get(key, {}).items():
            dd = {float(z): float(v) for z, v in s["Data"].items()}
            w["surveys"].setdefault(date, {})[var] = dict(
                data=dd, source=s.get("Source", ""), comment=s.get("Comment", ""),
                stable=s.get("StableFlag", None), units=s.get("Units", {}))
    # records of flows, enthalpies and transient pressures
    w["records"] = {}
    for key, v in d.items():
        if key in ("Production_Mass_Flow", "Injection_Mass_Flow", "Transient_Enthalpy", "Injection_Enthalpy"):
            w["records"][key] = dict(
                series={pdate(a): float(b) for a, b in v["Data"].items()}, units=v.get("Units", {}),
                source=v.get("Source", ""), comment=v.get("Comment", ""))
    w["transient_P"] = {}
    for z, v in d.get("Transient_Pressure", {}).items():
        w["transient_P"][float(z)] = dict(
            series={pdate(a): float(b) for a, b in v["Data"].items()}, units=v.get("Units", {}),
            source=v.get("Source", ""), comment=v.get("Comment", ""))
    return w


def track_xy_at(w, z):
    """Horizontal position of the well track at elevation z (linear along the track)."""
    tr = w["track"]
    zs = tr[:, 2]
    if len(tr) == 2 and np.allclose(tr[0, :2], tr[1, :2]):
        return float(tr[0, 0]), float(tr[0, 1])
    o = np.argsort(zs)
    x = np.interp(z, zs[o], tr[o, 0])
    y = np.interp(z, zs[o], tr[o, 1])
    return float(x), float(y)


def read_dat2(path):
    a = np.loadtxt(path, ndmin=2)
    return a[:, 0], a[:, 1]


def read_mapping(path):
    """Features as given in the file: name -> dict(type, segments=[ndarray...])."""
    m = json.load(open(path, encoding="utf-8"))
    feats = {}
    for name, v in m.items():
        data = v["Data"]
        if len(data) and isinstance(data[0][0], (list, tuple)):       # ragged: list of segments
            segs = [np.array(s, float) for s in data]
        else:
            segs = [np.array(data, float)]
        feats[name] = dict(type=v.get("type", ""), segments=segs,
                           meta={k: x for k, x in v.items() if k != "Data"})
    return feats


def seg_dist(px, py, a, b):
    ax, ay = a
    bx, by = b
    dx, dy = bx - ax, by - ay
    L2 = dx * dx + dy * dy
    t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((px - ax) * dx + (py - ay) * dy) / L2))
    return math.hypot(px - (ax + t * dx), py - (ay + t * dy))


def dist_to_feature(px, py, feat):
    best = float("inf")
    for s in feat["segments"]:
        for i in range(len(s) - 1):
            best = min(best, seg_dist(px, py, s[i], s[i + 1]))
    return best


def feat_key(name):
    return re.sub(r"[^A-Za-z0-9]+", "_", name).strip("_")


# ======================================================================================
# 5. per-well analysis
# ======================================================================================
def tsat_arr(P):
    return np.array([tsat_c(p) for p in P])


def make_profile(w):
    """Survey profile of a well (elevations descending). One survey date per well in this dataset;
    if several exist the earliest is used and the others are listed in w['other_dates']."""
    dates = sorted(w["surveys"], key=pdate)
    date = dates[0]
    w["survey_date"] = date
    w["other_dates"] = dates[1:]
    sv = w["surveys"][date]
    zs = sorted(set(sv["T"]["data"]) | set(sv["P"]["data"]), reverse=True)
    T = np.array([sv["T"]["data"].get(z, np.nan) for z in zs])
    P = np.array([sv["P"]["data"].get(z, np.nan) for z in zs])
    ts = tsat_arr(P)
    return dict(z=np.array(zs), T=T, P=P, Tsat=ts, margin=ts - T)


def in_interval(z, a, b, tol=1e-6):
    lo, hi = min(a, b), max(a, b)
    return lo - tol <= z <= hi + tol


def isothermal_sections(z, T):
    """Maximal sections whose temperature range is <= iso_range and whose length is >= iso_min_length.
    Works on the survey points (z descending). Returns a list of (z_top, z_bottom, Tmin, Tmax)."""
    n = len(z)
    wins = []
    for i in range(n):
        j = i
        while j + 1 < n and (np.nanmax(T[i:j + 2]) - np.nanmin(T[i:j + 2])) <= PRM["iso_range_degC"]:
            j += 1
        if j > i and z[i] - z[j] >= PRM["iso_min_length_m"]:
            wins.append((z[i], z[j], float(np.nanmin(T[i:j + 1])), float(np.nanmax(T[i:j + 1]))))
    # merge overlapping windows
    wins.sort(key=lambda a: -a[0])
    merged = []
    for a in wins:
        if merged and a[0] >= merged[-1][1]:
            m = merged[-1]
            merged[-1] = (m[0], min(m[1], a[1]), min(m[2], a[2]), max(m[3], a[3]))
        else:
            merged.append(a)
    return merged


def lsq_gradient(z, T, min_n=3):
    """Temperature gradient in degC/km (positive when T increases with depth), least squares."""
    if len(z) < min_n:
        return float("nan")
    s = np.polyfit(np.asarray(z, float), np.asarray(T, float), 1)[0]
    return float(-s * 1000.0)


def first_downward_crossing(z, T, thr):
    """Elevation of the first crossing of thr going downward (z descending), linear interpolation."""
    for i in range(len(z) - 1):
        if T[i] < thr <= T[i + 1] and T[i + 1] != T[i]:
            return float(z[i] + (thr - T[i]) / (T[i + 1] - T[i]) * (z[i + 1] - z[i]))
        if T[i] == thr and i == 0:
            return float(z[0])
    return float("nan")


def estimate_water_level(prof):
    """Water level (mRL): straight line P(z) through the shallowest `wl_points` survey points with
    P >= hydro_min_P_bar, solved for atmospheric pressure. May lie above the top survey point."""
    idx = [i for i in range(len(prof["z"])) if prof["P"][i] >= PRM["hydro_min_P_bar"]][:PRM["wl_points"]]
    if len(idx) < 2:
        return float("nan"), float("nan")
    z, P = prof["z"][idx], prof["P"][idx]
    b, a = np.polyfit(z, P, 1)
    if b >= 0:
        return float("nan"), float("nan")
    wl = (P_ATM_BAR - a) / b
    return float(wl), float(-b)


QUALITY_RANK_DEFAULT = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}


def quality_scheme(wells):
    """Map every quality tag found to a factor (CLAUDE.md section 5): highest tag 1.0, intermediate
    0.8, lowest tag or any tag that is not a quality grade 0.5."""
    tags = sorted({fz["quality"] for w in wells for fz in w["fz"]})
    grades = [t for t in tags if t.upper() in QUALITY_RANK_DEFAULT]
    mapping, notes = {}, {}
    ordered = sorted(grades, key=lambda t: -QUALITY_RANK_DEFAULT[t.upper()])
    for k, t in enumerate(ordered):
        if k == 0:
            mapping[t], notes[t] = 1.0, "highest quality grade in the dataset"
        elif k == len(ordered) - 1:
            mapping[t], notes[t] = 0.5, "lowest quality grade in the dataset"
        else:
            mapping[t], notes[t] = 0.8, "intermediate quality grade"
    for t in tags:
        if t not in mapping:
            mapping[t], notes[t] = 0.5, "not a quality grade (treated like the lowest grade)"
    return mapping, notes


def dominant_feedzone(w, qmap):
    """Pivot feed zone: highest quality factor, then the thicker zone, then the shallower zone."""
    if not w["fz"]:
        return None
    return sorted(w["fz"], key=lambda f: (-qmap[f["quality"]], -f["length"], -f["mid"]))[0]


def analyse_well(w, grid, wb, qmap, surf_info):
    """Profile-derived characteristics of one well (before the natural-state verdict)."""
    prof = make_profile(w)
    z, T, P = prof["z"], prof["T"], prof["P"]
    r = dict(well=w["name"])
    i = int(np.nanargmax(T))
    r["Tmax_degC"], r["Tmax_elev_mRL"] = float(T[i]), float(z[i])
    r["T_at_TD_degC"], r["TD_mRL"] = float(T[-1]), float(z[-1])
    r["track_TD_mRL"] = w["td_z"]
    # inversion: largest drop below the running maximum
    runmax = np.maximum.accumulate(T)
    zrun = []
    cur = 0
    for k in range(len(T)):
        if T[k] >= T[cur]:
            cur = k
        zrun.append(z[cur])
    drop = runmax - T
    k = int(np.argmax(drop))
    if drop[k] >= PRM["inversion_min_degC"]:
        r["inversion_degC"], r["inversion_top_mRL"], r["inversion_base_mRL"] = float(drop[k]), float(zrun[k]), float(z[k])
    else:
        r["inversion_degC"] = r["inversion_top_mRL"] = r["inversion_base_mRL"] = float("nan")
    wl, slope = estimate_water_level(prof)
    r["water_level_mRL"] = wl
    r["water_level_extrapolation_m"] = max(0.0, wl - z[0]) if not is_nan(wl) else float("nan")
    # shallow gradient: top 300 m of the liquid-filled survey
    zstart = min(z[0], wl) if not is_nan(wl) else z[0]
    sel = (z <= zstart + 1e-9) & (z >= zstart - PRM["shallow_fit_m"] - 1e-9)
    r["shallow_gradient_degC_per_km"] = lsq_gradient(z[sel], T[sel])
    r["shallow_fit_n"] = int(sel.sum())
    # deep gradient: below both the inversion base and the deepest feed zone
    lims = [min(f["bottom"], f["top"]) for f in w["fz"]]
    zdeep = min(lims) if lims else float("inf")
    if not is_nan(r["inversion_base_mRL"]):
        zdeep = min(zdeep, r["inversion_base_mRL"])
    if np.isfinite(zdeep) and zdeep - z[-1] >= PRM["deep_min_length_m"]:
        sel = z <= zdeep + 1e-9
        r["deep_gradient_degC_per_km"] = lsq_gradient(z[sel], T[sel], min_n=2)   # two points are a (exact) fit
        r["deep_fit_top_mRL"] = float(zdeep)
        r["deep_fit_n"] = int(sel.sum())
    else:
        r["deep_gradient_degC_per_km"] = float("nan")
        r["deep_fit_top_mRL"] = float(zdeep) if np.isfinite(zdeep) else float("nan")
        r["deep_fit_n"] = 0
    r["z_100C_mRL"] = first_downward_crossing(z, T, 100.0)
    r["z_200C_mRL"] = first_downward_crossing(z, T, 200.0)
    # boiling margin
    m = prof["margin"]
    ok = ~np.isnan(m)
    if ok.any():
        j = int(np.nanargmin(np.where(ok, m, np.inf)))
        r["min_boiling_margin_degC"], r["min_boiling_margin_elev_mRL"] = float(m[j]), float(z[j])
    else:
        r["min_boiling_margin_degC"] = r["min_boiling_margin_elev_mRL"] = float("nan")
    # pivot at the dominant feed zone
    dom = dominant_feedzone(w, qmap)
    if dom:
        r["pivot_elev_mRL"] = dom["mid"]
        r["pivot_P_bar"] = float(interp_z(z, P, dom["mid"]))
        r["pivot_feedzone"] = dom["n"]
    else:
        r["pivot_elev_mRL"] = r["pivot_P_bar"] = float("nan")
        r["pivot_feedzone"] = 0
    # isothermal sections
    iso = isothermal_sections(z, T)
    r["n_iso_sections"] = len(iso)
    # regime (rules in the report, section 3)
    flags_bpd = (m < PRM["bpd_margin_degC"]) & (P >= PRM["hydro_min_P_bar"])
    r["bpd_points"] = int(np.nansum(flags_bpd))
    if r["bpd_points"] > 0:
        regime = "two-phase"
    elif r["Tmax_degC"] < PRM["cold_well_Tmax_degC"]:
        regime = "recharge/downflow"
    elif not is_nan(r["inversion_degC"]) and r["inversion_degC"] >= 10.0:
        regime = "outflow"
    elif r["Tmax_degC"] >= 200.0 and r["Tmax_degC"] - r["T_at_TD_degC"] <= 5.0:
        regime = "upflow"
    else:
        regime = "conductive"
    r["regime"] = regime
    return prof, r, iso


# ======================================================================================
# 6. field-wide analysis
# ======================================================================================
def hydro_windows(prof, w):
    """Non-overlapping single-phase-liquid windows of >= hydro_min_length m (CLAUDE.md 4.3).
    Windows are cut at the casing shoe (they never straddle it) and at feed-zone mid-points where
    the minimum length allows. Returns a list of (z_top, z_bottom, band) and a list of notes."""
    z, T, P, m = prof["z"], prof["T"], prof["P"], prof["margin"]
    ok = (np.nan_to_num(m, nan=999.0) >= PRM["hydro_min_margin_degC"]) & (P >= PRM["hydro_min_P_bar"])
    runs, i = [], 0
    while i < len(z):
        if ok[i]:
            j = i
            while j + 1 < len(z) and ok[j + 1]:
                j += 1
            runs.append((z[i], z[j]))
            i = j + 1
        else:
            i += 1
    shoe = w["shoe_z"]
    mids = sorted([f["mid"] for f in w["fz"]], reverse=True)
    wins, notes = [], []
    Lmin = PRM["hydro_min_length_m"]
    for zt, zb in runs:
        pieces = []
        if not is_nan(shoe) and zb < shoe < zt:
            pieces = [(zt, shoe, "cased"), (shoe, zb, "open")]
        elif not is_nan(shoe) and zt <= shoe + 1e-9:
            pieces = [(zt, zb, "open")]
        else:
            pieces = [(zt, zb, "cased")]
        for a, b, band in pieces:
            if a - b < Lmin:
                notes.append(f"{w['name']}: {band} piece {a:.0f} to {b:.0f} mRL is {a - b:.0f} m (< {Lmin:.0f} m), not tested")
                continue
            start = a
            while True:
                q = next((s for s in mids if a > s > b and start - s >= Lmin and s - b >= Lmin), None)
                if q is None:
                    wins.append((start, b, band))
                    break
                wins.append((start, q, band))
                start = q
    return wins, notes


def window_stats(prof, za, zb):
    """Measured and IF97 pressure change between za (top) and zb (bottom)."""
    z, T, P = prof["z"], prof["T"], prof["P"]
    inside = [k for k in range(len(z)) if zb < z[k] < za]
    zz = [za] + [z[k] for k in inside] + [zb]
    TT = [float(interp_z(z, T, q)) for q in zz]
    PP = [float(interp_z(z, P, q)) for q in zz]
    rho = [rho_liq(t, p) for t, p in zip(TT, PP)]
    dP_if97 = sum(0.5 * (rho[k] + rho[k + 1]) * G * (zz[k] - zz[k + 1]) for k in range(len(zz) - 1)) / 1e5
    dP_meas = PP[-1] - PP[0]
    L = za - zb
    return dict(dP_measured=dP_meas, dP_IF97=dP_if97, ratio=dP_meas / dP_if97,
                rho_implied=dP_meas * 1e5 / (G * L), length=L,
                # each reading is rounded to 0.1 bar (error <= 0.05), so dP is uncertain by <= 0.1 bar
                rounding_ratio=PRM["survey_precision_bar"] / dP_if97,
                rounding_rho=PRM["survey_precision_bar"] * 1e5 / (G * L),
                mean_T=float(np.mean(TT)))


def first_flow(w):
    """Date of the first non-zero production or injection record and its kind."""
    best = None
    for key, kind in (("Production_Mass_Flow", "production"), ("Injection_Mass_Flow", "injection")):
        rec = w["records"].get(key)
        if rec:
            for d, v in sorted(rec["series"].items()):
                if abs(v) > 0:
                    if best is None or d < best[0]:
                        best = (d, kind)
                    break
    return best


def injectate_temperature(w):
    rec = w["records"].get("Injection_Enthalpy")
    if not rec:
        return float("nan"), float("nan")
    h = float(np.median(list(rec["series"].values())))
    return t_from_h_liq(h, 10.0), h


# ---- cross-checks between sources ------------------------------------------------------
DAT_NAME = re.compile(r"^(?P<pre>Wel+)(?P<n>\d+)_(?:(?P<elev>-?\d+(?:\.\d+)?)_)?(?P<kind>[a-z_]+)\.dat$")


def crosscheck_dat(W, datdir):
    """Compare every plain-column file with the JSON series of the same well (by name and by content)."""
    byn = {w["num"]: w for w in W}
    rows, refs = [], {}
    files = sorted(os.listdir(datdir))
    for f in files:
        m = DAT_NAME.match(f)
        row = dict(file=f, well_in_name="", content_well="", kind="", irregular_name="", n_dat=0, n_json=0,
                   max_abs_diff=float("nan"), max_axis_diff=float("nan"), ref_date="", sign="", status="")
        if not m:
            row.update(status="UNPARSED NAME", irregular_name="yes")
            rows.append(row)
            continue
        n, kind, pre = int(m.group("n")), m.group("kind"), m.group("pre")
        row.update(well_in_name=f"WELL-{n}", kind=kind, irregular_name=("yes: prefix '%s'" % pre) if pre != "Well" else "")
        a, b = read_dat2(os.path.join(datdir, f))
        row["n_dat"] = len(a)

        def series_for(w):
            """JSON series comparable to this file: (axis values, values, dates or None)."""
            if kind in ("downhole_temperature", "downhole_pressure"):
                v = "T" if kind.endswith("temperature") else "P"
                d = w["surveys"][w["survey_date"]][v]["data"]
                zs = sorted(d, reverse=True)
                return np.array(zs), np.array([d[z] for z in zs]), None
            if kind == "transient_pressure":
                if m.group("elev") is None:
                    return None
                zq = float(m.group("elev"))
                tp = w["transient_P"].get(zq)
                if tp is None:
                    return None
                ds = sorted(tp["series"])
                return None, np.array([tp["series"][d] for d in ds]), ds
            keymap = {"production_mass_flow": "Production_Mass_Flow", "injection_mass_flow": "Injection_Mass_Flow",
                      "transient_enthalpy": "Transient_Enthalpy", "injection_enthalpy": "Injection_Enthalpy"}
            rec = w["records"].get(keymap.get(kind, ""))
            if rec is None:
                return None
            ds = sorted(rec["series"])
            return None, np.array([rec["series"][d] for d in ds]), ds

        def compare(w):
            s = series_for(w)
            if s is None:
                return None
            ax, vals, ds = s
            if ds is None:                                        # survey: (value, elevation)
                if len(vals) != len(a):
                    return None
                o = np.argsort(-b)
                return dict(diff=float(np.max(np.abs(a[o] - vals))), axis=float(np.max(np.abs(b[o] - ax))), n=len(vals), sign="", ref="")
            if len(vals) != len(a):
                return None
            o = np.argsort(a)
            av, bv = a[o], b[o]
            sg = 1.0
            if np.max(np.abs(bv - vals)) > np.max(np.abs(-bv - vals)):
                sg = -1.0
            diff = float(np.max(np.abs(sg * bv - vals)))
            # reference date: first json date minus first time, then check all dates
            ref0 = ds[0] - dt.timedelta(days=round(av[0] * 365.25))
            cands = [ref0 + dt.timedelta(days=k) for k in range(-3, 4)]
            # reference date: among the candidates whose worst time error is within 0.003 yr of the best one,
            # take 1 January (the series are in years since a New-Year date), else the best match
            errs = {r: max(abs(years_since(d, r) - t) for d, t in zip(ds, av)) for r in cands}
            best = min(errs.values())
            ok_c = [r for r in cands if errs[r] <= best + 0.003]
            ref = next((r for r in ok_c if (r.month, r.day) == (1, 1)), min(ok_c, key=lambda r: errs[r]))
            axis = max(abs(years_since(d, ref) - t) for d, t in zip(ds, av))
            return dict(diff=diff, axis=float(axis), n=len(vals), sign=("dat value = -json" if sg < 0 else "same sign"),
                        ref=ref.isoformat())

        mine = compare(byn[n]) if n in byn else None
        if mine is None:
            row.update(status="NO MATCHING JSON SERIES")
        else:
            row.update(n_json=mine["n"], max_abs_diff=mine["diff"], max_axis_diff=mine["axis"], sign=mine["sign"],
                       ref_date=mine["ref"], content_well=f"WELL-{n}")
            tol_axis = 0.02 if kind not in ("downhole_temperature", "downhole_pressure") else 0.06
            ok = mine["diff"] < 1e-6 + 0.051 and mine["axis"] < tol_axis
            row["status"] = "OK" if ok else "MISMATCH"
        # content search among all wells to confirm the well
        hits = [x["name"] for x in W if compare(x) is not None and compare(x)["diff"] < 1e-6]
        if hits and row["content_well"] and f"WELL-{n}" not in hits:
            row["status"] = "MISMATCH (content matches " + ",".join(hits) + ")"
        row["content_matches"] = ",".join(hits)
        if mine and mine["ref"]:
            refs[f] = mine["ref"]
        rows.append(row)
    return pd.DataFrame(rows), refs


def integrity_checks(W, grid, wb, groups):
    """Track, survey, feed-zone, design and block-file consistency checks. One row per check/well."""
    rows = []

    def add(well, check, ok, detail, level=None):
        rows.append(dict(well=well, check=check, status="PASS" if ok else (level or "ISSUE"), detail=detail))

    gwells = sorted({x for a in groups.get("Areas", {}).values() for g in a.values() for x in g})
    add("ALL", "well list (groups file) equals the well data files", sorted(gwells) == sorted(w["name"] for w in W),
        f"{len(gwells)} wells in groups file, {len(W)} well JSON files")
    for w in W:
        n = w["name"]
        zt, zb = w["wellhead_z"], w["td_z"]
        sv = w["surveys"][w["survey_date"]]
        zs = sorted(set(sv["T"]["data"]) | set(sv["P"]["data"]), reverse=True)
        same = set(sv["T"]["data"]) == set(sv["P"]["data"])
        add(n, "T and P surveys share the same elevations", same, f"{len(sv['T']['data'])} T points, {len(sv['P']['data'])} P points")
        out = [z for z in zs if not (zb - 1e-6 <= z <= zt + 1e-6)]
        add(n, "every survey point lies within the well track", not out, f"outside track: {out}" if out else f"{zs[0]:.1f} to {zs[-1]:.1f} mRL inside {zt:.1f} to {zb:.1f} mRL")
        add(n, "survey reaches the wellhead", abs(zs[0] - zt) < 1e-6,
            f"top survey point {zs[0]:.1f} mRL is {zt - zs[0]:.1f} m below the wellhead ({zt:.1f} mRL)", level="NOTE")
        add(n, "survey reaches track bottom", abs(zs[-1] - zb) < 1e-6, f"deepest survey {zs[-1]:.1f} vs track bottom {zb:.1f} mRL", level="NOTE")
        for f in w["fz"]:
            inside = in_interval(f["top"], zb, zt) and in_interval(f["bottom"], zb, zt)
            add(n, f"feed zone {f['n']} within the well track", inside, f"{f['top']:.0f} to {f['bottom']:.0f} mRL")
            slot = in_interval(f["top"], w["liner_bottom_z"], w["liner_top_z"]) and in_interval(f["bottom"], w["liner_bottom_z"], w["liner_top_z"])
            add(n, f"feed zone {f['n']} within the slotted/open interval", slot,
                f"{f['top']:.0f} to {f['bottom']:.0f} mRL vs slotted {w['liner_top_z']:.0f} to {w['liner_bottom_z']:.0f} mRL")
            add(n, f"feed zone {f['n']} quality tag is a quality grade", f["quality"].upper() in QUALITY_RANK_DEFAULT,
                f"tag '{f['quality']}'", level="ISSUE")
        add(n, "casing top at wellhead", abs(w["casing_top_z"] - zt) < 1e-6, f"casing top {w['casing_top_z']:.0f} vs wellhead {zt:.0f} mRL")
        add(n, "liner bottom at track bottom", abs(w["liner_bottom_z"] - zb) < 1.0,
            f"liner bottom {w['liner_bottom_z']:.1f} vs track bottom {zb:.1f} mRL ({w['liner_bottom_z'] - zb:+.1f} m)")
        add(n, "liner overlaps the casing (hanger above shoe)", w["liner_top_z"] >= w["shoe_z"],
            f"liner top {w['liner_top_z']:.0f}, casing shoe {w['shoe_z']:.0f} mRL")
        # grid-file well track
        key = "WEL%d" % w["num"] if w["num"] >= 10 else "WELL%d" % w["num"]
        gt = grid["wells"].get(key)
        if gt is None:
            add(n, "track equals the grid-file well track", False, "no such well in the grid file WELLS section")
        else:
            gt = np.array(gt)
            d = np.max(np.abs(gt - w["track"])) if gt.shape == w["track"].shape else float("inf")
            add(n, "track equals the grid-file well track", d < 0.1, f"max coordinate difference {d:.2f} m ({len(gt)} points)")
        # block file
        bl = wb.get(n)
        if bl is None:
            add(n, "well is in the well-to-block file", False, "missing")
            continue
        ivs = [b["interval"] for b in bl["Blocks"]]
        add(n, "block file covers the track without gaps",
            abs(max(ivs[0]) - zt) < 1e-6 and abs(min(ivs[-1]) - zb) < 0.5 and all(abs(ivs[i][1] - ivs[i + 1][0]) < 0.5 for i in range(len(ivs) - 1)),
            f"{len(ivs)} blocks, {max(ivs[0]):.1f} to {min(ivs[-1]):.1f} mRL")
        bfz = list(bl["Feed_Zones"].values())
        okfz = len(bfz) == len(w["fz"])
        if okfz:
            for f, bf in zip(w["fz"], bfz):
                ints = [v["Interval"] for v in bf.values()]
                okfz = okfz and all(abs(i[0] - f["top"]) < 0.01 and abs(i[1] - f["bottom"]) < 0.01 for i in ints)
        add(n, "feed zones equal the block-file feed zones (count, order, interval)", okfz, f"{len(w['fz'])} in JSON, {len(bfz)} in block file")
        mism = []
        for z in zs:
            x, y = track_xy_at(w, z)
            c, l = column_at(grid, x, y), layer_at(grid, z)
            b = block_at(wb, n, z)
            if b is None or split_block(b) != (c, l):
                mism.append(z)
        add(n, "own point-in-polygon and layer search agrees with the block file", not mism,
            f"{len(zs) - len(mism)} of {len(zs)} survey points agree" + (f"; disagreeing: {mism}" if mism else ""))
    return pd.DataFrame(rows)


def provenance_table(W, groups):
    rows = []
    for a, g in groups.get("Areas", {}).items():
        for gname, wells in g.items():
            rows.append(dict(well="(group)", field=f"Area '{a}' / group '{gname}'", label="; ".join(wells)))
    for w in W:
        n = w["name"]
        rows.append(dict(well=n, field="Well_Design.Source", label=w["design_source"]))
        rows.append(dict(well=n, field="Feed_Zones.Units (Top/Bottom/Quality/Source)",
                         label="'m' (values are elevations in mRL) / 'Text' / 'Text' (no Source value given)"))
        for date, sv in w["surveys"].items():
            for v, s in sv.items():
                rows.append(dict(well=n, field=f"Downhole_{'Temperature' if v == 'T' else 'Pressure'}[{date}]",
                                 label=f"Source='{s['source']}'; Comment='{s['comment']}'; StableFlag={s['stable']}"))
        for k, rec in w["records"].items():
            rows.append(dict(well=n, field=k, label=f"Source='{rec['source']}'; Comment='{rec['comment']}'"))
        for z, tp in w["transient_P"].items():
            rows.append(dict(well=n, field=f"Transient_Pressure[{z:g}]", label=f"Source='{tp['source']}'; Comment='{tp['comment']}'"))
        for f in w["fz"]:
            rows.append(dict(well=n, field=f"Feed zone {f['n']} quality tag", label=f["quality"]))
    return pd.DataFrame(rows)


# ---- natural-state timing ----------------------------------------------------------------
def natural_state_assessment(W, profs, summ, qmap, isos):
    """Gauge test, trend test and injectate test (CLAUDE.md 4.1). Returns (rows, trend rows)."""
    flows = {w["name"]: first_flow(w) for w in W}
    field_first = min(v[0] for v in flows.values() if v)
    info = {}
    for w in W:
        n, prof = w["name"], profs[w["name"]]
        d = dict(well=n, role="none", well_first_flow=None, gauge_elev=float("nan"), survey_P_at_gauge=float("nan"),
                 first_gauge_date=None, first_gauge_P=float("nan"), diff=float("nan"))
        if flows[n]:
            d["well_first_flow"], d["role"] = flows[n]
        if w["transient_P"]:
            zq = sorted(w["transient_P"])[0]
            ser = w["transient_P"][zq]["series"]
            d0 = min(ser)
            d.update(gauge_elev=zq, first_gauge_date=d0, first_gauge_P=ser[d0])
            d["survey_P_at_gauge"] = float(interp_z(prof["z"], prof["P"], zq))
            d["diff"] = d["survey_P_at_gauge"] - ser[d0]
            ds = sorted(ser)
            sd = pdate(w["survey_date"])
            # gauge series interpolated to the survey date (linear in time)
            tt = [years_since(x, ds[0]) for x in ds]
            d["gauge_at_survey_date"] = float(np.interp(years_since(sd, ds[0]), tt, [ser[x] for x in ds])) if ds[0] <= sd <= ds[-1] else float("nan")
            d["survey_minus_gauge_at_survey_date"] = d["survey_P_at_gauge"] - d["gauge_at_survey_date"]
            d["days_gauge_start_after_well_first_flow"] = (d0 - d["well_first_flow"]).days if d["well_first_flow"] else None
            d["days_gauge_start_after_field_first_flow"] = (d0 - field_first).days
        info[n] = d
    tol = PRM["gauge_match_bar"]
    passed = [n for n, d in info.items() if d["first_gauge_date"] and abs(d["diff"]) <= tol
              and d["first_gauge_date"] <= field_first]
    # trend of passed wells at the feed zones of every well
    trend_rows = []
    for w in W:
        n = w["name"]
        for f in w["fz"]:
            p = float(interp_z(profs[n]["z"], profs[n]["P"], f["mid"]))
            vals = []
            for q in passed:
                if q == n:
                    continue
                v = float(interp_z(profs[q]["z"], profs[q]["P"], f["mid"]))
                if not np.isnan(v):
                    vals.append(v)
            tr = float(np.mean(vals)) if vals else float("nan")
            trend_rows.append(dict(well=n, feedzone=f["n"], mid_mRL=f["mid"], quality=f["quality"], P_bar=p,
                                   n_trend_wells=len(vals), trend_mean_bar=tr,
                                   trend_min_bar=float(np.min(vals)) if vals else float("nan"),
                                   trend_max_bar=float(np.max(vals)) if vals else float("nan"),
                                   residual_bar=(p - tr) if vals else float("nan")))
    trend = pd.DataFrame(trend_rows)
    out = []
    for w in W:
        n = w["name"]
        d = info[n]
        t = trend[(trend.well == n) & trend.residual_bar.notna()]
        maxres = float(t.residual_bar.abs().max()) if len(t) else float("nan")
        meanres = float(t.residual_bar.mean()) if len(t) else float("nan")
        Tmax = summ[n]["Tmax_degC"]
        Tinj, hinj = injectate_temperature(w)
        reasons = []
        verdict, conf = None, None
        if d["first_gauge_date"]:
            pre = d["first_gauge_date"] <= field_first
            reasons.append(
                f"gauge at {d['gauge_elev']:.0f} mRL starts {d['first_gauge_date'].isoformat()} "
                f"({d['days_gauge_start_after_field_first_flow']:+d} d vs first production in the field, "
                f"{d['days_gauge_start_after_well_first_flow']:+d} d vs this well's first flow); "
                f"survey - first gauge = {d['diff']:+.2f} bar")
            if abs(d["diff"]) <= tol and pre:
                verdict, conf = "natural state", "high"
                reasons.append("survey equals the first (pre-/start-of-production) gauge reading within the reading precision")
            elif d["diff"] > tol and pre:
                verdict, conf = "probably natural state", "medium"
                reasons.append("survey is ABOVE the first gauge reading, so it is not production-drawn-down; the first gauge value already sits "
                               "at the start of production at this well's flow rate, so an exact match cannot be demanded")
            elif d["diff"] < -tol and pre:
                verdict, conf = "possibly disturbed", "low"
                reasons.append("survey is BELOW the first gauge reading (drawdown before the survey)")
            else:
                reasons.append("gauge series starts after production began in the field, so it cannot serve as a natural-state reference; "
                               "decided from the feed-zone trend test")
        if verdict is None:
            if np.isnan(maxres):
                verdict, conf = "possibly disturbed", "low"
                reasons.append("no feed zone can be compared with the trend wells")
            elif maxres <= PRM["fz_trend_ok_bar"]:
                verdict, conf = "probably natural state", "medium"
                reasons.append(f"feed-zone pressure within {maxres:.1f} bar of the trend of the wells that passed the gauge test (limit {PRM['fz_trend_ok_bar']:.0f} bar)")
            else:
                verdict, conf = "possibly disturbed", "low"
                reasons.append(f"feed-zone pressure differs from the trend of the passed wells by {maxres:.1f} bar (limit {PRM['fz_trend_ok_bar']:.0f} bar)")
        if not np.isnan(maxres) and d["first_gauge_date"] and maxres > PRM["fz_trend_ok_bar"] and verdict != "possibly disturbed":
            reasons.append(f"(feed-zone trend residual {maxres:.1f} bar)")
        # injectors: isothermal section near the injectate temperature?
        if d["role"] == "injection":
            reasons.append(f"injectate {hinj:.0f} kJ/kg (median of the record) = {Tinj:.1f} degC (IF97 Region 1 liquid at 10 bar)")
            Tfz = [float(interp_z(profs[n]["z"], profs[n]["T"], f["mid"])) for f in w["fz"]]
            has_iso_near = any(abs(a - Tinj) <= 5 and b - a <= PRM["iso_range_degC"] for a, b in
                               [(s_[2], s_[3]) for s_ in isos[n]])
            if max(Tfz) < Tinj - 5 and not has_iso_near:
                reasons.append(f"feed-zone temperature ({', '.join(f'{t:.1f}' for t in Tfz)} degC) is colder than the injectate and the profile has no "
                               f"isothermal section near {Tinj:.0f} degC (Tmax {Tmax:.1f} degC): evidence against an injection-affected survey")
            else:
                if conf == "medium":
                    conf = "low"
                reasons.append(f"formation at the feed zone ({', '.join(f'{t:.1f}' for t in Tfz)} degC) is hotter than the injectate, so injection "
                               f"cooling could be hidden: confidence lowered")
        elif d["role"] == "none":
            reasons.append("no production or injection record for this well")
        if d["role"] == "production":
            reasons.append("producer: record starts " + d["well_first_flow"].isoformat())
        d.update(verdict=verdict, confidence=conf, reason="; ".join(reasons), trend_max_abs_resid=maxres,
                 trend_mean_resid=meanres, Tinj=Tinj, hinj=hinj)
        out.append(d)
    return out, trend, field_first, passed


# ---- tables -------------------------------------------------------------------------------
def pick_fz_block(wb, well, k, zmid):
    """Block of feed zone k (1-based, JSON order) from the well-to-block file: the block holding the
    mid-point, else the block with the greatest feed length."""
    fzs = list(wb[well]["Feed_Zones"].values())
    d = fzs[k - 1]
    inside = [b for b, _, _ in block_candidates(wb, well, zmid, tol=0.0) if b in d]
    if inside:
        return inside[0], sorted(d)
    best = max(d.items(), key=lambda kv: kv[1].get("Length", 0))[0]
    return best, sorted(d)


def location_weight(r):
    """Temperature location weight and reason (CLAUDE.md section 5 defaults)."""
    if r["above_model_surface"]:
        return W_ABOVE, "above the model surface"
    if r["near_top"]:
        return W_TOP, "within 30 m of the wellhead or in the top block of its column"
    if r["in_feedzone"]:
        return W_FZ, "feed zone"
    if r["behind_casing"]:
        return W_CASED, "behind casing"
    if r["in_iso"]:
        return W_ISO, "open hole within an isothermal section >= 300 m long (possible internal flow)"
    return W_OPEN, "open hole"


def build_targets(w, prof, summ, assess, iso, grid, wb, qmap):
    n = w["name"]
    cf = CONF_FACTOR[assess["confidence"]]
    rows = []
    for k, z in enumerate(prof["z"]):
        x, y = track_xy_at(w, z)
        blk = block_at(wb, n, z)
        col, lay = split_block(blk)
        if col is None:
            col, lay = column_at(grid, x, y), layer_at(grid, z)
            blk = f"{col}{lay}"
        surf = grid["surf"].get(col, float("nan"))
        top_layer = layer_at(grid, surf)
        fzs = [f for f in w["fz"] if in_interval(z, f["top"], f["bottom"])]
        r = dict(well=n, x=round(x, 1), y=round(y, 1), elevation_mRL=float(z), T_degC=float(prof["T"][k]),
                 P_bar=float(prof["P"][k]), Tsat_degC=float(prof["Tsat"][k]),
                 boiling_margin_degC=float(prof["margin"][k]), column=col, layer=lay, block=blk,
                 in_feedzone=bool(fzs), behind_casing=bool(z >= w["shoe_z"] - 1e-9),
                 above_model_surface=bool(z > surf + 1e-9))
        r["near_top"] = bool(z >= w["wellhead_z"] - PRM["top_block_m"] - 1e-9 or lay == top_layer)
        r["in_iso"] = any(a + 1e-9 >= z >= b - 1e-9 for a, b, _, _ in iso)
        lw, why = location_weight(r)
        r["ns_confidence"] = assess["confidence"]
        r["weight"] = round(lw * cf, 4)
        r["weight_reason"] = f"{why} ({lw:g}) x natural-state confidence {assess['confidence']} ({cf:g})"
        r["weight_P"] = round(max(qmap[f["quality"]] for f in fzs) * cf, 4) if fzs else 0.0
        r["feedzone"] = ",".join(str(f["n"]) for f in fzs)
        r.pop("near_top"), r.pop("in_iso")
        rows.append(r)
    return rows


def build_feedzones(w, prof, assess, wb, qmap):
    n = w["name"]
    cf = CONF_FACTOR[assess["confidence"]]
    rows = []
    for f in w["fz"]:
        T = float(interp_z(prof["z"], prof["T"], f["mid"]))
        P = float(interp_z(prof["z"], prof["P"], f["mid"]))
        ts = tsat_c(P) if not np.isnan(P) else float("nan")
        blk, allb = pick_fz_block(wb, n, f["n"], f["mid"])
        rows.append(dict(well=n, feedzone=f["n"], top_mRL=f["top"], bottom_mRL=f["bottom"], mid_mRL=f["mid"],
                         quality=f["quality"], T_degC=T, P_bar=P, Tsat_degC=ts, boiling_margin_degC=ts - T,
                         block=blk, ns_confidence=assess["confidence"], weight_T=round(W_FZ * cf, 4),
                         weight_P=round(qmap[f["quality"]] * cf, 4), quality_factor=qmap[f["quality"]],
                         blocks_in_block_file=" ".join(allb), length_m=f["length"]))
    return rows


def build_layer_centres(w, prof, grid, wb):
    n = w["name"]
    zmin, zmax = float(prof["z"].min()), float(prof["z"].max())
    rows = []
    for i, L in grid["layers"].items():
        c = L["centre"]
        if not (zmin - 1e-9 <= c <= zmax + 1e-9):
            continue
        blk = block_at(wb, n, c)
        if blk is None:
            continue
        T = float(interp_z(prof["z"], prof["T"], c))
        P = float(interp_z(prof["z"], prof["P"], c))
        ts = tsat_c(P)
        rows.append(dict(well=n, layer=i, layer_centre_mRL=c, block=blk, T_degC=T, P_bar=P, Tsat_degC=ts,
                         boiling_margin_degC=ts - T))
    return rows


COMMON_ELEV = [0.0, -250.0, -500.0, -750.0, -1000.0, -1500.0, -2000.0, -2500.0]


def common_elevation_tables(W, profs):
    rows = []
    for w in W:
        p = profs[w["name"]]
        for z in COMMON_ELEV:
            T = float(interp_z(p["z"], p["T"], z))
            P = float(interp_z(p["z"], p["P"], z))
            if np.isnan(T):
                continue
            ts = tsat_c(P)
            x, y = track_xy_at(w, z)
            rows.append(dict(well=w["name"], elevation_mRL=z, x=round(x, 1), y=round(y, 1), T_degC=T, P_bar=P,
                             Tsat_degC=ts, boiling_margin_degC=ts - T))
    long = pd.DataFrame(rows)
    wide_T = long.pivot(index="well", columns="elevation_mRL", values="T_degC").reindex([w["name"] for w in W])
    wide_P = long.pivot(index="well", columns="elevation_mRL", values="P_bar").reindex([w["name"] for w in W])
    out = []
    for wd in (wide_T, wide_P):
        wd.columns = [f"z_{int(c)}_mRL" for c in wd.columns]
        out.append(wd.reset_index())
    return long, out[0], out[1]


def geometry_table(W, grid, feats):
    faults = [k for k, v in feats.items() if str(v["meta"].get("legend", "")).lower().startswith("fault") or k.lower().startswith("fault")]
    waters = [k for k in feats if k not in faults]
    rows = []
    for w in W:
        x0, y0 = w["wellhead_xy"]
        col = column_at(grid, x0, y0)
        surf = grid["surf"].get(col, float("nan"))
        dev = float(np.max(np.hypot(w["track"][:, 0] - x0, w["track"][:, 1] - y0)))
        d = {k: dist_to_feature(x0, y0, feats[k]) for k in feats}
        nf = min(faults, key=lambda k: d[k]) if faults else ""
        nw = min(waters, key=lambda k: d[k]) if waters else ""
        r = dict(well=w["name"], wellhead_mRL=w["wellhead_z"], model_surface_mRL=surf,
                 wellhead_minus_surface_m=w["wellhead_z"] - surf, max_horizontal_deviation_m=dev,
                 nearest_fault=nf, dist_nearest_fault_m=d.get(nf, float("nan")),
                 nearest_water_body=nw, dist_nearest_water_body_m=d.get(nw, float("nan")))
        for k in feats:
            r[f"dist_{feat_key(k)}_m"] = d[k]
        r["model_column"] = col
        r["x_wellhead"], r["y_wellhead"] = x0, y0
        rows.append(r)
    return pd.DataFrame(rows), faults, waters


def hydro_tables(W, profs, summ):
    rows, notes = [], []
    for w in W:
        wins, nt = hydro_windows(profs[w["name"]], w)
        notes += nt
        for a, b, band in wins:
            s = window_stats(profs[w["name"]], a, b)
            rows.append(dict(well=w["name"], z_top_mRL=a, z_bottom_mRL=b, length_m=s["length"], band=band,
                             dP_measured_bar=s["dP_measured"], dP_IF97_bar=s["dP_IF97"],
                             ratio_measured_to_IF97=s["ratio"], rho_implied_kgm3=s["rho_implied"],
                             rounding_unc_ratio=s["rounding_ratio"], rounding_unc_rho_kgm3=s["rounding_rho"],
                             mean_T_degC=s["mean_T"], casing_shoe_mRL=w["shoe_z"],
                             cold_well=bool(summ[w["name"]]["Tmax_degC"] < PRM["cold_well_Tmax_degC"])))
    hy = pd.DataFrame(rows)
    # summaries per well and band, pooled, and the cold-well (datum/convention) component
    srows = []
    comp = {}
    for band in ("cased", "open"):
        c = hy[(hy.band == band) & hy.cold_well]
        comp[band] = float(c.ratio_measured_to_IF97.median()) if len(c) else float("nan")
    for band in ("cased", "open"):
        for scope, sub in [("all wells", hy[hy.band == band]), ("cold wells (Tmax<100 degC)", hy[(hy.band == band) & hy.cold_well]),
                           ("warm/hot wells", hy[(hy.band == band) & ~hy.cold_well])]:
            if len(sub):
                srows.append(dict(scope=scope, band=band, n_windows=len(sub), length_m=sub.length_m.sum(),
                                  ratio_median=sub.ratio_measured_to_IF97.median(), ratio_min=sub.ratio_measured_to_IF97.min(),
                                  ratio_max=sub.ratio_measured_to_IF97.max(),
                                  ratio_length_weighted=sub.dP_measured_bar.sum() / sub.dP_IF97_bar.sum(),
                                  cold_component=comp[band], excess_over_cold_component=float("nan")))
        for w in W:
            sub = hy[(hy.well == w["name"]) & (hy.band == band)]
            if len(sub):
                rw = sub.dP_measured_bar.sum() / sub.dP_IF97_bar.sum()
                srows.append(dict(scope=w["name"], band=band, n_windows=len(sub), length_m=sub.length_m.sum(),
                                  ratio_median=sub.ratio_measured_to_IF97.median(), ratio_min=sub.ratio_measured_to_IF97.min(),
                                  ratio_max=sub.ratio_measured_to_IF97.max(), ratio_length_weighted=rw,
                                  cold_component=comp[band], excess_over_cold_component=rw - comp[band]))
    hs = pd.DataFrame(srows)
    # integer column last: the acceptance check reads "8,5861.0" as one number
    hs = hs[[c for c in hs.columns if c != "n_windows"] + ["n_windows"]]
    return hy, hs, notes, comp


# ======================================================================================
# 7a. figures
# ======================================================================================
WELL_COLORS = plt.get_cmap("tab10")
FAULT_KW = dict(color="k", lw=1.6, zorder=3)
WATER_C = "#4aa3df"


def short(n):
    return n.replace("WELL-", "W")


def draw_features(ax, feats, labels=True, xlim=None, ylim=None):
    """Faults, rivers and the lake outline from the mapping file."""
    for name, f in feats.items():
        is_fault = name.lower().startswith("fault")
        for s in f["segments"]:
            if f["type"] == "fill":
                ax.fill(s[:, 0], s[:, 1], color=WATER_C, alpha=0.5, zorder=1, lw=0)
                ax.plot(s[:, 0], s[:, 1], color=WATER_C, lw=0.8, zorder=1)
            elif is_fault:
                ax.plot(s[:, 0], s[:, 1], **FAULT_KW)
            else:
                ax.plot(s[:, 0], s[:, 1], color=WATER_C, lw=1.8, zorder=2)
        if labels and is_fault and xlim is not None:
            s = f["segments"][0]
            # label at the point of the fault inside the window that is closest to its centre
            ts = np.linspace(0, 1, 400)
            xs = np.interp(ts, [0, 1], [s[0, 0], s[-1, 0]])
            ys = np.interp(ts, [0, 1], [s[0, 1], s[-1, 1]])
            inside = (xs > xlim[0]) & (xs < xlim[1]) & (ys > ylim[0]) & (ys < ylim[1])
            if inside.any():
                k = np.where(inside)[0]
                vertical = abs(s[-1, 1] - s[0, 1]) > abs(s[-1, 0] - s[0, 0])
                # label near the window edge so it does not sit on the wells
                kk = k[int(len(k) * 0.04)] if vertical else k[int(len(k) * 0.93)]
                ax.text(xs[kk], ys[kk], " " + name, fontsize=8, fontweight="bold", zorder=6,
                        bbox=dict(fc="w", ec="none", alpha=0.7, pad=0.5))


def savefig(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name), dpi=130, bbox_inches="tight")
    plt.close(fig)


def fig_well_map(W, grid, feats, summ):
    fig, axs = plt.subplots(1, 2, figsize=(15, 7.2), gridspec_kw=dict(width_ratios=[1, 1.15]))
    xs = np.concatenate([w["track"][:, 0] for w in W])
    ys = np.concatenate([w["track"][:, 1] for w in W])
    for k, ax in enumerate(axs):
        if k == 0:
            xlim, ylim = (2809000, 2830000), (6326000, 6343500)
        else:
            xlim, ylim = (xs.min() - 1500, xs.max() + 1500), (ys.min() - 1500, ys.max() + 1500)
        for c, pl in grid["polys"].items():
            p = np.array(pl + [pl[0]])
            ax.plot(p[:, 0], p[:, 1], color="0.82", lw=0.4, zorder=0)
        draw_features(ax, feats, xlim=xlim, ylim=ylim)
        for i, w in enumerate(W):
            tr = w["track"]
            if np.ptp(tr[:, 0]) > 1 or np.ptp(tr[:, 1]) > 1:
                ax.plot(tr[:, 0], tr[:, 1], color="C3", lw=1.8, zorder=4)
                ax.plot(tr[-1, 0], tr[-1, 1], "x", color="C3", ms=6 if k else 4, zorder=5)
            ax.plot(w["wellhead_xy"][0], w["wellhead_xy"][1], "o", color="k", mfc=WELL_COLORS(i), ms=9 if k else 5, mew=1.2, zorder=6)
            if k == 1 or (xlim[0] < w["wellhead_xy"][0] < xlim[1]):
                ax.annotate(short(w["name"]), w["wellhead_xy"], xytext=(7, 6), textcoords="offset points", fontsize=11 if k else 8,
                            fontweight="bold", zorder=7)
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal")
        ax.set_xlabel("Easting (m)")
        ax.set_ylabel("Northing (m)")
        ax.ticklabel_format(style="plain")
        ax.tick_params(axis="x", labelrotation=30, labelsize=8)
        ax.tick_params(axis="y", labelsize=8)
    axs[0].set_title("Model area: grid columns (grey), faults (black), lake and rivers (blue)")
    axs[1].set_title("Wells: wellhead (circle), deviated track and bottom hole (red line, x)")
    axs[0].add_patch(plt.Rectangle((xs.min() - 1500, ys.min() - 1500), np.ptp(xs) + 3000, np.ptp(ys) + 3000, fill=False, ec="C1", lw=1.5, zorder=5))
    savefig(fig, "well_map.png")


def fig_profile(w, prof, r, assess, iso, qmap):
    z, T, P = prof["z"], prof["T"], prof["P"]
    fig, axs = plt.subplots(1, 2, figsize=(11, 8), sharey=True)
    ztop = max(w["wellhead_z"], z[0], r["water_level_mRL"] if not is_nan(r["water_level_mRL"]) else -1e9) + 60
    zbot = min(w["td_z"], z[-1]) - 40
    ax = axs[0]
    ax.plot(T, z, "o-", color="C3", ms=4, label="T survey")
    ts = prof["Tsat"]
    ax.plot(ts, z, "--", color="C0", label="Tsat(P) (IF97 Region 4)")
    ax.set_xlabel("Temperature (degC)")
    ax.set_ylabel("Elevation (mRL)")
    ax2 = axs[1]
    ax2.plot(P, z, "o-", color="C2", ms=4, label="P survey")
    ax2.set_xlabel("Pressure (bar)")
    ax2.set_xlim(right=float(np.nanmax(P)) * 1.15)
    for a in axs:
        a.set_ylim(zbot, ztop)
        a.grid(alpha=0.3)
        for f in w["fz"]:
            a.axhspan(f["bottom"], f["top"], color="orange", alpha=0.35, lw=0)
        a.axhline(w["shoe_z"], color="k", ls="--", lw=1)
        a.axhline(w["liner_top_z"], color="k", ls=":", lw=1)
        a.axhline(w["wellhead_z"], color="0.4", ls="-", lw=0.8)
        for zt, zb, _, _ in iso:
            a.axhspan(zb, zt, color="0.6", alpha=0.18, lw=0)
    x0 = axs[0].get_xlim()[0]
    axs[0].text(axs[0].get_xlim()[1], w["shoe_z"], " casing shoe", ha="right", va="bottom", fontsize=8)
    axs[0].text(axs[0].get_xlim()[1], w["liner_top_z"], " liner top", ha="right", va="bottom", fontsize=8)
    axs[0].text(axs[0].get_xlim()[1], w["wellhead_z"], " wellhead", ha="right", va="bottom", fontsize=8)
    for f in w["fz"]:
        axs[0].text(axs[0].get_xlim()[0], f["mid"], f" FZ{f['n']} {f['quality']}", va="center", fontsize=8, color="darkorange", fontweight="bold")
        pf = float(interp_z(z, P, f["mid"]))
        ax2.plot(pf, f["mid"], "s", color="darkorange", mec="k", ms=8, zorder=5)
        ax2.annotate(f"{pf:.1f} bar", (pf, f["mid"]), xytext=(8, 4), textcoords="offset points", fontsize=8)
    if not is_nan(r["pivot_elev_mRL"]):
        ax2.plot(r["pivot_P_bar"], r["pivot_elev_mRL"], "*", color="r", ms=15, zorder=6, label=f"pivot (FZ{r['pivot_feedzone']})")
    if not is_nan(r["water_level_mRL"]):
        for a in axs:
            a.axhline(r["water_level_mRL"], color="C0", lw=1, ls="-.")
        axs[1].text(axs[1].get_xlim()[1], r["water_level_mRL"], " water level (extrapolated)" if r["water_level_extrapolation_m"] > 0 else " water level",
                    ha="right", va="bottom", fontsize=8, color="C0")
    axs[0].legend(loc="lower left", fontsize=8)
    ax2.legend(loc="lower left", fontsize=8)
    fig.suptitle(f"{w['name']}: survey {w['survey_date']}, regime {r['regime']}, natural-state verdict: {assess['verdict']} ({assess['confidence']})\n"
                 f"orange = feed zones, dashed = casing shoe, dotted = liner top, grey band = isothermal section >= 300 m", fontsize=10)
    savefig(fig, f"profile_{w['name']}.png")


def tri_field(ax, xs, ys, vs, cmap, label, levels=None):
    """Filled contours by linear interpolation inside the hull of the wells only (no extrapolation)."""
    xs, ys, vs = map(np.asarray, (xs, ys, vs))
    if len(xs) < 4:
        return None
    try:
        tri = mtri.Triangulation(xs, ys)
        # drop long thin triangles on the hull edge
        cs = ax.tricontourf(tri, vs, levels=levels or 12, cmap=cmap, alpha=0.75, zorder=1)
        ax.tricontour(tri, vs, levels=levels or 12, colors="k", linewidths=0.4, alpha=0.5, zorder=2)
        return cs
    except Exception:
        return None


def fig_maps(W, long, feats, var, name_prefix, cmap, unit, label):
    xs_all = np.array([w["wellhead_xy"][0] for w in W])
    ys_all = np.array([w["wellhead_xy"][1] for w in W])
    xlim = (xs_all.min() - 1500, xs_all.max() + 1500)
    ylim = (ys_all.min() - 1500, ys_all.max() + 1500)
    for z in [0.0, -500.0, -1000.0, -1500.0]:
        sub = long[long.elevation_mRL == z]
        fig, ax = plt.subplots(figsize=(8.5, 7.4))
        draw_features(ax, feats, xlim=xlim, ylim=ylim)
        cs = tri_field(ax, sub.x, sub.y, sub[var], cmap, label)
        have = set(sub.well)
        for i, w in enumerate(W):
            if w["name"] in have:
                row = sub[sub.well == w["name"]].iloc[0]
                sc = ax.scatter(row.x, row.y, c=[row[var]], cmap=cmap, vmin=sub[var].min(), vmax=sub[var].max(), s=160, ec="k", lw=1.2, zorder=5)
                ax.annotate(f"{short(w['name'])}\n{row[var]:.1f}", (row.x, row.y), xytext=(9, 5), textcoords="offset points", fontsize=9,
                            fontweight="bold", zorder=7, bbox=dict(fc="w", ec="none", alpha=0.7, pad=0.5))
            else:
                x, y = w["wellhead_xy"]
                ax.plot(x, y, "o", mfc="none", mec="0.4", ms=9, zorder=5)
                ax.annotate(f"{short(w['name'])}\nno data", (x, y), xytext=(9, 5), textcoords="offset points", fontsize=8, color="0.4", zorder=7)
        if cs is not None:
            cb = fig.colorbar(cs, ax=ax, shrink=0.8)
            cb.set_label(f"{label} ({unit})")
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
        ax.set_aspect("equal")
        ax.ticklabel_format(style="plain")
        ax.tick_params(axis="x", labelrotation=30, labelsize=8)
        ax.set_xlabel("Easting (m)")
        ax.set_ylabel("Northing (m)")
        ax.set_title(f"{label} at {z:.0f} mRL ({len(sub)} wells with data; interpolation inside the well hull only)")
        savefig(fig, f"{name_prefix}_{'minus' if z < 0 else ''}{abs(int(z))}m.png")


def section_geometry(W, names):
    """Line through the given wells (principal axis) and the distance of each well along it."""
    pts = np.array([W[[w["name"] for w in W].index(n)]["wellhead_xy"] for n in names])
    c = pts.mean(axis=0)
    u, s, vt = np.linalg.svd(pts - c)
    d = vt[0]
    if d[1] < 0:        # run south to north? keep west-to-east or south-to-north orientation
        d = -d
    along = (pts - c) @ d
    off = (pts - c) @ np.array([-d[1], d[0]])
    o = np.argsort(along)
    return [names[i] for i in o], along[o] - along[o].min(), off[o], d


def fig_section(W, profs, summ, names, var, fname, cmap, label, unit, levels):
    names, dist, off, d = section_geometry(W, names)
    az = (math.degrees(math.atan2(d[0], d[1])) + 360) % 360
    fig, ax = plt.subplots(figsize=(12, 7.5))
    zgrid = np.arange(-3000, 300, 10.0)
    F = np.full((len(zgrid), len(dist)), np.nan)
    for j, n in enumerate(names):
        p = profs[n]
        F[:, j] = interp_z(p["z"], p[var], zgrid)
    # interpolate between wells at each elevation only where both neighbours have data
    xg = np.linspace(dist.min(), dist.max(), 200)
    Z = np.full((len(zgrid), len(xg)), np.nan)
    for i in range(len(zgrid)):
        for j in range(len(dist) - 1):
            if not (np.isnan(F[i, j]) or np.isnan(F[i, j + 1])):
                m = (xg >= dist[j]) & (xg <= dist[j + 1])
                Z[i, m] = np.interp(xg[m], dist[[j, j + 1]], F[i, [j, j + 1]])
    cs = ax.contourf(xg, zgrid, np.ma.masked_invalid(Z), levels=levels, cmap=cmap, alpha=0.8)
    ax.contour(xg, zgrid, np.ma.masked_invalid(Z), levels=levels, colors="k", linewidths=0.4)
    cb = fig.colorbar(cs, ax=ax)
    cb.set_label(f"{label} ({unit})")
    w_by = {w["name"]: w for w in W}
    for j, n in enumerate(names):
        w = w_by[n]
        p = profs[n]
        ax.plot([dist[j]] * len(p["z"]), p["z"], "-", color="k", lw=1)
        ax.scatter([dist[j]] * len(p["z"]), p["z"], c=p[var], cmap=cmap, norm=cs.norm, s=14, ec="k", lw=0.3, zorder=5)
        for f in w["fz"]:
            ax.plot([dist[j] - 40, dist[j] + 40], [f["mid"]] * 2, color="darkorange", lw=4, zorder=6)
        ax.plot(dist[j], w["shoe_z"], "v", color="k", ms=7, zorder=7)
        ax.text(dist[j], w["wellhead_z"] + 60, f"{short(n)}\n(offset {off[j]:+.0f} m)", ha="center", fontsize=9, fontweight="bold")
        ax.plot(dist[j], w["wellhead_z"], "s", color="k", ms=5)
    ax.set_xlabel(f"Distance along the section line (m); the line runs at azimuth {az:.0f} deg (north = 0, east = 90), distance increases along that direction")
    ax.set_ylabel("Elevation (mRL)")
    ax.set_ylim(-3000, 450)
    ax.set_title(f"{label} section through {', '.join(short(n) for n in names)} (hottest wells); orange = feed zones, triangle = casing shoe;\n"
                 "interpolation between neighbouring wells only where both have data; offsets are distances of wells from the section line")
    ax.grid(alpha=0.3)
    savefig(fig, fname)


def fig_pt(W, profs, feats):
    fig, axs = plt.subplots(1, 2, figsize=(14, 6.5))
    Tsat_line = np.linspace(1, 373.9, 400)
    Psat_line = np.array([psat_bar(t) for t in Tsat_line])
    for k, ax in enumerate(axs):
        ax.plot(Tsat_line, Psat_line, "k-", lw=1.8, label="saturation curve (IF97 Region 4)")
        ax.plot(T_CRIT_C, P_CRIT_BAR, "k*", ms=12, label="critical point")
        for i, w in enumerate(W):
            p = profs[w["name"]]
            ax.plot(p["T"], p["P"], "-o", ms=3.5, lw=0.8, color=WELL_COLORS(i), label=short(w["name"]), alpha=0.9)
            for f in w["fz"]:
                T = float(interp_z(p["z"], p["T"], f["mid"]))
                P = float(interp_z(p["z"], p["P"], f["mid"]))
                ax.plot(T, P, "s", mfc=WELL_COLORS(i), mec="k", ms=8, zorder=5)
        if k == 0:
            ax.set_xlim(0, 400)
            ax.set_ylim(0, 300)
            ax.set_title("P-T state of every survey point (squares = feed zones)")
        else:
            ax.set_xlim(150, 285)
            ax.set_ylim(0, 130)
            ax.set_title("Zoom on the hot wells: distance to the saturation curve = boiling margin")
            w10 = profs["WELL-10"]
            ax.annotate("W10: within 5 degC of Tsat", (w10["T"][-4], w10["P"][-4]), xytext=(252, 12), textcoords="data", arrowprops=dict(arrowstyle="->"), fontsize=9)
        ax.set_xlabel("Temperature (degC)")
        ax.set_ylabel("Pressure (bar)")
        ax.grid(alpha=0.3)
    h, lab = axs[0].get_legend_handles_labels()
    fig.legend(h, lab, loc="center left", bbox_to_anchor=(0.99, 0.5), fontsize=9)
    savefig(fig, "PT_state.png")


# ======================================================================================
# 7b. derived tables, conceptual figure
# ======================================================================================
def hot_zone(prof, drop=None):
    """Elevation span around the temperature maximum where T >= Tmax - drop (interpolated)."""
    drop = PRM["hot_layer_drop_degC"] if drop is None else drop
    z, T = prof["z"], prof["T"]
    i = int(np.nanargmax(T))
    thr = T[i] - drop
    top = z[0]
    for k in range(i, 0, -1):
        if T[k - 1] < thr <= T[k]:
            top = z[k - 1] + (thr - T[k - 1]) / (T[k] - T[k - 1]) * (z[k] - z[k - 1])
            break
    bot = z[-1]
    for k in range(i, len(z) - 1):
        if T[k + 1] < thr <= T[k]:
            bot = z[k] + (T[k] - thr) / (T[k] - T[k + 1]) * (z[k + 1] - z[k])
            break
    return float(top), float(bot)


def augment_summary(S, W, grid, profs):
    """Extra columns after the required ones."""
    by = {w["name"]: w for w in W}
    ex = []
    for _, r in S.iterrows():
        w = by[r.well]
        col = column_at(grid, *w["wellhead_xy"])
        surf = grid["surf"][col]
        zt, zb = hot_zone(profs[r.well])
        d100 = w["wellhead_z"] - r.z_100C_mRL if not is_nan(r.z_100C_mRL) else float("nan")
        ex.append(dict(wellhead_mRL=w["wellhead_z"], model_surface_mRL=surf,
                       water_level_extrapolation_m=r.water_level_extrapolation_m,
                       water_level_below_model_surface_m=surf - r.water_level_mRL,
                       depth_100C_below_wellhead_m=d100,
                       hot_zone_top_mRL=zt, hot_zone_bottom_mRL=zb, hot_zone_thickness_m=zt - zb,
                       shallow_heat_flux_mW_m2_assumed_k=r.shallow_gradient_degC_per_km * PRM["k_rock_W_mK"],
                       deep_heat_flux_mW_m2_assumed_k=r.deep_gradient_degC_per_km * PRM["k_rock_W_mK"],
                       casing_shoe_mRL=w["shoe_z"], pivot_feedzone=r.pivot_feedzone, bpd_points=r.bpd_points,
                       deep_fit_top_mRL=r.deep_fit_top_mRL, shallow_fit_n=r.shallow_fit_n,
                       n_isothermal_sections=r.n_iso_sections))
    ex = pd.DataFrame(ex)
    drop = ["water_level_extrapolation_m", "pivot_feedzone", "bpd_points", "deep_fit_top_mRL", "shallow_fit_n", "n_iso_sections", "track_TD_mRL"]
    base = S.drop(columns=[c for c in drop if c in S.columns])
    return pd.concat([base.reset_index(drop=True), ex.reset_index(drop=True)], axis=1)


def iso_table(W, isos):
    rows = []
    for w in W:
        for a, b, t0, t1 in isos[w["name"]]:
            fz = [f["n"] for f in w["fz"] if a + 1e-9 >= f["mid"] >= b - 1e-9]
            rel_ = ("behind casing" if b >= w["shoe_z"] else "open hole" if a <= w["shoe_z"] else "spans the casing shoe")
            rows.append(dict(well=w["name"], z_top_mRL=a, z_bottom_mRL=b, length_m=a - b, T_min_degC=t0, T_max_degC=t1,
                             T_range_degC=t1 - t0, location=rel_, feedzones_inside=",".join(map(str, fz))))
    return pd.DataFrame(rows)


def boiling_flags(W, profs):
    rows = []
    for w in W:
        p = profs[w["name"]]
        m = p["margin"]
        flag = (m < PRM["bpd_margin_degC"]) & ~np.isnan(m)
        i = 0
        while i < len(flag):
            if flag[i]:
                j = i
                while j + 1 < len(flag) and flag[j + 1]:
                    j += 1
                k = i + int(np.argmin(m[i:j + 1]))
                rows.append(dict(well=w["name"], z_top_mRL=p["z"][i], z_bottom_mRL=p["z"][j], n_points=j - i + 1,
                                 min_margin_degC=m[k], min_margin_elev_mRL=p["z"][k], T_at_min_degC=p["T"][k], P_at_min_bar=p["P"][k]))
                i = j + 1
            else:
                i += 1
    return pd.DataFrame(rows, columns=["well", "z_top_mRL", "z_bottom_mRL", "n_points", "min_margin_degC",
                                       "min_margin_elev_mRL", "T_at_min_degC", "P_at_min_bar"])


def reference_model_summary(inv):
    """Short description of the reference model (labelled [reference model] in the report; never evidence)."""
    out = {}
    try:
        bc = json.load(open(role_paths(inv, "reference_bc")[0], encoding="utf-8"))
        out["heat_flux_default_W_m2"] = bc["Heat"]["Default"]["Heat_Flux"]
        out["atmosphere_T_degC"] = bc["Atmosphere"]["Default"]["Temperature"]
        out["rain_default_mm_yr"] = bc["Rain"]["Default"]["Rainfall"]
    except Exception as e:                                   # pragma: no cover
        out["error_bc"] = str(e)
    try:
        ns = [p for p in role_paths(inv, "reference_rock") if p.lower().endswith(".dat")][0]
        txt = open(ns, errors="ignore").read().splitlines()
        i = txt.index("GENER") if "GENER" in txt else None
        heat = []
        if i is not None:
            for ln in txt[i + 1:]:
                if re.match(r"^[A-Z]{4,}\s*$", ln):
                    break
                p = ln.split()
                if len(p) >= 4 and p[2] == "HEAT":
                    heat.append(float(p[3]))
        out["n_heat_sources"] = len(heat)
        out["heat_total_MW"] = sum(heat) / 1e6
        out["heat_per_block_min_kW"] = min(heat) / 1e3 if heat else float("nan")
        out["heat_per_block_max_kW"] = max(heat) / 1e3 if heat else float("nan")
    except Exception as e:                                   # pragma: no cover
        out["error_ns"] = str(e)
    return out


def fig_conceptual(W, profs, summ, S, grid, feats, assess):
    """Evidence sketch: plan view (left) and a south-west to north-east schematic section (right)."""
    by = {w["name"]: w for w in W}
    hot = [n for n in summ if summ[n]["Tmax_degC"] >= PRM["hot_core_Tmax_degC"]]
    cxy = np.mean([by[n]["wellhead_xy"] for n in hot], axis=0)
    fig, axs = plt.subplots(1, 2, figsize=(17, 7.8), gridspec_kw=dict(width_ratios=[1, 1.25]))
    ax = axs[0]
    xs = [w["wellhead_xy"][0] for w in W]
    ys = [w["wellhead_xy"][1] for w in W]
    xlim, ylim = (min(xs) - 1400, max(xs) + 1400), (min(ys) - 1300, max(ys) + 1300)
    draw_features(ax, feats, xlim=xlim, ylim=ylim)
    for c, pl in grid["polys"].items():
        p = np.array(pl + [pl[0]])
        ax.plot(p[:, 0], p[:, 1], color="0.9", lw=0.3, zorder=0)
    from matplotlib.patches import Ellipse, FancyArrowPatch
    # inferred hot core (wells with Tmax >= 255 degC)
    hx = np.array([by[n]["wellhead_xy"] for n in hot])
    ax.add_patch(Ellipse(hx.mean(axis=0), np.ptp(hx[:, 0]) + 1300, np.ptp(hx[:, 1]) + 1900, angle=0, fc="red", alpha=0.18, ec="r", lw=2, zorder=2))
    ax.text(0.02, 0.62, f"inferred upflow / hot core\n(Tmax >= {PRM['hot_core_Tmax_degC']:.0f} degC, shallow)", color="darkred", fontsize=10, fontweight="bold",
            ha="left", transform=ax.transAxes)
    cmap = plt.get_cmap("inferno")
    for n, r in summ.items():
        x, y = by[n]["wellhead_xy"]
        col = {"two-phase": "magenta", "outflow": "darkorange", "recharge/downflow": "dodgerblue", "conductive": "0.5", "upflow": "red"}[r["regime"]]
        ax.scatter(x, y, s=260, c=col, ec="k", lw=1.5, zorder=6)
        ax.annotate(f"{short(n)}  Tmax {r['Tmax_degC']:.0f}", (x, y), xytext=(10, -14), textcoords="offset points", fontsize=9, fontweight="bold",
                    bbox=dict(fc="w", ec="none", alpha=0.75, pad=0.5), zorder=8)
        if r["regime"] == "outflow":
            dx, dy = x - cxy[0], y - cxy[1]
            L = math.hypot(dx, dy)
            ax.add_patch(FancyArrowPatch((cxy[0] + dx * 0.35, cxy[1] + dy * 0.35), (x - dx / L * 120, y - dy / L * 120), arrowstyle="-|>", mutation_scale=22,
                                         color="darkorange", lw=3, zorder=4, alpha=0.9))
    for n in ("WELL-4", "WELL-8"):
        x, y = by[n]["wellhead_xy"]
        ax.add_patch(plt.Circle((x, y), 650, fc="dodgerblue", alpha=0.18, ec="b", lw=1.5, zorder=2))
    ax.text(by["WELL-4"]["wellhead_xy"][0] + 350, by["WELL-4"]["wellhead_xy"][1] - 1100, "cold downflow /\nrecharge (NE)", color="navy", fontsize=10, fontweight="bold")
    ax.text(by["WELL-8"]["wellhead_xy"][0] - 1700, by["WELL-8"]["wellhead_xy"][1] - 1000, "cold downflow /\nrecharge (SW)", color="navy", fontsize=10, fontweight="bold")
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    ax.ticklabel_format(style="plain")
    ax.tick_params(axis="x", labelrotation=30, labelsize=8)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    # right panel: wells ordered by distance from the hot-core centre (a projection onto one line would
    # put wells up to 3 km off the line); x positions are spread so that wells at similar distances do not overlap
    dist = {n: float(np.hypot(*(np.array(by[n]["wellhead_xy"]) - cxy))) for n in summ}
    order = sorted(summ, key=lambda n: dist[n])
    xpos, last = {}, -1e9
    for n in order:
        xpos[n] = max(dist[n], last + 330.0)
        last = xpos[n]
    ax = axs[1]
    xs_ = [xpos[n] for n in order]
    Sx = S.set_index("well")
    ax.plot(xs_, [by[n]["wellhead_z"] for n in order], "k-", lw=1.2, label="wellhead elevation")
    ax.plot(xs_, [Sx.loc[n, "water_level_mRL"] for n in order], "-.", color="C0", label="estimated water level")
    for k, n in enumerate(order):
        w, r = by[n], summ[n]
        ax.plot([xpos[n]] * 2, [w["wellhead_z"], w["td_z"]], color="0.35", lw=1.2, zorder=2)
        for f in w["fz"]:
            ax.plot([xpos[n] - 70, xpos[n] + 70], [f["mid"]] * 2, color="darkorange", lw=5, zorder=3)
        zt, zb = hot_zone(profs[n])
        if r["Tmax_degC"] >= 150:
            ax.add_patch(plt.Rectangle((xpos[n] - 130, zb), 260, zt - zb, fc="red", alpha=0.35, ec="none", zorder=1))
        ax.plot(xpos[n], r["Tmax_elev_mRL"], "*", color="r", ms=14, mec="k", zorder=5)
        ax.text(xpos[n], 330 + 120 * (k % 2), f"{short(n)} ({dist[n] / 1000:.1f} km)", ha="center", fontsize=9, fontweight="bold")
        ax.text(xpos[n] + 150, w["td_z"], f"TD {r['T_at_TD_degC']:.0f} degC", fontsize=7, va="center", bbox=dict(fc="w", ec="none", alpha=0.8, pad=0.5))
    for key, c, lab in (("z_100C_mRL", "royalblue", "100 degC isotherm (first downward crossing)"), ("z_200C_mRL", "firebrick", "200 degC isotherm")):
        pts = [(xpos[n], summ[n][key]) for n in order if not is_nan(summ[n][key])]
        if pts:
            ax.plot(*zip(*pts), "o-", color=c, lw=2, ms=7, mec="k", label=lab, zorder=4)
    for n in order:
        if is_nan(summ[n]["z_100C_mRL"]):
            ax.annotate("", (xpos[n] + 110, by[n]["td_z"]), (xpos[n] + 110, by[n]["td_z"] + 260), arrowprops=dict(arrowstyle="->", color="royalblue", lw=2))
    ax.plot([], [], "*", color="r", ms=12, mec="k", label="Tmax")
    ax.add_patch(plt.Rectangle((0, 0), 0, 0, fc="red", alpha=0.35, label="within 15 degC of Tmax (hot layer)"))
    ax.plot([], [], color="darkorange", lw=5, label="feed zones")
    ax.plot([], [], color="royalblue", marker=r"$\downarrow$", ls="", ms=10, label="100 degC not reached by TD")
    ax.set_xlabel("Distance of the wellhead from the hot-core centre, mean of W10, W2, W3 (m); positions spread slightly to avoid overlap")
    ax.set_ylabel("Elevation (mRL)")
    ax.set_ylim(-3050, 650)
    ax.set_title("Evidence by distance from the hot core: 100/200 degC isotherms, Tmax, hot layer", fontsize=11)
    ax.grid(alpha=0.3)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=4, fontsize=8)
    axs[0].set_title("Plan view (interpretation): wells coloured by regime", fontsize=11)
    fig.suptitle("Conceptual model sketch built from the well data only (orange arrows: outflow away from the hot core; blue: cold downflow; magenta: two-phase)", fontsize=11)
    savefig(fig, "conceptual_model.png")


# ======================================================================================
# 8. assembly, writing, main
# ======================================================================================
SUMMARY_REQ = ["well", "Tmax_degC", "Tmax_elev_mRL", "T_at_TD_degC", "TD_mRL", "inversion_degC", "inversion_top_mRL",
               "inversion_base_mRL", "shallow_gradient_degC_per_km", "deep_gradient_degC_per_km", "z_100C_mRL", "z_200C_mRL",
               "min_boiling_margin_degC", "min_boiling_margin_elev_mRL", "water_level_mRL", "pivot_elev_mRL", "pivot_P_bar",
               "regime", "ns_verdict", "ns_confidence"]


def write_csv(df, name, nd=3):
    """Write a table with floats rounded to nd decimals; every number in the report comes from such a file."""
    d = df.copy()
    for c in d.columns:
        if d[c].dtype.kind == "f":
            d[c] = d[c].round(1 if c in ("x", "y", "x_wellhead", "y_wellhead") else nd)
    d.to_csv(os.path.join(OUT, name), index=False, na_rep="", lineterminator="\n")
    return d


def find_json_files(d):
    fs = [os.path.join(d, f) for f in os.listdir(d) if f.lower().endswith(".json")]
    def key(p):
        m = re.findall(r"\d+", os.path.basename(p))
        return (int(m[0]) if m else 0, os.path.basename(p))
    return sorted(fs, key=key)


def build_context():
    ctx = {}
    log("== IF97 validation (every equivalent error must be <= 0.1 K)")
    val = validate_if97()
    for r in val:
        log(f"   {r[0]:42s} expected {r[1]:<14.9g} computed {r[2]:<14.9g} equivalent error {r[3]:.2e}")
    ctx["val"] = val
    log("== input inventory")
    inv = build_inventory()
    write_inventory(inv)
    ctx["inv"] = inv
    for r in inv:
        log(f"   {r[0]:18s} {r[1]}")
    grid_p = role_paths(inv, "geometry_dat")[0]
    wb_p = role_paths(inv, "well_blocks_json")[0]
    grid, wb = read_grid(grid_p), read_well_blocks(wb_p)
    log(f"== grid: {len(grid['verts'])} vertices, {len(grid['cols'])} columns, {len(grid['layers'])} layers, "
        f"top surface {grid['top_surface']:.0f} mRL; well-to-block file: {len(wb)} wells")
    # layers tile the grid with no gaps
    lay = grid["layers"]
    gaps = [i for i in range(2, len(lay) + 1) if abs(lay[i]["top"] - lay[i - 1]["bottom"]) > 1e-9]
    assert not gaps, f"layer gaps at {gaps}"
    ctx["grid"], ctx["wb"] = grid, wb
    jd = role_paths(inv, "well_data_json")[0]
    W = [read_well_json(p) for p in find_json_files(jd)]
    log(f"== {len(W)} well JSON files read: " + ", ".join(w["name"] for w in W))
    feats = read_mapping(find_json_files(role_paths(inv, "mapping")[0])[0])
    groups = json.load(open(role_paths(inv, "well_groups")[0], encoding="utf-8"))
    ctx.update(W=W, feats=feats, groups=groups)
    qmap, qnotes = quality_scheme(W)
    ctx["qmap"], ctx["qnotes"] = qmap, qnotes
    profs, summ, isos = {}, {}, {}
    for w in W:
        p, r, iso = analyse_well(w, grid, wb, qmap, None)
        profs[w["name"]], summ[w["name"]], isos[w["name"]] = p, r, iso
    ctx.update(profs=profs, summ=summ, isos=isos)
    cc, refs = crosscheck_dat(W, role_paths(inv, "well_data_dat")[0])
    ctx["CC"] = cc
    ic = integrity_checks(W, grid, wb, groups)
    ctx["IC"] = ic
    ctx["PV"] = provenance_table(W, groups)
    out, trend, field_first, passed = natural_state_assessment(W, profs, summ, qmap, isos)
    A = {d["well"]: d for d in out}
    ctx.update(A=A, TR=trend, field_first=field_first, passed=passed)
    # per-well tables
    T, F, L = [], [], []
    for w in W:
        n = w["name"]
        T += build_targets(w, profs[n], summ[n], A[n], isos[n], grid, wb, qmap)
        F += build_feedzones(w, profs[n], A[n], wb, qmap)
        L += build_layer_centres(w, profs[n], grid, wb)
    ctx["T"], ctx["FZ"], ctx["L"] = pd.DataFrame(T), pd.DataFrame(F), pd.DataFrame(L)
    rows = []
    for w in W:
        r = dict(summ[w["name"]])
        r["ns_verdict"], r["ns_confidence"] = A[w["name"]]["verdict"], A[w["name"]]["confidence"]
        rows.append(r)
    S = pd.DataFrame(rows)
    S = S[SUMMARY_REQ + [c for c in S.columns if c not in SUMMARY_REQ]]
    ctx["S"] = augment_summary(S, W, grid, profs)
    ctx["HY"], ctx["HS"], ctx["hnotes"], ctx["comp"] = hydro_tables(W, profs, summ)
    log(f"== hydrostatic windows: {len(ctx['HY'])} tested; pieces not tested because shorter than {PRM['hydro_min_length_m']:.0f} m:")
    for x in ctx["hnotes"]:
        log("   " + x)
    ctx["GT"], ctx["faults"], ctx["waters"] = geometry_table(W, grid, feats)
    ctx["long"], ctx["wT"], ctx["wP"] = common_elevation_tables(W, profs)
    ctx["ISO"] = iso_table(W, isos)
    ctx["BF"] = boiling_flags(W, profs)
    ctx["REF"] = reference_model_summary(inv)
    return ctx


def pressure_consistency_df(ctx):
    rows = []
    for w in ctx["W"]:
        d = ctx["A"][w["name"]]
        rows.append(dict(
            well=w["name"], gauge_elev_mRL=d["gauge_elev"], survey_P_at_gauge_bar=d["survey_P_at_gauge"],
            first_gauge_date=d["first_gauge_date"].isoformat() if d["first_gauge_date"] else "",
            first_gauge_P_bar=d["first_gauge_P"], diff_survey_minus_first_bar=d["diff"],
            verdict=d["verdict"], confidence=d["confidence"],
            well_role=d["role"], well_first_flow_date=d["well_first_flow"].isoformat() if d["well_first_flow"] else "",
            survey_date_label=w["survey_date"],
            gauge_start_minus_field_first_production_days=d.get("days_gauge_start_after_field_first_flow"),
            gauge_start_minus_well_first_flow_days=d.get("days_gauge_start_after_well_first_flow"),
            gauge_interpolated_to_survey_date_bar=d.get("gauge_at_survey_date"),
            survey_minus_gauge_at_survey_date_bar=d.get("survey_minus_gauge_at_survey_date"),
            feedzone_trend_max_abs_residual_bar=d["trend_max_abs_resid"], feedzone_trend_mean_residual_bar=d["trend_mean_resid"],
            injectate_T_degC=d["Tinj"], injectate_enthalpy_kJ_kg=d["hinj"], reason=d["reason"]))
    return pd.DataFrame(rows)


def write_all_tables(ctx):
    os.makedirs(FIG, exist_ok=True)
    W = ctx["W"]
    out = {}
    out["targets"] = write_csv(ctx["T"], "NS_calibration_targets.csv", 4)
    out["fz"] = write_csv(ctx["FZ"], "feedzones.csv", 4)
    out["layers"] = write_csv(ctx["L"], "T_P_at_layer_centres.csv", 3)
    out["summary"] = write_csv(ctx["S"], "well_summary.csv", 3)
    out["pc"] = write_csv(pressure_consistency_df(ctx), "pressure_consistency_check.csv", 3)
    out["hydro"] = write_csv(ctx["HY"], "hydrostatic_check.csv", 3)
    out["hydro_sum"] = write_csv(ctx["HS"], "hydrostatic_summary.csv", 3)
    out["geom"] = write_csv(ctx["GT"], "well_geometry_check.csv", 1)
    out["trend"] = write_csv(ctx["TR"], "pressure_trend_check.csv", 3)
    out["iso"] = write_csv(ctx["ISO"], "isothermal_sections.csv", 2)
    out["bf"] = write_csv(ctx["BF"], "boiling_margin_flags.csv", 3)
    out["long"] = write_csv(ctx["long"], "TP_at_common_elevations.csv", 3)
    out["wT"] = write_csv(ctx["wT"], "temperature_at_common_elevations.csv", 3)
    out["wP"] = write_csv(ctx["wP"], "pressure_at_common_elevations.csv", 3)
    out["cc"] = write_csv(ctx["CC"].drop(columns=["content_matches"]), "dat_json_crosscheck.csv", 4)
    out["ic"] = write_csv(ctx["IC"], "integrity_checks.csv", 3)
    out["pv"] = write_csv(ctx["PV"], "provenance_labels.csv", 3)
    wrows = []
    for w in W:
        d = ctx["T"][ctx["T"].well == w["name"]]
        wrows.append(dict(well=w["name"], ns_confidence=ctx["A"][w["name"]]["confidence"], sum_T_weights=d.weight.sum(),
                          mean_T_weight=d.weight.mean(), sum_P_weights=d.weight_P.sum(),
                          n_points=len(d), n_points_in_feedzone=int(d.in_feedzone.sum())))
    out["wbw"] = write_csv(pd.DataFrame(wrows), "weights_by_well.csv", 3)
    write_csv(pd.DataFrame([dict(check=r[0], expected=r[1], computed=r[2], equivalent_error_K_or_relative=r[3]) for r in ctx["val"]]),
              "if97_validation.csv", 9)
    prm = [dict(parameter=k, value=v) for k, v in PRM.items()]
    prm += [dict(parameter="g_m_s2", value=G), dict(parameter="P_critical_bar", value=P_CRIT_BAR),
            dict(parameter="weights_feed_zone", value=W_FZ), dict(parameter="weights_open_hole", value=W_OPEN),
            dict(parameter="weights_open_hole_isothermal", value=W_ISO), dict(parameter="weights_behind_casing", value=W_CASED),
            dict(parameter="weights_top_30m_or_top_block", value=W_TOP), dict(parameter="weights_above_surface", value=W_ABOVE)]
    prm += [dict(parameter=f"confidence_factor_{k}", value=v) for k, v in CONF_FACTOR.items()]
    prm += [dict(parameter=f"quality_factor_{k}", value=v) for k, v in ctx["qmap"].items()]
    write_csv(pd.DataFrame(prm), "analysis_parameters.csv", 5)
    write_csv(pd.DataFrame([dict(item=k, value=v) for k, v in ctx["REF"].items()]), "reference_model_summary.csv", 4)
    return out


def make_all_figures(ctx):
    W, profs, summ = ctx["W"], ctx["profs"], ctx["summ"]
    for f in glob.glob(os.path.join(FIG, "*.png")):
        os.remove(f)
    fig_well_map(W, ctx["grid"], ctx["feats"], summ)
    for w in W:
        n = w["name"]
        fig_profile(w, profs[n], summ[n], ctx["A"][n], ctx["isos"][n], ctx["qmap"])
    fig_maps(W, ctx["long"], ctx["feats"], "T_degC", "maps_temperature", "inferno", "degC", "Temperature")
    fig_maps(W, ctx["long"], ctx["feats"], "P_bar", "maps_pressure", "viridis", "bar", "Pressure")
    hot = sorted(summ, key=lambda k: -summ[k]["Tmax_degC"])[:4]
    ctx["hot4"] = hot
    fig_section(W, profs, summ, hot, "T", "section_hottest_wells.png", "inferno", "Temperature", "degC", np.arange(0, 281, 20))
    fig_section(W, profs, summ, hot, "P", "section_pressure_hottest_wells.png", "viridis", "Pressure", "bar", np.arange(0, 301, 20))
    fig_pt(W, profs, ctx["feats"])
    fig_conceptual(W, profs, summ, ctx["S"], ctx["grid"], ctx["feats"], ctx["A"])


class Tee:
    def __init__(self, stream):
        self.stream, self.lines = stream, []

    def write(self, s):
        self.stream.write(s)
        self.lines.append(s)

    def flush(self):
        self.stream.flush()


def main():
    tee = Tee(sys.stdout)
    sys.stdout = tee
    try:
        os.makedirs(FIG, exist_ok=True)
        ctx = build_context()
        tabs = write_all_tables(ctx)
        ctx["tabs"] = tabs
        make_all_figures(ctx)
        text, derived = write_report(ctx)
        with open(os.path.join(OUT, "NS_data_interpretation.md"), "w", encoding="utf-8", newline="\n") as fh:
            fh.write(text)
        write_csv(pd.DataFrame(derived, columns=["name", "value", "unit", "note"]), "report_derived_numbers.csv", 5)
        log("== outputs written")
        for f in sorted(os.listdir(OUT)):
            log("   Output/" + f)
        log(f"   Output/figures: {len(os.listdir(FIG))} figures")
    finally:
        sys.stdout = tee.stream
    with open(os.path.join(OUT, "run_console_log.txt"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write("".join(tee.lines))



# ======================================================================================
# 7c. report
# ======================================================================================
def md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |", "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(c) for c in r) + " |")
    return "\n".join(out)


def write_report(ctx):
    tb = ctx["tabs"]
    S = tb["summary"].set_index("well")
    PC = tb["pc"].set_index("well")
    FZ = tb["fz"]
    HY, HS = tb["hydro"], tb["hydro_sum"]
    GT = tb["geom"].set_index("well")
    TG = tb["targets"]
    TR = tb["trend"]
    wP = tb["wP"].set_index("well")
    wT = tb["wT"].set_index("well")
    ISO, BF = tb["iso"], tb["bf"]
    IC, CC = tb["ic"], tb["cc"]
    W = ctx["W"]
    wells = [w["name"] for w in W]
    derived = []

    def dv(name, val, unit="", note="", nd=1):
        v = round(float(val), nd)
        derived.append((name, v, unit, note))
        return f"{v:.{nd}f}"

    def f(x, nd=1):
        return "" if (x is None or (isinstance(x, float) and math.isnan(x))) else f"{x:.{nd}f}"

    def lst(names):
        return ", ".join(short(n) for n in names)

    def andlist(names):
        t = [short(n) for n in names]
        return t[0] if len(t) == 1 else ", ".join(t[:-1]) + " and " + t[-1]

    L = []
    add = L.append

    # ---------- numbers used in several places ----------
    n_w = len(wells)
    nat = [n for n in wells if S.loc[n, "ns_verdict"] == "natural state"]
    prob = [n for n in wells if S.loc[n, "ns_verdict"] == "probably natural state"]
    poss = [n for n in wells if S.loc[n, "ns_verdict"] == "possibly disturbed"]
    conf = {c: [n for n in wells if S.loc[n, "ns_confidence"] == c] for c in ("high", "medium", "low")}
    reg = {r: [n for n in wells if S.loc[n, "regime"] == r] for r in ("upflow", "outflow", "conductive", "recharge/downflow", "two-phase")}
    hot3 = [n for n in wells if S.loc[n, "Tmax_degC"] >= PRM["hot_core_Tmax_degC"]]
    cold = [n for n in wells if S.loc[n, "Tmax_degC"] < PRM["cold_well_Tmax_degC"]]
    cc_ok = int((CC.status == "OK").sum())
    irregular = CC[CC.irregular_name.notna() & (CC.irregular_name.astype(str) != "")].file.tolist()
    n_pass, n_issue, n_note = int((IC.status == "PASS").sum()), int((IC.status == "ISSUE").sum()), int((IC.status == "NOTE").sum())
    comp_c, comp_o = ctx["comp"]["cased"], ctx["comp"]["open"]
    hy_c = HY[HY.band == "cased"]
    hy_o = HY[HY.band == "open"]
    hs = HS.set_index(["scope", "band"])
    pw = lambda n, col, nd=1: f(S.loc[n, col], nd)
    tmax = lambda n: f(S.loc[n, "Tmax_degC"], 1)
    hz = lambda n: f(S.loc[n, "hot_zone_thickness_m"], 0)
    topT = lambda n: f(TG[TG.well == n].sort_values("elevation_mRL").T_degC.iloc[-1], 1)
    topdrop = lambda n: S.loc[n, "wellhead_mRL"] - TG[TG.well == n].elevation_mRL.max()
    w4d = [w for w in W if w["name"] == "WELL-4"][0]["liner_bottom_z"]
    p1500 = wP["z_-1500_mRL"].dropna()
    p2000 = wP["z_-2000_mRL"].dropna()
    p2500 = wP["z_-2500_mRL"].dropna()
    p1500x = p1500.drop("WELL-8", errors="ignore")
    tr_max = TR.residual_bar.abs().max()
    wl_ext = S[S.water_level_extrapolation_m > 5]
    wh_below = GT[GT.wellhead_minus_surface_m < -5]

    # ================= title =================
    add("# Natural-state data interpretation for model calibration")
    add("")
    add("Study area: the ten-well synthetic teaching dataset in `Baseline/` (all well files are labelled 'Synthetic data'). "
        "Everything below is generated by `Scripts/natural_state_analysis.py`; every number is in a file of `Output/`. "
        "Claims are labelled **[observation]** (read directly from the data), **[interpretation]** (inference from the data), "
        "**[assumption]** (a value not in the data) or **[reference model]** (a property of the starting model; never evidence). "
        "Units: mRL, bar, degC, kJ/kg, kg/s.")
    add("")

    # ================= 1 key conclusions =================
    add("## 1. Key conclusions for calibration")
    add("")
    gd = PC.loc["WELL-2", "diff_survey_minus_first_bar"]
    add(f"1. **The surveys can be used as natural-state targets.** [observation] The T and P surveys of all {n_w} wells carry one date label (29/07/2015). "
        f"The gauge test gives {andlist(nat)} a match to the first gauge reading within the reading precision "
        f"(survey minus first gauge {f(PC.loc['WELL-2', 'diff_survey_minus_first_bar'], 2)} and {f(PC.loc['WELL-3', 'diff_survey_minus_first_bar'], 2)} bar; both gauge series start on the first day of production in the field). "
        f"The other {len(prob)} wells are 'probably natural state' (feed-zone pressures within {dv('max_trend_resid', tr_max, 'bar', 'largest feed-zone pressure residual against the passed wells', 1)} bar of the trend; W9 survey is {f(PC.loc['WELL-9', 'diff_survey_minus_first_bar'], 2)} bar above its first gauge value and W10's gauge starts years after field production began). "
        f"**Action:** use all wells, with the confidence factors in `NS_calibration_targets.csv` (high {lst(conf['high'])}; low {lst(conf['low'])}; the rest medium).")
    add(f"2. **Pressure targets belong at feed zones only.** [observation] {len(FZ)} feed zones carry a formation-pressure estimate (`feedzones.csv`); they agree with each other across the field within {dv('tr_res_range', tr_max, 'bar', 'max abs residual', 1)} bar at equal elevation. "
        f"Between feed zones the surveyed pressure follows the well column (pivot at the dominant feed zone, `well_summary.csv`). **Action:** `weight_P` is non-zero only at feed-zone survey points; do not fit shallow wellbore pressures (about 1 bar at the top of several wells).")
    add(f"3. **The pressure gradients are steeper than an IF97 liquid column and this is the largest unresolved data problem.** [observation] "
        f"Cold wells (Tmax < {PRM['cold_well_Tmax_degC']:.0f} degC; {lst(cold)}) give a measured-to-IF97 pressure-change ratio of {f(comp_c, 3)} in the cased band and {f(comp_o, 3)} in the open band (medians). "
        f"The hot wells exceed this: the cased windows reach {f(hy_c.ratio_measured_to_IF97.max(), 3)} ({short(HY.loc[hy_c.ratio_measured_to_IF97.idxmax(), 'well'])}) and the open windows reach {f(hy_o.ratio_measured_to_IF97.max(), 3)} ({short(HY.loc[hy_o.ratio_measured_to_IF97.idxmax(), 'well'])}); a ratio above 1 means an implied density above that of water at the measured temperature. "
        f"Datum offsets, a gravity/unit convention, two-phase fluid and rounding cannot explain the hot-well excess (section 4c). **Action:** compare model and data at feed zones with a pressure tolerance of several bar, "
        f"and do not trust boiling margins at the shallow hot wells until the pressure convention is confirmed.")
    z200 = {n: S.loc[n, "z_200C_mRL"] for n in hot3}
    add(f"4. **A hot core with a thin cap sits at {lst(hot3)}.** [observation] Tmax is {', '.join(f'{tmax(n)} degC ({short(n)})' for n in hot3)} at only {f(S.loc['WELL-2', 'Tmax_elev_mRL'], 0)} to {f(S.loc['WELL-10', 'Tmax_elev_mRL'], 0)} mRL; the 200 degC isotherm is at "
        f"{', '.join(f'{f(z200[n], 1)} mRL ({short(n)})' for n in hot3)}, and the 100 degC isotherm is only {', '.join(f'{f(S.loc[n, 'depth_100C_below_wellhead_m'], 1)} m' for n in hot3)} below the wellheads. "
        f"[interpretation] This is the upflow zone, closest to Fault 3 (W2 {f(GT.loc['WELL-2', 'dist_nearest_fault_m'], 0)} m, W10 {f(GT.loc['WELL-10', 'dist_nearest_fault_m'], 0)} m). "
        f"**Action:** the model needs a deep hot source (>= {f(S.loc['WELL-2', 'T_at_TD_degC'], 1)} degC, still rising at the bottom of W2) under or near these wells and a low-permeability cap above it.")
    add(f"5. **Lateral outflow leaves the hot core at about -450 to -750 mRL.** [observation] {lst(reg['outflow'])} have their temperature maximum between {f(min(S.loc[n, 'Tmax_elev_mRL'] for n in reg['outflow']), 0)} and {f(max(S.loc[n, 'Tmax_elev_mRL'] for n in reg['outflow']), 0)} mRL "
        f"and a temperature inversion of {f(min(S.loc[n, 'inversion_degC'] for n in reg['outflow']), 1)} to {f(max(S.loc[n, 'inversion_degC'] for n in reg['outflow']), 1)} degC below it; "
        f"Tmax falls with distance from the core ({', '.join(f'{short(n)} {tmax(n)}' for n in sorted(reg['outflow'], key=lambda k: -S.loc[k, 'Tmax_degC']))} degC). "
        f"[interpretation] A permeable lateral layer carries hot water south (W6, W7), east (W5) and north-east (W9). **Action:** give the model a laterally permeable horizon at these elevations; the feed zones of the wells (`feedzones.csv`) mark where it is permeable.")
    bf = BF.iloc[0] if len(BF) else None
    if bf is not None:
        add(f"6. **Boiling-point-for-depth conditions are seen only at W10.** [observation] The boiling margin falls to {f(bf.min_margin_degC, 1)} degC at {f(bf.min_margin_elev_mRL, 0)} mRL ({bf.n_points} survey points below {PRM['bpd_margin_degC']:.0f} degC between {f(bf.z_top_mRL, 0)} and {f(bf.z_bottom_mRL, 1)} mRL), at the feed zone the margin is {f(FZ[FZ.well == 'WELL-10'].boiling_margin_degC.iloc[0], 1)} degC. "
            f"[interpretation] A two-phase or boiling zone is possible above the upflow; it depends on the pressure convention in conclusion 3: if the surveyed pressures are too high, the true Tsat is lower and every margin in the hot sections is smaller than listed. **Action:** allow a boiling zone near W10 and treat the smallest margins of W3 ({f(S.loc['WELL-3', 'min_boiling_margin_degC'], 1)} degC) and W2 ({f(S.loc['WELL-2', 'min_boiling_margin_degC'], 1)} degC) as upper bounds until the pressure convention is confirmed.")
    add(f"7. **Cold downflow or recharge at both margins.** [observation] W4 (north-east) reaches only {pw('WELL-4', 'Tmax_degC')} degC at {pw('WELL-4', 'TD_mRL', 0)} mRL (shallow gradient {pw('WELL-4', 'shallow_gradient_degC_per_km')} degC/km) and W8 (south-west) stays at {f(min(wT.loc['WELL-8', 'z_-1500_mRL'], wT.loc['WELL-8', 'z_-500_mRL'], wT.loc['WELL-8', 'z_-1000_mRL']), 1)} to {f(max(wT.loc['WELL-8', 'z_-1500_mRL'], wT.loc['WELL-8', 'z_-500_mRL'], wT.loc['WELL-8', 'z_-1000_mRL']), 1)} degC between -500 and -1500 mRL. "
        f"W8's shallow gradient ({pw('WELL-8', 'shallow_gradient_degC_per_km')} degC/km) is close to W1's ({pw('WELL-1', 'shallow_gradient_degC_per_km')} degC/km), so its cooling appears only below -550 mRL, whereas W4 is cold from the surface. [interpretation] Cold or cool water descends or enters laterally at the margins. **Action:** model cold recharge on the NE and SW flanks and a low basal heat input there; with an assumed rock conductivity of {PRM['k_rock_W_mK']:.1f} W/m/K [assumption] the shallow gradients of W4 and W1 correspond to {pw('WELL-4', 'shallow_heat_flux_mW_m2_assumed_k')} and {pw('WELL-1', 'shallow_heat_flux_mW_m2_assumed_k')} mW/m2.")
    add(f"8. **Deep pressures are connected across the field.** [observation] At -1500 mRL nine wells lie within {dv('spread_1500', p1500.max() - p1500.min(), 'bar', 'spread of P at -1500 mRL, all wells', 1)} bar ({f(p1500.min(), 1)} to {f(p1500.max(), 1)} bar; {dv('spread_1500_x8', p1500x.max() - p1500x.min(), 'bar', 'spread of P at -1500 mRL without W8', 1)} bar without the cold W8). "
        f"[interpretation] No lateral pressure barrier separates the wells at depth; W8's higher deep pressure is what a colder, denser column gives. **Action:** give the deep model high lateral permeability or a common deep pressure datum; fit the {len(FZ)} feed-zone pressures jointly.")
    add(f"9. **Geometry issues to fix before calibration.** [observation] {len(wh_below)} wellheads lie below the model surface of their column ({', '.join(f'{short(n)} by {f(-GT.loc[n, 'wellhead_minus_surface_m'], 0)} m' for n in wh_below.index)}); the W4 liner is designed to {f(w4d, 0)} mRL but the track ends at {f(S.loc['WELL-4', 'TD_mRL'], 0)} mRL; "
        f"the W1 water level extrapolates {f(S.loc['WELL-1', 'water_level_extrapolation_m'], 1)} m above its top survey point. All survey points and feed zones fall in the block given by the well-to-block file (own point-in-polygon search agrees at every survey point).")
    add("")

    # ================= 2 the data =================
    add("## 2. The data")
    add("")
    add("### 2a. What exists (`input_inventory.csv`)")
    add("")
    inv = pd.read_csv(os.path.join(OUT, "input_inventory.csv"))
    rows = []
    for role, g in inv.groupby("role"):
        rows.append((role, len(g), "; ".join(sorted(set(x.split("/")[-1] for x in g.path)))[:150]))
    add(md_table(["role", "rows", "items"], rows))
    add("")
    add("The ten per-well JSON files hold the track, casing/liner design, feed zones, one T and one P survey, and (for producers and injectors) flow, enthalpy and transient-pressure records. "
        "`Data_dat` holds plain-column copies, `Mapping.json` three faults ('Fault 1', 'Fault 2', 'Fault 3'), 'Lake Tarawera', 'Tarawera River' and 'river'. "
        "The reference model files (rock types, boundary conditions, incons, set-up and plotting files) are not evidence. The production history, forecast set-up and field-total flow file are out of scope. "
        "No role in section 2.2 of the instructions is empty; `Data/Make_up_json` is an empty directory and `Plotting` holds hard-coded paths of other machines.")
    add("")
    add("### 2b. How the data were checked")
    add("")
    add(f"- **JSON versus plain-column files:** {cc_ok} of {len(CC)} files match their JSON series exactly (largest value difference {f(CC.max_abs_diff.max(), 3)}). "
        f"Irregular names: {len(irregular)} files use the prefix 'Wel' instead of 'Well' ({irregular[0].split('_')[0]}; the others have the regular form 'WellN_...'); they were matched to the well by content as well as by name. "
        f"The plain-column time axes are years since 01/01/2014 (reference date found by matching the series against the JSON dates); flow files in `Data_dat` are stored with a negative sign for production and positive for injection. "
        f"[observation]")
    add(f"- **Track, survey, design and block file (`integrity_checks.csv`):** {n_pass} checks pass, {n_note} notes, {n_issue} issues. Every survey point and feed zone lies inside its track; every feed zone lies in the slotted/open interval; the JSON tracks equal the grid-file tracks. "
        f"Issues: W4 liner bottom ({f(w4d, 0)} mRL) versus track bottom ({f(S.loc['WELL-4', 'TD_mRL'], 0)} mRL); W10's feed-zone tag 'Synthetic' is not a quality grade. Notes: the W4 and W10 surveys start {f(topdrop('WELL-4'), 0)} m and {f(topdrop('WELL-10'), 0)} m below the wellhead.")
    add("- **Model blocks:** the grid file was parsed (vertex lines have the two coordinates run together, so they were read with a fixed pattern instead of whitespace; the other sections by whitespace). "
        "Layers tile the grid without gaps. For every survey point my own point-in-polygon and layer search gives the same block as the well-to-block file, which is used as the authority.")
    add("- **Thermodynamics:** IAPWS-IF97 Region 4 (Tsat, Psat) and Region 1 (liquid density, enthalpy) are implemented in the script and validated against the verification tables (`if97_validation.csv`; the script stops if any equivalent error exceeds 0.1 K). "
        f"Tsat is left blank above the critical pressure ({P_CRIT_BAR:.2f} bar); the deepest survey points of {andlist([w['name'] for w in W if np.nanmax(ctx['profs'][w['name']]['P']) > P_CRIT_BAR])} exceed it.")
    add("")
    add("### 2c. Provenance and quality labels (`provenance_labels.csv`)")
    add("")
    add("- Casing/liner design: Source = 'Synthetic data' in all ten wells. The group file has one area, 'Synthetic', with all ten wells.")
    add("- Survey `Source` and `Comment` fields are empty; every survey has `StableFlag` = True (a label only: the data do not show how stability was judged, so it is not treated as evidence of thermal recovery). Flow, enthalpy and transient-pressure `Source`/`Comment` fields are empty. Feed-zone `Source` is only a unit label ('Text').")
    add("- Survey date label in all wells: 29/07/2015. Feed-zone quality tags found and their factors for `weight_P`:")
    add("")
    add(md_table(["tag", "wells", "factor", "mapping"], [(t, lst([w['name'] for w in W if any(x['quality'] == t for x in w['fz'])]), f(ctx['qmap'][t], 1), ctx['qnotes'][t]) for t in sorted(ctx['qmap'], key=lambda t: -ctx['qmap'][t])]))
    add("")
    add("### 2d. Do the surveys represent the natural state? (`pressure_consistency_check.csv`)")
    add("")
    add("**Rule used.** The survey pressure interpolated to the gauge elevation is compared with the first transient value. A match within the reading precision decides 'natural state', whatever the date label says. "
        f"Wells without a gauge are compared at their feed zones with the mean pressure of the wells that passed (here {lst(ctx['passed'])}) at the same elevation (limit {PRM['fz_trend_ok_bar']:.0f} bar). "
        "Injectors: the injectate temperature is derived from the median injection enthalpy with IF97 Region 1 (liquid at 10 bar); a feed-zone temperature colder than the injectate with no isothermal section near it argues against an injection-affected survey, a hotter formation could hide cooling and lowers the confidence one step. "
        "Confidence: high = match to the pre-/start-of-production gauge; medium = probable; low = probable but injection cooling could be hidden.")
    add("")
    rows = []
    for n in wells:
        r = PC.loc[n]
        rows.append((short(n), r.well_role, r.well_first_flow_date if isinstance(r.well_first_flow_date, str) else "-", f(r.gauge_elev_mRL, 0) or "-", f(r.survey_P_at_gauge_bar, 2) or "-",
                     (r.first_gauge_date if isinstance(r.first_gauge_date, str) else "-"), f(r.first_gauge_P_bar, 1) or "-", f(r.diff_survey_minus_first_bar, 2) or "-",
                     f(r.feedzone_trend_max_abs_residual_bar, 2), r.verdict, r.confidence))
    add(md_table(["well", "role", "first flow", "gauge elev", "survey P at gauge", "first gauge", "first gauge P", "survey - first", "max FZ residual", "verdict", "confidence"], rows))
    add("")
    add(f"[observation] The gauge series of W2, W3 and W9 start on the day the first production record starts (also the first production anywhere in the field). W10's gauge starts {int(PC.loc['WELL-10', 'gauge_start_minus_field_first_production_days'])} days after field production began "
        f"and {int(PC.loc['WELL-10', 'gauge_start_minus_well_first_flow_days'])} days after the well's own first flow, so it cannot be a pre-production reference. W10's own first flow ({PC.loc['WELL-10', 'well_first_flow_date']}) is years later than the survey label.")
    add(f"**Date label versus gauge (data question).** Interpolating each gauge series to the survey date (29/07/2015) gives {f(PC.loc['WELL-2', 'gauge_interpolated_to_survey_date_bar'], 1)}, {f(PC.loc['WELL-3', 'gauge_interpolated_to_survey_date_bar'], 1)} and {f(PC.loc['WELL-9', 'gauge_interpolated_to_survey_date_bar'], 1)} bar for W2, W3 and W9, "
        f"so the gauge values are {f(PC.loc['WELL-2', 'survey_minus_gauge_at_survey_date_bar'], 1)}, {f(PC.loc['WELL-3', 'survey_minus_gauge_at_survey_date_bar'], 1)} and {f(PC.loc['WELL-9', 'survey_minus_gauge_at_survey_date_bar'], 1)} bar below the surveyed pressures. "
        "A survey taken on that date under production would sit at the gauge level, so the date label conflicts with the pre-production match. As instructed this is reported as a data question and does not lower the verdict.")
    add(f"**W9.** The survey exceeds its first gauge value by {f(PC.loc['WELL-9', 'diff_survey_minus_first_bar'], 2)} bar. Production or injection would lower, not raise, the survey pressure relative to the gauge, "
        "and the first gauge value is already at the production start at the largest flow rate in the field; the verdict is therefore 'probably natural state', medium.")
    add(f"**Injectors.** The injectate is {f(PC.loc['WELL-4', 'injectate_enthalpy_kJ_kg'], 0)} kJ/kg, i.e. {f(PC.loc['WELL-4', 'injectate_T_degC'], 1)} degC. "
        f"The feed-zone temperatures of W4 ({f(FZ[FZ.well == 'WELL-4'].T_degC.iloc[0], 1)} degC) and W8 ({f(FZ[FZ.well == 'WELL-8'].T_degC.iloc[0], 1)} degC) are colder than the injectate and neither profile has a section near it: evidence against an injection-affected survey. "
        f"W5's feed zones ({', '.join(f(x, 1) for x in FZ[FZ.well == 'WELL-5'].T_degC)} degC) are hotter than the injectate, so cooling could be hidden: confidence low.")
    add("")
    add("**Effect on the weights.** The confidence differs between wells (factors 1.0, 0.75, 0.5), so it does change relative weights: the temperature weights of W2 and W3 are twice those of W5 for the same location class. "
        "Location weights follow the defaults of the instructions without departure.")
    add("")

    # ================= 3 temperature =================
    add("## 3. Temperature structure")
    add("")
    add("Definitions (reproducible): inversion = largest drop below a running maximum if >= 3 degC; shallow gradient = least-squares fit over the top 300 m of the liquid-filled survey (liquid begins at the estimated water level); "
        "deep gradient = least-squares fit below both the inversion base and the deepest feed zone, blank if that interval is shorter than 300 m (two-point fits are exact slopes); isotherm crossings are the first downward crossings (linear interpolation); "
        "regime: two-phase if any boiling margin < 5 degC at P >= 1.5 bar; recharge/downflow if Tmax < 100 degC; outflow if the inversion is >= 10 degC; upflow if Tmax >= 200 degC at the bottom of the well; otherwise conductive.")
    add("")
    rows = []
    for n in wells:
        r = S.loc[n]
        inv_txt = f"{f(r.inversion_degC, 1)} ({f(r.inversion_top_mRL, 0)} to {f(r.inversion_base_mRL, 0)})" if not math.isnan(r.inversion_degC) else "none"
        rows.append((short(n), f(r.Tmax_degC, 1), f(r.Tmax_elev_mRL, 0), f(r.T_at_TD_degC, 1), f(r.TD_mRL, 0), inv_txt, f(r.shallow_gradient_degC_per_km, 0),
                     f(r.deep_gradient_degC_per_km, 1) or "-", f(r.z_100C_mRL, 0) or "-", f(r.z_200C_mRL, 0) or "-", f(r.min_boiling_margin_degC, 1), r.regime))
    add(md_table(["well", "Tmax degC", "at mRL", "T at TD degC", "TD mRL", "inversion degC (top to base mRL)", "shallow grad degC/km", "deep grad degC/km", "z 100C mRL", "z 200C mRL", "min margin degC", "regime"], rows))
    add("")
    add("**Observations.**")
    add("")
    add(f"- The three hottest wells ({lst(hot3)}) reach {f(min(S.loc[n, 'Tmax_degC'] for n in hot3), 1)} to {f(max(S.loc[n, 'Tmax_degC'] for n in hot3), 1)} degC between {f(max(S.loc[n, 'Tmax_elev_mRL'] for n in hot3), 0)} and {f(min(S.loc[n, 'Tmax_elev_mRL'] for n in hot3), 0)} mRL; "
        f"the layer within {PRM['hot_layer_drop_degC']:.0f} degC of Tmax is {f(S.loc['WELL-3', 'hot_zone_thickness_m'], 0)} m thick at W3 and {f(S.loc['WELL-2', 'hot_zone_thickness_m'], 0)} m at W2 (`well_summary.csv`).")
    add(f"- Depth of the 100 degC isotherm below the wellhead: {', '.join(f'{short(n)} {f(S.loc[n, 'depth_100C_below_wellhead_m'], 0)} m' for n in wells if not math.isnan(S.loc[n, 'depth_100C_below_wellhead_m']))}. "
        f"{lst([n for n in wells if math.isnan(S.loc[n, 'depth_100C_below_wellhead_m'])])} never reach 100 degC within their surveys. The cap (to 100 degC) thickens from {f(S.depth_100C_below_wellhead_m.min(), 0)} m over the hot core to {f(S.depth_100C_below_wellhead_m.max(), 0)} m at W1.")
    add(f"- Shallow gradients range from {f(S.shallow_gradient_degC_per_km.min(), 0)} (W4) to {f(S.shallow_gradient_degC_per_km.max(), 0)} degC/km (W3). With an assumed conductivity of {PRM['k_rock_W_mK']:.1f} W/m/K [assumption, typical rock value; not in the data] these are {f(S.shallow_heat_flux_mW_m2_assumed_k.min(), 0)} to {f(S.shallow_heat_flux_mW_m2_assumed_k.max(), 0)} mW/m2 "
        f"(column `shallow_heat_flux_mW_m2_assumed_k`). [interpretation] Values above a few hundred mW/m2 in the hot wells are probably not a basal conductive flux but the steep cap above a convective zone. The deep gradients are only {f(S.deep_gradient_degC_per_km.min(), 1)} to {f(S.deep_gradient_degC_per_km.max(), 1)} degC/km, i.e. nearly isothermal at depth.")
    add("- Isothermal sections of 300 m or more (`isothermal_sections.csv`; range <= 5 degC):")
    add("")
    add(md_table(["well", "top mRL", "bottom mRL", "length m", "T range degC", "location", "feed zones inside"],
                 [(short(r.well), f(r.z_top_mRL, 0), f(r.z_bottom_mRL, 0), f(r.length_m, 0), f(r.T_min_degC, 1) + " to " + f(r.T_max_degC, 1), r.location,
                   (str(r.feedzones_inside) if isinstance(r.feedzones_inside, str) or not math.isnan(r.feedzones_inside) else "none")) for r in ISO.itertuples()]))
    add("")
    add("- Inversions against the casing shoe and feed zones:")
    add("")
    rows = []
    for w in W:
        r = S.loc[w["name"]]
        if math.isnan(r.inversion_degC):
            continue
        fzin = [str(x["n"]) for x in w["fz"] if r.inversion_top_mRL + 1e-9 >= x["mid"] >= r.inversion_base_mRL - 1e-9]
        pos = ("below the shoe (open hole)" if r.inversion_top_mRL <= w["shoe_z"] + 1e-9
               else "starts behind casing, continues in open hole" if r.inversion_base_mRL < w["shoe_z"] else "behind casing")
        rows.append((short(w["name"]), f(r.inversion_degC, 1), f(r.inversion_top_mRL, 0), f(r.inversion_base_mRL, 0), f(w["shoe_z"], 0), pos,
                     ",".join(fzin) if fzin else "none"))
    add(md_table(["well", "drop degC", "top mRL", "base mRL", "casing shoe mRL", "position", "feed zones inside"], rows))
    add("")
    add("[observation] Each inversion starts below the casing shoe or straddles it, so the temperatures that define the maximum and the inversion are open-hole values that could be affected by flow in the well. "
        "The inversions of W2, W3, W5, W6, W7 and W9 are large and span feed zones, which is what outflow past the feed zones would produce [interpretation], but the data cannot exclude internal flow.")
    add("")
    add("**Wellbore check before interpreting.** [observation] Most of these sections are in open hole; the W3 section (-1100 to -2250 mRL, 174.5 to 177.2 degC) lies below the deepest feed zone and contains none, so flow between feed zones cannot produce it in the usual way. "
        "W4's cased section is a cold conductive profile. [interpretation] The deep isothermal sections of W2, W5, W6, W7, W8 are probably formation temperatures homogenised by deep convection or lateral flow, but internal wellbore flow cannot be excluded; "
        "open-hole points inside them get location weight 0.5 instead of 0.8. A repeat survey after a long shut-in, or a flow-spinner log, would decide.")
    add("")
    add(f"**Boiling margin (`boiling_margin_flags.csv`).** Only W10 has intervals below {PRM['bpd_margin_degC']:.0f} degC: {', '.join(f'{short(r.well)} {f(r.z_top_mRL, 0)} to {f(r.z_bottom_mRL, 1)} mRL (minimum {f(r.min_margin_degC, 1)} degC at {f(r.min_margin_elev_mRL, 0)} mRL, T {f(r.T_at_min_degC, 1)} degC, P {f(r.P_at_min_bar, 1)} bar)' for r in BF.itertuples())}. "
        f"The next smallest minimum is {f(S.drop('WELL-10').min_boiling_margin_degC.min(), 1)} degC ({short(S.drop('WELL-10').min_boiling_margin_degC.idxmin())}). The W10 profile is isothermal at {f(S.loc['WELL-10', 'T_at_TD_degC'], 1)} to {f(S.loc['WELL-10', 'Tmax_degC'], 1)} degC from its feed zone to the bottom, "
        "and its P-T path runs along the saturation curve between {f(BF.z_top_mRL.iloc[0], 0)} and {f(BF.z_bottom_mRL.iloc[0], 1)} mRL (`figures/PT_state.png`, right panel): a boiling-point-for-depth section **[interpretation]** that could be two-phase or a hot liquid at the edge of boiling.")
    add("")
    add("![Temperature maps](figures/maps_temperature_minus500m.png)")
    add("")
    add("Maps at 0, -500, -1000 and -1500 mRL (`figures/maps_temperature_*.png`) and sections (`figures/section_hottest_wells.png`) show the hot core near W10, W2 and W3 with temperature falling toward W6, W7, W5, W9 and the cold flanks W4 and W8; interpolation is inside the hull of wells with data only.")
    add("")

    # ================= 4 pressure =================
    add("## 4. Pressure structure")
    add("")
    add("### 4a. Pressure at common elevations (`pressure_at_common_elevations.csv`, bar; interpolation only)")
    add("")
    cols = [c for c in tb["wP"].columns if c != "well"]
    cols = sorted(cols, key=lambda c: -int(c.split("_")[1]))
    hdr = ["well"] + [c.replace("z_", "").replace("_mRL", "") + " mRL" for c in cols]
    add(md_table(hdr, [[short(n)] + [f(wP.loc[n, c], 1) or "-" for c in cols] for n in wells]))
    add("")
    add(f"[observation] Shallow pressures differ because the water levels differ: estimated levels (straight line through the shallowest three points with P >= {PRM['hydro_min_P_bar']} bar, solved for atmospheric pressure) range from {f(S.water_level_mRL.min(), 0)} (W1) to {f(S.water_level_mRL.max(), 0)} mRL (W4) "
        f"and lie {f(S.water_level_below_model_surface_m.min(), 0)} to {f(S.water_level_below_model_surface_m.max(), 0)} m below the model surface of their columns (`well_summary.csv`). "
        f"The W1 level is extrapolated {f(S.loc['WELL-1', 'water_level_extrapolation_m'], 0)} m above its top survey point and W5's {f(S.loc['WELL-5', 'water_level_extrapolation_m'], 0)} m, so these two are poor. "
        f"At -1500 mRL the wells lie within {dv('spread_1500b', p1500.max() - p1500.min(), 'bar', 'spread of P at -1500 mRL', 1)} bar, at -2000 mRL within {dv('spread_2000', p2000.max() - p2000.min(), 'bar', 'spread of P at -2000 mRL', 1)} bar and at -2500 mRL within {dv('spread_2500', p2500.max() - p2500.min(), 'bar', 'spread of P at -2500 mRL', 1)} bar; W8 (cold) has the highest deep pressure.")
    add("")
    add("### 4b. Feed-zone pressures and pivots (`feedzones.csv`)")
    add("")
    add("Pivot rule: the dominant feed zone is the one with the highest quality factor, then the thicker zone, then the shallower zone; the pivot pressure is the survey pressure at its mid-point.")
    add("")
    rows = []
    for r in FZ.itertuples():
        piv = "pivot" if abs(r.mid_mRL - S.loc[r.well, "pivot_elev_mRL"]) < 1e-6 else ""
        rows.append((short(r.well), r.feedzone, f(r.top_mRL, 0), f(r.bottom_mRL, 0), r.quality, f(r.T_degC, 1), f(r.P_bar, 1), f(r.boiling_margin_degC, 1), r.block, f(r.weight_T, 2), f(r.weight_P, 2), piv))
    add(md_table(["well", "fz", "top mRL", "bottom mRL", "quality", "T degC", "P bar", "boiling margin degC", "block", "weight_T", "weight_P", ""], rows))
    add("")
    add(f"[observation] Feed-zone pressures agree with the trend of the passed wells within {f(tr_max, 2)} bar (`pressure_trend_check.csv`); the largest residuals are W1 ({f(TR[TR.well == 'WELL-1'].residual_bar.min(), 2)} bar, the only well below the trend by more than 2.5 bar and close to the {PRM['fz_trend_ok_bar']:.0f} bar limit) and W4 ({f(TR[TR.well == 'WELL-4'].residual_bar.iloc[0], 2)} bar, an injector).")
    add("")
    add("### 4c. Hydrostatic checks (`hydrostatic_check.csv`, `hydrostatic_summary.csv`)")
    add("")
    add(f"Windows are single-phase liquid (boiling margin >= {PRM['hydro_min_margin_degC']:.0f} degC and P >= {PRM['hydro_min_P_bar']} bar throughout), non-overlapping, at least {PRM['hydro_min_length_m']:.0f} m long, cut at feed-zone mid-points where the minimum length allows, and never straddling the casing shoe "
        f"(a piece shorter than {PRM['hydro_min_length_m']:.0f} m on one side of the shoe is not tested: {len(ctx['hnotes'])} such pieces, listed in `run_console_log.txt`). g = {dv('g', G, 'm/s2', 'gravity used', 5)} m/s2. "
        f"Each reading is rounded to 0.1 bar, so a window's pressure change is uncertain by up to 0.1 bar: at most {dv('max_round_ratio', HY.rounding_unc_ratio.max(), '', 'largest rounding uncertainty of a ratio', 3)} in the ratio and {dv('max_round_rho', HY.rounding_unc_rho_kgm3.max(), 'kg/m3', 'largest rounding uncertainty of implied density', 1)} kg/m3 in the implied density.")
    add("")
    rows = []
    for (scope, band), r in hs.iterrows():
        if scope in ("all wells", "cold wells (Tmax<100 degC)", "warm/hot wells"):
            rows.append((scope, band, int(r.n_windows), f(r.length_m, 0), f(r.ratio_median, 3), f(r.ratio_min, 3), f(r.ratio_max, 3), f(r.ratio_length_weighted, 3), "", ""))
    for n in wells:
        for band in ("cased", "open"):
            if (n, band) in hs.index:
                r = hs.loc[(n, band)]
                rows.append((short(n), band, int(r.n_windows), f(r.length_m, 0), f(r.ratio_median, 3), f(r.ratio_min, 3), f(r.ratio_max, 3), f(r.ratio_length_weighted, 3), f(r.cold_component, 3), f(r.excess_over_cold_component, 3)))
    add(md_table(["scope", "band", "windows", "length m", "median ratio", "min", "max", "length-weighted ratio", "cold-well component", "excess over it"], rows))
    add("")
    add(f"[observation] The field-wide component, the median ratio of the cold wells, is **{f(comp_c, 3)} (cased)** and **{f(comp_o, 3)} (open)**. The open-band component rests on W8 alone (W4's open piece is only 250 m long). "
        f"The implied density of the cold-well cased windows is {f(HY[(HY.band == 'cased') & (HY.cold_well)].rho_implied_kgm3.min(), 0)} to {f(HY[(HY.band == 'cased') & (HY.cold_well)].rho_implied_kgm3.max(), 0)} kg/m3 against about {dv('if97_density_cold_cased', HY[(HY.band == 'cased') & (HY.cold_well)].dP_IF97_bar.sum() * 1e5 / (G * HY[(HY.band == 'cased') & (HY.cold_well)].length_m.sum()), 'kg/m3', 'length-weighted IF97 density of the cold-well cased windows', 0)} kg/m3 for IF97 water at those temperatures. "
        f"Per-well excess above the component: cased {', '.join(f'{short(r.scope)} {f(r.excess_over_cold_component, 3)}' for r in HS[(HS.band == 'cased') & HS.scope.str.startswith('WELL')].itertuples())}; open {', '.join(f'{short(r.scope)} {f(r.excess_over_cold_component, 3)}' for r in HS[(HS.band == 'open') & HS.scope.str.startswith('WELL')].itertuples())}.")
    add("")
    add("**Checks before calling the gradient impossible.**")
    add("")
    add("- *Datum offset:* adds the same pressure at every depth, cannot change a gradient. The shallowest survey pressures are 1.0 bar (atmospheric) at W2, W7 and W10, so there is no sign of an absolute offset.")
    add(f"- *Gravity or unit convention:* would scale every gradient by one constant at every depth in every well. The open-hole ratios of W5, W6 and W7 are {f(hs.loc[('WELL-5', 'open'), 'ratio_length_weighted'], 3)}, {f(hs.loc[('WELL-6', 'open'), 'ratio_length_weighted'], 3)} and {f(hs.loc[('WELL-7', 'open'), 'ratio_length_weighted'], 3)} while the cased ratios of the same wells are {f(hs.loc[('WELL-5', 'cased'), 'ratio_length_weighted'], 3)}, {f(hs.loc[('WELL-6', 'cased'), 'ratio_length_weighted'], 3)} and {f(hs.loc[('WELL-7', 'cased'), 'ratio_length_weighted'], 3)}, so one constant cannot explain the data; it could at most explain the small cold-well component.")
    add("- *Two-phase fluid:* makes the column lighter, so it lowers the ratio; it cannot explain a ratio above 1.")
    add("- *Rounding:* at most the uncertainties quoted above; far smaller than the excess.")
    add("")
    add(f"[interpretation] What remains: (a) a field-wide component of {dv('comp_o_minus1', comp_o - 1, '', 'open-band cold-well component minus 1', 3)} to {dv('comp_c_minus1', comp_c - 1, '', 'cased-band cold-well component minus 1', 3)} that may be a convention or fluid-density difference (e.g. dissolved solids; not testable here), and (b) an excess of up to {dv('max_cased_minus1', hy_c.ratio_measured_to_IF97.max() - 1, '', 'largest cased ratio minus 1', 3)} in the hot cased sections that no listed convention explains. "
        "In the cased band it rises broadly with the shallow temperature gradient (W3 836, W6 490, W9 319, W7 301, W5 295 degC/km, against W1 73, W8 85, W4 14 degC/km) and it is nearly absent in the deep open hole of W5, W6, W7, W8 (ratios within about 1.5 per cent of 1). A static single-phase column cannot be denser than pure water at the measured temperature (ratio 1) apart from a few per cent for dissolved solids, so the shallow hot-well pressures are probably formation-type values (vertical flow, or values taken from a model) rather than static well-column pressures; "
        "the data cannot decide. A shut-in, thermally recovered pressure survey in W3 or W6 would decide.")
    add("")
    add("![Pressure map at -1000 mRL](figures/maps_pressure_minus1000m.png)")
    add("")

    # ================= 5 conceptual =================
    add("## 5. Conceptual model")
    add("")
    add("![Conceptual model](figures/conceptual_model.png)")
    add("")
    add("Sketch built from the well data only (`figures/conceptual_model.png`). Each feature lists its supporting observation, the relation to the mapped faults and surface water (`well_geometry_check.csv`), and what would decide an alternative.")
    add("")
    nearest = lambda n: f"{GT.loc[n, 'nearest_fault']} {f(GT.loc[n, 'dist_nearest_fault_m'], 0)} m; {GT.loc[n, 'nearest_water_body']} {f(GT.loc[n, 'dist_nearest_water_body_m'], 0)} m"
    add(md_table(["well", "nearest fault and distance", "nearest water body and distance", "wellhead minus model surface m", "max deviation m"],
                 [(short(n), nearest(n).split("; ")[0], nearest(n).split("; ")[1], f(GT.loc[n, "wellhead_minus_surface_m"], 0), f(GT.loc[n, "max_horizontal_deviation_m"], 0)) for n in wells]))
    add("")
    add(f"1. **Upflow** (interpretation, medium confidence). Supported by the highest temperatures at the shallowest depths (W10 {pw('WELL-10', 'Tmax_degC')}, W2 {pw('WELL-2', 'Tmax_degC')}, W3 {pw('WELL-3', 'Tmax_degC')} degC), the shallowest 200 degC isotherm (W10 {pw('WELL-10', 'z_200C_mRL')} mRL), a cap of only {f(S.loc[hot3, 'depth_100C_below_wellhead_m'].min(), 0)} to {f(S.loc[hot3, 'depth_100C_below_wellhead_m'].max(), 0)} m above 100 degC, "
        f"and a temperature at the bottom of W2 that is still rising ({pw('WELL-2', 'T_at_TD_degC')} degC at {pw('WELL-2', 'TD_mRL', 0)} mRL; deep gradient {pw('WELL-2', 'deep_gradient_degC_per_km')} degC/km). W2 and W10 lie {f(GT.loc['WELL-2', 'dist_nearest_fault_m'], 0)} m and {f(GT.loc['WELL-10', 'dist_nearest_fault_m'], 0)} m from Fault 3. "
        "*Alternatives:* upflow in the matrix between Fault 3 and Fault 1, or on Fault 1 (W3 at "
        f"{f(GT.loc['WELL-3', 'dist_Fault_1_m'], 0)} m, W6 at {f(GT.loc['WELL-6', 'dist_Fault_1_m'], 0)} m). The data cannot say whether the faults are permeable, because no well is closer than {f(GT.dist_nearest_fault_m.min(), 0)} m to one. *Decisive:* a deep temperature log in W10 (it stops at -350 mRL) and in W2 below -2348 mRL, and a fault-zone well.")
    add(f"2. **Lateral outflow** (interpretation). {lst(reg['outflow'])} show a temperature maximum at {f(min(S.loc[n, 'Tmax_elev_mRL'] for n in reg['outflow']), 0)} to {f(max(S.loc[n, 'Tmax_elev_mRL'] for n in reg['outflow']), 0)} mRL, deepening away from the core "
        f"(W3 {pw('WELL-3', 'Tmax_elev_mRL', 0)}, W2 {pw('WELL-2', 'Tmax_elev_mRL', 0)}, W6 {pw('WELL-6', 'Tmax_elev_mRL', 0)}, W7 {pw('WELL-7', 'Tmax_elev_mRL', 0)}, W5 {pw('WELL-5', 'Tmax_elev_mRL', 0)}, W9 {pw('WELL-9', 'Tmax_elev_mRL', 0)} mRL), with cooler water below ({f(S.loc['WELL-3', 'inversion_degC'], 1)} degC inversion at W3). "
        f"Directions: south to W6 and W7, east to W5, north-east to W9 (plan view arrows). The layer within {PRM['hot_layer_drop_degC']:.0f} degC of Tmax is {', '.join(f'{short(n)} {hz(n)} m' for n in reg['outflow'] if n != 'WELL-5')} thick (W5 is nearly isothermal over {hz('WELL-5')} m). "
        "*Alternative:* the inversions are cold water entering from below rather than hot water leaving laterally; temperature alone cannot separate them. *Decisive:* chemistry or isotopes along the outflow path, or a deeper temperature log in W1 and W5. W1 (the most distant hot well) shows no inversion, so the outflow does not reach it at the surveyed depths.")
    add(f"3. **Cap.** [interpretation] The conductive cap is {f(S.depth_100C_below_wellhead_m.min(), 0)} to {f(S.depth_100C_below_wellhead_m.max(), 0)} m thick to 100 degC and thickens away from the core; the near-surface temperature of {topT('WELL-4')} degC at the top of W4 matches ordinary ground temperature. The cap is the zone with steep gradients ({f(S.shallow_gradient_degC_per_km.min(), 0)} to {f(S.shallow_gradient_degC_per_km.max(), 0)} degC/km) "
        "and sits above the first feed zones (the shallowest feed-zone top is at -250 mRL in W2 and W10).")
    add(f"4. **Cold downflow / recharge** (interpretation; W4 and W8 differ in detail). W4 ({GT.loc['WELL-4', 'nearest_fault']} at {f(GT.loc['WELL-4', 'dist_nearest_fault_m'], 0)} m, the highest ground, wellhead {pw('WELL-4', 'wellhead_mRL', 0)} mRL) and W8 (Tarawera River at {f(GT.loc['WELL-8', 'dist_nearest_water_body_m'], 0)} m, wellhead {pw('WELL-8', 'wellhead_mRL', 0)} mRL) reach no more than {pw('WELL-4', 'Tmax_degC')} and {pw('WELL-8', 'Tmax_degC')} degC. "
        f"W8 is nearly isothermal at {f(ISO[ISO.well == 'WELL-8'].T_min_degC.iloc[0], 1)} to {f(ISO[ISO.well == 'WELL-8'].T_max_degC.iloc[0], 1)} degC from -550 to -1700 mRL, with a {pw('WELL-8', 'inversion_degC')} degC inversion and a rise to {pw('WELL-8', 'T_at_TD_degC')} degC at the bottom, which fits lateral cool inflow at that depth as well as downflow. *Alternative for W4:* a cold, low-permeability flank with conduction only (W4's gradient is {pw('WELL-4', 'shallow_gradient_degC_per_km')} degC/km, far below the conductive values of the other wells). *Decisive:* a temperature log of W4 below -750 mRL, or the chemistry of its water.")
    add(f"5. **Two-phase zone** (interpretation, low confidence). Only W10, in the depth range {f(BF.z_top_mRL.iloc[0], 0)} to {f(BF.z_bottom_mRL.iloc[0], 1)} mRL, with the pressure caveat of section 4c.")
    add(f"6. **Pressure connection** (interpretation). Pressures at -1500 mRL agree within {f(p1500.max() - p1500.min(), 1)} bar and the feed-zone pressures lie on one trend (residuals <= {f(tr_max, 1)} bar), so the feed zones, from -250 to -1920 mRL, appear to be in one pressure system; the shallow water levels differ by up to {f(S.water_level_mRL.max() - S.water_level_mRL.min(), 0)} m, broadly following the model surface elevation ({f(S.model_surface_mRL.min(), 0)} to {f(S.model_surface_mRL.max(), 0)} mRL).")
    add("7. **Surface water.** [observation] The Tarawera River is within a kilometre of W3, W6, W7, W8 and the other 'river' is within 350 m of W5 and about 100 m of W1; Lake Tarawera is more than three kilometres from every well. [interpretation] The data do not show river-water cooling in the hot wells (W3 and W6 are hot). Only W8, next to the river, is cold, but W4, far from any water, is cold too, so surface water is not the only cause.")
    add("")

    # ================= 6 set-up =================
    add("## 6. Recommended calibration set-up")
    add("")
    add("### 6a. Targets and weights")
    add("")
    WBW = tb["wbw"].set_index("well")
    rows = [(short(n), f(WBW.loc[n, "sum_T_weights"], 2), f(WBW.loc[n, "mean_T_weight"], 3), f(WBW.loc[n, "sum_P_weights"], 2),
             WBW.loc[n, "ns_confidence"], int(WBW.loc[n, "n_points"]), int(WBW.loc[n, "n_points_in_feedzone"])) for n in wells]
    add(md_table(["well", "sum of T weights", "mean T weight", "sum of P weights", "confidence", "points", "points in feed zones"], rows))
    add("")
    add(f"- Temperature: all {len(TG)} survey points (`NS_calibration_targets.csv`) with location weight times confidence factor: feed zone 1.0, open hole 0.8, open hole in an isothermal section >= 300 m 0.5, behind casing 0.5, within 30 m of the wellhead or in the top block of the column 0.3, above the model surface 0 (no point is above the surface because every wellhead is at or below it). "
        f"The temperature weights add to {dv('sum_T_weights_all', TG.weight.sum(), '', 'sum of all temperature weights', 2)}.")
    add(f"- Pressure: only the {int((TG.weight_P > 0).sum())} survey points inside feed zones have `weight_P` (quality factor times confidence); {len(FZ)} feed-zone mid-point pressures are in `feedzones.csv`. Use the pivot pressures as the primary constraint on each well. No shallow wellbore pressure is a target.")
    add("- Pressure tolerance: use several bar (section 4c) and treat the cased-band ratios as a data question, not as overpressure.")
    add("- Layer-centre values for comparison with block temperatures and pressures are in `T_P_at_layer_centres.csv` (no extrapolation; blocks from the well-to-block file).")
    add("")
    add("### 6b. Heat and mass")
    add("")
    add(f"- **Heat in:** a deep source under the W10-W2-W3 area able to supply liquid at {f(S.loc['WELL-2', 'T_at_TD_degC'], 0)} to {f(S.loc['WELL-10', 'Tmax_degC'], 0)} degC [interpretation]. The margins need a low basal heat input: with the assumed conductivity the W1, W4 and W8 shallow gradients are {pw('WELL-1', 'shallow_heat_flux_mW_m2_assumed_k')}, {pw('WELL-4', 'shallow_heat_flux_mW_m2_assumed_k')} and {pw('WELL-8', 'shallow_heat_flux_mW_m2_assumed_k')} mW/m2 [assumption].")
    ref = ctx["REF"]
    add(f"- **[reference model]** (for comparison only): uniform basal heat flux {ref.get('heat_flux_default_W_m2', float('nan')):.1f} W/m2 in the boundary-condition file, {int(ref.get('n_heat_sources', 0))} base heat sources totalling {f(ref.get('heat_total_MW', float('nan')), 1)} MW "
        f"({f(ref.get('heat_per_block_min_kW', float('nan')), 1)} to {f(ref.get('heat_per_block_max_kW', float('nan')), 1)} kW per block), atmosphere {f(ref.get('atmosphere_T_degC', float('nan')), 1)} degC, rain {f(ref.get('rain_default_mm_yr', float('nan')), 0)} mm/year. "
        "These are the starting values of the model, not data. The shallow gradients of the margin wells ({pw('WELL-4', 'shallow_heat_flux_mW_m2_assumed_k', 0)} to {pw('WELL-8', 'shallow_heat_flux_mW_m2_assumed_k', 0)} mW/m2 with the assumed conductivity) bracket the uniform value, so the margins do not demand a change; the hot core needs a much larger, localised input.")
    add("- **Mass out:** lateral outflow from the hot core at about -450 to -750 mRL toward the south, east and north-east. **Mass in:** meteoric recharge, the cold flanks W4 and W8, and a deep hot upflow; the mass rates are not constrained by the natural-state data.")
    add("")
    add("### 6c. Permeability structure implied")
    add("")
    add("- A low-permeability cap over the hot core, 60 to 100 m thick to 100 degC and thickening outward [interpretation].")
    shal = [n for n in wells if (FZ[FZ.well == n].mid_mRL >= -900).any()]
    deep = [n for n in wells if ((FZ[FZ.well == n].mid_mRL <= -1100) & (FZ[FZ.well == n].mid_mRL > -1800)).any()]
    fs, fd = FZ[FZ.mid_mRL >= -900], FZ[(FZ.mid_mRL <= -1100) & (FZ.mid_mRL > -1800)]
    add(f"- Permeable feed zones: shallow ones at {f(fs.top_mRL.max(), 0)} to {f(fs.bottom_mRL.min(), 0)} mRL ({andlist(shal)}) and deep ones at {f(fd.top_mRL.max(), 0)} to {f(fd.bottom_mRL.min(), 0)} mRL ({andlist(deep)}), plus the W2 zone at {f(FZ[FZ.mid_mRL <= -1800].mid_mRL.iloc[0], 0)} mRL. Lateral permeability at these elevations is what the model must provide (the elevations are observations; their permeability is an inference).")
    add("- Vertical permeability high enough for upflow under the hot core and for cold downflow at the NE and SW margins [interpretation].")
    add("- High lateral permeability at depth, from the common deep pressure (conclusion 8).")
    add("")
    add("### 6d. Boundary requirements")
    add("")
    add(f"- Top: atmospheric pressure and ground temperature near 12 degC (the shallowest temperatures are {topT('WELL-4')} degC at W4 and {topT('WELL-8')} degC at W8); water table {f(S.water_level_below_model_surface_m.min(), 0)} to {f(S.water_level_below_model_surface_m.max(), 0)} m below the model surface.")
    add("- Sides: open to cold recharge at the NE and SW margins; no lateral pressure barrier at depth within the well field.")
    add("- Base: concentrated heat and mass input under W10/W2/W3; low heat input on the margins.")
    add("")

    # ================= 7 data quality =================
    add("## 7. Data quality issues and open questions (ranked by effect on calibration)")
    add("")
    iss = [
        ("Pressure gradients exceed the IF97 column", f"Cold-well component {f(comp_c, 3)} (cased) and {f(comp_o, 3)} (open); hot cased windows up to {f(hy_c.ratio_measured_to_IF97.max(), 3)}. Affects every pressure comparison and the boiling margins of the hot wells. Ask the data provider how the pressures were measured or produced; a shut-in survey in W3 or W6 would settle it."),
        ("Date label versus pre-production match", f"Surveys are labelled 29/07/2015. W2 and W3 match the pre-production gauge within 0.1 bar, yet the surveys of W2, W3 and W9 sit {f(PC.loc['WELL-3', 'survey_minus_gauge_at_survey_date_bar'], 1)} to {f(PC.loc['WELL-9', 'survey_minus_gauge_at_survey_date_bar'], 1)} bar above the gauge interpolated to the labelled date, and W10's first flow is years after the label. The natural-state use does not depend on the label; confirm the survey date."),
        ("W10 is the hottest well and the most uncertain", f"Its survey stops at {f(S.loc['WELL-10', 'TD_mRL'], 0)} mRL, its feed-zone tag is 'Synthetic', its gauge starts {int(PC.loc['WELL-10', 'gauge_start_minus_field_first_production_days'])} days after field production began, and its margin to boiling is {f(S.loc['WELL-10', 'min_boiling_margin_degC'], 1)} degC. Weighted as medium confidence and quality factor 0.5."),
        ("W9 survey above its first gauge value", f"{f(PC.loc['WELL-9', 'diff_survey_minus_first_bar'], 2)} bar; verdict probably natural state, medium."),
        ("Wellhead below the model surface", f"{', '.join(f'{short(n)} {f(-GT.loc[n, 'wellhead_minus_surface_m'], 0)} m' for n in wh_below.index)}. The block file starts at the wellhead, so the model cells above it have no data; the top-block weight rule uses the model surface."),
        ("Extrapolated water levels", f"W1 {f(S.loc['WELL-1', 'water_level_extrapolation_m'], 0)} m and W5 {f(S.loc['WELL-5', 'water_level_extrapolation_m'], 0)} m above the top survey point; W1 has {f(TG[TG.well == 'WELL-1'].sort_values('elevation_mRL').P_bar.iloc[-1], 1)} bar at its top point."),
        ("Isothermal open-hole sections", "Up to 1150 m long (W3, W8); wellbore flow versus formation convection not decidable from one survey. Weighted 0.5."),
        ("W4 design versus track", f"Liner bottom {f(w4d, 0)} mRL against the track bottom {f(S.loc['WELL-4', 'TD_mRL'], 0)} mRL; the track (and the grid file) were used."),
        ("Injector surveys", "W5's formation is hotter than the injectate, so injection cooling could be hidden (confidence low); W4 and W8 show no injection signature."),
        ("W1 feed-zone residual", f"{f(TR[TR.well == 'WELL-1'].residual_bar.min(), 2)} bar against the trend, close to the {PRM['fz_trend_ok_bar']:.0f} bar limit; lateral gradient or minor disturbance cannot be separated."),
        ("Empty provenance fields and irregular names", "Source/Comment fields are empty, StableFlag is True everywhere, and {len(irregular)} plain-column files use the prefix 'Wel10'; no effect on the numbers."),
    ]
    for i, (a, b) in enumerate(iss, 1):
        add(f"{i}. **{a}.** {b}")
    add("")
    add("**Open questions that observations could answer:** a shut-in pressure and temperature survey (pressure convention, wellbore effects); a deep log in W10 and below -2348 mRL in W2 (upflow depth and temperature); a W4 log below -750 mRL (cold downflow or conductive flank); a fault-crossing well (permeability of Faults 1 to 3).")
    add("")

    # ================= 8 index =================
    add("## 8. Output file index")
    add("")
    idx = [
        ("Scripts/natural_state_analysis.py", "regenerates everything, including this report and the input inventory"),
        ("Output/input_inventory.csv", "role of every file under Baseline/"),
        ("Output/NS_calibration_targets.csv", "one row per survey point: block, weights, weight_P"),
        ("Output/feedzones.csv", "T, P, margin, block, weights at each feed zone"),
        ("Output/T_P_at_layer_centres.csv", "interpolated T, P at model layer centres"),
        ("Output/well_summary.csv", "per-well characteristics, regime, verdict; extra columns for heat flux, cap depth, hot layer"),
        ("Output/pressure_consistency_check.csv", "natural-state gauge test and reasons"),
        ("Output/pressure_trend_check.csv", "feed-zone pressure against the trend of the passed wells"),
        ("Output/hydrostatic_check.csv, hydrostatic_summary.csv", "windows and ratio statistics by band"),
        ("Output/well_geometry_check.csv", "wellhead, surface, deviation, distances to faults and water"),
        ("Output/isothermal_sections.csv, boiling_margin_flags.csv", "isothermal sections and intervals with margin < 5 degC"),
        ("Output/TP_at_common_elevations.csv, temperature_at_common_elevations.csv, pressure_at_common_elevations.csv", "interpolated values at common elevations"),
        ("Output/dat_json_crosscheck.csv, integrity_checks.csv, provenance_labels.csv", "source cross-checks and labels"),
        ("Output/if97_validation.csv, analysis_parameters.csv, reference_model_summary.csv, report_derived_numbers.csv", "validation, parameters, reference-model description, numbers derived for this report"),
        ("Output/run_console_log.txt, acceptance_check.txt, process_log.md", "console log, acceptance check, process log"),
        ("Output/figures/well_map.png", "grid, faults, water bodies, wells"),
        ("Output/figures/profile_WELL-*.png", "T and P with Tsat, casing, feed zones, water level, pivot (one per well)"),
        ("Output/figures/maps_temperature_*.png, maps_pressure_*.png", "maps at 0, -500, -1000, -1500 mRL"),
        ("Output/figures/section_hottest_wells.png, section_pressure_hottest_wells.png", "sections through the four hottest wells"),
        ("Output/figures/PT_state.png", "P-T state against the saturation curve"),
        ("Output/figures/conceptual_model.png", "conceptual sketch"),
    ]
    add(md_table(["file", "content"], idx))
    add("")
    return "\n".join(L) + "\n", derived


if __name__ == "__main__":
    main()
