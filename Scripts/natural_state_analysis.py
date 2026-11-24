#!/usr/bin/env python
"""
Natural-state data interpretation for geothermal model calibration (TE14515).

Run from the repository root:
    python Scripts/natural_state_analysis.py

Reads only from Baseline/ (read-only). Writes every table, figure and the main
report to Output/. Rules followed are in CLAUDE.md.

Data used as evidence : Data/Data_json/WELL-*_data.json, Data/Data_dat/*.dat (cross-check),
                        Data/Mapping/Mapping.json, Data/Groups_data.json
Geometry only         : Models/TE14515/gTE14515.dat (VERTICES, GRID, LAYERS, SURFA, WELLS),
                        Models/TE14515/gTE14515_wells.json
Reference model only  : gTE14515_boundary_conditions.json, NaturalState/TE14515_NS_000.dat
                        (described, never used as evidence)
"""
import ast
import json
import math
import os
import sys
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.patches import Polygon as MplPolygon, FancyArrowPatch
from scipy.interpolate import griddata

# ----------------------------------------------------------------------------
# Paths and console logging
# ----------------------------------------------------------------------------
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BASE = os.path.join(ROOT, "Baseline", "Teaching_Model_3D")
DATA = os.path.join(BASE, "Data")
MODEL = os.path.join(BASE, "Models", "TE14515")
OUT = os.path.join(ROOT, "Output")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)


class Tee:
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
        self.stdout = sys.stdout

    def write(self, s):
        self.stdout.write(s)
        self.f.write(s)

    def flush(self):
        self.stdout.flush()
        self.f.flush()


sys.stdout = Tee(os.path.join(OUT, "run_console_log.txt"))
print(f"natural_state_analysis.py run at {datetime.now().isoformat(timespec='seconds')}")
print(f"Repository root: {ROOT}")

G = 9.80665  # m/s2, standard gravity (used for all hydrostatic calculations)
K_ROCK = 2.0  # W/m/K, ASSUMED rock thermal conductivity (not in data; typical value
# for altered volcaniclastic/ignimbrite sequences, e.g. Grant & Bixley 2011 range 1.5-2.5)
BM_FLAG = 5.0  # degC, boiling-margin flag threshold
WELLS = [f"WELL-{i}" for i in range(1, 11)]
DAT_NAME = {w: ("Wel10" if w == "WELL-10" else "Well" + w.split("-")[1]) for w in WELLS}

# ----------------------------------------------------------------------------
# 1. IAPWS-IF97 (Region 4 saturation, Region 1 compressed liquid)
# ----------------------------------------------------------------------------
_N4 = [0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
       0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
       -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
       0.65017534844798e3]
PC_MPA = 22.064


def psat_mpa(T_K):
    n = _N4
    th = T_K + n[8] / (T_K - n[9])
    A = th * th + n[0] * th + n[1]
    B = n[2] * th * th + n[3] * th + n[4]
    C = n[5] * th * th + n[6] * th + n[7]
    return (2 * C / (-B + math.sqrt(B * B - 4 * A * C))) ** 4


def tsat_K(p_mpa):
    if not (611.213e-6 <= p_mpa <= PC_MPA):
        return float("nan")
    n = _N4
    b = p_mpa ** 0.25
    E = b * b + n[2] * b + n[5]
    F = n[0] * b * b + n[3] * b + n[6]
    Gg = n[1] * b * b + n[4] * b + n[7]
    D = 2 * Gg / (-F - math.sqrt(F * F - 4 * E * Gg))
    return (n[9] + D - math.sqrt((n[9] + D) ** 2 - 4 * (n[8] + n[9] * D))) / 2


def tsat_c(p_bar):
    return tsat_K(p_bar / 10.0) - 273.15


_I1 = [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8,
       21, 23, 29, 30, 31, 32]
_J1 = [-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2,
       10, -8, -11, -6, -29, -31, -38, -39, -40, -41]
_n1 = [0.14632971213167, -0.84548187169114, -0.37563603672040e1, 0.33855169168385e1,
       -0.95791963387872, 0.15772038513228, -0.16616417199501e-1, 0.81214629983568e-3,
       0.28319080123804e-3, -0.60706301565874e-3, -0.18990068218419e-1,
       -0.32529748770505e-1, -0.21841717175414e-1, -0.52838357969930e-4,
       -0.47184321073267e-3, -0.30001780793026e-3, 0.47661393906987e-4,
       -0.44141845330846e-5, -0.72694996297594e-15, -0.31679644845054e-4,
       -0.28270797985312e-5, -0.85205128120103e-9, -0.22425281908000e-5,
       -0.65171222895601e-6, -0.14341729937924e-12, -0.40516996860117e-6,
       -0.12734301741641e-8, -0.17424871230634e-9, -0.68762131295531e-18,
       0.14478307828521e-19, 0.26335781662795e-22, -0.11947622640071e-22,
       0.18228094581404e-23, -0.93537087292458e-25]
R_KJ = 0.461526


def v_region1(T_K, p_mpa):
    pi = p_mpa / 16.53
    tau = 1386.0 / T_K
    gp = 0.0
    for I, J, n in zip(_I1, _J1, _n1):
        gp += -n * I * (7.1 - pi) ** (I - 1) * (tau - 1.222) ** J
    return R_KJ * T_K / (p_mpa * 1000.0) * pi * gp


def h_region1(T_K, p_mpa):
    """Specific enthalpy (kJ/kg), IF97 Region 1."""
    pi = p_mpa / 16.53
    tau = 1386.0 / T_K
    gt = 0.0
    for I, J, n in zip(_I1, _J1, _n1):
        gt += n * (7.1 - pi) ** I * J * (tau - 1.222) ** (J - 1)
    return R_KJ * T_K * tau * gt


def h_liquid(T_c, p_bar):
    return h_region1(T_c + 273.15, p_bar / 10.0)


def rho_liquid(T_c, p_bar):
    return 1.0 / v_region1(T_c + 273.15, p_bar / 10.0)


def validate_if97():
    print("\n=== IAPWS-IF97 validation against published values (IAPWS R7-97 tables) ===")
    checks = []
    # Tsat (K) at p (MPa): IF97 Table 35 + CLAUDE.md values
    for p, ref in [(0.1, 372.755919), (1.0, 453.035632), (10.0, 584.149488),
                   (0.101325, 373.124)]:
        val = tsat_K(p)
        checks.append(("Tsat", f"p={p} MPa", val, ref, abs(val - ref), 0.1))
    # CLAUDE.md degC values
    for pbar, ref in [(1.01325, 99.97), (100.0, 311.0)]:
        val = tsat_c(pbar)
        checks.append(("Tsat", f"p={pbar} bar", val, ref, abs(val - ref), 0.1))
    # Psat (MPa) at T (K): Table 33 -> convert to Tsat error equivalence via inverse
    for T, ref in [(300.0, 0.353658941e-2), (500.0, 0.263889776e1), (600.0, 0.123443146e2)]:
        val = psat_mpa(T)
        terr = abs(tsat_K(val) - T)  # round-trip
        checks.append(("Psat", f"T={T} K", val, ref, abs(val - ref) / ref, 1e-6))
        checks.append(("Tsat(Psat)", f"T={T} K", tsat_K(val), T, terr, 0.1))
    # Region 1 specific volume (m3/kg): Table 5
    for T, p, ref in [(300.0, 3.0, 0.100215168e-2), (300.0, 80.0, 0.971180894e-3),
                      (500.0, 3.0, 0.120241800e-2)]:
        val = v_region1(T, p)
        checks.append(("v_R1", f"T={T} K p={p} MPa", val, ref, abs(val - ref) / ref, 1e-6))
    for T, p, ref in [(300.0, 3.0, 0.115331273e3), (300.0, 80.0, 0.184142828e3), (500.0, 3.0, 0.975542239e3)]:
        val = h_region1(T, p)
        checks.append(("h_R1", f"T={T} K p={p} MPa", val, ref, abs(val - ref) / ref, 1e-6))
    ok = True
    for c in checks:
        status = "OK" if c[4] <= c[5] else "FAIL"
        if status == "FAIL":
            ok = False
        print(f"  {c[0]:10s} {c[1]:22s} calc={c[2]:.9g} ref={c[3]:.9g} err={c[4]:.3g} tol={c[5]} {status}")
    if not ok:
        raise SystemExit("IAPWS-IF97 validation FAILED (error above tolerance). Stopping.")
    print("  All IF97 checks passed (Tsat errors < 0.1 K, v, h and Psat relative errors < 1e-6).")
    return checks


IF97_CHECKS = validate_if97()

# ----------------------------------------------------------------------------
# 2. Load well data
# ----------------------------------------------------------------------------
def parse_list(s):
    return ast.literal_eval(s) if isinstance(s, str) else s


def load_wells():
    wells = {}
    for w in WELLS:
        d = json.load(open(os.path.join(DATA, "Data_json", f"{w}_data.json")))
        rec = {"name": w, "raw": d}
        rec["track"] = np.array([parse_list(s) for s in d["Well_Track"]["Data"]], float)
        rec["design"] = [parse_list(s) for s in d["Well_Design"]["Data"]]
        fz = []
        for s in d["Feed_Zones"]["Data"]:
            t, b, q = parse_list(s)
            fz.append({"top": float(t), "bot": float(b), "quality": str(q)})
        rec["fz"] = fz
        (tdate, tblk), = d["Downhole_Temperature"].items()
        (pdate, pblk), = d["Downhole_Pressure"].items()
        zT = np.array([float(k) for k in tblk["Data"].keys()])
        T = np.array(list(tblk["Data"].values()), float)
        zP = np.array([float(k) for k in pblk["Data"].keys()])
        P = np.array(list(pblk["Data"].values()), float)
        rec.update(survey_date_T=tdate, survey_date_P=pdate, zT=zT, T=T, zP=zP, P=P,
                   stable_T=tblk.get("StableFlag"), stable_P=pblk.get("StableFlag"))
        rec["transient_P"] = {}
        for el, blk in d.get("Transient_Pressure", {}).items():
            dates = list(blk["Data"].keys())
            rec["transient_P"][float(el)] = (dates, np.array(list(blk["Data"].values()), float))
        rec["prod"] = d.get("Production_Mass_Flow", {}).get("Data", {})
        rec["inj"] = d.get("Injection_Mass_Flow", {}).get("Data", {})
        rec["inj_h"] = d.get("Injection_Enthalpy", {}).get("Data", {})
        rec["prod_h"] = d.get("Transient_Enthalpy", {}).get("Data", {})
        casing = [r for r in rec["design"] if not r[3]]
        liner = [r for r in rec["design"] if r[3]]
        rec["casing_shoe"] = min(float(r[1]) for r in casing) if casing else None
        rec["casing_top"] = max(float(r[0]) for r in casing) if casing else None
        rec["liner_top"] = max(float(r[0]) for r in liner) if liner else None
        rec["liner_bot"] = min(float(r[1]) for r in liner) if liner else None
        rec["wellhead"] = rec["track"][0, 2]
        rec["td"] = rec["track"][-1, 2]
        wells[w] = rec
    return wells


W = load_wells()
groups = json.load(open(os.path.join(DATA, "Groups_data.json")))
group_list = groups["Areas"]["Synthetic"]["All"]
print(f"\nLoaded {len(W)} wells; Groups_data.json lists {len(group_list)} wells: "
      f"{'match' if sorted(group_list) == sorted(WELLS) else 'MISMATCH'}")


def track_xy(rec, z):
    """x, y on the well track at elevation z (track z decreasing)."""
    tr = rec["track"]
    zz = tr[::-1, 2]
    return (float(np.interp(z, zz, tr[::-1, 0])), float(np.interp(z, zz, tr[::-1, 1])))


def interp_profile(zs, vals, z):
    """Linear interpolation in elevation; NaN outside the surveyed range (no extrapolation)."""
    o = np.argsort(zs)
    zs, vals = np.asarray(zs)[o], np.asarray(vals)[o]
    z = np.atleast_1d(z).astype(float)
    out = np.interp(z, zs, vals)
    out[(z < zs[0] - 1e-6) | (z > zs[-1] + 1e-6)] = np.nan
    return out


# ----------------------------------------------------------------------------
# 3. Model geometry (geometry sections only)
# ----------------------------------------------------------------------------
def load_geometry():
    lines = open(os.path.join(MODEL, "gTE14515.dat")).read().splitlines()
    sec = None
    verts, cols, layers, surf, gwells = {}, {}, [], {}, {}
    cur = None
    for ln in lines:
        key = ln.strip()
        if key in ("VERTICES", "GRID", "CONNECTIONS", "LAYERS", "SURFA", "WELLS"):
            sec = key
            continue
        if ln.startswith("GENER"):
            continue
        if not key:
            continue
        if sec == "VERTICES":
            verts[ln[0:3]] = (float(ln[3:13]), float(ln[13:23]))
        elif sec == "GRID":
            if len(ln.rstrip()) > 3 and ln[3:].strip():
                cur = ln[0:3]
                cols[cur] = []
            else:
                cols[cur].append(ln[0:3])
        elif sec == "LAYERS":
            name = ln[0:2]
            parts = ln[2:].split()
            layers.append((name, float(parts[0]), float(parts[1])))
        elif sec == "SURFA":
            surf[ln[0:3]] = float(ln[3:].split()[0])
        elif sec == "WELLS":
            nm = ln[0:5]
            x, y, z = map(float, ln[5:].split())
            gwells.setdefault(nm, []).append((x, y, z))
    # layer tops
    lay = []
    top_surface = layers[0][1]
    prev_bot = top_surface
    for name, bot, cen in layers[1:]:
        lay.append({"name": name, "top": prev_bot, "bot": bot, "centre": cen})
        prev_bot = bot
    polys = {c: np.array([verts[v] for v in vs]) for c, vs in cols.items()}
    return verts, polys, lay, surf, gwells, top_surface


VERTS, COLS, LAYERS, SURF, GWELLS, TOPSURF = load_geometry()
print(f"Geometry: {len(VERTS)} vertices, {len(COLS)} columns, {len(LAYERS)} layers "
      f"(top {LAYERS[0]['top']} to bottom {LAYERS[-1]['bot']} mRL), "
      f"{len(SURF)} SURFA entries, {len(GWELLS)} wells in WELLS section")


def point_in_poly(x, y, poly):
    n = len(poly)
    inside = False
    j = n - 1
    for i in range(n):
        xi, yi = poly[i]
        xj, yj = poly[j]
        if ((yi > y) != (yj > y)) and (x < (xj - xi) * (y - yi) / (yj - yi + 1e-30) + xi):
            inside = not inside
        j = i
    return inside


def find_column(x, y):
    for c, p in COLS.items():
        if point_in_poly(x, y, p):
            return c
    return None


def col_surface(c):
    return SURF.get(c, TOPSURF)


def find_layer(z):
    for L in LAYERS:
        if L["bot"] < z <= L["top"] + 1e-9:
            return L
    if abs(z - LAYERS[-1]["bot"]) < 1e-6:
        return LAYERS[-1]
    return None


def block_of(x, y, z):
    c = find_column(x, y)
    if c is None:
        return None, None, None, None
    s = col_surface(c)
    if z > s + 1e-6:
        return c, None, None, True
    L = find_layer(z)
    if L is None:
        return c, None, None, False
    return c, L["name"].strip(), c + L["name"], False


def top_layer_of_column(c):
    s = col_surface(c)
    for L in LAYERS:
        if L["bot"] < s <= L["top"] + 1e-9:
            return L
    return None


# ----------------------------------------------------------------------------
# 4. Mapping
# ----------------------------------------------------------------------------
MAP = json.load(open(os.path.join(DATA, "Mapping", "Mapping.json")))
FAULTS = {k: np.array(v["Data"], float) for k, v in MAP.items() if k.startswith("Fault")}
LAKE = np.array(MAP["Lake Tarawera"]["Data"], float)
RIVERS = {k: [np.array(seg, float) for seg in MAP[k]["Data"]] for k in ["Tarawera River", "river"]}  # multiplot


def dist_point_polyline(px, py, pl):
    d = np.inf
    for (x1, y1), (x2, y2) in zip(pl[:-1], pl[1:]):
        dx, dy = x2 - x1, y2 - y1
        L2 = dx * dx + dy * dy
        t = 0 if L2 == 0 else max(0, min(1, ((px - x1) * dx + (py - y1) * dy) / L2))
        d = min(d, math.hypot(px - (x1 + t * dx), py - (y1 + t * dy)))
    return d


def seg_intersect(p1, p2, p3, p4):
    d = (p2[0] - p1[0]) * (p4[1] - p3[1]) - (p2[1] - p1[1]) * (p4[0] - p3[0])
    if abs(d) < 1e-12:
        return None
    t = ((p3[0] - p1[0]) * (p4[1] - p3[1]) - (p3[1] - p1[1]) * (p4[0] - p3[0])) / d
    u = ((p3[0] - p1[0]) * (p2[1] - p1[1]) - (p3[1] - p1[1]) * (p2[0] - p1[0])) / d
    if 0 <= t <= 1 and 0 <= u <= 1:
        return (p1[0] + t * (p2[0] - p1[0]), p1[1] + t * (p2[1] - p1[1]))
    return None


FAULT_X = {}
fk = list(FAULTS)
for i in range(len(fk)):
    for j in range(i + 1, len(fk)):
        p = seg_intersect(FAULTS[fk[i]][0], FAULTS[fk[i]][1], FAULTS[fk[j]][0], FAULTS[fk[j]][1])
        if p is not None:
            FAULT_X[f"{fk[i]} x {fk[j]}"] = p
print("Fault intersections (from Mapping.json traces):",
      {k: (round(v[0]), round(v[1])) for k, v in FAULT_X.items()})

# ============================================================================
# 5. Data integrity: JSON vs dat; points within track / slotted interval
# ============================================================================
print("\n=== 3.1 Data integrity ===")
integ_rows = []
for w, r in W.items():
    for kind, zs, vs in [("temperature", r["zT"], r["T"]), ("pressure", r["zP"], r["P"])]:
        fn = os.path.join(DATA, "Data_dat", f"{DAT_NAME[w]}_downhole_{kind}.dat")
        dat = np.loadtxt(fn)
        n_json, n_dat = len(zs), len(dat)
        dz = dv = np.nan
        if n_json == n_dat:
            dv = float(np.max(np.abs(dat[:, 0] - vs)))
            dz = float(np.max(np.abs(dat[:, 1] - zs)))
        integ_rows.append({"well": w, "check": f"json_vs_dat_{kind}", "n_json": n_json,
                           "n_dat": n_dat, "max_abs_diff_value": dv, "max_abs_diff_elev": dz,
                           "result": "OK" if (n_json == n_dat and dv < 1e-6 and dz < 1e-6) else "MISMATCH",
                           "detail": os.path.basename(fn)})
    for el, (dates, vals) in r["transient_P"].items():
        fn = os.path.join(DATA, "Data_dat", f"{DAT_NAME[w]}_{el}_transient_pressure.dat")
        dat = np.loadtxt(fn)
        dv = float(np.max(np.abs(dat[:, 1] - vals))) if len(dat) == len(vals) else np.nan
        integ_rows.append({"well": w, "check": "json_vs_dat_transient_pressure",
                           "n_json": len(vals), "n_dat": len(dat), "max_abs_diff_value": dv,
                           "max_abs_diff_elev": np.nan,
                           "result": "OK" if dv == 0 else "MISMATCH",
                           "detail": f"{os.path.basename(fn)}; dat time column starts at "
                                     f"{dat[0, 0]:g} (years after 2014-01-01 inferred), first date {dates[0]}"})
    # same elevations in T and P
    same = len(r["zT"]) == len(r["zP"]) and np.allclose(r["zT"], r["zP"])
    integ_rows.append({"well": w, "check": "T_and_P_elevations_identical", "n_json": len(r["zT"]),
                       "n_dat": len(r["zP"]), "max_abs_diff_value": np.nan, "max_abs_diff_elev": np.nan,
                       "result": "OK" if same else "MISMATCH", "detail": ""})
    wh, td = r["wellhead"], r["td"]
    zmax, zmin = r["zT"].max(), r["zT"].min()
    inside = (zmax <= wh + 1e-6) and (zmin >= td - 1e-6)
    integ_rows.append({"well": w, "check": "survey_within_track", "n_json": len(r["zT"]),
                       "n_dat": np.nan, "max_abs_diff_value": np.nan, "max_abs_diff_elev": np.nan,
                       "result": "OK" if inside else "OUTSIDE",
                       "detail": f"track {wh:g} to {td:g} mRL; survey {zmax:g} to {zmin:g} mRL; "
                                 f"first survey point {wh - zmax:g} m below wellhead"})
    # design vs track
    dtop = max(float(d[0]) for d in r["design"])
    dbot = min(float(d[1]) for d in r["design"])
    ok = abs(dtop - wh) < 1e-6 and abs(dbot - td) < 1e-6
    integ_rows.append({"well": w, "check": "design_matches_track", "n_json": np.nan, "n_dat": np.nan,
                       "max_abs_diff_value": np.nan, "max_abs_diff_elev": np.nan,
                       "result": "OK" if ok else "INCONSISTENT",
                       "detail": f"design {dtop:g} to {dbot:g} mRL vs track {wh:g} to {td:g} mRL"})
    for k, f in enumerate(r["fz"]):
        in_track = f["top"] <= wh and f["bot"] >= td
        in_slot = (r["liner_top"] is not None and f["top"] <= r["liner_top"] and f["bot"] >= r["liner_bot"])
        below_shoe = f["top"] <= r["casing_shoe"]
        res = "OK" if (in_track and in_slot and below_shoe) else "PROBLEM"
        integ_rows.append({"well": w, "check": f"feedzone_{k + 1}_location", "n_json": np.nan,
                           "n_dat": np.nan, "max_abs_diff_value": np.nan, "max_abs_diff_elev": np.nan,
                           "result": res,
                           "detail": f"FZ {f['top']:g} to {f['bot']:g} mRL ({f['quality']}); in track={in_track}; "
                                     f"in slotted liner ({r['liner_top']:g} to {r['liner_bot']:g})={in_slot}; "
                                     f"below casing shoe ({r['casing_shoe']:g})={below_shoe}"})
    # geometry file WELLS section vs json track
    gname = r["raw"]["AUT2_Name"] if r["raw"]["AUT2_Name"] in GWELLS else DAT_NAME[w].upper()[:5]
    gname = {"WELL-10": "WEL10"}.get(w, f"WELL{w.split('-')[1]}")
    gt = np.array(GWELLS.get(gname, []))
    if len(gt):
        dmax = float(np.max(np.abs(gt - r["track"][:len(gt)]))) if len(gt) == len(r["track"]) else np.nan
        integ_rows.append({"well": w, "check": "geometry_WELLS_vs_json_track", "n_json": len(r["track"]),
                           "n_dat": len(gt), "max_abs_diff_value": dmax, "max_abs_diff_elev": np.nan,
                           "result": "OK" if (not np.isnan(dmax) and dmax < 0.1) else "DIFF",
                           "detail": "geometry file rounds elevations to 0.1 m"})
INTEG = pd.DataFrame(integ_rows)
INTEG.to_csv(os.path.join(OUT, "data_integrity_check.csv"), index=False)
bad = INTEG[~INTEG["result"].isin(["OK"])]
print(f"  {len(INTEG)} integrity checks; {len(bad)} not OK:")
for _, b in bad.iterrows():
    print(f"   - {b['well']} {b['check']}: {b['result']} ({b['detail']})")

# ============================================================================
# 6. Natural-state timing test (survey P at gauge elevation vs first transient P)
# ============================================================================
print("\n=== 3.1 Natural-state timing ===")
pc_rows = []
for w, r in W.items():
    for el, (dates, vals) in r["transient_P"].items():
        ps = float(interp_profile(r["zP"], r["P"], el)[0])
        # value of the transient series one year later, to show how fast P falls
        pc_rows.append({"well": w, "gauge_elevation_mRL": el, "survey_date": r["survey_date_P"],
                        "survey_P_at_gauge_bar": round(ps, 2),
                        "first_transient_date": dates[0], "first_transient_P_bar": vals[0],
                        "diff_survey_minus_first_bar": round(ps - vals[0], 2),
                        "second_transient_date": dates[1], "second_transient_P_bar": vals[1],
                        "diff_survey_minus_second_bar": round(ps - vals[1], 2),
                        "transient_drop_first_year_bar": round(vals[0] - vals[1], 2),
                        "transient_last_date": dates[-1], "transient_last_P_bar": vals[-1]})
PC = pd.DataFrame(pc_rows)
# Field drawdown 2014 -> 2020 at the three 2014 gauges (used to judge WELL-10)
dd_2020 = {}
for w in ["WELL-2", "WELL-3", "WELL-9"]:
    dates, vals = list(W[w]["transient_P"].values())[0]
    i20 = dates.index("01/01/2020")
    dd_2020[w] = vals[0] - vals[i20]

NS = {}  # per-well verdict
for w, r in W.items():
    role = "producer" if r["prod"] else ("injector" if r["inj"] else "observation/unused")
    first_flow = None
    flows = r["prod"] or r["inj"]
    if flows:
        first_flow = next(d for d, v in flows.items() if v > 0)
    sub = PC[PC.well == w]
    if len(sub):
        diff = float(sub.diff_survey_minus_first_bar.iloc[0])
        first_date = sub.first_transient_date.iloc[0]
        if first_date.endswith("2014"):
            if abs(diff) <= 1.0:
                v, c, why = "natural state", "high", (
                    f"survey P at gauge ({sub.survey_P_at_gauge_bar.iloc[0]:.1f} bar) is within "
                    f"{abs(diff):.1f} bar of the first (01/01/2014, pre-production) gauge value; the 2015 gauge value is "
                    f"{sub.diff_survey_minus_second_bar.iloc[0]:.1f} bar lower than the survey")
            elif abs(diff) <= 3.0:
                v, c, why = "probably natural state", "medium", (
                    f"survey P at gauge differs by {diff:+.1f} bar from the first (2014) gauge value; much closer to 2014 "
                    f"than to the 2015 value ({sub.diff_survey_minus_second_bar.iloc[0]:+.1f} bar), but the offset "
                    f"exceeds the 0.1 bar reading precision (deviated well: gauge depth reference uncertain)")
            else:
                v, c, why = "possibly disturbed", "low", f"survey differs from first gauge value by {diff:+.1f} bar"
        else:
            mean_dd = np.mean(list(dd_2020.values()))
            v, c, why = "probably natural state", "medium", (
                f"gauge series starts {first_date}, after 6 years of field production; survey is {diff:+.1f} bar above "
                f"the 2020 gauge value, which is the right sign and the same order as the 2014-2020 drawdown at the "
                f"other gauges ({min(dd_2020.values()):.1f}-{max(dd_2020.values()):.1f} bar). Cannot be confirmed directly; "
                f"near-boiling profile adds uncertainty")
    else:
        if role == "injector":
            v, c, why = "probably natural state", "medium", (
                f"no pressure gauge. Injection of {list(r['inj_h'].values())[0]:.1f} kJ/kg water (about 80 degC) since "
                f"{first_flow}; survey dated {r['survey_date_T']}. Profile shows no isothermal ~80 degC injectate column; "
                f"producer surveys with the same date match pre-production gauges, so the survey date is probably not "
                f"the measurement date")
        else:
            v, c, why = "probably natural state", "medium", (
                f"no pressure gauge and no flow history. Survey dated {r['survey_date_T']} (after production start "
                f"01/01/2014) but producer surveys with the same date match the pre-production gauges")
    NS[w] = {"role": role, "verdict": v, "confidence": c, "reason": why, "first_flow": first_flow}
CONF_W = {"high": 1.0, "medium": 0.8, "low": 0.5}
PC.to_csv(os.path.join(OUT, "pressure_consistency_check.csv"), index=False)
for w in WELLS:
    print(f"  {w:8s} {NS[w]['role']:18s} {NS[w]['verdict']:24s} ({NS[w]['confidence']}): {NS[w]['reason'][:110]}...")

# ============================================================================
# 7. Per-point table: geometry, Tsat, boiling margin, location class, weights
# ============================================================================
print("\n=== Calibration-target table ===")
WCLASS = {"feed zone": 1.0, "open hole": 0.8, "open hole, possible internal flow": 0.4,
          "behind casing": 0.5, "near wellhead / top layer": 0.3, "above model surface": 0.0}


def isothermal_internal_flow_mask(r):
    """Open-hole points between the shallowest and deepest feed zone where the profile is
    near-isothermal (|dT/dz| < 0.03 degC/m over >= 200 m): candidate wellbore internal flow."""
    z, T = r["zT"], r["T"]
    m = np.zeros(len(z), bool)
    if len(r["fz"]) < 2:
        return m, []
    ztop = min(f["bot"] for f in r["fz"][:1])
    zbot = max(f["top"] for f in r["fz"][-1:])
    idx = [i for i in range(len(z)) if zbot <= z[i] <= ztop]
    segs = []
    # sliding search for sections with small gradient
    i = 0
    while i < len(idx) - 1:
        j = i + 1
        while j < len(idx) and abs(T[idx[j]] - T[idx[i]]) / max(1e-9, abs(z[idx[i]] - z[idx[j]])) < 0.03 \
                and abs(T[idx[j]] - T[idx[j - 1]]) / max(1e-9, abs(z[idx[j - 1]] - z[idx[j]])) < 0.03:
            j += 1
        if j - 1 > i and z[idx[i]] - z[idx[j - 1]] >= 200:
            for k in idx[i:j]:
                m[k] = True
            segs.append((z[idx[i]], z[idx[j - 1]], T[idx[i]], T[idx[j - 1]]))
            i = j
        else:
            i += 1
    return m, segs


rows = []
INTFLOW = {}
for w, r in W.items():
    mflow, segs = isothermal_internal_flow_mask(r)
    INTFLOW[w] = segs
    for i, (z, T, P) in enumerate(zip(r["zT"], r["T"], r["P"])):
        x, y = track_xy(r, z)
        c, lay, blk, above = block_of(x, y, z)
        ts = tsat_c(P)
        bm = ts - T if not np.isnan(ts) else np.nan
        in_fz = any(f["bot"] <= z <= f["top"] for f in r["fz"])
        behind = z >= r["casing_shoe"]
        toplay = top_layer_of_column(c) if c else None
        in_top_layer = (lay is not None and toplay is not None and lay == toplay["name"].strip())
        near_wh = (r["wellhead"] - z) <= 30.0
        above_wl = P < 1.2  # air/vapour-filled part of the wellbore at survey time
        flags = []
        if above:
            cls = "above model surface"
        elif in_fz:
            cls = "feed zone"
        elif near_wh or in_top_layer:
            cls = "near wellhead / top layer"
        elif behind:
            cls = "behind casing"
        elif mflow[i]:
            cls = "open hole, possible internal flow"
        else:
            cls = "open hole"
        if above_wl:
            flags.append("above wellbore water level (P<1.2 bar)")
        if not np.isnan(bm) and bm < BM_FLAG:
            flags.append(f"boiling margin < {BM_FLAG:g} degC")
        if np.isnan(ts):
            flags.append("P above critical pressure: no Tsat")
        wcls = WCLASS[cls]
        conf = NS[w]["confidence"]
        weight = round(wcls * CONF_W[conf], 3)
        reason = f"{cls} ({wcls}) x NS confidence {conf} ({CONF_W[conf]})"
        # pressure weight: formation pressure only at feed zones
        wP = round(1.0 * CONF_W[conf], 3) if (in_fz and not above) else 0.0
        rows.append({"well": w, "x": round(x, 1), "y": round(y, 1), "elevation_mRL": z,
                     "T_degC": T, "P_bar": P, "Tsat_degC": round(ts, 2) if not np.isnan(ts) else np.nan,
                     "boiling_margin_degC": round(bm, 2) if not np.isnan(bm) else np.nan,
                     "column": c, "layer": lay, "block": blk, "in_feedzone": in_fz,
                     "behind_casing": behind, "above_model_surface": bool(above),
                     "ns_confidence": conf, "weight": weight, "weight_reason": reason,
                     "location_class": cls, "weight_P": wP,
                     "rho_IF97_kgm3": round(rho_liquid(T, P), 1) if (np.isnan(bm) or bm > 0) else np.nan,
                     "flags": "; ".join(flags)})
TGT = pd.DataFrame(rows)
TGT.to_csv(os.path.join(OUT, "NS_calibration_targets.csv"), index=False)
print(f"  {len(TGT)} survey points written; class counts:")
print(TGT.location_class.value_counts().to_string())
print("  Points with boiling margin < 5 degC:")
print(TGT[TGT.boiling_margin_degC < BM_FLAG][["well", "elevation_mRL", "T_degC", "P_bar",
                                               "Tsat_degC", "boiling_margin_degC",
                                               "location_class"]].to_string(index=False))
for w, s in INTFLOW.items():
    if s:
        print(f"  {w}: near-isothermal open-hole section(s) between feed zones (possible internal flow): "
              + ", ".join(f"{a:g} to {b:g} mRL ({ta:.1f}-{tb:.1f} degC)" for a, b, ta, tb in s))

# ============================================================================
# 8. Well geometry check
# ============================================================================
print("\n=== Well geometry ===")
geo_rows = []
for w, r in W.items():
    x0, y0, z0 = r["track"][0]
    xt, yt, zt = r["track"][-1]
    c = find_column(x0, y0)
    s = col_surface(c)
    tl = top_layer_of_column(c)
    gt = r["raw"]
    row = {"well": w, "x_wellhead": x0, "y_wellhead": y0, "wellhead_mRL": z0, "column": c,
           "model_surface_mRL": s, "wellhead_minus_surface_m": z0 - s,
           "top_layer": tl["name"].strip() if tl else None,
           "top_block_centre_mRL": (s + tl["bot"]) / 2 if tl else None,
           "TD_mRL": zt, "horizontal_offset_at_TD_m": round(math.hypot(xt - x0, yt - y0), 1),
           "deviated": math.hypot(xt - x0, yt - y0) > 1,
           "casing_shoe_mRL": r["casing_shoe"], "liner_top_mRL": r["liner_top"],
           "liner_bottom_mRL": r["liner_bot"], "first_survey_mRL": r["zT"].max(),
           "survey_top_below_wellhead_m": z0 - r["zT"].max()}
    for fname, pl in FAULTS.items():
        row[f"dist_{fname.replace(' ', '')}_m"] = round(dist_point_polyline(x0, y0, pl), 0)
    for k, v in FAULT_X.items():
        row[f"dist_{k.replace(' ', '')}_m"] = round(math.hypot(x0 - v[0], y0 - v[1]), 0)
    row["dist_river_m"] = round(min(dist_point_polyline(x0, y0, sg) for sg in RIVERS["river"]), 0)
    row["dist_TaraweraRiver_m"] = round(min(dist_point_polyline(x0, y0, sg) for sg in RIVERS["Tarawera River"]), 0)
    row["dist_lake_shore_m"] = round(dist_point_polyline(x0, y0, LAKE), 0)
    row["inside_lake_polygon"] = point_in_poly(x0, y0, LAKE)
    # FZ blocks vs gTE14515_wells.json
    geo_rows.append(row)
GEO = pd.DataFrame(geo_rows)
GEO.to_csv(os.path.join(OUT, "well_geometry_check.csv"), index=False)
print(GEO[["well", "wellhead_mRL", "model_surface_mRL", "wellhead_minus_surface_m",
           "horizontal_offset_at_TD_m", "dist_Fault1_m", "dist_Fault2_m", "dist_Fault3_m",
           "dist_river_m", "dist_lake_shore_m"]].to_string(index=False))

# ============================================================================
# 9. Feed zones
# ============================================================================
print("\n=== Feed zones ===")
MWELLS = json.load(open(os.path.join(MODEL, "gTE14515_wells.json")))
fz_rows = []
for w, r in W.items():
    mfz = MWELLS.get(w, {}).get("Feed_Zones", {})
    qrank = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}
    best = max(range(len(r["fz"])), key=lambda k: (qrank.get(r["fz"][k]["quality"].upper(), 0), -k))
    for k, f in enumerate(r["fz"]):
        zm = 0.5 * (f["top"] + f["bot"])
        x, y = track_xy(r, zm)
        T = float(interp_profile(r["zT"], r["T"], zm)[0])
        P = float(interp_profile(r["zP"], r["P"], zm)[0])
        ts = tsat_c(P)
        c, lay, blk, above = block_of(x, y, zm)
        blocks_span = sorted({block_of(*track_xy(r, zz), zz)[2] for zz in np.arange(f["top"], f["bot"] - 0.1, -5.0)}
                             | {block_of(*track_xy(r, f["bot"]), f["bot"])[2]})
        mblk = list(mfz.get(f"Feed {k + 1}", {}).keys())
        npts = int(((r["zT"] <= f["top"]) & (r["zT"] >= f["bot"])).sum())
        fz_rows.append({"well": w, "feedzone": k + 1, "top_mRL": f["top"], "bottom_mRL": f["bot"],
                        "mid_mRL": zm, "quality": f["quality"], "dominant_feedzone": k == best,
                        "x": round(x, 1), "y": round(y, 1), "T_degC": round(T, 1), "P_bar": round(P, 2),
                        "Tsat_degC": round(ts, 1) if not np.isnan(ts) else np.nan,
                        "boiling_margin_degC": round(ts - T, 1) if not np.isnan(ts) else np.nan,
                        "h_liquid_kJkg": round(h_liquid(T, P), 1),
                        "column": c, "layer": lay, "block": blk,
                        "blocks_spanned": " ".join(b for b in blocks_span if b),
                        "model_wells_json_block": " ".join(mblk),
                        "block_matches_wells_json": (blk in mblk) if mblk else None,
                        "survey_points_in_fz": npts,
                        "in_slotted_liner": r["liner_top"] >= f["top"] and r["liner_bot"] <= f["bot"],
                        "ns_confidence": NS[w]["confidence"],
                        "quality_weight": {"HIGH": 1.0, "MEDIUM": 0.8, "LOW": 0.5}.get(f["quality"].upper(), 0.5),
                        "target_weight": round({"HIGH": 1.0, "MEDIUM": 0.8, "LOW": 0.5}.get(f["quality"].upper(), 0.5)
                                               * CONF_W[NS[w]["confidence"]], 3)})
FZ = pd.DataFrame(fz_rows)
FZ.to_csv(os.path.join(OUT, "feedzones.csv"), index=False)
print(FZ[["well", "feedzone", "mid_mRL", "quality", "T_degC", "P_bar", "boiling_margin_degC", "block",
          "model_wells_json_block", "target_weight"]].to_string(index=False))

# ============================================================================
# 10. Hydrostatic checks (single-phase liquid, long intervals, IF97 density)
# ============================================================================
print("\n=== 3.3 Hydrostatic checks ===")


def column_dp(r, z1, z2, dz=1.0):
    """Predicted pressure increase (bar) from z1 (top) to z2 (bottom) for a static liquid
    column with the measured temperature profile, integrated with IF97 Region 1 density."""
    zs = np.arange(z1, z2 - 1e-9, -dz)
    if zs[-1] > z2:
        zs = np.append(zs, z2)
    Tz = interp_profile(r["zT"], r["T"], zs)
    p = float(interp_profile(r["zP"], r["P"], z1)[0])
    rhos = []
    for i in range(len(zs)):
        rho = rho_liquid(Tz[i], p)
        rhos.append(rho)
        if i < len(zs) - 1:
            h = zs[i] - zs[i + 1]
            rho2 = rho_liquid(Tz[i + 1], p + rho * G * h / 1e5)
            p += 0.5 * (rho + rho2) * G * h / 1e5
    return p - float(interp_profile(r["zP"], r["P"], z1)[0]), float(np.mean(rhos))


hy_rows = []
for w, r in W.items():
    sub = TGT[TGT.well == w].sort_values("elevation_mRL", ascending=False)
    ok = sub[(sub.P_bar >= 1.5) & ((sub.boiling_margin_degC >= 10) | sub.boiling_margin_degC.isna())]
    # contiguous single-phase run(s)
    zs = ok.elevation_mRL.values
    allz = sub.elevation_mRL.values
    okset = set(zs)
    runs, cur = [], []
    for z in allz:
        if z in okset:
            cur.append(z)
        else:
            if cur:
                runs.append(cur)
            cur = []
    if cur:
        runs.append(cur)
    windows = []
    for run in runs:
        if run[0] - run[-1] >= 300:
            windows.append(("full single-phase run", run[0], run[-1]))
        i = 0
        while i < len(run) - 1:
            j = i + 1
            while j < len(run) and run[i] - run[j] < 300:
                j += 1
            if j < len(run):
                windows.append(("window >=300 m", run[i], run[j]))
            i = j
    for kind, za, zb in windows:
        pa = float(interp_profile(r["zP"], r["P"], za)[0])
        pb = float(interp_profile(r["zP"], r["P"], zb)[0])
        dp_pred, rho_mean = column_dp(r, za, zb)
        dp_meas = pb - pa
        rho_impl = dp_meas * 1e5 / (G * (za - zb))
        Tmin = float(sub[(sub.elevation_mRL <= za) & (sub.elevation_mRL >= zb)].T_degC.min())
        Tmax = float(sub[(sub.elevation_mRL <= za) & (sub.elevation_mRL >= zb)].T_degC.max())
        hy_rows.append({"well": w, "window_type": kind, "z_top_mRL": za, "z_bottom_mRL": zb,
                        "length_m": za - zb, "T_min_degC": Tmin, "T_max_degC": Tmax,
                        "P_top_bar": pa, "P_bottom_bar": pb, "dP_measured_bar": round(dp_meas, 2),
                        "dP_IF97_static_bar": round(dp_pred, 2),
                        "dP_excess_bar": round(dp_meas - dp_pred, 2),
                        "rho_implied_kgm3": round(rho_impl, 1), "rho_IF97_mean_kgm3": round(rho_mean, 1),
                        "ratio_measured_to_IF97": round(dp_meas / dp_pred, 4),
                        "implied_g_ms2": round(G * dp_meas / dp_pred, 3),
                        "steeper_than_cold_water": bool(rho_impl > 1000.0),
                        "rounding_uncertainty_ratio": round(0.1 / dp_pred, 4)})
HY = pd.DataFrame(hy_rows)
full = HY[HY.window_type == "full single-phase run"].copy()
cold = ["WELL-1", "WELL-4", "WELL-8"]  # Tmax < 120 degC: no hot column, no near-boiling zone
FIELD_RATIO = float(full[full.well.isin(cold)].ratio_measured_to_IF97.median())
FIELD_RATIO_ALL = float(full.ratio_measured_to_IF97.median())
# Depth bands: the excess is not uniform with depth, so the field-wide component is
# estimated per band from the cold wells, and each window is compared with its own band.
SHALLOW_CUT = -450.0  # windows whose bottom is at or above this elevation are "shallow"
HY["band"] = np.where(HY.z_bottom_mRL >= SHALLOW_CUT, "shallow (bottom >= -450 mRL)",
                      np.where(HY.z_top_mRL <= SHALLOW_CUT, "deep (top <= -450 mRL)", "spanning"))
HY.loc[HY.window_type.str.startswith("full"), "band"] = "full run"
BAND_RATIO = {}
for b in HY.band.unique():
    m = (HY.band == b) & HY.well.isin(cold)
    BAND_RATIO[b] = float(HY[m].ratio_measured_to_IF97.median()) if m.any() else np.nan
HY["field_ratio_cold_wells_same_band"] = HY.band.map(BAND_RATIO).round(4)
HY["excess_vs_field_ratio_pct"] = ((HY.ratio_measured_to_IF97 / HY.field_ratio_cold_wells_same_band - 1)
                                   * 100).round(2)
HY["excess_vs_IF97_pct"] = ((HY.ratio_measured_to_IF97 - 1) * 100).round(2)
print("  Cold-well median ratio by depth band:", {k: round(v, 4) for k, v in BAND_RATIO.items()})
BANDSUM = HY[HY.window_type != "full single-phase run"].groupby(["well", "band"]).agg(
    n=("ratio_measured_to_IF97", "size"), ratio_median=("ratio_measured_to_IF97", "median"),
    ratio_min=("ratio_measured_to_IF97", "min"), ratio_max=("ratio_measured_to_IF97", "max")).round(4).reset_index()
print(BANDSUM.to_string(index=False))
BANDSUM.to_csv(os.path.join(OUT, "hydrostatic_band_summary.csv"), index=False)
HY.to_csv(os.path.join(OUT, "hydrostatic_check.csv"), index=False)
print(f"  Field-wide ratio measured/IF97 (median of full runs in cold wells {cold}): {FIELD_RATIO:.4f}"
      f" -> implied g = {G * FIELD_RATIO:.3f} m/s2; median all wells {FIELD_RATIO_ALL:.4f}")
print(full[["well", "z_top_mRL", "z_bottom_mRL", "dP_measured_bar", "dP_IF97_static_bar", "ratio_measured_to_IF97",
            "implied_g_ms2"]].to_string(index=False))
print("  Per-window excess relative to field-wide ratio (|x| > 1%):")
big = HY[HY.excess_vs_field_ratio_pct.abs() > 1.0]
print(big[["well", "window_type", "z_top_mRL", "z_bottom_mRL", "ratio_measured_to_IF97",
           "excess_vs_field_ratio_pct"]].to_string(index=False) if len(big) else "   none")

# ============================================================================
# 11. Per-well summary
# ============================================================================
print("\n=== Per-well summary ===")


def first_crossing(z, T, iso):
    for i in range(len(z) - 1):
        if (T[i] - iso) * (T[i + 1] - iso) <= 0 and T[i] != T[i + 1] and T[i] < iso <= T[i + 1]:
            return z[i] + (iso - T[i]) * (z[i + 1] - z[i]) / (T[i + 1] - T[i])
    return np.nan


def linfit_grad(z, T):
    if len(z) < 2:
        return np.nan
    return -np.polyfit(z, T, 1)[0] * 1000.0  # degC/km, positive = increasing with depth


sum_rows = []
for w, r in W.items():
    z, T, P = r["zT"], r["T"], r["P"]
    sub = TGT[TGT.well == w]
    i_max = int(np.argmax(T))
    Tmax, zTmax = T[i_max], z[i_max]
    below = slice(i_max, len(z))
    i_min_b = i_max + int(np.argmin(T[below]))
    inv = Tmax - T[i_min_b]
    has_inv = inv > 2.0 and i_min_b != i_max
    reheat = T[-1] - T[i_min_b] if has_inv else np.nan
    # shallow gradient: liquid-filled top 250 m of the survey (behind casing)
    ztop_l = sub[sub.P_bar >= 1.2].elevation_mRL.max()
    msh = (z <= ztop_l) & (z >= ztop_l - 250)
    g_sh = linfit_grad(z[msh], T[msh])
    # deep gradient: below deepest feed zone (static part of the well) if >= 200 m, else deepest 400 m
    # deep gradient: static part of the well below the deepest feed zone and below any
    # inversion minimum; needs >= 3 points over >= 200 m, otherwise not reported
    zfz = min(f["bot"] for f in r["fz"])
    cands = [(min(zfz, z[i_min_b]) if has_inv else zfz, "below deepest FZ" + (" and inversion base" if has_inv else ""))]
    if has_inv:
        cands.append((z[i_min_b], "below inversion base (includes open hole between/at FZ)"))
    else:
        cands.append((z.min() + 500, "deepest 500 m"))
    g_dp, dp_int = np.nan, "no static interval >= 200 m with >= 3 points below the inversion"
    for zstart, lab in cands:
        mdp = z <= zstart
        if mdp.sum() >= 3 and (z[mdp].max() - z[mdp].min()) >= 200:
            g_dp = linfit_grad(z[mdp], T[mdp])
            dp_int = f"{z[mdp].max():g} to {z[mdp].min():g} ({lab})"
            break
    if sub.boiling_margin_degC.min() < BM_FLAG:
        g_dp, dp_int = np.nan, ("not reported: near-boiling convective well, TD only "
                                f"{r['wellhead'] - r['td']:.0f} m below wellhead; isothermal below feed zone")
    z100 = first_crossing(z, T, 100.0)
    z200 = first_crossing(z, T, 200.0)
    bm = sub.boiling_margin_degC.values
    ib = int(np.nanargmin(bm))
    # water level: from first two liquid points (P >= 1.5 bar) extrapolate up to 1.01325 bar
    liq = sub[sub.P_bar >= 1.5].sort_values("elevation_mRL", ascending=False)
    z1, p1, rho1 = liq.elevation_mRL.iloc[0], liq.P_bar.iloc[0], liq.rho_IF97_kgm3.iloc[0]
    wl = z1 + (p1 - 1.01325) * 1e5 / (rho1 * G)
    if sub.P_bar.iloc[0] >= 1.2:
        wl_note = ("pressurised at first survey point "
                   f"({sub.P_bar.iloc[0]:.1f} bar at {sub.elevation_mRL.iloc[0]:g} mRL)")
        if wl > r["wellhead"]:
            wl_note += (f"; equivalent level {wl - r['wellhead']:.0f} m above wellhead "
                        "(artesian or shut-in wellhead pressure)")
    else:
        wl_note = "air/vapour above level (P approx. 1 bar at top)"
    fzr = FZ[(FZ.well == w) & FZ.dominant_feedzone].iloc[0]
    # regime classification (rules, see report)
    min_bm_open = sub[~sub.behind_casing].boiling_margin_degC.min()
    if min_bm_open < BM_FLAG or bm[ib] < BM_FLAG:
        regime = "upflow / two-phase (near-boiling)"
    elif has_inv and inv >= 20:
        regime = ("outflow (hot tongue over cooler water)" if Tmax >= 200 else
                  "distal outflow (cooler tongue over cooler water)") + ("; deep reheating" if reheat > 10 else "")
    elif has_inv and inv >= 5:
        regime = "distal outflow (weak inversion)"
    elif Tmax < 100 or (g_sh < 40 and Tmax < 130):
        regime = "recharge / cold downflow"
    else:
        regime = "conductive"
    q = K_ROCK * g_dp / 1000.0 if not np.isnan(g_dp) else np.nan
    q_sh = K_ROCK * g_sh / 1000.0 if not np.isnan(g_sh) else np.nan
    sum_rows.append({"well": w, "role": NS[w]["role"], "wellhead_mRL": r["wellhead"], "TD_mRL": r["td"],
                     "Tmax_degC": Tmax, "Tmax_elev_mRL": zTmax, "T_at_TD_degC": T[-1],
                     "inversion_degC": round(inv, 1) if has_inv else 0.0,
                     "inversion_top_mRL": zTmax if has_inv else np.nan,
                     "inversion_base_mRL": z[i_min_b] if has_inv else np.nan,
                     "T_inversion_min_degC": T[i_min_b] if has_inv else np.nan,
                     "reheat_below_inversion_degC": round(reheat, 1) if has_inv else np.nan,
                     "shallow_gradient_degC_per_km": round(g_sh, 0),
                     "shallow_gradient_interval_mRL": f"{z[msh].max():g} to {z[msh].min():g}",
                     "deep_gradient_degC_per_km": round(g_dp, 0), "deep_gradient_interval_mRL": dp_int,
                     "conductive_flux_deep_Wm2_assumed_k2": round(q, 3),
                     "conductive_flux_cap_Wm2_assumed_k2": round(q_sh, 3),
                     "z_100C_mRL": round(z100, 0) if not np.isnan(z100) else np.nan,
                     "z_200C_mRL": round(z200, 0) if not np.isnan(z200) else np.nan,
                     "min_boiling_margin_degC": round(bm[ib], 1), "min_bm_elev_mRL": z[ib],
                     "n_points_bm_lt5": int((bm < BM_FLAG).sum()),
                     "water_level_mRL": round(wl, 0), "water_level_note": wl_note,
                     "pivot_fz_mid_mRL": fzr.mid_mRL, "pivot_P_bar": fzr.P_bar, "pivot_fz_quality": fzr.quality,
                     "regime": regime, "ns_verdict": NS[w]["verdict"], "ns_confidence": NS[w]["confidence"],
                     "casing_shoe_mRL": r["casing_shoe"],
                     "internal_flow_sections": "; ".join(f"{a:g} to {b:g} mRL" for a, b, _, _ in INTFLOW[w])})
SUM = pd.DataFrame(sum_rows)
SUM.to_csv(os.path.join(OUT, "well_summary.csv"), index=False)
print(SUM[["well", "Tmax_degC", "Tmax_elev_mRL", "T_at_TD_degC", "inversion_degC", "shallow_gradient_degC_per_km",
           "deep_gradient_degC_per_km", "z_100C_mRL", "z_200C_mRL", "min_boiling_margin_degC", "water_level_mRL",
           "regime"]].to_string(index=False))

# ============================================================================
# 12. T and P at layer centres and at common elevations (interpolation only)
# ============================================================================
lc_rows = []
for w, r in W.items():
    for L in LAYERS:
        zc = L["centre"]
        if not (r["zT"].min() <= zc <= r["zT"].max()):
            continue
        x, y = track_xy(r, zc)
        c, lay, blk, above = block_of(x, y, zc)
        tl = top_layer_of_column(c)
        is_top = tl is not None and tl["name"] == L["name"]
        if above:
            continue
        T = float(interp_profile(r["zT"], r["T"], zc)[0])
        P = float(interp_profile(r["zP"], r["P"], zc)[0])
        ts = tsat_c(P)
        sub = TGT[TGT.well == w]
        lc_rows.append({"well": w, "layer": L["name"].strip(), "layer_centre_mRL": zc,
                        "layer_top_mRL": L["top"], "layer_bottom_mRL": L["bot"], "column": c, "block": blk,
                        "x": round(x, 1), "y": round(y, 1), "T_degC": round(T, 1), "P_bar": round(P, 2),
                        "Tsat_degC": round(ts, 1) if not np.isnan(ts) else np.nan,
                        "boiling_margin_degC": round(ts - T, 1) if not np.isnan(ts) else np.nan,
                        "top_layer_of_column": is_top,
                        "behind_casing": zc >= r["casing_shoe"],
                        "in_feedzone": any(f["bot"] <= zc <= f["top"] for f in r["fz"]),
                        "note": "nominal layer centre lies in truncated top block" if is_top else ""})
LC = pd.DataFrame(lc_rows)
LC.to_csv(os.path.join(OUT, "T_P_at_layer_centres.csv"), index=False)
print(f"\nT/P at layer centres: {len(LC)} rows")

COMMON = [100, 0, -250, -500, -750, -1000, -1500, -2000, -2500]
ce_rows = []
for zc in COMMON:
    for w, r in W.items():
        T = float(interp_profile(r["zT"], r["T"], zc)[0])
        if np.isnan(T):
            continue
        P = float(interp_profile(r["zP"], r["P"], zc)[0])
        x, y = track_xy(r, zc)
        ts = tsat_c(P)
        ce_rows.append({"elevation_mRL": zc, "well": w, "x": round(x, 1), "y": round(y, 1), "T_degC": round(T, 1),
                        "P_bar": round(P, 2), "Tsat_degC": round(ts, 1) if not np.isnan(ts) else np.nan,
                        "boiling_margin_degC": round(ts - T, 1) if not np.isnan(ts) else np.nan,
                        "behind_casing": zc >= r["casing_shoe"],
                        "P_minus_field_median_bar": np.nan})
CE = pd.DataFrame(ce_rows)
for zc in COMMON:
    m = CE.elevation_mRL == zc
    CE.loc[m, "P_minus_field_median_bar"] = (CE.loc[m, "P_bar"] - CE.loc[m, "P_bar"].median()).round(2)
CE.to_csv(os.path.join(OUT, "T_P_at_common_elevations.csv"), index=False)

# Feed-zone pressures reduced to a common datum, to test deep pressure connection.
# Each FZ pressure is carried to -1000 mRL along a hypothetical static column with the
# field-median temperature-vs-elevation profile? No: we use a single reference liquid density
# so the comparison does not mix wellbore fluids. Report as "equivalent head".
fz_head = []
for _, f in FZ.iterrows():
    rho_fz = rho_liquid(f.T_degC, f.P_bar)
    head = f.mid_mRL + (f.P_bar - 1.01325) * 1e5 / (rho_fz * G)
    head_field = f.mid_mRL + (f.P_bar / FIELD_RATIO - 1.01325) * 1e5 / (rho_fz * G)
    fz_head.append({"well": f.well, "feedzone": f.feedzone, "mid_mRL": f.mid_mRL, "T_degC": f.T_degC,
                    "P_bar": f.P_bar, "rho_fz_IF97": round(rho_fz, 1),
                    "equiv_head_mRL_local_density": round(head, 1)})
FZH = pd.DataFrame(fz_head)
FZ = FZ.merge(FZH[["well", "feedzone", "rho_fz_IF97", "equiv_head_mRL_local_density"]], on=["well", "feedzone"])
FZ.to_csv(os.path.join(OUT, "feedzones.csv"), index=False)

# Fit a single field pressure-elevation curve through all feed zones; residuals show connection
fzfit = np.polyfit(FZ.mid_mRL, FZ.P_bar, 2)
FZ["P_residual_to_field_fit_bar"] = (FZ.P_bar - np.polyval(fzfit, FZ.mid_mRL)).round(2)
FZ.to_csv(os.path.join(OUT, "feedzones.csv"), index=False)
print("Feed-zone pressure residuals to a single field quadratic P(z) fit (bar):")
print(FZ[["well", "feedzone", "mid_mRL", "P_bar", "P_residual_to_field_fit_bar"]].to_string(index=False))

# ============================================================================
# 13. Reference model description (labelled; NOT evidence)
# ============================================================================
REF = {}
try:
    bc = json.load(open(os.path.join(MODEL, "gTE14515_boundary_conditions.json")))
    REF["heat_flux_default_Wm2"] = bc["Heat"]["Default"]["Heat_Flux"]
    REF["rain_default_mm"] = bc["Rain"]["Default"]["Rainfall"]
    REF["infiltration"] = bc["Rain"]["Default"]["Infiltration_Fraction"]
    REF["lake_columns"] = list(bc["Wet_Columns"]["Zones"]["Lake"]["Columns"].keys())
    gen_heat, gen_mass, deep_mass = [], 0, 0
    with open(os.path.join(MODEL, "NaturalState", "TE14515_NS_000.dat")) as f:
        insec = False
        for ln in f:
            if ln.startswith("GENER"):
                insec = True
                continue
            if insec:
                if not ln.strip() or ln[:5].strip() == "" and len(ln.strip()) == 0:
                    break
                if "HEAT" in ln:
                    gen_heat.append(float(ln.split()[-4]))
                elif "MASS" in ln or "COM1" in ln:
                    gen_mass += 1
                    if ln[3:5] in ("69", "68", "67"):
                        deep_mass += 1
    REF["n_heat_gen"] = len(gen_heat)
    REF["total_heat_MW"] = round(sum(gen_heat) / 1e6, 2)
    REF["n_mass_gen"] = gen_mass
    REF["n_deep_mass_gen"] = deep_mass
    print(f"\nReference model (NOT evidence): {REF['n_heat_gen']} basal HEAT generators, total "
          f"{REF['total_heat_MW']} MW; {REF['n_mass_gen']} MASS generators (rain, near surface), "
          f"{REF['n_deep_mass_gen']} in the deepest three layers")
except Exception as e:  # reference model description is optional
    print("Reference model description skipped:", e)

# ============================================================================
# 14. Figures
# ============================================================================
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3})
CMAP_T = "inferno"
xs_all = np.array([W[w]["track"][0, 0] for w in WELLS])
ys_all = np.array([W[w]["track"][0, 1] for w in WELLS])
XMIN, XMAX = xs_all.min() - 2500, xs_all.max() + 2500
YMIN, YMAX = ys_all.min() - 2500, ys_all.max() + 2500


def draw_base(ax, grid=True, zoom=True):
    if grid:
        for c, p in COLS.items():
            ax.add_patch(MplPolygon(p, closed=True, fill=False, lw=0.3, ec="0.75"))
    ax.add_patch(MplPolygon(LAKE, closed=True, fc="#9ecae1", ec="#3182bd", lw=0.8, alpha=0.6, label="Lake"))
    for k, pl in RIVERS.items():
        for seg in pl:
            ax.plot(seg[:, 0], seg[:, 1], color="#3182bd", lw=1.2)
    cx, cy = (XMIN + XMAX) / 2, (YMIN + YMAX) / 2
    for k, pl in FAULTS.items():
        ax.plot(pl[:, 0], pl[:, 1], "k-", lw=1.6)
        # label at the point of the trace nearest the view centre, shifted along the trace
        d = pl[1] - pl[0]
        t = np.clip(((cx - pl[0][0]) * d[0] + (cy - pl[0][1]) * d[1]) / (d @ d), 0, 1)
        t = np.clip(t + {"Fault 1": -0.25, "Fault 2": -0.3, "Fault 3": 0.12}.get(k, 0), 0, 1)
        q = pl[0] + t * d
        ax.annotate(k, q, fontsize=7, color="k", xytext=(3, 0), textcoords="offset points",
                    bbox=dict(fc="w", ec="none", alpha=0.7, pad=0.5))
    for k, v in FAULT_X.items():
        ax.plot(*v, "kx", ms=7)
    if zoom:
        ax.set_xlim(XMIN, XMAX)
        ax.set_ylim(YMIN, YMAX)
    ax.set_aspect("equal")
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.ticklabel_format(useOffset=False, style="plain")


def draw_wells(ax, label=True, color="k"):
    for w, r in W.items():
        tr = r["track"]
        ax.plot(tr[:, 0], tr[:, 1], "-", color=color, lw=1)
        ax.plot(tr[0, 0], tr[0, 1], "o", mfc="w", mec=color, ms=5)
        if label:
            ax.annotate(w.replace("WELL-", "W"), (tr[0, 0], tr[0, 1]), xytext=(4, 4),
                        textcoords="offset points", fontsize=8, weight="bold")


figs = []
# --- well map
fig, axs = plt.subplots(1, 2, figsize=(15, 7.5))
ax = axs[0]
draw_base(ax, zoom=False)
draw_wells(ax)
allp = np.vstack(list(COLS.values()))
ax.set_xlim(allp[:, 0].min(), allp[:, 0].max())
ax.set_ylim(allp[:, 1].min(), allp[:, 1].max())
ax.set_title("Model grid (gTE14515), faults, water bodies and wells")
ax = axs[1]
draw_base(ax)
sc = ax.scatter(SUM.merge(GEO, on="well").x_wellhead, SUM.merge(GEO, on="well").y_wellhead,
                c=SUM.Tmax_degC, cmap=CMAP_T, s=150, vmin=0, vmax=280, zorder=5, ec="k")
draw_wells(ax)
plt.colorbar(sc, ax=ax, label="Measured Tmax (degC)", shrink=0.8)
ax.set_title("Wells coloured by measured maximum temperature (x = fault intersection)")
fig.tight_layout()
fn = os.path.join(FIG, "well_map_grid_faults_water.png")
fig.savefig(fn, dpi=130)
plt.close(fig)
figs.append(fn)

# --- per-well profiles
for w, r in W.items():
    sub = TGT[TGT.well == w]
    fig, axs = plt.subplots(1, 3, figsize=(13, 7), sharey=True)
    z, T, P = r["zT"], r["T"], r["P"]
    ts = np.array([tsat_c(p) for p in P])
    ax = axs[0]
    ax.plot(T, z, "o-", color="tab:red", ms=3, label="Measured T")
    ax.plot(ts, z, "--", color="tab:blue", label="Tsat(P measured) IF97")
    ax.set_xlabel("Temperature (degC)")
    ax.set_ylabel("Elevation (mRL)")
    ax2 = axs[1]
    ax2.plot(P, z, "o-", color="tab:green", ms=3, label="Measured P")
    for _, f in FZ[FZ.well == w].iterrows():
        ax2.plot(f.P_bar, f.mid_mRL, "s", ms=8, mfc="none", mec="k")
    ax2.set_xlabel("Pressure (bar)")
    ax3 = axs[2]
    ax3.plot(sub.boiling_margin_degC, sub.elevation_mRL, "o-", color="tab:purple", ms=3)
    ax3.axvline(BM_FLAG, color="r", ls=":", label=f"{BM_FLAG:g} degC")
    ax3.set_xlabel("Boiling margin Tsat(P) - T (degC)")
    for a in axs:
        for f in r["fz"]:
            a.axhspan(f["bot"], f["top"], color="orange", alpha=0.35)
        a.axhline(r["casing_shoe"], color="k", lw=1.5, ls="-")
        a.axhline(r["wellhead"], color="0.4", lw=0.8, ls="--")
        for L in LAYERS:
            if r["td"] - 50 <= L["bot"] <= r["wellhead"]:
                a.axhline(L["bot"], color="0.85", lw=0.5, zorder=0)
        for s in INTFLOW[w]:
            a.axhspan(s[1], s[0], color="cyan", alpha=0.15)
    axs[0].plot([], [], color="orange", lw=8, alpha=0.35, label="Feed zone")
    axs[0].plot([], [], color="k", lw=1.5, label=f"Casing shoe {r['casing_shoe']:g} mRL")
    axs[0].plot([], [], color="0.4", lw=0.8, ls="--", label="Wellhead")
    if INTFLOW[w]:
        axs[0].plot([], [], color="cyan", lw=8, alpha=0.3, label="Isothermal open hole (possible internal flow)")
    axs[0].legend(loc="lower left", fontsize=7)
    ax3.legend(loc="lower right", fontsize=7)
    s = SUM[SUM.well == w].iloc[0]
    fig.suptitle(f"{w} ({NS[w]['role']}): Tmax {s.Tmax_degC:.1f} degC at {s.Tmax_elev_mRL:g} mRL; "
                 f"regime: {s.regime}; NS verdict: {NS[w]['verdict']} ({NS[w]['confidence']})\n"
                 f"grey lines = model layer boundaries", fontsize=9)
    fig.tight_layout()
    fn = os.path.join(FIG, f"profile_{w}.png")
    fig.savefig(fn, dpi=110)
    plt.close(fig)
    figs.append(fn)

# --- all-well overview profiles
fig, axs = plt.subplots(1, 2, figsize=(13, 8), sharey=True)
cm = plt.get_cmap("tab10")
for i, (w, r) in enumerate(W.items()):
    axs[0].plot(r["T"], r["zT"], "-", color=cm(i), label=w)
    axs[1].plot(r["P"], r["zP"], "-", color=cm(i), label=w)
zz = np.linspace(-2900, 300, 200)
axs[0].set_xlabel("Temperature (degC)")
axs[0].set_ylabel("Elevation (mRL)")
axs[1].set_xlabel("Pressure (bar)")
axs[0].legend(fontsize=8)
axs[0].set_title("All temperature profiles")
axs[1].set_title("All pressure profiles")
fig.tight_layout()
fn = os.path.join(FIG, "profiles_all_wells.png")
fig.savefig(fn, dpi=110)
plt.close(fig)
figs.append(fn)

# --- P-T state plot
fig, ax = plt.subplots(figsize=(9, 7))
Ts = np.linspace(20, 373, 300)
ax.plot(Ts, [psat_mpa(t + 273.15) * 10 for t in Ts], "k-", lw=2, label="Saturation curve (IF97 Region 4)")
for i, w in enumerate(WELLS):
    sub = TGT[TGT.well == w]
    ax.plot(sub.T_degC, sub.P_bar, "o", ms=4, color=cm(i), label=w)
    fz = FZ[FZ.well == w]
    ax.plot(fz.T_degC, fz.P_bar, "s", ms=9, mfc="none", mec=cm(i))
ax.set_xlabel("Temperature (degC)")
ax.set_ylabel("Pressure (bar)")
ax.set_ylim(0, 300)
ax.set_xlim(0, 380)
ax.legend(fontsize=7, ncol=2)
ax.set_title("P-T state of all survey points (squares = feed zones). Points left of the curve are liquid")
fig.tight_layout()
fn = os.path.join(FIG, "PT_state_vs_saturation.png")
fig.savefig(fn, dpi=120)
plt.close(fig)
figs.append(fn)

# --- maps by elevation
MAP_ELEVS = [0, -250, -500, -750, -1000, -1500, -2000]


def map_panel(ax, zc, var, vmin, vmax, cmap, unit):
    sub = CE[CE.elevation_mRL == zc]
    draw_base(ax, grid=False)
    if len(sub) >= 3:
        tri = mtri.Triangulation(sub.x.values, sub.y.values)
        cs = ax.tricontourf(tri, sub[var].values, levels=np.linspace(vmin, vmax, 15), cmap=cmap, alpha=0.8,
                            extend="both")
    sc = ax.scatter(sub.x, sub.y, c=sub[var], cmap=cmap, vmin=vmin, vmax=vmax, s=60, ec="k", zorder=5)
    for _, rr in sub.iterrows():
        ax.annotate(f"{rr.well.replace('WELL-', 'W')}\n{rr[var]:.0f}", (rr.x, rr.y), xytext=(4, 2),
                    textcoords="offset points", fontsize=7)
    draw_wells(ax, label=False, color="0.4")
    ax.set_title(f"{zc} mRL: {var.split('_')[0]} ({unit}); n={len(sub)} wells", fontsize=9)
    return sc


for var, vmin, vmax, cmap, unit, name in [("T_degC", 0, 280, CMAP_T, "degC", "maps_temperature_by_elevation.png"),
                                          ("P_minus_field_median_bar", -6, 6, "RdBu_r", "bar rel. to median at that elevation",
                                           "maps_pressure_by_elevation.png")]:
    fig, axs = plt.subplots(2, 4, figsize=(20, 10.5))
    sc = None
    for ax, zc in zip(axs.flat, MAP_ELEVS):
        sc = map_panel(ax, zc, var, vmin, vmax, cmap, unit)
    axs.flat[-1].axis("off")
    fig.colorbar(sc, ax=axs.flat[-1], fraction=0.8, label=unit)
    fig.suptitle(("Measured temperature" if var == "T_degC" else
                  "Measured pressure minus the median of the wells at the same elevation") +
                 " - linear interpolation inside the well hull only (no extrapolation). "
                 "Black lines: faults; blue: lake and rivers", fontsize=10)
    fig.tight_layout()
    fn = os.path.join(FIG, name)
    fig.savefig(fn, dpi=100)
    plt.close(fig)
    figs.append(fn)


# --- cross-sections
def section(wlist, name, title):
    pts = np.array([[W[w]["track"][0, 0], W[w]["track"][0, 1]] for w in wlist])
    seglen = np.r_[0, np.cumsum(np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1])))]

    def proj(x, y):
        best = (np.inf, 0, 0)
        for i in range(len(pts) - 1):
            p1, p2 = pts[i], pts[i + 1]
            d = p2 - p1
            L2 = d @ d
            t = max(0, min(1, ((x - p1[0]) * d[0] + (y - p1[1]) * d[1]) / L2))
            q = p1 + t * d
            dist = math.hypot(x - q[0], y - q[1])
            if dist < best[0]:
                best = (dist, seglen[i] + t * math.sqrt(L2), dist)
        return best[1], best[2]

    S, Z, T, P = [], [], [], []
    offs = {}
    for w in wlist:
        sub = TGT[TGT.well == w]
        for _, rr in sub.iterrows():
            s, d = proj(rr.x, rr.y)
            S.append(s), Z.append(rr.elevation_mRL), T.append(rr.T_degC), P.append(rr.P_bar)
            offs[w] = max(offs.get(w, 0), d)
    S, Z, T = map(np.array, (S, Z, T))
    si = np.linspace(S.min(), S.max(), 300)
    zi = np.linspace(Z.min(), Z.max(), 300)
    SI, ZI = np.meshgrid(si, zi)
    TI = griddata((S, Z * 1.0), T, (SI, ZI), method="linear")
    fig, ax = plt.subplots(figsize=(14, 7))
    cs = ax.contourf(SI, ZI, TI, levels=np.arange(0, 290, 20), cmap=CMAP_T, alpha=0.85)
    cl = ax.contour(SI, ZI, TI, levels=[100, 150, 200, 250], colors="w", linewidths=0.8)
    ax.clabel(cl, fmt="%d", fontsize=8)
    plt.colorbar(cs, ax=ax, label="T (degC), linear interpolation between wells, no extrapolation")
    for w in wlist:
        r = W[w]
        sub = TGT[TGT.well == w]
        ss = [proj(x, y)[0] for x, y in zip(sub.x, sub.y)]
        ax.plot(ss, sub.elevation_mRL, "k-", lw=1)
        sh = [proj(*track_xy(r, zz))[0] for zz in [r["wellhead"], r["casing_shoe"]]]
        ax.plot(sh, [r["wellhead"], r["casing_shoe"]], "k-", lw=3)
        for f in r["fz"]:
            sf = [proj(*track_xy(r, zz))[0] for zz in [f["top"], f["bot"]]]
            ax.plot(sf, [f["top"], f["bot"]], "-", color="cyan", lw=6)
        bm = sub[sub.boiling_margin_degC < BM_FLAG]
        if len(bm):
            ax.plot([proj(x, y)[0] for x, y in zip(bm.x, bm.y)], bm.elevation_mRL, "o", mfc="none", mec="lime",
                    ms=8)
        ax.annotate(f"{w}\n(off-line {offs[w]:.0f} m)", (proj(r['track'][0, 0], r['track'][0, 1])[0],
                                                          r["wellhead"]), xytext=(0, 8),
                    textcoords="offset points", ha="center", fontsize=8, weight="bold")
    # faults crossing the section line
    for k, pl in FAULTS.items():
        for i in range(len(pts) - 1):
            p = seg_intersect(pts[i], pts[i + 1], pl[0], pl[1])
            if p is not None:
                sf = proj(*p)[0]
                ax.axvline(sf, color="k", ls="--", lw=1)
                ax.annotate(f"{k} (trace)", (sf, Z.max() + 30), rotation=90, fontsize=7, va="bottom")
    ax.plot([], [], "k-", lw=3, label="casing")
    ax.plot([], [], "-", color="cyan", lw=6, label="feed zone")
    ax.plot([], [], "o", mfc="none", mec="lime", label=f"boiling margin < {BM_FLAG:g} degC")
    ax.legend(loc="lower left", fontsize=8)
    ax.set_xlabel("Distance along section (m)")
    ax.set_ylabel("Elevation (mRL)")
    ax.set_title(title)
    ax.set_ylim(Z.min() - 50, Z.max() + 250)
    ax.set_xlim(S.min() - 200, S.max() + 200)
    fig.tight_layout()
    fn = os.path.join(FIG, name)
    fig.savefig(fn, dpi=110)
    plt.close(fig)
    figs.append(fn)


hot = SUM.sort_values("Tmax_degC", ascending=False).well.tolist()
section(["WELL-2", "WELL-10", "WELL-3", "WELL-6", "WELL-7"], "section_hot_axis_W2-W10-W3-W6-W7.png",
        "Section A (NW-SE) through the hottest wells: WELL-2, WELL-10, WELL-3, WELL-6, WELL-7")
section(["WELL-8", "WELL-6", "WELL-3", "WELL-5", "WELL-1"], "section_west_east_W8-W6-W3-W5-W1.png",
        "Section B (W-E): WELL-8, WELL-6, WELL-3, WELL-5, WELL-1")
section(["WELL-2", "WELL-9", "WELL-4"], "section_north_W2-W9-W4.png",
        "Section C (towards the north-east): WELL-2, WELL-9, WELL-4")

# --- hydrostatic excess figure
fig, ax = plt.subplots(figsize=(9, 6))
for i, w in enumerate(WELLS):
    h = HY[HY.well == w]
    for _, rr in h.iterrows():
        ax.plot([rr.ratio_measured_to_IF97] * 2, [rr.z_top_mRL, rr.z_bottom_mRL], "-", color=cm(i), lw=3 if
                rr.window_type.startswith("full") else 1.5, alpha=0.9)
    ax.plot([], [], color=cm(i), label=w)
ax.axvline(1.0, color="k", lw=1, label="IF97 static column, g=9.80665")
ax.axvline(10.0 / G, color="r", ls="--", lw=1, label="ratio for g=10.0 m/s2")
ax.axvline(FIELD_RATIO, color="gray", ls=":", lw=1.5, label=f"cold-well median {FIELD_RATIO:.4f}")
ax.set_xlabel("Measured dP / IF97-predicted dP")
ax.set_ylabel("Elevation (mRL)")
ax.legend(fontsize=7)
ax.set_title("Hydrostatic test, single-phase liquid windows >= 300 m (thick = full run)")
fig.tight_layout()
fn = os.path.join(FIG, "hydrostatic_ratio_by_window.png")
fig.savefig(fn, dpi=120)
plt.close(fig)
figs.append(fn)

# --- natural state pressure test figure
fig, axs = plt.subplots(1, 4, figsize=(16, 4.5))
for ax, (_, rr) in zip(axs, PC.iterrows()):
    w = rr.well
    dates, vals = W[w]["transient_P"][rr.gauge_elevation_mRL]
    yrs = [datetime.strptime(d, "%d/%m/%Y") for d in dates]
    ax.plot(yrs, vals, "o-", label="gauge")
    ax.axhline(rr.survey_P_at_gauge_bar, color="r", ls="--", label="survey P at gauge elev.")
    ax.axvline(datetime.strptime(W[w]["survey_date_P"], "%d/%m/%Y"), color="gray", ls=":", label="survey date")
    ax.set_title(f"{w} at {rr.gauge_elevation_mRL:g} mRL: diff {rr.diff_survey_minus_first_bar:+.1f} bar", fontsize=9)
    ax.set_ylabel("P (bar)")
    ax.tick_params(axis="x", rotation=45)
    ax.legend(fontsize=7)
fig.tight_layout()
fn = os.path.join(FIG, "natural_state_pressure_test.png")
fig.savefig(fn, dpi=110)
plt.close(fig)
figs.append(fn)

# --- feed-zone pressure vs elevation
fig, ax = plt.subplots(figsize=(8, 7))
for i, w in enumerate(WELLS):
    f = FZ[FZ.well == w]
    ax.plot(f.P_bar, f.mid_mRL, "s", color=cm(i), ms=8, label=w)
zz = np.linspace(FZ.mid_mRL.min(), FZ.mid_mRL.max(), 50)
ax.plot(np.polyval(fzfit, zz), zz, "k--", lw=1, label="single field quadratic fit")
ax.set_xlabel("Feed-zone pressure (bar)")
ax.set_ylabel("Feed-zone mid elevation (mRL)")
ax.legend(fontsize=7)
ax.set_title("Feed-zone pressures: all on one curve within a few bar -> connected deep pressure")
fig.tight_layout()
fn = os.path.join(FIG, "feedzone_pressure_vs_elevation.png")
fig.savefig(fn, dpi=120)
plt.close(fig)
figs.append(fn)

# --- conceptual model sketch (plan view, built from the data tables)
fig, axs = plt.subplots(1, 2, figsize=(17, 8))
ax = axs[0]
draw_base(ax, grid=False)
s500 = CE[CE.elevation_mRL == -500]
tri = mtri.Triangulation(s500.x.values, s500.y.values)
cs = ax.tricontourf(tri, s500.T_degC.values, levels=np.arange(0, 290, 20), cmap=CMAP_T, alpha=0.6)
plt.colorbar(cs, ax=ax, label="T at -500 mRL (degC), interpolated inside well hull", shrink=0.7)
draw_wells(ax)
GS = SUM.merge(GEO[["well", "x_wellhead", "y_wellhead"]], on="well")
up = GS[GS.regime.str.startswith("upflow")]
for _, rr in up.iterrows():
    ax.add_patch(plt.Circle((rr.x_wellhead, rr.y_wellhead), 450, fc="none", ec="red", lw=2.5, ls="-"))
    ax.annotate("UPFLOW\n(near boiling)", (rr.x_wellhead, rr.y_wellhead), xytext=(-110, 30),
                textcoords="offset points", color="red", weight="bold", fontsize=9)
src = up[["x_wellhead", "y_wellhead"]].mean().values if len(up) else GS.loc[GS.Tmax_degC.idxmax(),
                                                                              ["x_wellhead", "y_wellhead"]].values
for _, rr in GS[GS.regime.str.contains("outflow")].iterrows():
    ax.add_patch(FancyArrowPatch(src, (rr.x_wellhead, rr.y_wellhead), arrowstyle="-|>", mutation_scale=18,
                                 color="darkorange", lw=2, shrinkA=10, shrinkB=10))
    ax.annotate(f"Tmax {rr.Tmax_degC:.0f} @ {rr.Tmax_elev_mRL:.0f}", (rr.x_wellhead, rr.y_wellhead),
                xytext=(6, -14), textcoords="offset points", fontsize=7, color="darkorange")
for _, rr in GS[GS.regime.str.startswith("recharge")].iterrows():
    ax.annotate("RECHARGE /\ncold downflow", (rr.x_wellhead, rr.y_wellhead), xytext=(10, -30),
                textcoords="offset points", color="blue", weight="bold", fontsize=8)
    ax.plot(rr.x_wellhead, rr.y_wellhead, "v", ms=16, color="blue")
for _, rr in GS[GS.regime == "conductive"].iterrows():
    ax.annotate("conductive\nmargin", (rr.x_wellhead, rr.y_wellhead), xytext=(10, -30),
                textcoords="offset points", color="green", fontsize=8)
ax.set_title("Conceptual model, plan view (data-derived; arrows = inferred outflow paths)")
# schematic section
ax = axs[1]
order = ["WELL-8", "WELL-6", "WELL-10", "WELL-2", "WELL-3", "WELL-9", "WELL-7", "WELL-5", "WELL-1", "WELL-4"]
order = [w for w in SUM.sort_values("Tmax_degC").well if w] if False else order
for i, w in enumerate(order):
    r = W[w]
    sc = ax.scatter(np.full(len(r["zT"]), i), r["zT"], c=r["T"], cmap=CMAP_T, vmin=0, vmax=280, s=40)
    s = SUM[SUM.well == w].iloc[0]
    ax.annotate(f"{w.replace('WELL-', 'W')}\n{s.regime.split(' (')[0].split(';')[0]}", (i, r["wellhead"] + 30),
                ha="center", fontsize=7, rotation=0)
    if s.inversion_degC > 0:
        ax.annotate("", xy=(i + 0.35, s.Tmax_elev_mRL), xytext=(i - 0.35, s.Tmax_elev_mRL),
                    arrowprops=dict(arrowstyle="->", color="darkorange", lw=2))
    for f in r["fz"]:
        ax.plot([i, i], [f["top"], f["bot"]], "-", color="cyan", lw=7, alpha=0.6)
    ax.plot([i - 0.25, i + 0.25], [r["casing_shoe"]] * 2, "k-", lw=2)
plt.colorbar(sc, ax=ax, label="Measured T (degC)", shrink=0.7)
ax.set_xticks(range(len(order)))
ax.set_xticklabels([w.replace("WELL-", "W") for w in order])
ax.set_ylabel("Elevation (mRL)")
ax.set_title("Wells ordered W->E-ish; orange arrows = lateral outflow at Tmax elevation;\n"
             "cyan = feed zones, black bar = casing shoe")
fig.tight_layout()
fn = os.path.join(FIG, "conceptual_model_sketch.png")
fig.savefig(fn, dpi=110)
plt.close(fig)
figs.append(fn)
print(f"\n{len(figs)} figures written to Output/figures/")

# ============================================================================
# 16. Report (defined below, called at the end)
# ============================================================================
def write_report():
    """Write Output/NS_data_interpretation.md. Every number is taken from the tables written above."""
    S = SUM.set_index("well")
    Gd = GEO.set_index("well")
    F = FZ.copy()
    K = lambda w: w.replace("WELL-", "W")

    def n0(x):
        return "n/a" if pd.isna(x) else f"{x:.0f}"

    def n1(x):
        return "n/a" if pd.isna(x) else f"{x:.1f}"

    # ---- derived quantities used in the text (all from output tables)
    up = S[S.regime.str.startswith("upflow")].index.tolist()
    outf = S[S.regime.str.contains("outflow")].index.tolist()
    rech = S[S.regime.str.startswith("recharge")].index.tolist()
    cond = S[S.regime == "conductive"].index.tolist()
    w10 = S.loc["WELL-10"]
    dFX = Gd["dist_Fault1xFault3_m"]
    x10, y10 = Gd.loc["WELL-10", ["x_wellhead", "y_wellhead"]]
    dist_from_up = {w: math.hypot(Gd.loc[w, "x_wellhead"] - x10, Gd.loc[w, "y_wellhead"] - y10) for w in WELLS}
    bm10 = TGT[(TGT.well == "WELL-10") & (TGT.boiling_margin_degC < BM_FLAG)]
    fz10 = F[F.well == "WELL-10"].iloc[0]
    hb = BANDSUM.set_index(["well", "band"])
    sh_band = "shallow (bottom >= -450 mRL)"
    dp_band = "deep (top <= -450 mRL)"
    shallow_ratio = {w: hb.loc[(w, sh_band), "ratio_median"] for w in WELLS if (w, sh_band) in hb.index}
    deep_ratio = {w: hb.loc[(w, dp_band), "ratio_median"] for w in WELLS if (w, dp_band) in hb.index}
    hot_sh = {w: v for w, v in shallow_ratio.items() if w not in cold}
    cold_sh = {w: v for w, v in shallow_ratio.items() if w in cold}
    steep = HY[HY.steeper_than_cold_water & (HY.window_type != "full single-phase run")]
    pcw = PC.set_index("well")
    res = F.P_residual_to_field_fit_bar
    ce = CE.set_index(["elevation_mRL", "well"])
    shallow_dev = CE[CE.elevation_mRL.isin([0, -250])].groupby("well").P_minus_field_median_bar.mean()
    lay_n = LC.groupby("well").size()
    above_n = int(TGT.above_model_surface.sum())
    n_pts = len(TGT)
    integ_bad = INTEG[INTEG.result != "OK"]
    fzq = F.groupby("quality").size().to_dict()
    hi_fz = F[F.quality == "HIGH"].sort_values("mid_mRL", ascending=False)
    cap_hot = S.loc[[w for w in WELLS if S.loc[w, "Tmax_degC"] >= 200]]
    w2 = S.loc["WELL-2"]
    els_all = [0, -250, -500, -750, -1000, -1500, -2000]
    deepfz = F[(F.mid_mRL <= -1000) & (F.well != "WELL-8")]
    z240 = {w: first_crossing(W[w]["zT"], W[w]["T"], 240.0) for w in ["WELL-3", "WELL-10", "WELL-2"]}
    capth = {w: W[w]["wellhead"] - z240[w] for w in z240}
    all_blocks_match = bool(F.block_matches_wells_json.all())
    # isotherms vs casing shoe
    iso_behind = {w: (S.loc[w, "z_100C_mRL"] >= W[w]["casing_shoe"], S.loc[w, "z_200C_mRL"] >= W[w]["casing_shoe"])
                  for w in cap_hot.index}
    iso_open = [f"{K(w)} 200 degC at {S.loc[w, 'z_200C_mRL']:.0f} mRL (shoe {W[w]['casing_shoe']:.0f})"
                for w, (a, b) in iso_behind.items() if not b]
    w8iso = W["WELL-8"]["T"][(W["WELL-8"]["zT"] <= -550) & (W["WELL-8"]["zT"] >= -950)]
    within2 = {e: int((CE[CE.elevation_mRL == e].P_minus_field_median_bar.abs() <= 2.0).sum()) for e in els_all}
    ntot = {e: int((CE.elevation_mRL == e).sum()) for e in els_all}
    tmax_vs_fz = {}
    for w in ["WELL-3", "WELL-6", "WELL-7", "WELL-9", "WELL-5"]:
        zt = S.loc[w, "Tmax_elev_mRL"]
        d = min((0 if f["bot"] <= zt <= f["top"] else min(abs(zt - f["top"]), abs(zt - f["bot"]))) for f in W[w]["fz"])
        tmax_vs_fz[w] = (d, zt >= W[w]["casing_shoe"])
    ddrop = PC[PC.first_transient_date.str.endswith("2014")].transient_drop_first_year_bar
    dat_ok = not INTEG[INTEG.check.str.startswith("json_vs_dat")].result.ne("OK").any()

    L = []
    A = L.append
    A("# Natural-state data interpretation: TE14515 (Teaching_Model_3D)")
    A("")
    A(f"Generated by `Scripts/natural_state_analysis.py` on {datetime.now().strftime('%Y-%m-%d %H:%M')}. "
      "Every number below is read from a file in `Output/` (file named in brackets). "
      "**Observation** = what the data show; **Interpretation** = an inference, with its confidence. "
      "Units: mRL (metres relative to sea level, positive up), bar (absolute), degC, kJ/kg, kg/s. "
      "Content about the reference model TE14515 is labelled **[reference model]** and is never used as evidence.")
    A("")
    # ------------------------------------------------------------------ 1
    A("## 1. Key conclusions for calibration")
    A("")
    out_rows = ", ".join(f"{K(w)} {S.loc[w, 'Tmax_degC']:.1f} degC at {S.loc[w, 'Tmax_elev_mRL']:.0f} mRL "
                         f"({dist_from_up[w] / 1000:.1f} km)" for w in sorted(outf, key=lambda w: dist_from_up[w]))
    A(f"1. **The upflow is at WELL-10 / WELL-2, next to the Fault 1 x Fault 3 intersection.** Observation: WELL-10 "
      f"is the hottest well ({w10.Tmax_degC:.1f} degC at {w10.Tmax_elev_mRL:.0f} mRL) and the only one near boiling: "
      f"boiling margin {bm10.boiling_margin_degC.min():.1f}-{bm10.boiling_margin_degC.max():.1f} degC from "
      f"{bm10.elevation_mRL.max():.0f} to {bm10.elevation_mRL.min():.1f} mRL; the profile is isothermal below its feed zone. "
      f"WELL-2 is {dFX['WELL-2']:.0f} m from the mapped Fault 1 x Fault 3 intersection; WELL-10 is "
      f"{Gd.loc['WELL-10', 'dist_Fault1_m']:.0f} m from Fault 1 and {Gd.loc['WELL-10', 'dist_Fault3_m']:.0f} m from Fault 3 "
      f"({dFX['WELL-10']:.0f} m from the intersection). WELL-2 reheats with depth, to {w2.T_at_TD_degC:.1f} degC at TD ({w2.TD_mRL:.0f} mRL), "
      f"the hottest deep temperature in the field. "
      f"Action: put the deep hot mass input under the F1 x F3 intersection (WELL-2/WELL-10 blocks), with enthalpy of at "
      f"least that of liquid at {w10.Tmax_degC:.0f} degC (about {h_liquid(w10.Tmax_degC, fz10.P_bar):.0f} kJ/kg from "
      f"IF97 at the WELL-10 feed-zone pressure). [well_summary.csv, well_geometry_check.csv, NS_calibration_targets.csv]")
    A(f"2. **A hot outflow layer spreads east, south-east and south at about -250 to -750 mRL.** Observation: every "
      f"well except WELL-10 (upflow), the recharge wells and the eastern margin well WELL-1 has a temperature inversion, and Tmax falls with distance from WELL-10: "
      f"{out_rows}. Below the inversions the temperature drops to "
      f"{S.loc[outf, 'T_inversion_min_degC'].min():.1f}-{S.loc[outf, 'T_inversion_min_degC'].max():.1f} degC. "
      f"Interpretation (high confidence): lateral outflow of hot fluid above cooler water. Action: the model needs a "
      f"laterally permeable layer at about -250 to -750 mRL that connects the upflow to WELL-3, WELL-6, WELL-9, WELL-7 "
      f"and WELL-5, and lower vertical permeability below it, so that the cool water under the outflow stays cool. "
      f"[well_summary.csv, figures/section_hot_axis_*.png]")
    A(f"3. **Cold recharge or downflow comes in at WELL-4 (north-east, on Fault 3) and WELL-8 (south-west).** "
      f"Observation: WELL-4 reaches only {S.loc['WELL-4', 'Tmax_degC']:.1f} degC at TD ({S.loc['WELL-4', 'TD_mRL']:.0f} mRL), "
      f"with a gradient of {S.loc['WELL-4', 'shallow_gradient_degC_per_km']:.0f}-{S.loc['WELL-4', 'deep_gradient_degC_per_km']:.0f} "
      f"degC/km, and is {Gd.loc['WELL-4', 'dist_Fault3_m']:.0f} m from Fault 3. "
      f"WELL-8 reaches only {S.loc['WELL-8', 'Tmax_degC']:.1f} degC at {S.loc['WELL-8', 'TD_mRL']:.0f} mRL. Action: the model needs cold "
      f"water ({W['WELL-4']['T'].min():.0f}-{W['WELL-4']['T'].max():.0f} degC in WELL-4, {w8iso.min():.0f}-{w8iso.max():.0f} degC "
      f"at -550 to -950 mRL in WELL-8) entering at depth from the north-east along Fault 3 and from the south-west. "
      f"A purely conductive base would make these wells far too hot. [well_summary.csv]")
    A(f"4. **The cap is thin and highly conductive over the hot zone.** Observation: in the wells with Tmax >= 200 degC "
      f"the 100 degC isotherm is at {cap_hot.z_100C_mRL.min():.0f} to {cap_hot.z_100C_mRL.max():.0f} mRL, "
      f"and the shallow gradients are {cap_hot.shallow_gradient_degC_per_km.min():.0f}-"
      f"{cap_hot.shallow_gradient_degC_per_km.max():.0f} degC/km. That is a conductive heat loss of "
      f"{cap_hot.conductive_flux_cap_Wm2_assumed_k2.min():.1f}-{cap_hot.conductive_flux_cap_Wm2_assumed_k2.max():.1f} W/m2, "
      f"assuming k = {K_ROCK} W/m/K. Action: the model needs a low-permeability layer from the surface down to about "
      f"-100 to -250 mRL above the outflow. The shallow temperatures must be matched by conduction through it, not by "
      f"through-flow. [well_summary.csv]")
    A(f"5. **The deep pressures are connected across the field.** Observation: all {len(F)} feed-zone pressures lie "
      f"within {res.abs().max():.1f} bar of a single pressure-elevation curve (RMS {np.sqrt((res ** 2).mean()):.1f} bar). "
      f"The largest positive residuals are in the recharge wells (WELL-4 {F[F.well == 'WELL-4'].P_residual_to_field_fit_bar.iloc[0]:+.1f}, "
      f"WELL-8 {F[F.well == 'WELL-8'].P_residual_to_field_fit_bar.iloc[0]:+.1f} bar), and the most negative are in the "
      f"eastern margin well WELL-1 ({F[F.well == 'WELL-1'].P_residual_to_field_fit_bar.min():+.1f} bar). Interpretation "
      f"(medium confidence): a single connected reservoir, with a small head gradient from the recharge sides towards the "
      f"east, where the river is. Action: use the feed-zone pressures as the pressure targets "
      f"(feedzones.csv), with a tolerance of about 1-2 bar. [feedzones.csv, figures/feedzone_pressure_vs_elevation.png]")
    A(f"6. **Above about -450 mRL the survey pressures are not a static column (data problem, high impact).** "
      f"Observation: the measured pressure gradient is steeper than an IF97 static column of the measured temperature "
      f"profile in every well. The median ratio is {min(cold_sh.values()):.3f}-{max(cold_sh.values()):.3f} in the cold wells "
      f"and {min(hot_sh.values()):.2f}-{max(hot_sh.values()):.2f} in the hot wells, and in {len(steep)} windows the "
      f"gradient is steeper than a column of cold water (implied density up to {steep.rho_implied_kgm3.max():.0f} kg/m3). "
      f"Below -450 mRL most windows are within a few percent of the IF97 static column. Checks: a datum offset cannot "
      f"change a gradient; a gravity or unit convention would act at every depth; two-phase fluid would make the column "
      f"lighter, not heavier; the 0.1 bar rounding is below 0.5 % over 300 m. "
      f"Action: give zero weight to pressures that are not at feed zones (column weight_P), and ask the data provider "
      f"how the shallow pressures were obtained. [hydrostatic_check.csv, hydrostatic_band_summary.csv]")
    A(f"7. **The surveys represent the pre-production natural state, even though every survey is dated 29/07/2015.** "
      f"Observation: at the gauge elevation, the survey pressure differs from the first (01/01/2014) gauge reading by "
      f"{pcw.loc['WELL-2', 'diff_survey_minus_first_bar']:+.1f} bar at WELL-2 and "
      f"{pcw.loc['WELL-3', 'diff_survey_minus_first_bar']:+.1f} bar at WELL-3. The 01/01/2015 readings are already "
      f"{pcw.loc['WELL-2', 'transient_drop_first_year_bar']:.1f} bar and {pcw.loc['WELL-3', 'transient_drop_first_year_bar']:.1f} bar lower than the first. "
      f"WELL-9 differs by {pcw.loc['WELL-9', 'diff_survey_minus_first_bar']:+.1f} bar. "
      f"Action: calibrate the natural state against all surveys, with the confidence multipliers in section 6. "
      f"[pressure_consistency_check.csv]")
    A(f"8. **The model surface is above several wellheads.** Observation: the wellhead is below the model surface of its "
      f"column at WELL-1 ({Gd.loc['WELL-1', 'wellhead_minus_surface_m']:.0f} m), WELL-8 "
      f"({Gd.loc['WELL-8', 'wellhead_minus_surface_m']:.0f} m), WELL-7 ({Gd.loc['WELL-7', 'wellhead_minus_surface_m']:.0f} m), "
      f"WELL-5 ({Gd.loc['WELL-5', 'wellhead_minus_surface_m']:.0f} m) and WELL-9 "
      f"({Gd.loc['WELL-9', 'wellhead_minus_surface_m']:.0f} m). No survey point is above the model surface. Action: correct "
      f"SURFA for these columns, or accept that the top-layer model temperatures sit tens of metres above the real ground "
      f"surface. [well_geometry_check.csv]")
    A(f"9. **[reference model] The starting model has no deep mass upflow.** It has {REF.get('n_heat_gen')} basal HEAT "
      f"generators totalling {REF.get('total_heat_MW')} MW and {REF.get('n_mass_gen')} near-surface rain MASS "
      f"generators. None of its MASS generators is in the deepest three layers. Conclusions 1 to 3 imply a hot mass "
      f"upflow and cold lateral inflow. The starting model includes neither. [key_numbers.json]")
    A("")
    # ------------------------------------------------------------------ 2
    A("## 2. The data")
    A("")
    A(f"**Inventory.** There are 10 wells (Groups_data.json). Each has one temperature survey and one pressure survey "
      f"at the same {n_pts} elevations in total, plus a well track, the casing and liner design, and feed zones "
      f"({', '.join(f'{k}: {v}' for k, v in fzq.items())}). WELL-2, 3, 9 and 10 have production rates, enthalpy and "
      f"pressure gauge series. WELL-4, 5 and 8 have injection rates at {list(W['WELL-4']['inj_h'].values())[0]:.1f} kJ/kg. "
      f"WELL-1, 6 and 7 have no flow data. WELL-8 and WELL-9 are deviated, with offsets at TD of "
      f"{Gd.loc['WELL-8', 'horizontal_offset_at_TD_m']:.0f} m and {Gd.loc['WELL-9', 'horizontal_offset_at_TD_m']:.0f} m.")
    A("")
    A(f"**Checks** (data_integrity_check.csv, {len(INTEG)} checks):")
    A("- JSON and dat files " + ("match exactly for every temperature, pressure and transient-pressure value (difference 0), "
      "and the T and P elevations are identical in each well." if dat_ok else "do NOT all match: see data_integrity_check.csv."))
    A("- Every survey point lies within the well track. Every feed zone lies within the track and the slotted liner, "
      "below the casing shoe. " + ("Every feed-zone block computed here matches the block in `gTE14515_wells.json` (feedzones.csv)."
      if all_blocks_match else "Some feed-zone blocks differ from `gTE14515_wells.json` (feedzones.csv)."))
    for _, b in integ_bad.iterrows():
        A(f"- **{b.well}: {b.check} = {b.result}.** {b.detail}. Treat the track TD as correct, since the survey ends there.")
    A(f"- The first survey point is well below the wellhead at WELL-4 "
      f"({Gd.loc['WELL-4', 'survey_top_below_wellhead_m']:.0f} m) and WELL-10 "
      f"({Gd.loc['WELL-10', 'survey_top_below_wellhead_m']:.0f} m).")
    A(f"- IAPWS-IF97 (Region 4 saturation, Region 1 density and enthalpy) is implemented in the script and validated "
      f"against the IF97 verification tables. All {len(IF97_CHECKS)} checks passed, including Tsat(1.01325 bar) = "
      f"{tsat_c(1.01325):.2f} degC and Tsat(100 bar) = {tsat_c(100):.2f} degC. The script stops if any Tsat error "
      f"exceeds 0.1 K (run_console_log.txt).")
    A("")
    A("**Natural-state timing** (pressure_consistency_check.csv, figures/natural_state_pressure_test.png). "
      "All surveys are dated 29/07/2015. Production and injection began on 01/01/2014.")
    A("")
    A("| Well | Role | Gauge elev. (mRL) | Survey P at gauge (bar) | First gauge P (bar, date) | Survey minus first (bar) | Verdict | Confidence | Weight multiplier |")
    A("|---|---|---|---|---|---|---|---|---|")
    for w in WELLS:
        if w in pcw.index:
            r = pcw.loc[w]
            A(f"| {w} | {NS[w]['role']} | {r.gauge_elevation_mRL:.0f} | {r.survey_P_at_gauge_bar:.1f} | "
              f"{r.first_transient_P_bar:.1f} ({r.first_transient_date}) | {r.diff_survey_minus_first_bar:+.1f} | "
              f"{NS[w]['verdict']} | {NS[w]['confidence']} | {CONF_W[NS[w]['confidence']]} |")
        else:
            A(f"| {w} | {NS[w]['role']} | - | - | - | - | {NS[w]['verdict']} | {NS[w]['confidence']} | "
              f"{CONF_W[NS[w]['confidence']]} |")
    A("")
    A("Reasons for each verdict (well_summary.csv `ns_verdict`):")
    for w in WELLS:
        A(f"- {w}: {NS[w]['reason']}.")
    A(f"- **Injectors.** If WELL-4, 5 or 8 had been surveyed while injecting, the open hole would show an isothermal "
      f"column near the injectate temperature (about 80 degC). WELL-4 reads {S.loc['WELL-4', 'Tmax_degC']:.0f} degC or "
      f"less, far colder than the injectate, so it is not an injection profile. WELL-5 is isothermal at about 155-158 "
      f"degC behind casing, above any injectate temperature. WELL-8 has an isothermal section at 67-68 degC from -650 to "
      f"-950 mRL, around the casing shoe at {W['WELL-8']['casing_shoe']:.0f} mRL. This is the only feature that could "
      f"be cooled injectate. Its points remain in the targets at the WELL-8 confidence, and it is listed in section 7.")
    A("")
    # ------------------------------------------------------------------ 3
    A("## 3. Temperature structure")
    A("")
    A("Per-well characteristics (well_summary.csv; profiles in figures/profile_WELL-*.png). Gradients are least-squares "
      "fits. The shallow gradient uses the top 250 m of the liquid-filled survey; the deep gradient uses the interval "
      "given in the CSV. Isotherms are the first downward crossing.")
    A("")
    A("| Well | Tmax (degC) @ mRL | T at TD (degC) @ mRL | Inversion (degC), range mRL | Shallow grad. (degC/km) | Deep grad. (degC/km) | z100 (mRL) | z200 (mRL) | Min boiling margin (degC) @ mRL | Water level / pivot | Regime |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for w in WELLS:
        r = S.loc[w]
        inv = (f"{r.inversion_degC:.1f}, {r.inversion_top_mRL:.0f} to {r.inversion_base_mRL:.0f}"
               if r.inversion_degC > 0 else "none")
        A(f"| {w} | {r.Tmax_degC:.1f} @ {r.Tmax_elev_mRL:.0f} | {r.T_at_TD_degC:.1f} @ {r.TD_mRL:.0f} | {inv} | "
          f"{n0(r.shallow_gradient_degC_per_km)} | {n0(r.deep_gradient_degC_per_km)} | {n0(r.z_100C_mRL)} | "
          f"{n0(r.z_200C_mRL)} | {r.min_boiling_margin_degC:.1f} @ {r.min_bm_elev_mRL:.0f} | "
          f"WL {r.water_level_mRL:.0f} mRL; pivot {r.pivot_P_bar:.1f} bar @ {r.pivot_fz_mid_mRL:.0f} ({r.pivot_fz_quality}) | "
          f"{r.regime} |")
    A("")
    A("**Reading the profiles against the wellbore** (all wells; see section 7 for the open questions):")
    A(f"- **Behind casing** ({int(TGT.behind_casing.sum())} of {n_pts} points), no internal flow is possible. The steep "
      f"shallow gradients and the 100 degC isotherms are behind casing in every well with Tmax >= 200 degC, and so is the "
      f"200 degC isotherm" + (f" except {'; '.join(iso_open)}" if iso_open else "") + ". These are taken as "
      f"conductive formation temperatures.")
    A("- **Inversions.** Distance from the inversion maximum to the nearest feed zone: "
      + ", ".join(f"{K(w)} {d:.0f} m" + (" (behind casing)" if bc else "") for w, (d, bc) in tmax_vs_fz.items())
      + ". In WELL-5 the maximum is behind casing, where no internal flow is possible. In all five the cooling "
      f"continues to or below the deepest feed zone, where a static well has no internal flow. In WELL-3 the minimum ({S.loc['WELL-3', 'T_inversion_min_degC']:.1f} degC at "
      f"{S.loc['WELL-3', 'inversion_base_mRL']:.0f} mRL) is about 800 m below its deepest feed zone. Interpretation "
      f"(high confidence): the inversions are formation features (lateral outflow), not wellbore artefacts.")
    segs = INTFLOW["WELL-2"]
    if segs:
        a, b, ta, tb = segs[0]
        A(f"- **WELL-2 isothermal section.** Between {a:.0f} and {b:.0f} mRL, in open hole between feed zones 2 "
          f"(-720 to -770 mRL) and 3 (-1900 to -1920 mRL), the temperature is {ta:.1f}-{tb:.1f} degC. This "
          f"is the signature of internal upflow from FZ3 to FZ2 inside the wellbore. Interpretation (medium confidence): "
          f"the formation between these depths may be cooler than measured. The points are down-weighted to 0.4.")
    A(f"- **Two-phase or boiling-point-for-depth.** Only WELL-10 comes within {BM_FLAG:g} degC of boiling "
      f"({len(bm10)} points, {bm10.elevation_mRL.max():.0f} to {bm10.elevation_mRL.min():.1f} mRL, behind casing and at "
      f"the feed zone). Interpretation (medium confidence): a two-phase or near-boiling zone at the top of the upflow, "
      f"from about -75 to -265 mRL. The next-closest well is WELL-3, with a minimum margin of "
      f"{S.loc['WELL-3', 'min_boiling_margin_degC']:.1f} degC at {S.loc['WELL-3', 'min_bm_elev_mRL']:.0f} mRL. "
      f"Boiling margins at every point are in NS_calibration_targets.csv and figures/PT_state_vs_saturation.png.")
    A(f"- **Recharge wells.** WELL-4 stays at {W['WELL-4']['T'].min():.1f}-{W['WELL-4']['T'].max():.1f} degC over its "
      f"{W['WELL-4']['zT'].max() - W['WELL-4']['zT'].min():.0f} m survey. WELL-8 stays at 64-68 degC from -550 to "
      f"-950 mRL, then warms conductively at {S.loc['WELL-8', 'deep_gradient_degC_per_km']:.0f} degC/km below its feed "
      f"zone. Both show cold water moving down or sideways. Neither is a conductive profile.")
    A(f"- **WELL-1 (east)** is monotonic, with a shallow gradient of {S.loc['WELL-1', 'shallow_gradient_degC_per_km']:.0f} degC/km and "
      f"{S.loc['WELL-1', 'T_at_TD_degC']:.1f} degC at TD. It is conductive, at the edge of the field. Its wellhead pressure is "
      f"{W['WELL-1']['P'][0]:.1f} bar at {W['WELL-1']['wellhead']:.0f} mRL, which makes it artesian (equivalent water "
      f"level {S.loc['WELL-1', 'water_level_mRL']:.0f} mRL). It is {Gd.loc['WELL-1', 'dist_river_m']:.0f} m from the river.")
    A("")
    A("**Temperature at common elevations** (T_P_at_common_elevations.csv; maps in figures/maps_temperature_by_elevation.png; "
      "interpolation only, no extrapolation):")
    A("")
    els = [0, -250, -500, -750, -1000, -1500, -2000]
    A("| Well | " + " | ".join(f"{e} mRL" for e in els) + " |")
    A("|---|" + "---|" * len(els))
    for w in WELLS:
        A(f"| {w} | " + " | ".join(n1(ce.loc[(e, w), 'T_degC']) if (e, w) in ce.index else "-" for e in els) + " |")
    A("")
    # ------------------------------------------------------------------ 4
    A("## 4. Pressure structure")
    A("")
    A("Measured pressure is a formation pressure only at feed zones. The feed-zone values below are the pressure "
      "targets. Elsewhere, P follows whatever is in the wellbore (see conclusion 6).")
    A("")
    A("| Well | FZ | Mid (mRL) | Quality | T (degC) | P (bar) | Tsat (degC) | Boiling margin (degC) | h liquid (kJ/kg) | Block | Residual to field curve (bar) |")
    A("|---|---|---|---|---|---|---|---|---|---|---|")
    for _, f in F.iterrows():
        A(f"| {f.well} | {f.feedzone} | {f.mid_mRL:.1f} | {f.quality} | {f.T_degC:.1f} | {f.P_bar:.1f} | {n1(f.Tsat_degC)} | "
          f"{n1(f.boiling_margin_degC)} | {f.h_liquid_kJkg:.0f} | {f.block} | {f.P_residual_to_field_fit_bar:+.2f} |")
    A("")
    A("**Pressure at common elevations**, as the deviation from the median of the wells at that elevation "
      "(T_P_at_common_elevations.csv; maps in figures/maps_pressure_by_elevation.png):")
    A("")
    A("| Well | " + " | ".join(f"{e} mRL" for e in els) + " |")
    A("|---|" + "---|" * len(els))
    for w in WELLS:
        A(f"| {w} | " + " | ".join(
            f"{ce.loc[(e, w), 'P_bar']:.1f} ({ce.loc[(e, w), 'P_minus_field_median_bar']:+.1f})" if (e, w) in ce.index else "-"
            for e in els) + " |")
    A("")
    A(f"- Observation: at 0 and -250 mRL the upflow wells have the highest pressures relative to the median "
      f"(mean deviation WELL-10 {shallow_dev['WELL-10']:+.1f}, WELL-2 {shallow_dev['WELL-2']:+.1f}, "
      f"WELL-3 {shallow_dev['WELL-3']:+.1f} bar), and WELL-1 the lowest ({shallow_dev['WELL-1']:+.1f} bar). "
      f"Wells within 2 bar of the median: " + ", ".join(f"{within2[e]}/{ntot[e]} at {e} mRL" for e in els_all)
      + f". WELL-1 is the most persistent low ({ce.loc[(-750, 'WELL-1'), 'P_minus_field_median_bar']:+.1f} bar at -750 mRL). "
      f"WELL-8 is high below its feed zone. Interpretation: the shallow contrast "
      f"partly reflects conclusion 6, so do not use it as a target.")
    A("- The pivot pressure in section 3 is the pressure at the dominant feed zone. Each well has only one survey, so a "
      "pivot point cannot be located by intersecting profiles.")
    A(f"- Observation: below the feed zone of WELL-8 the pressure rises above the median "
      f"({ce.loc[(-2000, 'WELL-8'), 'P_minus_field_median_bar']:+.1f} bar at -2000 mRL). Below its only feed zone "
      f"(-1100 to -1150 mRL) the pressure follows WELL-8's own cold, dense column, so this is a wellbore effect, not "
      f"formation pressure.")
    A("")
    A("**Hydrostatic checks** (hydrostatic_check.csv, hydrostatic_band_summary.csv, figures/hydrostatic_ratio_by_window.png). "
      "The tests use only single-phase liquid intervals: P >= 1.5 bar, boiling margin >= 10 degC, windows >= 300 m. "
      "The measured dP is compared with dP integrated through IF97 Region 1 density over the measured T profile, with "
      f"g = {G}.")
    A("")
    A("| Well | Shallow windows: median ratio | Deep windows: median ratio (min-max) |")
    A("|---|---|---|")
    for w in WELLS:
        if w in shallow_ratio or w in deep_ratio:
            dr = (f"{deep_ratio[w]:.3f} ({hb.loc[(w, dp_band), 'ratio_min']:.3f}-{hb.loc[(w, dp_band), 'ratio_max']:.3f})"
                  if w in deep_ratio else "-")
            A(f"| {w} | {shallow_ratio.get(w, float('nan')):.3f} | {dr} |")
    A("")
    A(f"- **Field-wide component.** The cold wells WELL-1, WELL-4 and WELL-8 have no hot column and no boiling. Their "
      f"full-run median ratio is {FIELD_RATIO:.4f}, equivalent to g = {G * FIELD_RATIO:.2f} m/s2. A uniform "
      f"g = 10 m/s2 convention would give {10 / G:.4f}. Split by depth, the cold-well medians are "
      f"{BAND_RATIO[sh_band]:.4f} in the shallow windows and {BAND_RATIO[dp_band]:.4f} in the deep windows. "
      f"Report: a possible datum or convention issue of about {(FIELD_RATIO - 1) * 100:.1f} % in the cold wells, not physical overpressure. It "
      f"is not uniform with depth, so a single gravity constant does not explain it.")
    A(f"- **Well-specific component.** The shallow excess grows with temperature, from "
      f"{min(cold_sh.values()):.3f}-{max(cold_sh.values()):.3f} in the cold wells to {max(hot_sh.values()):.2f} in "
      f"{max(hot_sh, key=hot_sh.get)}. WELL-2 also has a deep excess "
      f"({deep_ratio.get('WELL-2', float('nan')):.3f}), possibly from flow inside the wellbore through the isothermal "
      f"section (section 3). WELL-5, 6 and 7 below -450 mRL agree with the IF97 static column to within "
      f"{max(abs(deep_ratio[w] - 1) for w in ['WELL-5', 'WELL-6', 'WELL-7']) * 100:.1f} %. This confirms that the "
      f"density calculation and the units are right where the column is static.")
    A("- Interpretation (low confidence): the shallow pressures in the hot wells behave like formation pressures with "
      "upward flow through a low-permeability cap, or like synthetic block pressures. They are not readings in a static "
      "wellbore. Either way, they are not reliable targets.")
    A("")
    # ------------------------------------------------------------------ 5
    A("## 5. Conceptual model (from the data alone)")
    A("")
    A("![conceptual model](figures/conceptual_model_sketch.png)")
    A("")
    A("Cross-sections: figures/section_hot_axis_W2-W10-W3-W6-W7.png (NW-SE through the hottest wells), "
      "figures/section_west_east_W8-W6-W3-W5-W1.png and figures/section_north_W2-W9-W4.png. "
      "Well map: figures/well_map_grid_faults_water.png.")
    A("")
    A("| Feature | Where / depth | Supporting observation | Relation to mapped features | Confidence |")
    A("|---|---|---|---|---|")
    A(f"| Upflow | WELL-10 and WELL-2, from depth up to about -250 mRL | WELL-10 {w10.Tmax_degC:.1f} degC, boiling margin "
      f"{bm10.boiling_margin_degC.min():.1f} degC, isothermal below FZ; WELL-2 reheats to {w2.T_at_TD_degC:.1f} degC at "
      f"{w2.TD_mRL:.0f} mRL | W2 {dFX['WELL-2']:.0f} m from the Fault 1 x Fault 3 intersection; W10 "
      f"{Gd.loc['WELL-10', 'dist_Fault1_m']:.0f} m from F1 and {Gd.loc['WELL-10', 'dist_Fault3_m']:.0f} m from F3 | high for location; medium for depth extent (W10 TD only {w10.TD_mRL:.0f} mRL) |")
    A(f"| Two-phase zone | WELL-10, about -75 to -265 mRL | boiling margin < {BM_FLAG:g} degC at {len(bm10)} points | "
      f"top of the upflow | medium |")
    A(f"| Outflow | -250 to -750 mRL towards E/SE/S (W3, W6, W9, W7, W5) | inversions; Tmax falls and deepens with distance "
      f"from W10 (W3 {S.loc['WELL-3', 'Tmax_elev_mRL']:.0f}, W6 {S.loc['WELL-6', 'Tmax_elev_mRL']:.0f}, "
      f"W7 {S.loc['WELL-7', 'Tmax_elev_mRL']:.0f}, W9 {S.loc['WELL-9', 'Tmax_elev_mRL']:.0f} mRL) | flows away from F1 x F3 "
      f"and reaches WELL-5, east of Fault 2 ({S.loc['WELL-5', 'Tmax_degC']:.0f} degC at {S.loc['WELL-5', 'Tmax_elev_mRL']:.0f} mRL), "
      f"so Fault 2 is not a barrier to it; the shallow HIGH feed zones lie in this layer | high |")
    A(f"| Cool water below outflow | below about -1000 mRL under W3, W6, W7, W9 | T at -1500 mRL: "
      + ", ".join(f"{K(w)} {ce.loc[(-1500, w), 'T_degC']:.0f}" for w in ['WELL-3', 'WELL-6', 'WELL-7', 'WELL-9'] if (-1500, w) in ce.index)
      + f" degC; feed zones below -1000 mRL (excluding WELL-8) are at {deepfz.T_degC.min():.0f}-{deepfz.T_degC.max():.0f} degC | not tied to a mapped feature | medium |")
    A(f"| Recharge / downflow | NE (WELL-4) and SW (WELL-8) | W4 at most {S.loc['WELL-4', 'Tmax_degC']:.0f} degC; W8 at most "
      f"{S.loc['WELL-8', 'Tmax_degC']:.0f} degC, with an isothermal section; both have positive FZ pressure residuals | "
      f"W4 is {Gd.loc['WELL-4', 'dist_Fault3_m']:.0f} m from Fault 3 (cold water down or along F3 to the NE); W8 is "
      f"{Gd.loc['WELL-8', 'dist_lake_shore_m']:.0f} m from the lake shore, the nearest well to the lake | medium |")
    A(f"| Conductive margin / discharge | E (WELL-1; WELL-5 is distal outflow) | W1 conductive, artesian ({W['WELL-1']['P'][0]:.1f} bar at the "
      f"wellhead), lowest FZ pressure residuals | W1 is {Gd.loc['WELL-1', 'dist_river_m']:.0f} m and W5 "
      f"{Gd.loc['WELL-5', 'dist_river_m']:.0f} m from the river: likely the discharge area | medium |")
    A(f"| Cap | surface to about -100/-250 mRL over the hot zone | gradients {cap_hot.shallow_gradient_degC_per_km.min():.0f}-"
      f"{cap_hot.shallow_gradient_degC_per_km.max():.0f} degC/km behind casing; 100 degC isotherm at "
      f"{cap_hot.z_100C_mRL.min():.0f} to {cap_hot.z_100C_mRL.max():.0f} mRL | thickness from wellhead to the 240 degC isotherm: "
      + ", ".join(f"{K(w)} {v:.0f} m" for w, v in capth.items()) + " | medium |")
    A(f"| Deep pressure connection | whole field, -260 to -1910 mRL | FZ residuals within {res.abs().max():.1f} bar of one curve | "
      f"no pressure offset > 2.2 bar across Fault 1 or Fault 2 | medium |")
    A("")
    # ------------------------------------------------------------------ 6
    A("## 6. Recommended calibration set-up")
    A("")
    A("**Targets** (NS_calibration_targets.csv, feedzones.csv, T_P_at_layer_centres.csv).")
    A(f"- Temperature: use all {n_pts} survey points with the `weight` column. If the model compares at block "
      f"centres, use T_P_at_layer_centres.csv instead ({len(LC)} rows, interpolated, never extrapolated).")
    A(f"- Pressure: use only the {len(F)} feed-zone pressures (`weight_P` > 0 in NS_calibration_targets.csv, or "
      f"feedzones.csv `target_weight`). Give all other pressures zero weight (see conclusion 6).")
    A("- Enthalpy: no natural-state enthalpy data exist. Production enthalpy is out of scope.")
    A("")
    A("**Weighting scheme.** Weight = location weight x natural-state confidence multiplier.")
    A("")
    A("| Location class | Weight | Reason | Points |")
    A("|---|---|---|---|")
    reasons = {"feed zone": "formation T and P are measured directly",
               "open hole": "close to formation, but internal flow is possible",
               "open hole, possible internal flow": "isothermal open hole between feed zones (WELL-2)",
               "behind casing": "conductive formation T; no formation pressure",
               "near wellhead / top layer": "within 30 m of the wellhead or in the top model block; atmospheric/air effects",
               "above model surface": "no model block"}
    vc = TGT.location_class.value_counts()
    for k, v in WCLASS.items():
        A(f"| {k} | {v} | {reasons[k]} | {int(vc.get(k, 0))} |")
    A("")
    A("Confidence multipliers: high 1.0 (WELL-2, WELL-3), medium 0.8 (all others), from section 2. For feed-zone "
      "pressures (feedzones.csv `target_weight`), the quality tag (HIGH 1.0, MEDIUM 0.8, LOW or Synthetic 0.5) is "
      "multiplied by the same confidence.")
    A("")
    A("**Where heat and mass enter and leave** (interpretation from sections 3-5):")
    A(f"- **In:** a hot liquid upflow at the base, under the WELL-2/WELL-10 columns "
      f"({Gd.loc['WELL-2', 'column']}, {Gd.loc['WELL-10', 'column']}), next to the F1 x F3 intersection. Enthalpy at least "
      f"{h_liquid(w10.Tmax_degC, fz10.P_bar):.0f} kJ/kg. The mass rate is not constrained by these data; tune it to match "
      f"the outflow temperatures. Cold water enters laterally or from above at the NE (WELL-4 / Fault 3) and SW (WELL-8) "
      f"margins.")
    A(f"- **Out:** conduction through the cap over the hot zone (about {cap_hot.conductive_flux_cap_Wm2_assumed_k2.min():.1f}-"
      f"{cap_hot.conductive_flux_cap_Wm2_assumed_k2.max():.1f} W/m2, assuming k = {K_ROCK} W/m/K), and lateral outflow to "
      f"the east and south-east, towards the river near WELL-5 and WELL-1. WELL-1's artesian head suggests the river "
      f"reach is a discharge boundary.")
    A("")
    A("**Permeability structure the data imply** (interpretation; these are not values from the reference model):")
    A(f"- A vertical high-permeability conduit at the F1 x F3 intersection (WELL-10, WELL-2), from the base up to about -250 mRL.")
    A(f"- A laterally continuous permeable layer at about -250 to -750 mRL, which holds the HIGH feed zones at "
      + ", ".join(f"{K(r.well)} {r.mid_mRL:.0f}" for _, r in hi_fz[hi_fz.mid_mRL > -800].iterrows())
      + " mRL.")
    A("- Low vertical permeability between the outflow layer and the cooler water below. Without it the inversions "
      "(up to 98 degC in WELL-6) would not survive.")
    A("- A second permeable horizon at about -1100 to -1550 mRL that carries cooler water (HIGH feed zones at "
      + ", ".join(f"{K(r.well)} {r.mid_mRL:.0f}" for _, r in hi_fz[hi_fz.mid_mRL <= -800].iterrows())
      + " mRL).")
    A("- A low-permeability cap from the surface to about -100/-250 mRL ("
      + ", ".join(f"240 degC reached {v:.0f} m below the wellhead at {K(w)}" for w, v in capth.items()) + ").")
    A("- Fault 3 to the NE of the intersection is a path for cold water, not hot (WELL-4).")
    A("")
    A("**Boundary requirements.**")
    A("- Top: atmospheric, with rain recharge. The model surface must be lowered at WELL-1 and WELL-8 (section 1, item 8).")
    A("- Base: localised hot mass input at F1 x F3, plus background conductive heat elsewhere. A uniform basal heat flux "
      "cannot make a 98 degC inversion or a near-boiling well 600 m deep.")
    A("- Sides: recharge or open boundaries at the NE and SW, and a discharge or open boundary to the east (river).")
    A(f"- **[reference model]** Default heat flux {REF.get('heat_flux_default_Wm2')} W/m2; rain {REF.get('rain_default_mm')} "
      f"mm/yr with infiltration {REF.get('infiltration')}; lake wet columns {', '.join(REF.get('lake_columns', []))}. "
      f"The deep mass upflow is missing.")
    A("")
    # ------------------------------------------------------------------ 7
    A("## 7. Data quality issues and open questions (ranked by impact on calibration)")
    A("")
    A(f"1. **Shallow super-hydrostatic pressures** (section 4). These decide whether shallow P can be a target and how "
      f"permeable the cap must be. Ask the provider whether these are measured wellbore pressures or derived formation "
      f"pressures. Until then, only feed-zone P are used.")
    A(f"2. **Survey date 29/07/2015 for all wells.** The pressures say pre-production. If the surveys really are from "
      f"mid-2015, the producer-area pressures would be {ddrop.min():.0f}-"
      f"{ddrop.max():.0f} bar too low for the natural state. Confirm the date.")
    A(f"3. **Model surface above wellheads** by up to {-Gd.wellhead_minus_surface_m.min():.0f} m (WELL-1). This "
      f"biases the top-layer temperature and pressure comparisons. Fix SURFA.")
    A("4. **WELL-2 isothermal section** (-1100 to -1900 mRL), a possible wellbore internal flow. The deep temperature "
      "there may not be formation temperature; it is down-weighted.")
    A(f"5. **WELL-10 depth and tag.** WELL-10 ends at {w10.TD_mRL:.0f} mRL and its feed-zone quality is "
      f"'{fz10.quality}'. The depth extent of the upflow is constrained only by WELL-2's deep reheating.")
    A(f"6. **WELL-9 offset of {pcw.loc['WELL-9', 'diff_survey_minus_first_bar']:+.1f} bar** at the gauge. Check the "
      f"gauge depth reference (measured depth or elevation) on this deviated well.")
    A("7. **WELL-8 isothermal 67-68 degC** around its casing shoe. Is it injection or natural downflow? Medium confidence "
      "is kept.")
    A(f"8. **WELL-4 design** lists the liner to -3251 mRL, but the track ends at {W['WELL-4']['td']:.0f} mRL. The survey "
      f"starts {Gd.loc['WELL-4', 'survey_top_below_wellhead_m']:.0f} m below the wellhead.")
    A(f"9. **No data near the lake.** The nearest well is WELL-8, {Gd.dist_lake_shore_m.min():.0f} m from the shore. "
      f"The western boundary is not constrained.")
    A(f"10. **Rock conductivity is assumed** ({K_ROCK} W/m/K) in every heat-flux number. Scale the fluxes linearly if "
      f"another value is used.")
    A("")
    # ------------------------------------------------------------------ 8
    A("## 8. Output file index")
    A("")
    idx = [("Scripts/natural_state_analysis.py", "regenerates every output below"),
           ("NS_data_interpretation.md", "this report"),
           ("NS_calibration_targets.csv", "one row per survey point: block, Tsat, boiling margin, class, weight, weight_P, flags"),
           ("feedzones.csv", "feed-zone targets with blocks, T, P, boiling margin, enthalpy, residuals"),
           ("T_P_at_layer_centres.csv", "T and P at each model layer centre along each well (interpolated only)"),
           ("T_P_at_common_elevations.csv", "T and P at fixed elevations, used for the maps and tables"),
           ("well_summary.csv", "per-well characteristics and regime"),
           ("pressure_consistency_check.csv", "survey P vs first gauge P"),
           ("hydrostatic_check.csv", "column density tests per window"),
           ("hydrostatic_band_summary.csv", "hydrostatic ratios by well and depth band"),
           ("well_geometry_check.csv", "wellhead vs model surface, deviation, distances to faults, river and lake"),
           ("data_integrity_check.csv", "JSON vs dat, track, design and feed-zone checks"),
           ("key_numbers.json", "scalar results (field ratios, reference-model generator summary)"),
           ("run_console_log.txt", "console output of the final run"),
           ("process_log.md", "prompt, rules applied, steps and departures"),
           ("figures/", "well map, profiles, P-T plot, maps by elevation, cross-sections, hydrostatic, "
                        "natural-state test, feed-zone pressure, conceptual sketch")]
    A("| File (in Output/ unless stated) | Contents |")
    A("|---|---|")
    for a, b in idx:
        A(f"| {a} | {b} |")
    A("")
    A("Figures: " + ", ".join(f"`{os.path.basename(f)}`" for f in figs))
    txt = "\n".join(L) + "\n"
    # degC -> °C in prose (CSV column names such as T_degC are left unchanged)
    import re
    txt = re.sub(r"(?<![_A-Za-z])degC", "°C", txt)
    open(os.path.join(OUT, "NS_data_interpretation.md"), "w", encoding="utf-8").write(txt)
    print(f"Report written: Output/NS_data_interpretation.md ({len(L)} lines)")


# ============================================================================
# 15. Save key scalars for the report
# ============================================================================
KEY = {"field_ratio_cold": FIELD_RATIO, "field_ratio_all": FIELD_RATIO_ALL, "band_ratio_cold": BAND_RATIO,
       "implied_g": G * FIELD_RATIO, "ref": REF, "fault_x": {k: list(v) for k, v in FAULT_X.items()},
       "fz_fit": list(fzfit), "dd_2020": dd_2020}
json.dump(KEY, open(os.path.join(OUT, "key_numbers.json"), "w"), indent=2, default=float)
print("Key numbers:", json.dumps(KEY, default=float)[:600])

write_report()
print("Done.")
