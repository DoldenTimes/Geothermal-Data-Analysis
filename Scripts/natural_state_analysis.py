#!/usr/bin/env python
"""
Natural-state data interpretation for geothermal model calibration.

Run from the repository root:   python Scripts/natural_state_analysis.py

Reads ONLY from Baseline/Teaching_Model_3D (read-only), writes ONLY to Output/ (and Output/figures/).
Sections of the run:
  0  IAPWS-IF97 implementation + validation (aborts if any error > 0.1 K)
  1  Load data (JSON = data, .dat = cross-check), grid geometry, mapping polylines
  2  Integrity checks (JSON vs dat, survey/feed-zone/track/design consistency)
  3  Natural-state timing tests and per-well verdicts / confidence
  4  Hydrostatic checks, boiling margin
  5  Per-well characteristics, regime classification
  6  Lateral structure tables (common elevations, layer centres), calibration targets
  7  Figures
  8  Report (Output/NS_data_interpretation.md)
"""
import ast
import json
import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.patches import Polygon as MplPolygon

ROOT = Path(__file__).resolve().parent.parent
BASE = ROOT / "Baseline" / "Teaching_Model_3D"
OUT = ROOT / "Output"
FIG = OUT / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)


class Tee:
    """Duplicate console output into Output/run_console_log.txt"""
    def __init__(self, path):
        self.f = open(path, "w", encoding="utf-8")
        self.o = sys.stdout
    def write(self, s):
        self.o.write(s); self.f.write(s)
    def flush(self):
        self.o.flush(); self.f.flush()


sys.stdout = Tee(OUT / "run_console_log.txt")


def hdr(s):
    print("\n" + "=" * 100 + f"\n{s}\n" + "=" * 100)


G = 9.80665          # m/s2 (standard gravity)
P_ATM = 1.01325      # bar
K_ROCK_ASSUMED = 2.0  # W/m/K  ASSUMPTION (typical volcanic/sedimentary rock; not in the data)

# ----------------------------------------------------------------------------------------------
# 0. IAPWS-IF97  (Region 4 saturation, Region 1 liquid)
# ----------------------------------------------------------------------------------------------
_N4 = [0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5,
       -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6,
       -0.23855557567849, 0.65017534844798e3]


def psat_MPa(T_K):
    n = _N4; T = np.asarray(T_K, float)
    th = T + n[8] / (T - n[9])
    A = th ** 2 + n[0] * th + n[1]; B = n[2] * th ** 2 + n[3] * th + n[4]; C = n[5] * th ** 2 + n[6] * th + n[7]
    return (2 * C / (-B + np.sqrt(B ** 2 - 4 * A * C))) ** 4


def tsat_K(P_MPa):
    n = _N4; b = np.asarray(P_MPa, float) ** 0.25
    E = b ** 2 + n[2] * b + n[5]; F = n[0] * b ** 2 + n[3] * b + n[6]; Gg = n[1] * b ** 2 + n[4] * b + n[7]
    D = 2 * Gg / (-F - np.sqrt(F ** 2 - 4 * E * Gg))
    return (n[9] + D - np.sqrt((n[9] + D) ** 2 - 4 * (n[8] + n[9] * D))) / 2


def tsat_C(P_bar):
    return tsat_K(np.asarray(P_bar, float) / 10.0) - 273.15


_R1 = [(0, -2, 0.14632971213167), (0, -1, -0.84548187169114), (0, 0, -0.37563603672040e1), (0, 1, 0.33855169168385e1),
       (0, 2, -0.95791963387872), (0, 3, 0.15772038513228), (0, 4, -0.16616417199501e-1), (0, 5, 0.81214629983568e-3),
       (1, -9, 0.28319080123804e-3), (1, -7, -0.60706301565874e-3), (1, -1, -0.18990068218419e-1),
       (1, 0, -0.32529748770505e-1), (1, 1, -0.21841717175414e-1), (1, 3, -0.52838357969930e-4),
       (2, -3, -0.47184321073267e-3), (2, 0, -0.30001780793026e-3), (2, 1, 0.47661393906987e-4),
       (2, 3, -0.44141845330846e-5), (2, 17, -0.72694996297594e-15), (3, -4, -0.31679644845054e-4),
       (3, 0, -0.28270797985312e-5), (3, 6, -0.85205128120103e-9), (4, -5, -0.22425281908000e-5),
       (4, -2, -0.65171222895601e-6), (4, 10, -0.14341729937924e-12), (5, -8, -0.40516996860117e-6),
       (8, -11, -0.12734301741641e-8), (8, -6, -0.17424871230634e-9), (21, -29, -0.68762131295531e-18),
       (23, -31, 0.14478307828521e-19), (29, -38, 0.26335781662795e-22), (30, -39, -0.11947622640071e-22),
       (31, -40, 0.18228094581404e-23), (32, -41, -0.93537087292458e-25)]


def rho_liq(T_K, P_MPa):
    """Region 1 liquid density, kg/m3 (valid 273.15-623.15 K, up to 100 MPa)."""
    pi = np.asarray(P_MPa, float) / 16.53; tau = 1386.0 / np.asarray(T_K, float)
    gp = sum(-n * I * (7.1 - pi) ** (I - 1) * (tau - 1.222) ** J for I, J, n in _R1)
    v = 0.461526e-3 * np.asarray(T_K, float) / np.asarray(P_MPa, float) * pi * gp  # m3/kg
    return 1.0 / v


hdr("0. IAPWS-IF97 validation (script aborts if any T error > 0.1 K)")
checks = [
    ("Tsat(1.01325 bar)", float(tsat_C(1.01325)), 99.97, 0.1, "degC"),
    ("Tsat(100 bar)", float(tsat_C(100.0)), 311.0, 0.1, "degC"),
    ("Tsat(0.1 MPa) IAPWS table 35", float(tsat_K(0.1)), 372.755919, 0.01, "K"),
    ("Tsat(1 MPa)   IAPWS table 35", float(tsat_K(1.0)), 453.035632, 0.01, "K"),
    ("Tsat(10 MPa)  IAPWS table 35", float(tsat_K(10.0)), 584.149488, 0.01, "K"),
    ("Psat(300 K)   IAPWS table 35 [kPa]", float(psat_MPa(300.0)) * 1e3, 3.53658941, 1e-4, "kPa"),
    ("Psat(500 K)   IAPWS table 35 [MPa]", float(psat_MPa(500.0)), 2.63889776, 1e-6, "MPa"),
    ("Psat(600 K)   IAPWS table 35 [MPa]", float(psat_MPa(600.0)), 12.3443146, 1e-5, "MPa"),
    ("rho(3 MPa,300 K) IAPWS table 5", float(rho_liq(300, 3)), 1 / 0.100215168e-2, 0.01, "kg/m3"),
    ("rho(80 MPa,300 K) IAPWS table 5", float(rho_liq(300, 80)), 1 / 0.971180894e-3, 0.01, "kg/m3"),
    ("rho(3 MPa,500 K) IAPWS table 5", float(rho_liq(500, 3)), 1 / 0.120241800e-2, 0.01, "kg/m3"),
]
for lab, c, r, tol, u in checks:
    ok = abs(c - r) <= tol
    print(f"  {lab:40s} computed={c:14.6f} reference={r:14.6f} error={c - r:+.2e} {u:8s} {'OK' if ok else 'FAIL'}")
    if not ok:
        sys.exit(f"IF97 validation failed for {lab}")
print("  All IAPWS-IF97 validation checks passed (T errors << 0.1 K).")


# ----------------------------------------------------------------------------------------------
# 1. Load
# ----------------------------------------------------------------------------------------------
def lit(x):
    """JSON lists are stored as a mixture of real lists and strings like '[-310.0, -360.0, 'HIGH']'."""
    if isinstance(x, str):
        return ast.literal_eval(x)
    return x


hdr("1. Loading data")
groups = json.load(open(BASE / "Data" / "Groups_data.json"))
WELLS = groups["Areas"]["Synthetic"]["All"]
print("Wells (Groups_data.json):", WELLS)

W = {}
for w in WELLS:
    n = int(w.split("-")[1])
    d = json.load(open(BASE / "Data" / "Data_json" / f"{w}_data.json"))
    rec = {"name": w, "n": n, "raw": d}
    trk = np.array([lit(t)[:3] for t in d["Well_Track"]["Data"]], float)
    rec["track"] = trk
    rec["fz"] = [lit(f) for f in d["Feed_Zones"]["Data"]]            # [top, bottom, quality] in mRL
    rec["design"] = [lit(x) for x in d["Well_Design"]["Data"]]       # [top, bottom, dia, slotted, spec]
    Td = d["Downhole_Temperature"]; Pd = d["Downhole_Pressure"]
    assert len(Td) == 1 and len(Pd) == 1, "more than one survey date: not handled"
    (tdate, tv), = Td.items(); (pdate, pv), = Pd.items()
    rec["survey_date_T"], rec["survey_date_P"] = tdate, pdate
    rec["stable_T"], rec["stable_P"] = tv["StableFlag"], pv["StableFlag"]
    zs = sorted(tv["Data"].keys(), key=lambda s: -float(s))
    rec["T_json"] = {float(k): tv["Data"][k] for k in zs}
    rec["P_json"] = {float(k): pv["Data"][k] for k in sorted(pv["Data"].keys(), key=lambda s: -float(s))}
    rec["z"] = np.array([float(k) for k in zs])
    rec["T"] = np.array([tv["Data"][k] for k in zs], float)
    rec["P"] = np.array([pv["Data"][k] for k in zs], float)
    rec["wh"] = trk[0]
    # production / injection
    rec["role"] = "producer" if "Production_Mass_Flow" in d else ("injector" if "Injection_Mass_Flow" in d else "monitor (no flow data)")
    rec["flow"] = d.get("Production_Mass_Flow", d.get("Injection_Mass_Flow"))
    rec["enth"] = d.get("Transient_Enthalpy", d.get("Injection_Enthalpy"))
    rec["tp"] = d.get("Transient_Pressure")
    W[w] = rec
    print(f"  {w}: role={rec['role']:24s} n_survey={len(rec['z']):2d} survey date={tdate} track_pts={len(trk)} "
          f"wellhead z={trk[0][2]:.0f} TD z={trk[-1][2]:.1f}  FZ={[(f[0], f[1], f[2]) for f in rec['fz']]}")

# mapping
mapping = json.load(open(BASE / "Data" / "Mapping" / "Mapping.json"))
faults = {k: np.array(v["Data"], float) for k, v in mapping.items() if k.startswith("Fault")}
rivers = {}
for k in ("Tarawera River", "river"):
    rivers[k] = [np.array(seg, float) for seg in mapping[k]["Data"]]
lake = np.array(mapping["Lake Tarawera"]["Data"], float)
for k, v in faults.items():
    print(f"  {k}: {len(v)} pts  from ({v[0][0]:.0f},{v[0][1]:.0f}) to ({v[-1][0]:.0f},{v[-1][1]:.0f})")
print(f"  Lake polygon: {len(lake)} pts, extent x {lake[:,0].min():.0f}-{lake[:,0].max():.0f} y {lake[:,1].min():.0f}-{lake[:,1].max():.0f}")

# grid (geometry only: VERTICES, GRID, LAYERS, SURFA, WELLS)
Lines = [l.rstrip("\r") for l in open(BASE / "Models" / "TE14515" / "gTE14515.dat").read().split("\n")]
iv, ig, ic, il, isf, iw = (Lines.index(k) for k in ("VERTICES", "GRID", "CONNECTIONS", "LAYERS", "SURFA", "WELLS"))
VERT = {l[:3]: (float(l[3:13]), float(l[13:23])) for l in Lines[iv + 1:ig] if l.strip()}
colv = {}
cur = None
for l in Lines[ig + 1:ic]:
    if not l.strip():
        continue
    if len(l.strip()) > 3:
        cur = l[:3]; colv[cur] = []
    else:
        colv[cur].append(l.strip())
COLPOLY = {c: np.array([VERT[v] for v in vs]) for c, vs in colv.items()}
COLCEN = {c: p.mean(axis=0) for c, p in COLPOLY.items()}
SURF = {l[:3]: float(l[3:]) for l in Lines[isf + 1:iw] if l.strip()}
LAY = []
for l in Lines[il + 1:isf]:
    if l.strip():
        a = l.split(); LAY.append((int(a[0]), float(a[1]), float(a[2])))   # name, bottom, centre
LAY = pd.DataFrame(LAY, columns=["layer", "bottom", "centre"])
LAY["top"] = LAY["bottom"].shift(1)
print(f"  Grid: {len(VERT)} vertices, {len(COLPOLY)} columns, {len(LAY)} layers (layer 0 = atmosphere), "
      f"surface elevation range {min(SURF.values()):.0f}..{max(SURF.values()):.0f} mRL")

wellsjson = json.load(open(BASE / "Models" / "TE14515" / "gTE14515_wells.json"))


def inpoly(x, y, P):
    c = False; n = len(P)
    for i in range(n):
        x1, y1 = P[i]; x2, y2 = P[(i + 1) % n]
        if (y1 > y) != (y2 > y) and x < (x2 - x1) * (y - y1) / (y2 - y1) + x1:
            c = not c
    return c


_poly_cache = {}


def find_column(x, y):
    for c, P in COLPOLY.items():
        mn = P.min(axis=0); mx = P.max(axis=0)
        if mn[0] <= x <= mx[0] and mn[1] <= y <= mx[1] and inpoly(x, y, P):
            return c
    return None


def find_layer(z):
    for _, r in LAY.iloc[1:].iterrows():
        if r.bottom < z <= r.top:
            return int(r.layer)
    return None


def track_xy(w, z):
    t = W[w]["track"]
    if len(t) == 1:
        return t[0][:2]
    zz = t[::-1, 2]
    return np.array([np.interp(z, zz, t[::-1, 0]), np.interp(z, zz, t[::-1, 1])])


def block_at(w, z):
    x, y = track_xy(w, z)
    c = find_column(x, y); l = find_layer(z)
    return c, l, (f"{c}{l:02d}" if c and l is not None else None)


def seg_dist(p, a, b):
    ab = b - a; t = np.clip(np.dot(p - a, ab) / np.dot(ab, ab), 0, 1)
    return np.hypot(*(p - (a + t * ab)))


def dist_polyline(p, pts):
    return min(seg_dist(p, pts[i], pts[i + 1]) for i in range(len(pts) - 1))


def dist_polys(p, polys):
    return min(dist_polyline(p, s) for s in polys)


def interp_profile(w, z, arr="T"):
    """Linear interpolation within survey range only (NaN outside)."""
    zz = W[w]["z"][::-1]; vv = W[w][arr][::-1]
    z = np.asarray(z, float)
    out = np.interp(z, zz, vv)
    out = np.where((z < zz.min() - 1e-9) | (z > zz.max() + 1e-9), np.nan, out)
    return out


# ----------------------------------------------------------------------------------------------
# 2. Integrity checks
# ----------------------------------------------------------------------------------------------
hdr("2. Data integrity: JSON vs .dat, survey/feed zones/track/design")
integ = []


def add(check, well, detail, result, value=np.nan):
    integ.append(dict(check=check, well=well, detail=detail, result=result, value=value))


def datfile(w, kind):
    n = W[w]["n"]
    pre = "Wel10" if n == 10 else f"Well{n}"
    return BASE / "Data" / "Data_dat" / f"{pre}_downhole_{kind}.dat"


for w, r in W.items():
    for kind, key in (("temperature", "T"), ("pressure", "P")):
        a = np.loadtxt(datfile(w, kind)); a = a.reshape(-1, 2)
        dat_z = a[:, 1]; dat_v = a[:, 0]
        js = r[key + "_json"]
        same_n = len(a) == len(js)
        maxdiff = 0.0; missing = 0
        for z, v in zip(dat_z, dat_v):
            if z in js:
                maxdiff = max(maxdiff, abs(js[z] - v))
            else:
                missing += 1
        ok = same_n and missing == 0 and maxdiff < 1e-6
        add(f"JSON vs dat {kind}", w, f"n_json={len(js)} n_dat={len(a)} missing_in_json={missing} max|diff|={maxdiff:.2e}",
            "PASS" if ok else "FLAG", maxdiff)
    # transient pressure
    if r["tp"]:
        for zk, tpv in r["tp"].items():
            pre = "Wel10" if r["n"] == 10 else f"Well{r['n']}"
            a = np.loadtxt(BASE / "Data" / "Data_dat" / f"{pre}_{float(zk)}_transient_pressure.dat").reshape(-1, 2)
            js = list(tpv["Data"].values())
            ok = len(a) == len(js) and np.allclose(a[:, 1], js)
            add("JSON vs dat transient pressure", w, f"gauge {zk} mRL; n_json={len(js)} n_dat={len(a)}; dat time column = years since 2014 "
                f"({a[0,0]:.0f}..{a[-1,0]:.0f})", "PASS" if ok else "FLAG")
    # survey elevations within track
    tz_top, tz_bot = r["track"][0][2], r["track"][-1][2]
    out = [z for z in r["z"] if z > tz_top + 1e-6 or z < tz_bot - 1e-6]
    add("Survey points inside track z-range", w, f"track {tz_top:.1f}..{tz_bot:.1f} mRL; survey {r['z'][0]:.1f}..{r['z'][-1]:.1f}; outside={len(out)}",
        "PASS" if not out else "FLAG", len(out))
    # T & P survey duplicates / monotonic
    add("Survey elevations strictly descending", w, "", "PASS" if np.all(np.diff(r["z"]) < 0) else "FLAG")
    # stable flags
    add("StableFlag", w, f"T={r['stable_T']} P={r['stable_P']} (date {r['survey_date_T']})", "PASS" if r["stable_T"] and r["stable_P"] else "FLAG")
    # design continuity
    des = r["design"]
    for i, d in enumerate(des):
        add("Design interval inside track z-range", w, f"{d[4]}: {d[0]}..{d[1]} mRL, slotted={d[3]}; track bottom {tz_bot:.1f}",
            "PASS" if d[1] >= tz_bot - 1e-6 else "FLAG (design extends below track end)", d[1] - tz_bot)
    cas = [d for d in des if not d[3]]
    liners = [d for d in des if d[3]]
    shoe = min(c[1] for c in cas) if cas else np.nan
    r["shoe"] = shoe
    if liners:
        lap = max(l[0] for l in liners) - shoe   # >=0: liner top above casing shoe -> liner laps inside casing
        add("Casing/liner overlap or gap", w, f"casing shoe {shoe:.0f}, liner top {max(l[0] for l in liners):.0f}", "PASS (liner laps inside casing)" if lap >= 0 else "FLAG (gap)", lap)
    # feed zones
    for k, f in enumerate(r["fz"], 1):
        top, bot, q = f
        inside = tz_bot - 1e-6 <= bot and top <= tz_top + 1e-6
        below_shoe = top <= shoe + 1e-6
        in_liner = any(l[3] and l[1] <= bot and top <= l[0] for l in liners)
        status = "slotted liner" if in_liner else ("open hole below shoe" if below_shoe else "behind blank casing")
        add("Feed zone within track and open/slotted interval", w, f"FZ{k} {top:.0f}..{bot:.0f} mRL [{q}] -> {status}",
            "PASS" if inside and below_shoe else "FLAG")
        # survey point at FZ?
        has_pt = any(bot - 1e-6 <= z <= top + 1e-6 for z in r["z"])
        add("Survey point inside feed zone", w, f"FZ{k} {top:.0f}..{bot:.0f}", "PASS" if has_pt else "FLAG (interpolated only)")
    add("Feed-zone units", w, "Feed_Zones.Units say 'm' but values are elevations (mRL, negative below sea level) - treated as mRL", "NOTE")
    # wells.json cross-check: blocks
    wj = wellsjson[w]
    blocknames = {b["block"] for b in wj["Blocks"]}
    miss = 0
    for z in r["z"]:
        c, l, b = block_at(w, z)
        if b not in blocknames:
            miss += 1
    add("Survey point blocks (own grid lookup) in gTE14515_wells.json", w, f"{len(r['z'])} points, not found={miss}", "PASS" if miss == 0 else "FLAG", miss)

integ = pd.DataFrame(integ)
integ.to_csv(OUT / "data_integrity_check.csv", index=False)
print(integ.groupby(["check", "result"]).size().to_string())
print("FLAGS:\n", integ[integ.result.str.startswith("FLAG")][["check", "well", "detail"]].to_string())

# ----------------------------------------------------------------------------------------------
# 3. Natural-state timing
# ----------------------------------------------------------------------------------------------
hdr("3. Natural-state timing tests")


def parse_date(s):
    return pd.to_datetime(s, format="%d/%m/%Y")


rows = []
for w, r in W.items():
    row = dict(well=w, role=r["role"], survey_date=r["survey_date_P"])
    # first flow date
    if r["flow"]:
        fl = {parse_date(k): v for k, v in r["flow"]["Data"].items()}
        d0 = min(fl); row["first_flow_date"] = d0.strftime("%d/%m/%Y"); row["first_flow_kg_s"] = fl[d0]
        prior = [k for k in fl if k <= parse_date(r["survey_date_P"])]
        row["flow_at_survey_date_kg_s"] = fl[max(prior)] if prior else 0.0   # 0 = no flow recorded before survey
    if r["role"] == "injector":
        h = list(r["enth"]["Data"].values())[0]
        row["injectate_enthalpy_kJ_kg"] = h
        row["injectate_T_approx_degC"] = round(h / 4.19, 1)   # ASSUMPTION cp=4.19 kJ/kg/K (liquid water at ~80 C)
    if r["tp"]:
        (zk, tpv), = r["tp"].items()
        dates = [parse_date(k) for k in tpv["Data"]]; vals = list(tpv["Data"].values())
        zg = float(zk)
        row["gauge_elev_mRL"] = zg
        row["first_transient_date"] = dates[0].strftime("%d/%m/%Y"); row["first_transient_P_bar"] = vals[0]
        row["survey_P_at_gauge_bar"] = float(np.round(interp_profile(w, zg, "P"), 2))
        row["delta_survey_minus_first_bar"] = round(row["survey_P_at_gauge_bar"] - vals[0], 2)
        row["last_transient_P_bar"] = vals[-1]
        row["first_transient_after_first_flow"] = bool(dates[0] > parse_date(row["first_flow_date"]))
        row["survey_before_first_transient_yrs"] = round((dates[0] - parse_date(r["survey_date_P"])).days / 365.25, 2)
        row["drawdown_first_to_last_bar"] = round(vals[0] - vals[-1], 2)
        # drawdown of this well by 01/01/2015 (first anniversary) if series long enough
        for d, v in zip(dates, vals):
            if d.year == 2015:
                row["P_at_2015-01-01_bar"] = v
            if d.year == 2020:
                row["P_at_2020-01-01_bar"] = v
    rows.append(row)
pc = pd.DataFrame(rows)

# neighbour-referenced consistency (verified wells only)
# step 1: verified wells = tested wells with valid test and |delta|<=0.5 bar
verified = [r.well for r in pc.itertuples() if hasattr(r, "delta_survey_minus_first_bar") and
            pd.notna(r.delta_survey_minus_first_bar) and r.first_transient_after_first_flow is False and abs(r.delta_survey_minus_first_bar) <= 0.5]
print("Wells verified as natural state by gauge test (valid test, |delta| <= 0.5 bar):", verified)
# Neighbour comparison at FEED-ZONE elevations (formation pressure), NOT over the whole column: column pressures at depth also depend
# on each well's own fluid density, so only feed-zone / pivot pressures are comparable between wells.
def dominant_fz_idx(w):
    rank = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}
    best = None
    for k, f in enumerate(W[w]["fz"]):
        mid = 0.5 * (f[0] + f[1])
        key = (rank.get(f[2], 2), float(interp_profile(w, mid, "T")))
        if best is None or key > best[0]:
            best = (key, k)
    return best[1]


for w in W:
    W[w]["dom"] = dominant_fz_idx(w)
FZ_OFFSETS = {}
for w in W:
    offs = []
    for k, f in enumerate(W[w]["fz"]):
        mid = 0.5 * (f[0] + f[1])
        ref = np.nanmean([interp_profile(v, mid, "P") for v in verified if v != w])
        offs.append((k, mid, float(interp_profile(w, mid, "P")), float(ref)))
    FZ_OFFSETS[w] = offs
    dom = W[w]["dom"]
    pc.loc[pc.well == w, "dP_at_feedzones_vs_verified_bar"] = "; ".join(f"FZ{k+1}@{m:.0f}:{p - r:+.1f}{'*' if k == dom else ''}" for k, m, p, r in offs)
    pc.loc[pc.well == w, "dP_pivot_fz_vs_verified_bar"] = offs[dom][2] - offs[dom][3]
    pc.loc[pc.well == w, "max_abs_dP_fz_vs_verified_bar"] = max(abs(p - r) for k, m, p, r in offs)
pc["mean_dP_vs_verified_bar"] = pc["dP_pivot_fz_vs_verified_bar"]
pc["max_abs_dP_vs_verified_bar"] = pc["max_abs_dP_fz_vs_verified_bar"]

verdict, conf_label, reason = {}, {}, {}
CONF_NUM = {("natural state", "high"): 1.0, ("probably natural state", "medium"): 0.75,
            ("probably natural state", "low"): 0.6, ("possibly disturbed", "low"): 0.4,
            ("possibly disturbed", "medium"): 0.4}
for r in pc.itertuples():
    w = r.well
    has_tp = pd.notna(getattr(r, "first_transient_P_bar", np.nan))
    mdp = r.mean_dP_vs_verified_bar if pd.notna(r.mean_dP_vs_verified_bar) else np.nan
    if has_tp and not r.first_transient_after_first_flow and abs(r.delta_survey_minus_first_bar) <= 0.5:
        v, c = "natural state", "high"
        why = (f"survey P at gauge ({r.gauge_elev_mRL:.0f} mRL) = {r.survey_P_at_gauge_bar:.1f} bar vs first transient value {r.first_transient_P_bar:.1f} bar "
               f"on {r.first_transient_date} (day production began): delta {r.delta_survey_minus_first_bar:+.2f} bar")
    elif has_tp and not r.first_transient_after_first_flow and r.delta_survey_minus_first_bar > 0.5:
        v, c = "probably natural state", "medium"
        why = (f"survey P at gauge ({r.gauge_elev_mRL:.0f} mRL) = {r.survey_P_at_gauge_bar:.1f} bar is {r.delta_survey_minus_first_bar:+.2f} bar ABOVE first transient "
               f"{r.first_transient_P_bar:.1f} bar ({r.first_transient_date}); sign is opposite to drawdown, so not production-disturbed, but the offset exceeds gauge rounding (0.1 bar)")
    elif has_tp and not r.first_transient_after_first_flow and r.delta_survey_minus_first_bar < -0.5:
        v, c = "possibly disturbed", "low"
        why = f"survey P is {-r.delta_survey_minus_first_bar:.1f} bar BELOW first transient value: consistent with drawdown before survey"
    elif has_tp and r.first_transient_after_first_flow:
        v, c = "probably natural state", "low"
        why = (f"transient series starts {r.first_transient_date}, AFTER first flow ({r.first_flow_date}) and {r.survey_before_first_transient_yrs:.1f} yr after the survey: "
               f"test is not valid; survey P at gauge {r.survey_P_at_gauge_bar:.1f} bar vs first transient {r.first_transient_P_bar:.1f} bar ({r.delta_survey_minus_first_bar:+.1f} bar) is "
               f"the sign expected for later drawdown; the survey P at its feed zone is {mdp:+.1f} bar from the verified wells (W2, W3) at the same elevation")
    elif r.role == "injector":
        inj_t = r.injectate_T_approx_degC
        t_fz = W[w]["T"][0] * 0 + float(interp_profile(w, 0.5 * sum(W[w]["fz"][W[w]["dom"]][:2]), "T"))
        hotter = t_fz > inj_t
        pos = pd.notna(mdp) and mdp > 2.0
        if pos or hotter:
            v, c = "possibly disturbed", "low"
        else:
            v, c = "probably natural state", "low"
        why = (f"injector (injection from {r.first_flow_date}, {r.first_flow_kg_s:.0f} kg/s at ~{inj_t:.0f} degC, i.e. {r.first_flow_date} to survey ~1.6 yr); no pressure gauge. "
               f"Survey T at pivot feed zone {t_fz:.0f} degC is " + ("ABOVE the injectate temperature, so cooling by injection would be invisible as a warm signature and cannot be excluded"
               if hotter else "BELOW the injectate temperature, so the profile shows no warm-injection signature") +
               f"; P vs verified wells (W2, W3) at the pivot feed zone: {mdp:+.1f} bar" + ("; positive offset >2 bar is the sign of injection pressurisation" if pos else ""))
    else:
        if pd.notna(mdp) and abs(mdp) <= 3.0 and r.max_abs_dP_vs_verified_bar <= 4.0:
            v, c = "probably natural state", "medium"
        else:
            v, c = "probably natural state", "low"
        why = (f"no flow data and no transient gauge: cannot be tested directly; P offset vs verified wells (W2, W3) at the pivot feed zone: {mdp:+.1f} bar "
               f"(max |offset| over its feed zones {r.max_abs_dP_vs_verified_bar:.1f} bar)" + ("" if abs(mdp) <= 3 else "; offset is large and unexplained"))
    verdict[w], conf_label[w], reason[w] = v, c, why
pc["verdict"] = pc.well.map(verdict); pc["confidence_level"] = pc.well.map(conf_label)
pc["ns_confidence"] = [CONF_NUM[(verdict[w], conf_label[w])] for w in pc.well]
pc["reason"] = pc.well.map(reason)
pc.to_csv(OUT / "pressure_consistency_check.csv", index=False)
for r in pc.itertuples():
    print(f"  {r.well:8s} {r.role:24s} -> {r.verdict:24s} [{r.confidence_level}] ns_conf={r.ns_confidence}\n      {r.reason}")
NSCONF = dict(zip(pc.well, pc.ns_confidence))

# ----------------------------------------------------------------------------------------------
# 4. Survey point table: boiling margin, grid, location class, weights
# ----------------------------------------------------------------------------------------------
hdr("4. Survey-point table (boiling margin, block, location class)")


def dominant_fz(w):
    """Dominant feed zone: highest quality tag; ties -> the zone with the higher survey T."""
    rank = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}
    best = None
    for k, f in enumerate(W[w]["fz"]):
        mid = 0.5 * (f[0] + f[1])
        key = (rank.get(f[2], 2), float(interp_profile(w, mid, "T")))
        if best is None or key > best[0]:
            best = (key, k)
    return best[1]


NEAR_WH = 50.0  # m below wellhead counted as 'near the wellhead / top layer'


def isothermal_sections(z, T, min_len=300.0, tol=3.0):
    """Maximal depth intervals (>= min_len) in which T varies by <= tol degC. Returns list of (ztop, zbot)."""
    out = []; n = len(z); i = 0
    while i < n - 1:
        j = i + 1
        while j < n and (T[i:j + 1].max() - T[i:j + 1].min()) <= tol:
            j += 1
        j -= 1
        if j > i and z[i] - z[j] >= min_len:
            out.append((z[i], z[j])); i = j
        else:
            i += 1
    return out


for w in W:
    W[w]["isoth"] = isothermal_sections(W[w]["z"], W[w]["T"])

pts = []
for w, r in W.items():
    wh_z = r["wh"][2]
    for z, T, P in zip(r["z"], r["T"], r["P"]):
        x, y = track_xy(w, z)
        c, l, b = block_at(w, z)
        ts = float(tsat_C(P))
        infz = any(f[1] - 1e-6 <= z <= f[0] + 1e-6 for f in r["fz"])
        casing = [d for d in r["design"] if not d[3]]
        behind = any(d[1] < z <= d[0] + 1e-6 for d in casing)
        above = (c is not None) and (z > SURF[c] + 1e-6)
        iso = any(a + 1e-6 >= z >= b_ - 1e-6 for a, b_ in r["isoth"])
        dom_fz = r["fz"][r["dom"]]
        is_pivot = dom_fz[1] - 1e-6 <= z <= dom_fz[0] + 1e-6
        # location class and weight (T)
        if above:
            cls, wt, why = "above model surface", 0.0, "above model surface"
        elif infz:
            cls, wt, why = "feed zone", 1.0, "feed zone (formation connection)"
        elif z >= wh_z - NEAR_WH:
            cls, wt, why = "near wellhead", 0.3, f"within {NEAR_WH:.0f} m of wellhead / top layer"
        elif behind:
            cls, wt, why = "behind casing", 0.5, "behind blank casing (conductive formation T, not formation P)"
        else:
            if iso:
                cls, wt, why = "open hole (isothermal section)", 0.5, "open hole/slotted liner inside an isothermal section >=300 m: possible internal wellbore flow; downgraded 0.8 -> 0.5"
            else:
                cls, wt, why = "open hole", 0.8, "open hole / slotted liner"
        conf = NSCONF[w]
        wP = wt if (infz or above) else 0.25 * wt
        if infz:
            prole = "pivot feed zone (formation P)" if is_pivot else "feed zone (formation P)"
        else:
            prole = "wellbore column (not formation P)"
        pts.append(dict(well=w, x=round(x, 1), y=round(y, 1), elevation_mRL=z, T_degC=T, P_bar=P, Tsat_degC=round(ts, 2),
                        boiling_margin_degC=round(ts - T, 2), column=c, layer=l, block=b, in_feedzone=infz, behind_casing=behind,
                        above_model_surface=above, ns_confidence=conf, weight=round(wt * conf, 3),
                        weight_reason=f"{why}; base {wt:g} x ns_confidence {conf:g}",
                        weight_P=round(wP * conf, 3), P_role=prole, location_class=cls, in_isothermal_section=iso))
PT = pd.DataFrame(pts)
cols_req = ["well", "x", "y", "elevation_mRL", "T_degC", "P_bar", "Tsat_degC", "boiling_margin_degC", "column", "layer", "block",
            "in_feedzone", "behind_casing", "above_model_surface", "ns_confidence", "weight", "weight_reason"]
PT[cols_req + ["weight_P", "P_role", "location_class"]].to_csv(OUT / "NS_calibration_targets.csv", index=False)
print(PT.groupby("location_class").size().to_string())
print("Points with boiling margin < 5 degC:\n", PT[PT.boiling_margin_degC < 5][["well", "elevation_mRL", "T_degC", "P_bar", "Tsat_degC", "boiling_margin_degC"]].to_string())

# ----------------------------------------------------------------------------------------------
# 5. Hydrostatic check
# ----------------------------------------------------------------------------------------------
hdr("5. Hydrostatic column-density tests (single-phase intervals >= 300 m only)")
hs = []
MIN_LEN, TARGET_LEN, MARGIN_MIN, COLD_T = 300.0, 400.0, 10.0, 100.0
for w, r in W.items():
    z, T, P = r["z"], r["T"], r["P"]
    ts = tsat_C(P); marg = ts - T
    rho = rho_liq(T + 273.15, P / 10.0)
    k = 0; n = len(z)
    while k < n - 1:
        m = k + 1
        while m < n - 1 and z[k] - z[m] < TARGET_LEN:
            m += 1
        if z[m] - z[-1] > -MIN_LEN and m < n - 1 and (z[m] - z[-1]) < MIN_LEN:   # remainder too short: absorb
            m = n - 1
        sl = slice(k, m + 1)
        L = z[k] - z[m]
        if L >= MIN_LEN:
            dz = -np.diff(z[sl])
            exp = float(np.sum(0.5 * (rho[sl][1:] + rho[sl][:-1]) * G * dz) / 1e5)
            obs = float(P[m] - P[k])
            ok = bool(marg[sl].min() >= MARGIN_MIN)
            rm = float(rho[sl].mean())
            hs.append(dict(well=w, z_top_mRL=z[k], z_bottom_mRL=z[m], length_m=L, n_points=m - k + 1, T_mean_degC=round(T[sl].mean(), 1),
                           T_range_degC=f"{T[sl].min():.0f}-{T[sl].max():.0f}", min_boiling_margin_degC=round(marg[sl].min(), 1),
                           single_phase_ok=ok, dP_obs_bar=round(obs, 2), dP_expected_IF97_bar=round(exp, 2),
                           excess_pct=round((obs / exp - 1) * 100, 2), excess_kg_m3=round((obs / exp - 1) * rm, 1),
                           apparent_density_kg_m3=round(obs / exp * rm, 1), model_density_kg_m3=round(rm, 1),
                           rounding_uncertainty_pct=round(0.1 / obs * 100 * 1.41, 2),   # +-0.1 bar on two readings
                           cold=bool(T[sl].max() < COLD_T),
                           max_liquid_density_kg_m3=round(float(rho_liq(277.15, P[sl].mean() / 10.0)), 1)))
        k = m
HS = pd.DataFrame(hs)
HS["exceeds_max_liquid_gradient"] = (HS.apparent_density_kg_m3 * (1 - HS.rounding_uncertainty_pct / 100 * 2)) > HS.max_liquid_density_kg_m3
HS.loc[~HS.single_phase_ok, ["excess_pct", "excess_kg_m3"]] = np.nan   # excluded from statistics (near-boiling)
HS.to_csv(OUT / "hydrostatic_check.csv", index=False)
ok = HS[HS.single_phase_ok]
cold_ex = ok[ok.cold].excess_pct
FIELD_EXCESS_COLD = float(cold_ex.median())
cold_shallow = ok[ok.cold & (ok.z_top_mRL > -250)].excess_pct
cold_deeper = ok[ok.cold & (ok.z_top_mRL <= -250)].excess_pct
cold_shallow = ok[ok.cold & (ok.z_top_mRL > -250)].excess_pct
cold_deeper = ok[ok.cold & (ok.z_top_mRL <= -250)].excess_pct
DEEP_EXCESS = ok[ok.z_top_mRL <= -1000].excess_pct
print(HS.to_string())
print(f"\nCold windows (max T < {COLD_T:.0f} degC): shallow (top > -250 mRL) n={len(cold_shallow)} median {cold_shallow.median():.2f} %, range {cold_shallow.min():.2f}..{cold_shallow.max():.2f}; "
      f"deeper n={len(cold_deeper)} median {cold_deeper.median():.2f} %, range {cold_deeper.min():.2f}..{cold_deeper.max():.2f}")
print(f"Field-wide excess estimate: median of cold windows (max T<{COLD_T:.0f} degC, n={len(cold_ex)}): {FIELD_EXCESS_COLD:.2f} %  "
      f"(range {cold_ex.min():.2f}..{cold_ex.max():.2f})")
print(f"Deep windows (top <= -1000 mRL, n={len(DEEP_EXCESS)}): median {DEEP_EXCESS.median():.2f} %, range {DEEP_EXCESS.min():.2f}..{DEEP_EXCESS.max():.2f}")
print(f"For comparison, using g = 10 m/s2 instead of 9.80665 would give +{(10 / G - 1) * 100:.2f} % on every window.")
HS["excess_minus_field_pct"] = (HS.excess_pct - FIELD_EXCESS_COLD).round(2)
HS.to_csv(OUT / "hydrostatic_check.csv", index=False)

# ----------------------------------------------------------------------------------------------
# 6. Per-well characteristics
# ----------------------------------------------------------------------------------------------
hdr("6. Per-well characteristics")


def crossings(z, T, level):
    """All elevations where profile crosses level (linear interp), top to bottom, with direction."""
    out = []
    for i in range(len(z) - 1):
        a, b = T[i] - level, T[i + 1] - level
        if a == 0:
            out.append((z[i], "at"))
        elif a * b < 0:
            zc = z[i] + (z[i + 1] - z[i]) * (level - T[i]) / (T[i + 1] - T[i])
            out.append((zc, "down" if b > 0 else "up"))   # 'down' = T rising with depth through level
    return out


def regress(zz, vv):
    return np.polyfit(zz, vv, 1)[0]


def water_level(w):
    """Elevation where the wellbore column pressure extrapolates to 1 atm. Uses the first 3 points with P>=2 bar."""
    r = W[w]; m = r["P"] >= 2.0
    zz, pp = r["z"][m][:3], r["P"][m][:3]
    s, b = np.polyfit(zz, pp, 1)         # bar per m (negative slope in z)
    return (P_ATM - b) / s, -s * 100.0   # elevation, gradient in bar/100 m


summ = []
for w, r in W.items():
    z, T, P = r["z"], r["T"], r["P"]
    i = int(np.argmax(T)); tmax, zmax = T[i], z[i]
    ts = tsat_C(P); marg = ts - T; j = int(np.argmin(marg))
    runmax = np.maximum.accumulate(T)
    dd = runmax - T                       # drop below the highest temperature above
    jb = int(np.argmax(dd)); drop = float(dd[jb]); tminb = T[jb]
    i_pk = int(np.where(T[:jb + 1] == runmax[jb])[0][0])   # elevation of the peak that this inversion hangs from
    ipeak_z = z[i_pk]
    inv = drop >= 3.0   # operational threshold for an inversion
    # cap & deep gradients
    msh = z >= z[0] - 300.0
    g_sh = regress(z[msh], T[msh]) * -1000.0
    mdp = z <= z[-1] + 500.0
    g_dp = regress(z[mdp], T[mdp]) * -1000.0
    # isotherms
    cr100 = crossings(z, T, 100.0); cr200 = crossings(z, T, 200.0)
    first = lambda cr: cr[0][0] if cr else np.nan
    last_up = lambda cr: [c[0] for c in cr if c[1] == "up"][-1] if any(c[1] == "up" for c in cr) else np.nan
    cr90 = crossings(z, T, 0.9 * tmax) if tmax >= 150 else []
    wl, wlg = water_level(w)
    dfz = r["fz"][r["dom"]]; dmid = 0.5 * (dfz[0] + dfz[1])
    c0 = find_column(*r["wh"][:2])
    # near-boiling thickness
    nb = z[marg < 5.0]
    # regime
    tdr = T[-1] / tmax
    nb_pts = int(np.sum(marg < 5.0))
    if nb_pts >= 2:
        regime = "two-phase / near-boiling upflow"
    elif tmax >= 200 and tdr >= 0.9:
        regime = "upflow (hot to TD, with shallow outflow peak)" if inv else "upflow"
    elif inv and drop >= 10:
        regime = "outflow over cooler deep water (inversion)"
    elif tmax < 100 and g_sh < 40 and abs(g_dp) < 30:
        regime = "recharge / downflow (cold, sub-conductive gradient)"
    elif inv:
        regime = "conductive-cold with weak inversion (lateral cold flow)"
    else:
        regime = "conductive"
    q_cap = K_ROCK_ASSUMED * g_sh   # mW/m2   (W/m/K * K/km = mW/m2)
    iso_txt = "; ".join(f"{a:.0f} to {b_:.0f} mRL" for a, b_ in r["isoth"]) if r["isoth"] else ""
    summ.append(dict(
        well=w, role=r["role"], x=r["wh"][0], y=r["wh"][1], wellhead_mRL=r["wh"][2], td_mRL=z[-1],
        T_max_degC=tmax, T_max_elev_mRL=zmax, T_at_TD_degC=T[-1],
        inversion=inv, inversion_size_degC=round(drop, 1) if inv else 0.0, inversion_top_mRL=ipeak_z if inv else np.nan,
        inversion_bottom_mRL=z[jb] if inv else np.nan, rebound_below_inversion_degC=round(T[-1] - tminb, 1) if inv else np.nan,
        casing_shoe_mRL=r["shoe"], inversion_vs_shoe=("peak below casing shoe" if (inv and ipeak_z < r["shoe"]) else ("peak behind casing" if inv else "")),
        cap_gradient_K_km=round(g_sh, 1), deep_gradient_K_km=round(g_dp, 1),
        cond_heat_flux_mW_m2_ASSUMED_k2=round(q_cap, 0),
        elev_100C_first_mRL=first(cr100), elev_100C_deepest_up_mRL=last_up(cr100),
        elev_200C_first_mRL=first(cr200), elev_200C_deepest_up_mRL=last_up(cr200),
        elev_T90pct_Tmax_mRL=first(cr90),   # operational 'cap base' (hot wells only)
        
        min_boiling_margin_degC=round(marg[j], 1), min_margin_elev_mRL=z[j], n_points_margin_lt5=nb_pts,
        n_points_margin_lt10=int(np.sum(marg < 10)),
        water_level_est_mRL=round(wl, 0), wellhead_minus_water_level_m=round(r["wh"][2] - wl, 0), shallow_gradient_bar_100m=round(wlg, 2),
        pivot_fz_mid_mRL=dmid, pivot_P_bar=float(interp_profile(w, dmid, "P")), pivot_T_degC=float(interp_profile(w, dmid, "T")),
        pivot_quality=dfz[2], isothermal_sections=iso_txt,
        model_surface_mRL=SURF[c0], regime=regime))
SUM = pd.DataFrame(summ)
SUM.to_csv(OUT / "well_summary.csv", index=False)
print(SUM.T.to_string())

# ----------------------------------------------------------------------------------------------
# Feed zones
# ----------------------------------------------------------------------------------------------
fzr = []
for w, r in W.items():
    wj = wellsjson[w]["Feed_Zones"]
    for k, f in enumerate(r["fz"]):
        top, bot, q = f; mid = 0.5 * (top + bot)
        T = float(interp_profile(w, mid, "T")); P = float(interp_profile(w, mid, "P"))
        ts = float(tsat_C(P))
        c, l, b = block_at(w, mid)
        # blocks spanned
        zs = np.linspace(top, bot, 11)
        bl = []
        for zz in zs:
            bb = block_at(w, zz)[2]
            if bb not in bl:
                bl.append(bb)
        ref = list(wj[f"Feed {k + 1}"].keys())
        mk = (r["z"] >= bot - 1e-6) & (r["z"] <= top + 1e-6)
        fzr.append(dict(well=w, feed=f"Feed {k + 1}", top_mRL=top, bottom_mRL=bot, mid_mRL=mid, thickness_m=top - bot, quality=q,
                        T_degC=round(T, 1), P_bar=round(P, 1), Tsat_degC=round(ts, 1), boiling_margin_degC=round(ts - T, 1),
                        survey_points_in_interval=int(mk.sum()), column=c, layer=l, block_at_mid=b,
                        blocks_spanned=" ".join(bl), blocks_gTE14515_wells_json=" ".join(ref),
                        block_match=bool(b in ref), is_pivot=(k == r["dom"]), ns_confidence=NSCONF[w],
                        weight_T=round(1.0 * NSCONF[w], 3), weight_P=round(1.0 * NSCONF[w], 3),
                        completion=("slotted liner" if any(d[3] and d[1] <= bot and top <= d[0] for d in r["design"]) else "open hole below shoe")))
FZ = pd.DataFrame(fzr)
FZ.to_csv(OUT / "feedzones.csv", index=False)
print(FZ.to_string())

# ----------------------------------------------------------------------------------------------
# Geometry check
# ----------------------------------------------------------------------------------------------
geo = []
for w, r in W.items():
    wh = r["wh"]; td = r["track"][-1]
    c0 = find_column(wh[0], wh[1]); c1 = find_column(td[0], td[1])
    p = wh[:2]
    df = {k: dist_polyline(p, v) for k, v in faults.items()}
    dr = {k: dist_polys(p, v) for k, v in rivers.items()}
    dlake_edge = min(dist_polyline(p, lake[i:i + 2]) for i in range(len(lake) - 1))
    in_lake = inpoly(p[0], p[1], lake)
    kf = min(df, key=df.get); kr = min(dr, key=dr.get)
    # max horizontal deviation anywhere on track
    dev = np.hypot(r["track"][:, 0] - wh[0], r["track"][:, 1] - wh[1]).max()
    fz_blocks_ok = bool(FZ[FZ.well == w].block_match.all())
    geo.append(dict(well=w, x=wh[0], y=wh[1], wellhead_mRL=wh[2], column=c0, model_surface_mRL=SURF[c0],
                    wellhead_minus_model_surface_m=wh[2] - SURF[c0], td_x=td[0], td_y=td[1], td_mRL=td[2], td_column=c1,
                    max_horizontal_deviation_m=round(dev, 1), offset_from_column_centre_m=round(float(np.hypot(*(COLCEN[c0] - wh[:2]))), 0),
                    nearest_fault=kf, dist_nearest_fault_m=round(df[kf], 0), dist_fault1_m=round(df["Fault 1"], 0),
                    dist_fault2_m=round(df["Fault 2"], 0), dist_fault3_m=round(df["Fault 3"], 0),
                    nearest_river=kr, dist_nearest_river_m=round(dr[kr], 0), dist_lake_edge_m=round(dlake_edge, 0),
                    wellhead_inside_lake=in_lake, td_dist_nearest_fault_m=round(min(dist_polyline(td[:2], v) for v in faults.values()), 0),
                    first_block_wells_json=wellsjson[w]["Blocks"][0]["block"], feed_zone_blocks_match_wells_json=fz_blocks_ok))
GEO = pd.DataFrame(geo)


def seg_intersect(a, b):
    """Intersection of the infinite lines through two fault segments (returns None if parallel)."""
    p, r = a[0], a[-1] - a[0]; q, s_ = b[0], b[-1] - b[0]
    den = r[0] * s_[1] - r[1] * s_[0]
    if abs(den) < 1e-9:
        return None
    t = ((q - p)[0] * s_[1] - (q - p)[1] * s_[0]) / den
    u = ((q - p)[0] * r[1] - (q - p)[1] * r[0]) / den
    return p + t * r, (0 <= t <= 1 and 0 <= u <= 1)


FINT = []
names = list(faults)
for i in range(len(names)):
    for j in range(i + 1, len(names)):
        res = seg_intersect(faults[names[i]], faults[names[j]])
        if res is not None:
            xy, within = res
            near = {w: round(float(np.hypot(*(W[w]["wh"][:2] - xy))), 0) for w in W}
            FINT.append(dict(fault_a=names[i], fault_b=names[j], x=round(float(xy[0]), 0), y=round(float(xy[1]), 0), within_both_segments=bool(within),
                             nearest_well=min(near, key=near.get), nearest_well_dist_m=min(near.values()),
                             dist_to_WELL_2_m=near["WELL-2"], dist_to_WELL_10_m=near["WELL-10"], dist_to_WELL_3_m=near["WELL-3"]))
FINT = pd.DataFrame(FINT)
FINT.to_csv(OUT / "fault_intersections.csv", index=False)
print(FINT.to_string())
GEO.to_csv(OUT / "well_geometry_check.csv", index=False)
print(GEO.T.to_string())

# ----------------------------------------------------------------------------------------------
# 7. Lateral structure tables
# ----------------------------------------------------------------------------------------------
hdr("7. Lateral structure tables (interpolation only)")
COMMON = [100, 0, -250, -500, -750, -1000, -1500, -2000, -2500]
ce = []
for e in COMMON:
    for w, r in W.items():
        T = float(interp_profile(w, e, "T")); P = float(interp_profile(w, e, "P"))
        if np.isnan(T):
            continue
        x, y = track_xy(w, e)
        ce.append(dict(elevation_mRL=e, well=w, x=round(x, 1), y=round(y, 1), T_degC=round(T, 1), P_bar=round(P, 1),
                       Tsat_degC=round(float(tsat_C(P)), 1), boiling_margin_degC=round(float(tsat_C(P)) - T, 1)))
CE = pd.DataFrame(ce)
CE.to_csv(OUT / "T_P_at_common_elevations.csv", index=False)
print(CE.pivot(index="well", columns="elevation_mRL", values="T_degC").to_string())
print(CE.pivot(index="well", columns="elevation_mRL", values="P_bar").to_string())

lc = []
for w, r in W.items():
    for rr in LAY.iloc[1:].itertuples():
        zc = rr.centre
        T = float(interp_profile(w, zc, "T"))
        if np.isnan(T):
            continue
        c, l, b = block_at(w, zc)
        if c is None or zc > SURF[c]:
            continue
        P = float(interp_profile(w, zc, "P"))
        infz = any(f[1] - 1e-6 <= zc <= f[0] + 1e-6 for f in r["fz"])
        behind = any(d[1] < zc <= d[0] + 1e-6 for d in r["design"] if not d[3])
        x, y = track_xy(w, zc)
        lc.append(dict(well=w, layer=rr.layer, layer_centre_mRL=zc, column=c, block=f"{c}{rr.layer:02d}", x=round(x, 1), y=round(y, 1),
                       T_degC=round(T, 1), P_bar=round(P, 1), Tsat_degC=round(float(tsat_C(P)), 1),
                       boiling_margin_degC=round(float(tsat_C(P)) - T, 1), in_feedzone=infz, behind_casing=behind))
LC = pd.DataFrame(lc)
LC.to_csv(OUT / "T_P_at_layer_centres.csv", index=False)
print(f"Layer-centre rows: {len(LC)} ({LC.groupby('well').size().to_dict()})")

# save state for plotting/report stage
STATE = dict(PT=PT, HS=HS, SUM=SUM, FZ=FZ, GEO=GEO, CE=CE, LC=LC, pc=pc, integ=integ)


# ----------------------------------------------------------------------------------------------
# 8. Figures
# ----------------------------------------------------------------------------------------------
hdr("8. Figures")
# colour-blind-safe palette (Okabe-Ito)
OI = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000", "#882255", "#44AA99"]
WCOL = {w: OI[i] for i, w in enumerate(WELLS)}
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 110})


def save(fig, name):
    fig.savefig(FIG / name, dpi=130, bbox_inches="tight")
    plt.close(fig)
    print("  wrote", name)


# --- well map ---------------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(11, 9))
for c, P_ in COLPOLY.items():
    ax.add_patch(MplPolygon(P_, closed=True, fill=False, ec="#bbbbbb", lw=0.4))
ax.add_patch(MplPolygon(lake, closed=True, fc="#9ecae1", ec="#3182bd", lw=0.6, alpha=0.8, label="Lake Tarawera"))
first = True
for k, segs in rivers.items():
    for sg in segs:
        ax.plot(sg[:, 0], sg[:, 1], color="#3182bd", lw=1.4, label="River" if first else None); first = False
for i, (k, v) in enumerate(faults.items()):
    ax.plot(v[:, 0], v[:, 1], "k-", lw=2, label="Fault" if i == 0 else None)
    mid = v.mean(axis=0); ax.text(v[0][0], v[0][1], " " + k, fontsize=8, va="bottom")
for _, r_ in FINT[FINT.within_both_segments].iterrows():
    ax.plot(r_.x, r_.y, marker="X", color="#D55E00", ms=11, mec="k", ls="", label="Fault intersection")
for w, r in W.items():
    ax.plot(r["track"][:, 0], r["track"][:, 1], "-", color=WCOL[w], lw=2)
    ax.plot(*r["wh"][:2], "o", color=WCOL[w], mec="k", ms=8)
    ax.annotate(w.replace("WELL-", "W"), r["wh"][:2], xytext=(6, 6), textcoords="offset points", fontweight="bold", fontsize=9)
h, l = ax.get_legend_handles_labels()
by = dict(zip(l, h)); ax.legend(by.values(), by.keys(), loc="lower right")
ax.set_xlim(2805000, 2829000); ax.set_ylim(6322000, 6343000); ax.set_aspect("equal")
ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
ax.set_title("Well locations, model grid columns (grey), faults, rivers and Lake Tarawera\n(deviated tracks W8, W9 drawn as lines from wellhead dot)")
save(fig, "well_map_grid_faults_water.png")

fig, ax = plt.subplots(figsize=(9, 8))   # zoomed map on the wells
for c, P_ in COLPOLY.items():
    ax.add_patch(MplPolygon(P_, closed=True, fill=False, ec="#cccccc", lw=0.5))
    if c in SURF:
        pass
for i, (k, v) in enumerate(faults.items()):
    ax.plot(v[:, 0], v[:, 1], "k-", lw=2, label="Fault" if i == 0 else None); ax.text(v[0][0] if k != "Fault 3" else 2817000, v[0][1] if k != "Fault 3" else 6333800, " " + k, fontsize=8)
first = True
for k, segs in rivers.items():
    for sg in segs:
        ax.plot(sg[:, 0], sg[:, 1], color="#3182bd", lw=1.4, label="River" if first else None); first = False
for _, r_ in FINT[FINT.within_both_segments].iterrows():
    ax.plot(r_.x, r_.y, marker="X", color="#D55E00", ms=11, mec="k", ls="", label="Fault intersection")
for w, r in W.items():
    row = SUM[SUM.well == w].iloc[0]
    ax.plot(*r["wh"][:2], "o", ms=14, mec="k", color=plt.cm.inferno(min(row.T_max_degC, 270) / 280))
    ax.annotate(f"{w.replace('WELL-', 'W')}\n{row.T_max_degC:.0f}°C", r["wh"][:2], xytext=(9, 4), textcoords="offset points", fontsize=8, fontweight="bold")
ax.set_xlim(2816500, 2825500); ax.set_ylim(6330500, 6339500); ax.set_aspect("equal")
h, l = ax.get_legend_handles_labels(); by = dict(zip(l, h)); ax.legend(by.values(), by.keys(), loc="lower right")
ax.set_title("Wells coloured by maximum surveyed temperature (label: well, Tmax)")
ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
save(fig, "well_map_zoom_tmax.png")

# --- profiles -----------------------------------------------------------------------------------
QCOL = {"HIGH": "#009E73", "MEDIUM": "#E69F00", "LOW": "#D55E00", "Synthetic": "#999999"}
for w, r in W.items():
    z, T, P = r["z"], r["T"], r["P"]
    srow = SUM[SUM.well == w].iloc[0]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 8), sharey=True)
    ts = tsat_C(P)
    a1.plot(T, z, "-o", color="#0072B2", ms=3, label="Survey T")
    a1.plot(ts, z, "--", color="#D55E00", label="Tsat(P) (IF97)")
    for k, f in enumerate(r["fz"]):
        a1.axhspan(f[1], f[0], color=QCOL.get(f[2], "#999"), alpha=0.35, label=f"Feed zone [{f[2]}]" + (" – pivot" if k == r["dom"] else ""))
        a2.axhspan(f[1], f[0], color=QCOL.get(f[2], "#999"), alpha=0.35)
    for d in r["design"]:
        col = "#7f7f7f" if not d[3] else "#c7c7c7"
        for a in (a1, a2):
            a.axhspan(d[1], d[0], xmin=0.0, xmax=0.02, color=col, alpha=0.9)
    a1.axhline(r["shoe"], color="k", ls=":", lw=1, label=f"Casing shoe {r['shoe']:.0f}")
    a1.axhline(r["wh"][2], color="#555", lw=0.8); a2.axhline(r["wh"][2], color="#555", lw=0.8)
    cs = SURF[find_column(*r["wh"][:2])]
    a1.axhline(cs, color="#999", lw=0.8, ls="-."); a2.axhline(cs, color="#999", lw=0.8, ls="-.")
    a1.text(a1.get_xlim()[1], cs, " model surface", va="bottom", ha="right", fontsize=7, color="#666")
    a1.plot(srow.T_max_degC, srow.T_max_elev_mRL, "k*", ms=12, label=f"Tmax {srow.T_max_degC:.0f}°C @ {srow.T_max_elev_mRL:.0f}")
    for lev, c_ in ((100, "#56B4E9"), (200, "#CC79A7")):
        a1.axvline(lev, color=c_, lw=0.8, ls=":")
    for a, b_ in r["isoth"]:
        a1.annotate("", xy=(T.min() - 5, b_), xytext=(T.min() - 5, a), arrowprops=dict(arrowstyle="-", color="#CC79A7", lw=3))
    m = PT[PT.well == w]
    a1.scatter(m[m.boiling_margin_degC < 10].T_degC, m[m.boiling_margin_degC < 10].elevation_mRL, s=70, facecolors="none", edgecolors="red", label="margin < 10 °C")
    a1.set_xlabel("Temperature (°C)"); a1.set_ylabel("Elevation (mRL)")
    a1.set_title(f"{w}: {srow.regime}", fontsize=9)
    a1.legend(fontsize=7, loc="lower left")
    a2.plot(P, z, "-o", color="#0072B2", ms=3, label="Survey P (wellbore column)")
    dz_ = r["fz"][r["dom"]]; mid = 0.5 * (dz_[0] + dz_[1])
    a2.plot(srow.pivot_P_bar, mid, "D", color="#009E73", mec="k", ms=9, label=f"Pivot (dominant FZ) {srow.pivot_P_bar:.1f} bar")
    for k, f in enumerate(r["fz"]):
        mm = 0.5 * (f[0] + f[1]); a2.plot(float(interp_profile(w, mm, "P")), mm, "s", color="#E69F00", mec="k", ms=6, label="Feed-zone P" if k == 0 else None)
    a2.plot(P_ATM, srow.water_level_est_mRL, "v", color="#D55E00", ms=9, label=f"Water level est. {srow.water_level_est_mRL:.0f} mRL (extrapolated)")
    ph = P[0] + (r["wh"][2] - z) * 0.0980665 * 1.0   # 1000 kg/m3 water-column reference slope from the top point
    a2.plot(P[0] + (z[0] - z) * 0.0980665, z, ":", color="#777", label="Cold-water gradient (0.098 bar/m) from top point")
    a2.set_xlabel("Pressure (bar)"); a2.set_title("Pressure: only feed-zone / pivot values are formation targets", fontsize=9)
    a2.legend(fontsize=7, loc="lower left")
    a2.set_xlim(left=0)
    fig.suptitle(f"{w}  NS verdict: {pc[pc.well == w].verdict.iloc[0]} ({pc[pc.well == w].confidence_level.iloc[0]}); grey strip = blank casing (dark) / slotted liner (light)", fontsize=9)
    fig.tight_layout()
    save(fig, f"profile_{w}.png")

# --- maps by elevation ---------------------------------------------------------------------------
MAP_ELEV = [0, -250, -500, -750, -1000, -1500, -2000, -2500]


def map_panel(ax, e, col, cmap, label, anomaly=False):
    d = CE[CE.elevation_mRL == e]
    if len(d) == 0:
        ax.set_visible(False); return
    vals = d[col].values.astype(float)
    ref = np.median(vals) if anomaly else 0
    vv = vals - ref
    for k, v in faults.items():
        ax.plot(v[:, 0], v[:, 1], "k-", lw=1)
    for sg in rivers["Tarawera River"]:
        ax.plot(sg[:, 0], sg[:, 1], color="#3182bd", lw=1)
    if len(d) >= 3:
        tri = mtri.Triangulation(d.x.values, d.y.values)
        lin = mtri.LinearTriInterpolator(tri, vv)
        xi, yi = np.meshgrid(np.linspace(2816500, 2825500, 150), np.linspace(6330500, 6339500, 150))
        zi = lin(xi, yi)
        vmin, vmax = (-np.nanmax(np.abs(vv)) - 1e-6, np.nanmax(np.abs(vv)) + 1e-6) if anomaly else (None, None)
        cf = ax.contourf(xi, yi, zi, 14, cmap=cmap, alpha=0.75, vmin=vmin, vmax=vmax)
        plt.colorbar(cf, ax=ax, shrink=0.8, label=label)
    sc = ax.scatter(d.x, d.y, c=vv, cmap=cmap, s=70, edgecolor="k", zorder=5,
                    vmin=(-np.nanmax(np.abs(vv)) - 1e-6 if anomaly else None), vmax=(np.nanmax(np.abs(vv)) + 1e-6 if anomaly else None))
    for _, r_ in d.iterrows():
        ax.annotate(f"{r_.well.replace('WELL-', 'W')}\n{r_[col]:.0f}" if col == "T_degC" else f"{r_.well.replace('WELL-', 'W')}\n{r_[col]:.1f}", (r_.x, r_.y), xytext=(5, 4), textcoords="offset points", fontsize=7)
    ax.set_xlim(2816500, 2825500); ax.set_ylim(6330500, 6339500); ax.set_aspect("equal")
    ax.set_title(f"{e:+d} mRL  ({len(d)} wells" + (f"; P median {ref:.1f} bar)" if anomaly else ")"), fontsize=9)
    ax.tick_params(labelsize=6)


for col, cmap, label, fname, anomaly in (("T_degC", "inferno", "T (°C)", "maps_temperature_by_elevation.png", False),
                                          ("P_bar", "RdBu_r", "P − median (bar)", "maps_pressure_by_elevation.png", True)):
    fig, axs = plt.subplots(2, 4, figsize=(18, 9))
    for ax, e in zip(axs.ravel(), MAP_ELEV):
        map_panel(ax, e, col, cmap, label, anomaly)
    fig.suptitle(("Temperature" if col == "T_degC" else "Pressure anomaly (value − median of wells at that elevation; labels are absolute P)") +
                 " at common elevations: linear interpolation inside the well hull only, no extrapolation; black = faults, blue = Tarawera River", fontsize=10)
    fig.tight_layout()
    save(fig, fname)

# --- cross-sections ------------------------------------------------------------------------------
SECTIONS = {
    "hot_axis_W2-W10-W3": (("WELL-10", "WELL-3"), 800.0, "Section through the hottest wells (W2, W10, W3; W2 projected ~570 m onto the W10–W3 line)"),
    "regional_W8-W1": (("WELL-8", "WELL-1"), 900.0, "Regional SW–NE section W8 → W1 (wells within 900 m of the line)"),
    "north_south_W4-W7": (("WELL-4", "WELL-7"), 1500.0, "North–south section W4 → W7 (wells within 1500 m of the line)"),
}
ZG = np.arange(300, -2900, -25.0)


def section_data(a, b, corridor):
    pa = W[a]["wh"][:2]; pb = W[b]["wh"][:2]
    v = pb - pa; L = np.hypot(*v); u = v / L
    sel = []
    for w, r in W.items():
        d = r["wh"][:2] - pa
        along = float(np.dot(d, u)); off = float(u[0] * d[1] - u[1] * d[0])
        if abs(off) <= corridor and -300 <= along <= L + 300:
            sel.append((w, along, off))
    sel.sort(key=lambda t: t[1])
    return sel, L


def interp_section(sel, arr):
    Zs = np.full((len(ZG), 1), np.nan)
    ds = np.array([s[1] for s in sel])
    out = np.full((len(ZG), 200), np.nan)
    xs = np.linspace(ds.min(), ds.max(), 200)
    for i, zg in enumerate(ZG):
        xv, yv = [], []
        for w, al, off in sel:
            val = float(interp_profile(w, zg, arr))
            if not np.isnan(val):
                xv.append(al); yv.append(val)
        if len(xv) >= 2:
            xv = np.array(xv); yv = np.array(yv)
            m = (xs >= xv.min()) & (xs <= xv.max())
            out[i, m] = np.interp(xs[m], xv, yv)
    return xs, out


SECINFO = {}
for name, (ends, corridor, title) in SECTIONS.items():
    sel, L = section_data(ends[0], ends[1], corridor)
    if len(sel) < 2:
        continue
    SECINFO[name] = sel
    fig, axs = plt.subplots(2, 1, figsize=(12, 10), sharex=True)
    for ax, arr, cmap, lab in ((axs[0], "T", "inferno", "T (°C)"), (axs[1], "P", "viridis", "P (bar)")):
        xs, F = interp_section(sel, arr)
        cf = ax.contourf(xs, ZG, F, 18, cmap=cmap)
        plt.colorbar(cf, ax=ax, label=lab)
        if arr == "T":
            cs = ax.contour(xs, ZG, F, levels=[100, 200, 250], colors=["#56B4E9", "#CC79A7", "white"], linewidths=1.2)
            ax.clabel(cs, fmt="%d°C", fontsize=7)
        for w, al, off in sel:
            r = W[w]
            ax.plot([al, al], [r["wh"][2], r["z"][-1]], "-", color="white", lw=2.5)
            ax.plot([al, al], [r["wh"][2], r["z"][-1]], "-", color="k", lw=0.8)
            for f in r["fz"]:
                ax.plot([al, al], [f[1], f[0]], "-", color=QCOL.get(f[2], "#999"), lw=6, solid_capstyle="butt")
            ax.plot(al, r["shoe"], "kv", ms=6)
            ax.text(al, r["wh"][2] + 40, f"{w.replace('WELL-','W')}\n({off:+.0f} m)", ha="center", fontsize=8, fontweight="bold")
            if arr == "T":
                ax.plot(al, SUM[SUM.well == w].T_max_elev_mRL.iloc[0], "w*", mec="k", ms=10)
        ax.set_ylim(-2900, 500); ax.set_ylabel("Elevation (mRL)")
    axs[1].set_xlabel("Distance along section (m); offset of each well from the line in brackets")
    axs[0].set_title(title + "\nInterpolation between wells at the same elevation only (blank = not covered by ≥2 wells). Coloured bars = feed zones (green HIGH, amber MEDIUM, red LOW), ▼ = casing shoe, ★ = Tmax", fontsize=9)
    fig.tight_layout()
    save(fig, f"section_{name}.png")
    print("   section", name, [(s[0], round(s[1]), round(s[2])) for s in sel])

# --- P-T state plot -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 8))
Tc = np.linspace(20, 370, 300)
ax.plot(Tc, psat_MPa(Tc + 273.15) * 10, "k-", lw=2, label="Saturation curve (IF97 Region 4)")
for w, r in W.items():
    ax.plot(r["T"], r["P"], "-", color=WCOL[w], lw=1.2, alpha=0.9)
    ax.plot(r["T"], r["P"], ".", color=WCOL[w], ms=4, label=w)
    for f in r["fz"]:
        mm = 0.5 * (f[0] + f[1])
        ax.plot(float(interp_profile(w, mm, "T")), float(interp_profile(w, mm, "P")), "s", color=WCOL[w], mec="k", ms=8)
ax.plot([], [], "s", color="w", mec="k", label="Feed-zone midpoints")
ax.set_xlabel("Temperature (°C)"); ax.set_ylabel("Pressure (bar)"); ax.set_yscale("log"); ax.set_ylim(0.9, 320)
ax.set_title("P–T state of all survey points vs the saturation curve\n(liquid is to the left/above; W10 is the only well touching the curve)")
ax.legend(fontsize=7, ncol=2, loc="lower right")
save(fig, "PT_state_vs_saturation.png")

# --- hydrostatic excess --------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(8, 7))
for w in WELLS:
    d = HS[(HS.well == w) & HS.single_phase_ok]
    ax.plot(d.excess_pct, 0.5 * (d.z_top_mRL + d.z_bottom_mRL), "-o", color=WCOL[w], label=w, ms=5)
ax.axvline(0, color="k", lw=1); ax.axvline(FIELD_EXCESS_COLD, color="#D55E00", ls="--", label=f"cold-window median {FIELD_EXCESS_COLD:.1f} %")
ax.axvline((10 / G - 1) * 100, color="#999", ls=":", label="g = 10 m/s² would give +2.0 %")
ax.set_xlabel("Gradient excess over IF97 liquid column (%)"); ax.set_ylabel("Window mid-elevation (mRL)")
ax.set_title("Column-density test: observed ΔP / IF97 ΔP − 1 (windows 345–990 m, single-phase, margin ≥ 10 °C)")
ax.legend(fontsize=7, ncol=2)
save(fig, "hydrostatic_excess_by_window.png")

# --- natural-state test --------------------------------------------------------------------------
fig, axs = plt.subplots(1, 4, figsize=(15, 4), sharey=False)
for ax, w in zip(axs, ["WELL-2", "WELL-3", "WELL-9", "WELL-10"]):
    r = W[w]; (zk, tpv), = r["tp"].items()
    yrs = [parse_date(k).year + (parse_date(k).dayofyear - 1) / 365.25 for k in tpv["Data"]]
    ax.plot(yrs, list(tpv["Data"].values()), "-o", color="#0072B2", label=f"transient gauge {zk} mRL")
    sy = parse_date(r["survey_date_P"]); sy = sy.year + (sy.dayofyear - 1) / 365.25
    ax.plot(sy, float(interp_profile(w, float(zk), "P")), "D", color="#D55E00", ms=10, label="survey P at gauge elev. (29/07/2015)")
    ax.set_title(f"{w}: {pc[pc.well == w].verdict.iloc[0]} [{pc[pc.well == w].confidence_level.iloc[0]}]", fontsize=9)
    ax.set_xlabel("Year"); ax.set_ylabel("P (bar)"); ax.legend(fontsize=7)
fig.suptitle("Survey pressure vs transient gauge: a natural-state survey should match the FIRST transient value (production start)", fontsize=10)
fig.tight_layout()
save(fig, "natural_state_pressure_test.png")


# --- conceptual sketch (INTERPRETATION, drawn from well data only) ---------------------------------
S_ = SUM.set_index("well")
fig, ax = plt.subplots(figsize=(10, 9))
for k, v in faults.items():
    ax.plot(v[:, 0], v[:, 1], "k-", lw=2.2)
for sg in rivers["Tarawera River"]:
    ax.plot(sg[:, 0], sg[:, 1], color="#3182bd", lw=1.5)
hot = ["WELL-10", "WELL-2", "WELL-3"]
cx = np.mean([W[w]["wh"][0] for w in hot]); cy = np.mean([W[w]["wh"][1] for w in hot])
ax.add_patch(plt.Circle((cx, cy), 900, color="#D55E00", alpha=0.25, label="Inferred upflow / hottest zone (Tmax >= 260 °C)"))
ax.plot(cx, cy, "*", color="#D55E00", ms=22, mec="k")
for w in ["WELL-6", "WELL-7", "WELL-9", "WELL-5"]:
    ax.annotate("", xy=W[w]["wh"][:2], xytext=(cx, cy), arrowprops=dict(arrowstyle="-|>", color="#E69F00", lw=3, alpha=0.9))
for w in ["WELL-8", "WELL-4"]:
    ax.annotate("", xy=(cx + (W[w]["wh"][0] - cx) * 0.55, cy + (W[w]["wh"][1] - cy) * 0.55), xytext=W[w]["wh"][:2],
                arrowprops=dict(arrowstyle="-|>", color="#0072B2", lw=3, alpha=0.9))
for w, r in W.items():
    row = S_.loc[w]
    ax.plot(*r["wh"][:2], "o", ms=12, mec="k", color=plt.cm.inferno(min(row.T_max_degC, 270) / 280))
    ax.annotate(f"{w.replace('WELL-', 'W')}: {row.T_max_degC:.0f}°C @ {row.T_max_elev_mRL:.0f}\n{row.regime.split(' (')[0]}", r["wh"][:2], xytext={"WELL-2": (-150, 22), "WELL-10": (-165, -34), "WELL-3": (10, -26), "WELL-1": (10, -26), "WELL-5": (8, 14)}.get(w, (8, -4)), textcoords="offset points", fontsize=7)
ax.plot([], [], color="#E69F00", lw=3, label="Inferred lateral outflow (-450...-750 mRL)")
ax.plot([], [], color="#0072B2", lw=3, label="Inferred cold recharge / downflow")
ax.plot([], [], "k-", lw=2, label="Faults (Mapping.json)"); ax.plot([], [], color="#3182bd", label="Tarawera River")
ax.legend(loc="lower right", fontsize=8)
ax.set_xlim(2816500, 2825500); ax.set_ylim(6330500, 6339500); ax.set_aspect("equal")
ax.set_title("Conceptual sketch from the well data only (plan view; INTERPRETATION - arrows are inferred, not measured)")
ax.set_xlabel("Easting (m)"); ax.set_ylabel("Northing (m)")
save(fig, "conceptual_model_sketch.png")

# ----------------------------------------------------------------------------------------------
# 9. Report (all numbers come from the tables above)
# ----------------------------------------------------------------------------------------------
hdr("9. Writing report")


def md(df, floatfmt=1):
    d = df.copy()
    out = ["| " + " | ".join(str(c) for c in d.columns) + " |", "|" + "|".join("---" for _ in d.columns) + "|"]
    for _, r in d.iterrows():
        cells = []
        for v in r:
            if isinstance(v, (float, np.floating)):
                cells.append("–" if np.isnan(v) else f"{v:.{floatfmt}f}")
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


pcs = pc.set_index("well"); GE = GEO.set_index("well")
f1x = FINT[(FINT.fault_a == "Fault 1") & (FINT.fault_b == "Fault 3")].iloc[0]
g = lambda w, c: S_.loc[w, c]
hot3 = ["WELL-10", "WELL-2", "WELL-3"]
HSo = HS[HS.single_phase_ok]
shallow = HS[(HS.z_top_mRL >= 0) & HS.single_phase_ok]
nexc = int(shallow.exceeds_max_liquid_gradient.sum())
deepw = HSo[HSo.z_top_mRL <= -1300]
zero_deep = sorted(set(HSo[(HSo.z_top_mRL <= -750) & (HSo.excess_pct.abs() <= 0.5)].well))
w2_2015 = pcs.loc["WELL-2", "P_at_2015-01-01_bar"]; w3_2015 = pcs.loc["WELL-3", "P_at_2015-01-01_bar"]
CEi = CE.set_index(["well", "elevation_mRL"])

t_temp = SUM[["well", "regime", "T_max_degC", "T_max_elev_mRL", "T_at_TD_degC", "inversion_size_degC", "inversion_top_mRL", "inversion_bottom_mRL",
              "cap_gradient_K_km", "deep_gradient_K_km", "elev_100C_first_mRL", "elev_200C_first_mRL", "elev_T90pct_Tmax_mRL",
              "min_boiling_margin_degC", "min_margin_elev_mRL"]].copy()
t_temp.columns = ["well", "regime", "Tmax °C", "Tmax mRL", "T@TD °C", "inversion K", "inv. top mRL", "inv. bottom mRL", "cap grad K/km", "deep grad K/km",
                  "100 °C mRL", "200 °C mRL", "0.9·Tmax mRL", "min margin °C", "at mRL"]


def pivot_tab(col):
    t = CE.pivot(index="well", columns="elevation_mRL", values=col)
    t = t[[c for c in COMMON if c in t.columns]].reset_index()
    t.columns = ["well"] + [f"{c:+d} mRL" for c in t.columns[1:]]
    return t


t_pres_e, t_temp_e = pivot_tab("P_bar"), pivot_tab("T_degC")
t_piv = SUM[["well", "water_level_est_mRL", "model_surface_mRL", "pivot_fz_mid_mRL", "pivot_quality", "pivot_P_bar", "pivot_T_degC"]].copy()
t_piv.columns = ["well", "water level est. mRL (extrapolated)", "model surface mRL", "pivot FZ mRL", "quality", "pivot P bar", "pivot T °C"]
t_ns = pc[["well", "role", "verdict", "confidence_level", "ns_confidence", "gauge_elev_mRL", "survey_P_at_gauge_bar", "first_transient_P_bar",
           "delta_survey_minus_first_bar", "dP_at_feedzones_vs_verified_bar"]].copy()
t_ns.columns = ["well", "role", "verdict", "confidence", "factor", "gauge mRL", "survey P at gauge bar", "first transient P bar", "Δ bar", "FZ P offset vs W2/W3 (bar; * pivot)"]
t_fz = FZ[["well", "feed", "mid_mRL", "quality", "T_degC", "P_bar", "Tsat_degC", "boiling_margin_degC", "block_at_mid", "is_pivot", "weight_T"]].copy()
t_fz.columns = ["well", "feed", "mid mRL", "quality", "T °C", "P bar", "Tsat °C", "margin °C", "block", "pivot", "weight"]
t_cls = PT.groupby("location_class").agg(points=("weight", "size"), mean_weight=("weight", "mean")).reset_index().round(2)
t_geo = GEO[["well", "wellhead_mRL", "model_surface_mRL", "wellhead_minus_model_surface_m", "max_horizontal_deviation_m", "dist_fault1_m", "dist_fault2_m",
             "dist_fault3_m", "dist_nearest_river_m", "dist_lake_edge_m"]].copy()
t_geo.columns = ["well", "wellhead mRL", "model surface mRL", "wellhead − surface m", "max deviation m", "→Fault 1 m", "→Fault 2 m", "→Fault 3 m", "→river m", "→lake edge m"]
t_hs = HS[["well", "z_top_mRL", "z_bottom_mRL", "T_range_degC", "min_boiling_margin_degC", "dP_obs_bar", "dP_expected_IF97_bar", "excess_pct",
           "apparent_density_kg_m3", "max_liquid_density_kg_m3", "exceeds_max_liquid_gradient"]].copy()
t_hs.columns = ["well", "top mRL", "bottom mRL", "T range °C", "min margin °C", "ΔP obs bar", "ΔP IF97 bar", "excess %", "apparent ρ", "max liquid ρ", "> max?"]
iso_rows = SUM[SUM.isothermal_sections.notna() & (SUM.isothermal_sections != "")][["well", "isothermal_sections"]]
flags = "; ".join(integ[integ.result.str.startswith("FLAG")].apply(lambda r: r.well + ": " + r.detail, axis=1))
wells_excess_list = ", ".join(sorted(shallow[shallow.exceeds_max_liquid_gradient].well))
n_pts = len(PT)
p15 = CE[CE.elevation_mRL == -1500].P_bar; p20 = CE[CE.elevation_mRL == -2000].P_bar
inv6 = SUM[SUM.well.isin(["WELL-2", "WELL-3", "WELL-5", "WELL-6", "WELL-7", "WELL-9"])].inversion_size_degC
nonw10_minmargin = SUM[SUM.well != "WELL-10"].min_boiling_margin_degC.min()
w2deep = HS[(HS.well == "WELL-2") & (HS.z_top_mRL <= -1100)].excess_pct.mean()
wh_top_T = PT.groupby("well").T_degC.first()
wh_top_P1 = PT[PT.well == "WELL-1"].P_bar.iloc[0]

R = f"""# Natural-state data interpretation – Teaching_Model_3D (TE14515 grid)

*Scope: natural state only. Inputs: `Baseline/Teaching_Model_3D` (read-only). Every number below is taken from a table in `Output/` written by `Scripts/natural_state_analysis.py`. Claims are tagged **[Obs]** observation, **[Int]** interpretation, **[Assump]** assumption. Elevations mRL, pressures bar (absolute), temperatures °C.*

## 1. Key conclusions for calibration

1. **[Obs] The surveys dated 29/07/2015 reproduce the pressure at the *start* of production, not at the survey date.** WELL-2 and WELL-3 survey pressures at the gauge elevation ({pcs.loc['WELL-2','survey_P_at_gauge_bar']:.1f} and {pcs.loc['WELL-3','survey_P_at_gauge_bar']:.1f} bar) equal the first transient values ({pcs.loc['WELL-2','first_transient_P_bar']:.1f} and {pcs.loc['WELL-3','first_transient_P_bar']:.1f} bar, 01/01/2014; Δ = {pcs.loc['WELL-2','delta_survey_minus_first_bar']:+.2f} and {pcs.loc['WELL-3','delta_survey_minus_first_bar']:+.2f} bar), whereas the gauges had already fallen to {w2_2015:.1f} and {w3_2015:.1f} bar by 01/01/2015 (`pressure_consistency_check.csv`). **Action:** treat the surveys as natural-state (pre-production) data and resolve the date label (§7.1). Verdicts and confidence (§2) set the weights.
2. **[Obs] The hottest, shallowest water is a compact zone around WELL-10, WELL-2 and WELL-3** (Tmax {g('WELL-10','T_max_degC'):.0f}, {g('WELL-2','T_max_degC'):.0f}, {g('WELL-3','T_max_degC'):.0f} °C at {g('WELL-10','T_max_elev_mRL'):.0f}, {g('WELL-2','T_max_elev_mRL'):.0f}, {g('WELL-3','T_max_elev_mRL'):.0f} mRL), falling to {g('WELL-6','T_max_degC'):.0f} (W6), {g('WELL-9','T_max_degC'):.0f} (W9), {g('WELL-7','T_max_degC'):.0f} (W7), {g('WELL-5','T_max_degC'):.0f} (W5) and only {g('WELL-8','T_max_degC'):.0f} (W8), {g('WELL-4','T_max_degC'):.0f} (W4). **[Int]** the upflow is under or next to this zone, which is {f1x.dist_to_WELL_2_m:.0f} m (WELL-2) from the Fault 1 × Fault 3 intersection. **Action:** put the deep heat/mass input here. WELL-2 is still {g('WELL-2','T_at_TD_degC'):.0f} °C at {g('WELL-2','td_mRL'):.0f} mRL, so no well constrains the source depth.
3. **[Obs] WELL-10 follows the boiling-point-for-depth curve** (minimum margin {g('WELL-10','min_boiling_margin_degC'):.1f} °C at {g('WELL-10','min_margin_elev_mRL'):.0f} mRL; {int(g('WELL-10','n_points_margin_lt5'))} points < 5 °C, {int(g('WELL-10','n_points_margin_lt10'))} < 10 °C). No other well is below {nonw10_minmargin:.1f} °C. **Action:** the model must reproduce a near-saturated or two-phase zone between about −75 and −350 mRL near WELL-10; match T together with Tsat(P).
4. **[Obs] Six wells (W2, W3, W5, W6, W7, W9) have temperature inversions of {inv6.min():.0f}–{inv6.max():.0f} K** with peaks at {SUM[SUM.well.isin(['WELL-2','WELL-3','WELL-5','WELL-6','WELL-7','WELL-9'])].T_max_elev_mRL.max():.0f} to {SUM[SUM.well.isin(['WELL-2','WELL-3','WELL-5','WELL-6','WELL-7','WELL-9'])].T_max_elev_mRL.min():.0f} mRL (W3 {g('WELL-3','inversion_size_degC'):.0f} K, W6 {g('WELL-6','inversion_size_degC'):.0f} K). **[Int]** a lateral outflow layer near −450…−750 mRL over cooler water. **Action:** the model needs a permeable layer there; deep isothermal sections get reduced weight (§3).
5. **[Obs] Pressure in the upper ~250–650 m of every well rises with depth faster than any liquid column can** (the first window exceeds the maximum liquid gradient in {nexc} of the {len(shallow)} wells that can be tested, W10 excluded as near-boiling: {wells_excess_list}). Cold wells W1, W4, W8 show {cold_shallow.min():.1f}–{cold_shallow.max():.1f} %, hot wells up to {HS.excess_pct.max():.0f} %. **Action:** do not use shallow pressure gradients as a target (weight 0.25 × location weight); use pressure at feed zones/pivots. Cause unresolved (§7.2).
6. **[Obs] Below about −750 mRL the pressure profiles of {', '.join(zero_deep)} equal the IAPWS-IF97 liquid column within ±0.5 %**, so units and gravity are consistent. Deep excess elsewhere is {deepw.excess_pct.min():.1f} to {deepw.excess_pct.max():.1f} % (W2 largest): well-specific, not a single field-wide offset.
7. **[Obs] Deep pressures are nearly uniform across the field** (−1500 mRL: {p15.min():.1f}–{p15.max():.1f} bar; −2000 mRL: {p20.min():.1f}–{p20.max():.1f} bar; W8 highest). **[Int]** laterally connected deep reservoir with near-hydrostatic head over ≥ 5 km; W8's high value is consistent with its cold, dense column rather than a lateral head.
8. **[Obs] Wellheads lie below the model surface by up to {abs(GEO.wellhead_minus_model_surface_m.min()):.0f} m** (WELL-1); estimated water levels span {SUM.water_level_est_mRL.min():.0f}–{SUM.water_level_est_mRL.max():.0f} mRL. **Action:** use an atmospheric-pressure water-table boundary; do not take the wellhead elevation as the grid top.
9. **[Obs] WELL-1 (east), WELL-4 (north) and WELL-8 (south-west) are cold:** W1 gradient {g('WELL-1','cap_gradient_K_km'):.0f} K/km (top) and {g('WELL-1','deep_gradient_K_km'):.0f} (deep); W4 only {g('WELL-4','cap_gradient_K_km'):.0f} K/km and {g('WELL-4','T_at_TD_degC'):.0f} °C at {g('WELL-4','td_mRL'):.0f} mRL; W8 is 64–68 °C from about −650 to −1500 mRL. **[Int]** cold recharge/downflow north and south-west, conductive east.

## 2. The data

**Inventory.** 10 wells (`Groups_data.json`), one T and one P survey each, all dated 29/07/2015 and flagged stable. Producers W2, W3, W9, W10 have transient gauges; W4, W5, W8 are injectors; W1, W6, W7 have no flow data. Faults, river, lake from `Mapping.json`; grid columns, {len(LAY)-1} active layers and surface from `gTE14515.dat`.

**Checks** (`data_integrity_check.csv`, {len(integ)} rows): JSON vs `.dat` for every T and P point and every transient series — no difference ({int((integ.result=='PASS').sum())} PASS rows in total); all survey points lie within the track; every feed zone is below the casing shoe in slotted/open hole and contains a survey point; all survey and feed-zone blocks agree with `gTE14515_wells.json`. Flag: {flags} (design/track inconsistency; no survey affected). Notes: feed-zone values are labelled "m" but are elevations (mRL); JSON track/design/feed-zone entries mix lists and strings (parsed). IAPWS-IF97 Regions 4 and 1 were validated against published values (errors ≪ 0.1 K; `run_console_log.txt`).

**Natural-state timing and verdicts** (`pressure_consistency_check.csv`):

{md(t_ns, 2)}

Method: survey P interpolated to the gauge elevation vs the **first** transient value. Valid test and |Δ| ≤ 0.5 bar ⇒ *natural state* (high, factor 1.0). WELL-9 is {pcs.loc['WELL-9','delta_survey_minus_first_bar']:+.2f} bar *above* the first value: opposite to drawdown, so not production-disturbed but unexplained ⇒ *probably natural state* (medium, 0.75). WELL-10's series starts {pcs.loc['WELL-10','first_transient_date']}, after production began ({pcs.loc['WELL-10','first_flow_date']}), so the test is invalid; the survey is {pcs.loc['WELL-10','delta_survey_minus_first_bar']:+.1f} bar above the 2020 value (sign expected from later drawdown) and its feed-zone P is {pcs.loc['WELL-10','dP_pivot_fz_vs_verified_bar']:+.1f} bar from W2/W3 ⇒ *probably natural state* (low, 0.6). Wells without gauges are compared at feed-zone elevations with W2/W3. **Injectors** ran at {pcs.loc['WELL-4','first_flow_kg_s']:.0f}/{pcs.loc['WELL-5','first_flow_kg_s']:.0f}/{pcs.loc['WELL-8','first_flow_kg_s']:.0f} kg/s (W4/W5/W8) of water at ≈ {pcs.loc['WELL-4','injectate_T_approx_degC']:.0f} °C [Assump: cp 4.19 kJ/kg/K] from 01/01/2014. W4 and W8 profiles are colder than the injectate (no warm-injection signature) ⇒ probably natural (low, 0.6); W5's feed-zone T ({g('WELL-5','pivot_T_degC'):.0f} °C) exceeds the injectate, so cooling could be hidden ⇒ possibly disturbed (0.4). No injector has a gauge. WELL-1 is {pcs.loc['WELL-1','dP_pivot_fz_vs_verified_bar']:+.1f} bar from W2/W3 at its pivot (probably natural, medium, 0.75), as are W6 and W7.

![natural-state test](figures/natural_state_pressure_test.png)

![well map](figures/well_map_grid_faults_water.png)

**Geometry** (`well_geometry_check.csv`):

{md(t_geo, 0)}

W8 and W9 are deviated ({GE.loc['WELL-8','max_horizontal_deviation_m']:.0f} and {GE.loc['WELL-9','max_horizontal_deviation_m']:.0f} m at TD); blocks follow the track. No wellhead is inside the lake; the nearest lake edge is {GEO.dist_lake_edge_m.min():.0f} m away.

## 3. Temperature structure

Per-well characteristics (`well_summary.csv`). Cap gradient = top 300 m, deep = bottom 500 m (K/km, positive = T rises with depth). Inversion = largest drop below an earlier maximum (threshold 3 K). 0.9·Tmax is an operational marker of the cap base (hot wells only).

{md(t_temp, 0)}

**Isothermal sections** (≥ 300 m within 3 K):

{md(iso_rows, 0)}

**Wellbore check.** Inversion peaks lie below the casing shoe (open hole) in W2, W3, W6, W7 and W9, and behind casing in W5 (peak {g('WELL-5','T_max_elev_mRL'):.0f}, shoe {g('WELL-5','casing_shoe_mRL'):.0f}); internal flow could alter the open-hole parts. **[Int]** The long isothermal sections below or between feed zones (e.g. W3, W6, W7) may be partly wellbore convection; open-hole points inside them get weight 0.5 not 0.8 and are not used to infer reservoir uniformity. The peaks lie close to a feed zone in W3 (peak {g('WELL-3','T_max_elev_mRL'):.0f}, feed zone −335), W6 ({g('WELL-6','T_max_elev_mRL'):.0f} vs −475), W7 ({g('WELL-7','T_max_elev_mRL'):.0f} vs −630) and W2 ({g('WELL-2','T_max_elev_mRL'):.0f} vs −745 and −260), which favours a reservoir origin there. In W5 the peak ({g('WELL-5','T_max_elev_mRL'):.0f} mRL) lies ~650 m above its shallowest feed zone (−1325), but it is behind casing, where T is closer to conductive formation T, so it is a formation feature without a nearby feed zone (an open question for the outflow depth there); W9's peak ({g('WELL-9','T_max_elev_mRL'):.0f}) lies between its feed zones (−535, −838).

**Temperature at common elevations** (`T_P_at_common_elevations.csv`, interpolation only):

{md(t_temp_e, 0)}

![T maps](figures/maps_temperature_by_elevation.png)

![hot-axis section](figures/section_hot_axis_W2-W10-W3.png)

![regional section](figures/section_regional_W8-W1.png)

The N–S section (W4→W7) is also in `figures/`. Sections interpolate between wells at the same elevation only: they are schematic between wells, not constrained images.

![P–T state](figures/PT_state_vs_saturation.png)

**Heat-flux contrast (assumption-based).** With an assumed {K_ROCK_ASSUMED:.1f} W/m/K [Assump, not in the data] the cap gradients give {g('WELL-4','cond_heat_flux_mW_m2_ASSUMED_k2'):.0f} (W4), {g('WELL-1','cond_heat_flux_mW_m2_ASSUMED_k2'):.0f} (W1), {g('WELL-8','cond_heat_flux_mW_m2_ASSUMED_k2'):.0f} (W8) mW/m² at the cold margins and {g('WELL-2','cond_heat_flux_mW_m2_ASSUMED_k2'):.0f} (W2) to {g('WELL-3','cond_heat_flux_mW_m2_ASSUMED_k2'):.0f} (W3) mW/m² in the core. The core values include a near-surface hydrological gradient and are not pure conduction; use only as a margin-to-core contrast.

## 4. Pressure structure

**Rule applied:** wellbore pressure equals formation pressure only at feed zones; elsewhere it follows the fluid column, pivoting at the dominant feed zone (highest quality tag, ties by higher survey T). Wellhead pressures of ≈ 1–11 bar are not reservoir pressure.

Pressure at common elevations (bar):

{md(t_pres_e, 1)}

Pivot pressures and water levels (`well_summary.csv`):

{md(t_piv, 1)}

![P maps](figures/maps_pressure_by_elevation.png)

**Hydrostatic checks** (`hydrostatic_check.csv`; windows ≥ 300 m, single-phase only, margin ≥ {MARGIN_MIN:.0f} °C; W10 excluded). Rounding (0.1 bar) limits each window to about ±0.4 %.

{md(t_hs, 1)}

**[Obs]** Windows starting near the surface need an apparent density of {shallow.apparent_density_kg_m3.min():.0f}–{shallow.apparent_density_kg_m3.max():.0f} kg/m³ against a maximum liquid density of {shallow.max_liquid_density_kg_m3.min():.0f}–{shallow.max_liquid_density_kg_m3.max():.0f}, so these profiles are not static single-phase columns. Cold windows (max T < {COLD_T:.0f} °C): shallow median {cold_shallow.median():.1f} % (range {cold_shallow.min():.1f}–{cold_shallow.max():.1f}), deeper median {cold_deeper.median():.1f} %; overall median {FIELD_EXCESS_COLD:.1f} %. For scale, g = 10 m/s² would give +{(10/G-1)*100:.1f} %. **[Int]** The shallow excess common to the cold wells may be a datum/convention issue (not physical overpressure), but cannot explain the larger hot-well excess or why W5–W7 are exact at depth.

![hydrostatic excess](figures/hydrostatic_excess_by_window.png)

**Explanations checked:** units (bar and mRL agree between JSON and `.dat`); density (IF97 validated, Region 1 used on liquid windows only); datum (windows are elevation differences, so a constant offset cancels); wellbore (a static column cannot exceed the cold-water gradient). None resolves the hot-well excess.

## 5. Conceptual model

![conceptual sketch](figures/conceptual_model_sketch.png)

| Feature | Interpretation (strength) | Supporting observation |
|---|---|---|
| Upflow | Under/near W10–W2–W3, close to Fault 1 × Fault 3 ({f1x.dist_to_WELL_2_m:.0f} m from W2) (moderate) | Tmax {g('WELL-10','T_max_degC'):.0f}/{g('WELL-2','T_max_degC'):.0f}/{g('WELL-3','T_max_degC'):.0f} °C; W2 T@TD {g('WELL-2','T_at_TD_degC'):.0f} °C and rising; W10 on the boiling curve; 200 °C at {g('WELL-10','elev_200C_first_mRL'):.0f}/{g('WELL-2','elev_200C_first_mRL'):.0f}/{g('WELL-3','elev_200C_first_mRL'):.0f} mRL |
| Outflow | Lateral layer ≈ −450…−750 mRL towards S, SE, E (W6, W7, W9, W5) (moderate; direction inferred from the T pattern only) | inversion peaks {g('WELL-6','T_max_elev_mRL'):.0f} (W6), {g('WELL-7','T_max_elev_mRL'):.0f} (W7), {g('WELL-9','T_max_elev_mRL'):.0f} (W9), {g('WELL-5','T_max_elev_mRL'):.0f} (W5) mRL (near feed zones except W5); peak T {g('WELL-6','T_max_degC'):.0f}/{g('WELL-9','T_max_degC'):.0f}/{g('WELL-7','T_max_degC'):.0f}/{g('WELL-5','T_max_degC'):.0f} °C |
| Cap | Conductive top; thin over W10, deeper away (moderate) | cap gradients {SUM[SUM.well.isin(hot3)].cap_gradient_K_km.min():.0f}–{SUM[SUM.well.isin(hot3)].cap_gradient_K_km.max():.0f} K/km (core); 0.9·Tmax at {g('WELL-10','elev_T90pct_Tmax_mRL'):.0f} (W10), {g('WELL-3','elev_T90pct_Tmax_mRL'):.0f} (W3), {g('WELL-2','elev_T90pct_Tmax_mRL'):.0f} (W2), {g('WELL-6','elev_T90pct_Tmax_mRL'):.0f} (W6), {g('WELL-9','elev_T90pct_Tmax_mRL'):.0f} (W9) mRL |
| Two-phase zone | Small, near W10, ≈ −75…−350 mRL (one well) | margin {g('WELL-10','min_boiling_margin_degC'):.1f} °C |
| Cold recharge/downflow | North (W4), south-west (W8) (weak–moderate) | W4 {g('WELL-4','cap_gradient_K_km'):.0f} K/km; W8 64–68 °C at −650…−1500 mRL; river {GE.loc['WELL-8','dist_nearest_river_m']:.0f} m from W8; lake ≥ {GEO.dist_lake_edge_m.min():.0f} m away, so no direct lake influence is seen |
| Conductive margin | East (W1) (moderate) | monotonic T, {g('WELL-1','T_at_TD_degC'):.0f} °C at {g('WELL-1','td_mRL'):.0f} mRL |
| Deep pressure connection | Laterally connected, near-hydrostatic (moderate) | P at −1500 mRL {p15.min():.1f}–{p15.max():.1f} bar |
| Faults | Hot wells lie within {GE.loc['WELL-6','dist_fault1_m']:.0f}–{GE.loc['WELL-3','dist_fault1_m']:.0f} m of Fault 1 (W6, W2, W10, W3) and within {GE.loc['WELL-2','dist_fault3_m']:.0f}/{GE.loc['WELL-10','dist_fault3_m']:.0f} m of Fault 3 (W2/W10); W9/W7 are {GE.loc['WELL-9','dist_fault2_m']:.0f}/{GE.loc['WELL-7','dist_fault2_m']:.0f} m from Fault 2 (spatial association only; dip and permeability unknown) | `well_geometry_check.csv`, `fault_intersections.csv` |

## 6. Recommended calibration set-up

**Targets and weights** (`NS_calibration_targets.csv`, {n_pts} survey points; weight = location weight × natural-state factor; extra columns `weight_P`, `P_role`):

{md(t_cls, 2)}

Scheme: feed zone 1.0; near wellhead (top 50 m) 0.3; behind casing 0.5; open hole 0.8, **reduced to 0.5 inside isothermal sections ≥ 300 m** (possible wellbore flow); above model surface 0 (none occur); multiplied by the factor 1.0 / 0.75 / 0.6 / 0.4. Pressure outside feed zones gets 0.25 × the location weight because it reflects the wellbore column. **Feed-zone targets** (`feedzones.csv`; all blocks match `gTE14515_wells.json`):

{md(t_fz, 1)}

**Heat and mass.** Deep hot input near W10–W2–W3 (T ≥ {g('WELL-2','T_at_TD_degC'):.0f} °C at depth in W2); the data give no flux, only the temperatures to match. Cooling/recharge from the north (W4) and south-west (W8). Conductive heat loss at the top (§3 contrast). **Permeability structure implied [Int]:** low-permeability cap above roughly −150…−450 mRL (steep T gradient; steeper-than-hydrostatic P); a permeable lateral layer near −450…−750 mRL (inversion peaks near feed zones in W2, W3, W6, W7; HIGH-quality feed zones at W2 −745, W3 −335, W6 −475/−1125, W7 −630/−1325, W8 −1125); a laterally connected, near-hydrostatic deep column; vertical connection in the hot core. **Boundaries:** atmospheric pressure at the water table (estimated {SUM.water_level_est_mRL.min():.0f}–{SUM.water_level_est_mRL.max():.0f} mRL); near-surface T about {wh_top_T.min():.0f}–{wh_top_T.max():.0f} °C at the wellheads; lateral sides open to deep pressure support; no flux can be derived for the base.

`T_P_at_layer_centres.csv` gives {len(LC)} well/layer pairs for block-wise comparison.

## 7. Data quality issues and open questions (ranked by effect on calibration)

1. **Survey date vs content.** Surveys dated 29/07/2015 match the 01/01/2014 pressure ({pcs.loc['WELL-2','survey_P_at_gauge_bar']:.1f} bar vs {w2_2015:.1f} bar at W2 on 01/01/2015). Either the date is wrong or the data are reconstructed pre-production values. Confirm: the whole natural-state calibration depends on it.
2. **Pressure gradients steeper than a water column** in the upper part of all wells (hot wells up to {HS.excess_pct.max():.0f} %; cold wells {cold_shallow.min():.1f}–{cold_shallow.max():.1f} %) and a well-dependent excess at depth (W2 mean {w2deep:.1f} % below −1100 mRL; W5–W7 ≈ 0). Obtain survey conditions (static or flowing, gauge reference) before using non-feed-zone pressure.
3. **WELL-10:** no valid natural-state test; near-boiling, so P and T errors are coupled; feed-zone quality is "Synthetic".
4. **Injectors W4, W5, W8:** no gauges; W5 possibly disturbed; W4's survey ends at {g('WELL-4','td_mRL'):.0f} mRL while its liner is designed to −3251 mRL.
5. **WELL-1 pressure** is {pcs.loc['WELL-1','dP_pivot_fz_vs_verified_bar']:+.1f} bar from W2/W3 at its pivot; its wellhead is {abs(GE.loc['WELL-1','wellhead_minus_model_surface_m']):.0f} m below the grid surface and reads {wh_top_P1:.1f} bar at the wellhead: check the datum.
6. **Grid surface vs wellhead elevation** differs by up to {abs(GEO.wellhead_minus_model_surface_m.min()):.0f} m; confirm topography.
7. Feed-zone units labelled "m" but are mRL; mixed JSON list/string formats; dominant feed zone not given (chosen by quality, then temperature).
8. Faults have no dip or permeability; no lake data (lake assumed unimportant only because it is far).
9. Rock conductivity for heat flux is assumed ({K_ROCK_ASSUMED:.1f} W/m/K).

## 8. Output file index

| File | Content |
|---|---|
| `Scripts/natural_state_analysis.py` | regenerates everything except the hand-written `process_log.md` |
| `Output/NS_data_interpretation.md` | this report |
| `Output/NS_calibration_targets.csv` | {n_pts} survey points, weights, blocks (+ `weight_P`, `P_role`, `location_class`) |
| `Output/feedzones.csv` | feed-zone T/P/margin, blocks, weights |
| `Output/T_P_at_layer_centres.csv` | T, P at layer centres, no extrapolation |
| `Output/T_P_at_common_elevations.csv` | T, P at common elevations |
| `Output/well_summary.csv` | per-well characteristics and regime |
| `Output/pressure_consistency_check.csv` | natural-state timing tests and verdicts |
| `Output/hydrostatic_check.csv` | column-density windows |
| `Output/well_geometry_check.csv`, `Output/fault_intersections.csv` | geometry, distances, fault intersections |
| `Output/data_integrity_check.csv` | all integrity checks |
| `Output/figures/` | well maps, 10 profiles, T/P maps, 3 sections, P–T plot, hydrostatic excess, natural-state test, conceptual sketch |
| `Output/run_console_log.txt` | console output of the final run |
| `Output/process_log.md` | prompt, rules applied, steps |
"""
open(OUT / "NS_data_interpretation.md", "w", encoding="utf-8").write(R)
print("  wrote NS_data_interpretation.md", len(R), "chars")
print("\nDONE")
