"""Natural-state data interpretation for geothermal model calibration.

Run from the repository root:

    python Scripts/natural_state_analysis.py

Reads only Baseline/ (read-only) and regenerates every output in Output/ and Output/figures/,
including the input inventory and the report Output/NS_data_interpretation.md.
Production / injection records are used only to test whether the surveys represent the
natural state (CLAUDE.md section 4.1).
"""
import ast
import json
import math
import os
import re
import sys
from datetime import datetime

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402
from matplotlib.tri import Triangulation, LinearTriInterpolator  # noqa: E402

BASE = "Baseline"
OUT = "Output"
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

LOG = []
REF_TOP = {}


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.append(s)


# ======================================================================================
# Analysis parameters (all written to Output/analysis_parameters.csv)
# ======================================================================================
G = 9.80665                 # m/s2
P_ATM = 1.01325             # bar
P_CRIT = 220.64             # bar
WIN_MIN = 300.0             # m, minimum hydrostatic window / gradient interval
MARGIN_LIQ = 10.0           # degC, boiling margin for single-phase liquid windows
P_LIQ = 1.5                 # bar, minimum pressure for liquid-filled points
MARGIN_FLAG = 5.0           # degC, low boiling-margin flag
INV_MIN = 3.0               # degC, minimum inversion size
ISO_GRAD = 10.0             # degC/km, |dT/dz| below which consecutive points are 'isothermal'
ISO_LEN = 300.0             # m, isothermal section length that reduces weights
WH_DIST = 30.0              # m, distance below wellhead for weight 0.3
COLD_TMAX = 100.0           # degC, cold-well definition
GAUGE_TOL = 0.1             # bar, reading precision used for the gauge test (two values at 0.1 bar)
PCONN_TOL = 2.0             # bar, wells further than this below the median P at -1000 mRL are reported apart
TREND_TOL = 5.0             # bar, negative residual against the passed-well trend that flags disturbance
CAP_GRAD = 100.0            # degC/km, conductive cap gradient threshold
K_ROCK = 2.5                # W/m/K, assumed rock thermal conductivity (typical volcanic rock)
P_INJ_ASSUMED = 10.0        # bar, pressure used to turn injectate enthalpy into temperature
COLD_COL_T = 20.0           # degC, temperature of the reference cold water column at feed zones
INJ_ISO_TOL = 10.0          # degC, isothermal section within this of injectate T = injection signature
COMMON_ELEVS = [200, 0, -250, -500, -750, -1000, -1250, -1500, -2000, -2500]
MAP_ELEVS = [0, -250, -500, -750, -1000, -1500]
LOC_W = {"above_surface": 0.0, "wellhead_or_top_block": 0.3, "feed_zone": 1.0,
         "behind_casing": 0.5, "open_hole_isothermal": 0.5, "open_hole": 0.8}
CONF_F = {"high": 1.0, "medium": 0.75, "low": 0.5}

PARAMS = [
    ("g_m_s2", G, "gravity used for IF97 hydrostatic integration"),
    ("P_atm_bar", P_ATM, "atmospheric pressure (water-level extrapolation target)"),
    ("P_crit_bar", P_CRIT, "critical pressure; no Tsat above"),
    ("window_min_m", WIN_MIN, "minimum hydrostatic window and gradient interval"),
    ("liquid_margin_degC", MARGIN_LIQ, "boiling margin required throughout a hydrostatic window"),
    ("liquid_P_bar", P_LIQ, "minimum P for liquid-filled points"),
    ("low_margin_flag_degC", MARGIN_FLAG, "boiling margin flag threshold"),
    ("inversion_min_degC", INV_MIN, "minimum inversion reported"),
    ("isothermal_grad_degC_per_km", ISO_GRAD, "|dT/dz| below which adjacent points are isothermal"),
    ("isothermal_len_m", ISO_LEN, "isothermal length that triggers weight 0.5"),
    ("wellhead_dist_m", WH_DIST, "points within this of the wellhead get weight 0.3"),
    ("cold_Tmax_degC", COLD_TMAX, "cold-well definition"),
    ("gauge_tol_bar", GAUGE_TOL, "gauge-test match tolerance (reading precision)"),
    ("pconn_tol_bar", PCONN_TOL, "P at -1000 mRL further below the median than this is reported separately"),
    ("trend_tol_bar", TREND_TOL, "negative trend residual that flags possible drawdown"),
    ("cap_gradient_degC_per_km", CAP_GRAD, "local gradient above which an interval is counted as conductive cap"),
    ("k_rock_W_per_mK", K_ROCK, "ASSUMED rock conductivity (typical volcanic rock; not from the data)"),
    ("P_injectate_bar", P_INJ_ASSUMED, "pressure used to convert injectate enthalpy to temperature"),
    ("cold_column_T_degC", COLD_COL_T, "reference cold column for feed-zone overpressure"),
    ("inj_iso_tol_degC", INJ_ISO_TOL, "isothermal section this close to injectate T = injection signature"),
] + [(f"weight_{k}", v, "location weight") for k, v in LOC_W.items()] + \
    [(f"conf_factor_{k}", v, "natural-state confidence factor") for k, v in CONF_F.items()]


# ======================================================================================
# IAPWS-IF97 (Region 1 and Region 4), validated against the IF97 verification tables
# ======================================================================================
_R = 0.461526  # kJ/kg/K
_N4 = [None, 0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5,
       -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6,
       -0.23855557567849, 0.65017534844798e3]
_I1 = np.array([0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8,
                21, 23, 29, 30, 31, 32], dtype=float)
_J1 = np.array([-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10,
                -8, -11, -6, -29, -31, -38, -39, -40, -41], dtype=float)
_n1 = np.array([0.14632971213167, -0.84548187169114, -0.37563603672040e1, 0.33855169168385e1,
                -0.95791963387872, 0.15772038513228, -0.16616417199501e-1, 0.81214629983568e-3,
                0.28319080123804e-3, -0.60706301565874e-3, -0.18990068218419e-1,
                -0.32529748770505e-1, -0.21841717175414e-1, -0.52838357969930e-4,
                -0.47184321073267e-3, -0.30001780793026e-3, 0.47661393906987e-4,
                -0.44141845330846e-5, -0.72694996297594e-15, -0.31679644845054e-4,
                -0.28270797985312e-5, -0.85205128120103e-9, -0.22425281908000e-5,
                -0.65171222895601e-6, -0.14341729937924e-12, -0.40516996860117e-6,
                -0.12734301741641e-8, -0.17424871230634e-9, -0.68762131295531e-18,
                0.14478307828521e-19, 0.26335781662795e-22, -0.11947622640071e-22,
                0.18228094581404e-23, -0.93537087292458e-25])


def tsat_K_MPa(p):
    n = _N4
    b = p ** 0.25
    e = b * b + n[3] * b + n[6]
    f = n[1] * b * b + n[4] * b + n[7]
    g = n[2] * b * b + n[5] * b + n[8]
    d = 2 * g / (-f - math.sqrt(f * f - 4 * e * g))
    return (n[10] + d - math.sqrt((n[10] + d) ** 2 - 4 * (n[9] + n[10] * d))) / 2


def psat_MPa_K(T):
    n = _N4
    th = T + n[9] / (T - n[10])
    A = th * th + n[1] * th + n[2]
    B = n[3] * th * th + n[4] * th + n[5]
    C = n[6] * th * th + n[7] * th + n[8]
    return (2 * C / (-B + math.sqrt(B * B - 4 * A * C))) ** 4


def tsat_c(p_bar):
    """Saturation temperature in degC for P in bar; NaN above the critical pressure or below triple."""
    if p_bar is None or not np.isfinite(p_bar):
        return float("nan")
    p = p_bar / 10.0
    if p < 611.213e-6 or p > 22.064:
        return float("nan")
    return tsat_K_MPa(p) - 273.15


def _r1(T_K, p_MPa):
    pi = p_MPa / 16.53
    tau = 1386.0 / T_K
    a = 7.1 - pi
    b = tau - 1.222
    gpi = np.sum(-_n1 * _I1 * a ** (_I1 - 1) * b ** _J1)
    gtau = np.sum(_n1 * a ** _I1 * _J1 * b ** (_J1 - 1))
    v = _R * T_K / p_MPa * pi * gpi / 1000.0  # m3/kg (R in kJ -> /1000 with p in MPa*1000 kPa)
    h = _R * T_K * tau * gtau
    return v, h


def v1(T_K, p_MPa):
    pi = p_MPa / 16.53
    tau = 1386.0 / T_K
    gpi = np.sum(-_n1 * _I1 * (7.1 - pi) ** (_I1 - 1) * (tau - 1.222) ** _J1)
    return _R * 1000.0 * T_K / (p_MPa * 1e6) * pi * gpi  # m3/kg


def h1(T_K, p_MPa):
    pi = p_MPa / 16.53
    tau = 1386.0 / T_K
    gtau = np.sum(_n1 * (7.1 - pi) ** _I1 * _J1 * (tau - 1.222) ** (_J1 - 1))
    return _R * T_K * tau * gtau  # kJ/kg


def rho_liq(T_c, p_bar):
    return 1.0 / v1(T_c + 273.15, p_bar / 10.0)


def validate_if97():
    rows = []
    for T, p, v_ref, h_ref in [(300, 3, 0.100215168e-2, 0.115331273e3),
                               (300, 80, 0.971180894e-3, 0.184142828e3),
                               (500, 3, 0.120241800e-2, 0.975542239e3)]:
        v = v1(T, p)
        h = h1(T, p)
        rows.append(("Region1 v", f"T={T}K p={p}MPa", v_ref, v, abs(v / v_ref - 1)))
        rows.append(("Region1 h", f"T={T}K p={p}MPa", h_ref, h, abs(h / h_ref - 1)))
    for p, T_ref in [(0.1, 372.755919), (1.0, 453.035632), (10.0, 584.149488)]:
        T = tsat_K_MPa(p)
        rows.append(("Region4 Tsat", f"p={p}MPa", T_ref, T, abs(T - T_ref)))
    for T, p_ref in [(300, 0.353658941e-2), (500, 0.263889776e1), (600, 0.123443146e2)]:
        p = psat_MPa_K(T)
        rows.append(("Region4 Psat", f"T={T}K", p_ref, p, abs(p / p_ref - 1)))
    for p_bar, T_ref in [(1.01325, 99.97), (100.0, 311.0)]:
        T = tsat_c(p_bar)
        rows.append(("Tsat reference", f"p={p_bar}bar", T_ref, T, abs(T - T_ref)))
    df = pd.DataFrame(rows, columns=["quantity", "state", "reference", "computed", "error"])
    df["error_K_equiv"] = np.nan
    for i, r in df.iterrows():
        if r.quantity.startswith("Region4 Tsat") or r.quantity == "Tsat reference":
            df.loc[i, "error_K_equiv"] = r.error
    bad = df[((df.quantity.isin(["Region4 Tsat", "Tsat reference"])) & (df.error > 0.1)) |
             ((df.quantity.str.startswith("Region1")) & (df.error > 1e-6)) |
             ((df.quantity == "Region4 Psat") & (df.error > 1e-6))]
    if len(bad):
        raise SystemExit(f"IF97 validation failed:\n{bad}")
    df.to_csv(os.path.join(OUT, "if97_validation.csv"), index=False, float_format="%.9g")
    log(f"IF97 validation passed ({len(df)} checks; max Tsat error "
        f"{df.loc[df.error_K_equiv.notna(), 'error'].max():.4f} K)")


def T_from_h(h_kJkg, p_bar):
    lo, hi = 274.0, 620.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if h1(mid, p_bar / 10.0) > h_kJkg:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi) - 273.15


def integrate_P(zt, zb, Pt, zs, Ts):
    """Integrate dP/dz = rho_IF97(T(z), P) g downward from zt to zb; zs descending, Ts matching."""
    n_ = max(int(math.ceil((zt - zb) / 1.0)), 1)
    h = (zt - zb) / n_
    P = Pt
    z = zt
    zs = np.asarray(zs, float)
    Ts = np.asarray(Ts, float)
    for _ in range(n_):
        T1 = float(np.interp(z, zs[::-1], Ts[::-1]))
        k1 = rho_liq(T1, P) * G / 1e5
        T2 = float(np.interp(z - h, zs[::-1], Ts[::-1]))
        k2 = rho_liq(T2, P + h * k1) * G / 1e5
        P += 0.5 * h * (k1 + k2)
        z -= h
    return P


# ======================================================================================
# 1. Inventory
# ======================================================================================
ROLE_NOTES = {
    "well_data_json": "Per-well data (track, design, feed zones, T/P surveys, transients, flows). DATA.",
    "well_data_dat": "Plain-column copies of well data. DATA, cross-check only.",
    "mapping": "Fault, lake and river polylines. DATA.",
    "well_groups": "Well list and groupings. DATA.",
    "geometry_dat": "MULgraph grid file. GEOMETRY ONLY (vertices, grid, layers, surface, wells).",
    "well_blocks_json": "Well-to-block intersections. GEOMETRY ONLY; authority for block assignment.",
    "reference_rock": "Rock-type dictionary of the reference model. REFERENCE ONLY.",
    "reference_bc": "Boundary conditions of the reference model. REFERENCE ONLY.",
    "reference_results": "Initial conditions / results of the reference model. REFERENCE ONLY.",
    "reference_other": "Model set-up, plotting or auxiliary file. REFERENCE ONLY.",
    "out_of_scope": "Not read for the natural-state analysis.",
}


def classify(path):
    rel = path.replace("\\", "/")
    low = rel.lower()
    ext = os.path.splitext(low)[1]
    if any(k in low for k in ("/futurescenarios/", "/production/", "productionhistory")):
        return "out_of_scope", "production / forecast set-up or history"
    if ext == ".json":
        try:
            d = json.load(open(path, encoding="utf-8"))
        except Exception:
            return "reference_other", "unparseable json"
        if isinstance(d, list):
            return "reference_other", "plotting definition (contains hard-coded paths from other machines)"
        if "Well_Track" in d:
            return "well_data_json", ROLE_NOTES["well_data_json"]
        vals = [v for v in d.values() if isinstance(v, dict)]
        if vals and all("Blocks" in v for v in vals):
            return "well_blocks_json", ROLE_NOTES["well_blocks_json"]
        if "Areas" in d:
            return "well_groups", ROLE_NOTES["well_groups"]
        if vals and all("Data" in v and "type" in v for v in vals):
            return "mapping", ROLE_NOTES["mapping"]
        if "lithology" in d:
            return "reference_rock", ROLE_NOTES["reference_rock"]
        if any(k in d for k in ("Heat", "Rain", "Atmosphere")):
            return "reference_bc", ROLE_NOTES["reference_bc"]
        if "setup" in d:
            return "reference_other", "model set-up file"
        if "Total" in d:
            return "out_of_scope", "field-total steam flow and power output (production)"
        return "reference_other", "other json"
    if ext == ".auxjson":
        return "reference_other", "auxiliary geometry scene file (copy of mapping with hard-coded source path)"
    if ext == ".incon":
        return "reference_results", "initial conditions of the reference natural-state model"
    if ext == ".dat":
        txt = open(path, encoding="utf-8", errors="ignore").read(200000)
        if re.search(r"^VERTICES", txt, re.M) and re.search(r"^LAYERS", txt, re.M):
            return "geometry_dat", ROLE_NOTES["geometry_dat"]
        if re.search(r"^ROCKS", txt, re.M):
            return "reference_other", ("reference natural-state simulator input (ROCKS = reference rock "
                                       "properties, GENER = reference heat/mass sources). REFERENCE ONLY")
        lines = [ln.split() for ln in txt.splitlines() if ln.strip()]
        if lines and all(len(x) == 2 for x in lines):
            return "well_data_dat", ROLE_NOTES["well_data_dat"]
        return "reference_other", "other dat"
    if ext in (".xlsx", ".xls"):
        return "out_of_scope", "production / reinjection / monitoring spreadsheets"
    return "out_of_scope", "system or unrelated file"


def build_inventory():
    files = []
    for dp, dns, fns in os.walk(BASE):
        dns.sort()
        for f in sorted(fns):
            files.append(os.path.join(dp, f).replace("\\", "/"))
    cls = {f: classify(f) for f in files}
    # collapse directories whose files all share one role (top-most such directory)
    dirs = sorted({os.path.dirname(f) for f in files})
    rows, covered = [], set()

    def files_under(d):
        return [f for f in files if f.startswith(d + "/")]

    for d in sorted(dirs, key=lambda s: s.count("/")):
        if d in covered or any(d.startswith(c + "/") for c in covered):
            continue
        fs = files_under(d)
        roles = {cls[f][0] for f in fs}
        if len(roles) == 1 and d.count("/") >= 3 and len(fs) > 1:
            role = roles.pop()
            note = cls[fs[0]][1]
            if role == "out_of_scope":
                note = "; ".join(sorted({cls[f][1] for f in fs}))
            rows.append((role, d, f"directory ({len(fs)} files): {note}"))
            covered.add(d)
    # also cover empty sub-directories implicitly through parent rows; add remaining files
    for f in files:
        if any(f.startswith(c + "/") for c in covered):
            continue
        rows.append((cls[f][0], f, cls[f][1]))
    inv = pd.DataFrame(rows, columns=["role", "path", "note"]).sort_values(["role", "path"])
    inv.to_csv(os.path.join(OUT, "input_inventory.csv"), index=False)
    log(f"Inventory: {len(files)} files under {BASE}/ in {len(inv)} rows")
    for r in inv.itertuples():
        log(f"  {r.role:18s} {r.path}")
    return inv


# ======================================================================================
# 2. Readers
# ======================================================================================
def lit(s):
    return ast.literal_eval(s) if isinstance(s, str) else s


def ddate(s):
    return datetime.strptime(s, "%d/%m/%Y")


def natkey(s):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def read_geometry(path):
    lines = open(path, encoding="utf-8").read().splitlines()
    sec, verts, grid, layers, surf, wells = None, {}, {}, [], {}, {}
    cur = None
    for ln in lines:
        head = ln.strip().split()[0] if ln.strip() else ""
        if head in ("VERTICES", "GRID", "CONNECTIONS", "LAYERS", "SURFA", "WELLS") and ln[:1] != " ":
            sec = head
            continue
        if not ln.strip():
            sec = None if sec != "GRID" else sec
            if sec == "GRID":
                continue
            continue
        if sec == "VERTICES":
            # fixed width: 3-char name, x (10), y (10) -- this section has no whitespace separators
            verts[ln[:3]] = (float(ln[3:13]), float(ln[13:23]))
        elif sec == "GRID":
            if len(ln.rstrip()) > 3 and ln[3:].strip():
                cur = ln[:3]
                grid[cur] = []
            else:
                grid[cur].append(ln[:3])
        elif sec == "LAYERS":
            p = ln.split()
            layers.append((int(p[0]), float(p[1]), float(p[2])))
        elif sec == "SURFA":
            p = ln.split()
            surf[p[0]] = float(p[1])
        elif sec == "WELLS":
            p = ln.split()
            wells.setdefault(p[0], []).append(tuple(float(x) for x in p[1:4]))
    lay = {}
    for i, b, c in layers:
        if i == 0:
            continue
        top = [x for x in layers if x[0] == i - 1][0][1]
        lay[i] = dict(top=top, bottom=b, centre=c)
    return dict(verts=verts, grid=grid, layers=lay, layer0=layers[0][1], surf=surf, wells=wells)


def col_polys(geo):
    return {c: MPath(np.array([geo["verts"][v] for v in vs])) for c, vs in geo["grid"].items()}


def find_col(polys, x, y):
    hits = [c for c, p in polys.items() if p.contains_point((x, y), radius=1e-6)]
    return hits[0] if hits else None


def col_top_layer(geo, col):
    s = geo["surf"].get(col, geo["layer0"])
    cands = [i for i, L in geo["layers"].items() if L["bottom"] < s - 1e-6]
    return min(cands)


def read_wells(json_dir, groups_path):
    order = None
    if groups_path:
        g = json.load(open(groups_path))
        for area in g.get("Areas", {}).values():
            for lst in area.values():
                order = list(lst)
    W = {}
    for f in sorted(os.listdir(json_dir)):
        if not f.endswith(".json"):
            continue
        d = json.load(open(os.path.join(json_dir, f)))
        name = d["Well_Name"]
        track = np.array([lit(r) for r in d["Well_Track"]["Data"]], dtype=float)
        design = []
        for r in d["Well_Design"]["Data"]:
            t, b, dia, slot, spec = lit(r)
            design.append(dict(top=float(t), bottom=float(b), dia=dia, slotted=bool(slot), spec=spec))
        fz = []
        for r in d["Feed_Zones"]["Data"]:
            t, b, q = lit(r)
            fz.append(dict(top=float(max(t, b)), bottom=float(min(t, b)), quality=str(q),
                           mid=0.5 * (float(t) + float(b))))
        surveys = {}
        for kind in ("Downhole_Temperature", "Downhole_Pressure"):
            for dt, s in d.get(kind, {}).items():
                surveys.setdefault(dt, {})[kind] = s
        sv = []
        for dt, s in surveys.items():
            T = s["Downhole_Temperature"]["Data"]
            P = s["Downhole_Pressure"]["Data"]
            zs = sorted({float(k) for k in T} | {float(k) for k in P}, reverse=True)
            Tm = {float(k): v for k, v in T.items()}
            Pm = {float(k): v for k, v in P.items()}
            df = pd.DataFrame(dict(z=zs, T=[Tm.get(z, np.nan) for z in zs], P=[Pm.get(z, np.nan) for z in zs]))
            meta = {k: {kk: vv for kk, vv in s[k].items() if kk not in ("Data", "Units")} for k in s}
            sv.append(dict(date=dt, df=df, meta=meta, units={k: s[k].get("Units") for k in s}))
        gauges = {}
        for el, s in d.get("Transient_Pressure", {}).items():
            gauges[float(el)] = dict(series=[(ddate(k), v) for k, v in s["Data"].items()],
                                     meta={k: v for k, v in s.items() if k not in ("Data", "Units")})
        flows = {}
        for k in ("Production_Mass_Flow", "Injection_Mass_Flow", "Transient_Enthalpy", "Injection_Enthalpy"):
            if k in d:
                flows[k] = dict(series=[(ddate(a), b) for a, b in d[k]["Data"].items()],
                                meta={kk: vv for kk, vv in d[k].items() if kk not in ("Data", "Units")})
        W[name] = dict(name=name, aut2=d.get("AUT2_Name"), track=track, design=design, fz=fz,
                       surveys=sv, gauges=gauges, flows=flows, raw=d, file=f,
                       fz_units=d["Feed_Zones"].get("Units"), track_units=d["Well_Track"].get("Units"),
                       design_source=d["Well_Design"].get("Source"))
    names = sorted(W, key=natkey)
    if order:
        names = [n for n in order if n in W] + [n for n in names if n not in order]
    return {n: W[n] for n in names}


def read_mapping(path):
    m = json.load(open(path))
    feats = {}
    for name, v in m.items():
        data = v["Data"]
        arr = np.array(data, dtype=object)
        segs = []
        if v.get("type") == "multiplot" or (len(data) and isinstance(data[0][0], list)):
            segs = [np.array(s, dtype=float) for s in data]
        else:
            segs = [np.array(data, dtype=float)]
        legend = str(v.get("legend", ""))
        kind = "fault" if ("fault" in name.lower() or legend.lower() == "fault") else "water"
        feats[name] = dict(segs=segs, kind=kind, fill=v.get("type") == "fill", legend=legend,
                           type=v.get("type"))
        del arr
    return feats


def seg_dist(px, py, seg):
    best = np.inf
    for i in range(len(seg) - 1):
        (x1, y1), (x2, y2) = seg[i], seg[i + 1]
        dx, dy = x2 - x1, y2 - y1
        L2 = dx * dx + dy * dy
        t = 0.0 if L2 == 0 else max(0.0, min(1.0, ((px - x1) * dx + (py - y1) * dy) / L2))
        best = min(best, math.hypot(px - x1 - t * dx, py - y1 - t * dy))
    return best


def feat_dist(px, py, f):
    if f["fill"]:
        if MPath(f["segs"][0]).contains_point((px, py)):
            return 0.0
    return min(seg_dist(px, py, s) for s in f["segs"])


# ======================================================================================
# helpers
# ======================================================================================
def interp_z(zq, z, v):
    """Linear interpolation in elevation; NaN outside [min z, max z] (never extrapolate)."""
    z = np.asarray(z, float)
    v = np.asarray(v, float)
    o = np.argsort(z)
    z, v = z[o], v[o]
    ok = np.isfinite(v)
    z, v = z[ok], v[ok]
    zq = np.atleast_1d(np.asarray(zq, float))
    out = np.interp(zq, z, v)
    out[(zq < z[0] - 1e-6) | (zq > z[-1] + 1e-6)] = np.nan
    return out


def track_xy(w, z):
    tr = w["track"]
    o = np.argsort(tr[:, 2])
    zz = tr[o, 2]
    return float(np.interp(z, zz, tr[o, 0])), float(np.interp(z, zz, tr[o, 1]))


def shoe(w):
    cs = [d for d in w["design"] if not d["slotted"]]
    return min(d["bottom"] for d in cs) if cs else w["track"][:, 2].max()


def rnd(x, n=2):
    return None if x is None or (isinstance(x, float) and not np.isfinite(x)) else round(float(x), n)


DERIVED = []


def F(v, nd=1, label=""):
    """Format a number for the report and record it in report_derived_numbers.csv."""
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "n/a"
    s = f"{v:.{nd}f}"
    DERIVED.append((label, s))
    return s


def md_table(df, fmts=None):
    fmts = fmts or {}
    cols = list(df.columns)
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, (float, np.floating)):
                if not np.isfinite(v):
                    cells.append("")
                else:
                    nd = fmts.get(c, 1)
                    cells.append(F(float(v), nd, f"table:{c}"))
            elif v is None:
                cells.append("")
            else:
                cells.append(str(v))
        out.append("| " + " | ".join(cells) + " |")
    return "\n".join(out)


# ======================================================================================
# MAIN
# ======================================================================================
def main():
    log("Natural-state analysis -- run from", os.path.basename(os.getcwd()) or ".")
    validate_if97()
    pd.DataFrame(PARAMS, columns=["parameter", "value", "meaning"]).to_csv(
        os.path.join(OUT, "analysis_parameters.csv"), index=False)

    inv = build_inventory()

    def role_paths(role):
        return list(inv.loc[inv.role == role, "path"])

    json_dir = role_paths("well_data_json")[0]
    dat_dir = role_paths("well_data_dat")[0] if role_paths("well_data_dat") else None
    map_path = role_paths("mapping")[0]
    grp = role_paths("well_groups")
    geo_path = role_paths("geometry_dat")[0]
    wb_path = role_paths("well_blocks_json")[0]
    W = read_wells(json_dir, grp[0] if grp else None)
    wells = list(W)
    log(f"Wells: {wells}")
    geo = read_geometry(geo_path)
    log(f"Geometry: {len(geo['verts'])} vertices, {len(geo['grid'])} columns, {len(geo['layers'])} layers, "
        f"{len(geo['surf'])} surface entries, {len(geo['wells'])} wells")
    WB = json.load(open(wb_path))
    MAP = read_mapping(map_path)
    polys = col_polys(geo)

    # ---------------- reference model summary ([reference model]) ----------------
    ref_rows = []
    for p in role_paths("reference_bc"):
        bc = json.load(open(p))
        for k, v in bc.items():
            if isinstance(v, dict):
                ref_rows.append(("[reference model] boundary condition", k,
                                 json.dumps(v.get("Default", ""))[:120], os.path.basename(p)))
    for p in role_paths("reference_other"):
        if p.endswith(".dat") and os.path.isfile(p):
            txt = open(p, encoding="utf-8", errors="ignore").read()
            rocks = re.search(r"^ROCKS\n(.*?)\n\s*\n", txt, re.S | re.M)
            if rocks:
                perms = []
                for ln in rocks.group(1).splitlines():
                    try:
                        perms.append(float(ln[30:40]))
                    except ValueError:
                        pass
                ref_rows.append(("[reference model] rock types", "count", len(perms), os.path.basename(p)))
                if perms:
                    ref_rows.append(("[reference model] rock types", "k_x range m2",
                                     f"{min(perms):.1e} to {max(perms):.1e}", os.path.basename(p)))
            gen = re.search(r"^GENER\n(.*?)(\n\s*\n|\Z)", txt, re.S | re.M)
            if gen:
                types = {}
                for ln in gen.group(1).splitlines():
                    toks = ln[10:].split()
                    t = next((x for x in toks if re.fullmatch(r"[A-Z][A-Z0-9]*", x)), None)
                    if t is None:
                        continue
                    i = toks.index(t)
                    try:
                        q = float(toks[i + 1])
                    except (IndexError, ValueError):
                        q = 0.0
                    try:
                        hx = float(toks[i + 2])
                    except (IndexError, ValueError):
                        hx = float("nan")
                    types.setdefault(t, []).append((ln[:5], q, hx))
                for t, lst in sorted(types.items()):
                    lays = sorted({int(b[3:5]) for b, _, _ in lst if b[3:5].strip().isdigit()})
                    ref_rows.append(("[reference model] GENER " + t, "count / total rate (kg/s or W)",
                                     f"{len(lst)} / {sum(q for _, q, _ in lst):.4g}", os.path.basename(p)))
                    ref_rows.append(("[reference model] GENER " + t, "layers of source blocks",
                                     f"{lays[0]} to {lays[-1]}" if lays else "", os.path.basename(p)))
                    if t != "HEAT":
                        hs_ = sorted({round(h) for _, _, h in lst if np.isfinite(h)})
                        ref_rows.append(("[reference model] GENER " + t, "injected enthalpy J/kg",
                                         ", ".join(f"{h:.4g}" for h in hs_[:5]), os.path.basename(p)))
                    qs = np.array([q for _, q, _ in lst])
                    ref_rows.append(("[reference model] GENER " + t, "rate min / max per source",
                                     f"{qs.min():.4g} / {qs.max():.4g}", os.path.basename(p)))
                    top = sorted(lst, key=lambda r: -r[1])[:10]
                    REF_TOP[t] = [b[:3] for b, _, _ in top]
    pd.DataFrame(ref_rows, columns=["item", "quantity", "value", "source"]).to_csv(
        os.path.join(OUT, "reference_model_summary.csv"), index=False)

    # ---------------- provenance ----------------
    prov = []
    for n, w in W.items():
        d = w["raw"]
        prov.append((n, "Well_Design", "Source", w["design_source"]))
        prov.append((n, "Well_Track", "Units", json.dumps(w["track_units"])))
        prov.append((n, "Feed_Zones", "Units", json.dumps(w["fz_units"])))
        prov.append((n, "Feed_Zones", "quality tags", ",".join(f["quality"] for f in w["fz"])))
        for s in w["surveys"]:
            for k, m in s["meta"].items():
                for kk, vv in m.items():
                    prov.append((n, f"{k} {s['date']}", kk, repr(vv)))
        for el, g in w["gauges"].items():
            for kk, vv in g["meta"].items():
                prov.append((n, f"Transient_Pressure {el}", kk, repr(vv)))
        for k, fl in w["flows"].items():
            for kk, vv in fl["meta"].items():
                prov.append((n, k, kk, repr(vv)))
        for k in d:
            if k not in ("Well_Name", "AUT2_Name", "Feed_Zones", "Well_Track", "Well_Design", "Downhole_Temperature",
                         "Downhole_Pressure", "Transient_Pressure", "Production_Mass_Flow", "Injection_Mass_Flow",
                         "Transient_Enthalpy", "Injection_Enthalpy") and not k.endswith("Generator_Name"):
                prov.append((n, k, "unhandled key", json.dumps(d[k])[:100]))
    prov = pd.DataFrame(prov, columns=["well", "item", "field", "value"])
    prov.to_csv(os.path.join(OUT, "provenance_labels.csv"), index=False)

    # ---------------- integrity checks ----------------
    integ = []

    def chk(well, check, ok, detail):
        integ.append((well, check, "OK" if ok else "ISSUE", detail))

    survey_dates = sorted({s["date"] for w in W.values() for s in w["surveys"]})
    for n, w in W.items():
        tr = w["track"]
        ztop, zbot = tr[:, 2].max(), tr[:, 2].min()
        if len(w["surveys"]) != 1:
            chk(n, "one survey", False, f"{len(w['surveys'])} surveys; the first is used")
        sdf = w["surveys"][0]["df"]
        chk(n, "T and P at same elevations", sdf[["T", "P"]].notna().all().all(),
            f"{int(sdf['T'].isna().sum())} missing T, {int(sdf['P'].isna().sum())} missing P")
        inside = sdf.z.between(zbot - 0.5, ztop + 0.5)
        chk(n, "survey points within track", inside.all(),
            f"survey {sdf.z.max():.1f} to {sdf.z.min():.1f} mRL; track {ztop:.1f} to {zbot:.1f} mRL")
        if sdf.z.max() < ztop - 1:
            chk(n, "survey starts at wellhead", False,
                f"first survey point {sdf.z.max():.1f} mRL is {ztop - sdf.z.max():.1f} m below the wellhead {ztop:.1f}")
        slotted = [d for d in w["design"] if d["slotted"]]
        for i, f in enumerate(w["fz"]):
            chk(n, f"feed zone {i + 1} within track", zbot - 0.5 <= f["bottom"] and f["top"] <= ztop + 0.5,
                f"{f['top']:.1f} to {f['bottom']:.1f}")
            ins = any(d["bottom"] - 0.5 <= f["bottom"] and f["top"] <= d["top"] + 0.5 for d in slotted)
            chk(n, f"feed zone {i + 1} within slotted interval", ins,
                f"{f['top']:.1f} to {f['bottom']:.1f}; slotted " +
                ",".join(f"{d['top']:.0f}/{d['bottom']:.0f}" for d in slotted))
        dtop = max(d["top"] for d in w["design"])
        dbot = min(d["bottom"] for d in w["design"])
        chk(n, "design top = wellhead", abs(dtop - ztop) < 1, f"design top {dtop:.1f}, wellhead {ztop:.1f}")
        chk(n, "design bottom = track bottom", abs(dbot - zbot) < 1, f"design bottom {dbot:.1f}, track bottom {zbot:.1f}")
        ds = sorted(w["design"], key=lambda d: -d["top"])
        for a, b in zip(ds, ds[1:]):
            chk(n, "strings overlap (liner hung inside casing)", b["top"] >= a["bottom"] - 0.5,
                f"{a['spec']} to {a['bottom']:.1f}; {b['spec']} from {b['top']:.1f}")
        gw = geo["wells"].get(w["aut2"])
        if gw:
            g = np.array(gw)
            same = len(g) == len(tr) and np.allclose(g, tr, atol=0.11)
            chk(n, "grid-file WELLS track = JSON track", same, f"{len(g)} vs {len(tr)} points")
        else:
            chk(n, "grid-file WELLS track = JSON track", False, f"{w['aut2']} not in grid file")
        if n not in WB:
            chk(n, "well in well-to-block file", False, "missing")
        else:
            fzj = list(WB[n]["Feed_Zones"].values())
            chk(n, "feed-zone count matches block file", len(fzj) == len(w["fz"]), f"{len(w['fz'])} vs {len(fzj)}")
            for i, (f, fj) in enumerate(zip(w["fz"], fzj)):
                ints = [tuple(sorted(v["Interval"])) for v in fj.values()]
                ok = any(abs(a - f["bottom"]) < 0.6 and abs(b - f["top"]) < 0.6 for a, b in ints) or \
                    (abs(min(a for a, _ in ints) - f["bottom"]) < 0.6 and abs(max(b for _, b in ints) - f["top"]) < 0.6)
                chk(n, f"feed zone {i + 1} interval = block file", ok, f"{f['top']}/{f['bottom']} vs {ints}")
            bl = WB[n]["Blocks"]
            chk(n, "block-file top = wellhead", abs(max(bl[0]["interval"]) - ztop) < 1,
                f"block-file top {max(bl[0]['interval']):.1f}, wellhead {ztop:.1f}")
            chk(n, "block-file bottom = track bottom", abs(min(bl[-1]["interval"]) - zbot) < 1,
                f"block-file bottom {min(bl[-1]['interval']):.1f}, track bottom {zbot:.1f}")

    # ---------------- dat cross-check ----------------
    xc = []
    irregular = []
    if dat_dir:
        for f in sorted(os.listdir(dat_dir)):
            m = re.match(r"^(Wel+)(\d+)_(.*)\.dat$", f)
            if not m:
                irregular.append(f + " (unmatched)")
                continue
            num = int(m.group(2))
            if m.group(1) != "Well":
                irregular.append(f)
            well = next((n for n in W if re.search(r"(\d+)$", n) and int(re.search(r"(\d+)$", n).group(1)) == num), None)
            arr = np.loadtxt(os.path.join(dat_dir, f), ndmin=2)
            kind = m.group(3)
            w = W[well]
            if kind in ("downhole_temperature", "downhole_pressure"):
                col = "T" if "temperature" in kind else "P"
                sdf = w["surveys"][0]["df"]
                js = dict(zip(sdf.z, sdf[col]))
                diffs = [abs(js.get(z, np.nan) - v) for v, z in arr]
                nmiss = sum(1 for _, z in arr if z not in js)
                xc.append((well, f, col, len(arr), len(sdf), nmiss, np.nanmax(diffs), "", ""))
            else:
                gm = re.match(r"^(-?[\d.]+)_transient_pressure$", kind)
                if gm:
                    series = w["gauges"][float(gm.group(1))]["series"]
                    label = f"transient P at {gm.group(1)}"
                elif kind == "production_mass_flow":
                    series = w["flows"]["Production_Mass_Flow"]["series"]
                    label = "production mass flow"
                elif kind == "transient_enthalpy":
                    series = w["flows"]["Transient_Enthalpy"]["series"]
                    label = "production enthalpy"
                elif kind == "injection_mass_flow":
                    series = w["flows"]["Injection_Mass_Flow"]["series"]
                    label = "injection mass flow"
                elif kind == "injection_enthalpy":
                    series = w["flows"]["Injection_Enthalpy"]["series"]
                    label = "injection enthalpy"
                else:
                    irregular.append(f + " (unknown kind)")
                    continue
                n_ = min(len(series), len(arr))
                refs = [series[i][0] - pd.Timedelta(days=arr[i, 0] * 365.25) for i in range(n_)]
                ref = pd.Series(refs).median()
                ref = pd.Timestamp(ref.year, ref.month, 1) if ref.day <= 15 else                     (pd.Timestamp(ref.year, ref.month, 1) + pd.offsets.MonthBegin(1))
                tj = np.array([(d - ref).days / 365.25 for d, _ in series[:n_]])
                vj = np.array([v for _, v in series[:n_]])
                sign = np.sign(np.nansum(arr[:n_, 1] * vj)) if np.any(vj) else 1
                dv = np.abs(np.abs(arr[:n_, 1]) - np.abs(vj))
                note = "dat value has opposite sign (generation convention: production negative)" if sign < 0 else ""
                xc.append((well, f, label, len(arr), len(series), abs(len(arr) - len(series)),
                           float(dv.max()), ref.strftime("%d/%m/%Y") + f" (max time misfit {np.max(np.abs(tj - arr[:n_, 0])):.3f} yr)", note))
    xc = pd.DataFrame(xc, columns=["well", "dat_file", "quantity", "n_dat", "n_json", "n_unmatched",
                                   "max_abs_diff", "reference_date", "note"])
    xc.to_csv(os.path.join(OUT, "dat_json_crosscheck.csv"), index=False)
    for r in xc.itertuples():
        chk(r.well, f"dat cross-check {r.quantity}", r.n_unmatched == 0 and r.max_abs_diff < 1e-6,
            f"{r.dat_file}: max diff {r.max_abs_diff:.3g}, unmatched {r.n_unmatched}")
    wells_with_dat = set(xc.well)
    for n in W:
        if n not in wells_with_dat:
            chk(n, "dat files present", False, "no plain-column files")
    log(f"dat cross-check: {len(xc)} files; irregular names: {irregular}")

    # ---------------- block assignment and own point-in-polygon check ----------------
    def block_for(n, z):
        bl = WB[n]["Blocks"]
        for b in bl:
            top, bot = max(b["interval"]), min(b["interval"])
            if bot < z <= top + 1e-6:
                return b["block"]
        for b in bl:
            top, bot = max(b["interval"]), min(b["interval"])
            if bot - 0.5 <= z <= top + 0.5:
                return b["block"]
        return None

    # ---------------- per-well survey tables ----------------
    first_prod_field = min((d for w in W.values() for d, v in w["flows"].get("Production_Mass_Flow", {}).get("series", [])
                            if v > 0), default=None)
    first_inj_field = min((d for w in W.values() for d, v in w["flows"].get("Injection_Mass_Flow", {}).get("series", [])
                           if v > 0), default=None)
    pip_rows = []
    rows = []
    for n, w in W.items():
        sdf = w["surveys"][0]["df"].copy()
        sdf["Tsat"] = [tsat_c(p) for p in sdf.P]
        sdf["margin"] = sdf.Tsat - sdf["T"]
        w["sdf"] = sdf
        w["shoe"] = shoe(w)
        w["zwh"] = w["track"][:, 2].max()
        w["x0"], w["y0"] = w["track"][0, 0], w["track"][0, 1]
        bl0 = WB[n]["Blocks"][0]["block"]
        w["wh_col"] = bl0[:3]
        w["surface"] = geo["surf"].get(w["wh_col"], geo["layer0"])
        own = find_col(polys, w["x0"], w["y0"])
        pip_rows.append((n, "wellhead", w["zwh"], w["wh_col"], own, own == w["wh_col"]))
        for z in sdf.z:
            x, y = track_xy(w, z)
            b = block_for(n, z)
            if b is not None:
                own = find_col(polys, x, y)
                pip_rows.append((n, "survey", z, b[:3], own, own == b[:3]))
    pip = pd.DataFrame(pip_rows, columns=["well", "point", "elevation_mRL", "column_block_file", "column_own_pip", "agree"])
    pip.to_csv(os.path.join(OUT, "column_assignment_check.csv"), index=False)
    log(f"Own point-in-polygon column vs block file: {int(pip.agree.sum())}/{len(pip)} agree")

    # ---------------- injectate temperature ----------------
    inj_T = {}
    for n, w in W.items():
        if "Injection_Enthalpy" in w["flows"]:
            hs = np.array([v for _, v in w["flows"]["Injection_Enthalpy"]["series"]])
            inj_T[n] = dict(h_mean=float(hs.mean()), T=T_from_h(float(hs.mean()), P_INJ_ASSUMED))
    for n, v in inj_T.items():
        log(f"{n}: mean injectate enthalpy {v['h_mean']:.1f} kJ/kg -> {v['T']:.1f} degC (IF97 Region 1 at {P_INJ_ASSUMED} bar)")

    # ---------------- isothermal sections ----------------
    iso_rows = []
    for n, w in W.items():
        sdf = w["sdf"]
        z, T = sdf.z.values, sdf["T"].values
        w["iso"] = []
        i = 0
        while i < len(z) - 1:
            j = i
            while j < len(z) - 1 and abs(T[j + 1] - T[j]) / (z[j] - z[j + 1]) * 1000 < ISO_GRAD:
                j += 1
            if j > i and z[i] - z[j] >= ISO_LEN:
                fzin = [k + 1 for k, f in enumerate(w["fz"]) if f["bottom"] <= z[i] and f["top"] >= z[j]]
                loc = ("open hole" if z[i] < w["shoe"] else "behind casing" if z[j] >= w["shoe"]
                       else "straddles casing shoe")
                w["iso"].append((z[i], z[j]))
                # internal wellbore flow needs an entry and an exit: a feed-zone pair bracketing the section
                Tm = float(np.mean(T[i:j + 1]))
                fzT = [float(interp_z(f["mid"], z, T)[0]) for f in w["fz"]]
                pairs = []
                for a, fa in enumerate(w["fz"]):
                    for b, fb in enumerate(w["fz"]):
                        if fa["mid"] <= fb["mid"]:
                            continue
                        ov = min(fa["mid"], z[i]) - max(fb["mid"], z[j])
                        if ov >= 0.5 * (z[i] - z[j]):   # at least half the section lies between the two feed zones
                            pairs.append((a + 1, b + 1))
                match = [k + 1 for k, t_ in enumerate(fzT) if abs(t_ - Tm) <= 5.0]
                pairs.sort(key=lambda ab: w["fz"][ab[0] - 1]["mid"] - w["fz"][ab[1] - 1]["mid"])
                if pairs:
                    verdict_if = (f"internal flow possible between feed zones {pairs[0][0]} and {pairs[0][1]}"
                                  + (f"; section T matches feed zone {','.join(map(str, match))}" if match else ""))
                else:
                    verdict_if = "no feed-zone pair brackets it: internal flow unlikely, probably a formation feature"
                iso_rows.append((n, z[i], z[j], z[i] - z[j], T[i:j + 1].min(), T[i:j + 1].max(), loc,
                                 ",".join(map(str, fzin)), w["shoe"], verdict_if))
            i = max(j, i + 1)
    iso = pd.DataFrame(iso_rows, columns=["well", "top_mRL", "bottom_mRL", "length_m", "Tmin_degC", "Tmax_degC",
                                          "location", "feedzones_inside", "shoe_mRL", "internal_flow_test"])
    iso.to_csv(os.path.join(OUT, "isothermal_sections.csv"), index=False, float_format="%.2f")

    # ---------------- natural-state timing test (4.1.2) ----------------
    pc_rows = []
    verdict = {}
    passed = []
    gauge_info = {}
    for n, w in W.items():
        sdf = w["sdf"]
        own_prod = next((d for d, v in w["flows"].get("Production_Mass_Flow", {}).get("series", []) if v > 0), None)
        own_inj = next((d for d, v in w["flows"].get("Injection_Mass_Flow", {}).get("series", []) if v > 0), None)
        own_first = min([d for d in (own_prod, own_inj) if d is not None], default=None)
        w["own_first"] = own_first
        w["role"] = "producer" if own_prod else "injector" if own_inj else "monitoring/unused"
        if w["gauges"]:
            el = sorted(w["gauges"])[0]
            ser = sorted(w["gauges"][el]["series"])
            d0, p0 = ser[0]
            ps = float(interp_z(el, sdf.z, sdf.P)[0])
            diff = ps - p0
            gauge_info[n] = dict(el=el, d0=d0, p0=p0, ps=ps, diff=diff, ser=ser)
            pc_rows.append(dict(well=n, gauge_elev_mRL=el, survey_P_at_gauge_bar=round(ps, 3),
                                first_gauge_date=d0.strftime("%d/%m/%Y"), first_gauge_P_bar=p0,
                                diff_survey_minus_first_bar=round(diff, 3),
                                gauge_start_vs_own_first_flow_days=(d0 - own_first).days if own_first else None,
                                gauge_start_vs_field_first_production_days=(d0 - first_prod_field).days if first_prod_field else None,
                                own_first_flow_date=own_first.strftime("%d/%m/%Y") if own_first else "",
                                field_first_production_date=first_prod_field.strftime("%d/%m/%Y") if first_prod_field else "",
                                second_gauge_date=ser[1][0].strftime("%d/%m/%Y") if len(ser) > 1 else "",
                                second_gauge_P_bar=ser[1][1] if len(ser) > 1 else None))
        else:
            pc_rows.append(dict(well=n, own_first_flow_date=own_first.strftime("%d/%m/%Y") if own_first else "",
                                field_first_production_date=first_prod_field.strftime("%d/%m/%Y") if first_prod_field else ""))
    # gauge verdicts
    for n, gi in gauge_info.items():
        w = W[n]
        pre = first_prod_field is None or gi["d0"] <= first_prod_field
        if pre and abs(gi["diff"]) <= GAUGE_TOL + 1e-9:
            verdict[n] = ("natural state", "high",
                          f"survey P at the gauge ({gi['ps']:.2f} bar) matches the first gauge reading "
                          f"({gi['p0']:.1f} bar, {gi['d0']:%d/%m/%Y}, the date production started) within "
                          f"the {GAUGE_TOL} bar reading precision; the next reading is "
                          f"{gi['ser'][1][1]:.1f} bar, so a survey made after production started could not match")
            passed.append(n)
    # drawdown of passed gauges for comparison
    for n, gi in gauge_info.items():
        if n in verdict:
            continue
        pre = first_prod_field is None or gi["d0"] <= first_prod_field
        if pre:
            later = gi["ser"][1][1] if len(gi["ser"]) > 1 else None
            if gi["diff"] > 0 and later is not None and gi["ps"] > later:
                verdict[n] = ("probably natural state", "medium",
                              f"survey P at the gauge ({gi['ps']:.2f} bar) is {gi['diff']:.2f} bar ABOVE the first "
                              f"gauge reading ({gi['p0']:.1f} bar on {gi['d0']:%d/%m/%Y}, the date production started) "
                              f"and far above the next reading ({later:.1f} bar); production lowers pressure, so the "
                              f"survey pre-dates the drawdown, but the mismatch exceeds the reading precision")
            else:
                verdict[n] = ("possibly disturbed", "low",
                              f"survey P at the gauge differs from the first reading by {gi['diff']:.2f} bar")
        else:
            # gauge starts after field production: compare with drawdown at passed gauges over the same period
            dd = []
            for m in passed:
                s = gauge_info[m]["ser"]
                pt = [p for d, p in s if d <= gi["d0"]]
                if pt:
                    dd.append(s[0][1] - pt[-1])
            gi["passed_drawdown"] = dd
            gi["deferred"] = True
            continue
    # pressure trend from passed wells' feed zones
    trend_pts = []
    for m in passed:
        w = W[m]
        for f in w["fz"]:
            trend_pts.append((f["mid"], float(interp_z(f["mid"], w["sdf"].z, w["sdf"].P)[0])))
    trend_pts = np.array(trend_pts)
    if len(trend_pts) >= 2:
        tcoef = np.polyfit(trend_pts[:, 0], trend_pts[:, 1], 1)
    else:
        tcoef = None
    trend_rows = []
    for n, w in W.items():
        for i, f in enumerate(w["fz"]):
            p = float(interp_z(f["mid"], w["sdf"].z, w["sdf"].P)[0])
            T = float(interp_z(f["mid"], w["sdf"].z, w["sdf"]["T"])[0])
            ptr = float(np.polyval(tcoef, f["mid"])) if tcoef is not None else np.nan
            inside = trend_pts[:, 0].min() - 1 <= f["mid"] <= trend_pts[:, 0].max() + 1
            trend_rows.append((n, i + 1, f["mid"], round(p, 2), round(T, 1), round(ptr, 2), round(p - ptr, 2),
                               "interpolated" if inside else "outside trend range (extrapolated trend)",
                               n in passed))
    trend = pd.DataFrame(trend_rows, columns=["well", "feedzone", "mid_mRL", "P_bar", "T_degC", "P_trend_bar",
                                              "residual_bar", "trend_use", "trend_well"])
    trend.to_csv(os.path.join(OUT, "pressure_trend_check.csv"), index=False)
    for n, gi in gauge_info.items():
        if not gi.get("deferred"):
            continue
        res = trend[trend.well == n]
        dd = gi["passed_drawdown"]
        rr = (res.residual_bar.min(), res.residual_bar.max())
        yrs = (gi["d0"] - first_prod_field).days / 365.25
        if gi["diff"] > 0 and rr[0] >= -TREND_TOL:
            verdict[n] = ("probably natural state", "medium",
                          f"gauge starts {gi['d0']:%d/%m/%Y}, {yrs:.1f} yr after field production began, so the gauge "
                          f"test cannot be applied; the survey is {gi['diff']:.1f} bar above the first gauge reading, "
                          f"less than the {min(dd):.1f}-{max(dd):.1f} bar drawdown at the passed gauges over the same "
                          f"period, but the feed-zone pressure lies {rr[0]:+.1f} bar from the pre-production trend of "
                          f"the passed wells, whereas a survey taken after production began would sit several bar "
                          f"below it")
        else:
            verdict[n] = ("possibly disturbed", "low",
                          f"gauge starts {yrs:.1f} yr after field production; survey-minus-gauge {gi['diff']:.1f} bar; "
                          f"feed-zone residual against the passed-well trend {rr[0]:+.1f} bar")
    for n, w in W.items():
        if n in verdict:
            continue
        res = trend[trend.well == n]
        rmin = res.residual_bar.min()
        if n in inj_T:
            Ti = inj_T[n]["T"]
            fzT = res.T_degC.values
            iso_near = [s for s in iso_rows if s[0] == n and abs(0.5 * (s[4] + s[5]) - Ti) <= INJ_ISO_TOL]
            hotter = np.nanmax(fzT) > Ti
            trend_txt = f"feed-zone P residuals vs the passed-well trend {res.residual_bar.min():+.1f} to {res.residual_bar.max():+.1f} bar"
            if iso_near:
                verdict[n] = ("possibly disturbed", "low",
                              f"injector; isothermal section near the injectate temperature ({Ti:.0f} degC); {trend_txt}")
            elif hotter:
                verdict[n] = ("probably natural state", "low",
                              f"injector with no isothermal section at the injectate temperature ({Ti:.0f} degC), but the "
                              f"formation at the feed zones ({np.nanmin(fzT):.0f}-{np.nanmax(fzT):.0f} degC) is hotter than "
                              f"the injectate, so injection cooling could be hidden; {trend_txt}")
            else:
                verdict[n] = ("probably natural state", "medium",
                              f"injector; profile is colder than the injectate ({Ti:.0f} degC) at every feed zone "
                              f"(max {np.nanmax(fzT):.0f} degC) with no isothermal section near the injectate "
                              f"temperature: evidence against an injection-affected survey; no gauge; {trend_txt}")
        else:
            if rmin >= -TREND_TOL:
                verdict[n] = ("probably natural state", "medium",
                              f"no gauge and no production/injection record; feed-zone pressures lie "
                              f"{res.residual_bar.min():+.1f} to {res.residual_bar.max():+.1f} bar from the trend of the "
                              f"gauge-passed wells (no drawdown signature)")
            else:
                verdict[n] = ("possibly disturbed", "low",
                              f"no gauge; feed-zone pressure {rmin:+.1f} bar below the passed-well trend")
    for r in pc_rows:
        r["verdict"], r["confidence"] = verdict[r["well"]][:2]
        r["reason"] = verdict[r["well"]][2]
        r["role_in_field"] = W[r["well"]]["role"]
        if r["well"] in gauge_info and "passed_drawdown" in gauge_info[r["well"]]:
            r["passed_gauge_drawdown_same_period_bar"] = ";".join(f"{x:.1f}" for x in gauge_info[r["well"]]["passed_drawdown"])
    pcols = ["well", "gauge_elev_mRL", "survey_P_at_gauge_bar", "first_gauge_date", "first_gauge_P_bar",
             "diff_survey_minus_first_bar", "verdict", "confidence"]
    pc = pd.DataFrame(pc_rows)
    pc = pc[pcols + [c for c in pc.columns if c not in pcols]]
    pc.to_csv(os.path.join(OUT, "pressure_consistency_check.csv"), index=False)
    for n in W:
        W[n]["verdict"], W[n]["conf"], W[n]["vreason"] = verdict[n]
        log(f"{n}: {verdict[n][0]} / {verdict[n][1]} -- {verdict[n][2]}")

    # ---------------- feed-zone quality mapping ----------------
    grade_order = ["VERY HIGH", "HIGH", "GOOD", "MEDIUM", "MODERATE", "FAIR", "LOW", "POOR", "VERY LOW"]
    tags = sorted({f["quality"] for w in W.values() for f in w["fz"]})
    graded = sorted([t for t in tags if t.upper() in grade_order], key=lambda t: grade_order.index(t.upper()))
    qmap = {}
    for t in tags:
        if t not in graded:
            qmap[t] = (0.5, "not a quality grade -> lowest factor")
        elif t == graded[0]:
            qmap[t] = (1.0, "highest grade in dataset")
        elif t == graded[-1]:
            qmap[t] = (0.5, "lowest grade in dataset")
        else:
            qmap[t] = (0.8, "intermediate grade")
    log(f"Feed-zone quality tags: {qmap}")

    # ---------------- calibration targets ----------------
    trows = []
    for n, w in W.items():
        sdf = w["sdf"]
        conf = w["conf"]
        cf = CONF_F[conf]
        top_layer_by_col = {}
        for z, T, P, ts, mg in sdf[["z", "T", "P", "Tsat", "margin"]].itertuples(index=False):
            x, y = track_xy(w, z)
            b = block_for(n, z)
            above = z > w["surface"] + 0.5
            col = b[:3] if b else w["wh_col"]
            lay = int(b[3:]) if b else None
            if col not in top_layer_by_col:
                top_layer_by_col[col] = col_top_layer(geo, col)
            infz = [k for k, f in enumerate(w["fz"]) if f["bottom"] - 1e-6 <= z <= f["top"] + 1e-6]
            behind = z >= w["shoe"]
            in_iso = any(a >= z >= c for a, c in w["iso"])
            if above:
                reason = "above_surface"
            elif w["zwh"] - z <= WH_DIST or (lay is not None and lay == top_layer_by_col[col]):
                reason = "wellhead_or_top_block"
            elif infz:
                reason = "feed_zone"
            elif behind:
                reason = "behind_casing"
            elif in_iso:
                reason = "open_hole_isothermal"
            else:
                reason = "open_hole"
            wt = LOC_W[reason] * cf
            wp = qmap[w["fz"][infz[0]]["quality"]][0] * cf if (infz and not above) else 0.0
            trows.append(dict(well=n, x=round(x, 1), y=round(y, 1), elevation_mRL=z, T_degC=T, P_bar=P,
                              Tsat_degC=rnd(ts, 3), boiling_margin_degC=rnd(ts - T, 3) if np.isfinite(ts) else None,
                              column=col, layer=lay, block=b, in_feedzone=bool(infz), behind_casing=bool(behind),
                              above_model_surface=bool(above), ns_confidence=conf, weight=round(wt, 4),
                              weight_reason=f"{reason} ({LOC_W[reason]}) x {conf} ({cf})",
                              weight_P=round(wp, 4),
                              feedzone=infz[0] + 1 if infz else None, in_isothermal_section=in_iso,
                              low_margin_flag=bool(np.isfinite(ts) and ts - T < MARGIN_FLAG and P >= P_LIQ)))
    tgt = pd.DataFrame(trows)
    tgt["layer"] = tgt["layer"].astype("Int64")
    tgt["feedzone"] = tgt["feedzone"].astype("Int64")
    tgt.to_csv(os.path.join(OUT, "NS_calibration_targets.csv"), index=False)
    log(f"Calibration targets: {len(tgt)} survey points")

    # weights by well summary
    wbw = tgt.groupby("well", sort=False).agg(n_points=("weight", "size"), sum_weight_T=("weight", "sum"),
                                              n_P_targets=("weight_P", lambda s: int((s > 0).sum())),
                                              ns_confidence=("ns_confidence", "first")).reset_index()
    wbw.to_csv(os.path.join(OUT, "weights_by_well.csv"), index=False, float_format="%.3f")

    # boiling-margin flags
    bm = tgt[tgt.low_margin_flag][["well", "elevation_mRL", "T_degC", "P_bar", "Tsat_degC", "boiling_margin_degC",
                                   "behind_casing", "in_feedzone"]]
    bm.to_csv(os.path.join(OUT, "boiling_margin_flags.csv"), index=False)

    # ---------------- feed zones ----------------
    frows = []
    for n, w in W.items():
        fzj = list(WB[n]["Feed_Zones"].values())
        for i, f in enumerate(w["fz"]):
            T = float(interp_z(f["mid"], w["sdf"].z, w["sdf"]["T"])[0])
            P = float(interp_z(f["mid"], w["sdf"].z, w["sdf"].P)[0])
            ts = tsat_c(P)
            blk = max(fzj[i].items(), key=lambda kv: kv[1].get("Fraction", 0))[0] if i < len(fzj) else None
            allb = ";".join(f"{k}:{v.get('Fraction', '')}" for k, v in fzj[i].items()) if i < len(fzj) else ""
            cf = CONF_F[w["conf"]]
            frows.append(dict(well=n, feedzone=i + 1, top_mRL=f["top"], bottom_mRL=f["bottom"], mid_mRL=f["mid"],
                              quality=f["quality"], T_degC=round(T, 2), P_bar=round(P, 2), Tsat_degC=rnd(ts, 3),
                              boiling_margin_degC=rnd(ts - T, 3), block=blk, ns_confidence=w["conf"],
                              weight_T=round(LOC_W["feed_zone"] * cf, 4),
                              weight_P=round(qmap[f["quality"]][0] * cf, 4),
                              quality_factor=qmap[f["quality"]][0], blocks_and_fractions=allb,
                              n_survey_points_inside=int(((w["sdf"].z <= f["top"]) & (w["sdf"].z >= f["bottom"])).sum()),
                              behind_casing=f["top"] > w["shoe"]))
    fzdf = pd.DataFrame(frows)
    # pressure excess over a cold (20 degC, IF97) water column standing from the wellhead at 1 atm
    fzdf["P_cold_column_from_wellhead_bar"] = [
        round(integrate_P(W[r.well]["zwh"], r.mid_mRL, P_ATM, [W[r.well]["zwh"], r.mid_mRL], [COLD_COL_T, COLD_COL_T]), 2)
        for r in fzdf.itertuples()]
    fzdf["P_excess_over_cold_column_bar"] = (fzdf.P_bar - fzdf.P_cold_column_from_wellhead_bar).round(2)
    fzdf["P_cold_column_from_model_surface_bar"] = [
        round(integrate_P(W[r.well]["surface"], r.mid_mRL, P_ATM, [W[r.well]["surface"], r.mid_mRL], [COLD_COL_T, COLD_COL_T]), 2)
        for r in fzdf.itertuples()]
    fzdf["P_excess_over_cold_column_from_surface_bar"] = (fzdf.P_bar - fzdf.P_cold_column_from_model_surface_bar).round(2)
    fzdf.to_csv(os.path.join(OUT, "feedzones.csv"), index=False)

    # ---------------- layer centres ----------------
    lrows = []
    for n, w in W.items():
        sdf = w["sdf"]
        zmin, zmax = sdf.z.min(), sdf.z.max()
        for i, L in sorted(geo["layers"].items()):
            c = L["centre"]
            if not (zmin - 1e-6 <= c <= zmax + 1e-6):
                continue
            T = float(interp_z(c, sdf.z, sdf["T"])[0])
            P = float(interp_z(c, sdf.z, sdf.P)[0])
            ts = tsat_c(P)
            b = block_for(n, c)
            if c > w["surface"] + 0.5:
                b = None
            lrows.append(dict(well=n, layer=i, layer_centre_mRL=c, block=b, T_degC=round(T, 2), P_bar=round(P, 2),
                              Tsat_degC=rnd(ts, 3), boiling_margin_degC=rnd(ts - T, 3)))
    lc = pd.DataFrame(lrows)
    lc.to_csv(os.path.join(OUT, "T_P_at_layer_centres.csv"), index=False)

    # ---------------- common elevations ----------------
    crow = []
    for n, w in W.items():
        sdf = w["sdf"]
        for z in COMMON_ELEVS:
            T = float(interp_z(z, sdf.z, sdf["T"])[0])
            P = float(interp_z(z, sdf.z, sdf.P)[0])
            if np.isfinite(T):
                x, y = track_xy(w, z)
                crow.append(dict(well=n, elevation_mRL=z, x=round(x, 1), y=round(y, 1), T_degC=round(T, 2),
                                 P_bar=round(P, 2), Tsat_degC=rnd(tsat_c(P), 3),
                                 boiling_margin_degC=rnd(tsat_c(P) - T, 3)))
    ce = pd.DataFrame(crow)
    ce.to_csv(os.path.join(OUT, "TP_at_common_elevations.csv"), index=False)
    ceT = ce.pivot(index="well", columns="elevation_mRL", values="T_degC").reindex(wells)
    ceP = ce.pivot(index="well", columns="elevation_mRL", values="P_bar").reindex(wells)
    ceT = ceT[sorted(ceT.columns, reverse=True)]
    ceP = ceP[sorted(ceP.columns, reverse=True)]
    ceT.to_csv(os.path.join(OUT, "temperature_at_common_elevations.csv"))
    ceP.to_csv(os.path.join(OUT, "pressure_at_common_elevations.csv"))

    # ---------------- hydrostatic check ----------------
    hrows = []
    skipped = []
    for n, w in W.items():
        sdf = w["sdf"]
        z, T, P, mg = sdf.z.values, sdf["T"].values, sdf.P.values, sdf.margin.values
        ok = (mg >= MARGIN_LIQ) & (P >= P_LIQ) & np.isfinite(mg)
        runs, cur = [], []
        for i in range(len(z)):
            if ok[i]:
                cur.append(i)
            elif cur:
                runs.append(cur)
                cur = []
        if cur:
            runs.append(cur)
        for run in runs:
            for band, sel in (("cased", [i for i in run if z[i] >= w["shoe"]]),
                              ("open", [i for i in run if z[i] < w["shoe"]])):
                if len(sel) < 2:
                    continue
                zt, zb = z[sel[0]], z[sel[-1]]
                if zt - zb < WIN_MIN:
                    skipped.append((n, band, zt, zb, "band interval shorter than minimum window"))
                    continue
                cands = []
                for f in w["fz"]:
                    for ze in (f["top"], f["bottom"]):
                        k = min(sel, key=lambda i: abs(z[i] - ze))
                        if k not in cands and k not in (sel[0], sel[-1]):
                            cands.append(k)
                cands = sorted(cands)
                start = sel[0]
                cuts = [start]
                for c in cands:
                    if z[start] - z[c] >= WIN_MIN and z[c] - zb >= WIN_MIN:
                        cuts.append(c)
                        start = c
                cuts.append(sel[-1])
                for a, b in zip(cuts, cuts[1:]):
                    L = z[a] - z[b]
                    dPm = P[b] - P[a]
                    dPi = integrate_P(z[a], z[b], P[a], z, T) - P[a]
                    hrows.append(dict(well=n, z_top_mRL=z[a], z_bottom_mRL=z[b], length_m=L, band=band,
                                      dP_measured_bar=round(dPm, 3), dP_IF97_bar=round(dPi, 3),
                                      ratio_measured_to_IF97=round(dPm / dPi, 4),
                                      rho_implied_kgm3=round(dPm * 1e5 / (G * L), 1),
                                      rho_IF97_mean_kgm3=round(dPi * 1e5 / (G * L), 1),
                                      ratio_rounding_uncertainty=round(0.1 / dPi, 4),
                                      T_top_degC=T[a], T_bottom_degC=T[b],
                                      shoe_mRL=w["shoe"], well_Tmax_degC=float(np.nanmax(T))))
    hy = pd.DataFrame(hrows)
    cold = [n for n, w in W.items() if np.nanmax(w["sdf"]["T"]) < COLD_TMAX]
    cold_med = {b: float(hy[hy.well.isin(cold) & (hy.band == b)].ratio_measured_to_IF97.median())
                for b in ("cased", "open")}
    hy["cold_well_median_ratio_band"] = hy.band.map(cold_med).round(4)
    hy["excess_over_cold_component"] = (hy.ratio_measured_to_IF97 - hy.cold_well_median_ratio_band).round(4)
    hy.to_csv(os.path.join(OUT, "hydrostatic_check.csv"), index=False)
    hs = []
    for n in wells + ["ALL"]:
        for b in ("cased", "open", "all"):
            s = hy if n == "ALL" else hy[hy.well == n]
            if b != "all":
                s = s[s.band == b]
            if len(s) == 0:
                continue
            hs.append(dict(well=n, band=b, n_windows=len(s), total_length_m=s.length_m.sum(),
                           median_ratio=round(s.ratio_measured_to_IF97.median(), 4),
                           min_ratio=round(s.ratio_measured_to_IF97.min(), 4),
                           max_ratio=round(s.ratio_measured_to_IF97.max(), 4),
                           cold_well_component=round(cold_med.get(b, np.nan), 4) if b != "all" else None,
                           excess_over_cold_component=round(s.ratio_measured_to_IF97.median() - cold_med[b], 4)
                           if b != "all" and np.isfinite(cold_med.get(b, np.nan)) else None,
                           max_rounding_uncertainty=round(s.ratio_rounding_uncertainty.max(), 4)))
    hsum = pd.DataFrame(hs)
    hsum.to_csv(os.path.join(OUT, "hydrostatic_summary.csv"), index=False)
    pd.DataFrame(skipped, columns=["well", "band", "z_top_mRL", "z_bottom_mRL", "reason"]).to_csv(
        os.path.join(OUT, "hydrostatic_skipped.csv"), index=False)
    log(f"Hydrostatic: {len(hy)} windows; cold wells {cold}; cold-well median ratio {cold_med}")

    # ---------------- well summary ----------------
    srows = []
    for n, w in W.items():
        sdf = w["sdf"]
        z, T, P, mg = sdf.z.values, sdf["T"].values, sdf.P.values, sdf.margin.values
        imax = int(np.nanargmax(T))
        # inversion: largest drop below a running maximum
        best = (0.0, None, None)
        rmax, rz = -np.inf, None
        for i in range(len(z)):
            if T[i] > rmax:
                rmax, rz = T[i], z[i]
            if rmax - T[i] > best[0]:
                best = (rmax - T[i], rz, z[i])
        inv_sz, inv_top, inv_base = best if best[0] >= INV_MIN else (None, None, None)
        # shallow gradient over top 300 m of the liquid-filled survey
        liq = np.where(P >= P_LIQ)[0]
        zl0 = z[liq[0]]
        sel = [i for i in liq if zl0 - z[i] <= WIN_MIN + 1e-6]
        sg = -np.polyfit(z[sel], T[sel], 1)[0] * 1000 if len(sel) >= 2 else None
        # deep gradient below inversion base and deepest feed zone
        zlim = min([f["bottom"] for f in w["fz"]] + ([inv_base] if inv_base is not None else []))
        dsel = [i for i in range(len(z)) if z[i] <= zlim + 1e-6]
        dg = None
        dspan = z[dsel[0]] - z[dsel[-1]] if len(dsel) >= 2 else 0.0
        if len(dsel) >= 2 and zlim - z[-1] >= WIN_MIN:
            dg = -np.polyfit(z[dsel], T[dsel], 1)[0] * 1000
        # isotherm crossings

        def crossing(Tc):
            for i in range(len(z) - 1):
                if T[i] < Tc <= T[i + 1]:
                    return z[i] + (Tc - T[i]) / (T[i + 1] - T[i]) * (z[i + 1] - z[i])
            return None

        liqm = (P >= P_LIQ) & np.isfinite(mg)
        im = int(np.nanargmin(np.where(liqm, mg, np.inf)))
        # water level: linear fit of the shallowest liquid points (within 100 m), extrapolated to P_atm
        wsel = [i for i in liq if zl0 - z[i] <= 100 + 1e-6][:5]
        if len(wsel) < 2:
            wsel = list(liq[:2])
        pf = np.polyfit(z[wsel], P[wsel], 1)
        wl = (P_ATM - pf[1]) / pf[0]
        # pivot: dominant feed zone = highest quality grade, then thickest, then shallowest
        rank = {t: (graded.index(t) if t in graded else 99) for t in tags}
        dom = sorted(range(len(w["fz"])), key=lambda k: (rank[w["fz"][k]["quality"]],
                                                         -(w["fz"][k]["top"] - w["fz"][k]["bottom"]),
                                                         -w["fz"][k]["top"]))[0]
        pz = w["fz"][dom]["mid"]
        pP = float(interp_z(pz, z, P)[0])
        # cap: from wellhead down to the first point below which the local gradient < CAP_GRAD
        lg = [(T[i + 1] - T[i]) / (z[i] - z[i + 1]) * 1000 for i in range(len(z) - 1)]
        cap_base = None
        if sg is not None and sg >= CAP_GRAD:
            for i in range(liq[0], len(z) - 1):
                if lg[i] < CAP_GRAD and z[i] < zl0 - 50:
                    cap_base = z[i]
                    break
        Tmax = float(T[imax])
        mean_grad = (T[-1] - T[0]) / (z[0] - z[-1]) * 1000
        if mg[im] < 2.0:
            regime = "two-phase"
        elif inv_sz is not None and inv_sz >= 10:
            regime = "outflow"
        elif Tmax >= 200:
            regime = "upflow"
        elif mean_grad < 30:
            regime = "recharge/downflow"
        else:
            regime = "conductive"
        notes = []
        if regime == "upflow" and mg[im] < MARGIN_FLAG:
            notes.append(f"near-boiling (min margin {mg[im]:.1f} degC)")
        if regime == "outflow" and T[-1] > Tmax - 20 and T[-1] - np.nanmin(T[imax:]) >= 5:
            notes.append(f"reheats below the inversion (T at TD {T[-1]:.0f} degC, within 20 degC of Tmax)")
        if regime == "outflow":
            notes.append(f"outflow tongue centred near {z[imax]:.0f} mRL")
        if w["iso"]:
            notes.append("isothermal section(s): " + ", ".join(f"{a:.0f} to {b:.0f}" for a, b in w["iso"]))
        heat = K_ROCK * dg / 1000 if dg is not None and dg > 0 else None
        srows.append(dict(well=n, Tmax_degC=Tmax, Tmax_elev_mRL=z[imax], T_at_TD_degC=T[-1], TD_mRL=z[-1],
                          inversion_degC=rnd(inv_sz, 2), inversion_top_mRL=inv_top, inversion_base_mRL=inv_base,
                          shallow_gradient_degC_per_km=rnd(sg, 1), deep_gradient_degC_per_km=rnd(dg, 1),
                          z_100C_mRL=rnd(crossing(100.0), 1), z_200C_mRL=rnd(crossing(200.0), 1),
                          min_boiling_margin_degC=round(float(mg[im]), 2), min_boiling_margin_elev_mRL=z[im],
                          water_level_mRL=round(wl, 1), pivot_elev_mRL=pz, pivot_P_bar=round(pP, 2),
                          regime=regime, ns_verdict=w["verdict"], ns_confidence=w["conf"],
                          regime_note="; ".join(notes), wellhead_mRL=w["zwh"], shoe_mRL=w["shoe"],
                          dominant_feedzone=dom + 1, dominant_feedzone_quality=w["fz"][dom]["quality"],
                          mean_gradient_degC_per_km=round(mean_grad, 1),
                          cap_base_mRL=cap_base, cap_thickness_m=rnd(w["zwh"] - cap_base, 0) if cap_base is not None else None,
                          deep_conductive_heat_flux_W_m2=rnd(heat, 3), role_in_field=w["role"],
                          deep_interval_top_mRL=zlim, deep_gradient_points_span_m=dspan,
                          water_level_above_wellhead=bool(wl > w["zwh"])))
        w["summary"] = srows[-1]
    ws = pd.DataFrame(srows)
    ws.to_csv(os.path.join(OUT, "well_summary.csv"), index=False)

    # ---------------- geometry check ----------------
    grows = []
    feats = list(MAP)
    for n, w in W.items():
        dev = float(np.max(np.hypot(w["track"][:, 0] - w["x0"], w["track"][:, 1] - w["y0"])))
        d = {f: feat_dist(w["x0"], w["y0"], MAP[f]) for f in feats}
        faults = [f for f in feats if MAP[f]["kind"] == "fault"]
        water = [f for f in feats if MAP[f]["kind"] == "water"]
        nf = min(faults, key=lambda f: d[f]) if faults else None
        nw = min(water, key=lambda f: d[f]) if water else None
        row = dict(well=n, wellhead_mRL=w["zwh"], model_surface_mRL=w["surface"],
                   wellhead_minus_surface_m=round(w["zwh"] - w["surface"], 1),
                   max_horizontal_deviation_m=round(dev, 1), nearest_fault=nf,
                   dist_nearest_fault_m=round(d[nf], 1) if nf else None, nearest_water_body=nw,
                   dist_nearest_water_body_m=round(d[nw], 1) if nw else None)
        for f in feats:
            row[f"dist_{re.sub(r'[^A-Za-z0-9]+', '_', f)}_m"] = round(d[f], 1)
        row["wellhead_column"] = w["wh_col"]
        row["column_top_layer"] = col_top_layer(geo, w["wh_col"])
        grows.append(row)
    gc = pd.DataFrame(grows)
    gc.to_csv(os.path.join(OUT, "well_geometry_check.csv"), index=False)
    for r in gc.itertuples():
        chk(r.well, "wellhead vs model surface", abs(r.wellhead_minus_surface_m) <= 10,
            f"wellhead {r.wellhead_mRL:.1f}, model surface {r.model_surface_mRL:.1f} mRL")
    ig = pd.DataFrame(integ, columns=["well", "check", "status", "detail"])
    ig.to_csv(os.path.join(OUT, "integrity_checks.csv"), index=False)
    log(f"Integrity checks: {len(ig)} checks, {int((ig.status == 'ISSUE').sum())} issues")
    for r in ig[ig.status == "ISSUE"].itertuples():
        log(f"  ISSUE {r.well}: {r.check} -- {r.detail}")

    # ======================================================================================
    # FIGURES
    # ======================================================================================
    colors = plt.cm.tab10(np.linspace(0, 1, 10))
    wcol = {n: colors[i % 10] for i, n in enumerate(wells)}
    xs = np.array([W[n]["x0"] for n in wells])
    ys = np.array([W[n]["y0"] for n in wells])
    pad = 1500
    xlim = (xs.min() - pad, xs.max() + pad)
    ylim = (ys.min() - pad, ys.max() + pad)

    def draw_base(ax, grid=True, labels=True, lim=True):
        if grid:
            for c, p in polys.items():
                v = p.vertices
                ax.plot(np.r_[v[:, 0], v[0, 0]], np.r_[v[:, 1], v[0, 1]], color="0.85", lw=0.4, zorder=1)
        for name, f in MAP.items():
            for s in f["segs"]:
                if f["kind"] == "fault":
                    ax.plot(s[:, 0], s[:, 1], "k-", lw=1.6, zorder=3)
                elif f["fill"]:
                    ax.fill(s[:, 0], s[:, 1], color="lightskyblue", alpha=0.8, zorder=2)
                else:
                    ax.plot(s[:, 0], s[:, 1], color="deepskyblue", lw=1.2, zorder=2)
        for name, f in MAP.items():
            if f["kind"] == "fault":
                s = f["segs"][0]
                # label at the point of the fault nearest to the map centre
                cx, cy = np.mean(xlim), np.mean(ylim)
                t = np.linspace(0, 1, 200)
                px = s[0, 0] + t * (s[-1, 0] - s[0, 0])
                py = s[0, 1] + t * (s[-1, 1] - s[0, 1])
                inb = (px > xlim[0]) & (px < xlim[1]) & (py > ylim[0]) & (py < ylim[1])
                if inb.any():
                    k = np.where(inb)[0][int(0.85 * inb.sum()) - 1]
                    ax.text(px[k], py[k], name, fontsize=7, color="k", rotation=0, zorder=6,
                            bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.5))
        for n in wells:
            w = W[n]
            ax.plot(w["track"][:, 0], w["track"][:, 1], "-", color="0.3", lw=1, zorder=4)
            ax.plot(w["x0"], w["y0"], "o", ms=5, mfc="white", mec="k", zorder=5)
            if labels:
                ax.annotate(n, (w["x0"], w["y0"]), xytext=(4, 4), textcoords="offset points", fontsize=7, zorder=6)
        if lim:
            ax.set_xlim(*xlim)
            ax.set_ylim(*ylim)
        ax.set_aspect("equal")
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.tick_params(labelsize=7)

    # well map (full grid + zoom)
    fig, axs = plt.subplots(1, 2, figsize=(15, 7.5))
    draw_base(axs[0], lim=False)
    allv = np.array(list(geo["verts"].values()))
    axs[0].set_xlim(allv[:, 0].min(), allv[:, 0].max())
    axs[0].set_ylim(allv[:, 1].min(), allv[:, 1].max())
    axs[0].plot([xlim[0], xlim[1], xlim[1], xlim[0], xlim[0]], [ylim[0], ylim[0], ylim[1], ylim[1], ylim[0]], "r--", lw=1)
    axs[0].set_title("Model grid, mapped features and wells (red box = zoom)", fontsize=10)
    draw_base(axs[1])
    for n in wells:
        w = W[n]
        if w["summary"]["max_horizontal_deviation_m"] > 50 if "max_horizontal_deviation_m" in w["summary"] else False:
            pass
    axs[1].set_title("Well field (deviated tracks shown as lines)", fontsize=10)
    for ax in axs:
        ax.set_xlabel("Easting (m)", fontsize=8)
        ax.set_ylabel("Northing (m)", fontsize=8)
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    axs[1].legend(handles=[Line2D([], [], color="k", lw=1.6, label="fault"),
                           Patch(color="lightskyblue", label="lake"),
                           Line2D([], [], color="deepskyblue", label="river"),
                           Line2D([], [], color="0.85", label="model columns"),
                           Line2D([], [], marker="o", mfc="white", mec="k", ls="", label="wellhead")],
                  fontsize=7, loc="lower left")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "well_map.png"), dpi=130)
    plt.close(fig)

    # profiles
    for n, w in W.items():
        sdf = w["sdf"]
        fig, axs = plt.subplots(1, 3, figsize=(13, 7.5), sharey=True, gridspec_kw=dict(width_ratios=[3, 3, 2]))
        a, b, c = axs
        a.plot(sdf["T"], sdf.z, "o-", color="tab:red", ms=3, label="measured T")
        ok = sdf.P < P_CRIT
        a.plot(sdf.Tsat[ok], sdf.z[ok], "--", color="tab:purple", label="Tsat(measured P), IF97")
        if n in inj_T:
            a.axvline(inj_T[n]["T"], color="tab:cyan", ls=":", label=f"injectate T ({inj_T[n]['T']:.0f} °C)")
        b.plot(sdf.P, sdf.z, "o-", color="tab:blue", ms=3, label="measured P")
        for (zg, gi) in [(gauge_info[n]["el"], gauge_info[n])] if n in gauge_info else []:
            b.plot(gi["p0"], zg, "k*", ms=12, label=f"first gauge reading {gi['d0']:%Y}")
        c.plot(sdf.margin, sdf.z, "o-", color="tab:green", ms=3, label="boiling margin")
        c.axvline(MARGIN_FLAG, color="tab:orange", ls=":", label=f"{MARGIN_FLAG:.0f} °C flag")
        for ax in axs:
            for f in w["fz"]:
                ax.axhspan(f["bottom"], f["top"], color="gold", alpha=0.5)
            ax.axhline(w["shoe"], color="k", lw=1.2, ls="-.")
            ax.axhline(w["zwh"], color="saddlebrown", lw=1)
            ax.axhline(w["surface"], color="green", lw=0.8, ls=":")
            for a_, b_ in w["iso"]:
                ax.axhspan(b_, a_, color="0.6", alpha=0.15)
            ax.grid(alpha=0.3)
            ax.tick_params(labelsize=8)
        xl = a.get_xlim()
        a.text(xl[0], w["shoe"], f" casing shoe {w['shoe']:.0f}", fontsize=7, va="bottom")
        a.text(xl[0], w["zwh"], f" wellhead {w['zwh']:.0f}", fontsize=7, va="bottom", color="saddlebrown")
        for k, f in enumerate(w["fz"]):
            b.text(b.get_xlim()[1], f["mid"], f"FZ{k + 1} {f['quality']} ", fontsize=7, ha="right", va="center")
        a.set_xlabel("Temperature (°C)")
        b.set_xlabel("Pressure (bar)")
        c.set_xlabel("Tsat − T (°C)")
        a.set_ylabel("Elevation (mRL)")
        handles = [Patch(color="gold", alpha=0.5, label="feed zone"), Line2D([], [], color="k", ls="-.", label="casing shoe"),
                   Line2D([], [], color="saddlebrown", label="wellhead"),
                   Line2D([], [], color="green", ls=":", label="model surface")]
        if w["iso"]:
            handles.append(Patch(color="0.6", alpha=0.3, label="isothermal ≥300 m"))
        a.legend(handles=a.get_legend_handles_labels()[0] + handles, fontsize=7, loc="best", framealpha=0.85)
        b.legend(fontsize=7, loc="upper right")
        c.legend(fontsize=7, loc="lower right")
        c.set_xlim(left=min(-5, np.nanmin(sdf.margin) - 5))
        fig.suptitle(f"{n} — survey labelled {w['surveys'][0]['date']} — NS verdict: {w['verdict']} ({w['conf']}) — "
                     f"regime: {w['summary']['regime']}", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, f"profile_{n}.png"), dpi=110)
        plt.close(fig)

    # P-T state
    fig, ax = plt.subplots(figsize=(9, 7))
    Ts = np.linspace(10, 373.9, 300)
    ax.plot(Ts, [psat_MPa_K(t + 273.15) * 10 for t in Ts], "k-", lw=2, label="saturation (IF97 Region 4)")
    for n, w in W.items():
        s = w["sdf"]
        ax.plot(s["T"], s.P, ".-", color=wcol[n], lw=0.8, ms=4, label=n)
        fzr = fzdf[fzdf.well == n]
        ax.plot(fzr.T_degC, fzr.P_bar, "s", color=wcol[n], mec="k", ms=7)
    ax.plot([], [], "s", color="w", mec="k", label="feed zone")
    ax.set_xlabel("Temperature (°C)")
    ax.set_ylabel("Pressure (bar)")
    ax.set_xlim(0, 320)
    ax.set_ylim(0, 300)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, ncol=2)
    ax.set_title("Survey P–T state against the saturation curve (liquid lies right/below... left of curve = subcooled)", fontsize=9)
    ax.set_title("Survey P–T state against the saturation curve (points left of the curve are compressed liquid)", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "PT_state.png"), dpi=130)
    plt.close(fig)

    # maps at elevations
    for var, lab, cmap, fname in (("T_degC", "T (°C)", "inferno", "maps_temperature.png"),
                                  ("P_bar", "P (bar)", "viridis", "maps_pressure.png")):
        fig, axs = plt.subplots(2, 3, figsize=(16, 10.5))
        for ax, z in zip(axs.flat, MAP_ELEVS):
            draw_base(ax, grid=False, labels=False)
            s = ce[ce.elevation_mRL == z]
            vmin, vmax = s[var].min(), s[var].max()
            if len(s) >= 3:
                tri = Triangulation(s.x.values, s.y.values)
                cs = ax.tricontourf(tri, s[var].values, levels=12, cmap=cmap, alpha=0.75, zorder=2.5)
                fig.colorbar(cs, ax=ax, shrink=0.75, label=lab)
            sc = ax.scatter(s.x, s.y, c=s[var], cmap=cmap, vmin=vmin, vmax=vmax, edgecolor="k", s=50, zorder=6)
            for r in s.itertuples():
                ax.annotate(f"{r.well}\n{getattr(r, var):.0f}", (r.x, r.y), xytext=(4, 3), textcoords="offset points",
                            fontsize=6.5, zorder=7, bbox=dict(fc="white", ec="none", alpha=0.6, pad=0.3))
            ax.set_title(f"{lab} at {z} mRL ({len(s)} wells; interpolation inside well hull only)", fontsize=8.5)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, fname), dpi=110)
        plt.close(fig)

    # cross-section through hottest wells
    hot = ws.sort_values("Tmax_degC", ascending=False)
    hot_wells = list(hot[hot.Tmax_degC >= 200].well)
    if len(hot_wells) < 3:
        hot_wells = list(hot.well[:3])
    pts = np.array([[W[n]["x0"], W[n]["y0"]] for n in hot_wells])
    cen = pts.mean(axis=0)
    u, sv_, vt = np.linalg.svd(pts - cen)
    dirv = vt[0]
    if dirv[1] < 0:
        dirv = -dirv  # S -> N
    normal = np.array([-dirv[1], dirv[0]])
    SEC_HALF = 1200.0
    sec_rows = []
    for n, w in W.items():
        for z, T, P in w["sdf"][["z", "T", "P"]].itertuples(index=False):
            x, y = track_xy(w, z)
            d = np.dot([x - cen[0], y - cen[1]], dirv)
            off = np.dot([x - cen[0], y - cen[1]], normal)
            sec_rows.append((n, d, off, z, T, P))
    sec = pd.DataFrame(sec_rows, columns=["well", "dist_m", "offset_m", "z", "T", "P"])
    in_sec = sec.groupby("well").offset_m.apply(lambda s: s.abs().min() <= SEC_HALF)
    sec_wells = [n for n in wells if in_sec[n]]
    secs = sec[sec.well.isin(sec_wells)]
    secs.round(2).to_csv(os.path.join(OUT, "section_points.csv"), index=False)
    sec_order = sorted(sec_wells, key=lambda n: secs[secs.well == n].dist_m.iloc[0])
    log(f"Section through hottest wells {hot_wells}: wells within {SEC_HALF:.0f} m: {sec_order}")
    for var, lab, cmap, fname in (("T", "Temperature (°C)", "inferno", "section_T_hottest_wells.png"),
                                  ("P", "Pressure (bar)", "viridis", "section_P_hottest_wells.png")):
        fig, ax = plt.subplots(figsize=(13, 7.5))
        tri = Triangulation(secs.dist_m.values, secs.z.values)
        # mask long thin triangles that would bridge large gaps
        tr_x = secs.dist_m.values[tri.triangles]
        mask = (tr_x.max(axis=1) - tr_x.min(axis=1)) > 3000
        tri.set_mask(mask)
        cs = ax.tricontourf(tri, secs[var].values, levels=np.arange(0, 300, 20) if var == "T" else 15, cmap=cmap)
        if var == "T":
            cl = ax.tricontour(tri, secs[var].values, levels=[100, 150, 200, 250], colors="w", linewidths=0.8)
            ax.clabel(cl, fmt="%d", fontsize=7)
        fig.colorbar(cs, ax=ax, label=lab)
        for n in sec_order:
            s = secs[secs.well == n]
            w = W[n]
            ax.plot(s.dist_m, s.z, "k-", lw=1)
            d0 = s.dist_m.iloc[0]
            ax.text(d0, s.z.max() + 40, f"{n}\n(off {s.offset_m.abs().min():.0f} m)", ha="center", fontsize=7.5)
            for f in w["fz"]:
                dz = np.interp(f["mid"], s.z[::-1], s.dist_m[::-1])
                ax.plot(dz, f["mid"], "s", color="gold", mec="k", ms=6)
            dsh = np.interp(w["shoe"], s.z[::-1], s.dist_m[::-1])
            ax.plot(dsh, w["shoe"], "k_", ms=14, mew=2.5)
        ax.set_xlabel(f"Distance along section (m), direction ({dirv[0]:.2f}E, {dirv[1]:.2f}N) through {', '.join(hot_wells)}")
        ax.set_ylabel("Elevation (mRL)")
        ax.set_title(f"{lab} section through the hottest wells (wells within {SEC_HALF:.0f} m of the line; linear "
                     "interpolation between wells only). Squares = feed zones, bars = casing shoes", fontsize=9)
        ax.set_ylim(min(secs.z.min(), -3000) - 50, max(secs.z.max(), 300) + 150)
        fig.tight_layout()
        fig.savefig(os.path.join(FIG, fname), dpi=120)
        plt.close(fig)

    # ======================================================================================
    # Conceptual-model quantities (from data only)
    # ======================================================================================
    hottest = ws.sort_values("Tmax_degC", ascending=False).iloc[0]
    hw = W[hottest.well]
    cm_rows = []
    for r in ws.itertuples():
        w = W[r.well]
        d = math.hypot(w["x0"] - hw["x0"], w["y0"] - hw["y0"])
        az = (math.degrees(math.atan2(w["x0"] - hw["x0"], w["y0"] - hw["y0"])) + 360) % 360
        cm_rows.append(dict(well=r.well, dist_from_hottest_m=round(d, 0), azimuth_from_hottest_deg=round(az, 0),
                            Tmax_degC=r.Tmax_degC, Tmax_elev_mRL=r.Tmax_elev_mRL, inversion_degC=r.inversion_degC,
                            regime=r.regime, cap_base_mRL=r.cap_base_mRL, cap_thickness_m=r.cap_thickness_m,
                            T_at_minus1500=ceT.loc[r.well].get(-1500, np.nan),
                            P_at_minus1000=ceP.loc[r.well].get(-1000, np.nan),
                            nearest_fault=gc.set_index("well").loc[r.well, "nearest_fault"],
                            dist_nearest_fault_m=gc.set_index("well").loc[r.well, "dist_nearest_fault_m"]))
    cm = pd.DataFrame(cm_rows).sort_values("dist_from_hottest_m")
    cm.to_csv(os.path.join(OUT, "conceptual_model_evidence.csv"), index=False)
    outfl = cm[cm.regime == "outflow"]
    corr_T = float(np.corrcoef(cm.dist_from_hottest_m, cm.Tmax_degC)[0, 1])
    corr_z = float(np.corrcoef(outfl.dist_from_hottest_m, outfl.Tmax_elev_mRL)[0, 1]) if len(outfl) > 2 else np.nan

    # conceptual figure
    fig, (a, b) = plt.subplots(1, 2, figsize=(17, 8), gridspec_kw=dict(width_ratios=[1, 1.25]))
    draw_base(a, grid=False, labels=False)
    regcol = {"upflow": "red", "outflow": "orange", "two-phase": "magenta", "conductive": "tab:blue",
              "recharge/downflow": "navy"}
    for r in ws.itertuples():
        w = W[r.well]
        a.scatter(w["x0"], w["y0"], s=40 + 2.0 * r.Tmax_degC, color=regcol[r.regime], edgecolor="k", zorder=7)
        a.annotate(f"{r.well}\n{r.Tmax_degC:.0f}°C @ {r.Tmax_elev_mRL:.0f}", (w["x0"], w["y0"]), xytext=(7, -4),
                   textcoords="offset points", fontsize=7, zorder=8, bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.3))
    # upflow zone: circle around wells with Tmax within 15 degC of the hottest and no large inversion
    up = [r.well for r in ws.itertuples() if r.Tmax_degC >= hottest.Tmax_degC - 15]
    upx = np.mean([W[n]["x0"] for n in up])
    upy = np.mean([W[n]["y0"] for n in up])
    rad = max(400, max(math.hypot(W[n]["x0"] - upx, W[n]["y0"] - upy) for n in up) + 300)
    a.add_patch(plt.Circle((upx, upy), rad, fc="red", alpha=0.12, ec="red", ls="--", zorder=2))
    a.text(upx, upy + rad + 80, "inferred upflow zone\n(interpretation)", ha="center", color="darkred", fontsize=8, zorder=8)
    for r in outfl.itertuples():
        w = W[r.well]
        if r.well in up:
            continue
        a.annotate("", xy=(w["x0"], w["y0"]), xytext=(upx, upy),
                   arrowprops=dict(arrowstyle="-|>", color="darkorange", lw=1.5 + 0.02 * (r.inversion_degC or 0),
                                   alpha=0.85), zorder=6)
    a.legend(handles=[Line2D([], [], marker="o", ls="", color=c, mec="k", label=k) for k, c in regcol.items()] +
             [Line2D([], [], color="darkorange", marker=">", label="inferred outflow direction")],
             fontsize=7, loc="lower right", title="regime (size ~ Tmax)", title_fontsize=7)
    for r in ws.itertuples():
        if r.regime in ("recharge/downflow", "conductive"):
            w = W[r.well]
            a.text(w["x0"], w["y0"] - 260, "cool margin /\nrecharge?" if r.regime != "conductive" else "conductive\nmargin",
                   ha="center", va="top", fontsize=7, color="navy", style="italic", zorder=8)
    a.set_title("Plan: regimes, inferred upflow and outflow paths (interpretation)", fontsize=9)
    # section panel
    tri = Triangulation(secs.dist_m.values, secs.z.values)
    tr_x = secs.dist_m.values[tri.triangles]
    tri.set_mask((tr_x.max(axis=1) - tr_x.min(axis=1)) > 3000)
    cs = b.tricontourf(tri, secs["T"].values, levels=np.arange(0, 300, 20), cmap="inferno", alpha=0.85)
    fig.colorbar(cs, ax=b, label="T (°C)")
    for n in sec_order:
        s = secs[secs.well == n]
        b.plot(s.dist_m, s.z, "k-", lw=0.8)
        b.text(s.dist_m.iloc[0], s.z.max() + 40, n, ha="center", fontsize=7.5)
        cb = W[n]["summary"]["cap_base_mRL"]
        if cb is not None and np.isfinite(cb):
            b.plot(np.interp(cb, s.z[::-1], s.dist_m[::-1]), cb, "D", color="cyan", mec="k", ms=6, zorder=6)
        r = W[n]["summary"]
        if r["inversion_degC"]:
            b.plot(np.interp(r["Tmax_elev_mRL"], s.z[::-1], s.dist_m[::-1]), r["Tmax_elev_mRL"], ">", color="darkorange",
                   mec="k", ms=10, zorder=6)
    lowm = tgt[tgt.low_margin_flag & tgt.well.isin(sec_order)]
    for r in lowm.itertuples():
        s = secs[(secs.well == r.well) & (secs.z == r.elevation_mRL)]
        b.plot(s.dist_m, s.z, "o", color="magenta", mec="k", ms=6, zorder=7)
    hs_ = secs[secs.well == hottest.well]
    deep_hot = ws[(ws.regime_note.fillna("").str.contains("reheats")) & (ws.T_at_TD_degC >= 200)]
    if len(hs_):
        src = secs[secs.well == deep_hot.well.iloc[0]] if len(deep_hot) and deep_hot.well.iloc[0] in sec_order else hs_
        b.annotate("", xy=(hs_.dist_m.iloc[-1], hs_.z.min()), xytext=(src.dist_m.iloc[-1], src.z.min()),
                   arrowprops=dict(arrowstyle="-|>", color="darkred", lw=2.5, connectionstyle="arc3,rad=-0.2"))
        b.text(src.dist_m.iloc[-1], src.z.min() - 150, "inferred upflow (interpretation)", ha="center", fontsize=8,
               color="darkred")
        b.text(hs_.dist_m.iloc[0] + 60, -110, "near-boiling (margin < 5 °C)", ha="left", va="center", fontsize=7,
               color="magenta", bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.3), zorder=8)
    b.legend(handles=[Line2D([], [], marker="D", color="cyan", mec="k", ls="", label="cap base (gradient < 100 °C/km)"),
                      Line2D([], [], marker=">", color="darkorange", mec="k", ls="", label="outflow tongue (T max above an inversion)"),
                      Line2D([], [], marker="o", color="magenta", mec="k", ls="", label="boiling margin < 5 °C")],
             fontsize=7, loc="lower left")
    b.set_xlabel("Distance along section (m)")
    b.set_ylabel("Elevation (mRL)")
    b.set_title("Section through hottest wells with interpreted features", fontsize=9)
    b.set_ylim(min(secs.z.min(), -3000) - 600, max(secs.z.max(), 300) + 150)
    fig.suptitle("Conceptual model built from the well data alone (symbols = observations; arrows/zones = interpretation)",
                 fontsize=10)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "conceptual_model.png"), dpi=120)
    plt.close(fig)

    # ======================================================================================
    # REPORT
    # ======================================================================================
    ctx = dict(W=W, ws=ws, fzdf=fzdf, hy=hy, hsum=hsum, cold=cold, cold_med=cold_med, pc=pc, trend=trend, tcoef=tcoef,
               iso=iso, ig=ig, xc=xc, irregular=irregular, inv=inv, qmap=qmap, tgt=tgt, ceT=ceT, ceP=ceP, gc=gc,
               cm=cm, corr_T=corr_T, corr_z=corr_z, hottest=hottest, up=up, inj_T=inj_T, passed=passed,
               gauge_info=gauge_info, first_prod_field=first_prod_field, first_inj_field=first_inj_field,
               survey_dates=survey_dates, pip=pip, prov=prov, sec_order=sec_order, hot_wells=hot_wells,
               wbw=wbw, bm=bm, lc=lc, geo=geo, polys=polys, skipped=skipped, ref=pd.read_csv(os.path.join(OUT, "reference_model_summary.csv")))
    write_report(ctx)

    pd.DataFrame(DERIVED, columns=["label", "value"]).drop_duplicates().to_csv(
        os.path.join(OUT, "report_derived_numbers.csv"), index=False)
    log("Done.")
    with open(os.path.join(OUT, "run_console_log.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG) + "\n")


# ======================================================================================
# Report writer
# ======================================================================================
def write_report(c):
    W, ws, fz, hy, hsum = c["W"], c["ws"], c["fzdf"], c["hy"], c["hsum"]
    pc, trend, iso, ig, xc, inv = c["pc"], c["trend"], c["iso"], c["ig"], c["xc"], c["inv"]
    ceT, ceP, gc, cm, tgt = c["ceT"], c["ceP"], c["gc"], c["cm"], c["tgt"]
    wells = list(W)
    wsi = ws.set_index("well")
    gci = gc.set_index("well")
    L = []
    A = L.append

    def wl(n):
        return wsi.loc[n]

    hot = c["hottest"].well
    conf_counts = ws.ns_confidence.value_counts()
    outfl = ws[ws.regime == "outflow"]
    coldw = c["cold"]
    p1000 = ceP[-1000].dropna()
    p1000_main = p1000.drop([n for n in p1000.index if p1000[n] < p1000.median() - PCONN_TOL])
    p1000_low = p1000[p1000 < p1000.median() - PCONN_TOL]
    deep_open = hy[(hy.band == "open") & (hy.z_top_mRL <= -1000)]
    rnd_mult = (c["cold_med"]["open"] - 1) / hy[hy.well.isin(c["cold"]) & (hy.band == "open")].ratio_rounding_uncertainty.max()
    dd1 = {n: gi_[ "ser"][0][1] - gi_["ser"][1][1] for n, gi_ in c["gauge_info"].items() if n in c["passed"]}
    rho_cold_max = max(rho_liq(float(W[n]["sdf"]["T"].min()), float(W[n]["sdf"].P.max())) for n in W)

    def pl(nn, word):
        return f"{nn} {word}" + ("" if nn == 1 else "s")
    t1000 = ceT[-1000].dropna()
    cased = hy[hy.band == "cased"]
    opn = hy[hy.band == "open"]
    bm = c["bm"]
    deep_hot = ws[(ws.regime_note.fillna("").str.contains("reheats")) & (ws.T_at_TD_degC >= 200)]
    dh = deep_hot.well.iloc[0] if len(deep_hot) else None
    capw = ws[ws.cap_base_mRL.notna()]
    fi = c["first_prod_field"]
    gi = c["gauge_info"]
    passed = c["passed"]
    fz_exc = fz[fz.well.isin(passed)]

    A("# Natural-state data interpretation for model calibration")
    A("")
    A("Generated by `Scripts/natural_state_analysis.py` from `Baseline/` only. Every number below is taken "
      "from a table in `Output/` (derived numbers are also listed in `report_derived_numbers.csv`). "
      "Statements are labelled **[obs]** (observation in the data) or **[interp]** (interpretation). "
      "Content describing the existing model is labelled **[reference model]** and is never used as evidence.")
    A("")

    # ------------------------------------------------------------------ 1
    A("## 1. Key conclusions for calibration")
    A("")
    k = 1
    w10 = wl(hot)
    A(f"{k}. **Upflow beneath the {hot} / {dh} area, at the junction of Fault 1 and Fault 3.** [obs] {hot} has the highest "
      f"temperature in the field, {F(w10.Tmax_degC, 1)} °C at {F(w10.Tmax_elev_mRL, 0)} mRL, and sits on the boiling curve "
      f"(minimum boiling margin {F(w10.min_boiling_margin_degC, 2)} °C at {F(w10.min_boiling_margin_elev_mRL, 0)} mRL). "
      f"{dh} reheats below its inversion to {F(wl(dh).T_at_TD_degC, 1)} °C at its total depth {F(wl(dh).TD_mRL, 0)} mRL, the "
      f"hottest deep temperature measured. {hot} and {dh} are {F(gci.loc[hot, 'dist_Fault_3_m'], 1)} m and "
      f"{F(gci.loc[dh, 'dist_Fault_3_m'], 1)} m from Fault 3, near its junction with the north end of Fault 1. "
      f"[interp] Hot fluid rises at depth beneath {dh} and reaches its shallowest level at {hot}. "
      f"**Action:** put the deep heat and mass input under the {dh}–{hot} columns, with a permeable vertical path "
      f"along the fault junction up to about {F(w10.Tmax_elev_mRL, 0)} mRL.")
    k += 1
    A(f"{k}. **A lateral outflow tongue at about −450 to −750 mRL spreads south, south-east, east and north-east from the upflow.** "
      f"[obs] {len(outfl)} wells ({', '.join(outfl.well)}) have temperature inversions of "
      f"{F(outfl.inversion_degC.min(), 1)}–{F(outfl.inversion_degC.max(), 1)} °C, with the temperature peak at "
      f"{F(outfl.Tmax_elev_mRL.max(), 0)} to {F(outfl.Tmax_elev_mRL.min(), 0)} mRL. The peak temperature falls with distance "
      f"from {hot} (correlation coefficient {F(c['corr_T'], 2)}; see `conceptual_model_evidence.csv`). "
      f"[interp] This is a hot, laterally permeable aquifer fed from the upflow, overlying cooler rock. "
      f"**Action:** make a laterally permeable layer at roughly −400 to −800 mRL connecting the upflow to "
      f"{', '.join(outfl.well[~outfl.well.isin([dh])])}. Make the rock below it, down to the inversion bases at "
      f"{F(outfl.inversion_base_mRL.max(), 0)} to {F(outfl.inversion_base_mRL.min(), 0)} mRL, less permeable or recharged "
      f"from the side so the inversions survive.")
    k += 1
    A(f"{k}. **A conductive cap covers the hot area.** [obs] The shallow gradients in the hot wells are "
      f"{F(capw.shallow_gradient_degC_per_km.min(), 1)}–{F(capw.shallow_gradient_degC_per_km.max(), 1)} °C/km. "
      f"The cap base, where the local gradient first falls below {F(CAP_GRAD, 0)} °C/km, lies at "
      f"{F(capw.cap_base_mRL.max(), 0)} to {F(capw.cap_base_mRL.min(), 0)} mRL, giving {F(capw.cap_thickness_m.min(), 0)}–"
      f"{F(capw.cap_thickness_m.max(), 0)} m of cap below the wellhead. The cap is thinnest at "
      f"{capw.sort_values('cap_thickness_m').well.iloc[0]} and {capw.sort_values('cap_thickness_m').well.iloc[1]}. "
      f"**Action:** use low-permeability cap rock from the surface down to about −250 mRL over the upflow, deepening to about "
      f"−550 mRL toward the outflow wells.")
    k += 1
    A(f"{k}. **No well is two-phase, but {hot} is at boiling point for depth.** [obs] {len(bm)} points, all in {hot} "
      f"between {F(bm.elevation_mRL.max(), 0)} and {F(bm.elevation_mRL.min(), 1)} mRL, have a boiling margin below "
      f"{F(MARGIN_FLAG, 0)} °C. The smallest is {F(bm.boiling_margin_degC.min(), 2)} °C. Every other well keeps at least "
      f"{F(ws[ws.well != hot].min_boiling_margin_degC.min(), 2)} °C of margin. **Action:** the model must hold the "
      f"{hot} column within a few °C of saturation without boiling it dry. If the model boils there, it has too much heat "
      f"or too little pressure at the top of the upflow.")
    k += 1
    A(f"{k}. **Deep pressures are connected across the field.** [obs] At −1000 mRL, {len(p1000_main)} wells read "
      f"{F(p1000_main.min(), 2)}–{F(p1000_main.max(), 2)} bar, although their temperatures there range from "
      f"{F(t1000[p1000_main.index].min(), 2)} to {F(t1000[p1000_main.index].max(), 2)} °C. Every feed-zone pressure lies within "
      f"{F(trend.residual_bar.min(), 2)} to +{F(trend.residual_bar.max(), 2)} bar of one linear trend through the feed zones "
      f"of the gauge-verified wells ({', '.join(passed)}). "
      + (f"The exception is {', '.join(p1000_low.index)} on the eastern margin, at {F(p1000_low.min(), 2)} bar, "
         f"{F(p1000_main.min() - p1000_low.min(), 2)} bar below the others. " if len(p1000_low) else "")
      + "[interp] The reservoir below the cap is a single, well-connected pressure system. **Action:** set high lateral "
      f"permeability at depth, and treat feed-zone pressures as one connected data set; their spread is only a few bar.")
    k += 1
    A(f"{k}. **The surveys represent the natural state, despite their date label.** [obs] Every survey is labelled "
      f"{', '.join(c['survey_dates'])}, but field production began on {fi:%d/%m/%Y}. In {', '.join(passed)} the survey pressure at the "
      f"gauge matches the first gauge reading to within {F(GAUGE_TOL, 1)} bar, while the next annual readings are "
      + ", ".join(f"{F(gi[n]['ser'][1][1], 1)} bar" for n in passed) +
      ". A survey made in 2015 could not match. Verdicts: "
      + ", ".join(f"{pl(conf_counts.get(x, 0), 'well')} {x}" for x in ("high", "medium", "low")) +
      " confidence (§2.3). **Action:** use all surveys as natural-state targets with the weights in "
      "`NS_calibration_targets.csv`, and ask the data owner for the true survey dates.")
    k += 1
    A(f"{k}. **Pressures above the feed zones are not usable as a hydrostatic column.** [obs] Behind casing, the measured ΔP is "
      f"{F(cased.ratio_measured_to_IF97.min(), 4)}–{F(cased.ratio_measured_to_IF97.max(), 4)} times the IF97 liquid column for "
      f"the measured temperatures, with implied densities up to {F(cased.rho_implied_kgm3.max(), 1)} kg/m³. In open hole "
      f"below −1000 mRL the ratio is only {F(deep_open.ratio_measured_to_IF97.min(), 4)}–"
      f"{F(deep_open.ratio_measured_to_IF97.max(), 4)}. A datum offset, a unit or gravity convention, two-phase fluid and "
      f"rounding cannot explain the shallow excess (§4.1). "
      f"**Action:** use pressure only at feed zones (`weight_P` is zero elsewhere), and do not tune shallow permeability to fit "
      f"the cased pressure profile.")
    k += 1
    A(f"{k}. **Cold margins on the north and south-west.** [obs] {', '.join(coldw)} never exceed {F(COLD_TMAX, 0)} °C. "
      f"{coldw[0]} has a mean gradient of {F(wl(coldw[0]).mean_gradient_degC_per_km, 1)} °C/km and {coldw[1]} of "
      f"{F(wl(coldw[1]).mean_gradient_degC_per_km, 1)} °C/km. WELL-1 in the east is conductive, reaching "
      f"{F(wl('WELL-1').T_at_TD_degC, 1)} °C at {F(wl('WELL-1').TD_mRL, 0)} mRL. [interp] These are recharge or cold-margin "
      f"areas. **Action:** use open, cold side boundaries or rain recharge near {', '.join(coldw)}, and no heat input there "
      f"beyond the background conductive flux.")
    k += 1
    nb = gc[gc.wellhead_minus_surface_m < -10]
    A(f"{k}. **The model topography disagrees with the wellheads.** [obs] In {len(nb)} wells the wellhead is "
      f"{F(-nb.wellhead_minus_surface_m.max(), 1)}–{F(-nb.wellhead_minus_surface_m.min(), 1)} m **below** the model surface "
      f"of its column; WELL-1 is the extreme case. **Action:** check the wellhead datums or the grid surface before "
      f"calibrating near-surface temperatures. Weights treat these points as below the surface.")
    A("")

    # ------------------------------------------------------------------ 2
    A("## 2. The data")
    A("")
    A("### 2.1 Inventory")
    A("")
    A(md_table(inv[["role", "path", "note"]]))
    A("")
    A("Roles that could not be filled: none. `reference_results` is the initial-conditions file only; no model results "
      "(listing) files are present. Production and injection records inside the well files were read only for the "
      "timing test in §2.3.")
    A("")
    A("### 2.2 Checks and provenance")
    A("")
    nis = ig[ig.status == "ISSUE"]
    A(f"[obs] {len(ig)} integrity checks were run (`integrity_checks.csv`) and {len(nis)} raised an issue:")
    A("")
    A(md_table(nis[["well", "check", "detail"]]))
    A("")
    A(f"- **JSON vs plain-column files.** All {len(xc)} plain-column files match the JSON exactly, value by value "
      f"(largest difference {F(float(xc.max_abs_diff.max()), 1)}; `dat_json_crosscheck.csv`). Transient times are years "
      f"since {xc.reference_date.dropna().iloc[0].split(' ')[0]}, quantised to months. Production mass flows have the opposite "
      f"sign in the plain-column files (simulator convention). Irregular file names: {', '.join(c['irregular'])} "
      f"(prefix `Wel10` instead of `Well10`; matched by number and content).")
    A(f"- **Every survey point and feed zone lies inside its track. Every feed zone lies inside a slotted liner.** "
      f"Casing strings overlap the liner tops correctly. In WELL-4 the liner is listed down to "
      f"{F(min(d['bottom'] for d in W['WELL-4']['design']), 0)} mRL, but the track and the grid file end at "
      f"{F(W['WELL-4']['track'][:, 2].min(), 0)} mRL.")
    A(f"- **Block assignment.** Survey points and feed zones take their blocks from the well-to-block file. An independent "
      f"point-in-polygon search on the grid file agrees for {int(c['pip'].agree.sum())} of {len(c['pip'])} points "
      f"(`column_assignment_check.csv`). The grid-file WELLS section matches the JSON tracks.")
    A("- **Provenance labels** (`provenance_labels.csv`). Every `Well_Design` block has `Source: 'Synthetic data'`. All "
      "survey, transient and flow `Source` and `Comment` fields are empty strings. Every survey has `StableFlag: True`. "
      "Feed-zone quality tags are HIGH, MEDIUM and LOW, except WELL-10, whose only tag is `'Synthetic'`, which is not a "
      "grade. Track elevations are labelled `masl` and feed zones `m`, but both are used as mRL, as instructed. "
      "**The dataset declares itself synthetic.** The tags give no basis for ranking wells by measurement quality, apart "
      "from the feed-zone grades.")
    A("")
    A("### 2.3 Natural-state timing verdict for each well")
    A("")
    A(f"Field production began on {fi:%d/%m/%Y} (first non-zero production record in any well); the first injection was on "
      f"{c['first_inj_field']:%d/%m/%Y}. Every survey carries the date label {', '.join(c['survey_dates'])}. "
      f"Injectate enthalpy averages {F(c['inj_T']['WELL-4']['h_mean'], 1)} kJ/kg, which IF97 Region 1 converts to "
      f"{F(c['inj_T']['WELL-4']['T'], 1)} °C at an assumed {F(P_INJ_ASSUMED, 1)} bar. Pressure has a negligible effect on "
      f"this conversion.")
    A("")
    t = pc[["well", "role_in_field", "gauge_elev_mRL", "survey_P_at_gauge_bar", "first_gauge_date", "first_gauge_P_bar",
            "diff_survey_minus_first_bar", "verdict", "confidence"]].copy()
    A(md_table(t, {"gauge_elev_mRL": 1, "survey_P_at_gauge_bar": 2, "first_gauge_P_bar": 1, "diff_survey_minus_first_bar": 2}))
    A("")
    A("Reasons (full text in `pressure_consistency_check.csv`):")
    A("")
    for r in pc.itertuples():
        A(f"- **{r.well}** ({r.verdict}, {r.confidence}): {r.reason}.")
    A("")
    A(f"For wells without a gauge, the feed-zone pressures were compared with a linear P(z) trend through the "
      f"feed zones of {', '.join(passed)}: P = {F(c['tcoef'][1], 2)} {'−' if c['tcoef'][0] < 0 else '+'} "
      f"{F(abs(c['tcoef'][0]), 4)} × z bar (`pressure_trend_check.csv`). All residuals lie between "
      f"{F(trend.residual_bar.min(), 2)} and +{F(trend.residual_bar.max(), 2)} bar. A survey made during field production "
      f"would sit about as far below the trend as the drawdown seen in the passed gauges after one year "
      f"({', '.join(f'{n} {F(v, 1)} bar' for n, v in dd1.items())}).")
    A("")
    A("**How confidence is used.** The confidence factor multiplies the location weight (high 1.0, medium 0.75, low 0.5). "
      f"Confidence differs between wells ({', '.join(f'{v} {k}' for k, v in conf_counts.items())}), so it changes the "
      "relative weights. WELL-5 is low because it is an injector whose formation is hotter than the injectate, so injection "
      "cooling could be hidden.")
    A("")

    # ------------------------------------------------------------------ 3
    A("## 3. Temperature structure")
    A("")
    t = ws[["well", "Tmax_degC", "Tmax_elev_mRL", "T_at_TD_degC", "TD_mRL", "inversion_degC", "inversion_top_mRL",
            "inversion_base_mRL", "shallow_gradient_degC_per_km", "deep_gradient_degC_per_km", "z_100C_mRL", "z_200C_mRL",
            "min_boiling_margin_degC", "regime"]].copy()
    A(md_table(t, {"Tmax_elev_mRL": 0, "TD_mRL": 1, "inversion_degC": 1, "inversion_top_mRL": 0, "inversion_base_mRL": 0,
                   "z_100C_mRL": 1, "z_200C_mRL": 1, "min_boiling_margin_degC": 2}))
    A("")
    A("Definitions (see `analysis_parameters.csv`):")
    A("- The inversion is the largest drop below the running maximum, reported if it is at least 3 °C.")
    A("- The shallow gradient is a least-squares fit over the top 300 m of the liquid-filled survey (P ≥ 1.5 bar).")
    A("- The deep gradient is a fit below both the inversion base and the deepest feed zone. It is blank if that interval is "
      "under 300 m. "
      f"WELL-2's interval is {F(wl('WELL-2').deep_interval_top_mRL - wl('WELL-2').TD_mRL, 0)} m, but it holds only 2 survey "
      f"points {F(wl('WELL-2').deep_gradient_points_span_m, 0)} m apart, so treat that value as weak.")
    A("- Regime rules, applied in order: two-phase (margin < 2 °C), outflow (inversion ≥ 10 °C), upflow (Tmax ≥ 200 °C), "
      "recharge/downflow (mean gradient < 30 °C/km), otherwise conductive.")
    A("")
    A("**Well by well**, with each feature located against the casing shoe and feed zones. Features are [obs]; readings "
      "are [interp]:")
    A("")
    for n in wells:
        r = wl(n)
        w = W[n]
        txt = (f"- **{n}** ({r.regime}; shoe {F(w['shoe'], 0)} mRL). Tmax {F(r.Tmax_degC, 1)} °C at {F(r.Tmax_elev_mRL, 0)} mRL")
        if pd.notna(r.inversion_degC):
            if r.inversion_top_mRL < w["shoe"]:
                loc = "entirely below the shoe, in open hole"
            elif r.inversion_base_mRL >= w["shoe"]:
                loc = "entirely behind casing"
            else:
                loc = "top behind casing, base in open hole"
            txt += f", then an inversion of {F(r.inversion_degC, 1)} °C down to {F(r.inversion_base_mRL, 0)} mRL ({loc})"
        txt += "."
        isos = iso[iso.well == n]
        for s in isos.itertuples():
            fzs = f"; contains feed zone {s.feedzones_inside}" if pd.notna(s.feedzones_inside) and str(s.feedzones_inside) != "" else ""
            txt += (f" Isothermal {F(s.Tmin_degC, 1)}–{F(s.Tmax_degC, 1)} °C from {F(s.top_mRL, 0)} to {F(s.bottom_mRL, 0)} mRL "
                    f"({s.location}{fzs}); feed-zone test: {s.internal_flow_test}.")
        if n == hot:
            txt += (" The near-boiling section starts behind casing, above the shoe, and continues into the feed zone. Behind "
                    "casing the well exchanges no fluid with the formation, so those temperatures are most likely formation "
                    "temperatures. [interp] Single-phase upflow at saturation. The flat base of the profile below the only "
                    "feed zone, 268.7 °C down to TD, rests on 3 points and cannot be interpreted further.")
        elif r.regime == "outflow":
            txt += " [interp] A hot lateral outflow at the Tmax level, over cooler rock."
            if len(isos):
                if isos.internal_flow_test.str.startswith("internal flow possible").any():
                    txt += (" The inversion below it may be partly a wellbore effect: internal flow between the bracketing "
                            "feed zones would flatten the profile. A static survey after long shut-in would decide. Weights "
                            "are reduced to 0.5 there.")
                else:
                    txt += (" No feed-zone pair brackets the deep isothermal section, so internal wellbore flow is "
                            "unlikely. It probably reflects cool, well-mixed or laterally recharged formation water "
                            "below the outflow.")
        elif n in coldw:
            txt += " [interp] Cold margin or recharge: low gradient and no hot fluid."
        elif r.regime == "conductive":
            txt += (f" Smooth profile with no inversion; the gradient falls from {F(r.shallow_gradient_degC_per_km, 1)} °C/km "
                    "near the surface. [interp] Conductive margin; the falling gradient may reflect a warm shallow lateral "
                    "flow or deeper cool recharge.")
        A(txt)
    A("")
    A("Temperatures at common elevations (°C; linear interpolation within each survey, no extrapolation; "
      "`temperature_at_common_elevations.csv`):")
    A("")
    tt = ceT.reset_index()
    tt.columns = ["well"] + [f"{int(x)} mRL" for x in ceT.columns]
    A(md_table(tt, {cc: 1 for cc in tt.columns}))
    A("")
    A("Maps: `figures/maps_temperature.png`. Section: `figures/section_T_hottest_wells.png`.")
    A("")

    # ------------------------------------------------------------------ 4
    A("## 4. Pressure structure")
    A("")
    A("Pressures at common elevations (bar; `pressure_at_common_elevations.csv`):")
    A("")
    tp = ceP.reset_index()
    tp.columns = ["well"] + [f"{int(x)} mRL" for x in ceP.columns]
    A(md_table(tp, {cc: 2 for cc in tp.columns}))
    A("")
    A(f"[obs] Below the cap, pressure at a given elevation hardly varies between wells. At −1000 mRL the spread is "
      f"{F(p1000_main.max() - p1000_main.min(), 2)} bar across {len(p1000_main)} wells. The cold wells are slightly higher "
      f"deeper down, as heavier columns would make them. [interp] This is a hydraulically connected liquid reservoir.")
    A("")
    A("**Feed-zone pressures (the formation-pressure targets)** and the pivot of each well (`feedzones.csv`, "
      "`well_summary.csv`). The dominant feed zone is the one with the highest quality grade, then the thickest, then the "
      "shallowest. The last column is the pressure in excess of a 20 °C water column (IF97) standing from the wellhead at "
      "1 atm:")
    A("")
    t = fz[["well", "feedzone", "mid_mRL", "quality", "T_degC", "P_bar", "boiling_margin_degC", "block",
            "P_excess_over_cold_column_bar"]]
    A(md_table(t, {"mid_mRL": 1, "T_degC": 2, "P_bar": 2, "boiling_margin_degC": 3, "P_excess_over_cold_column_bar": 2}))
    A("")
    A(f"[obs] Feed-zone pressures are {F(fz.P_excess_over_cold_column_bar.min(), 2)} to "
      f"+{F(fz.P_excess_over_cold_column_bar.max(), 2)} bar relative to a cold column from the wellhead. "
      f"Relative to a cold column from the model surface they are {F(fz.P_excess_over_cold_column_from_surface_bar.min(), 2)} "
      f"to +{F(fz.P_excess_over_cold_column_from_surface_bar.max(), 2)} bar (`feedzones.csv`). "
      f"[interp] Shallow feed zones in the hot wells are overpressured relative to any cold column from the ground; "
      f"WELL-3's −335 mRL feed zone is the clearest case. This can be physical, an upflow driven under a low-permeability "
      f"cap, and it agrees with the super-hydrostatic shallow gradients in §4.1. Deep feed zones are close to or below a "
      f"cold column, as connected hot water would be. A cold water table at the model surface explains much of WELL-1's "
      f"excess, because its wellhead lies 105 m below that surface.")
    A("")
    A("### 4.1 Hydrostatic check (single-phase liquid windows, IF97 Region 1, g = 9.80665 m/s²)")
    A("")
    A("Windows are non-overlapping and at least 300 m long, with boiling margin ≥ 10 °C and P ≥ 1.5 bar throughout. "
      "They are always split at the casing shoe and at feed zones where the minimum length allows. A band shorter than "
      "300 m is not tested (`hydrostatic_skipped.csv`). The rounding uncertainty of the ratio is ±0.1 bar / ΔP_IF97. "
      "Results by band (`hydrostatic_summary.csv`):")
    A("")
    t = hsum[hsum.band != "all"][["well", "band", "n_windows", "total_length_m", "median_ratio", "min_ratio", "max_ratio",
                                  "excess_over_cold_component", "max_rounding_uncertainty"]]
    A(md_table(t, {"total_length_m": 1, "median_ratio": 4, "min_ratio": 4, "max_ratio": 4,
                   "excess_over_cold_component": 4, "max_rounding_uncertainty": 4}))
    A("")
    cm_ = c["cold_med"]
    A(f"- **Field-wide (cold-well) component.** The cold wells ({', '.join(coldw)}, Tmax < 100 °C) give median ratios of "
      f"{F(cm_['cased'], 4)} (cased) and {F(cm_['open'], 4)} (open). [interp] This may be a datum or convention issue "
      f"rather than physical overpressure. A single gravity or unit convention would give the same ratio in both bands. The "
      f"open-hole value, {F((cm_['open'] - 1) * 100, 2)} % above IF97, is {F(rnd_mult, 1)} times the largest rounding "
      f"uncertainty of the cold-well open windows and could be a small convention or salinity effect. The cased value cannot be a convention, because it differs from the open-hole "
      f"value.")
    A(f"- **Excess above that component.** Cased windows in the hot wells exceed it by up to "
      f"{F(hsum[(hsum.band == 'cased') & (hsum.well != 'ALL')].excess_over_cold_component.max(), 4)} "
      f"(WELL-3, implied density {F(cased.rho_implied_kgm3.max(), 1)} kg/m³). In open hole the excess is largest in the "
      f"window through the outflow of WELL-2 ({F(opn.ratio_measured_to_IF97.max(), 4)}). Deep open-hole windows are within "
      f"{F(opn[opn.z_top_mRL <= -1000].ratio_measured_to_IF97.min(), 4)}–{F(opn[opn.z_top_mRL <= -1000].ratio_measured_to_IF97.max(), 4)} "
      f"of IF97.")
    A("- **Checked before calling this anomalous:**")
    A("  - A **datum offset** shifts P but cannot change ΔP over a window.")
    A("  - A **gravity or unit convention** would act at every depth. Deep open-hole ratios near 1.00 rule it out as the main cause.")
    A("  - **Two-phase fluid** makes a column lighter, not heavier.")
    A(f"  - **Rounding** is at most {F(hy.ratio_rounding_uncertainty.max(), 4)} in the ratio.")
    A(f"  - **A different thermal state** when P was measured: pure water at the coldest survey temperature and the highest "
      f"survey pressure has a density of {F(rho_cold_max, 1)} kg/m³ (IF97), below the implied densities of the hot-well "
      f"cased windows.")
    A("- [interp] **Result.** The pressure profile above the feed zones is not a static wellbore column in equilibrium with "
      "the measured temperatures. It behaves like a formation pressure profile with an upward, super-hydrostatic gradient "
      "through the cap. That is what synthetic, model-generated block pressures would show, and it is consistent with the "
      "dataset's own `Synthetic data` label. A real static survey could not show it. **Treat the cased pressures as a data "
      "question; they are not calibration targets.** Ask the data owner how the pressures were generated. If they are "
      "formation pressures, the cap needs low vertical permeability and upward flow to hold the excess pressure.")
    A("")

    # ------------------------------------------------------------------ 5
    A("## 5. Conceptual model (from the data alone)")
    A("")
    A("![conceptual model](figures/conceptual_model.png)")
    A("")
    cap_up = wl(hot).cap_base_mRL
    A(f"- **Upflow.** [interp] Upflow lies beneath {dh} and {hot}, close to the junction of Fault 1 and Fault 3. Supporting observations [obs]:")
    A(f"  - {hot} has the field maximum, {F(wl(hot).Tmax_degC, 1)} °C, at only {F(wl(hot).Tmax_elev_mRL, 0)} mRL, at boiling point for depth.")
    A(f"  - {dh} reheats to {F(wl(dh).T_at_TD_degC, 1)} °C at {F(wl(dh).TD_mRL, 0)} mRL.")
    A(f"  - At −250 mRL the temperature falls in every direction away from {hot} (`maps_temperature.png`).")
    A(f"  - *Undecided:* whether the deep upflow enters under {dh} (to the north, near Fault 3) or further west. No well near "
      f"the upflow goes deeper than {dh}'s {F(wl(dh).TD_mRL, 0)} mRL. A deep well, or a deeper survey in {hot}, would decide.")
    A(f"- **Outflow.** [interp] Hot water spreads at about −450 to −750 mRL in three directions:")
    A("  - south along Fault 1 (WELL-3, WELL-6);")
    A("  - south-east toward Fault 2 (WELL-7, WELL-5);")
    A("  - north-east (WELL-9).")
    A(f"  [obs] Each outflow well has an inversion below its peak. The peak cools with distance ({', '.join(f'{r.well} {F(r.Tmax_degC, 1)} °C' for r in cm[cm.regime == 'outflow'].itertuples())}). "
      f"The peak elevation varies from −450 to −750 mRL with no clear trend in distance (correlation coefficient "
      f"{F(c['corr_z'], 2)}). East of Fault 2, WELL-5 shows only a weak outflow ({F(wl('WELL-5').inversion_degC, 1)} °C "
      "inversion) and WELL-1 none. *Undecided:* whether Fault 1 and Fault 2 guide the outflow, or it is a stratigraphic "
      "aquifer that thins eastward. More wells between and beyond the faults at −500 mRL would decide.")
    A("- **Recharge and downflow.** [interp] Cold-margin recharge lies on the north (WELL-4, on Fault 3, far from the "
      "upflow) and the south-west (WELL-8, near the river and the lake). [obs] WELL-4 never exceeds 29.0 °C, and WELL-8 reaches only "
      "93.1 °C at −2716 mRL. Deep cool, nearly isothermal sections lie under the outflow (for example WELL-3 at 174.5–177.2 °C, "
      "−1100 to −2250 mRL). No feed-zone pair brackets them in WELL-3, WELL-5, WELL-6 or WELL-7 (`isothermal_sections.csv`). "
      "[interp] Cool water moves in laterally beneath the outflow tongue. *Alternative:* unrecognised minor feed zones driving "
      "internal flow. A long-shut-in survey would decide. In WELL-2 a feed-zone pair does bracket the section, so internal "
      "flow there is possible.")
    A(f"- **Cap.** [obs] The cap base lies at {F(capw.cap_base_mRL.max(), 0)} to {F(capw.cap_base_mRL.min(), 0)} mRL, "
      f"{F(capw.cap_thickness_m.min(), 0)}–{F(capw.cap_thickness_m.max(), 0)} m thick below the wellheads. It is shallowest "
      f"over the upflow (−250 mRL at {hot} and WELL-3) and deepens outward. No cap is resolved in the cold wells or in "
      "WELL-1.")
    A(f"- **Two-phase zone.** [obs] None is measured. {hot} is within 2.6–4.3 °C of boiling from −75 to −262.5 mRL. "
      "Above −75 mRL the margin grows steadily toward the top of the survey. [interp] The well therefore shows no steam "
      "zone. A two-phase zone could exist laterally, between WELL-10 and the fault junction where no well is drilled, or in "
      "the top 50 m. The survey stops 50 m below "
      "the wellhead.")
    A(f"- **Pressure connection.** [obs] Pressures at depth agree within a few bar across the field, and feed-zone "
      f"pressures follow one trend to within {F(trend.residual_bar.abs().max(), 2)} bar. WELL-1, in the east beyond Fault 2, "
      f"is about 4–5 bar lower at −1000 mRL. [interp] There is one connected reservoir. WELL-1 may sit beyond a partial "
      "barrier (Fault 2), or be drained toward the river; *undecided*. A gauge or interference test in WELL-1 would decide.")
    A("- **Surface water.** The river crosses the field between WELL-3/WELL-9 and WELL-6/WELL-7, and passes close to WELL-1 "
      f"({F(gci.loc['WELL-1', 'dist_river_m'], 1)} m) and WELL-8 ({F(gci.loc['WELL-8', 'dist_Tarawera_River_m'], 1)} m). "
      "No outflow reaches the surface in the data. WELL-1's water level is above its wellhead (estimated "
      f"{F(wl('WELL-1').water_level_mRL, 1)} mRL). [interp] Artesian conditions near the river; or the wellhead datum is wrong.")
    A("")

    # ------------------------------------------------------------------ 6
    A("## 6. Recommended calibration set-up")
    A("")
    A(f"**Targets** (`NS_calibration_targets.csv`, {len(tgt)} survey points; `feedzones.csv`, {len(fz)} feed zones):")
    A("- Temperature at every survey point, with weight = location weight × confidence factor:")
    A("  - feed zone 1.0;")
    A("  - open hole 0.8;")
    A("  - open hole in an isothermal section ≥ 300 m 0.5;")
    A("  - behind casing 0.5;")
    A("  - within 30 m of the wellhead or in the column's top block 0.3;")
    A("  - above the model surface 0.")
    A("  No point is above the model surface, because every wellhead is at or below it.")
    A("- Pressure at feed zones only: `weight_P` = feed-zone quality factor × confidence factor. Quality tags found and mapped:")
    for t_, (f_, why) in c["qmap"].items():
        A(f"  - `{t_}` → {f_} ({why})")
    A("- Match the boiling margin in WELL-10, not just its temperature: the model should stay liquid at margins of 2.6–4.3 °C.")
    A("")
    t = c["wbw"]
    A(md_table(t, {"sum_weight_T": 3}))
    A("")
    A("**Where heat and mass enter and leave** [interp]:")
    A("- **In.** A deep hot mass and heat source under the WELL-2–WELL-10 columns (fault junction). It must supply at least "
      "about 270 °C fluid; its enthalpy cannot be fixed from the data. Elsewhere at the base, background conductive heat: "
      f"deep gradients of {F(ws.deep_gradient_degC_per_km.min(), 1)}–{F(ws.deep_gradient_degC_per_km.max(), 1)} °C/km give "
      f"{F(ws.deep_conductive_heat_flux_W_m2.min(), 3)}–{F(ws.deep_conductive_heat_flux_W_m2.max(), 3)} W/m², using an "
      f"**assumed** conductivity of {F(K_ROCK, 1)} W/m/K (typical volcanic rock, not from the data). Under the outflow these "
      "gradients are lower bounds, because the inversions there are advective.")
    A("- **Out.** Lateral outflow at −450 to −750 mRL toward the south, south-east and north-east boundaries, or discharge to "
      "the river. Some vertical leakage through the cap over the upflow, needed to keep WELL-10 at boiling point for depth "
      "under a thin cap.")
    A("- **Recharge.** Cold water from the north (WELL-4) and south-west (WELL-8), and beneath the outflow tongue.")
    A("")
    A("**Permeability structure the data imply** [interp]:")
    A(f"- A vertical high-permeability path along the Fault 1 / Fault 3 junction, down to at least {F(wl(dh).TD_mRL, 0)} mRL "
      f"(the deepest hot measurement, in {dh}).")
    A("- A laterally permeable outflow layer at about −400 to −800 mRL.")
    A("- A low-permeability cap down to −250 mRL (centre) to −550 mRL (flanks).")
    A("- Moderate permeability below −1000 mRL, connected laterally: pressures there agree within a few bar.")
    A("- A possible partial barrier at Fault 2 toward WELL-1 (lower pressure).")
    A("")
    A("**Boundary requirements:**")
    A("- A top boundary that holds the water table near the surface. Check the wellhead-vs-surface mismatch first.")
    A("- Open (constant-pressure) lateral boundaries for the outflow and cold recharge on the north, south-west and east.")
    A("- The base boundary has a localised hot source under the fault junction and conductive flux elsewhere.")
    A("")
    ref = c["ref"]
    A("[reference model] For comparison only, not evidence. The starting model applies:")
    for r in ref.itertuples():
        A(f"  - {r.item.replace('[reference model] ', '')}: {r.quantity} = {r.value}")
    if REF_TOP.get("HEAT"):
        pts_ = [c["polys"][cc].vertices.mean(axis=0) for cc in REF_TOP["HEAT"] if cc in c["polys"]]
        cx_, cy_ = np.mean(pts_, axis=0)
        dref = math.hypot(cx_ - W[hot]["x0"], cy_ - W[hot]["y0"])
        DERIVED.append(("reference top-10 HEAT centroid distance to hottest well", f"{dref:.0f}"))
        pd.DataFrame([("[reference model] top-10 HEAT generator columns", ";".join(REF_TOP["HEAT"]),
                       round(cx_, 1), round(cy_, 1), round(dref, 0))],
                     columns=["item", "columns", "x_centroid", "y_centroid", "dist_to_hottest_well_m"]).to_csv(
            os.path.join(OUT, "reference_model_heat_centroid.csv"), index=False)
        A(f"  [reference model] The model's MASS generators all sit in the top layers and carry low-enthalpy water; they "
          f"look like rain recharge, not a hot upflow. Its heat enters as basal HEAT generators. The 10 largest of these are "
          f"centred {F(dref, 0)} m from {hot} (`reference_model_heat_centroid.csv`). The data call for a hot **mass** upflow "
          f"under the {dh}–{hot} columns (conclusion 1), which the starting model, as read here, does not have.")
    A("")

    # ------------------------------------------------------------------ 7
    A("## 7. Data quality issues and open questions (ranked by impact on calibration)")
    A("")
    A("1. **Pressure profiles above the feed zones are heavier than any water column** (§4.1). Up to "
      f"{F(cased.rho_implied_kgm3.max(), 1)} kg/m³ implied. Treat this as synthetic or formation pressure, not a wellbore "
      "column. Feed-zone pressures remain usable, but the generating convention should be confirmed.")
    A("2. **Survey date labels conflict with the gauges.** Every survey is labelled "
      f"{', '.join(c['survey_dates'])}, after production started on {fi:%d/%m/%Y}. The gauge test shows the surveys are pre-production "
      "(WELL-2, WELL-3 match within 0.1 bar). Ask for the real survey dates. This does not change the verdicts.")
    A(f"3. **Wellheads lie below the model surface** in {len(nb)} wells, by up to {F(-nb.wellhead_minus_surface_m.min(), 1)} m "
      "(`well_geometry_check.csv`). This affects shallow temperature targets and water levels.")
    A(f"4. **WELL-9's survey disagrees with its first gauge reading** by {F(gi['WELL-9']['diff'], 2)} bar, more than the "
      "reading precision. Its confidence is medium. A gauge depth or datum error is possible.")
    A(f"5. **WELL-10's gauge starts in {gi['WELL-10']['d0']:%Y}**, so no pre-production reading exists. The survey is "
      f"{F(gi['WELL-10']['diff'], 1)} bar above the first reading, which is less than the drawdown at the other gauges "
      f"({', '.join(F(x, 1) for x in gi['WELL-10'].get('passed_drawdown', []))} bar). The difference may be boiling buffering "
      "at the near-saturated feed zone. Its feed-zone tag `'Synthetic'` is not a grade.")
    A("6. **WELL-5 (injector, formation hotter than the injectate)**: injection cooling could be hidden, so its confidence "
      "is low.")
    A("7. **WELL-4 design**: the liner is listed to −3251 mRL, but the track ends at −750 mRL. The survey starts 110 m below "
      "the wellhead.")
    A("8. **Deep isothermal sections** (`isothermal_sections.csv`) could be wellbore internal flow. Weights are already "
      "reduced; a long-shut-in survey would decide.")
    A("9. **Deep coverage is thin near the upflow.** WELL-10 stops at −350 mRL, and WELL-2's deep gradient rests on 2 points.")
    A("10. **Irregular file names** (`Wel10_*`), and production flows with opposite signs in the plain-column files. Both "
      "are cosmetic.")
    A("")

    # ------------------------------------------------------------------ 8
    A("## 8. Output file index")
    A("")
    idx = [
        ("input_inventory.csv", "Role of every Baseline file/directory"),
        ("NS_calibration_targets.csv", "All survey points with block, flags and weights (main target file)"),
        ("feedzones.csv", "Feed-zone T, P, Tsat, margin, block, weights, overpressure"),
        ("T_P_at_layer_centres.csv", "T and P interpolated to model layer centres (no extrapolation)"),
        ("well_summary.csv", "Per-well characteristics, regime, verdict"),
        ("pressure_consistency_check.csv", "Gauge test and timing per well"),
        ("pressure_trend_check.csv", "Feed-zone pressures vs passed-well trend"),
        ("hydrostatic_check.csv / hydrostatic_summary.csv / hydrostatic_skipped.csv", "IF97 hydrostatic windows"),
        ("well_geometry_check.csv", "Wellhead vs surface, deviation, distances to faults and water"),
        ("TP_at_common_elevations.csv, temperature_/pressure_at_common_elevations.csv", "Lateral structure"),
        ("integrity_checks.csv, dat_json_crosscheck.csv, column_assignment_check.csv, provenance_labels.csv", "Data checks"),
        ("isothermal_sections.csv, boiling_margin_flags.csv, weights_by_well.csv", "Supporting tables"),
        ("conceptual_model_evidence.csv, section_points.csv", "Evidence for the conceptual model"),
        ("if97_validation.csv, analysis_parameters.csv, reference_model_summary.csv, report_derived_numbers.csv", "Method, parameters, reference model, report numbers"),
        ("figures/well_map.png", "Grid, faults, water bodies, wells"),
        ("figures/profile_<well>.png", "T, P, Tsat, boiling margin, casing, feed zones per well"),
        ("figures/maps_temperature.png, figures/maps_pressure.png", "Maps at 0 to −1500 mRL"),
        ("figures/section_T_hottest_wells.png, figures/section_P_hottest_wells.png", "Sections through the hottest wells"),
        ("figures/PT_state.png", "P–T against saturation"),
        ("figures/conceptual_model.png", "Conceptual model sketch"),
        ("run_console_log.txt, acceptance_check.txt, process_log.md", "Run log, acceptance check, process log"),
    ]
    A(md_table(pd.DataFrame(idx, columns=["file (in Output/)", "contents"])))
    A("")
    with open(os.path.join(OUT, "NS_data_interpretation.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    log("Report written: Output/NS_data_interpretation.md")


if __name__ == "__main__":
    main()
