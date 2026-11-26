"""Natural-state analysis for geothermal model calibration (CLAUDE.md v4).

Run from the repository root:

    python Scripts/natural_state_analysis.py

Regenerates every output under Output/ from Baseline/ (read-only): the input inventory,
the calibration targets, the interpretation report, the supporting tables and all figures.
Output is deterministic (no timestamps). Process log and console log are written separately.
"""
import ast
import datetime as dtm
import glob
import json
import math
import os
import re
import textwrap

import iapws
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402
from matplotlib.path import Path  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402
from scipy.interpolate import LinearNDInterpolator  # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BASE = os.path.join(ROOT, "Baseline")
CAL = os.path.join(ROOT, "Output", "Calibration_Input")
ADD = os.path.join(ROOT, "Output", "Additional_Output")
CALF = os.path.join(CAL, "figures")
ADDF = os.path.join(ADD, "figures")
for _d in (CAL, ADD, CALF, ADDF):
    os.makedirs(_d, exist_ok=True)

# ----------------------------------------------------------------------------------------
# Parameters (every value here is also written to Additional_Output/analysis_parameters.csv)
# ----------------------------------------------------------------------------------------
G = 9.80665              # m/s2, CLAUDE.md 4.3
P_ATM = 1.01325         # bar
P_CRIT = 220.64         # bar, IF97 critical pressure
T_CRIT = 373.946        # degC, IF97 critical temperature
WIN_MIN = 300.0         # m, minimum hydrostatic window length (§4.3)
ISO_RANGE = 5.0         # degC, isothermal-section range (§4.2)
ISO_SPAN = 300.0        # m, isothermal-section span (§4.2)
HYDRO_BM = 10.0         # degC, single-phase test boiling margin (§4.3)
HYDRO_P = 1.5           # bar, single-phase test pressure (§4.3)
NEAR_BOIL = 5.0         # degC, boiling-margin flag (4.3, 5.1)
TWO_PHASE_BM = 1.0      # degC, two-phase regime threshold (§4.4)
INVERSION_MIN = 3.0     # degC, inversion threshold (§4.4)
OUTFLOW_MIN = 10.0      # degC, outflow regime threshold (§4.4)
SHALLOW_WIN = 300.0     # m, shallow-gradient window (§4.4)
DEEP_SPAN = 300.0       # m, deep-gradient minimum span (§4.4)
COLD_T = 100.0          # degC, cold-well threshold (§4.3)
WELLHEAD_NEAR = 30.0    # m, weight rule (5.5)
GAUGE_PREC = 0.1        # bar, reading precision of the transient records [assumption]
TREND_TOL_MIN = 3.0     # bar, tolerance floor (§4.1)
TREND_TOL_RMS = 2.0     # x RMS of gauge-verified feed zones (§4.1)
TREND_NEAR = 500.0      # m, well-to-line distance for sections (8.2 F4)
BAND_M = 500.0          # m, elevation band thickness (8.2 F3)
BAND_MIN_WELLS = 3      # wells needed to end the band range (8.2 F3)
INJ_P = 10.0            # bar, pressure at which injectate enthalpy is converted to T [assumption]
ISO_T_BIG = (100.0, 200.0)
ISO_LEVELS = [50.0 * k for k in range(1, 7)]   # 50..300 degC
CONF_F = {"high": 1.0, "medium": 0.75, "low": 0.5}
LOC_F = {"above": 0.0, "top": 0.3, "feed": 1.0, "casing": 0.5, "iso": 0.5, "open": 0.8}
VE_NOTE = "vertical scale equal to horizontal unless stated"
T_MAP_STEP = 10.0       # degC colour step for sections (8.1)
P_DEV_STEP = 1.0        # bar colour step for pressure-deviation sections (8.1)
CLASS_T = [(10, 0.0, "within +/-10 degC"), (20, 1, "10-20"), (30, 2, "20-30"), (40, 3, "30-40"),
           (np.inf, 4, "> 40")]
CLASS_P = [(2, "within +/-2 bar"), (4, "2-4"), (6, "4-6"), (8, "6-8"), (np.inf, "> 8")]

PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b", "#e377c2",
           "#7f7f7f", "#bcbd22", "#17becf"]

FIG_DPI = 150

# ----------------------------------------------------------------------------------------
# IAPWS-IF97
# ----------------------------------------------------------------------------------------
_N = [None, 0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5,
      -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6,
      -0.23855557567849, 0.65017534844798e3]


def tsat_c(p_bar):
    """Saturation temperature (degC) at p_bar, IF97 Region 4 (eq. 31). NaN outside the range."""
    p = p_bar / 10.0
    if not (611.213e-6 <= p <= 22.064):
        return float("nan")
    n = _N
    b = p ** 0.25
    e = b * b + n[3] * b + n[6]
    f = n[1] * b * b + n[4] * b + n[7]
    g = n[2] * b * b + n[5] * b + n[8]
    d = 2 * g / (-f - math.sqrt(f * f - 4 * e * g))
    t = (n[10] + d - math.sqrt((n[10] + d) ** 2 - 4 * (n[9] + n[10] * d))) / 2
    return t - 273.15


def rho_liquid(t_c, p_bar):
    """Liquid density (kg/m3), IF97 Region 1 (via the iapws package)."""
    return iapws.IAPWS97(P=p_bar / 10.0, T=t_c + 273.15).rho


def injectate_t(h_kj):
    """Temperature (degC) of liquid at enthalpy h_kj (kJ/kg), IF97 Region 1, at INJ_P bar."""
    return iapws.IAPWS97(P=INJ_P / 10.0, h=h_kj).T - 273.15


def validate_thermo():
    """Make the script fail if the Region 4 implementation is off by more than 0.1 K."""
    errs = []
    for p, tk in [(1.0, 372.755919), (10.0, 453.035632), (100.0, 584.149488)]:
        errs.append(abs(tsat_c(p) + 273.15 - tk))           # IF97 verification table values
    for p in np.linspace(1.0, 220.0, 60):
        errs.append(abs(tsat_c(p) - (iapws.IAPWS97(P=p / 10.0, x=0).T - 273.15)))
    assert abs(tsat_c(P_ATM) - 99.97) < 0.1 and abs(tsat_c(100.0) - 311.0) < 0.1
    assert max(errs) < 0.1, "IF97 Region 4 implementation error"
    return max(errs)


# ----------------------------------------------------------------------------------------
# Small helpers
# ----------------------------------------------------------------------------------------
def lit(s):
    return ast.literal_eval(s) if isinstance(s, str) else s


def iso_date(s):
    return dtm.datetime.strptime(s, "%d/%m/%Y").date().isoformat()


def fnum(v, nd=2):
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return ""
    return round(float(v), nd)


def write_csv(df, path):
    df.to_csv(path, index=False, lineterminator="\n")


def rel(path):
    return os.path.relpath(path, ROOT).replace("\\", "/")


def dist_to_polyline(px, py, pts):
    """Distance (m) from point to polyline (array Nx2)."""
    a = pts[:-1]
    b = pts[1:]
    ab = b - a
    ap = np.array([px, py]) - a
    t = np.clip(np.sum(ap * ab, axis=1) / np.maximum(np.sum(ab * ab, axis=1), 1e-12), 0, 1)
    proj = a + ab * t[:, None]
    return float(np.min(np.hypot(*(proj - [px, py]).T)))


def seg_line_intersections(pts, A, u):
    """Parameters s (along the line A + s u) where the polyline pts crosses that line."""
    out = []
    for p, q in zip(pts[:-1], pts[1:]):
        d = q - p
        den = u[0] * d[1] - u[1] * d[0]
        if abs(den) < 1e-12:
            continue
        w = p - A
        t = (w[0] * d[1] - w[1] * d[0]) / den       # along the line
        k = (w[0] * u[1] - w[1] * u[0]) / den       # along the segment
        if 0.0 <= k <= 1.0:
            out.append(t)
    return out


# ----------------------------------------------------------------------------------------
# Loaders
# ----------------------------------------------------------------------------------------
def load_groups():
    g = json.load(open(os.path.join(BASE, "Teaching_Model_3D", "Data", "Groups_data.json"),
                       encoding="utf-8"))
    groups = {}
    for area, gd in g["Areas"].items():
        for grp, wells in gd.items():
            groups[grp] = wells
    order = []
    for grp, wells in groups.items():
        order += wells
    return order, groups


def load_wells(order):
    wells = {}
    for fn in sorted(glob.glob(os.path.join(BASE, "Teaching_Model_3D", "Data", "Data_json", "*.json"))):
        d = json.load(open(fn, encoding="utf-8"))
        name = d["Well_Name"]
        tr = np.array([lit(s) for s in d["Well_Track"]["Data"]], float)
        order_idx = np.argsort(-tr[:, 2], kind="stable")
        tr = tr[order_idx]
        casings = []
        for s in d["Well_Design"]["Data"]:
            top, bot, dia, slotted, spec = lit(s)
            casings.append(dict(top=float(top), bot=float(bot), diam=str(dia), slotted=bool(slotted),
                                spec=str(spec)))
        shoe = min(c["bot"] for c in casings if not c["slotted"])
        fz = []
        for s in d["Feed_Zones"]["Data"]:
            top, bot, q = lit(s)
            fz.append(dict(top=float(top), bot=float(bot), q=str(q)))
        # survey dictionaries: date -> {elev: value}
        tdate = list(d["Downhole_Temperature"].keys())
        pdate = list(d["Downhole_Pressure"].keys())
        assert len(tdate) == 1 and len(pdate) == 1, f"{name}: more than one survey date"
        tv = d["Downhole_Temperature"][tdate[0]]
        pv = d["Downhole_Pressure"][pdate[0]]
        tz = np.array([float(k) for k in tv["Data"]])
        tt = np.array([float(v) for v in tv["Data"].values()])
        pz = np.array([float(k) for k in pv["Data"]])
        pp = np.array([float(v) for v in pv["Data"].values()])
        it = np.argsort(-tz, kind="stable")
        ip = np.argsort(-pz, kind="stable")
        assert np.array_equal(tz[it], pz[ip]), f"{name}: T and P elevations differ"
        gauges = {}
        if "Transient_Pressure" in d:
            for el, v in d["Transient_Pressure"].items():
                gauges[float(el)] = sorted((iso_date(k), float(val)) for k, val in v["Data"].items())
        flows = {}
        for key in ("Production_Mass_Flow", "Injection_Mass_Flow"):
            if key in d:
                flows[key] = sorted((iso_date(k), float(val)) for k, val in d[key]["Data"].items())
        enth = None
        for key in ("Transient_Enthalpy", "Injection_Enthalpy"):
            if key in d:
                enth = sorted((iso_date(k), float(val)) for k, val in d[key]["Data"].items())
        kind = "injector" if "Injection_Mass_Flow" in flows else (
            "producer" if "Production_Mass_Flow" in flows else "none")
        flow_series = flows.get("Production_Mass_Flow") or flows.get("Injection_Mass_Flow") or []
        first_flow = next((dd for dd, v in flow_series if v > 0), None)
        wells[name] = dict(
            name=name, track=tr, wellhead_z=float(tr[0, 2]), casings=casings, shoe=float(shoe),
            fz=fz, survey_date=iso_date(tdate[0]), tz=tz[it], tt=tt[it], pz=pz[ip], pp=pp[ip],
            source_notes=dict(design_source=d["Well_Design"].get("Source", ""),
                              survey_source=tv.get("Source", ""), survey_comment=tv.get("Comment", ""),
                              survey_stable=tv.get("StableFlag", None)),
            gauges=gauges, flows=flows, kind=kind, first_flow=first_flow,
            enth=enth, fz_units=d["Feed_Zones"]["Units"], n_survey=len(tt),
            aut2=d.get("AUT2_Name", ""),
        )
    missing = [w for w in order if w not in wells]
    assert not missing, f"wells in Groups_data.json without a data file: {missing}"
    return {w: wells[w] for w in order}


def load_well_blocks(path, names):
    wj = json.load(open(path, encoding="utf-8"))
    out = {}
    for w in names:
        d = wj[w]
        blocks = []
        for b in d["Blocks"]:
            top, bot = max(b["interval"]), min(b["interval"])
            m = re.match(r"^([a-z]+)(\d+)$", b["block"])
            blocks.append(dict(block=b["block"], top=float(top), bot=float(bot),
                               col=m.group(1), layer=int(m.group(2))))
        fzb = []
        for k, v in d["Feed_Zones"].items():
            fzb.append(dict(feed=k, block=list(v.keys())[0]))
        out[w] = dict(blocks=blocks, fzb=fzb)
    return out


def load_geometry(path):
    sec = {}
    cur = None
    for ln in open(path, encoding="utf-8", errors="replace"):
        s = ln.rstrip("\n")
        if re.match(r"^[A-Z]{4,}\s*$", s):
            cur = s.strip()
            sec.setdefault(cur, [])
            continue
        if cur is not None:
            sec[cur].append(s)
    verts = {}
    for ln in sec.get("VERTICES", []):
        m = re.match(r"^([a-z]{3})(\d+\.\d\d)(\d+\.\d\d)\s*$", ln)
        if m:
            verts[m.group(1)] = (float(m.group(2)), float(m.group(3)))
    cols = []
    lines = sec.get("GRID", [])
    i = 0
    while i < len(lines):
        m = re.match(r"^([a-z]{3})\d*\s+(\d+)\s*$", lines[i])
        if m:
            n = int(m.group(2))
            names = [x.strip() for x in lines[i + 1:i + 1 + n]]
            cols.append((m.group(1), names))
            i += 1 + n
        else:
            i += 1
    layers = {}
    for ln in sec.get("LAYERS", []):
        p = ln.split()
        if len(p) == 3:
            layers[int(p[0])] = dict(bottom=float(p[1]), centre=float(p[2]))
    surf = {}
    for ln in sec.get("SURFA", []):
        p = ln.split()
        if len(p) == 2:
            surf[p[0]] = float(p[1])
    grid_top = layers[0]["bottom"]
    layer_info = {}
    for i in sorted(layers):
        if i == 0:
            continue
        layer_info[i] = dict(top=layers[i - 1]["bottom"], bottom=layers[i]["bottom"],
                             centre=layers[i]["centre"])
    return dict(verts=verts, cols=cols, layers=layer_info, surf=surf, grid_top=grid_top,
                n_vertex_lines=len(sec.get("VERTICES", [])), n_sections=sorted(sec.keys()))


def load_mapping(path):
    m = json.load(open(path, encoding="utf-8"))
    feats = []
    for name, v in m.items():
        typ = v.get("type", "")
        data = v["Data"]
        if typ == "multiplot":
            segs = [np.array(s, float) for s in data]
        elif typ == "fill":
            segs = [np.array(data, float)]
        else:
            segs = [np.array(data, float)]
        leg = v.get("legend", "")
        kind = "fault" if (v.get("c") == "k" or name.lower().startswith("fault")) else "water"
        feats.append(dict(name=name, type=typ, kind=kind, fill=(typ == "fill"), segs=segs,
                          legend=leg, colour=v.get("c", "")))
    return feats


def load_dat_pairs(folder):
    """Plain-column files: survey (value, elevation) and transient (time, value)."""
    out = []
    for fn in sorted(glob.glob(os.path.join(folder, "*.dat"))):
        data = np.loadtxt(fn, ndmin=2)
        out.append((os.path.basename(fn), data))
    return out


# ----------------------------------------------------------------------------------------
# Geometry helpers
# ----------------------------------------------------------------------------------------
def build_column_polys(geo):
    polys = []
    for name, vnames in geo["cols"]:
        pts = [geo["verts"][v] for v in vnames if v in geo["verts"]]
        if len(pts) >= 3:
            polys.append((name, Path(np.array(pts))))
    return polys


def column_surface(geo, polys, x, y):
    for name, path in polys:
        if path.contains_point((x, y)):
            return geo["surf"].get(name, geo["grid_top"])
    return float("nan")


def model_extent(geo):
    v = np.array(list(geo["verts"].values()))
    return v[:, 0].min(), v[:, 0].max(), v[:, 1].min(), v[:, 1].max()


def interp_track(track, z):
    """x, y at elevation z along the track (track sorted by descending z)."""
    zz = track[::-1, 2]
    return float(np.interp(z, zz, track[::-1, 0])), float(np.interp(z, zz, track[::-1, 1]))


def interp_profile(zs, vs, z):
    """Linear interpolation of a profile (zs descending) at z; NaN outside the survey range."""
    zs_a = zs[::-1]
    vs_a = vs[::-1]
    if z < zs_a[0] - 1e-9 or z > zs_a[-1] + 1e-9:
        return float("nan")
    return float(np.interp(z, zs_a, vs_a))


def block_for(blocks, z):
    cands = [b for b in blocks if b["bot"] - 0.5 <= z <= b["top"] + 0.5]
    if not cands:
        return None
    inside = [b for b in cands if b["bot"] <= z < b["top"]]
    return (inside or cands)[0]


def top_block_of_column(blocks, col):
    same = [b for b in blocks if b["col"] == col]
    if not same:
        return None
    return min(same, key=lambda b: b["layer"])["block"]


def isothermal_flags(z, t, open_mask):
    """Open-hole points inside a run of consecutive points spanning >= ISO_SPAN with T range <= ISO_RANGE."""
    flag = np.zeros(len(z), bool)
    n = len(z)
    for i in range(n):
        if not open_mask[i]:
            continue
        j = i
        lo = hi = t[i]
        last_ok = None
        while j < n and open_mask[j]:
            lo = min(lo, t[j])
            hi = max(hi, t[j])
            if hi - lo > ISO_RANGE:
                break
            if z[i] - z[j] >= ISO_SPAN:
                last_ok = j
            j += 1
        if last_ok is not None:
            flag[i:last_ok + 1] = True
    return flag


def quality_factors(all_tags):
    grades = {"HIGH": 3, "MEDIUM": 2, "LOW": 1}
    present = sorted({t for t in all_tags if t in grades}, key=lambda t: -grades[t])
    f = {}
    for t in all_tags:
        f[t] = 0.5
    if present:
        f[present[0]] = 1.0
    if len(present) >= 3:
        f[present[1]] = 0.8
    if len(present) >= 2:
        f[present[-1]] = 0.5
    return f, present


# ----------------------------------------------------------------------------------------
# Main analysis
# ----------------------------------------------------------------------------------------
def main():
    log = []

    def say(msg):
        print(msg)
        log.append(msg)

    err = validate_thermo()
    say(f"IF97 Region 4 check passed; max error vs IF97 verification and iapws = {err:.2e} K")

    order, groups = load_groups()
    wells = load_wells(order)
    wnames = list(wells.keys())
    well_colour = {w: PALETTE[i % len(PALETTE)] for i, w in enumerate(wnames)}
    wb = load_well_blocks(os.path.join(BASE, "Teaching_Model_3D", "Models", "TE14515",
                                       "gTE14515_wells.json"), wnames)
    geo = load_geometry(os.path.join(BASE, "Teaching_Model_3D", "Models", "TE14515", "gTE14515.dat"))
    feats = load_mapping(os.path.join(BASE, "Teaching_Model_3D", "Data", "Mapping", "Mapping.json"))
    polys = build_column_polys(geo)
    say(f"Grid: {len(geo['verts'])} vertices, {len(polys)} column polygons, "
        f"{len(geo['layers'])} layers, {len(geo['surf'])} surface entries")
    xmin, xmax, ymin, ymax = model_extent(geo)

    # --- inventory --------------------------------------------------------------------
    write_inventory()

    # --- cross-check plain-column files against JSON (4.1 step 1) ------------------------
    xc_rows = dat_cross_check(wells)
    write_csv(pd.DataFrame(xc_rows), os.path.join(ADD, "data_crosscheck.csv"))
    xc_bad = [r for r in xc_rows if r["status"] == "mismatch"]
    say(f"Cross-check: {len(xc_rows)} plain-column series compared with JSON; "
        f"{len(xc_bad)} not matching")

    # --- per-point targets -----------------------------------------------------------------
    tgt_rows = []
    for w in wnames:
        W = wells[w]
        B = wb[w]
        tr = W["track"]
        fz_all_top = [f["top"] for f in W["fz"]]
        fz_all_bot = [f["bot"] for f in W["fz"]]
        fz_hi = max(fz_all_top) if fz_all_top else -1e9
        for i, z in enumerate(W["tz"]):
            T = float(W["tt"][i])
            P = float(W["pp"][i])
            x, y = interp_track(tr, z)
            blk = block_for(B["blocks"], z)
            block = blk["block"] if blk else ""
            col = blk["col"] if blk else ""
            layer = blk["layer"] if blk else ""
            surf_col = geo["surf"].get(col, geo["grid_top"]) if col else float("nan")
            above = bool(z > surf_col)
            Ts = tsat_c(P)
            bm = (Ts - T) if np.isfinite(Ts) else float("nan")
            fz_idx = [k for k, f in enumerate(W["fz"]) if min(f["top"], f["bot"]) - 0.5 <= z <= max(f["top"], f["bot"]) + 0.5]
            in_fz = len(fz_idx) > 0
            behind = bool(z >= W["shoe"])
            qfl = W["fz"][fz_idx[0]]["q"] if in_fz else ""
            tops = top_block_of_column(B["blocks"], col) if col else None
            top_col = tops is not None and block == tops
            near_wh = (W["wellhead_z"] - z) <= WELLHEAD_NEAR
            tgt_rows.append(dict(well=w, x=x, y=y, elevation_mRL=float(z), T_degC=T, P_bar=P,
                                 Tsat_degC=Ts, boiling_margin_degC=bm, column=col, layer=layer,
                                 block=block, in_feedzone=in_fz, behind_casing=behind,
                                 above_model_surface=above, model_surface_here_mRL=surf_col,
                                 _fz=fz_idx, _top_col=top_col, _near_wh=near_wh, _iso=False,
                                 _open=not behind, _liquid=None))
    tgt = pd.DataFrame(tgt_rows)
    # liquid flag and isothermal sections in open hole, per well
    tgt["_liquid"] = (tgt.P_bar >= HYDRO_P) & ((tgt.boiling_margin_degC >= 0) |
                                               (tgt.boiling_margin_degC.isna() & (tgt.P_bar > P_CRIT)))
    for w in wnames:
        m = tgt.well == w
        sub = tgt[m]
        iso = isothermal_flags(sub.elevation_mRL.to_numpy(float), sub.T_degC.to_numpy(float),
                               sub._open.to_numpy(bool))
        tgt.loc[m, "_iso"] = iso

    # --- natural-state tests (§4.1) ----------------------------------------------------------
    ns = natural_state_tests(wells, tgt, geo)
    conf = {w: ns[w]["confidence"] for w in wnames}
    verdict = {w: ns[w]["verdict"] for w in wnames}
    tgt["ns_confidence"] = tgt.well.map(conf)
    fz_quality_tags = sorted({f["q"] for w in wnames for f in wells[w]["fz"]})
    qf, qpresent = quality_factors(fz_quality_tags)

    # --- weights (5.5) -----------------------------------------------------------------------
    wt, wr, wp = [], [], []
    for _, r in tgt.iterrows():
        cf = CONF_F[r.ns_confidence]
        if r.above_model_surface:
            loc, why = LOC_F["above"], "above model surface"
        elif r._top_col or r._near_wh:
            loc, why = LOC_F["top"], "within 30 m of wellhead or top block of column"
        elif r.in_feedzone:
            loc, why = LOC_F["feed"], "feed zone"
        elif r.behind_casing:
            loc, why = LOC_F["casing"], "behind casing"
        elif r._iso:
            loc, why = LOC_F["iso"], "open hole in isothermal section (possible internal flow)"
        else:
            loc, why = LOC_F["open"], "open hole"
        wt.append(round(loc * cf, 6))
        wr.append(why)
        if r.in_feedzone and not r.above_model_surface:
            q = wells[r.well]["fz"][r._fz[0]]["q"]
            wp.append(round(qf[q] * cf, 6))
        else:
            wp.append(0.0)
    tgt["weight"] = wt
    tgt["weight_reason"] = wr
    tgt["weight_P"] = wp
    tgt["quality_tag"] = [wells[r.well]["fz"][r._fz[0]]["q"] if r.in_feedzone else ""
                          for _, r in tgt.iterrows()]

    # --- thermo: layer centres ---------------------------------------------------------------
    lc_rows = []
    for w in wnames:
        W = wells[w]
        zmin, zmax = W["tz"].min(), W["tz"].max()
        for b in wb[w]["blocks"]:
            L = geo["layers"].get(b["layer"])
            if L is None:
                continue
            c = L["centre"]
            if not (zmin - 1e-6 <= c <= zmax + 1e-6):
                continue
            T = interp_profile(W["tz"], W["tt"], c)
            P = interp_profile(W["pz"], W["pp"], c)
            Ts = tsat_c(P)
            lc_rows.append(dict(well=w, layer=b["layer"], layer_centre_mRL=c, block=b["block"],
                                T_degC=T, P_bar=P, Tsat_degC=Ts,
                                boiling_margin_degC=(Ts - T) if np.isfinite(Ts) else float("nan")))
    lc = pd.DataFrame(lc_rows)

    # --- feed zones --------------------------------------------------------------------------
    fz_rows = []
    for w in wnames:
        W = wells[w]
        sub_fz = list(W["fz"])
        fb = wb[w]["fzb"]
        order_fz = sorted(range(len(sub_fz)), key=lambda k: -(sub_fz[k]["top"] + sub_fz[k]["bot"]))
        for k, f in enumerate(sub_fz):
            top, bot = max(f["top"], f["bot"]), min(f["top"], f["bot"])
            mid = 0.5 * (top + bot)
            T = interp_profile(W["tz"], W["tt"], mid)
            P = interp_profile(W["pz"], W["pp"], mid)
            Ts = tsat_c(P)
            bmv = (Ts - T) if np.isfinite(Ts) else float("nan")
            blk = fb[k]["block"] if k < len(fb) else ""
            if k < len(fb):
                bl = next((bb for bb in wb[w]["blocks"] if bb["block"] == blk), None)
                if bl and not (bl["bot"] - 0.5 <= mid <= bl["top"] + 0.5):
                    say(f"WARNING {w} feed zone {k + 1}: block {blk} does not contain mid {mid}")
            cf = CONF_F[conf[w]]
            fz_rows.append(dict(well=w, feedzone=k + 1, top_mRL=top, bottom_mRL=bot, mid_mRL=mid,
                                quality=f["q"], T_degC=T, P_bar=P, Tsat_degC=Ts,
                                boiling_margin_degC=bmv, block=blk, ns_confidence=conf[w],
                                weight_T=round(LOC_F["feed"] * cf, 6),
                                weight_P=round(qf[f["q"]] * cf, 6), _thick=top - bot))
    fz = pd.DataFrame(fz_rows)

    # --- well summary (§4.4) -------------------------------------------------------------------
    ws_rows = []
    for w in wnames:
        g = tgt[tgt.well == w]
        ws_rows.append(summary_for(w, g, fz[fz.well == w], wells[w], ns[w], conf, verdict))
    ws = pd.DataFrame(ws_rows)

    # --- pressure consistency, hydrostatic and bands ------------------------------------------
    pc = ns_table(ns, wells)
    hy, hy_meta = hydrostatic(wells, tgt, conf)
    bands_df, band_list = band_table(tgt, wells, ns, wnames)
    geom_df = geometry_table(wells, wb, geo, feats)

    # --- write tables -----------------------------------------------------------------------
    tgt_out = tgt.copy()
    tgt_out["ns_confidence"] = tgt_out.ns_confidence
    cols_targets = ["well", "x", "y", "elevation_mRL", "T_degC", "P_bar", "Tsat_degC",
                    "boiling_margin_degC", "column", "layer", "block", "in_feedzone", "behind_casing",
                    "above_model_surface", "ns_confidence", "weight", "weight_reason", "weight_P",
                    "quality_tag", "model_surface_here_mRL", "near_wellhead_or_top_block",
                    "isothermal_open_hole", "gauge_or_survey"]
    tgt_out["near_wellhead_or_top_block"] = tgt._near_wh | tgt._top_col
    tgt_out["isothermal_open_hole"] = tgt._iso
    tgt_out["gauge_or_survey"] = "survey " + tgt.well.map(lambda w: wells[w]["survey_date"])
    tgt_out = tgt_out[cols_targets]
    for c in ["x", "y", "elevation_mRL", "T_degC", "P_bar", "Tsat_degC", "boiling_margin_degC",
              "model_surface_here_mRL"]:
        tgt_out[c] = tgt_out[c].map(lambda v: fnum(v, 3) if np.isfinite(v) else "")
    write_csv(tgt_out, os.path.join(CAL, "NS_calibration_targets.csv"))

    fz_out = fz.drop(columns=["_thick"])
    write_csv(round_df(fz_out, 3), os.path.join(CAL, "feedzones.csv"))
    lc_out = round_df(lc, 3)
    write_csv(lc_out, os.path.join(CAL, "T_P_at_layer_centres.csv"))
    write_csv(round_df(ws, 3), os.path.join(CAL, "well_summary.csv"))
    write_csv(round_df(pc, 3), os.path.join(ADD, "pressure_consistency_check.csv"))
    write_csv(round_df(hy, 4), os.path.join(ADD, "hydrostatic_check.csv"))
    write_csv(round_df(geom_df, 2), os.path.join(ADD, "well_geometry_check.csv"))
    write_csv(round_df(bands_df, 3), os.path.join(ADD, "T_P_at_elevation_bands.csv"))
    params = parameters_table()
    write_csv(params, os.path.join(ADD, "analysis_parameters.csv"))
    write_csv(pd.DataFrame(hy_meta["component_rows"]), os.path.join(ADD, "hydrostatic_cold_component.csv"))

    # --- figures ------------------------------------------------------------------------------
    ctx = dict(wells=wells, wnames=wnames, wb=wb, geo=geo, polys=polys, feats=feats, tgt=tgt,
               fz=fz, ws=ws, pc=pc, ns=ns, conf=conf, verdict=verdict, well_colour=well_colour,
               extent=(xmin, xmax, ymin, ymax), lc=lc, hy=hy, bands=bands_df, band_list=band_list,
               trend=ns["_trend"], qf=qf, qpresent=qpresent, groups=groups)
    figures_all(ctx)

    # --- figure index, supporting tables and report -------------------------------------------
    figure_index(ctx)
    extra_tables(ctx)
    write_report(ctx)
    say(f"Wrote {len(tgt_out)} targets, {len(fz)} feed zones, {len(lc)} layer rows, {len(ws)} wells")
    return log


def round_df(df, nd):
    out = df.copy()
    for c in out.columns:
        if out[c].dtype.kind == "f":
            out[c] = out[c].round(nd)
    return out


def write_inventory():
    T = "Baseline/Teaching_Model_3D"
    rows = [
        ("well_data_json", f"{T}/Data/Data_json", "10 per-well JSON files: track, casing/liner design, feed zones, downhole T and P surveys, transient pressures, production/injection and enthalpy records."),
        ("well_data_dat", f"{T}/Data/Data_dat", "Plain-column copies: survey files (value, elevation), transient files (time in years, value), flow and enthalpy files. Used only for the T, P and transient cross-check. Irregular prefix Wel10 for WELL-10."),
        ("mapping", f"{T}/Data/Mapping/Mapping.json", "Faults (Fault 1-3, polylines) and water bodies (Lake Tarawera fill; Tarawera River and river as multi-segment lines)."),
        ("well_groups", f"{T}/Data/Groups_data.json", "One group 'All' (Areas/Synthetic) listing WELL-1 to WELL-10."),
        ("geometry_dat", f"{T}/Models/TE14515/gTE14515.dat", "MULgraph grid. Used: VERTICES, GRID (column polygons), LAYERS, SURFA. WELLS section not used (tracks come from the well JSON)."),
        ("well_blocks_json", f"{T}/Models/TE14515/gTE14515_wells.json", "Well-to-block intersections and interpolation weights. Authority for block assignment."),
        ("reference_rock", f"{T}/Models/TE14515/gTE14515_dict_rocktype.json", "Reference rock-type dictionary. Not used as evidence."),
        ("reference_bc", f"{T}/Models/TE14515/gTE14515_boundary_conditions.json", "Reference boundary conditions. Not used as evidence."),
        ("reference_results", f"{T}/Models/TE14515/NaturalState/TE14515_NS_000.dat", "Reference model natural-state input file (AUTOUGH2 format). Not used as evidence; not read beyond its header."),
        ("reference_results", f"{T}/Models/TE14515/NaturalState/TE14515_NS_000.incon", "Reference model initial conditions. Not used as evidence."),
        ("reference_other", f"{T}/Models/TE14515/gTE14515.auxjson", "Scene/GIS configuration of the reference model; holds hard-coded source paths from another machine. Not used."),
        ("reference_other", f"{T}/Models/TE14515/currentGeomSetup.json", "Model set-up file. Not used."),
        ("reference_other", f"{T}/Models/TE14515/NaturalState/currentNSSetup.json", "Natural-state set-up file. Not used."),
        ("reference_other", f"{T}/Data/Plotting", "Plotting configuration files; contain hard-coded paths from another machine. Not used."),
        ("out_of_scope", f"{T}/Models/TE14515/Production/currentPRSetup.json", "Production model set-up. Not read."),
        ("out_of_scope", f"{T}/Models/TE14515/FutureScenarios", "Future-scenario model set-up. Not read."),
        ("out_of_scope", f"{T}/Parameters/FutureScenarios", "Future-scenario parameters (empty folder tree in the baseline). Not read."),
        ("out_of_scope", f"{T}/Parameters/ProductionHistory", "Production and reinjection history spreadsheets and a .DS_Store file. Not read; production timing taken from the well JSON flow records."),
        ("out_of_scope", f"{T}/Data/Area_data.json", "Field-total steam and power history. Not read."),
        ("out_of_scope", f"{T}/Data/Make_up_json", "Empty directory (no files). Not read."),
    ]
    for r in rows:
        if not os.path.exists(os.path.join(ROOT, r[1])):
            raise FileNotFoundError(r[1])
    df = pd.DataFrame(rows, columns=["role", "path", "note"])
    # cover every file under Baseline
    covered = []
    for role, p, _ in rows:
        full = os.path.join(ROOT, p)
        covered.append(os.path.normpath(full))
    uncovered = []
    for dp, _, fs in os.walk(BASE):
        for f in fs:
            fp = os.path.normpath(os.path.join(dp, f))
            if not any(fp == c or fp.startswith(c + os.sep) for c in covered):
                uncovered.append(fp)
    assert not uncovered, f"uncovered baseline files: {uncovered[:5]}"
    write_csv(df, os.path.join(ADD, "input_inventory.csv"))


def dat_cross_check(wells):
    rows = []
    folder = os.path.join(BASE, "Teaching_Model_3D", "Data", "Data_dat")
    for fn, data in load_dat_pairs(folder):
        base = os.path.basename(fn)
        m = re.match(r"^(Well\d+|Wel10)_(.*)\.dat$", base)
        if not m:
            continue
        wname = "WELL-" + re.sub(r"\D", "", m.group(1))
        W = wells.get(wname)
        kind = m.group(2)
        if W is None:
            continue
        if kind.startswith("downhole_temperature") or kind.startswith("downhole_pressure"):
            is_t = "temperature" in kind
            vals, zs = data[:, 0], data[:, 1]
            ref_z = W["tz"] if is_t else W["pz"]
            ref_v = W["tt"] if is_t else W["pp"]
            d = []
            for z, v in zip(zs, vals):
                jj = np.where(np.isclose(ref_z, z, atol=1e-6))[0]
                d.append(abs(ref_v[jj[0]] - v) if len(jj) else np.inf)
            ok = len(zs) == len(ref_z) and max(d) < 1e-6
            rows.append(dict(file=f"Data/Data_dat/{base}", well=wname,
                             series="survey T" if is_t else "survey P", n_points=len(zs),
                             status="match" if ok else "mismatch",
                             max_abs_diff=round(float(max(d)), 6) if d else ""))
        elif kind.startswith(("transient_pressure", "-")) or "transient_pressure" in kind:
            el = float(re.search(r"(-?\d+\.\d+)", kind).group(1))
            ser = W["gauges"].get(el, [])
            t_years, vals = data[:, 0], data[:, 1]
            ref0 = dtm.date(2014, 1, 1)
            d = []
            for t, v in zip(t_years, vals):
                dd = ref0 + dtm.timedelta(days=round(t * 365.25))
                hit = [val for ds, val in ser if abs((dtm.date.fromisoformat(ds) - dd).days) <= 2]
                d.append(abs(hit[0] - v) if hit else np.inf)
            ok = len(ser) == len(t_years) and max(d) < 1e-6
            rows.append(dict(file=f"Data/Data_dat/{base}", well=wname, series=f"transient P at {el:g} mRL",
                             n_points=len(t_years), status="match" if ok else "mismatch",
                             max_abs_diff=round(float(max(d)), 6) if d else ""))
        else:
            rows.append(dict(file=f"Data/Data_dat/{base}", well=wname, series=kind, n_points=len(data),
                             status="not checked (flow or enthalpy series; not part of 4.1 step 1)",
                             max_abs_diff=""))
    return rows


def natural_state_tests(wells, tgt, geo):
    """4.1: gauge test, feed-zone trend test and verdict rubric."""
    out = {}
    prec = GAUGE_PREC
    # gauge tests
    for w, W in wells.items():
        g = {}
        if W["gauges"]:
            el, series = next(iter(W["gauges"].items()))
            first_date, first_val = series[0]
            sv = interp_profile(W["pz"], W["pp"], el)
            diff = sv - first_val
            own = W["first_flow"]
            pre = (own is None) or (first_date <= own)
            passed = pre and abs(diff) <= prec + 1e-9
            g = dict(gauge_elev=el, survey_P=sv, first_date=first_date, first_P=first_val,
                     diff=diff, own_first_flow=own, series_start=first_date, pre=pre,
                     passed=passed, n=len(series), last_date=series[-1][0], series=series)
        out[w] = dict(gauge=g)
    # gauge-verified wells -> trend from their feed zones
    verified = [w for w in wells if out[w]["gauge"] and out[w]["gauge"]["passed"]]
    pts = []
    for w in verified:
        W = wells[w]
        for f in W["fz"]:
            mid = 0.5 * (f["top"] + f["bot"])
            P = interp_profile(W["pz"], W["pp"], mid)
            pts.append((-mid, P))
    if len(pts) >= 2:
        A = np.array([[1.0, p[0]] for p in pts])
        y = np.array([p[1] for p in pts])
        coef = np.linalg.lstsq(A, y, rcond=None)[0]
        resid = y - A @ coef
        rms = float(np.sqrt(np.mean(resid ** 2)))
    else:
        coef = np.array([float("nan"), float("nan")])
        rms = float("nan")
    tol = max(TREND_TOL_MIN, TREND_TOL_RMS * rms) if np.isfinite(rms) else TREND_TOL_MIN
    trend = dict(a=float(coef[0]), b=float(coef[1]), rms=rms, tol=tol, verified=verified,
                 n_points=len(pts))

    def trend_p(z):
        return trend["a"] + trend["b"] * (-z)

    for w, W in wells.items():
        resids = []
        fzres = []
        for f in W["fz"]:
            mid = 0.5 * (f["top"] + f["bot"])
            P = interp_profile(W["pz"], W["pp"], mid)
            r = P - trend_p(mid)
            fzres.append(r)
            resids.append(r)
        trend_ok = all(abs(r) <= tol for r in resids) if resids else False
        out[w]["fz_resid"] = fzres
        out[w]["trend_ok"] = trend_ok
        out[w]["max_resid"] = max(abs(r) for r in resids) if resids else float("nan")
        # injector formation vs injectate
        injector = W["kind"] == "injector"
        h_inj = None
        if injector and W["enth"]:
            h_inj = float(np.median([v for _, v in W["enth"]]))
        inj_T = injectate_t(h_inj) if h_inj is not None else float("nan")
        fz_T = [interp_profile(W["tz"], W["tt"], 0.5 * (f["top"] + f["bot"])) for f in W["fz"]]
        form_T = float(np.median(fz_T)) if fz_T else float("nan")
        hotter = bool(injector and np.isfinite(inj_T) and np.isfinite(form_T) and form_T > inj_T)
        out[w]["injector"] = injector
        out[w]["inj_h"] = h_inj
        out[w]["inj_T"] = inj_T
        out[w]["form_T"] = form_T
        out[w]["hotter"] = hotter
        g = out[w]["gauge"]
        if g and g["passed"]:
            v, c = "natural state", "high"
            rule = "rule 1: survey matches a pre-production gauge reading within the reading precision"
        elif hotter:
            v, c = "probably natural state at best", "low"
            rule = "rule 2: injector whose formation is hotter than the injectate; cooling could be hidden"
        elif g and g["pre"] and trend_ok:
            v, c = "probably natural state", "medium"
            rule = "rule 3: gauge test misses by more than the precision; feed zones on the verified trend"
        elif trend_ok:
            v, c = "probably natural state", "medium"
            rule = "rule 4: no pre-production gauge reading; feed zones on the verified trend"
        else:
            v, c = "possibly disturbed", "low"
            rule = "rule 5: feed-zone pressure off the verified trend by more than the tolerance"
        out[w]["verdict"] = v
        out[w]["confidence"] = c
        out[w]["rule"] = rule
    out["_trend"] = trend
    return out


def ns_table(ns, wells):
    rows = []
    trend = ns["_trend"]
    for w in wells:
        g = ns[w]["gauge"]
        if g:
            rows.append(dict(well=w, gauge_elev_mRL=g["gauge_elev"], survey_P_at_gauge_bar=g["survey_P"],
                             first_gauge_date=g["first_date"], first_gauge_P_bar=g["first_P"],
                             diff_survey_minus_first_bar=g["diff"], verdict=ns[w]["verdict"],
                             confidence=ns[w]["confidence"],
                             own_first_flow_date=g["own_first_flow"] or "",
                             gauge_series_start_vs_own_flow=("before or at own first flow" if g["pre"]
                                                              else "after own first flow"),
                             pre_production_reading=g["pre"],
                             gauge_test=("pass" if g["passed"] else "miss"),
                             last_gauge_date=g["last_date"],
                             trend_max_abs_resid_bar=ns[w]["max_resid"], trend_tol_bar=trend["tol"],
                             rule=ns[w]["rule"]))
        else:
            rows.append(dict(well=w, gauge_elev_mRL=np.nan, survey_P_at_gauge_bar=np.nan,
                             first_gauge_date="", first_gauge_P_bar=np.nan,
                             diff_survey_minus_first_bar=np.nan, verdict=ns[w]["verdict"],
                             confidence=ns[w]["confidence"], own_first_flow_date=wells[w]["first_flow"] or "",
                             gauge_series_start_vs_own_flow="no gauge", pre_production_reading=False,
                             gauge_test="no gauge", last_gauge_date="",
                             trend_max_abs_resid_bar=ns[w]["max_resid"], trend_tol_bar=trend["tol"],
                             rule=ns[w]["rule"]))
    df = pd.DataFrame(rows)
    df["gauge_P_at_survey_date_bar"] = [
        gauge_at_date(ns[w]["gauge"]["series"], wells[w]["survey_date"]) if ns[w]["gauge"] else np.nan
        for w in df.well]
    df["gauge_P_at_survey_date_bar"] = df.gauge_P_at_survey_date_bar.map(lambda v: "" if pd.isna(v) else round(v, 3))
    df["gauge_elev_mRL"] = df.gauge_elev_mRL.map(lambda v: "" if pd.isna(v) else v)
    df["survey_P_at_gauge_bar"] = df.survey_P_at_gauge_bar.map(lambda v: "" if pd.isna(v) else round(v, 3))
    df["first_gauge_P_bar"] = df.first_gauge_P_bar.map(lambda v: "" if pd.isna(v) else v)
    df["diff_survey_minus_first_bar"] = df.diff_survey_minus_first_bar.map(lambda v: "" if pd.isna(v) else round(v, 3))
    return df


def hydrostatic(wells, tgt, conf):
    rows = []
    comp_rows = []
    per_cold = {"cased": [], "open": []}
    cold_wells = []
    for w in wells:
        W = wells[w]
        g = tgt[tgt.well == w].sort_values("elevation_mRL", ascending=False).reset_index(drop=True)
        z = g.elevation_mRL.to_numpy(float)
        t = g.T_degC.to_numpy(float)
        p = g.P_bar.to_numpy(float)
        bm = g.boiling_margin_degC.to_numpy(float)
        liq = (p >= HYDRO_P) & ((bm >= HYDRO_BM) | (np.isnan(bm) & (p > P_CRIT)))
        bounds = [f["top"] for f in W["fz"]] + [f["bot"] for f in W["fz"]] + [W["shoe"]]
        wins = windows(z, liq, bounds)
        is_cold = t.max() < COLD_T
        if is_cold:
            cold_wells.append(w)
        rho_cache = {}

        def rho_at(k):
            if k not in rho_cache:
                rho_cache[k] = rho_liquid(t[k], p[k])
            return rho_cache[k]

        for (s, e) in wins:
            zt, zb = z[s], z[e]
            length = zt - zb
            dP_meas = p[e] - p[s]
            dP_if = 0.0
            for k in range(s, e):
                dP_if += 0.5 * (rho_at(k) + rho_at(k + 1)) * G * (z[k] - z[k + 1])   # Pa
            dP_if = dP_if / 1e5
            straddle = bool(zt > W["shoe"] > zb)
            if straddle:
                band = "straddles shoe"          # neither band (§4.3): top above and bottom below the shoe
            else:
                band = "cased" if zb >= W["shoe"] else "open"
            ratio = dP_meas / dP_if
            rho_imp = dP_meas * 1e5 / (G * length)
            unc = 0.1 / dP_if          # rounding: +/-0.1 bar on the difference of two 0.1-bar readings
            rows.append(dict(well=w, z_top_mRL=zt, z_bottom_mRL=zb, length_m=length, band=band,
                             dP_measured_bar=dP_meas, dP_IF97_bar=dP_if, ratio_measured_to_IF97=ratio,
                             rho_implied_kgm3=rho_imp, ratio_rounding_uncertainty=unc,
                             n_points=e - s + 1, cold_well=is_cold,
                             straddles_shoe=straddle))
            if is_cold and band in per_cold:
                per_cold[band].append(ratio)
    comp = {b: (float(np.median(v)) if v else float("nan")) for b, v in per_cold.items()}
    for b in comp:
        comp_rows.append(dict(band=b, n_cold_windows=len(per_cold[b]), cold_wells=",".join(cold_wells),
                              median_ratio=comp[b]))
    df = pd.DataFrame(rows)
    df["excess_pct_over_cold_component"] = [
        (100.0 * (r.ratio_measured_to_IF97 / comp[r.band] - 1.0)) if (r.band in comp and np.isfinite(comp[r.band])) else np.nan
        for _, r in df.iterrows()]
    return df, dict(component=comp, cold_wells=cold_wells, component_rows=comp_rows)


def windows(z, liq, bounds):
    """4.3 window algorithm. z descending; liq boolean per point; bounds = feed-zone and shoe elevations."""
    wins = []
    n = len(z)
    i = 0
    while i < n:
        if not liq[i]:
            i += 1
            continue
        j = i
        while j + 1 < n and liq[j + 1]:
            j += 1
        wins += window_run(z, i, j, bounds)
        i = j + 1
    return wins


def window_run(z, i, j, bounds):
    if z[i] - z[j] < WIN_MIN:
        return []
    wins = []
    s = i
    k = i + 1
    while k <= j:
        reached = any(z[k] - 1e-9 <= B <= z[k - 1] + 1e-9 for B in bounds)
        if z[s] - z[k - 1] >= WIN_MIN and reached:
            wins.append((s, k - 1))
            s = k
        k += 1
    if z[s] - z[j] >= WIN_MIN:
        wins.append((s, j))
    elif wins:
        wins[-1] = (wins[-1][0], j)
    return wins


def summary_for(w, g, fzw, W, nsw, conf, verdict):
    """Per-well characteristics (§4.4); mirrors the definitions used by the acceptance check."""
    g = g.sort_values("elevation_mRL", ascending=False).groupby("elevation_mRL", sort=False).mean(numeric_only=True).reset_index()
    z = g.elevation_mRL.to_numpy(float)
    t = g.T_degC.to_numpy(float)
    p = g.P_bar.to_numpy(float)
    bm = g.boiling_margin_degC.to_numpy(float)
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
    inv = best if best[0] >= INVERSION_MIN else (np.nan, np.nan, np.nan)
    e["inversion_degC"], e["inversion_top_mRL"], e["inversion_base_mRL"] = inv
    liq = (p >= HYDRO_P) & ((bm >= 0) | (np.isnan(bm) & (p > P_CRIT)))
    zl, tl, pl, bml = z[liq], t[liq], p[liq], bm[liq]
    e["shallow_gradient_degC_per_km"] = np.nan
    n_sh = 0
    if len(zl):
        wsel = zl >= zl[0] - SHALLOW_WIN
        n_sh = int(wsel.sum())
        if n_sh >= 3:
            e["shallow_gradient_degC_per_km"] = -np.polyfit(zl[wsel], tl[wsel], 1)[0] * 1000
    fzb = fzw.bottom_mRL.min() if len(fzw) else np.nan
    bounds = [b for b in (inv[2], fzb) if not np.isnan(b)]
    sel = z <= min(bounds) if bounds else np.ones_like(z, bool)
    e["deep_gradient_degC_per_km"] = np.nan
    deep_span = float(z[sel].max() - z[sel].min()) if sel.sum() else 0.0
    if sel.sum() >= 3 and deep_span >= DEEP_SPAN:
        e["deep_gradient_degC_per_km"] = -np.polyfit(z[sel], t[sel], 1)[0] * 1000
    e["z_100C_mRL"] = crossing(z, t, 100.0)
    e["z_200C_mRL"] = crossing(z, t, 200.0)
    ok = ~np.isnan(bml)
    if ok.any():
        j = int(np.argmin(np.where(ok, bml, np.inf)))
        e["min_boiling_margin_degC"], e["min_boiling_margin_elev_mRL"] = bml[j], zl[j]
    else:
        e["min_boiling_margin_degC"] = e["min_boiling_margin_elev_mRL"] = np.nan
    e["water_level_mRL"] = np.nan
    wl_ext = np.nan
    if len(zl) >= 3:
        b, a = np.polyfit(zl[:3], pl[:3], 1)
        if b != 0:
            e["water_level_mRL"] = (P_ATM - a) / b
            wl_ext = e["water_level_mRL"] - zl[0]
    # pivot (5.x / 4.4)
    pv = fzw.copy()
    cfw = CONF_F[nsw["confidence"]] if False else 1.0
    pivot_fz = None
    if len(pv):
        qmax = pv.weight_P.max()
        cand = pv[pv.weight_P == qmax].copy()
        cand["th"] = cand._thick
        cand = cand.sort_values(["th", "mid_mRL"], ascending=[False, True])
        pivot_fz = cand.iloc[0]
        e["pivot_elev_mRL"] = pivot_fz.mid_mRL
        e["pivot_P_bar"] = pivot_fz.P_bar
        e["pivot_feedzone"] = int(pivot_fz.feedzone)
    else:
        e["pivot_elev_mRL"] = np.nan
        e["pivot_P_bar"] = np.nan
        e["pivot_feedzone"] = np.nan
    # regime (§4.4)
    if ok.any() and np.nanmin(bml) <= TWO_PHASE_BM:
        reg = "two-phase"
    elif e["Tmax_degC"] < COLD_T:
        reg = "recharge/downflow"
    elif not np.isnan(inv[0]) and inv[0] >= OUTFLOW_MIN:
        reg = "outflow"
    elif e["Tmax_degC"] >= 200 and np.isnan(inv[0]):
        reg = "upflow"
    else:
        reg = "conductive"
    near_boil = bool(ok.any() and np.nanmin(bml) < NEAR_BOIL)
    note = ""
    if reg == "upflow" and near_boil:
        note = "upflow by the Tmax rule; minimum boiling margin is below 5 degC, so the top of the upflow is near boiling"
    elif reg == "recharge/downflow":
        note = "cold profile (Tmax below 100 degC); the regime does not by itself show direction of flow"
    elif reg == "outflow":
        note = "inversion present, so the regime is outflow; the inversion is measured, its cause is not"
    return {
        "well": w,
        "Tmax_degC": e["Tmax_degC"], "Tmax_elev_mRL": e["Tmax_elev_mRL"],
        "T_at_TD_degC": e["T_at_TD_degC"], "TD_mRL": e["TD_mRL"],
        "inversion_degC": e["inversion_degC"], "inversion_top_mRL": e["inversion_top_mRL"],
        "inversion_base_mRL": e["inversion_base_mRL"],
        "shallow_gradient_degC_per_km": e["shallow_gradient_degC_per_km"],
        "deep_gradient_degC_per_km": e["deep_gradient_degC_per_km"],
        "z_100C_mRL": e["z_100C_mRL"], "z_200C_mRL": e["z_200C_mRL"],
        "min_boiling_margin_degC": e["min_boiling_margin_degC"],
        "min_boiling_margin_elev_mRL": e["min_boiling_margin_elev_mRL"],
        "water_level_mRL": e["water_level_mRL"], "pivot_elev_mRL": e["pivot_elev_mRL"],
        "pivot_P_bar": e["pivot_P_bar"], "regime": reg,
        "ns_verdict": verdict[w], "ns_confidence": conf[w],
        "deep_fit_n": int(sel.sum()) if sel.sum() >= 3 and deep_span >= DEEP_SPAN else "",
        "deep_fit_span_m": deep_span if sel.sum() >= 3 and deep_span >= DEEP_SPAN else "",
        "shallow_fit_n": n_sh,
        "water_level_extrapolation_m": wl_ext,
        "near_boiling": near_boil,
        "pivot_feedzone": e["pivot_feedzone"],
        "regime_note": note,
    }


def crossing(z, t, x):
    for i in range(len(t) - 1):
        if t[i] < x <= t[i + 1]:
            return z[i] + (x - t[i]) * (z[i + 1] - z[i]) / (t[i + 1] - t[i])
    return float("nan")


def band_edges(tgt):
    zmax = tgt.elevation_mRL.max()
    zmin = tgt.elevation_mRL.min()
    top0 = math.ceil(zmax / BAND_M) * BAND_M
    bots = []
    k = 0
    while top0 - (k + 1) * BAND_M >= math.floor(zmin / BAND_M) * BAND_M:
        bots.append((top0 - k * BAND_M, top0 - (k + 1) * BAND_M))
        k += 1
    return bots


def band_table(tgt, wells, ns, wnames):
    trend = ns["_trend"]
    edges = band_edges(tgt)
    counts = []
    for top, bot in edges:
        ws_in = tgt[(tgt.elevation_mRL <= top) & (tgt.elevation_mRL > bot)].well.nunique()
        counts.append(ws_in)
    idx = [i for i, c in enumerate(counts) if c >= BAND_MIN_WELLS]
    sel = edges[idx[0]:idx[-1] + 1] if idx else []
    rows = []
    band_list = []
    for top, bot in sel:
        sub = tgt[(tgt.elevation_mRL <= top) & (tgt.elevation_mRL > bot)]
        well_means = {}
        for w in wnames:
            ww = sub[sub.well == w]
            if len(ww):
                well_means[w] = ww.T_degC.mean()
        med = float(np.median(list(well_means.values()))) if well_means else float("nan")
        band_list.append((top, bot, med, len(well_means)))
        for w in wnames:
            ww = sub[sub.well == w]
            if not len(ww):
                continue
            tm = ww.T_degC.mean()
            pdev = np.mean([pp - (trend["a"] + trend["b"] * (-zz)) for pp, zz in zip(ww.P_bar, ww.elevation_mRL)])
            rows.append(dict(well=w, band_top_mRL=top, band_bottom_mRL=bot, n_points=len(ww),
                             T_mean_degC=tm, T_field_median_degC=med, T_minus_field_median_degC=tm - med,
                             P_mean_bar=ww.P_bar.mean(), P_minus_trend_bar=pdev))
    return pd.DataFrame(rows), band_list


def geometry_table(wells, wb, geo, feats):
    rows = []
    polys = build_column_polys(geo)
    faults = [f for f in feats if f["kind"] == "fault"]
    water = [f for f in feats if f["kind"] == "water"]
    for w in wells:
        W = wells[w]
        wh = W["track"][0]
        blocks = wb[w]["blocks"]
        top_block = max(blocks, key=lambda b: b["top"])["block"] if blocks else ""
        col = re.match(r"^([a-z]+)", top_block).group(1) if top_block else ""
        msurf = geo["surf"].get(col, geo["grid_top"])
        dev = float(np.max(np.hypot(W["track"][:, 0] - wh[0], W["track"][:, 1] - wh[1])))
        row = dict(well=w, wellhead_mRL=W["wellhead_z"], model_surface_mRL=msurf,
                   wellhead_minus_surface_m=W["wellhead_z"] - msurf,
                   max_horizontal_deviation_m=dev)
        dists = {}
        for f in feats:
            if f["kind"] == "fault":
                d = min(dist_to_polyline(wh[0], wh[1], s) for s in f["segs"])
            else:
                d = min(dist_to_polyline(wh[0], wh[1], s) for s in f["segs"])
                if f["fill"]:
                    pth = Path(f["segs"][0])
                    if pth.contains_point((wh[0], wh[1])):
                        d = 0.0
            dists[f["name"]] = d
        fn = min(faults, key=lambda f: dists[f["name"]])
        wn = min(water, key=lambda f: dists[f["name"]])
        row["nearest_fault"] = fn["name"]
        row["dist_nearest_fault_m"] = dists[fn["name"]]
        row["nearest_water_body"] = wn["name"]
        row["dist_nearest_water_body_m"] = dists[wn["name"]]
        for f in feats:
            row[f"dist_{re.sub(r'[^A-Za-z0-9]+', '_', f['name']).strip('_')}_m"] = dists[f["name"]]
        rows.append(row)
    return pd.DataFrame(rows)


def parameters_table():
    rows = [
        ("gravity", G, "m/s2", "CLAUDE 4.3"),
        ("atmospheric pressure (water-level reference)", P_ATM, "bar", "CLAUDE 4.4"),
        ("critical pressure (IF97)", P_CRIT, "bar", "CLAUDE 4.3"),
        ("critical temperature (IF97)", T_CRIT, "degC", "IF97"),
        ("minimum hydrostatic window length", WIN_MIN, "m", "CLAUDE 4.3"),
        ("isothermal section range", ISO_RANGE, "degC", "CLAUDE 4.2"),
        ("isothermal section span", ISO_SPAN, "m", "CLAUDE 4.2"),
        ("single-phase test boiling margin", HYDRO_BM, "degC", "CLAUDE 4.3"),
        ("single-phase test pressure", HYDRO_P, "bar", "CLAUDE 4.3"),
        ("boiling-margin flag (near boiling)", NEAR_BOIL, "degC", "CLAUDE 4.3, 5.1"),
        ("two-phase threshold", TWO_PHASE_BM, "degC", "CLAUDE 4.4"),
        ("inversion threshold", INVERSION_MIN, "degC", "CLAUDE 4.4"),
        ("outflow threshold", OUTFLOW_MIN, "degC", "CLAUDE 4.4"),
        ("shallow gradient window", SHALLOW_WIN, "m", "CLAUDE 4.4"),
        ("deep gradient minimum span", DEEP_SPAN, "m", "CLAUDE 4.4"),
        ("cold-well threshold (Tmax)", COLD_T, "degC", "CLAUDE 4.3"),
        ("well-head weight distance", WELLHEAD_NEAR, "m", "CLAUDE 5.5"),
        ("gauge reading precision", GAUGE_PREC, "bar", "assumption: one unit of the 0.1 bar record"),
        ("trend tolerance floor", TREND_TOL_MIN, "bar", "CLAUDE 4.1"),
        ("trend tolerance, multiple of RMS", TREND_TOL_RMS, "-", "CLAUDE 4.1"),
        ("well distance from section line", TREND_NEAR, "m", "CLAUDE 8.2 F4"),
        ("elevation band thickness", BAND_M, "m", "CLAUDE 8.2 F3"),
        ("wells needed to bound the band range", BAND_MIN_WELLS, "wells", "CLAUDE 8.2 F3"),
        ("pressure for injectate enthalpy-to-temperature", INJ_P, "bar", "assumption: subcooled liquid, IF97 Region 1"),
        ("weight: wellhead / top block", LOC_F["top"], "-", "CLAUDE 5.5"),
        ("weight: feed zone", LOC_F["feed"], "-", "CLAUDE 5.5"),
        ("weight: behind casing", LOC_F["casing"], "-", "CLAUDE 5.5"),
        ("weight: open hole", LOC_F["open"], "-", "CLAUDE 5.5"),
        ("weight: open hole, isothermal section", LOC_F["iso"], "-", "CLAUDE 5.5"),
        ("confidence factor high", CONF_F["high"], "-", "CLAUDE 5.5"),
        ("confidence factor medium", CONF_F["medium"], "-", "CLAUDE 5.5"),
        ("confidence factor low", CONF_F["low"], "-", "CLAUDE 5.5"),
        ("temperature colour step", T_MAP_STEP, "degC", "CLAUDE 8.1"),
        ("pressure-deviation colour step", P_DEV_STEP, "bar", "CLAUDE 8.1"),
        ("isotherm spacing", 50.0, "degC", "CLAUDE 8.1"),
        ("isotherm heavy levels", 100.0, "degC", "CLAUDE 8.1"),
        ("isotherm heavy levels", 200.0, "degC", "CLAUDE 8.1"),
        ("elevation axis rounding", 100.0, "m", "CLAUDE 8.1"),
    ]
    return pd.DataFrame(rows, columns=["parameter", "value", "unit", "source"])


# ----------------------------------------------------------------------------------------
# Summary table of wells (used by report and figures)
# ----------------------------------------------------------------------------------------
ELEV_LIM = (-2900.0, 400.0)      # highest wellhead / deepest survey, rounded outward to 100 m (8.1)
T_NORM = BoundaryNorm(np.arange(0.0, 300.0 + T_MAP_STEP, T_MAP_STEP), 30)
T_CMAP = plt.get_cmap("inferno", 30)
P_NORM = BoundaryNorm(np.arange(-10.0, 10.0 + P_DEV_STEP, P_DEV_STEP), 20)
P_CMAP = plt.get_cmap("RdBu_r", 20)
ISO_LW = {50.0: 0.7, 100.0: 1.8, 150.0: 0.7, 200.0: 1.8, 250.0: 0.7, 300.0: 0.7}


def elev_axis_lim(ctx):
    return ELEV_LIM


def trend_p(ctx, z):
    tr = ctx["trend"]
    return tr["a"] + tr["b"] * (-np.asarray(z, float))


def well_surface(ctx, w):
    blocks = ctx["wb"][w]["blocks"]
    top_block = max(blocks, key=lambda b: b["top"])
    return ctx["geo"]["surf"].get(top_block["col"], ctx["geo"]["grid_top"])


def plot_base(ax, ctx, highlight=None, zoom=None, labels=True, fs=6, lines=None):
    geo = ctx["geo"]
    segs = []
    for _, vn in geo["cols"]:
        pts = [geo["verts"][v] for v in vn if v in geo["verts"]]
        if len(pts) >= 3:
            segs.append(np.array(pts + [pts[0]]))
    ax.add_collection(LineCollection(segs, colors="#c4c4c4", linewidths=0.35, zorder=1))
    xmin, xmax, ymin, ymax = ctx["extent"]
    ax.add_patch(Rectangle((xmin, ymin), xmax - xmin, ymax - ymin, fill=False, ls="--",
                           ec="0.35", lw=0.8, zorder=2))
    for f in ctx["feats"]:
        if f["kind"] == "water":
            for s in f["segs"]:
                if f["fill"]:
                    ax.fill(s[:, 0], s[:, 1], color="#9ecae1", alpha=0.8, lw=0.6, ec="#2171b5", zorder=2)
                else:
                    ax.plot(s[:, 0], s[:, 1], color="#2171b5", lw=1.3, zorder=2)
            if labels:
                s = f["segs"][0]
                c = np.clip(s.mean(axis=0), [xmin + 200, ymin + 200], [xmax - 200, ymax - 200]) if f["fill"] else s[np.argmin(np.abs(s[:, 0] - (xmin + xmax) / 2) + np.abs(s[:, 1] - (ymin + ymax) / 2))]
                ax.text(c[0], c[1], f["name"], fontsize=fs, color="#08519c", zorder=6, ha="center")
        else:
            for s in f["segs"]:
                ax.plot(s[:, 0], s[:, 1], color="k", lw=1.2, zorder=3)
            if labels:
                s = f["segs"][0]
                c = s[np.argmin(np.abs(s[:, 0] - (xmin + xmax) / 2) + np.abs(s[:, 1] - (ymin + ymax) / 2))]
                ax.text(c[0], c[1], f["name"], fontsize=fs, color="k", zorder=6, ha="left",
                        va="bottom")
    for w in ctx["wnames"]:
        tr = ctx["wells"][w]["track"]
        if highlight is None:
            colr, lw = "0.55", 1.0
        elif w == highlight:
            colr, lw = ctx["well_colour"][w], 2.8
        else:
            colr, lw = "0.7", 0.8
        ax.plot(tr[:, 0], tr[:, 1], color=colr, lw=lw, zorder=4)
        ax.plot(tr[0, 0], tr[0, 1], "o", color="k", ms=3.2, zorder=7)
        if labels:
            ax.text(tr[0, 0] + 60, tr[0, 1] + 60, w, fontsize=fs, zorder=8,
                    color=ctx["well_colour"][w] if highlight in (None, w) else "0.3", fontweight="bold")
    if lines:
        for name, (p0, p1) in lines.items():
            ax.plot([p0[0], p1[0]], [p0[1], p1[1]], ls=(0, (5, 3)), color="#d62728", lw=1.2, zorder=5)
            ax.text(p0[0], p0[1], name.split("|")[0], color="#d62728", fontsize=fs + 1, zorder=9, ha="right")
            ax.text(p1[0], p1[1], name.split("|")[1], color="#d62728", fontsize=fs + 1, zorder=9, ha="left")
    if zoom is None:
        ax.set_xlim(xmin - 100, xmax + 100)
        ax.set_ylim(ymin - 100, ymax + 100)
    else:
        ax.set_xlim(zoom[0], zoom[1])
        ax.set_ylim(zoom[2], zoom[3])
    ax.set_aspect("equal")
    ax.set_xlabel("Easting (m, NZTM)", fontsize=8)
    ax.set_ylabel("Northing (m, NZTM)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.ticklabel_format(style="plain", useOffset=False)
    x0, y0 = ax.get_xlim()[0] + 0.05 * np.diff(ax.get_xlim())[0], ax.get_ylim()[0] + 0.05 * np.diff(ax.get_ylim())[0]
    ax.plot([x0, x0 + 1000], [y0, y0], "k-", lw=2.2, zorder=9)
    ax.text(x0 + 500, y0 + 0.012 * np.diff(ax.get_ylim())[0], "1 km", fontsize=6, ha="center", zorder=9)


def base_legend_handles():
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    return [Line2D([], [], color="#c4c4c4", lw=0.8, label="model grid columns"),
            Line2D([], [], color="0.35", ls="--", lw=0.8, label="model extent"),
            Line2D([], [], color="k", lw=1.2, label="fault (named)"),
            Patch(color="#9ecae1", label="water body"),
            Line2D([], [], color="k", marker="o", ls="", ms=3.2, label="wellhead"),
            Line2D([], [], color="0.55", lw=1.0, label="well track"),
            Line2D([], [], color="#d62728", ls=(0, (5, 3)), lw=1.2, label="section line")]


def section_setup(ctx):
    wells, ws, extent = ctx["wells"], ctx["ws"], ctx["extent"]
    hot = ws.sort_values(["Tmax_degC", "well"], ascending=[False, True]).well.tolist()
    wa, wb_ = hot[0], hot[1]
    A = wells[wa]["track"][0, :2].copy()
    B = wells[wb_]["track"][0, :2].copy()
    d = B - A
    u = d / np.hypot(*d)
    n = np.array([-u[1], u[0]])
    xmin, xmax, ymin, ymax = extent
    t0, t1 = -np.inf, np.inf
    for comp, lo, hi in ((0, xmin, xmax), (1, ymin, ymax)):
        ta = (lo - A[comp]) / u[comp]
        tb = (hi - A[comp]) / u[comp]
        t0 = max(t0, min(ta, tb))
        t1 = min(t1, max(ta, tb))
    sg = np.arange(np.ceil(t0 / 10.0) * 10.0, t1, 10.0)
    pts = A[None, :] + sg[:, None] * u[None, :]
    surf = np.full(len(sg), np.nan)
    for name, path in ctx["polys"]:
        m = path.contains_points(pts) & np.isnan(surf)
        if m.any():
            surf[m] = ctx["geo"]["surf"].get(name, ctx["geo"]["grid_top"])
    proj = {}
    offsets = {}
    for w in ctx["wnames"]:
        tr = wells[w]["track"]
        off = (tr[:, 0] - A[0]) * n[0] + (tr[:, 1] - A[1]) * n[1]
        s_tr = (tr[:, 0] - A[0]) * u[0] + (tr[:, 1] - A[1]) * u[1]
        mean_off = float(np.mean(np.abs(off)))
        offsets[w] = mean_off
        if mean_off > TREND_NEAR:
            continue
        g = ctx["tgt"][ctx["tgt"].well == w]
        xi = np.array([interp_track(tr, z)[0] for z in g.elevation_mRL])
        yi = np.array([interp_track(tr, z)[1] for z in g.elevation_mRL])
        si = (xi - A[0]) * u[0] + (yi - A[1]) * u[1]
        proj[w] = dict(s_med=float(np.median(si)), s_pts=si, z=g.elevation_mRL.to_numpy(float),
                       T=g.T_degC.to_numpy(float), P=g.P_bar.to_numpy(float),
                       bm=g.boiling_margin_degC.to_numpy(float), s_tr=s_tr, z_tr=tr[:, 2],
                       off=mean_off, shoe=wells[w]["shoe"], fz=wells[w]["fz"],
                       colour=ctx["well_colour"][w])
    # faults and water crossing the line
    cross = []
    for f in ctx["feats"]:
        for s in f["segs"]:
            for sv in seg_line_intersections(s, A, u):
                if t0 <= sv <= t1:
                    cross.append((sv, f["name"], f["kind"]))
    return dict(A=A, u=u, n=n, wa=wa, wb=wb_, s0=t0, s1=t1, sg=sg, surf=surf, proj=proj,
                offsets=offsets, cross=cross, B=B)


def prof_interp(z_desc, v_desc, zarr):
    zs = np.asarray(z_desc, float)[::-1]
    vs = np.asarray(v_desc, float)[::-1]
    out = np.interp(zarr, zs, vs)
    out[(zarr < zs[0]) | (zarr > zs[-1])] = np.nan
    return out


def envelope_field(sec, zg, kind, ctx):
    """Data envelope (4 section rules, 8.2 F4): linear between neighbouring wells where both have data."""
    sg = sec["sg"]
    F = np.full((len(zg), len(sg)), np.nan)
    order = sorted(sec["proj"], key=lambda w: sec["proj"][w]["s_med"])

    def value(w):
        p = sec["proj"][w]
        v = prof_interp(p["z"], p["T"] if kind == "T" else p["P"], zg)
        if kind == "P":
            v = v - trend_p(ctx, zg)
        return v

    for a, b in zip(order[:-1], order[1:]):
        sa, sb = sec["proj"][a]["s_med"], sec["proj"][b]["s_med"]
        if sb - sa < 1.0:
            continue
        cols = np.where((sg >= sa) & (sg <= sb))[0]
        if not len(cols):
            continue
        wgt = (sg[cols] - sa) / (sb - sa)
        va = value(a)
        vb = value(b)
        blk = (1 - wgt)[None, :] * va[:, None] + wgt[None, :] * vb[:, None]
        cur = F[:, cols]
        F[:, cols] = np.where(np.isnan(blk), cur, blk)
    return F, order


def draw_section_frame(ax, sec, ctx, kind):
    """Wells, faults, surface and labels common to F4 and its companions."""
    sg = sec["sg"]
    ax.plot(sg, sec["surf"], color="#2ca02c", lw=1.2, zorder=6)
    ax.text(sg[0] + 30, np.nanmax(sec["surf"]) + 8, "model surface", color="#2ca02c", fontsize=6, zorder=7)
    for w, p in sec["proj"].items():
        ax.plot(p["s_tr"], p["z_tr"], color="k", lw=1.0, zorder=7)
        ax.text(np.median(p["s_pts"]), ELEV_LIM[0] + 40 + 160 * (sorted(sec["proj"]).index(w) % 2), f"{w}\n{p['off']:.0f} m off-line",
                ha="center", va="bottom", fontsize=6, zorder=9, color=p["colour"],
                bbox=dict(fc="white", ec="none", alpha=0.7, pad=0.5))
        ax.hlines(p["shoe"], np.median(p["s_pts"]) - 35, np.median(p["s_pts"]) + 35, color="k", lw=1.6, zorder=8)
        ax.text(np.median(p["s_pts"]) + 38, p["shoe"], "shoe", fontsize=5, va="center", zorder=8)
        for f in p["fz"]:
            ax.vlines(np.median(p["s_pts"]), min(f["top"], f["bot"]), max(f["top"], f["bot"]),
                      color="#d62728", lw=4.0, zorder=8)
            ax.text(np.median(p["s_pts"]) - 40, 0.5 * (f["top"] + f["bot"]), f["q"][:3].lower(),
                    fontsize=5, ha="right", va="center", color="#d62728", zorder=9)
    for sv, name, kind_f in sec["cross"]:
        ax.plot([sv, sv], [ELEV_LIM[1] - 25, ELEV_LIM[1] - 5], color="k" if kind_f == "fault" else "#2171b5",
                lw=1.5, zorder=9)
        ax.text(sv, ELEV_LIM[1] - 30, name, rotation=90, fontsize=6, ha="center", va="top", zorder=9)
    ax.set_xlim(sec["s0"], sec["s1"])
    ax.set_ylim(*ELEV_LIM)
    ax.set_xlabel("Distance along section A-A' from WELL-A (m)", fontsize=8)
    ax.set_ylabel("Elevation (mRL)", fontsize=8)
    ax.tick_params(labelsize=7)
    ax.text(sec["s0"] + 60, ELEV_LIM[0] + 120, "A", fontsize=12, weight="bold", color="#d62728")
    ax.text(sec["s1"] - 150, ELEV_LIM[0] + 120, "A'", fontsize=12, weight="bold", color="#d62728")


def stipple_low_margin(ax, sec):
    for w, p in sec["proj"].items():
        z, bm = p["z"], p["bm"]
        for i in range(len(z) - 1):
            if np.isfinite(bm[i]) and np.isfinite(bm[i + 1]) and bm[i] < NEAR_BOIL and bm[i + 1] < NEAR_BOIL:
                zz = np.linspace(z[i + 1], z[i], 12)
                ax.plot(np.full_like(zz, p["s_med"]) + 6, zz, "k.", ms=2.2, zorder=8)


def fig_section(ctx, sec, kind, path, interp=False):
    if kind == "T":
        zg = np.arange(ELEV_LIM[1], ELEV_LIM[0] - 5, -5.0)
        F, order = envelope_field(sec, zg, "T", ctx)
        cmap, norm, unit = T_CMAP, T_NORM, "degC"
        iso = ISO_LEVELS
    else:
        zg = np.arange(ELEV_LIM[1], ELEV_LIM[0] - 5, -5.0)
        F, order = envelope_field(sec, zg, "P", ctx)
        cmap, norm, unit = P_CMAP, P_NORM, "bar (deviation from trend)"
        iso = [-6.0, -4.0, -2.0, 2.0, 4.0, 6.0]
    sg = sec["sg"]
    Sg, Zg = np.meshgrid(sg, zg)
    surf_grid = np.interp(sg, sg, np.where(np.isnan(sec["surf"]), ELEV_LIM[1], sec["surf"]))[None, :]
    above = Zg > surf_grid
    F = np.where(above, np.nan, F)
    nodata = np.isnan(F) & ~above
    if interp:
        F, nodata = interpolate_section(sec, sg, zg, kind, ctx, above)
    fig, ax = plt.subplots(figsize=(12, 7.4))
    fig.subplots_adjust(left=0.07, right=0.92, top=0.90, bottom=0.22)
    ax.set_facecolor("white")
    mesh = ax.pcolormesh(Sg, Zg, np.ma.masked_invalid(F), cmap=cmap, norm=norm, shading="auto", zorder=2)
    ax.pcolormesh(Sg, Zg, np.ma.masked_where(~nodata, np.zeros_like(Zg)), cmap=ListedColormap(["#ececec"]),
                  shading="auto", zorder=1, vmin=0, vmax=1)
    plt.rcParams["hatch.color"] = "#9a9a9a"
    ax.contourf(Sg, Zg, nodata.astype(float), levels=[0.5, 1.5], colors="none", hatches=["////"], zorder=3)
    plt.rcParams["hatch.color"] = "k"
    if kind == "T":
        cs = ax.contour(Sg, Zg, F, levels=iso, colors="k", linewidths=[ISO_LW[v] for v in iso], zorder=4)
        ax.clabel(cs, fmt="%d degC", fontsize=6, inline=True)
    else:
        cs = ax.contour(Sg, Zg, F, levels=iso, colors="0.25", linewidths=0.6, zorder=4)
        ax.clabel(cs, fmt="%+d bar", fontsize=6, inline=True)
    valid = np.isfinite(F).astype(float)
    if not interp:
        ax.contour(Sg, Zg, valid, levels=[0.5], colors="k", linestyles=":", linewidths=0.7, zorder=5)
    draw_section_frame(ax, sec, ctx, kind)
    for w, p in sec["proj"].items():
        if kind == "T":
            vals = p["T"]
            ax.scatter(np.full(len(p["z"]), p["s_med"]), p["z"], c=vals, cmap=cmap, norm=norm,
                       s=16, edgecolor="k", lw=0.3, zorder=9)
        else:
            vals = p["P"] - trend_p(ctx, p["z"])
            ax.scatter(np.full(len(p["z"]), p["s_med"]), p["z"], c=vals, cmap=cmap, norm=norm,
                       s=16, edgecolor="k", lw=0.3, zorder=9)
    if not interp:
        stipple_low_margin(ax, sec)
    cb = fig.colorbar(mesh, ax=ax, pad=0.01, fraction=0.025)
    cb.set_label("Temperature (degC)" if kind == "T" else "Pressure deviation from trend (bar)", fontsize=8)
    cb.ax.tick_params(labelsize=7)
    ve = ((ELEV_LIM[1] - ELEV_LIM[0]) / (fig.get_size_inches()[1] * 0.80)) / ((sec["s1"] - sec["s0"]) / (fig.get_size_inches()[0] * 0.73))
    if kind == "T":
        title = f"Section A-A' temperature (degC), elevation mRL; vertical exaggeration x{ve:.1f}"
    else:
        title = f"Section A-A' pressure deviation from the verified trend (bar), elevation mRL; vertical exaggeration x{ve:.1f}"
    if interp:
        title = "[interpolated, not data] " + title
    ax.set_title(title, fontsize=9)
    handles = [plt.Rectangle((0, 0), 1, 1, fc="#ececec", ec="0.5", hatch="////", label="no data (outside the envelope)"),
               plt.Line2D([], [], color="k", ls=":", label="data envelope (between neighbouring wells)"),
               plt.Line2D([], [], color="k", marker="o", ls="", ms=4, label="survey point (colour = value)"),
               plt.Line2D([], [], color="k", lw=1.0, label="model surface"),
               plt.Line2D([], [], color="#d62728", lw=4.0, label="feed zone (label = quality)"),
               plt.Line2D([], [], color="k", lw=1.6, label="casing shoe"),
               plt.Line2D([], [], color="k", ls="", marker=".", ms=3, label="boiling margin < 5 degC (stipple)")]
    if interp:
        handles.insert(1, plt.Line2D([], [], color="k", ls="--", label="F4 data envelope (dashed)"))
        handles.append(plt.Line2D([], [], color="0.4", ls="", label="colours only inside convex hull of data;\nnot constrained elsewhere"))
    if kind == "T":
        handles.append(plt.Line2D([], [], color="k", lw=1.8, label="100 and 200 degC isotherms"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=7, frameon=False)
    fig.savefig(path, dpi=FIG_DPI)
    plt.close(fig)
    return ve


def interpolate_section(sec, sg, zg, kind, ctx, above):
    pts, vals = [], []
    for w, p in sec["proj"].items():
        v = p["T"] if kind == "T" else (p["P"] - trend_p(ctx, p["z"]))
        pts += [(s, z) for s, z in zip(np.full(len(p["z"]), p["s_med"]), p["z"])]
        vals += list(v)
    pts = np.array(pts)
    vals = np.array(vals)
    interp = LinearNDInterpolator(pts, vals)
    Sg, Zg = np.meshgrid(sg, zg)
    F = interp(Sg, Zg)
    F = np.where(above, np.nan, F)
    nodata = np.isnan(F) & ~above
    return F, nodata


def pt_figure(ctx, path):
    fig, ax = plt.subplots(figsize=(7.5, 5.8))
    tp = np.linspace(1.0, P_CRIT, 300)
    tt = np.array([tsat_c(p) for p in tp])
    ax.plot(tt, tp, color="k", lw=1.2, label="saturation curve (IF97 Region 4)")
    ax.plot([T_CRIT], [P_CRIT], "k*", ms=9, label="critical point (IF97)")
    for w in ctx["wnames"]:
        g = ctx["tgt"][ctx["tgt"].well == w]
        ax.scatter(g.T_degC, g.P_bar, s=14, color=ctx["well_colour"][w], edgecolor="k", lw=0.3, label=w, zorder=3)
    fz = ctx["fz"]
    ax.scatter(fz.T_degC, fz.P_bar, s=70, facecolor="none", edgecolor="k", lw=1.2, zorder=4,
               label="feed zone (mid elevation)")
    ax.set_xlabel("Temperature (degC)")
    ax.set_ylabel("Pressure (bar)")
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 300)
    ax.set_title("P-T state of survey points against the saturation curve (natural state)", fontsize=9)
    ax.grid(alpha=0.25)
    ax.legend(fontsize=6, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    fig.subplots_adjust(right=0.74, left=0.1, bottom=0.12, top=0.92)
    fig.savefig(path, dpi=FIG_DPI)
    plt.close(fig)


def fig_profile_axes(axT, axP, ctx, w, shared=True):
    W = ctx["wells"][w]
    g = ctx["tgt"][ctx["tgt"].well == w]
    colour = ctx["well_colour"][w]
    z = g.elevation_mRL.to_numpy(float)
    T = g.T_degC.to_numpy(float)
    P = g.P_bar.to_numpy(float)
    bm = g.boiling_margin_degC.to_numpy(float)
    zz = np.arange(z.min(), z.max() + 0.01, 2.0)
    pz = np.interp(zz, z[::-1], P[::-1])
    ts = np.array([tsat_c(p) for p in pz])
    ok = np.isfinite(ts) & (pz <= P_CRIT)
    axT.plot(ts[ok], zz[ok], color="k", lw=0.8, ls="--", label="Tsat(P)")
    for i in range(len(z) - 1):
        if np.isfinite(bm[i]) and np.isfinite(bm[i + 1]) and bm[i] < NEAR_BOIL and bm[i + 1] < NEAR_BOIL:
            axT.axhspan(z[i + 1], z[i], color="orange", alpha=0.35, zorder=1)
    axT.scatter(T, z, s=16, color=colour, edgecolor="k", lw=0.3, zorder=4, label="survey T")
    for f in W["fz"]:
        top, bot = max(f["top"], f["bot"]), min(f["top"], f["bot"])
        axT.axhspan(bot, top, color="#d62728", alpha=0.18, zorder=2)
        axT.text(2, 0.5 * (top + bot), f"feed zone ({f['q']})", fontsize=6, va="center", color="#a50f15")
    axT.axhline(W["shoe"], color="0.25", ls="-.", lw=0.9)
    axT.text(298, W["shoe"], "shoe", fontsize=6, ha="right", va="bottom")
    surf = well_surface(ctx, w)
    axT.axhline(surf, color="#2ca02c", ls=":", lw=1.0)
    axT.text(298, surf, "model surface", fontsize=6, ha="right", va="bottom", color="#2ca02c")
    axT.axhline(W["wellhead_z"], color="0.4", ls=":", lw=0.8)
    axT.text(298, W["wellhead_z"], "wellhead", fontsize=6, ha="right", va="bottom", color="0.3")
    ws = ctx["ws"].set_index("well").loc[w]
    for k, v in (("z_100C_mRL", 100), ("z_200C_mRL", 200)):
        if np.isfinite(ws[k]):
            axT.axhline(ws[k], color="#6a51a3", lw=1.2 if v == 200 else 0.9, ls="-")
            axT.text(2, ws[k], f"{v} degC at {ws[k]:.0f} mRL", fontsize=6, va="bottom", color="#6a51a3")
    axT.set_xlim(0, 300)
    axT.set_ylim(*ELEV_LIM)
    axT.set_xlabel("Temperature (degC)", fontsize=8)
    axT.set_ylabel("Elevation (mRL)", fontsize=8)
    axT.tick_params(labelsize=7)
    axT.grid(alpha=0.2)
    if axP is not None:
        axP.scatter(P, z, s=16, color=colour, edgecolor="k", lw=0.3, zorder=4, label="survey P")
        fzmid = [0.5 * (f["top"] + f["bot"]) for f in W["fz"]]
        fzP = [interp_profile(W["pz"], W["pp"], m) for m in fzmid]
        axP.scatter(fzP, fzmid, s=70, color=colour, edgecolor="k", lw=1.2, zorder=5, label="feed-zone P")
        pv = ws.pivot_elev_mRL
        if np.isfinite(pv):
            axP.axhline(pv, color="0.2", ls="--", lw=0.8)
            axP.text(2, pv, "pivot", fontsize=6, va="bottom")
        wl = ws.water_level_mRL
        if np.isfinite(wl) and ELEV_LIM[0] <= wl <= ELEV_LIM[1]:
            axP.axhline(wl, color="#2171b5", ls=":", lw=1.0)
            axP.text(2, wl, "water level", fontsize=6, va="bottom", color="#2171b5")
        pmax = max(P.max(), 50.0)
        axP.set_xlim(0, np.ceil(pmax / 50.0) * 50.0)
        axP.set_xlabel("Pressure (bar)", fontsize=8)
        axP.tick_params(labelsize=7, labelleft=False)
        axP.grid(alpha=0.2)
        axP.set_ylim(*ELEV_LIM)


def fig_profiles(ctx):
    for w in ctx["wnames"]:
        fig = plt.figure(figsize=(9.6, 8.2))
        axT = fig.add_axes([0.085, 0.09, 0.30, 0.84])
        axP = fig.add_axes([0.42, 0.09, 0.21, 0.84], sharey=axT)
        axL = fig.add_axes([0.72, 0.60, 0.25, 0.30])
        axI = fig.add_axes([0.72, 0.08, 0.26, 0.40])
        fig_profile_axes(axT, axP, ctx, w)
        axT.set_title("Temperature", fontsize=9)
        axP.set_title("Pressure", fontsize=9)
        plot_base(axL, ctx, highlight=w, labels=False, zoom=None)
        axL.set_title("Locator", fontsize=8)
        axL.tick_params(labelsize=5)
        axL.set_xlabel("")
        axL.set_ylabel("")
        xs = [ctx["wells"][v]["track"][:, 0] for v in ctx["wnames"]]
        ys = [ctx["wells"][v]["track"][:, 1] for v in ctx["wnames"]]
        xc, yc = ctx["wells"][w]["track"][0, 0], ctx["wells"][w]["track"][0, 1]
        span = 1500.0
        axL.set_xlim(xc - span, xc + span)
        axL.set_ylim(yc - span, yc + span)
        axI.axis("off")
        r = ctx["ws"].set_index("well").loc[w]
        p = ctx["pc"].set_index("well").loc[w]
        g = ctx["ns"][w]
        lines = [
            f"{w}",
            f"survey date: {ctx['wells'][w]['survey_date']}",
            f"natural-state verdict: {r.ns_verdict}",
            f"confidence: {r.ns_confidence}",
            f"regime: {r.regime}",
            f"test: {g['rule'].split(':')[0]}",
            "rule: " + textwrap.fill(g['rule'].split(': ', 1)[1], 42),
            f"Tmax {r.Tmax_degC:.1f} degC at {r.Tmax_elev_mRL:.0f} mRL",
            f"gauge: {p.gauge_test if isinstance(p.gauge_test, str) else 'none'}",
        ]
        axI.text(0.0, 1.0, "\n".join(lines), fontsize=7.5, va="top", ha="left", family="DejaVu Sans")
        fig.suptitle(f"{w}: temperature and pressure profiles (survey {ctx['wells'][w]['survey_date']}), "
                     f"elevation mRL; Tsat line only below {P_CRIT:.2f} bar", fontsize=9)
        fig.savefig(os.path.join(CALF, f"profile_{w}.png"), dpi=FIG_DPI)
        plt.close(fig)
    # group panels
    for grp, members in ctx["groups"].items():
        members = [m for m in members if m in ctx["wells"]]
        n = len(members)
        ncol = min(5, n)
        nrow = int(math.ceil(n / ncol))
        fig = plt.figure(figsize=(3.0 * ncol + 2.8, 4.1 * nrow + 0.6))
        axs = []
        for i, w in enumerate(members):
            r_, c_ = divmod(i, ncol)
            ax = fig.add_axes([0.07 + c_ * (0.82 / ncol), 0.08 + (nrow - 1 - r_) * (0.86 / nrow),
                               0.82 / ncol * 0.80, 0.86 / nrow * 0.84])
            axs.append(ax)
            fig_profile_axes(ax, None, ctx, w)
            ax.set_title(w, fontsize=8, color=ctx["well_colour"][w])
            if c_ > 0:
                ax.tick_params(labelleft=False)
                ax.set_ylabel("")
        axL = fig.add_axes([0.875, 0.60, 0.12, 0.30])
        plot_base(axL, ctx, highlight=None, labels=False)
        for w in members:
            tr = ctx["wells"][w]["track"]
            axL.plot(tr[:, 0], tr[:, 1], color=ctx["well_colour"][w], lw=2.0)
        axL.set_title(f"Locator: {grp}", fontsize=7)
        axL.tick_params(labelsize=4)
        axL.set_xlabel("")
        axL.set_ylabel("")
        fig.suptitle(f"Temperature profiles, group '{grp}': same axes in every panel "
                     f"(0-300 degC; elevation mRL; orange = boiling margin < 5 degC)", fontsize=9)
        fig.savefig(os.path.join(CALF, f"profile_panels_{grp}.png"), dpi=FIG_DPI)
        plt.close(fig)


def band_symbol(ax, x, y, d, kind):
    if kind == "T":
        if abs(d) <= 10:
            ax.plot(x, y, marker="+", ms=9, mew=1.6, color="0.2", zorder=5)
            return
        size = 40 if abs(d) <= 20 else 80 if abs(d) <= 30 else 130 if abs(d) <= 40 else 200
    else:
        if abs(d) <= 2:
            ax.plot(x, y, marker="+", ms=9, mew=1.6, color="0.2", zorder=5)
            return
        size = 40 if abs(d) <= 4 else 80 if abs(d) <= 6 else 130 if abs(d) <= 8 else 200
    colr = "#d62728" if d > 0 else "#2171b5"
    ax.scatter([x], [y], s=size, color=colr, edgecolor="k", lw=0.4, zorder=5)


def fig_bands(ctx, kind):
    bl = ctx["band_list"]
    bands = ctx["bands"]
    n = len(bl)
    ncol = min(4, n)
    nrow = int(math.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.9 * ncol + 2.6, 3.9 * nrow + 0.8), squeeze=False)
    col_name = "T_minus_field_median_degC" if kind == "T" else "P_minus_trend_bar"
    for k, ax in enumerate(axes.ravel()):
        if k >= n:
            ax.axis("off")
            continue
        top, bot, med, nw = bl[k]
        mid = 0.5 * (top + bot)
        plot_base(ax, ctx, labels=False, zoom=(2817000, 2824500, 6331000, 6338500))
        sub = bands[(bands.band_top_mRL == top) & (bands.band_bottom_mRL == bot)]
        have = set(sub.well)
        for w in ctx["wnames"]:
            tr = ctx["wells"][w]["track"]
            zz = tr[:, 2]
            x = float(np.interp(mid, zz[::-1], tr[::-1, 0]))
            y = float(np.interp(mid, zz[::-1], tr[::-1, 1]))
            if w in have:
                d = float(sub[sub.well == w][col_name].iloc[0])
                band_symbol(ax, x, y, d, kind)
            else:
                ax.plot(x, y, "o", mfc="none", mec="0.6", ms=5, zorder=5)
            ax.text(x + 40, y + 40, w.replace("WELL-", "W"), fontsize=6, zorder=6)
        if kind == "T":
            ax.set_title(f"{top:g} to {bot:g} mRL\nfield median {med:.0f} degC ({nw} wells)", fontsize=8)
        else:
            ax.set_title(f"{top:g} to {bot:g} mRL\ndeviation from trend (field wells: {nw})", fontsize=8)
        ax.tick_params(labelsize=5)
        ax.set_xlabel("")
        ax.set_ylabel("")
    handles = [plt.Line2D([], [], marker="+", ls="", color="0.2", ms=9, label="within +/-10 degC" if kind == "T" else "within +/-2 bar"),
               plt.Line2D([], [], marker="o", ls="", mfc="#d62728", mec="k", ms=6, label="hotter" if kind == "T" else "above trend"),
               plt.Line2D([], [], marker="o", ls="", mfc="#2171b5", mec="k", ms=6, label="colder" if kind == "T" else "below trend"),
               plt.Line2D([], [], marker="o", ls="", mfc="none", mec="0.6", ms=5, label="no data in band")]
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=7, frameon=False)
    if kind == "T":
        fig.suptitle("Temperature minus field median by elevation band (degC): symbol size = class "
                     "(10-20, 20-30, 30-40, > 40 degC)", fontsize=9)
        path = os.path.join(CALF, "maps_temperature_bands.png")
    else:
        fig.suptitle("Pressure minus gauge-verified hydrostatic trend by elevation band (bar): symbol size = class "
                     "(2-4, 4-6, 6-8, > 8 bar)", fontsize=9)
        path = os.path.join(CALF, "maps_pressure_bands.png")
    fig.subplots_adjust(right=0.86, left=0.05, top=0.90, bottom=0.05, wspace=0.25, hspace=0.35)
    fig.savefig(path, dpi=FIG_DPI)
    plt.close(fig)


def interp_band_maps(ctx, kind):
    bl = ctx["band_list"]
    bands = ctx["bands"]
    n = len(bl)
    ncol = min(4, n)
    nrow = int(math.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.9 * ncol + 1.4, 3.9 * nrow + 0.9), squeeze=False)
    col_name = "T_minus_field_median_degC" if kind == "T" else "P_minus_trend_bar"
    if kind == "T":
        cmap, norm = plt.get_cmap("RdBu_r", 30), plt.Normalize(-40, 40)
    else:
        cmap, norm = P_CMAP, plt.Normalize(-10, 10)
    mesh = None
    for k, ax in enumerate(axes.ravel()):
        if k >= n:
            ax.axis("off")
            continue
        top, bot, med, nw = bl[k]
        mid = 0.5 * (top + bot)
        zoom = (2817000, 2824500, 6331000, 6338500)
        sub = bands[(bands.band_top_mRL == top) & (bands.band_bottom_mRL == bot)]
        pts, vals = [], []
        for w in sub.well:
            tr = ctx["wells"][w]["track"]
            zz = tr[:, 2]
            x = float(np.interp(mid, zz[::-1], tr[::-1, 0]))
            y = float(np.interp(mid, zz[::-1], tr[::-1, 1]))
            pts.append((x, y))
            vals.append(float(sub[sub.well == w][col_name].iloc[0]))
        xs = np.linspace(zoom[0], zoom[1], 160)
        ys = np.linspace(zoom[2], zoom[3], 160)
        Xg, Yg = np.meshgrid(xs, ys)
        if len(pts) >= 3:
            interpo = LinearNDInterpolator(np.array(pts), np.array(vals))
            Z = interpo(Xg, Yg)
        else:
            Z = np.full(Xg.shape, np.nan)
        cm_ = ax.pcolormesh(Xg, Yg, np.ma.masked_invalid(Z), cmap=cmap, norm=norm, shading="auto", zorder=2)
        mesh = cm_
        nd = np.isnan(Z).astype(float)
        ax.contourf(Xg, Yg, nd, levels=[0.5, 1.5], colors="#ececec", zorder=2.5)
        plot_base(ax, ctx, labels=False, zoom=zoom)
        for (x, y), v in zip(pts, vals):
            ax.plot(x, y, "k.", ms=4, zorder=6)
        if len(pts) >= 3:
            hull_pts = np.array(pts)
            from scipy.spatial import ConvexHull
            h = ConvexHull(hull_pts)
            hv = np.append(h.vertices, h.vertices[0])
            ax.plot(hull_pts[hv, 0], hull_pts[hv, 1], ls="--", color="k", lw=0.8, zorder=6)
        ax.set_title(f"{top:g} to {bot:g} mRL: {nw} wells", fontsize=8)
        ax.tick_params(labelsize=5)
        ax.set_xlabel("")
        ax.set_ylabel("")
    if mesh is not None:
        cb = fig.colorbar(mesh, ax=axes.ravel().tolist(), fraction=0.025, pad=0.02)
        cb.set_label("degC minus field median" if kind == "T" else "bar, deviation from trend", fontsize=8)
    title_kind = "temperature minus field median (degC)" if kind == "T" else "pressure minus trend (bar)"
    fig.suptitle(f"[interpolated, not data] {title_kind} by elevation band. Linear interpolation on a Delaunay triangulation "
                 f"of well values; grey = outside the convex hull of wells (not constrained)", fontsize=8)
    fig.subplots_adjust(left=0.05, right=0.88, top=0.88, bottom=0.05, wspace=0.25, hspace=0.35)
    fname = "maps_temperature_interpolated.png" if kind == "T" else "maps_pressure_interpolated.png"
    fig.savefig(os.path.join(ADDF, fname), dpi=FIG_DPI)
    plt.close(fig)


def fig_conceptual(ctx, sec, path):
    fig = plt.figure(figsize=(13.5, 6.8))
    axP = fig.add_axes([0.05, 0.22, 0.40, 0.72])
    axS = fig.add_axes([0.52, 0.22, 0.40, 0.72])
    plot_base(axP, ctx, labels=True, zoom=None)
    axP.set_title("1  Plan view (base layer; arrows as keyed in the report table)", fontsize=9)
    wnames = ctx["wnames"]
    ws = ctx["ws"].set_index("well")
    upflow = [w for w in wnames if ws.loc[w, "regime"] == "upflow"]
    outflow = [w for w in wnames if ws.loc[w, "regime"] == "outflow"]
    cold = [w for w in wnames if ws.loc[w, "regime"] == "recharge/downflow"]
    wa = sec["wa"]
    A = ctx["wells"][wa]["track"][0, :2]
    el = 1
    key = []
    for w in upflow:
        P0 = ctx["wells"][w]["track"][0, :2]
        axP.annotate("", xy=(P0[0] + 300, P0[1] + 300), xytext=(P0[0], P0[1]),
                     arrowprops=dict(arrowstyle="-|>", color="#d62728", lw=2.0, ls=(0, (4, 2))), zorder=10)
        axP.text(P0[0] + 320, P0[1] + 330, f"{el}", color="#d62728", fontsize=9, weight="bold", zorder=11)
        key.append((el, "Deep hot upflow (upward arrow)", w,
                    f"Tmax {ws.loc[w, 'Tmax_degC']:.1f} degC at {ws.loc[w, 'Tmax_elev_mRL']:.0f} mRL, no inversion; minimum boiling margin {ws.loc[w, 'min_boiling_margin_degC']:.1f} degC", "interp", "dashed (inferred from the highest temperature)"))
        el += 1
    for w in outflow:
        P0 = ctx["wells"][w]["track"][0, :2]
        if upflow:
            P_up = ctx["wells"][upflow[0]]["track"][0, :2]
            axP.annotate("", xy=P0, xytext=P_up,
                         arrowprops=dict(arrowstyle="-", color="#fd8d3c", lw=1.0, ls=(0, (4, 3))), zorder=6)
        d = P0 - (ctx["wells"][upflow[0]]["track"][0, :2] if upflow else P0)
        u = d / (np.hypot(*d) + 1e-9)
        axP.annotate("", xy=P0 + 260 * u, xytext=P0 + 20 * u,
                     arrowprops=dict(arrowstyle="-|>", color="#fd8d3c", lw=2.2), zorder=10)
        iv = ws.loc[w, "inversion_degC"]
        key.append((el, "Lateral outflow at inversion depth (solid arrow in plan)", w,
                    f"inversion {iv:.1f} degC between {ws.loc[w, 'inversion_top_mRL']:.0f} and {ws.loc[w, 'inversion_base_mRL']:.0f} mRL", "obs", "solid (measured inversion); path to the upflow dashed"))
        axP.text(P0[0] + 300 * u[0] + 40, P0[1] + 300 * u[1] + 40, f"{el}", color="#d94801", fontsize=9, weight="bold")
        el += 1
    for w in cold:
        P0 = ctx["wells"][w]["track"][0, :2]
        axP.plot(P0[0], P0[1], "o", mfc="none", mec="#2171b5", ms=14, mew=2.0, zorder=9)
        axP.annotate("", xy=P0 + np.array([0, -250]), xytext=P0 + np.array([0, 0]),
                     arrowprops=dict(arrowstyle="-|>", color="#2171b5", lw=2.0, ls=(0, (3, 2))), zorder=10)
        axP.text(P0[0] + 40, P0[1] - 290, f"{el}", color="#2171b5", fontsize=9, weight="bold")
        key.append((el, "Cold recharge or downflow (blue arrow)", w,
                    f"Tmax {ws.loc[w, 'Tmax_degC']:.1f} degC; regime recharge/downflow", "obs", "dashed (direction inferred from cold profile)"))
        el += 1
    # pressure connections between gauge-verified wells
    ver = ctx["trend"]["verified"]
    if len(ver) >= 2:
        a, b = ver[0], ver[1]
        za = ctx["wells"][a]
        zb = ctx["wells"][b]
        common = sorted(set(np.round(za["tz"], 1)) & set(np.round(zb["tz"], 1)), reverse=True)
        pick = None
        for z in common:
            if ELEV_LIM[0] < z < 0:
                pick = z
                break
        if pick is not None:
            dP = interp_profile(za["pz"], za["pp"], pick) - interp_profile(zb["pz"], zb["pp"], pick)
            pa, pb = za["track"][0, :2], zb["track"][0, :2]
            axP.plot([pa[0], pb[0]], [pa[1], pb[1]], ls=(0, (2, 2)), color="0.2", lw=1.0, zorder=7)
            key.append((el, "Pressure connection between wells (thin dashed line)", f"{a}-{b}",
                        f"pressure difference {dP:.2f} bar at {pick:.0f} mRL (common elevation)", "obs",
                        "dashed (connection inferred from pressures)"))
            mid = 0.5 * (pa + pb)
            axP.text(mid[0], mid[1], f"{el}  {dP:+.1f} bar", fontsize=8, color="0.1", zorder=11,
                     bbox=dict(fc="white", ec="none", alpha=0.7))
            el += 1
    # near-boiling stipple at upflow
    for w in upflow:
        P0 = ctx["wells"][w]["track"][0, :2]
        axP.plot(P0[0] - 60, P0[1] - 60, "k.", ms=5, zorder=11)
    key.append((el, "Two-phase or near-boiling zone (stipple)", ",".join(upflow) or "none",
                f"minimum boiling margin {ws.loc[upflow[0], 'min_boiling_margin_degC']:.1f} degC"
                if upflow else "no upflow well", "obs", "stipple (measured margin)"))
    el += 1
    key.append((el, "Cap (low-permeability, conductive gradient)", "none",
                "not resolved by the survey data; no cap band drawn", "interp", "absent; see report"))
    axP.set_xlabel("Easting (m, NZTM)", fontsize=8)
    # section panel
    fz_back = envelope_field(sec, np.arange(ELEV_LIM[1], ELEV_LIM[0] - 5, -5.0), "T", ctx)[0]
    zg = np.arange(ELEV_LIM[1], ELEV_LIM[0] - 5, -5.0)
    sg = sec["sg"]
    Sg, Zg = np.meshgrid(sg, zg)
    surf_grid = np.interp(sg, sg, np.where(np.isnan(sec["surf"]), ELEV_LIM[1], sec["surf"]))[None, :]
    F = np.where(Zg > surf_grid, np.nan, fz_back)
    nodata = np.isnan(F) & ~(Zg > surf_grid)
    axS.pcolormesh(Sg, Zg, np.ma.masked_invalid(F), cmap=T_CMAP, norm=T_NORM, shading="auto")
    axS.pcolormesh(Sg, Zg, np.ma.masked_where(~nodata, np.zeros_like(Zg)), cmap=ListedColormap(["#ececec"]),
                   shading="auto", vmin=0, vmax=1)
    cs = axS.contour(Sg, Zg, F, levels=ISO_LEVELS, colors="k", linewidths=[ISO_LW[v] for v in ISO_LEVELS])
    axS.clabel(cs, fmt="%d degC", fontsize=6, inline=True)
    draw_section_frame(axS, sec, ctx, "T")
    # arrows in section for on-line wells
    for w, p in sec["proj"].items():
        s_m = p["s_med"]
        if w in upflow:
            zt = p["z"][np.argmax(p["T"])]
            axS.annotate("", xy=(s_m, zt + 380), xytext=(s_m, zt - 20),
                         arrowprops=dict(arrowstyle="-|>", color="#d62728", lw=2.2, ls=(0, (4, 2))), zorder=12)
            axS.text(s_m + 40, zt + 300, "1", color="#d62728", fontsize=10, weight="bold", zorder=13)
            ztop = np.nanmax(p["bm"][p["bm"] < NEAR_BOIL]) if np.any(p["bm"] < NEAR_BOIL) else None
        if w in outflow:
            iv_top = ws.loc[w, "inversion_top_mRL"]
            iv_base = ws.loc[w, "inversion_base_mRL"]
            if np.isfinite(iv_top):
                zmid = 0.5 * (iv_top + iv_base)
                axS.annotate("", xy=(s_m + 260, zmid), xytext=(s_m + 20, zmid),
                             arrowprops=dict(arrowstyle="-|>", color="#fd8d3c", lw=2.4), zorder=12)
                num = [k[0] for k in key if k[2] == w and "Lateral" in k[1]]
                axS.text(s_m + 280, zmid, f"{num[0] if num else ''}", color="#d94801", fontsize=10, weight="bold", zorder=13)
        if w in cold:
            zc = p["z"].min() + 150
            axS.annotate("", xy=(s_m, zc - 300), xytext=(s_m, zc),
                         arrowprops=dict(arrowstyle="-|>", color="#2171b5", lw=2.0, ls=(0, (3, 2))), zorder=12)
    axS.set_title("2  Section A-A' (same line and colour scale as F4); dashed = inferred, solid = measured", fontsize=9)
    fig.subplots_adjust(right=0.97)
    fig.suptitle("Conceptual model built from survey data (numbers key to the report table)", fontsize=10)
    handles = [
        plt.Line2D([], [], color="#d62728", lw=2.2, ls=(0, (4, 2)), label="deep hot upflow (dashed = inferred)"),
        plt.Line2D([], [], color="#fd8d3c", lw=2.4, label="lateral outflow (solid = measured inversion)"),
        plt.Line2D([], [], color="#2171b5", lw=2.0, ls=(0, (3, 2)), label="cold recharge or downflow (dashed)"),
        plt.Rectangle((0, 0), 1, 1, fc="#ececec", ec="0.5", label="cap: absent (not resolved)"),
        plt.Line2D([], [], color="k", ls="", marker=".", ms=5, label="two-phase / near-boiling (stipple)"),
        plt.Line2D([], [], color="0.2", ls=(0, (2, 2)), lw=1.0, label="pressure connection (labelled, bar)"),
        plt.Line2D([], [], color="k", ls="", marker="o", mfc="none", ms=8, label="cold region marker"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=7, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.savefig(path, dpi=FIG_DPI)
    plt.close(fig)
    return key


def fig_gauge(ctx, path):
    pcs = ctx["pc"].set_index("well")
    gw = [w for w in ctx["wnames"] if ctx["ns"][w]["gauge"]]
    if not gw:
        return
    fig, axs = plt.subplots(1, len(gw), figsize=(4.2 * len(gw) + 0.5, 4.6), squeeze=False)
    field_first = min(ctx["wells"][w]["first_flow"] for w in ctx["wnames"]
                      if ctx["wells"][w]["kind"] == "producer" and ctx["wells"][w]["first_flow"])
    for ax, w in zip(axs[0], gw):
        g = ctx["ns"][w]["gauge"]
        ser = g["series"]
        t = [dtm.date.fromisoformat(d) for d, _ in ser]
        v = [val for _, val in ser]
        ax.plot(t, v, "o-", color=ctx["well_colour"][w], ms=3, lw=1.0, label="transient gauge")
        ax.axhspan(g["first_P"] - GAUGE_PREC, g["first_P"] + GAUGE_PREC, color="#fdae6b", alpha=0.6, label="first reading +/- precision")
        ax.axhline(g["survey_P"], color="k", ls="--", lw=1.0, label=f"survey P at gauge ({g['survey_P']:.1f} bar, {ctx['wells'][w]['survey_date']})")
        own = g["own_first_flow"]
        if own:
            ax.axvline(dtm.date.fromisoformat(own), color="#2171b5", ls="-", lw=1.0, label="own first flow")
        ax.axvline(dtm.date.fromisoformat(field_first), color="0.3", ls=":", lw=1.2, label="first production in field")
        flows = ctx["wells"][w]["flows"].get("Production_Mass_Flow") or ctx["wells"][w]["flows"].get("Injection_Mass_Flow") or []
        for i, (d, q) in enumerate(flows):
            if q > 0:
                t0 = dtm.date.fromisoformat(d)
                t1 = dtm.date.fromisoformat(flows[i + 1][0]) if i + 1 < len(flows) else t0 + dtm.timedelta(days=30)
                ax.axvspan(t0, t1, color="0.85", alpha=0.6, lw=0, zorder=0)
        ax.set_title(f"{w}: {pcs.loc[w, 'verdict']} ({pcs.loc[w, 'confidence']})\n"
                     f"gauge at {g['gauge_elev']:.0f} mRL; diff {g['diff']:+.2f} bar", fontsize=8)
        ax.set_ylabel("Pressure (bar)", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.2)
        ax.legend(fontsize=5.5, loc="upper right", frameon=False)
        import matplotlib.dates as mdates
        ax.xaxis.set_major_locator(mdates.YearLocator(2))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
    fig.suptitle("Gauge check (grey shading = flowing). Reading precision and survey value shown; "
                 "verdict in each panel title", fontsize=9)
    fig.subplots_adjust(left=0.06, right=0.99, top=0.80, bottom=0.12, wspace=0.28)
    fig.savefig(path, dpi=FIG_DPI)
    plt.close(fig)


def figures_all(ctx):
    fig_map(ctx)
    fig_profiles(ctx)
    fig_bands(ctx, "T")
    fig_bands(ctx, "P")
    sec = section_setup(ctx)
    ctx["sec"] = sec
    fig_section(ctx, sec, "T", os.path.join(CALF, "section_A.png"))
    fig_section(ctx, sec, "P", os.path.join(CALF, "section_A_P.png"))
    pt_figure(ctx, os.path.join(CALF, "PT_state.png"))
    ctx["concept_key"] = fig_conceptual(ctx, sec, os.path.join(CALF, "conceptual_model.png"))
    fig_gauge(ctx, os.path.join(CALF, "gauge_check.png"))
    # interpolated companions (8.4)
    fig_section(ctx, sec, "T", os.path.join(ADDF, "section_A_interpolated.png"), interp=True)
    fig_section(ctx, sec, "P", os.path.join(ADDF, "section_A_P_interpolated.png"), interp=True)
    interp_band_maps(ctx, "T")
    interp_band_maps(ctx, "P")


def fig_map(ctx):
    fig, ax = plt.subplots(figsize=(9.5, 8.0))
    sec = section_setup(ctx)
    lines = {"A|A'": (sec["A"] - sec["u"] * 3000 * 0, sec["A"])}
    s0, s1 = sec["s0"], sec["s1"]
    p0 = sec["A"] + s0 * sec["u"]
    p1 = sec["A"] + s1 * sec["u"]
    plot_base(ax, ctx, labels=True, zoom=None, lines={"A|A'": (p0, p1)})
    ax.set_title("Base map: wells, faults, water bodies, model grid and section line A-A' "
                 "(through the two hottest wells, extended to the model extent)", fontsize=8)
    handles = base_legend_handles()
    fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.80, 0.5), fontsize=7, frameon=False)
    fig.subplots_adjust(left=0.09, right=0.78, top=0.93, bottom=0.08)
    fig.savefig(os.path.join(CALF, "well_map.png"), dpi=FIG_DPI)
    plt.close(fig)
    ctx["section_line_ends"] = (p0, p1)


def index_rows(ctx):
    rows = []
    rows.append(("F1", rel(os.path.join(CALF, "well_map.png")), "wells, faults, water bodies, grid, section line A-A'",
                 ",".join(ctx["wnames"]), "map", "Baseline mapping, geometry, well tracks", "obs"))
    for w in ctx["wnames"]:
        rows.append(("F2", rel(os.path.join(CALF, f"profile_{w}.png")), f"T and P profiles, shoe, feed zones, Tsat, verdict for {w}",
                     w, "0-2350 mRL span shown on fixed axis", "survey T and P, well JSON", "obs"))
    for grp in ctx["groups"]:
        rows.append(("F2", rel(os.path.join(CALF, f"profile_panels_{grp}.png")), f"temperature profiles of group {grp} on common axes",
                     ",".join(w for w in ctx["groups"][grp] if w in ctx["wells"]), "elevation", "survey T", "obs"))
    rows.append(("F3", rel(os.path.join(CALF, "maps_temperature_bands.png")), "temperature minus field median in 500 m bands",
                 ",".join(ctx["wnames"]), "elevation bands", "T_P_at_elevation_bands.csv", "obs"))
    rows.append(("F3", rel(os.path.join(CALF, "maps_pressure_bands.png")), "pressure minus gauge-verified trend in 500 m bands",
                 ",".join(ctx["wnames"]), "elevation bands", "T_P_at_elevation_bands.csv", "obs"))
    sec = ctx["sec"]
    rows.append(("F4", rel(os.path.join(CALF, "section_A.png")), "temperature section through the two hottest wells; data envelope; no-data hatched",
                 ",".join(sorted(sec["proj"].keys())), f"section A-A' (end points listed in the report)", "NS_calibration_targets.csv", "mixed"))
    rows.append(("F4", rel(os.path.join(CALF, "section_A_P.png")), "pressure deviation from trend, same section and envelope",
                 ",".join(sorted(sec["proj"].keys())), "section A-A'", "NS_calibration_targets.csv", "mixed"))
    rows.append(("F5", rel(os.path.join(CALF, "PT_state.png")), "P-T state of survey points with the IF97 saturation curve",
                 ",".join(ctx["wnames"]), "all survey points", "NS_calibration_targets.csv, feedzones.csv", "obs"))
    rows.append(("F6", rel(os.path.join(CALF, "conceptual_model.png")), "conceptual model on F1 (plan) and F4 (section)",
                 ",".join(ctx["wnames"]), "plan and section A-A'", "well_summary.csv, feedzones.csv, F1, F4", "mixed"))
    gw = [w for w in ctx["wnames"] if ctx["ns"][w]["gauge"]]
    if gw:
        rows.append(("F7", rel(os.path.join(CALF, "gauge_check.png")), "transient gauge against survey pressure, production timing",
                     ",".join(gw), "gauge elevations as given in the report", "transient records, flow records, survey", "obs"))
    return rows


def figure_index(ctx):
    rows = index_rows(ctx)
    df = pd.DataFrame(rows, columns=["figure_id", "file", "shows", "wells", "elevation_or_line", "data_sources", "claim_type"])
    write_csv(df, os.path.join(CAL, "figures_index.csv"))


def gauge_at_date(series, date_iso):
    """Gauge value at a date, linear in time between readings; blank outside the record."""
    t = [dtm.date.fromisoformat(d).toordinal() for d, _ in series]
    v = [x for _, x in series]
    x = dtm.date.fromisoformat(date_iso).toordinal()
    if x < t[0] or x > t[-1]:
        return float("nan")
    return float(np.interp(x, t, v))


def extra_tables(ctx):
    """Supporting tables used by the report (all numbers quoted in the report come from these files)."""
    wells, wnames, tgt, ns = ctx["wells"], ctx["wnames"], ctx["tgt"], ctx["ns"]
    rows = []
    for z in [0.0, -250.0, -500.0, -1000.0, -1500.0, -2000.0]:
        for w in wnames:
            W = wells[w]
            if W["tz"].min() - 1e-9 <= z <= W["tz"].max() + 1e-9:
                T = interp_profile(W["tz"], W["tt"], z)
                P = interp_profile(W["pz"], W["pp"], z)
                rows.append(dict(elevation_mRL=z, well=w, T_degC=round(T, 2), P_bar=round(P, 2),
                                 P_minus_trend_bar=round(P - float(trend_p(ctx, z)), 2)))
    write_csv(pd.DataFrame(rows), os.path.join(ADD, "T_P_at_common_elevations.csv"))
    tr = ctx["trend"]
    par = pd.read_csv(os.path.join(ADD, "analysis_parameters.csv"))
    extra = pd.DataFrame([
        ("verified trend intercept a (P = a + b x depth)", tr["a"], "bar", "fit to gauge-verified feed zones (§4.1)"),
        ("verified trend slope b", tr["b"], "bar/m", "fit to gauge-verified feed zones (§4.1)"),
        ("verified trend RMS residual", tr["rms"], "bar", "gauge-verified feed zones about the trend"),
        ("verified trend tolerance applied", tr["tol"], "bar", "max(3, 2 x RMS) (§4.1)"),
    ], columns=par.columns)
    write_csv(round_df(pd.concat([par, extra], ignore_index=True), 6), os.path.join(ADD, "analysis_parameters.csv"))

    hy = ctx["hy"]
    if len(hy):
        g = hy.groupby(["well", "band"], sort=False).agg(
            n_windows=("ratio_measured_to_IF97", "size"),
            ratio_min=("ratio_measured_to_IF97", "min"),
            ratio_median=("ratio_measured_to_IF97", "median"),
            ratio_max=("ratio_measured_to_IF97", "max"),
            rho_implied_min_kgm3=("rho_implied_kgm3", "min"),
            rho_implied_max_kgm3=("rho_implied_kgm3", "max"),
            excess_median_pct=("excess_pct_over_cold_component", "median"),
        ).reset_index()
        g["well_order"] = g.well.map({w: i for i, w in enumerate(wnames)})
        g = g.sort_values(["well_order", "band"]).drop(columns="well_order")
        write_csv(round_df(g, 4), os.path.join(ADD, "hydrostatic_by_well_band.csv"))

    rows = []
    for w in wnames:
        W = wells[w]
        n = ns[w]
        rows.append(dict(well=w, kind=W["kind"], injectate_enthalpy_kJkg=n["inj_h"] if n["inj_h"] is not None else "",
                         injectate_T_degC=n["inj_T"], formation_median_T_at_feedzones_degC=n["form_T"],
                         formation_hotter_than_injectate=n["hotter"] if W["kind"] == "injector" else "",
                         shoe_mRL=W["shoe"], first_flow_date=W["first_flow"] or "",
                         open_hole_isothermal_points=int(tgt[(tgt.well == w) & tgt._iso].shape[0]),
                         open_hole_isothermal_top_mRL=(tgt[(tgt.well == w) & tgt._iso].elevation_mRL.max()
                                                       if tgt[(tgt.well == w) & tgt._iso].shape[0] else ""),
                         open_hole_isothermal_base_mRL=(tgt[(tgt.well == w) & tgt._iso].elevation_mRL.min()
                                                        if tgt[(tgt.well == w) & tgt._iso].shape[0] else ""),
                         open_hole_isothermal_T_min_degC=(tgt[(tgt.well == w) & tgt._iso].T_degC.min()
                                                          if tgt[(tgt.well == w) & tgt._iso].shape[0] else ""),
                         open_hole_isothermal_T_max_degC=(tgt[(tgt.well == w) & tgt._iso].T_degC.max()
                                                          if tgt[(tgt.well == w) & tgt._iso].shape[0] else "")))
    inj = pd.DataFrame(rows)
    for c in ["injectate_T_degC", "formation_median_T_at_feedzones_degC", "open_hole_isothermal_T_min_degC",
              "open_hole_isothermal_T_max_degC"]:
        inj[c] = inj[c].map(lambda v: round(float(v), 2) if v not in ("", None) and np.isfinite(float(v)) else "")
    write_csv(inj, os.path.join(ADD, "injector_and_isothermal_check.csv"))

    wr = tgt.groupby("weight_reason", sort=False).agg(n_points=("weight", "size"), weight_min=("weight", "min"),
                                                      weight_max=("weight", "max")).reset_index()
    write_csv(round_df(wr, 3), os.path.join(ADD, "weights_summary.csv"))


def write_report(ctx):
    """Generates Calibration_Input/NS_data_interpretation.md from the output tables."""
    wells, wnames = ctx["wells"], ctx["wnames"]
    ws = pd.read_csv(os.path.join(CAL, "well_summary.csv"))
    pc = pd.read_csv(os.path.join(ADD, "pressure_consistency_check.csv"))
    hy = pd.read_csv(os.path.join(ADD, "hydrostatic_check.csv"))
    hb = pd.read_csv(os.path.join(ADD, "hydrostatic_by_well_band.csv"))
    hc = pd.read_csv(os.path.join(ADD, "hydrostatic_cold_component.csv"))
    geo_chk = pd.read_csv(os.path.join(ADD, "well_geometry_check.csv"))
    fz = pd.read_csv(os.path.join(CAL, "feedzones.csv"))
    tgt = pd.read_csv(os.path.join(CAL, "NS_calibration_targets.csv"))
    inj = pd.read_csv(os.path.join(ADD, "injector_and_isothermal_check.csv"))
    wr = pd.read_csv(os.path.join(ADD, "weights_summary.csv"))
    cem = pd.read_csv(os.path.join(ADD, "T_P_at_common_elevations.csv"))
    xc = pd.read_csv(os.path.join(ADD, "data_crosscheck.csv"))
    inv = pd.read_csv(os.path.join(ADD, "input_inventory.csv"))
    par = pd.read_csv(os.path.join(ADD, "analysis_parameters.csv"))
    trend = ctx["trend"]
    ns = ctx["ns"]
    key = ctx["concept_key"]
    write_csv(pd.DataFrame(key, columns=["element_no", "element", "well_or_location", "observation", "claim_type", "drawing"]), os.path.join(ADD, "conceptual_model_key.csv"))
    sec = ctx["sec"]
    A_name, B_name = sec["wa"], sec["wb"]
    p0 = sec["A"] + sec["s0"] * sec["u"]
    p1 = sec["A"] + sec["s1"] * sec["u"]
    inj_df = inj[inj.kind == "injector"]
    inj_hot = inj_df[inj_df.formation_hotter_than_injectate == True].well.tolist()  # noqa: E712
    inj_cold = [w for w in inj_df.well if w not in inj_hot]
    ws_i = ws.set_index("well")
    pc_i = pc.set_index("well")
    g_wells = [w for w in wnames if ns[w]["gauge"]]
    verdict_groups = {}
    for w in wnames:
        verdict_groups.setdefault((ws_i.loc[w, "ns_verdict"], ws_i.loc[w, "ns_confidence"]), []).append(w)
    def wl(ws_list):
        return ", ".join(ws_list) if ws_list else "none"
    up = ws_i[ws_i.regime == "upflow"].index.tolist()
    outf = ws_i[ws_i.regime == "outflow"].index.tolist()
    cold = ws_i[ws_i.regime == "recharge/downflow"].index.tolist()
    conduc = ws_i[ws_i.regime == "conductive"].index.tolist()
    inv_wells = ws_i[ws_i.inversion_degC.notna()].index.tolist()
    comp = {r.band: r.median_ratio for _, r in hc.iterrows()}
    comp_cased = comp.get("cased", float("nan"))
    comp_open = comp.get("open", float("nan"))
    wh_below = geo_chk[geo_chk.wellhead_minus_surface_m < 0]
    wh_rng = (wh_below.wellhead_minus_surface_m.min(), wh_below.wellhead_minus_surface_m.max()) if len(wh_below) else (0, 0)
    w2 = pc_i.loc["WELL-2"] if "WELL-2" in pc_i.index else None
    w3 = pc_i.loc["WELL-3"] if "WELL-3" in pc_i.index else None
    gauge_w = pc_i[pc_i.gauge_test.isin(["pass", "miss"])]
    cased3 = hy[(hy.well == "WELL-3") & (hy.band == "cased")]
    worst_cased = hy[hy.band == "cased"].sort_values("rho_implied_kgm3", ascending=False).iloc[0]
    open_hb = hb[hb.band == "open"]
    cased_hb = hb[hb.band == "cased"]
    straddle = hy[hy.band == "straddles shoe"]
    fz_off = fz.copy()
    tlist = trend["verified"]
    first_prod = ctx.get("field_first", "")
    L = []
    a = L.append

    a("# Natural-state data interpretation: Teaching_Model_3D (TE14515)")
    a("")
    a("Generated by `Scripts/natural_state_analysis.py` from `Baseline/` (read-only). Every number below is read from a table in this folder or in `../Additional_Output/`. Tags: `[obs]` measured in the data, `[interp]` an interpretation of measured data, `[assumption]` an assumed value, `[reference model]` content of the starting model (never evidence).")
    a("")
    a("## 1. Key conclusions for calibration")
    a("")
    n = 1
    a(f"{n}. **Natural-state verdicts.** Of {len(wnames)} wells, {len(verdict_groups.get(('natural state','high'), []))} are natural state with high confidence ({wl(verdict_groups.get(('natural state','high'), []))}); {len(verdict_groups.get(('probably natural state','medium'), []))} are probably natural state with medium confidence ({wl(verdict_groups.get(('probably natural state','medium'), []))}); {wl(verdict_groups.get(('probably natural state at best','low'), []))} is probably natural state at best, low confidence; and {wl(verdict_groups.get(('possibly disturbed','low'), []))} is possibly disturbed, low confidence. [interp] Source: `pressure_consistency_check.csv`, rules in section 2.")
    a(f"   **Action:** use the natural-state targets at the weights in `NS_calibration_targets.csv`; treat {wl(verdict_groups.get(('possibly disturbed','low'), []))} and {wl(verdict_groups.get(('probably natural state at best','low'), []))} as low-weight targets only.")
    a("")
    n += 1
    if w2 is not None and w3 is not None:
        a(f"{n}. **The two gauge-verified wells match their first gauge reading, but the survey date conflicts with the gauge.** WELL-2 survey pressure at {w2.gauge_elev_mRL:.0f} mRL is {w2.survey_P_at_gauge_bar:.1f} bar, against a first gauge reading of {w2.first_gauge_P_bar:.1f} bar on {w2.first_gauge_date}; WELL-3 differs by {w3.diff_survey_minus_first_bar:.2f} bar. [obs] On the survey label date ({wells['WELL-2']['survey_date']}) the WELL-2 gauge read {w2.gauge_P_at_survey_date_bar:.1f} bar, so the survey agrees with the first reading and not with the gauge at its own date. [interp] The gauge test decides the verdict (§4.1), so both wells stay at high confidence; the date label is a data question.")
        a(f"   **Action:** confirm the survey date for WELL-2 and WELL-3. If the survey was taken at the start of production, these two surveys are the best natural-state anchors in the dataset.")
    a("")
    n += 1
    a(f"{n}. **Cased-section pressures disagree with the temperature profiles.** Above each production casing shoe, the measured pressure gradient is steeper than the IF97 column for the measured temperatures in every well that has a cased window: the cased-band ratio runs from {cased_hb.ratio_min.min():.2f} to {cased_hb.ratio_max.max():.2f}, against a cased component of {comp_cased:.2f} (one cold window, WELL-8). [obs] The steepest cased window ({worst_cased.well}, {worst_cased.z_top_mRL:.0f} to {worst_cased.z_bottom_mRL:.0f} mRL) implies a liquid density of {worst_cased.rho_implied_kgm3:.0f} kg/m3, which is not possible for water. [interp] A datum, gravity or unit convention would act at every depth and cannot explain this, and two-phase fluid would make the column lighter; the cased pressures are most likely the wellbore column (§4.2), not formation pressure.")
    a("   **Action:** do not use cased-section pressures as formation targets. Use feed-zone pressures only (`feedzones.csv`).")
    a("")
    n += 1
    big_open = open_hb[open_hb.excess_median_pct > 3.0].well.tolist()
    a(f"{n}. **Open-hole pressures agree with IF97 for most wells, apart from a field-wide offset.** The open-band cold-well component is {comp_open:.3f}, a field-wide offset of about {100 * (comp_open - 1):.1f} percent. [obs] Open-hole median excesses above 3 percent occur in {wl(big_open)} ([hydrostatic_by_well_band.csv](../Additional_Output/hydrostatic_by_well_band.csv)). [obs] A datum or convention offset is the first thing to check (§4.3). Flowing wells could also carry bore friction in the measured gradient. [interp]")
    a("   **Action:** confirm the gravity and unit convention and the elevation datum used for the gauges before calibrating to open-hole pressure gradients; carry the component as a known offset.")
    a("")
    n += 1
    if up:
        upw = up[0]
        a(f"{n}. **Hot upflow at {upw}.** {upw} has the highest temperature in the field ({ws_i.loc[upw,'Tmax_degC']:.1f} degC at {ws_i.loc[upw,'Tmax_elev_mRL']:.1f} mRL), no inversion, and a minimum boiling margin of {ws_i.loc[upw,'min_boiling_margin_degC']:.1f} degC at {ws_i.loc[upw,'min_boiling_margin_elev_mRL']:.0f} mRL, so the top of the upflow is near boiling. [obs] Regime: upflow by the Tmax rule (§4.4). [interp]")
        a(f"   **Action:** place the heat source beneath {upw}, with a near-boiling target at {ws_i.loc[upw,'min_boiling_margin_elev_mRL']:.0f} mRL and a feed-zone temperature target at {fz[fz.well==upw].mid_mRL.iloc[0]:.1f} mRL.")
    a("")
    n += 1
    a(f"{n}. **Lateral outflow at depth.** Measured inversions (temperature drop below the running maximum, at least 3 degC) occur in {wl(inv_wells)}, with tops between {ws_i.inversion_top_mRL.min():.0f} and {ws_i.inversion_top_mRL.max():.0f} mRL. [obs] Regime outflow for {wl(outf)} (inversion at least 10 degC). [interp] Lateral outflow away from the upflow at these depths.")
    a("   **Action:** put permeable horizons at the inversion depths listed in `well_summary.csv`, and calibrate lateral permeability to reproduce those inversions.")
    a("")
    n += 1
    a(f"{n}. **Cold recharge or downflow in the cold wells.** {wl(cold)} have Tmax below 100 degC ({', '.join(f'{ws_i.loc[w,'Tmax_degC']:.1f}' for w in cold)} degC). [obs] Both wells have Tmax below 100 degC. [interp] Cold downflow or recharge on that side; the regime does not show the direction of flow.")
    a("   **Action:** set cold recharge boundaries near these wells at their measured temperatures; the open-hole isothermal section of WELL-8 is listed in `injector_and_isothermal_check.csv`.")
    a("")
    n += 1
    a(f"{n}. **Wellheads sit below the model surface.** In {len(wh_below)} wells the wellhead is {abs(wh_rng[1]):.0f} to {abs(wh_rng[0]):.0f} m below the model surface of its column. [obs] {int(tgt.above_model_surface.astype(bool).sum())} survey points lie above the model surface. [obs] The grid topography and the wellhead elevations do not agree. [interp]")
    a("   **Action:** reconcile the model topography with the wellhead elevations before using near-surface temperatures.")
    a("")
    n += 1
    a(f"{n}. **The cap is not resolved by the survey data.** No interval of conductive gradient is identified above the upflow that can be separated from the outflow profile. [interp] The cap is therefore drawn as absent in F6. What would decide it: a shallow vertical profile in a well on the upflow with a linear conductive gradient above the 100 degC isotherm, or a lithology log of the caprock.")
    a("   **Action:** calibrate cap thickness and permeability contrast as free parameters, with the outflow inversions as the constraint.")
    a("")
    a("## 2. The data")
    a("")
    a(f"**Inventory.** {len(inv)} inventory rows classify every file under `Baseline/` (see [input_inventory.csv](../Additional_Output/input_inventory.csv)). Geometry and block assignment come from the grid file and the well-to-block file. Model properties (reference rock, boundary conditions, results, set-up) are reference only and were not used as evidence. Production and future-scenario files are out of scope. Field production timing was taken from the well flow records.")
    a("")
    a(f"**Checks.** Plain-column copies were compared with the JSON for every survey and transient point ([data_crosscheck.csv](../Additional_Output/data_crosscheck.csv)): {int((xc.status == 'match').sum())} series match exactly and {int((xc.status == 'mismatch').sum())} disagree. [obs] Flow and enthalpy series were not cross-checked. Irregular file names: `Wel10_*` for WELL-10. Survey dates in use: {', '.join(sorted(set(v['survey_date'] for v in wells.values())))}; each well has one survey date. [obs]")
    a("")
    a("**Provenance labels.**")
    srcs = sorted({wells[w]["source_notes"]["design_source"] for w in wnames if wells[w]["source_notes"]["design_source"]})
    a("- Design source field of every well reads " + ", ".join(f"'{x}'" for x in srcs) + ". [obs]")
    a("- Feed-zone quality tags found: " + ", ".join(sorted(fz.quality.unique())) + f". The WELL-10 tag 'Synthetic' is not a quality grade. [obs]")
    a("- Feed-zone top and bottom are labelled 'm' in the JSON but are elevations in mRL (§3). [obs]")
    a("- Survey StableFlag is true for every survey; transient records carry no flag. [obs]")
    a("")
    a("**Natural-state timing for each well with a gauge** ([pressure_consistency_check.csv](../Additional_Output/pressure_consistency_check.csv)):")
    a("")
    a("| Well | Gauge elevation (mRL) | First reading (date, bar) | Own first flow | Pre-production reading | Survey minus first reading (bar) | Verdict | Confidence |")
    a("|---|---|---|---|---|---|---|---|")
    for w in g_wells:
        r = pc_i.loc[w]
        a(f"| {w} | {r.gauge_elev_mRL:.0f} | {r.first_gauge_date}, {r.first_gauge_P_bar:.1f} | {r.own_first_flow_date} | {'yes' if r.pre_production_reading else 'no'} | {r.diff_survey_minus_first_bar:+.2f} | {r.verdict} | {r.confidence} |")
    a("")
    a(f"Wells without a gauge: {wl([w for w in wnames if not ns[w]['gauge']])}. For these, the feed-zone pressures were tested against the verified trend: P = {trend['a']:.2f} + {trend['b']:.5f} x depth bar, fitted to the gauge-verified feed zones of {wl(tlist)} (RMS {trend['rms']:.2f} bar). The tolerance is {trend['tol']:.1f} bar, the floor, because twice the RMS is smaller. [interp] Figure F7 shows the gauge test: ![F7](figures/gauge_check.png) [obs]")
    a("")
    a("Figure F7 (gauge check) [obs]: transient pressure against time, survey pressure at the gauge elevation, first production in the field, and flowing periods.")
    a("")
    a("**Verdict rules applied (§4.1).** Rule 1: survey matches a pre-production gauge reading within the precision. Rule 2: injector whose formation is hotter than its injectate. Rule 3: gauge misses by more than the precision, feed zones on the trend. Rule 4: no pre-production gauge reading, feed zones on the trend. Rule 5: feed-zone pressures off the trend by more than the tolerance. Rules applied per well: " + "; ".join(f"{w} {ns[w]['rule'].split(':')[0]}" for w in wnames) + ". [interp]")
    a("")
    a("Figure F1 (base map) [obs]: wells, faults, water bodies, model grid and section line A-A'. ![F1](figures/well_map.png)")
    a("")
    a("## 3. Temperature structure")
    a("")
    a("| Well | Tmax (degC) | at (mRL) | T at TD (degC) | Inversion (degC) | Regime | 100 degC at (mRL) | Min boiling margin (degC) | Verdict | Confidence |")
    a("|---|---|---|---|---|---|---|---|---|---|")
    for w in wnames:
        r = ws_i.loc[w]
        inv_s = f"{r.inversion_degC:.1f}" if np.isfinite(r.inversion_degC) else "none"
        z100 = f"{r.z_100C_mRL:.0f}" if np.isfinite(r.z_100C_mRL) else "not reached"
        a(f"| {w} | {r.Tmax_degC:.1f} | {r.Tmax_elev_mRL:.0f} | {r.T_at_TD_degC:.1f} | {inv_s} | {r.regime} | {z100} | {r.min_boiling_margin_degC:.1f} | {r.ns_verdict} | {r.ns_confidence} |")
    a("")
    a("[obs] Values are from `well_summary.csv`. Inversions, Tmax and isotherm crossings follow the §4.4 definitions.")
    a("")
    a(f"Profiles and their shoe, feed-zone and boiling-margin annotations: figures [F2](figures/profile_panels_All.png) (all wells on common axes), and one profile per well in `figures/profile_<well>.png`. [obs]")
    a("")
    a(f"Lateral temperature structure at depth bands without interpolation: [F3 temperature](figures/maps_temperature_bands.png). [obs] Section along the line through the two hottest wells: [F4](figures/section_A.png), the line runs through {A_name} and {B_name}, from E{p0[0]:.0f} N{p0[1]:.0f} to E{p1[0]:.0f} N{p1[1]:.0f} (NZTM, m). Offsets of wells from the section line are printed on the figure. [obs] Interpolated companion (not data, not used for targets): [section](../Additional_Output/figures/section_A_interpolated.png) and [plan](../Additional_Output/figures/maps_temperature_interpolated.png). [interp]")
    a("")
    a("## 4. Pressure structure")
    a("")
    a("Pressure at common elevations, interpolated within each survey and never extrapolated ([T_P_at_common_elevations.csv](../Additional_Output/T_P_at_common_elevations.csv)):")
    a("")
    a("| Elevation (mRL) | Well | T (degC) | P (bar) | P minus trend (bar) |")
    a("|---|---|---|---|---|")
    for _, r in cem.iterrows():
        a(f"| {r.elevation_mRL:.0f} | {r.well} | {r.T_degC:.1f} | {r.P_bar:.1f} | {r.P_minus_trend_bar:+.1f} |")
    a("")
    a("[obs] Hydrostatic checks by band (4.3; [hydrostatic_by_well_band.csv](../Additional_Output/hydrostatic_by_well_band.csv)):")
    a("")
    a("| Well | Band | Windows | Ratio min | Ratio median | Ratio max | Median excess over cold component (%) |")
    a("|---|---|---|---|---|---|---|")
    for _, r in hb.iterrows():
        exc = f"{r.excess_median_pct:+.1f}" if pd.notna(r.excess_median_pct) else ""
        a(f"| {r.well} | {r.band} | {int(r.n_windows)} | {r.ratio_min:.3f} | {r.ratio_median:.3f} | {r.ratio_max:.3f} | {exc} |")
    a("")
    a(f"Cold-well component (median of the cold wells' windows): cased {comp_cased:.3f}, open {comp_open:.3f} ([hydrostatic_cold_component.csv](../Additional_Output/hydrostatic_cold_component.csv)). [interp] Windows that straddle a casing shoe, because a remnant joined the window above it (§4.3), are reported as 'straddles shoe' and are in neither band; there are {len(straddle)} such windows. [obs] Each ratio has a rounding uncertainty of about 0.1 bar over the window's pressure difference ([hydrostatic_check.csv](../Additional_Output/hydrostatic_check.csv)). [obs]")
    a("")
    a("Pressure versions of the sections and band maps: [F3 pressure](figures/maps_pressure_bands.png), [F4 pressure](figures/section_A_P.png). [obs] Interpolated companions: [section](../Additional_Output/figures/section_A_P_interpolated.png) and [plan](../Additional_Output/figures/maps_pressure_interpolated.png). [interp]")
    a("")
    a("Figure F5 (P-T state) [obs]: survey points of pressure against temperature with the IF97 saturation curve, ending at the critical point; feed zones marked. ![F5](figures/PT_state.png)")
    a("")
    a("## 5. Conceptual model")
    a("")
    a("Figure F6 draws the model on the base map (F1) and on the section (F4). Each element is numbered in the figure and listed below with its supporting observation and claim type. Solid arrows are measured; dashed arrows are inferred ([conceptual_model.png](figures/conceptual_model.png), [F6](figures/conceptual_model.png)).")
    a("")
    a("| No. | Element | Well or location | Supporting observation | Claim type | Drawing |")
    a("|---|---|---|---|---|---|")
    for el, name, loc, obs, claim, symb in key:
        a(f"| {el} | {name} | {loc} | {obs} | [{claim}] | {symb} |")
    a("")
    if up:
        gu = geo_chk.set_index("well").loc[up[0]]
        a(f"Relation to the mapping: the upflow well {up[0]} has its nearest fault at {gu.nearest_fault} ({gu.dist_nearest_fault_m:.0f} m) and its nearest water body at {gu.nearest_water_body} ({gu.dist_nearest_water_body_m:.0f} m) ([well_geometry_check.csv](../Additional_Output/well_geometry_check.csv)). [obs] Whether the faults or the water bodies control the upflow or the outflow cannot be decided from the survey data. [interp] What would decide it: a cross-fault pressure test, or a temperature survey run across the fault.")
    a("")
    a("## 6. Recommended calibration set-up")
    a("")
    a("**Targets and weights.** Temperature targets are every survey point (`NS_calibration_targets.csv`, `weight`); pressure targets are the feed zones only (`feedzones.csv`, `weight_P`). The weight counts by rule are in [weights_summary.csv](../Additional_Output/weights_summary.csv). [obs] The confidence factor changes relative weights because the confidence differs between wells (verdicts above). [interp]")
    a("")
    a("**Where heat and mass enter and leave.** Heat enters beneath the upflow well (section 5, element 1). Mass and heat leave laterally at the inversion depths (elements 2 to 7). Cold recharge enters at the cold wells (elements 8 and 9). [interp]")
    a("")
    a("**Permeability structure the data imply.** A high-permeability path at each inversion depth connecting the upflow to the outflow wells; low-permeability elsewhere above the upflow. The cap is not resolved (section 1, conclusion 9). [interp]")
    a("")
    a(f"**Boundary requirements.** Cold recharge at the cold wells; a fixed-pressure or fixed-temperature condition at the model top, with the topography reconciled first (section 1, conclusion 8). Injectate temperatures are in [injector_and_isothermal_check.csv](../Additional_Output/injector_and_isothermal_check.csv). The formation is hotter than the injectate at {wl(inj_hot)} and colder at {wl(inj_cold)}. [obs] Where the formation is hotter, injection cooling could be hidden and the confidence is lowered (§4.1). [interp]")
    a("")
    a("## 7. Data quality issues and open questions (ranked by effect on calibration)")
    a("")
    a("1. Survey date for WELL-2 and WELL-3 conflicts with the gauge (section 1, conclusion 2). [obs] **Action:** confirm the date.")
    a("2. Cased-section pressures imply densities that are impossible for water (section 1, conclusion 3). [obs] **Action:** exclude them from formation targets.")
    a(f"3. Wellheads below the model surface in {len(wh_below)} wells (section 1, conclusion 8). [obs] **Action:** reconcile the topography.")
    if g_wells:
        r9 = pc_i.loc["WELL-9"] if "WELL-9" in pc_i.index else None
        if r9 is not None and r9.gauge_test == "miss":
            a(f"4. WELL-9 gauge misses its first reading by {r9.diff_survey_minus_first_bar:+.2f} bar, more than the {GAUGE_PREC:.1f} bar precision, although its feed zones lie on the trend. [obs] **Action:** check the gauge elevation and the survey depth.")
    r10 = pc_i.loc["WELL-10"] if "WELL-10" in pc_i.index else None
    if r10 is not None:
        a(f"5. WELL-10 gauge starts on {r10.first_gauge_date}, after its own first flow on {r10.own_first_flow_date}; its survey predates its production. [obs] The first readings of WELL-2, WELL-3 and WELL-9 are dated on their own first production day and are treated as pre-production. [assumption] **Action:** obtain an early gauge record if one exists.")
    r1 = ws_i.loc["WELL-1"]
    a(f"6. WELL-1 has a feed zone {fz[(fz.well=='WELL-1')].iloc[0].P_bar - float(trend['a'] + trend['b'] * (-fz[(fz.well=='WELL-1')].iloc[0].mid_mRL)):+.1f} bar from the trend, outside the tolerance, and its extrapolated water level is {r1.water_level_mRL:.0f} mRL, above its wellhead. [obs] **Action:** treat WELL-1 pressure as low confidence.")
    a(f"7. The gauge reading precision of {GAUGE_PREC:.1f} bar is an assumption. [assumption] All verdicts that depend on it should be rechecked when the instrument specification is known.")
    a("8. The cold-well component rests on WELL-8 alone: the only WELL-4 window straddles the shoe, so it is excluded. [obs] **Action:** treat the datum offset as provisional.")
    a("9. WELL-5 is an injector whose formation is hotter than the injectate; cooling could be hidden. [obs] **Action:** keep its weights low.")
    a("10. The WELL-10 feed-zone tag 'Synthetic' is not a quality grade and is weighted at the lowest factor. [obs] **Action:** confirm the quality of that zone.")
    a("11. Irregular file names (`Wel10_*`) and empty provenance fields. [obs] **Action:** none for calibration; note for data handling.")
    a("")
    a("## 8. Output file index")
    a("")
    a("**Output/Calibration_Input (the calibration set):**")
    a("")
    for f, d in [("NS_data_interpretation.md", "this report"), ("NS_calibration_targets.csv", "one row per survey point with weights"),
                 ("feedzones.csv", "feed zones with survey T, P and weights"), ("T_P_at_layer_centres.csv", "T and P at model layer centres (no extrapolation)"),
                 ("well_summary.csv", "per-well characteristics and verdicts"), ("figures_index.csv", "index of figures F1 to F7"),
                 ("figures/", "standard figures F1 to F7")]:
        a(f"- [{f}]({f if not f.endswith('/') else f}) : {d}")
    a("")
    a("**Output/Additional_Output (supporting evidence):**")
    a("")
    for f, d in [("input_inventory.csv", "role of every Baseline file"), ("data_crosscheck.csv", "plain-column against JSON"),
                 ("pressure_consistency_check.csv", "gauge tests, trend tolerance and rules"),
                 ("hydrostatic_check.csv", "hydrostatic windows"), ("hydrostatic_by_well_band.csv", "hydrostatic summary by well and band"),
                 ("hydrostatic_cold_component.csv", "cold-well component"), ("well_geometry_check.csv", "wellhead, surface and feature distances"),
                 ("T_P_at_elevation_bands.csv", "data for F3"), ("T_P_at_common_elevations.csv", "T and P at common elevations"),
                 ("injector_and_isothermal_check.csv", "injectate and isothermal sections"), ("weights_summary.csv", "weight counts by rule"),
                 ("analysis_parameters.csv", "every threshold used"), ("figures/", "interpolated companions and supporting figures"),
                 ("run_console_log.txt", "console output of the final run"), ("acceptance_check.txt", "acceptance check output"),
                 ("process_log.md", "prompt, steps, checks and departures")]:
        a(f"- [{f}](../Additional_Output/{f}) : {d}")
    a("")
    text = "\n".join(L) + "\n"
    with open(os.path.join(CAL, "NS_data_interpretation.md"), "w", encoding="utf-8", newline="\n") as fh:
        fh.write(text)


if __name__ == "__main__":
    main()
