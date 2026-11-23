#!/usr/bin/env python
"""
Natural-state data interpretation for geothermal model calibration (TE14515 / Teaching_Model_3D).

Run from the repository root:
    python Scripts/natural_state_analysis.py

Reads only from Baseline/ (read-only). Writes every deliverable to Output/ and Output/figures/.
Rules followed: CLAUDE.md sections 1-4 and 7. Reference-only files (NaturalState results,
boundary conditions, rock types, Plotting, Parameters) are NOT read at all. Geometry is read only
from gTE14515.dat (VERTICES, GRID, LAYERS, SURFA, WELLS) and gTE14515_wells.json (well-to-block).
"""
import ast
import csv
import json
import math
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.path import Path as MPath  # noqa: E402
from matplotlib.tri import Triangulation, LinearTriInterpolator  # noqa: E402
from iapws import IAPWS97  # noqa: E402  (IAPWS-IF97 implementation: Region 4 for Tsat, Region 1 for liquid rho)

# ----------------------------------------------------------------------------- paths and constants
ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "Baseline" / "Teaching_Model_3D"
DATA_JSON = BASE / "Data" / "Data_json"
DATA_DAT = BASE / "Data" / "Data_dat"
MAPPING = BASE / "Data" / "Mapping" / "Mapping.json"
GEOM_DAT = BASE / "Models" / "TE14515" / "gTE14515.dat"
WELLS_JSON = BASE / "Models" / "TE14515" / "gTE14515_wells.json"
OUT = ROOT / "Output"
FIG = OUT / "figures"

G = 9.80665            # m/s^2
P_CRIT = 220.64        # bar, critical pressure of water
ATM = 1.01325          # bar
MARGIN_FLAG = 5.0      # K, boiling-margin flag threshold (CLAUDE 3.3)
HYDRO_MIN_LEN = 200.0  # m, minimum interval length for hydrostatic test (CLAUDE 3.3)
NEAR_WH_M = 50.0       # m below the top of the track counted as "near wellhead / top layer"
ISO_TOL = 1.0          # K, isothermal tolerance
ISO_MIN_LEN = 50.0     # m
CONF_FACTOR = {"high": 1.0, "medium": 0.6, "low": 0.3}
MAP_ELEVS = [200.0, 0.0, -250.0, -500.0, -750.0, -1000.0, -1500.0]
TRANSIENT_ELEVS_EXPECTED = {2: -750.0, 3: -550.0, 9: -675.0, 10: -250.0}
SECTION_WELLS = (10, 3)  # section line from WELL-10 to WELL-3 (two hottest wells)
SECTION_CORRIDOR_M = 400.0

LOG = []


def log(msg=""):
    msg = str(msg)
    print(msg)
    LOG.append(msg)


# ----------------------------------------------------------------------------- thermodynamics (IAPWS-IF97)
def tsat_c(p_bar):
    """Saturation temperature, deg C (IF97 Region 4). NaN at or above the critical pressure."""
    if not np.isfinite(p_bar) or p_bar >= P_CRIT:
        return float("nan")
    return IAPWS97(P=p_bar / 10.0, x=0).T - 273.15


def psat_bar(t_c):
    """Saturation pressure, bar, from temperature (IF97 Region 4)."""
    return IAPWS97(T=t_c + 273.15, x=0).P * 10.0


def rho_liq(p_bar, t_c):
    """Liquid density, kg/m^3 (IF97 Region 1). Only call for compressed liquid."""
    return IAPWS97(P=p_bar / 10.0, T=t_c + 273.15).rho


def validate_thermo():
    log("IF97 validation (CLAUDE 3.3: fail if any error > 0.1 K):")
    ok = True
    for p, ref in [(ATM, 99.97), (100.0, 311.0)]:
        calc = tsat_c(p)
        err = calc - ref
        log(f"  Tsat({p:g} bar) = {calc:.3f} C | reference {ref:.2f} C | error {err:+.4f} K")
        if abs(err) > 0.1:
            ok = False
    rho = rho_liq(ATM, 25.0)
    log(f"  rho_liquid(1.01325 bar, 25 C) = {rho:.2f} kg/m3 | reference 997.05 | error {rho - 997.05:+.2f}")
    if abs(rho - 997.05) > 0.5:
        ok = False
    if not ok:
        raise SystemExit("IF97 validation FAILED - stopping.")
    log("  IF97 validation PASSED.")


# ----------------------------------------------------------------------------- small helpers
def fmt(x, nd=2):
    if x is None:
        return ""
    if isinstance(x, (bool, np.bool_)):
        return "True" if x else "False"
    if isinstance(x, (int, np.integer)):
        return str(int(x))
    try:
        xf = float(x)
    except (TypeError, ValueError):
        return str(x)
    if not math.isfinite(xf):
        return ""
    return f"{xf:.{nd}f}"


def write_csv(path, header, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        for r in rows:
            w.writerow(r)
    log(f"  wrote {path.relative_to(ROOT).as_posix()} ({len(rows)} rows)")


def parse_date(s):
    return datetime.strptime(s, "%d/%m/%Y")


def decimal_year(d):
    start = datetime(d.year, 1, 1)
    return d.year + (d - start).days / 365.25


def interp_no_extrap(zs, vs, zq):
    """Linear interpolation in elevation. NaN outside the surveyed range (no extrapolation)."""
    zs = np.asarray(zs, float)
    vs = np.asarray(vs, float)
    o = np.argsort(zs)
    zs, vs = zs[o], vs[o]
    if zq < zs[0] - 1e-9 or zq > zs[-1] + 1e-9:
        return float("nan")
    return float(np.interp(zq, zs, vs))


def track_xy(track, z):
    """x,y of a well track at elevation z (track = list of (x,y,z) from top to bottom)."""
    for (x1, y1, z1), (x2, y2, z2) in zip(track[:-1], track[1:]):
        lo, hi = min(z1, z2), max(z1, z2)
        if lo - 1e-9 <= z <= hi + 1e-9:
            if abs(z2 - z1) < 1e-9:
                return x1, y1
            f = (z - z1) / (z2 - z1)
            return x1 + f * (x2 - x1), y1 + f * (y2 - y1)
    return float("nan"), float("nan")


def dist_to_polyline(px, py, pts):
    best = float("inf")
    p = np.array([px, py], float)
    for a, b in zip(pts[:-1], pts[1:]):
        a = np.asarray(a, float)
        b = np.asarray(b, float)
        ab = b - a
        denom = float(np.dot(ab, ab))
        t = 0.0 if denom == 0 else float(np.clip(np.dot(p - a, ab) / denom, 0.0, 1.0))
        d = float(np.linalg.norm(a + t * ab - p))
        best = min(best, d)
    return best


# ----------------------------------------------------------------------------- geometry (dat + wells json)
def parse_geometry():
    lines = GEOM_DAT.read_text(encoding="utf-8", errors="replace").splitlines()
    sections = {"VERTICES", "GRID", "CONNECTIONS", "LAYERS", "SURFA", "WELLS"}
    cur = None
    verts, layers, surfa, wells_geom = [], {}, {}, {}
    for ln in lines:
        s = ln.strip()
        if s in sections:
            cur = s
            continue
        if cur == "VERTICES" and len(ln) >= 23 and s:
            verts.append((ln[:3], float(ln[3:13]), float(ln[13:23])))
        elif cur == "LAYERS" and s:
            p = s.split()
            if len(p) == 3:
                layers[int(p[0])] = (min(float(p[1]), float(p[2])), max(float(p[1]), float(p[2])))
        elif cur == "SURFA" and s:
            p = s.split()
            if len(p) == 2:
                surfa[p[0]] = float(p[1])
        elif cur == "WELLS" and s:
            p = s.split()
            if len(p) == 4:
                wells_geom.setdefault(p[0], []).append((float(p[1]), float(p[2]), float(p[3])))
    log(f"Geometry parsed: {len(verts)} vertices, {len(layers)} layers, {len(surfa)} surface columns, "
        f"{len(wells_geom)} well tracks in WELLS section.")
    return verts, layers, surfa, wells_geom


def parse_mapping():
    m = json.loads(MAPPING.read_text(encoding="utf-8"))
    out = {}
    for name, obj in m.items():
        d = obj["Data"]
        if isinstance(d[0][0], (list, tuple)):
            polylines = [np.asarray(pl, float) for pl in d]
        else:
            polylines = [np.asarray(d, float)]
        out[name] = polylines
    log("Mapping loaded: " + ", ".join(f"{k} ({len(v)} polyline(s))" for k, v in out.items()))
    return out


# ----------------------------------------------------------------------------- well JSON
def load_well(n):
    name = f"WELL-{n}"
    d = json.loads((DATA_JSON / f"{name}_data.json").read_text(encoding="utf-8"))
    w = {"name": name, "aut2": d.get("AUT2_Name", "")}
    w["track"] = [tuple(float(v) for v in ast.literal_eval(s)) for s in d["Well_Track"]["Data"]]
    w["design"] = []
    for s in d["Well_Design"]["Data"]:
        top, bot, size, slotted, spec = ast.literal_eval(s)
        w["design"].append({"top": float(top), "bot": float(bot), "size": size,
                            "slotted": bool(slotted), "spec": spec})
    w["design_source"] = d["Well_Design"].get("Source", "")
    w["feedzones"] = []
    for s in d.get("Feed_Zones", {}).get("Data", []):
        top, bot, q = ast.literal_eval(s)
        w["feedzones"].append({"top": float(top), "bot": float(bot), "quality": q})
    # survey temperature and pressure (one survey date per well in these files)
    tkey = list(d["Downhole_Temperature"].keys())
    pkey = list(d["Downhole_Pressure"].keys())
    w["survey_date"] = tkey[0]
    if len(tkey) > 1 or len(pkey) > 1:
        log(f"  NOTE {name}: more than one survey date present; using the first.")
    tdat = d["Downhole_Temperature"][tkey[0]]["Data"]
    pdat = d["Downhole_Pressure"][pkey[0]]["Data"]
    w["T_raw"] = {float(k): float(v) for k, v in tdat.items()}
    w["P_raw"] = {float(k): float(v) for k, v in pdat.items()}
    # transient pressure
    w["transient"] = {}
    for el, obj in d.get("Transient_Pressure", {}).items():
        series = sorted(((parse_date(k), float(v)) for k, v in obj["Data"].items()), key=lambda t: t[0])
        w["transient"][float(el)] = series
    # flow records
    w["flow"] = None
    for key in ("Production_Mass_Flow", "Injection_Mass_Flow"):
        if key in d:
            series = sorted(((parse_date(k), float(v)) for k, v in d[key]["Data"].items()), key=lambda t: t[0])
            w["flow"] = {"type": "production" if key.startswith("Prod") else "injection",
                         "start": series[0][0], "end": series[-1][0], "first": series[0][1],
                         "last": series[-1][1], "n": len(series)}
    if "Transient_Enthalpy" in d:
        ser = sorted(((parse_date(k), float(v)) for k, v in d["Transient_Enthalpy"]["Data"].items()),
                     key=lambda t: t[0])
        w["enth"] = {"start": ser[0][0], "first": ser[0][1], "n": len(ser)}
    if "Injection_Enthalpy" in d:
        ser = sorted(((parse_date(k), float(v)) for k, v in d["Injection_Enthalpy"]["Data"].items()),
                     key=lambda t: t[0])
        w["inj_enth"] = {"first": ser[0][1]}
    return w


def load_wells_geometry_json():
    wj = json.loads(WELLS_JSON.read_text(encoding="utf-8"))
    return wj


def block_parts(block):
    """'arz38' -> ('arz', 38). Column = first three letters, layer = trailing digits."""
    col, lay = block[:3], block[3:]
    return col, int(lay)


def block_for(wj_well, z):
    """Model block for elevation z: Interp_Weights (max weight) first, then the well's own block intervals."""
    iw = {float(k): v for k, v in wj_well["Interp_Weights"].items()}
    for k, v in iw.items():
        if abs(k - z) < 1e-6:
            blocks, weights = v["blocks"], v["weights"]
            return blocks[int(np.argmax(weights))], "interp_weights"
    for blk in wj_well["Blocks"]:
        lo, hi = sorted(blk["interval"])
        if lo - 1e-6 <= z <= hi + 1e-6:
            return blk["block"], "track_interval"
    return None, "none"


# ----------------------------------------------------------------------------- integrity checks
def integrity_checks(wells, wj_all, geom_wells, verts, layers, surfa):
    rows = []
    log("Data integrity (CLAUDE 3.1.1):")
    for w in wells:
        n = w["num"]
        name = w["name"]
        wj = wj_all[name]
        # dat cross-check for T and P
        for qty, raw, fname in [
            ("temperature", w["T_raw"], f"Well{n}_downhole_temperature.dat"),
            ("pressure", w["P_raw"], f"Well{n}_downhole_pressure.dat"),
        ]:
            f = DATA_DAT / fname
            if not f.exists():
                f = DATA_DAT / fname.replace(f"Well{n}_", f"Wel{n}_")  # WELL-10 files use 'Wel10_'
            if not f.exists():
                rows.append([name, qty, len(raw), "", "", "", "dat file missing"])
                log(f"  {name} {qty}: dat file {fname} not found")
                continue
            arr = np.loadtxt(f, ndmin=2)
            vals_dat, z_dat = arr[:, 0], arr[:, 1]
            matched, maxdiff = 0, 0.0
            for v, z in zip(vals_dat, z_dat):
                if z in raw:
                    matched += 1
                    maxdiff = max(maxdiff, abs(raw[z] - v))
            status = "OK" if (matched == len(raw) == len(vals_dat) and maxdiff < 1e-6) else "MISMATCH"
            rows.append([name, qty, len(raw), len(vals_dat), matched, f"{maxdiff:.6g}", status])
            log(f"  {name} {qty}: JSON {len(raw)} pts, dat {len(vals_dat)} rows, matched {matched}, "
                f"max |diff| {maxdiff:.3g} -> {status}")
        # transient pressure cross-check (dat index = year - 2014)
        for el, series in w["transient"].items():
            fname = [p.name for p in DATA_DAT.glob("*transient_pressure.dat")
                     if (p.name.startswith(f"Well{n}_") or p.name.startswith(f"Wel{n}_")) and f"_{el:.1f}_" in p.name]
            if not fname:
                rows.append([name, f"transient {el:.0f} mRL", len(series), "", "", "", "dat file missing"])
                continue
            arr = np.loadtxt(DATA_DAT / fname[0], ndmin=2)
            json_by_year = {d.year: v for d, v in series}
            matched, maxdiff = 0, 0.0
            for idx, v in arr:
                yr = 2014 + int(idx)
                if yr in json_by_year:
                    matched += 1
                    maxdiff = max(maxdiff, abs(json_by_year[yr] - v))
            status = "OK" if matched == len(arr) == len(series) and maxdiff < 1e-6 else "MISMATCH"
            rows.append([name, f"transient {el:.0f} mRL", len(series), len(arr), matched,
                         f"{maxdiff:.6g}", status])
            log(f"  {name} transient {el:.0f} mRL: JSON {len(series)}, dat {len(arr)}, matched {matched} -> {status}")
        # survey points inside the track and in the design / feed-zone intervals
        track_z = [t[2] for t in w["track"]]
        zmin, zmax = min(track_z), max(track_z)
        out_track = [z for z in w["T_raw"] if not (zmin - 1e-6 <= z <= zmax + 1e-6)]
        rows.append([name, "survey points inside track z-range", len(w["T_raw"]), "", "",
                     "", "OK" if not out_track else f"outside: {out_track}"])
        for fz in w["feedzones"]:
            lo, hi = sorted([fz["top"], fz["bot"]])
            in_slot = any(dsg["slotted"] and dsg["bot"] - 1e-6 <= lo and hi <= dsg["top"] + 1e-6
                          for dsg in w["design"])
            in_track = (zmin - 1e-6 <= lo) and (hi <= zmax + 1e-6)
            rows.append([name, f"feed zone {lo:.0f}..{hi:.0f} ({fz['quality']})", "", "", "", "",
                         "in track" if in_track else "NOT in track"
                         + ("; in slotted liner" if in_slot else "; NOT in slotted liner")])
        # feed-zone intervals in the well-to-block json match the JSON feed zones
        for fz in w["feedzones"]:
            found = False
            for fzk, blocks in wj["Feed_Zones"].items():
                for blk, info in blocks.items():
                    a, b = info["Interval"]
                    if abs(a - fz["top"]) < 1e-6 and abs(b - fz["bot"]) < 1e-6:
                        found = True
            if not found:
                rows.append([name, "feed zone not in gTE14515_wells.json", "", "", "", "", "CHECK"])
                log(f"  {name}: feed zone {fz['top']}..{fz['bot']} not matched in gTE14515_wells.json")
        # geometry WELLS section vs JSON track (top and bottom)
        aut2 = w["aut2"]
        if aut2 in geom_wells:
            gp = geom_wells[aut2]
            jt, jb = w["track"][0], w["track"][-1]
            gt, gb = gp[0], gp[-1]
            d_top = math.dist(jt[:2], gt[:2]) + abs(jt[2] - gt[2])
            d_bot = math.dist(jb[:2], gb[:2]) + abs(jb[2] - gb[2])
            st = "OK" if d_top < 0.5 and d_bot < 0.5 else "CHECK"
            rows.append([name, "track ends: JSON vs gTE14515.dat WELLS", "", "", "", f"{max(d_top, d_bot):.3g}", st])
            log(f"  {name} track ends vs WELLS section: max diff {max(d_top, d_bot):.3g} m -> {st}")
        # layer consistency of the block assigned to each survey point
        n_in, n_out = 0, 0
        for z in w["T_raw"]:
            blk, _ = block_for(wj, z)
            if blk is None:
                n_out += 1
                continue
            _, lay = block_parts(blk)
            if lay in layers and layers[lay][0] - 1e-6 <= z <= layers[lay][1] + 1e-6:
                n_in += 1
            else:
                n_out += 1
        rows.append([name, "survey points whose block layer contains z", len(w["T_raw"]), "", n_in, "",
                     "OK" if n_out == 0 else f"{n_out} point(s) fall in layer gaps or outside layers"])
    return rows


# ----------------------------------------------------------------------------- analysis
def build_points(w, wj, layers, surfa):
    """One record per survey point (T and P together)."""
    recs = []
    track = w["track"]
    zs = sorted(w["T_raw"].keys(), reverse=True)
    top_z = track[0][2]
    for z in zs:
        T = w["T_raw"][z]
        P = w["P_raw"].get(z, float("nan"))
        x, y = track_xy(track, z)
        ts = tsat_c(P)
        margin = ts - T if np.isfinite(ts) else float("nan")
        blk, blk_src = block_for(wj, z)
        col, lay = (block_parts(blk) if blk else ("", -1))
        in_fz = any(fz["bot"] - 1e-9 <= z <= fz["top"] + 1e-9 for fz in w["feedzones"]) \
            or any(fz["top"] - 1e-9 <= z <= fz["bot"] + 1e-9 for fz in w["feedzones"])
        slotted = any(dsg["slotted"] and dsg["bot"] - 1e-9 <= z <= dsg["top"] + 1e-9 for dsg in w["design"])
        casing = any((not dsg["slotted"]) and dsg["bot"] - 1e-9 <= z <= dsg["top"] + 1e-9 for dsg in w["design"])
        surf = surfa.get(col, float("nan"))
        above = bool(np.isfinite(surf) and z > surf)
        near_wh = z > top_z - NEAR_WH_M
        recs.append(dict(well=w["name"], num=w["num"], x=x, y=y, z=z, T=T, P=P, Tsat=ts, margin=margin,
                         block=blk or "", block_src=blk_src, column=col, layer=lay, in_fz=in_fz,
                         open_hole=slotted, casing=casing and not slotted, above=above, near_wh=near_wh,
                         surf=surf))
    return recs


def classify_point_weight(r, conf):
    """Weight scheme (CLAUDE 4, default table; reason recorded)."""
    if r["above"]:
        base, why = 0.0, "above model surface (SURFA)"
    elif r["near_wh"]:
        base, why = 0.3, f"within {NEAR_WH_M:.0f} m of top of track (wellhead / top layer)"
    elif r["in_fz"]:
        base, why = 1.0, "feed zone"
    elif r["open_hole"]:
        base, why = 0.8, "open hole (slotted liner)"
    elif r["casing"]:
        base, why = 0.5, "behind casing"
    else:
        base, why = 0.5, "no design interval matched (default)"
    weight = base * CONF_FACTOR[conf]
    return weight, f"{why}; x confidence {conf} ({CONF_FACTOR[conf]:.1f})"


def crossings(zs, vs, thr):
    """Elevations where the profile crosses threshold thr (linear interpolation), top to bottom."""
    out = []
    for i in range(len(zs) - 1):
        a, b = vs[i] - thr, vs[i + 1] - thr
        if a == 0:
            out.append(zs[i])
        elif a * b < 0:
            f = a / (a - b)
            out.append(zs[i] + f * (zs[i + 1] - zs[i]))
    return out


def classify_regime(zs, Ts, margins):
    """Rule-based regime label. zs and Ts ordered top to bottom."""
    Tmax = max(Ts)
    k = int(np.argmax(Ts))
    n_near_boil = int(np.sum(np.array([m for m in margins if np.isfinite(m)]) < MARGIN_FLAG))
    if Tmax < 100:
        return "cold (<100 C): recharge/downflow or injection-cooled"
    if n_near_boil >= 2 and Tmax >= 200:
        return "two-phase / near-boiling upflow (>=2 points with margin < 5 K)"
    if 0 < k < len(Ts) - 1 and Ts[-1] <= Tmax - 10:
        return "upflow with lateral outflow (interior T maximum; cooler below)"
    if 0 < k < len(Ts) - 1:
        return "interior T maximum with weak cooling below (<10 K): mixed"
    if Ts[-1] >= Ts[0]:
        return "conductive (T increases with depth, no inversion)"
    return "mixed / downflow"


def well_summary(w, pts, wj, surfa):
    pts = sorted(pts, key=lambda r: -r["z"])  # top to bottom
    zs = np.array([r["z"] for r in pts])
    Ts = np.array([r["T"] for r in pts])
    Ps = np.array([r["P"] for r in pts])
    margins = [r["margin"] for r in pts]
    kmax = int(np.argmax(Ts))
    z_maxT, Tmax = zs[kmax], Ts[kmax]
    zTD, T_TD = zs[-1], Ts[-1]
    # inversion: interior maximum, then the lowest temperature below it
    inv_top, inv_base, inv_drop = "", "", 0.0
    if 0 < kmax < len(Ts) - 1:
        below = Ts[kmax + 1:]
        j = kmax + 1 + int(np.argmin(below))
        inv_top, inv_base, inv_drop = zs[kmax], zs[j], Tmax - Ts[j]
    # isothermal sections (|dT| < 1 K over >= 50 m), reported as top:base
    iso = []
    i = 0
    while i < len(zs) - 1:
        j = i
        while j + 1 < len(zs) and abs(Ts[j + 1] - Ts[i]) <= ISO_TOL:
            j += 1
        if j > i and (zs[i] - zs[j]) >= ISO_MIN_LEN:
            iso.append(f"{zs[i]:.0f}:{zs[j]:.0f}")
        i = max(j, i + 1)
    # gradients in K/km (positive = warmer with depth)
    def grad(mask):
        if mask.sum() < 2:
            return float("nan")
        slope = np.polyfit(zs[mask], Ts[mask], 1)[0]  # dT/dz
        return -slope * 1000.0
    shallow = grad((zs >= 0) & (zs <= zs[0]))
    deep = grad(zs <= zTD + 500.0)
    t100 = crossings(list(zs), list(Ts), 100.0)
    t200 = crossings(list(zs), list(Ts), 200.0)
    finite_m = [(m, z) for m, z in zip(margins, zs) if np.isfinite(m)]
    min_m, z_min_m = (min(finite_m)[0], min(finite_m)[1]) if finite_m else (float("nan"), float("nan"))
    n_flag = sum(1 for m, _ in finite_m if m < MARGIN_FLAG)
    # estimated water level from the top survey point (liquid column to atmosphere)
    top = pts[0]
    wl = float("nan")
    wl_note = ""
    if np.isfinite(top["margin"]) and top["margin"] > MARGIN_FLAG:
        rho_top = rho_liq(top["P"], top["T"])
        wl = top["z"] + (top["P"] - ATM) * 1e5 / (rho_top * G)
        if wl > top["z"] + 5:
            wl_note = "water level above top survey point: top pressure exceeds 1 atm (check wellhead datum)"
    else:
        wl_note = "not estimated: top point is at or near saturation"
    feed_P = []
    for fz in w["feedzones"]:
        mid = 0.5 * (fz["top"] + fz["bot"])
        feed_P.append(interp_no_extrap(zs, Ps, mid))
    feed_T = [interp_no_extrap(zs, Ts, 0.5 * (fz["top"] + fz["bot"])) for fz in w["feedzones"]]
    surf_top = surfa.get(pts[0]["column"], float("nan"))
    regime = classify_regime(list(zs), list(Ts), margins)
    # natural-state verdict and confidence (CLAUDE 3.1)
    verdict, conf, vreason = natural_state_verdict(w)
    return {
        "well": w["name"], "n_survey_points": len(pts), "survey_date": w["survey_date"],
        "track_top_mRL": w["track"][0][2], "track_bottom_mRL": w["track"][-1][2],
        "max_T_degC": Tmax, "max_T_elev_mRL": z_maxT, "T_at_TD_degC": T_TD, "TD_mRL": zTD,
        "shallow_grad_K_per_km": shallow, "deep_grad_K_per_km": deep,
        "inversion_top_mRL": inv_top, "inversion_base_mRL": inv_base, "inversion_drop_K": inv_drop,
        "isothermal_sections_mRL": ";".join(iso),
        "T100_crossings_mRL": ";".join(f"{v:.0f}" for v in t100),
        "T200_crossings_mRL": ";".join(f"{v:.0f}" for v in t200),
        "min_boiling_margin_K": min_m, "min_margin_elev_mRL": z_min_m, "n_points_margin_below_5K": n_flag,
        "est_water_level_mRL": wl, "water_level_note": wl_note,
        "feedzone_P_bar": ";".join(f"{v:.1f}" for v in feed_P),
        "feedzone_T_degC": ";".join(f"{v:.1f}" for v in feed_T),
        "regime": regime, "ns_verdict": verdict, "ns_confidence": conf, "verdict_reason": vreason,
        "wellhead_surface_SURFA_mRL": surf_top,
    }


def natural_state_verdict(w):
    """Rule-based verdict from the timing records (CLAUDE 3.1.2). Returns (verdict, confidence, reason)."""
    survey = parse_date("29/07/2015") if w["survey_date"] == "29/07/2015" else parse_date(w["survey_date"])
    fl = w["flow"]
    if fl and fl["type"] == "injection":
        return ("possibly disturbed", "low",
                f"injection {fl['first']:.0f} kg/s from {fl['start']:%Y} (injectate enthalpy "
                f"{w['inj_enth']['first']:.1f} kJ/kg); cold-water injection since record start")
    if fl and fl["type"] == "production" and fl["start"] <= survey:
        return ("possibly disturbed", "low",
                f"production record starts {fl['start']:%d/%m/%Y} (before survey {survey:%d/%m/%Y}); "
                f"survey is not pre-production")
    if fl and fl["type"] == "production" and fl["start"] > survey:
        return ("probably natural state", "low",
                f"survey ({survey:%d/%m/%Y}) predates production record ({fl['start']:%d/%m/%Y}); "
                f"transient disagreement unresolved")
    return ("probably natural state", "low",
            "no production, injection or transient record; timing cannot be confirmed (all surveys carry the same date)")


def hydrostatic_tests(w, pts):
    """Column-density tests (CLAUDE 3.3): intervals >= 200 m, single-phase at every point."""
    pts = sorted(pts, key=lambda r: -r["z"])
    rows = []
    for i in range(len(pts)):
        j = None
        for k in range(i + 1, len(pts)):
            if pts[i]["z"] - pts[k]["z"] >= HYDRO_MIN_LEN:
                j = k
                break
        if j is None:
            continue
        seg = pts[i:j + 1]
        single = all((not np.isfinite(r["margin"]) and r["P"] >= P_CRIT) or (r["margin"] >= MARGIN_FLAG)
                     for r in seg)
        a, b = pts[i], pts[j]
        dz = a["z"] - b["z"]          # positive: b is deeper
        dP = b["P"] - a["P"]          # positive: pressure rises with depth
        Tm = 0.5 * (a["T"] + b["T"])
        Pm = 0.5 * (a["P"] + b["P"])
        rho_if = rho_liq(Pm, Tm) if single else float("nan")
        rho_obs = dP * 1e5 / (G * dz)
        excess = (rho_obs / rho_if - 1.0) * 100.0 if single else float("nan")
        resid = dP - rho_if * G * dz / 1e5 if single else float("nan")
        rows.append(dict(well=w["name"], z_top=a["z"], z_bot=b["z"], dz=dz, dP=dP, Tm=Tm, Pm=Pm,
                         rho_obs=rho_obs, rho_if=rho_if, excess=excess, resid=resid, single=single))
    return rows


def pressure_consistency(w, pts):
    """Survey pressure at transient elevation vs first transient value (CLAUDE 3.1.2)."""
    out = []
    if not w["transient"]:
        return out
    zs = [r["z"] for r in pts]
    Ps = [r["P"] for r in pts]
    survey = parse_date(w["survey_date"])
    sy = decimal_year(survey)
    for el, series in w["transient"].items():
        p_survey = interp_no_extrap(zs, Ps, el)
        first_d, first_p = series[0]
        d_first_p = p_survey - first_p
        # transient value at survey date (time interpolation); NaN outside the record
        t_series = [decimal_year(d) for d, _ in series]
        v_series = [v for _, v in series]
        if sy < t_series[0] or sy > t_series[-1]:
            p_time, d_time = float("nan"), float("nan")
            note = "survey date outside transient record"
        else:
            p_time = float(np.interp(sy, t_series, v_series))
            d_time = p_survey - p_time
            note = ""
        if abs(d_first_p) <= 1.0:
            verdict = "agrees with first transient (<=1 bar)"
        elif abs(d_first_p) <= 3.0:
            verdict = "marginal agreement with first transient (1-3 bar)"
        else:
            verdict = "disagrees with first transient (>3 bar)"
        if np.isfinite(d_time):
            if abs(d_time) > 3.0:
                verdict += "; disagrees with transient at survey date (>3 bar)"
            else:
                verdict += "; consistent with transient at survey date"
        out.append(dict(well=w["name"], transient_mRL=el, survey_P_at_elev_bar=p_survey,
                        first_transient_date=f"{first_d:%d/%m/%Y}", first_transient_bar=first_p,
                        diff_first_bar=d_first_p, transient_at_survey_bar=p_time, diff_at_survey_bar=d_time,
                        note=note, verdict=verdict))
    return out


def dist_fault_river(w, mapping, track):
    x, y = track[0][0], track[0][1]
    res = {}
    for nm in ["Fault 1", "Fault 2", "Fault 3"]:
        res[nm] = min(dist_to_polyline(x, y, pl) for pl in mapping[nm])
    res["Tarawera River"] = min(dist_to_polyline(x, y, pl) for pl in mapping["Tarawera River"])
    res["river"] = min(dist_to_polyline(x, y, pl) for pl in mapping["river"])
    lake = mapping["Lake Tarawera"][0]
    closed = np.vstack([lake, lake[:1]])
    res["Lake shore"] = dist_to_polyline(x, y, closed)
    res["Inside lake"] = bool(MPath(lake).contains_point((x, y)))
    return res


# ----------------------------------------------------------------------------- figures
def style():
    plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.3, "figure.dpi": 110})


def fig_map(wells, mapping, verts):
    fig, ax = plt.subplots(figsize=(8, 7.5))
    vx = [v[1] for v in verts]
    vy = [v[2] for v in verts]
    ax.scatter(vx, vy, s=1.5, c="0.6", marker=".", label="model column vertices (VERTICES)")
    lake = mapping["Lake Tarawera"][0]
    ax.fill(lake[:, 0], lake[:, 1], color="#9ecae1", alpha=0.6, label="Lake Tarawera (outline)")
    for pl in mapping["Tarawera River"]:
        ax.plot(pl[:, 0], pl[:, 1], color="#3182bd", lw=1.2)
    for pl in mapping["river"]:
        ax.plot(pl[:, 0], pl[:, 1], color="#3182bd", lw=1.2, ls="--")
    for nm in ["Fault 1", "Fault 2", "Fault 3"]:
        pl = mapping[nm][0]
        ax.plot(pl[:, 0], pl[:, 1], color="#cb181d", lw=2, ls="-.")
        ax.text(pl[0, 0], pl[0, 1], nm, color="#cb181d", fontsize=8)
    for w in wells:
        xs = [t[0] for t in w["track"]]
        ys = [t[1] for t in w["track"]]
        ax.plot(xs, ys, color="k", lw=1.5)
        ax.scatter([xs[0]], [ys[0]], c="k", s=18)
        ax.text(xs[0] + 40, ys[0] + 40, w["name"], fontsize=8, weight="bold")
    ax.set_xlabel("x (NZTM, m)")
    ax.set_ylabel("y (NZTM, m)")
    ax.set_aspect("equal")
    ax.set_title("Wells (black), faults (red), water (blue), grid vertices (grey)")
    ax.legend(loc="lower left", fontsize=7)
    fig.tight_layout()
    fig.savefig(FIG / "well_map_grid_faults_water.png")
    plt.close(fig)


def fig_profile(w, pts, tsat_curve):
    pts = sorted(pts, key=lambda r: -r["z"])
    zs = np.array([r["z"] for r in pts])
    Ts = np.array([r["T"] for r in pts])
    Ps = np.array([r["P"] for r in pts])
    tss = np.array([r["Tsat"] for r in pts])
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.5, 6), sharey=True)
    for d in w["design"]:
        color = "#bdbdbd" if not d["slotted"] else "#c7e9c0"
        label = "casing" if not d["slotted"] else "slotted liner (open)"
        for ax in (a1, a2):
            ax.axhspan(d["bot"], d["top"], color=color, alpha=0.35, lw=0)
    for fz in w["feedzones"]:
        for ax in (a1, a2):
            ax.axhspan(fz["bot"], fz["top"], color="#fb6a4a", alpha=0.6, lw=0)
    a1.plot(Ts, zs, "o-", color="#d62728", ms=3, label="T survey")
    m = np.isfinite(tss)
    a1.plot(tss[m], zs[m], "--", color="#1f77b4", lw=1, label="Tsat(P) (IF97)")
    a1.set_xlabel("T (deg C)")
    a1.set_ylabel("Elevation (mRL)")
    a1.set_title(f"{w['name']}: temperature")
    a1.legend(fontsize=7, loc="lower right")
    a2.plot(Ps, zs, "o-", color="#1f77b4", ms=3, label="P survey")
    if w["transient"]:
        for el, series in w["transient"].items():
            a2.scatter([series[0][1]], [el], marker="s", color="k", s=28, label="first transient")
            a2.scatter([series[-1][1]], [el], marker="^", color="#2ca02c", s=28, label="last transient")
    a2.set_xlabel("P (bar)")
    a2.set_title("pressure (survey, transient)")
    a2.legend(fontsize=7, loc="lower right")
    fig.suptitle("Grey: casing; green: slotted liner; red bands: feed zones", fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / f"profile_{w['name']}.png")
    plt.close(fig)


def fig_map_by_elev(field, what, cmap, unit):
    """field: {well: (track, zs, vals)}. Each well is placed at its own track position at the elevation."""
    fig, axes = plt.subplots(2, 4, figsize=(12, 6.5))
    axes = axes.ravel()
    for ax, e in zip(axes, MAP_ELEVS):
        xs, ys, vs, names = [], [], [], []
        for nm, (track, zs_, vals_) in field.items():
            v = interp_no_extrap(zs_, vals_, e)
            x, y = track_xy(track, e)
            if np.isfinite(v) and np.isfinite(x):
                xs.append(x)
                ys.append(y)
                vs.append(v)
                names.append(nm)
        ax.set_title(f"{what} at {e:.0f} mRL", fontsize=9)
        if len(xs) >= 3:
            tri = Triangulation(np.array(xs), np.array(ys))
            interp = LinearTriInterpolator(tri, np.array(vs))
            xg = np.linspace(min(xs), max(xs), 120)
            yg = np.linspace(min(ys), max(ys), 120)
            X, Y = np.meshgrid(xg, yg)
            Z = interp(X, Y)
            cf = ax.contourf(X, Y, Z, levels=12, cmap=cmap)
            fig.colorbar(cf, ax=ax, label=unit, shrink=0.8)
        ax.scatter(xs, ys, c=vs, cmap=cmap, edgecolors="k", s=40, zorder=3)
        for nm, x, y, v in zip(names, xs, ys, vs):
            ax.annotate(f"{nm.replace('WELL-', 'W')}\n{v:.0f}", (x, y), fontsize=6, xytext=(3, 3),
                        textcoords="offset points")
        ax.set_aspect("equal")
        ax.tick_params(labelsize=6)
    for ax in axes[len(MAP_ELEVS):]:
        ax.axis("off")
    fig.suptitle(f"{what}: interpolated only inside the convex hull of wells that have data at that elevation "
                 "(no extrapolation)", fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / f"maps_{what.split()[0].lower()}_by_elevation.png")
    plt.close(fig)


def fig_section(wells_pts, wells):
    """wells_pts: {name: (x, y)} at wellhead; wells: {name: [records]}."""
    x0, y0 = wells_pts[f"WELL-{SECTION_WELLS[0]}"]
    x1, y1 = wells_pts[f"WELL-{SECTION_WELLS[1]}"]
    u = np.array([x1 - x0, y1 - y0], float)
    L = float(np.linalg.norm(u))
    u /= L
    rows = []
    for nm, pts in wells.items():
        for r in pts:
            d = (r["x"] - x0) * u[0] + (r["y"] - y0) * u[1]
            off = -(r["x"] - x0) * u[1] + (r["y"] - y0) * u[0]
            rows.append((nm, d, off, r["z"], r["T"], r["P"]))
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5.5), sharey=True)
    inside = [r for r in rows if abs(r[2]) <= SECTION_CORRIDOR_M]
    if len(inside) >= 4:
        xs = np.array([r[1] for r in inside])
        zs_ = np.array([r[3] for r in inside])
        Ts = np.array([r[4] for r in inside])
        Ps = np.array([r[5] for r in inside])
        tri = Triangulation(xs, zs_)
        cf = a1.tricontourf(tri, Ts, levels=14, cmap="inferno")
        fig.colorbar(cf, ax=a1, label="T (deg C)", shrink=0.8)
        ps_ok = np.isfinite(Ps)
        if ps_ok.sum() >= 4:
            tri2 = Triangulation(xs[ps_ok], zs_[ps_ok])
            cf2 = a2.tricontourf(tri2, Ps[ps_ok], levels=14, cmap="viridis")
            fig.colorbar(cf2, ax=a2, label="P (bar)", shrink=0.8)
    for nm, d, off, z, T, P in rows:
        if abs(off) <= SECTION_CORRIDOR_M:
            a1.scatter([d], [z], c=[T], cmap="inferno", s=14, edgecolors="k", vmin=0, vmax=290, zorder=3)
        else:
            a1.scatter([d], [z], c="none", s=10, edgecolors="0.6", zorder=2)  # outside corridor: outline only
    for nm in wells:
        dd = [r[1] for r in rows if r[0] == nm and abs(r[2]) <= SECTION_CORRIDOR_M]
        oo = [r[2] for r in rows if r[0] == nm and abs(r[2]) <= SECTION_CORRIDOR_M]
        if dd:
            a1.text(np.mean(dd), 230, f"{nm}\n({np.mean(oo):.0f} m off-line)", fontsize=6, ha="center")
    if inside:
        xin = [r[1] for r in inside]
        a1.set_xlim(min(xin) - 150, max(xin) + 150)
        a2.set_xlim(min(xin) - 150, max(xin) + 150)
    a1.set_ylim(-3000, 350)
    a1.set_xlabel(f"distance along section from WELL-{SECTION_WELLS[0]} to WELL-{SECTION_WELLS[1]} (m)")
    a1.set_ylabel("Elevation (mRL)")
    a1.set_title("T: survey points within 400 m of the line (black edge)")
    a2.set_xlabel("distance along section (m)")
    a2.set_title("P: survey points within 400 m of the line")
    fig.tight_layout()
    fig.savefig(FIG / "section_WELL-10_to_WELL-3.png")
    plt.close(fig)
    return u, L, rows


def fig_PT(wells):
    fig, ax = plt.subplots(figsize=(7, 6))
    tt = np.linspace(100.0, 373.9, 300)
    pp = np.array([psat_bar(t) for t in tt])
    ax.plot(tt, pp, "k-", lw=1.5, label="saturation curve (IF97)")
    ax.axhline(P_CRIT, color="0.5", ls=":", lw=1)
    ax.text(105, P_CRIT + 3, "critical pressure 220.64 bar", fontsize=7, color="0.4")
    cmap = plt.get_cmap("tab10")
    for i, (nm, pts) in enumerate(wells.items()):
        ax.plot([r["T"] for r in pts], [r["P"] for r in pts], "o-", ms=3, lw=0.8,
                color=cmap(i % 10), label=nm)
    ax.set_xlabel("T (deg C)")
    ax.set_ylabel("P (bar)")
    ax.set_xlim(0, 380)
    ax.set_ylim(0, 300)
    ax.set_title("P-T state: points left/above the curve are liquid; points on it are at saturation")
    ax.legend(fontsize=6, ncol=2, loc="upper left")
    fig.tight_layout()
    fig.savefig(FIG / "PT_state_vs_saturation.png")
    plt.close(fig)


# ----------------------------------------------------------------------------- main
def main():
    OUT.mkdir(exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    style()
    log("=== natural_state_analysis.py run ===")
    log(f"Run date: {datetime.now():%Y-%m-%d %H:%M}")
    log("Reading inputs from Baseline/ only (read-only).")
    validate_thermo()

    verts, layers, surfa, geom_wells = parse_geometry()
    mapping = parse_mapping()
    wj_all = {f"WELL-{n}": None for n in range(1, 11)}
    wjson = load_wells_geometry_json()
    for k in wj_all:
        wj_all[k] = wjson[k]

    wells = []
    for n in range(1, 11):
        w = load_well(n)
        w["num"] = n
        wells.append(w)
        log(f"Loaded {w['name']}: {len(w['T_raw'])} T pts, {len(w['P_raw'])} P pts, "
            f"{len(w['feedzones'])} feed zones, survey {w['survey_date']}, flow: "
            f"{w['flow']['type'] + ' from ' + format(w['flow']['start'], '%d/%m/%Y') if w['flow'] else 'none'}")
    log(f"Survey dates present: {sorted(set(w['survey_date'] for w in wells))}")
    log(f"Design source labels: {sorted(set(w['design_source'] for w in wells))}")

    # ---- integrity
    integ_rows = integrity_checks(wells, wj_all, geom_wells, verts, layers, surfa)
    write_csv(OUT / "data_integrity_check.csv",
              ["well", "check", "n_json", "n_dat", "n_matched", "max_abs_diff", "status"], integ_rows)

    # ---- per-point records
    all_pts = {}
    cal_rows = []
    for w in wells:
        pts = build_points(w, wj_all[w["name"]], layers, surfa)
        all_pts[w["name"]] = pts
        for r in pts:
            conf = "low"
            verdict, conf, _ = natural_state_verdict(w)
            r["conf"] = conf
    # confidence per well (verdict) is attached to each point
    ver = {}
    for w in wells:
        ver[w["name"]] = natural_state_verdict(w)

    # ---- targets CSV
    for w in wells:
        verdict, conf, _ = ver[w["name"]]
        for r in sorted(all_pts[w["name"]], key=lambda q: -q["z"]):
            wgt, why = classify_point_weight(r, conf)
            cal_rows.append([
                r["well"], fmt(r["x"], 1), fmt(r["y"], 1), fmt(r["z"], 1), fmt(r["T"], 2), fmt(r["P"], 2),
                fmt(r["Tsat"], 2), fmt(r["margin"], 2), r["column"], r["layer"], r["block"],
                r["in_fz"], r["casing"], r["above"], conf, fmt(wgt, 3), why,
            ])
    cal_rows = cal_rows
    write_csv(OUT / "NS_calibration_targets.csv",
              ["well", "x", "y", "elevation_mRL", "T_degC", "P_bar", "Tsat_degC", "boiling_margin_degC",
               "column", "layer", "block", "in_feedzone", "behind_casing", "above_model_surface",
               "ns_confidence", "weight", "weight_reason"], cal_rows)

    # ---- feed zones
    fz_rows = []
    for w in wells:
        verdict, conf, _ = ver[w["name"]]
        pts = all_pts[w["name"]]
        zs = [r["z"] for r in pts]
        Ts = [r["T"] for r in pts]
        Ps = [r["P"] for r in pts]
        wj = wj_all[w["name"]]
        for k, fz in enumerate(w["feedzones"], start=1):
            mid = 0.5 * (fz["top"] + fz["bot"])
            T = interp_no_extrap(zs, Ts, mid)
            P = interp_no_extrap(zs, Ps, mid)
            ts = tsat_c(P)
            margin = ts - T if np.isfinite(ts) else float("nan")
            wj_blocks = wj["Feed_Zones"].get(f"Feed {k}", {})
            wj_block = list(wj_blocks.keys())[0] if wj_blocks else ""
            blk_mid, _ = block_for(wj, mid)
            blk_mid = blk_mid or ""
            n_inside = sum(1 for z in zs if min(fz["top"], fz["bot"]) - 1e-9 <= z <= max(fz["top"], fz["bot"]) + 1e-9)
            x, y = track_xy(w["track"], mid)
            wgt = 1.0 * CONF_FACTOR[conf]
            fz_rows.append([w["name"], f"Feed {k}", fmt(fz["top"], 1), fmt(fz["bot"], 1), fz["quality"],
                            fmt(mid, 1), fmt(x, 1), fmt(y, 1), fmt(T, 2), fmt(P, 2), fmt(ts, 2), fmt(margin, 2),
                            n_inside, wj_block, blk_mid, blk_mid[:3] if blk_mid else "",
                            blk_mid[3:] if blk_mid else "", conf, fmt(wgt, 3)])
    write_csv(OUT / "feedzones.csv",
              ["well", "feed_zone", "top_mRL", "bottom_mRL", "quality_tag", "midpoint_mRL", "x", "y",
               "T_degC_at_mid", "P_bar_at_mid", "Tsat_degC", "boiling_margin_degC", "n_survey_pts_inside",
               "block_from_wells_json", "block_at_midpoint", "column", "layer", "ns_confidence", "weight"],
              fz_rows)

    # ---- layer centres
    lc_rows = []
    for w in wells:
        pts = all_pts[w["name"]]
        zs = [r["z"] for r in pts]
        Ts = [r["T"] for r in pts]
        Ps = [r["P"] for r in pts]
        for lay in sorted(layers):
            if lay == 0:
                continue
            lo, hi = layers[lay]
            c = 0.5 * (lo + hi)
            T = interp_no_extrap(zs, Ts, c)
            if not np.isfinite(T):
                continue
            P = interp_no_extrap(zs, Ps, c)
            ts = tsat_c(P)
            lc_rows.append([w["name"], lay, fmt(lo, 1), fmt(hi, 1), fmt(c, 1), fmt(T, 2), fmt(P, 2),
                            fmt(ts, 2), fmt(ts - T if np.isfinite(ts) else float("nan"), 2)])
    write_csv(OUT / "T_P_at_layer_centres.csv",
              ["well", "layer", "layer_bottom_mRL", "layer_top_mRL", "layer_centre_mRL", "T_degC", "P_bar",
               "Tsat_degC", "boiling_margin_degC"], lc_rows)

    # ---- T and P at common elevations (interpolation only, blank outside survey range)
    ce_rows = []
    for w in wells:
        pts = all_pts[w["name"]]
        zs = [r["z"] for r in pts]
        Ts = [r["T"] for r in pts]
        Ps = [r["P"] for r in pts]
        row = [w["name"]]
        for e in MAP_ELEVS:
            row.append(fmt(interp_no_extrap(zs, Ts, e), 1))
            row.append(fmt(interp_no_extrap(zs, Ps, e), 1))
        ce_rows.append(row)
    ce_hdr = ["well"]
    for e in MAP_ELEVS:
        ce_hdr += [f"T_at_{e:.0f}_mRL_degC", f"P_at_{e:.0f}_mRL_bar"]
    write_csv(OUT / "T_P_at_common_elevations.csv", ce_hdr, ce_rows)

    # ---- well summary
    summ = []
    for w in wells:
        summ.append(well_summary(w, all_pts[w["name"]], wj_all[w["name"]], surfa))
    sum_cols = list(summ[0].keys())
    write_csv(OUT / "well_summary.csv", sum_cols,
              [[fmt(s[c], 2) if isinstance(s[c], (float, np.floating)) else s[c] for c in sum_cols] for s in summ])

    # ---- pressure consistency
    pc_rows = []
    for w in wells:
        for r in pressure_consistency(w, all_pts[w["name"]]):
            pc_rows.append([r["well"], fmt(r["transient_mRL"], 1), fmt(r["survey_P_at_elev_bar"], 2),
                            r["first_transient_date"], fmt(r["first_transient_bar"], 2),
                            fmt(r["diff_first_bar"], 2), fmt(r["transient_at_survey_bar"], 2),
                            fmt(r["diff_at_survey_bar"], 2), r["note"], r["verdict"]])
    write_csv(OUT / "pressure_consistency_check.csv",
              ["well", "transient_elev_mRL", "survey_P_at_elev_bar", "first_transient_date",
               "first_transient_bar", "diff_survey_minus_first_bar", "transient_at_survey_date_bar",
               "diff_survey_minus_transient_at_survey_bar", "note", "verdict"], pc_rows)
    for r in pc_rows:
        log(f"  pressure check {r[0]} @ {r[1]} mRL: survey {r[2]} bar vs first {r[3]} {r[4]} bar "
            f"(diff {r[5]}); at survey date diff {r[7]} -> {r[9]}")

    # ---- hydrostatic
    hy_rows = []
    all_excess = []
    per_well_ex = {}
    for w in wells:
        hs = hydrostatic_tests(w, all_pts[w["name"]])
        for h in hs:
            hy_rows.append([h["well"], fmt(h["z_top"], 1), fmt(h["z_bot"], 1), fmt(h["dz"], 1),
                            fmt(h["dP"], 2), fmt(h["Tm"], 1), fmt(h["Pm"], 1), fmt(h["rho_obs"], 1),
                            fmt(h["rho_if"], 1), fmt(h["excess"], 2), fmt(h["resid"], 2),
                            "single-phase" if h["single"] else "excluded (two-phase or near saturation)"])
            if h["single"] and np.isfinite(h["excess"]):
                all_excess.append(h["excess"])
                per_well_ex.setdefault(h["well"], []).append(h["excess"])
    field_med = float(np.median(all_excess)) if all_excess else float("nan")
    for nm, vals in per_well_ex.items():
        log(f"  hydrostatic {nm}: median excess {np.median(vals):+.2f} % over {len(vals)} interval(s)")
    log(f"  FIELD-WIDE median excess (all valid intervals): {field_med:+.2f} %, "
        f"range {min(all_excess):+.2f} to {max(all_excess):+.2f} %")
    hy_rows.append(["FIELD", "", "", "", "", "", "", "", "", fmt(field_med, 2), "",
                    f"median of {len(all_excess)} valid intervals"])
    write_csv(OUT / "hydrostatic_check.csv",
              ["well", "z_top_mRL", "z_bottom_mRL", "dz_m", "dP_bar", "T_mean_degC", "P_mean_bar",
               "rho_implied_kg_m3", "rho_IF97_kg_m3", "excess_pct", "residual_bar", "status"], hy_rows)

    # ---- geometry
    geo_rows = []
    for w in wells:
        g = dist_fault_river(w, mapping, w["track"])
        top_block = wj_all[w["name"]]["Blocks"][0]["block"]
        col = top_block[:3]
        surf = surfa.get(col, float("nan"))
        wh = w["track"][0][2]
        dev_h = math.dist(w["track"][0][:2], w["track"][-1][:2])
        nearest_fault = min(g["Fault 1"], g["Fault 2"], g["Fault 3"])
        geo_rows.append([w["name"], fmt(w["track"][0][0], 1), fmt(w["track"][0][1], 1), fmt(wh, 1),
                         fmt(w["track"][-1][2], 1), col, fmt(surf, 1), fmt(wh - surf, 1), fmt(dev_h, 1),
                         fmt(g["Fault 1"], 0), fmt(g["Fault 2"], 0), fmt(g["Fault 3"], 0),
                         fmt(nearest_fault, 0), fmt(g["Tarawera River"], 0), fmt(g["river"], 0),
                         fmt(g["Lake shore"], 0), str(g["Inside lake"])])
    write_csv(OUT / "well_geometry_check.csv",
              ["well", "wellhead_x", "wellhead_y", "wellhead_mRL", "track_bottom_mRL", "column_top",
               "model_surface_SURFA_mRL", "wellhead_minus_surface_m", "horizontal_offset_top_to_bottom_m",
               "dist_fault1_m", "dist_fault2_m", "dist_fault3_m", "nearest_fault_m",
               "dist_Tarawera_River_m", "dist_river_m", "dist_lake_shore_m", "inside_lake_outline"], geo_rows)

    # ---- figures
    log("Figures:")
    fig_map(wells, mapping, verts)
    log("  well map")
    for w in wells:
        fig_profile(w, all_pts[w["name"]], None)
    log("  10 well profiles")
    # maps by elevation: T and P for wells with data
    field_T, field_P = {}, {}
    for w in wells:
        pts = all_pts[w["name"]]
        zs = [r["z"] for r in pts]
        field_T[w["name"]] = (w["track"], zs, [r["T"] for r in pts])
        field_P[w["name"]] = (w["track"], zs, [r["P"] for r in pts])
    fig_map_by_elev(field_T, "Temperature T", "inferno", "deg C")
    fig_map_by_elev(field_P, "Pressure P", "viridis", "bar")
    log("  T and P maps by elevation")
    sec_heads = {w["name"]: (w["track"][0][0], w["track"][0][1]) for w in wells}
    sec_u, sec_L, sec_rows = fig_section(sec_heads, {w["name"]: all_pts[w["name"]] for w in wells})
    log(f"  section WELL-10 to WELL-3, length {sec_L:.0f} m")
    fig_PT({w["name"]: sorted(all_pts[w["name"]], key=lambda r: -r["z"]) for w in wells})
    log("  P-T state plot")

    # ---- console log
    log("Done.")
    (OUT / "run_console_log.txt").write_text("\n".join(LOG) + "\n", encoding="utf-8")


if __name__ == "__main__":
    sys.exit(main())
