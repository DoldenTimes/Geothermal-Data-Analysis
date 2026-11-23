"""
Natural-state data analysis for the Teaching_Model_3D geothermal dataset.

Purpose: interpret the downhole temperature / pressure data and well geometry
so that an engineer (or agent) can see what the natural-state calibration
targets are, and where the data are internally consistent or inconsistent.

Inputs  (read-only):
    Baseline/Teaching_Model_3D/Data/Data_json/WELL-n_data.json   (geometry + downhole P/T)
    Baseline/Teaching_Model_3D/Data/Data_dat/*.dat               (cross-check of same P/T)
    Baseline/Teaching_Model_3D/Data/Mapping/Mapping.json         (lake / faults for maps)

The model files under Baseline/.../Models are NOT used as data (per CLAUDE.md).

Outputs (written to Output/):
    well_geometry_check.csv
    downhole_profiles_long.csv
    hydrostatic_check.csv
    pressure_datum_reduction.csv
    well_pressure_consistency.csv
    feedzones_summary.csv
    T_P_at_common_elevations.csv
    dat_vs_json_check.csv
    NS_summary_stats.json
    figures/*.png

Run:  python Scripts/natural_state_analysis.py      (from the repository root)
"""

import glob
import json
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------------
ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "Baseline", "Teaching_Model_3D", "Data")
JSON_DIR = os.path.join(DATA, "Data_json")
DAT_DIR = os.path.join(DATA, "Data_dat")
MAPPING = os.path.join(DATA, "Mapping", "Mapping.json")
OUT = os.path.join(ROOT, "Output")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

G = 9.80665          # m/s2
COMMON_Z = [0, -250, -500, -750, -1000, -1250, -1500]   # masl, for comparison tables
LOG = []


def log(msg):
    print(msg)
    LOG.append(msg)


# ----------------------------------------------------------------------------
# Fluid properties
# ----------------------------------------------------------------------------
def rho_water(T_C):
    """Liquid water density (kg/m3), Kell (1975) polynomial. Valid ~0-150 C at
    saturation; used here for T < 150 C, which covers all measured data."""
    T = np.asarray(T_C, dtype=float)
    num = (999.83952 + 16.945176 * T - 7.9870401e-3 * T**2
           - 46.170461e-6 * T**3 + 105.56302e-9 * T**4 - 280.54253e-12 * T**5)
    return num / (1.0 + 16.879850e-3 * T)


def t_sat_C(p_bar):
    """Saturation temperature (C) of water at pressure p_bar, IAPWS-IF97 Region 4
    backward-form (Wagner-type) equation. Validated against tabulated values below."""
    p = np.asarray(p_bar, dtype=float) / 10.0          # MPa
    n1, n2, n3, n4, n5 = 0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5, -0.32325550322333e7
    n6, n7, n8, n9, n10 = 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849, 0.65017534844798e3
    beta = p ** 0.25
    E = beta**2 + n3 * beta + n6
    F = n1 * beta**2 + n4 * beta + n7
    Gc = n2 * beta**2 + n5 * beta + n8
    D = 2.0 * Gc / (-F - np.sqrt(F**2 - 4.0 * E * Gc))
    Tk = 0.5 * (n10 + D - np.sqrt((n10 + D)**2 - 4.0 * (n9 + n10 * D)))
    return Tk - 273.15


def validate_tsat():
    tests = {1.01325: 99.974, 1.0: 99.606, 10.0: 179.88, 100.0: 311.00}   # bar -> C (IAPWS-IF97 tables, approx.)
    rows = []
    for p, Tref in tests.items():
        Tc = float(t_sat_C(p))
        rows.append((p, Tref, Tc, Tc - Tref))
        log(f"  Tsat check p={p} bar: model={Tc:.2f} C, reference={Tref:.2f} C, diff={Tc - Tref:+.2f} K")
    worst = max(abs(r[3]) for r in rows)
    if worst > 0.5:
        raise RuntimeError(f"Tsat validation failed (max |diff| {worst:.2f} K)")
    return rows


# ----------------------------------------------------------------------------
# Hydrostatic integration helpers
# ----------------------------------------------------------------------------
def column_integrate(z_from, z_to, P_from_bar, z_prof, T_prof, n_steps=None):
    """Integrate hydrostatic pressure from z_from to z_to (masl) using the
    well's measured temperature profile (linear interp, clamped at ends).
    Returns pressure (bar) at z_to. Pressure increases downward."""
    z_prof = np.asarray(z_prof, float)
    T_prof = np.asarray(T_prof, float)
    order = np.argsort(z_prof)
    zp, Tp = z_prof[order], T_prof[order]
    if abs(z_to - z_from) < 1e-9:
        return P_from_bar
    if n_steps is None:
        n_steps = max(2, int(abs(z_to - z_from)) + 2)
    zs = np.linspace(z_from, z_to, n_steps)
    T = np.interp(zs, zp, Tp)
    integrand = rho_water(T) * G                        # Pa/m
    dP_Pa = np.trapezoid(integrand, zs)                 # signed
    return P_from_bar - dP_Pa / 1e5


# ----------------------------------------------------------------------------
# Load geometry and downhole data
# ----------------------------------------------------------------------------
def load_wells():
    wells = {}
    for f in glob.glob(os.path.join(JSON_DIR, "WELL-*_data.json")):
        d = json.load(open(f))
        name = d["Well_Name"]
        n = int(re.search(r"WELL-(\d+)", name).group(1))
        dh_T = d["Downhole_Temperature"]
        dh_P = d["Downhole_Pressure"]
        # Only one survey date exists per well; check and take it.
        assert len(dh_T) == 1 and len(dh_P) == 1, f"{name}: >1 survey date"
        date_T = list(dh_T)[0]
        date_P = list(dh_P)[0]
        T_rec = dh_T[date_T]
        P_rec = dh_P[date_P]
        z_T = np.array([float(k) for k in T_rec["Data"]])
        T_val = np.array(list(T_rec["Data"].values()), float)
        z_P = np.array([float(k) for k in P_rec["Data"]])
        P_val = np.array(list(P_rec["Data"].values()), float)

        # Track: list of "[x, y, z]" strings
        track = np.array([json.loads(s.replace("'", '"')) for s in d["Well_Track"]["Data"]], float)
        # Feed zones: "[top, bottom, 'quality']"
        fz = []
        for s in d["Feed_Zones"]["Data"]:
            parts = re.findall(r"\[\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)\s*,\s*'([^']*)'\s*\]", s)
            top, bot, q = parts[0]
            fz.append((float(top), float(bot), q))
        design = []
        for s in d["Well_Design"]["Data"]:
            m = re.match(r"\[\s*(-?[\d.]+)\s*,\s*(-?[\d.]+)\s*,\s*'([^']*)'\s*,\s*(True|False)\s*,\s*'([^']*)'\s*\]", s)
            design.append({"top": float(m.group(1)), "bottom": float(m.group(2)),
                           "diam": m.group(3), "slotted": m.group(4) == "True", "spec": m.group(5)})

        flow_keys = [k for k in ["Production_Mass_Flow", "Injection_Mass_Flow", "Transient_Enthalpy",
                                 "Injection_Enthalpy", "Transient_Pressure"] if k in d]
        flow_info = {}
        for k in flow_keys:
            dd = d[k].get("Data", {})
            flow_info[k] = (len(dd), list(dd.keys())[0] if dd else "")

        wells[name] = dict(
            n=n, name=name, AUT2=d["AUT2_Name"], date=date_T, date_P=date_P,
            z_T=z_T, T=T_val, z_P=z_P, P=P_val, track=track, feed=fz, design=design,
            flow_info=flow_info, json_raw=d,
            stable=(T_rec.get("StableFlag"), P_rec.get("StableFlag")),
        )
    return dict(sorted(wells.items(), key=lambda kv: kv[1]["n"]))


def well_xy_at_z(track, z):
    """Interpolate x,y along the well track at elevation z (masl). Track points are
    ordered top->bottom, z decreasing."""
    zt = track[:, 2]
    x, y = track[:, 0], track[:, 1]
    order = np.argsort(zt)
    zs, xs, ys = zt[order], x[order], y[order]
    if zs[0] == zs[-1]:
        return x[0], y[0]
    return float(np.interp(z, zs, xs)), float(np.interp(z, zs, ys))


# ----------------------------------------------------------------------------
# Main analysis
# ----------------------------------------------------------------------------
def main():
    import csv
    log("=== Natural-state analysis ===")
    log(f"Root: {ROOT}")
    log("Step 1: validate saturation-temperature and density functions")
    validate_tsat()
    log(f"  rho_water(20C)={float(rho_water(20)):.1f}, rho_water(100C)={float(rho_water(100)):.1f}, "
        f"rho_water(120C)={float(rho_water(120)):.1f} kg/m3")

    log("Step 2: load geometry and downhole P/T from Data_json")
    wells = load_wells()
    for w in wells.values():
        log(f"  {w['name']}: T n={len(w['z_T'])}, P n={len(w['z_P'])}, survey {w['date']}, "
            f"stable flags {w['stable']}, feed zones={len(w['feed'])}, flow={list(w['flow_info'].keys())}")

    # ---------------- Step 3: DAT vs JSON cross-check ----------------
    log("Step 3: cross-check Data_dat files against Data_json")
    dat_rows = []
    for w in wells.values():
        prefix = "Wel10" if w["n"] == 10 else f"Well{w['n']}"
        for kind, zkey, vkey in [("downhole_temperature", "z_T", "T"), ("downhole_pressure", "z_P", "P")]:
            fn = os.path.join(DAT_DIR, f"{prefix}_{kind}.dat")
            if not os.path.exists(fn):
                dat_rows.append([w["name"], kind, "missing file", "", "", ""])
                continue
            arr = np.loadtxt(fn)            # columns: value, elevation
            dat_z, dat_v = arr[:, 1], arr[:, 0]
            jz, jv = w[zkey], w[vkey]
            maxdiff = 0.0
            nmatch = 0
            for zz, vv in zip(dat_z, dat_v):
                idx = np.where(np.isclose(jz, zz, atol=1e-6))[0]
                if len(idx):
                    nmatch += 1
                    maxdiff = max(maxdiff, abs(jv[idx[0]] - vv))
            dat_rows.append([w["name"], kind, len(dat_z), nmatch, round(maxdiff, 6),
                             "OK" if (nmatch == len(dat_z) == len(jz) and maxdiff < 1e-6) else "CHECK"])
            log(f"  {w['name']} {kind}: dat n={len(dat_z)}, json n={len(jz)}, matched={nmatch}, max |diff|={maxdiff:.2e}")
    with open(os.path.join(OUT, "dat_vs_json_check.csv"), "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(["well", "variable", "n_dat", "n_matched_to_json", "max_abs_diff", "status"])
        cw.writerows(dat_rows)

    # ---------------- Step 4: geometry checks ----------------
    log("Step 4: geometry and feed-zone consistency checks")
    geo_rows = []
    for w in wells.values():
        tr = w["track"]
        z_head = tr[0, 2]
        z_bot = tr[-1, 2]
        liner = w["design"][-1]
        zmin_meas = min(w["z_T"].min(), w["z_P"].min())
        zmax_meas = max(w["z_T"].max(), w["z_P"].max())
        fz_inside = all((b >= z_bot - 1e-6) and (t <= z_head + 1e-6) for t, b, _ in w["feed"])
        fz_in_liner = all((b >= liner["bottom"] - 1e-6) and (t <= liner["top"] + 1e-6) for t, b, _ in w["feed"])
        meas_in_track = (zmin_meas >= z_bot - 1e-6) and (zmax_meas <= z_head + 1e-6)
        geo_rows.append([w["name"], w["AUT2"], len(tr), round(z_head, 1), round(z_bot, 1),
                         round(zmax_meas, 1), round(zmin_meas, 1), len(w["feed"]),
                         "; ".join(f"{t:.0f}/{b:.0f} {q}" for t, b, q in w["feed"]),
                         "yes" if meas_in_track else "NO", "yes" if fz_inside else "NO",
                         "yes" if fz_in_liner else "NO", w["design"][-1]["spec"],
                         "deviated" if len(tr) > 2 else "vertical (2-pt track)"])
        if not fz_inside:
            log(f"  WARNING {w['name']}: feed zone outside track range")
        if not fz_in_liner:
            log(f"  NOTE {w['name']}: feed zone(s) outside slotted liner interval")
    with open(os.path.join(OUT, "well_geometry_check.csv"), "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(["well", "AUT2_name", "n_track_pts", "wellhead_z_masl", "TD_z_masl",
                     "shallowest_meas_masl", "deepest_meas_masl", "n_feedzones",
                     "feedzones_top/bottom_masl_quality", "meas_within_track", "feedzones_in_track",
                     "feedzones_in_liner", "liner_spec", "track_type"])
        cw.writerows(geo_rows)

    # ---------------- Step 5: long table of profiles ----------------
    log("Step 5: build long table of T, P, Tsat margin, feed-zone flags")
    long_rows = []
    # per-well profiles keyed by elevation (merged on T elevations)
    prof = {}
    for w in wells.values():
        # Use the T elevations as the master list; P should share elevations.
        zs_T = w["z_T"]
        zs_P = w["z_P"]
        Tmap = dict(zip(np.round(zs_T, 3), w["T"]))
        Pmap = dict(zip(np.round(zs_P, 3), w["P"]))
        common = sorted(set(Tmap) & set(Pmap), reverse=True)
        z = np.array(common)
        T = np.array([Tmap[k] for k in common])
        P = np.array([Pmap[k] for k in common])
        Ts = t_sat_C(P)
        head = w["track"][0, 2]
        # gradient to next deeper point (C/100m) and feed zone flag
        grad = np.full_like(T, np.nan)
        for i in range(len(z) - 1):
            dz = z[i] - z[i + 1]                      # depth interval (>0)
            grad[i] = (T[i + 1] - T[i]) / dz * 100.0 if dz != 0 else np.nan   # C/100m, positive = warmer with depth
        in_fz, fq = [], []
        for zi in z:
            hit = [q for t, b, q in w["feed"] if b - 1e-6 <= zi <= t + 1e-6]
            in_fz.append("yes" if hit else "no")
            fq.append(hit[0] if hit else "")
        for i in range(len(z)):
            x, y = well_xy_at_z(w["track"], z[i])
            long_rows.append([w["name"], round(z[i], 2), round(head - z[i], 2), round(x, 1), round(y, 1),
                              round(T[i], 2), round(P[i], 3), round(float(Ts[i]), 2),
                              round(float(Ts[i] - T[i]), 2), (round(grad[i], 2) if np.isfinite(grad[i]) else ""),
                              in_fz[i], fq[i]])
        prof[w["name"]] = dict(z=z, T=T, P=P, Tsat=Ts, grad=grad)
    with open(os.path.join(OUT, "downhole_profiles_long.csv"), "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(["well", "z_masl", "depth_below_wellhead_m", "x_m", "y_m", "T_C", "P_bar",
                     "Tsat_at_P_C", "Tsat_minus_T_C", "dT_dz_to_next_C_per_100m", "in_feedzone", "feedzone_quality"])
        cw.writerows(long_rows)

    # ---------------- Step 6: hydrostatic check ----------------
    log("Step 6: hydrostatic consistency (pressure predicted from shallowest point with well's own T profile)")
    hyd_rows = []
    hyd_summary = {}
    for w in wells.values():
        p = prof[w["name"]]
        z, T, P = p["z"], p["T"], p["P"]
        z_anchor, P_anchor = z[0], P[0]
        resid_list = []
        for i in range(1, len(z)):
            P_pred = column_integrate(z_anchor, z[i], P_anchor, z, T)
            r = P[i] - P_pred
            resid_list.append(r)
            # segment check: implied mean density over this interval
            dz = z[i - 1] - z[i]
            dP = (P[i] - P[i - 1]) * 1e5
            rho_imp = dP / (G * dz) if dz > 0 else np.nan
            rho_model = float(rho_water(0.5 * (T[i] + T[i - 1])))
            hyd_rows.append([w["name"], round(z[i], 1), round(P[i], 3), round(P_pred, 3), round(r, 3),
                             round(z[i - 1], 1), round(z[i], 1), round(rho_imp, 1), round(rho_model, 1),
                             round(100 * (rho_imp / rho_model - 1), 2) if np.isfinite(rho_imp) else ""])
        r = np.array(resid_list)
        seg_vals = [row[9] for row in hyd_rows if row[0] == w["name"] and row[9] != ""]
        rel = float(np.nanmedian(seg_vals)) if seg_vals else float("nan")
        # whole-column comparison: implied mean density (measured dP / dz) vs Kell mean density
        P_pred_deep = column_integrate(z_anchor, z[-1], P_anchor, z, T)
        dz_col = z_anchor - z[-1]
        rho_col_meas = (P[-1] - P_anchor) * 1e5 / (G * dz_col)
        rho_col_kell = (P_pred_deep - P_anchor) * 1e5 / (G * dz_col)
        hyd_summary[w["name"]] = dict(resid_mean_bar=float(np.mean(r)), resid_max_abs_bar=float(np.max(np.abs(r))),
                                      resid_at_deepest_bar=float(r[-1]),
                                      seg_density_excess_pct_median=rel,
                                      column_rho_implied_kgm3=float(rho_col_meas),
                                      column_rho_kell_kgm3=float(rho_col_kell),
                                      column_excess_pct=float(100 * (rho_col_meas / rho_col_kell - 1)))
        log(f"  {w['name']}: anchor {z_anchor:.0f} masl P={P_anchor:.2f} bar; residual mean {np.mean(r):+.2f} bar, "
            f"deepest {r[-1]:+.2f} bar, max |res| {np.max(np.abs(r)):.2f} bar; column-mean density implied "
            f"{rho_col_meas:.0f} vs Kell {rho_col_kell:.0f} kg/m3 ({100 * (rho_col_meas / rho_col_kell - 1):+.1f}%); "
            f"median segment excess {rel:+.1f}%")
    with open(os.path.join(OUT, "hydrostatic_check.csv"), "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(["well", "z_masl", "P_meas_bar", "P_pred_from_shallowest_bar", "residual_bar",
                     "z_upper_masl", "z_lower_masl", "rho_implied_kgm3", "rho_Kell_kgm3", "implied_vs_Kell_pct"])
        cw.writerows(hyd_rows)

    # ---------------- Step 7: pressure reduced to datum 0 masl ----------------
    log("Step 7: reduce each pressure to datum 0 masl using the well's own column (cross-well consistency)")
    dat0_rows = []
    per_well_stats = []
    for w in wells.values():
        p = prof[w["name"]]
        z, T, P = p["z"], p["T"], p["P"]
        P0 = []
        for i in range(len(z)):
            P0_i = column_integrate(z[i], 0.0, P[i], z, T)
            P0.append(P0_i)
            dat0_rows.append([w["name"], round(z[i], 1), round(P[i], 3), round(P0_i, 3)])
        P0 = np.array(P0)
        # trend of datum pressure with depth: slope bar per 100 m
        slope = np.polyfit(z, P0, 1)[0] * 100 if len(z) > 1 else np.nan
        per_well_stats.append([w["name"], round(float(np.mean(P0)), 3), round(float(np.std(P0)), 3),
                               round(float(np.min(P0)), 3), round(float(np.max(P0)), 3), round(float(slope), 3), len(z)])
        log(f"  {w['name']}: P@0masl mean {np.mean(P0):.2f} bar, sd {np.std(P0):.2f}, range {np.min(P0):.2f}-{np.max(P0):.2f}, "
            f"trend with z {slope:+.3f} bar/100 m")
    with open(os.path.join(OUT, "pressure_datum_reduction.csv"), "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(["well", "z_masl", "P_meas_bar", "P_at_datum0_bar"])
        cw.writerows(dat0_rows)
    with open(os.path.join(OUT, "well_pressure_consistency.csv"), "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(["well", "P_datum0_mean_bar", "P_datum0_sd_bar", "P_datum0_min_bar", "P_datum0_max_bar",
                     "trend_bar_per_100m_vs_z", "n_points"])
        cw.writerows(per_well_stats)
    datum_means = {r[0]: r[1] for r in per_well_stats}
    spread = max(datum_means.values()) - min(datum_means.values())
    log(f"  Spread of well-mean datum pressures (0 masl): {spread:.2f} bar")

    # ---------------- Step 8: feed zones ----------------
    log("Step 8: feed-zone summary")
    fz_rows = []
    for w in wells.values():
        p = prof[w["name"]]
        z, T, P = p["z"], p["T"], p["P"]
        for top, bot, q in w["feed"]:
            # interpolate T and P at feed top and bottom (within measured range only)
            def at(zq, arr):
                zs = z[::-1]; a = arr[::-1]
                if zq < zs[0] or zq > zs[-1]:
                    return np.nan
                return float(np.interp(zq, zs, a))
            Tt, Tb = at(top, T), at(bot, T)
            Pt, Pb = at(top, P), at(bot, P)
            # local gradient of T across the zone
            gT = (Tb - Tt) / (top - bot) * 100 if np.isfinite(Tt) and np.isfinite(Tb) and top != bot else np.nan   # C/100m, positive = warmer with depth
            fz_rows.append([w["name"], top, bot, q, round(Tt, 2) if np.isfinite(Tt) else "", round(Tb, 2) if np.isfinite(Tb) else "",
                            round(Pt, 2) if np.isfinite(Pt) else "", round(Pb, 2) if np.isfinite(Pb) else "",
                            round(gT, 2) if np.isfinite(gT) else "", round(float(t_sat_C(Pt)) - Tt, 1) if np.isfinite(Pt) and np.isfinite(Tt) else ""])
    with open(os.path.join(OUT, "feedzones_summary.csv"), "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(["well", "top_masl", "bottom_masl", "quality", "T_top_C", "T_bottom_C", "P_top_bar", "P_bottom_bar",
                     "mean_dT_dz_C_per_100m", "Tsat_minus_T_at_top_C"])
        cw.writerows(fz_rows)

    # ---------------- Step 9: common elevations ----------------
    log("Step 9: interpolate T and P at common elevations (no extrapolation)")
    ce_rows = []
    for w in wells.values():
        p = prof[w["name"]]
        z, T, P = p["z"], p["T"], p["P"]
        zs, Ts_, Ps_ = z[::-1], T[::-1], P[::-1]
        row = [w["name"]]
        for zc in COMMON_Z:
            if zs[0] <= zc <= zs[-1]:
                row += [round(float(np.interp(zc, zs, Ts_)), 2), round(float(np.interp(zc, zs, Ps_)), 2)]
            else:
                row += ["", ""]
        ce_rows.append(row)
    hdr = ["well"]
    for zc in COMMON_Z:
        hdr += [f"T_{zc}masl_C", f"P_{zc}masl_bar"]
    with open(os.path.join(OUT, "T_P_at_common_elevations.csv"), "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(hdr)
        cw.writerows(ce_rows)

    # ---------------- Step 10: summary JSON ----------------
    summary = {}
    for w in wells.values():
        p = prof[w["name"]]
        summary[w["name"]] = dict(
            n_points=int(len(p["z"])), z_top=float(p["z"].max()), z_bottom=float(p["z"].min()),
            T_min=float(p["T"].min()), T_max=float(p["T"].max()),
            P_min=float(p["P"].min()), P_max=float(p["P"].max()),
            min_Tsat_minus_T=float(np.min(p["Tsat"] - p["T"])),
            mean_grad_upper_C_per_100m=float(np.nanmean(p["grad"][: max(1, len(p["grad"]) // 3)])),
            mean_grad_lower_C_per_100m=float(np.nanmean(p["grad"][-max(1, len(p["grad"]) // 3):])),
            hydrostatic=hyd_summary[w["name"]],
        )
    summary["_datum_spread_bar"] = float(spread)
    json.dump(summary, open(os.path.join(OUT, "NS_summary_stats.json"), "w"), indent=2)

    # ---------------- Step 11: figures ----------------
    log("Step 11: figures")
    colors = plt.cm.tab10(np.linspace(0, 1, 10))

    # Fig 1: well map
    fig, ax = plt.subplots(figsize=(7, 6))
    try:
        mp = json.load(open(MAPPING))
        for k, v in mp.items():
            if not (isinstance(v, dict) and isinstance(v.get("Data"), list) and v["Data"]):
                continue
            data = v["Data"]
            # a single polyline is a list of [x, y]; a multi-segment item is a list of such lists
            segs = [data] if isinstance(data[0][0], (int, float)) else data
            for seg in segs:
                arr = np.array(seg, float)
                if arr.ndim == 2 and arr.shape[1] >= 2:
                    ax.plot(arr[:, 0], arr[:, 1], color="0.6", lw=1, label=k)
    except Exception as e:
        log(f"  mapping overlay skipped: {e}")
    for i, (name, w) in enumerate(wells.items()):
        tr = w["track"]
        ax.plot(tr[:, 0], tr[:, 1], "-", color=colors[i % 10], lw=2)
        ax.plot(tr[0, 0], tr[0, 1], "o", color=colors[i % 10], ms=6)
        ax.annotate(name.replace("WELL-", "W"), (tr[0, 0], tr[0, 1]), fontsize=8, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_title("Well locations (wellhead dots; deviated tracks shown as lines)")
    ax.ticklabel_format(useOffset=False, style="plain")
    ax.tick_params(axis="x", rotation=30)
    ax.set_aspect("equal", adjustable="datalim")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "01_well_map.png"), dpi=130)
    plt.close(fig)

    # Fig 2: T and P vs elevation overlay, with feed zones
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 7), sharey=True)
    for i, (name, w) in enumerate(wells.items()):
        p = prof[name]
        a1.plot(p["T"], p["z"], "o-", ms=3, color=colors[i % 10], label=name.replace("WELL-", "W"))
        a2.plot(p["P"], p["z"], "o-", ms=3, color=colors[i % 10])
        for top, bot, q in w["feed"]:
            a1.axhspan(bot, top, color=colors[i % 10], alpha=0.12)
    a1.set_xlabel("Temperature (degC)")
    a2.set_xlabel("Pressure (bar)")
    a1.set_ylabel("Elevation (masl)")
    a1.set_title("Downhole temperature (2015-07-29)")
    a2.set_title("Downhole pressure (2015-07-29)")
    a1.grid(alpha=0.3); a2.grid(alpha=0.3)
    a1.legend(fontsize=7, ncol=2)
    fig.suptitle("Natural-state downhole profiles; shaded bands = feed zones")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "02_T_P_profiles.png"), dpi=130)
    plt.close(fig)

    # Fig 3: P-T state with saturation curve
    fig, ax = plt.subplots(figsize=(7, 6))
    pp = np.linspace(1, 260, 300)
    ax.plot(pp, t_sat_C(pp), "k--", lw=1.2, label="Boiling curve (Tsat)")
    for i, (name, w) in enumerate(wells.items()):
        p = prof[name]
        ax.plot(p["P"], p["T"], "o-", ms=3, color=colors[i % 10], label=name.replace("WELL-", "W"))
    ax.set_xlabel("Pressure (bar)")
    ax.set_ylabel("Temperature (degC)")
    ax.set_title("P-T state of measured points vs boiling curve")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "03_PT_state.png"), dpi=130)
    plt.close(fig)

    # Fig 4: pressure reduced to datum 0 masl, plotted against depth per well
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 6))
    for i, (name, w) in enumerate(wells.items()):
        rows = [r for r in dat0_rows if r[0] == name]
        z = np.array([r[1] for r in rows]); P0 = np.array([r[3] for r in rows])
        a1.plot(P0, z, "o-", ms=3, color=colors[i % 10], label=name.replace("WELL-", "W"))
    a1.set_xlabel("Pressure reduced to datum 0 masl (bar)")
    a1.set_ylabel("Measurement elevation (masl)")
    a1.set_title("Datum-reduced pressure (each well)")
    a1.grid(alpha=0.3); a1.legend(fontsize=7, ncol=2)
    names = list(datum_means.keys())
    means = [datum_means[n] for n in names]
    a2.bar([n.replace("WELL-", "W") for n in names], means, color="0.6")
    a2.set_ylabel("Mean datum pressure (bar)")
    a2.set_title("Well-mean datum pressure")
    a2.tick_params(axis="x", rotation=45)
    a2.set_ylim(min(means) - 2, max(means) + 2)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "04_datum_pressure_consistency.png"), dpi=130)
    plt.close(fig)

    # Fig 5: T and P at common elevations, plotted at well positions (scatter maps)
    fig, axes = plt.subplots(2, 4, figsize=(14, 7))
    show_z = [0, -500, -1000, -1500]
    for j, zc in enumerate(show_z):
        for row_i, var in enumerate(["T", "P"]):
            ax = axes[row_i, j]
            xs, ys, vs, labs = [], [], [], []
            for name, w in wells.items():
                p = prof[name]
                zs, arr = p["z"][::-1], (p["T"] if var == "T" else p["P"])[::-1]
                if zs[0] <= zc <= zs[-1]:
                    x, y = well_xy_at_z(w["track"], zc)
                    xs.append(x); ys.append(y); vs.append(float(np.interp(zc, zs, arr))); labs.append(name.replace("WELL-", "W"))
            if xs:
                sc = ax.scatter(xs, ys, c=vs, cmap="coolwarm" if var == "T" else "viridis", s=120, edgecolor="k")
                for x, y, l, v in zip(xs, ys, labs, vs):
                    ax.annotate(f"{l}\n{v:.0f}", (x, y), fontsize=6, xytext=(3, 3), textcoords="offset points")
                plt.colorbar(sc, ax=ax, shrink=0.8, label=("degC" if var == "T" else "bar"))
            ax.set_title(f"{var} at {zc} masl", fontsize=9)
            ax.ticklabel_format(useOffset=False, style="plain")
            ax.tick_params(labelsize=6)
            ax.set_aspect("equal", adjustable="datalim")
    fig.suptitle("Temperature (top) and pressure (bottom) at common elevations; wells with data only (no extrapolation)")
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "05_T_P_maps_by_elevation.png"), dpi=130)
    plt.close(fig)

    log("Done. Outputs written to " + OUT)
    with open(os.path.join(OUT, "run_console_log.txt"), "w") as f:
        f.write("\n".join(LOG) + "\n")


if __name__ == "__main__":
    main()
