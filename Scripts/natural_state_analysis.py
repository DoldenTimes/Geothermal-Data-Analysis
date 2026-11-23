"""
Natural-state data analysis for the TE14515 teaching geothermal model.

Reads ONLY from Baseline/ (never modified). Writes tables to Output/ and figures to Output/figures/.
Run:  python Scripts/natural_state_analysis.py

Data used (factual):  well downhole T/P surveys, feed zones, well tracks, well-to-block mapping,
                      geometry (layers, column vertices, surface elevations), boundary-condition
                      file, GIS mapping (lake / rivers / faults), initial production-test values.
Not used as fact:     the supplied natural-state model (TE14515_NS_000.*), reference only.
"""
import json, ast, os, re, glob, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
BASE = os.path.join(ROOT, "Baseline", "Teaching_Model_3D")
DATA = os.path.join(BASE, "Data")
MODEL = os.path.join(BASE, "Models", "TE14515")
OUT = os.path.join(ROOT, "Output")
FIG = os.path.join(OUT, "figures")
os.makedirs(FIG, exist_ok=True)

G = 9.80665
BAR = 1e5

# ----------------------------------------------------------------------------
# Steam-property helpers (IAPWS-IF97 region 4 saturation; simple liquid density)
# ----------------------------------------------------------------------------
_n = [0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5,
      -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6,
      -0.23855557567849, 0.65017534844798e3]


def psat_bar(T):
    """Saturation pressure (bar) at T degC (IF97 region 4)."""
    T = np.asarray(T, float) + 273.15
    th = T + _n[8] / (T - _n[9])
    A = th**2 + _n[0] * th + _n[1]
    B = _n[2] * th**2 + _n[3] * th + _n[4]
    C = _n[5] * th**2 + _n[6] * th + _n[7]
    p = (2 * C / (-B + np.sqrt(B * B - 4 * A * C)))**4  # MPa
    return p * 10.0


def tsat_c(P_bar):
    """Saturation temperature (degC) at P bar (IF97 region 4 inverse)."""
    P = np.asarray(P_bar, float) / 10.0  # MPa
    b = P**0.25
    E = b**2 + _n[2] * b + _n[5]
    F = _n[0] * b**2 + _n[3] * b + _n[6]
    Gc = _n[1] * b**2 + _n[4] * b + _n[7]
    D = 2 * Gc / (-F - np.sqrt(F * F - 4 * E * Gc))
    T = (_n[9] + D - np.sqrt((_n[9] + D)**2 - 4 * (_n[8] + _n[9] * D))) / 2
    return T - 273.15


# Saturated liquid enthalpy (kJ/kg), steam-table values
_hf_T = np.array([0, 25, 50, 75, 100, 125, 150, 175, 200, 220, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 340, 350])
_hf_h = np.array([0.0, 104.8, 209.3, 313.9, 419.1, 525.0, 632.2, 741.1, 852.4, 943.6, 1037.6, 1085.8, 1134.9, 1185.2,
                  1236.8, 1290.0, 1344.8, 1401.8, 1461.3, 1524.9, 1594.6, 1670.9])
_hg_h = np.array([2501, 2547, 2592, 2635, 2676, 2714, 2746, 2773, 2793, 2802, 2803, 2801, 2796, 2790, 2780, 2766, 2749,
                  2727, 2700, 2666, 2622, 2564], float)


def hf(T):
    return np.interp(T, _hf_T, _hf_h)


def T_from_hf(h):
    return np.interp(h, _hf_h, _hf_T)


def rho_liq(T, P_bar):
    """Compressed-liquid density (kg/m3): saturated-liquid correlation + constant compressibility."""
    T = np.asarray(T, float)
    # fit to steam-table saturated liquid densities (0-350C)
    Tt = np.array([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 350])
    rt = np.array([999.8, 998.2, 992.2, 983.2, 971.8, 958.4, 943.1, 926.1, 907.4, 886.9, 864.7, 840.3, 813.6, 783.7,
                   750.3, 712.1, 667.1, 610.7, 574.7])
    r0 = np.interp(T, Tt, rt)
    beta = (0.45e-9 + 2.5e-12 * np.clip(T, 0, 340)**1.6 / 100.0)  # 1/Pa, rough compressibility rise with T
    beta = np.clip(beta, 0.45e-9, 3e-9)
    return r0 * (1 + beta * np.asarray(P_bar, float) * BAR * 0.5)


# ----------------------------------------------------------------------------
# Load data
# ----------------------------------------------------------------------------
def wnum(name):
    return int(name.split("-")[1])


wells = {}
for f in sorted(glob.glob(os.path.join(DATA, "Data_json", "WELL-*_data.json")), key=lambda s: wnum(os.path.basename(s).split("_")[0])):
    d = json.load(open(f))
    nm = d["Well_Name"]
    track = np.array([ast.literal_eval(s) for s in d["Well_Track"]["Data"]], float)
    fz = [ast.literal_eval(s) for s in d["Feed_Zones"]["Data"]]
    Tdat = d["Downhole_Temperature"]
    Pdat = d["Downhole_Pressure"]
    dates = list(Tdat.keys())
    assert len(dates) == 1 and list(Pdat.keys()) == dates, "expected exactly one survey per well"
    Tz = np.array(sorted(((float(k), v) for k, v in Tdat[dates[0]]["Data"].items()), reverse=True))
    Pz = np.array(sorted(((float(k), v) for k, v in Pdat[dates[0]]["Data"].items()), reverse=True))
    assert np.allclose(Tz[:, 0], Pz[:, 0]), f"{nm}: T and P survey elevations differ"
    wells[nm] = dict(name=nm, aut2=d["AUT2_Name"], track=track, fz=fz, date=dates[0],
                     z=Tz[:, 0], T=Tz[:, 1], P=Pz[:, 1], raw=d,
                     design=[ast.literal_eval(s) for s in d["Well_Design"]["Data"]])
    wells[nm]["x"], wells[nm]["y"] = track[0, 0], track[0, 1]
    wells[nm]["ztop"], wells[nm]["zbot"] = track[0, 2], track[-1, 2]
    wells[nm]["deviated"] = len(track) > 2

wellnames = sorted(wells, key=wnum)

wj = json.load(open(os.path.join(MODEL, "gTE14515_wells.json")))

# geometry: layers + column centres + surface elevations
lines = open(os.path.join(MODEL, "gTE14515.dat")).read().splitlines()
iv = lines.index("VERTICES"); ic = lines.index("CONNECTIONS"); il = lines.index("LAYERS"); isf = [i for i, l in enumerate(lines) if l.startswith("SURFA")][0]
iw = lines.index("WELLS")
verts = {}
ig = lines.index("GRID")
for l in lines[iv + 1:ig]:
    if len(l) >= 23:
        verts[l[:3]] = (float(l[3:13]), float(l[13:23]))
cols = {}
i = ig + 1
while i < len(lines) and not lines[i].startswith(("CONNECTIONS", "LAYERS")):
    l = lines[i]
    if len(l) > 3 and l[3:].strip():
        nm_c = l[:3]; n = int(l[3:].split()[1]); vs = [lines[i + 1 + k][:3] for k in range(n)]
        xy = np.array([verts[v] for v in vs]); cols[nm_c] = list(xy.mean(axis=0)); i += n + 1
    else:
        i += 1
layers = []
for l in lines[il + 1:isf]:
    p = l.split()
    if len(p) == 3:
        layers.append((int(p[0]), float(p[1]), float(p[2])))
layers = np.array(layers)  # idx, bottom, centre
# layer tops: top of layer i = bottom of layer i-1 (layer 0 top == centre == 1000)
ltop = np.r_[1000.0, layers[:-1, 1]]
surf = {}
for l in lines[isf + 1:iw]:
    p = l.split()
    if len(p) == 2 and len(p[0]) == 3:
        surf[p[0]] = float(p[1])

# boundary conditions
bc = json.load(open(os.path.join(MODEL, "gTE14515_boundary_conditions.json")))
mapping = json.load(open(os.path.join(DATA, "Mapping", "Mapping.json")))
aux = json.load(open(os.path.join(MODEL, "gTE14515.auxjson")))


def layer_of_block(b):
    return int(b[3:])


def col_of_block(b):
    return b[:3]


# ----------------------------------------------------------------------------
# Per-well analysis
# ----------------------------------------------------------------------------
def interp_profile(z, v, zq):
    """Linear interpolation on elevation (z descending)."""
    return np.interp(zq, z[::-1], v[::-1], left=np.nan, right=np.nan)


rows = []
fzrows = []
hyd_rows = []
pcons_rows = []
for nm in wellnames:
    w = wells[nm]
    z, T, P = w["z"], w["T"], w["P"]
    # --- water level (elevation where gauge/abs pressure first exceeds ~1.2 bar): survey P is absolute-ish (1 bar at top)
    # pressure gradient between survey points (bar/100 m)
    dz = -np.diff(z)
    dPdz = np.diff(P) / dz * 100.0  # bar per 100 m (positive with depth)
    dTdz = np.diff(T) / dz * 100.0  # degC per 100 m
    # Boiling point with depth from measured pressure
    Tsat = tsat_c(np.maximum(P, 0.05))
    margin = Tsat - T  # >0 sub-boiling, ~0 boiling, <0 impossible for liquid -> steam
    # hydrostatic from measured T: integrate rho(T,P) g dz starting from the first survey point down
    Ph = np.zeros_like(P); Ph[0] = P[0]
    for i in range(1, len(z)):
        Tm = 0.5 * (T[i] + T[i - 1]); Pm = 0.5 * (P[i] + P[i - 1])
        Ph[i] = Ph[i - 1] + rho_liq(Tm, Pm) * G * (z[i - 1] - z[i]) / BAR
    resid = P - Ph
    # Water level: elevation where P = 1.013 bar (extrapolated from first points)
    # use first survey point: P0 at z0; pressure at which atmospheric; if P0 ~ 1 bar, level ~ z0
    # Level extrapolated upward / downward assuming liquid gradient of near-surface
    wl = z[0] - (1.013 - P[0]) * BAR / (rho_liq(T[0], 1) * G) if P[0] < 1.013 else z[0] + (P[0] - 1.013) * BAR / (rho_liq(T[0], 1) * G)
    # temperature features
    imax = int(np.argmax(T))
    Tmax, zmax = T[imax], z[imax]
    Tbot, zbot = T[-1], z[-1]
    # inversion: temperature drop below the max by >2C deeper than the max
    inv_drop = Tmax - T[imax:].min()
    zinv = z[imax:][int(np.argmin(T[imax:]))]
    # conductive/cap zone: from top down to where dT/dz drops below 5 C/100 m after exceeding it (cap base)
    # steepest gradient
    isteep = int(np.argmax(dTdz))
    steep = dTdz[isteep]
    zsteep = 0.5 * (z[isteep] + z[isteep + 1])
    # top of >=150C, >=200C
    def first_cross(thr):
        idx = np.where(T >= thr)[0]
        if len(idx) == 0:
            return np.nan
        i = idx[0]
        if i == 0:
            return z[0]
        t0, t1 = T[i - 1], T[i]; z0, z1 = z[i - 1], z[i]
        return z0 + (thr - t0) * (z1 - z0) / (t1 - t0)
    z100, z150, z200 = first_cross(100), first_cross(150), first_cross(200)
    # surface conductive gradient (uppermost 100 m of survey)
    top_dT = (T[np.searchsorted(-z, -(z[0] - 100))] - T[0]) if True else np.nan
    # min margin to boiling-point curve (sub-boiling margin), where P>1.5 bar
    mask = P > 1.5
    jm = int(np.argmin(np.where(mask, margin, 1e9)))
    two_phase = np.where(margin < 1.0)[0]
    rows.append(dict(well=nm, x=w["x"], y=w["y"], z_top=w["ztop"], z_bottom=w["zbot"], deviated=w["deviated"],
                     n_pts=len(z), survey_top=z[0], survey_bottom=z[-1], survey_date=w["date"],
                     T_top=T[0], P_top=P[0], water_level_est=round(wl, 1),
                     Tmax=Tmax, z_Tmax=zmax, T_bottom=Tbot, P_bottom=P[-1],
                     inversion_drop_C=round(inv_drop, 1), z_Tmin_below_max=zinv if inv_drop > 2 else np.nan,
                     steepest_dTdz_C_per_100m=round(steep, 1), z_steepest=round(zsteep, 0),
                     z_T100=z100, z_T150=z150, z_T200=z200,
                     min_boiling_margin_C=round(margin[jm], 1), z_min_margin=z[jm],
                     n_pts_within_1C_of_BPD=len(two_phase),
                     hydrostatic_residual_bottom_bar=round(resid[-1], 2),
                     hydrostatic_residual_max_abs_bar=round(np.abs(resid).max(), 2)))
    for zz, t, p, ts, mg, ph, rs in zip(z, T, P, Tsat, margin, Ph, resid):
        hyd_rows.append(dict(well=nm, elevation=zz, T=t, P=p, Tsat_at_P=round(float(ts), 1),
                             boiling_margin=round(float(mg), 1), P_hydrostatic_from_T=round(float(ph), 2),
                             P_resid=round(float(rs), 2)))
    # feed zones: T,P (interpolated over the interval; mean) and block
    mapfz = wj[nm]["Feed_Zones"]
    for k, (ztop, zbot_, q) in enumerate([(a[0], a[1], a[2]) if len(a) > 2 else (a[0], a[1], "") for a in w["fz"]], 1):
        zq = np.linspace(ztop, zbot_, 11)
        Tf = np.nanmean(interp_profile(z, T, zq)); Pf = np.nanmean(interp_profile(z, P, zq))
        blocks = [b for b in mapfz.get(f"Feed {k}", {}).keys()]
        fzrows.append(dict(well=nm, feed=k, top=ztop, bottom=zbot_, quality=q, mid=0.5 * (ztop + zbot_),
                           T_survey=round(float(Tf), 1), P_survey=round(float(Pf), 1),
                           Tsat_at_P=round(float(tsat_c(Pf)), 1),
                           blocks=";".join(blocks),
                           liquid_h_equiv=round(float(hf(Tf)), 0)))

import csv


def write_csv(path, rs):
    with open(path, "w", newline="") as f:
        wr = csv.DictWriter(f, fieldnames=list(rs[0].keys())); wr.writeheader(); wr.writerows(rs)


write_csv(os.path.join(OUT, "well_summary.csv"), rows)
write_csv(os.path.join(OUT, "feedzones.csv"), fzrows)
write_csv(os.path.join(OUT, "hydrostatic_check.csv"), hyd_rows)

# ----------------------------------------------------------------------------
# Production-test (initial) values vs surveys: pressure and enthalpy consistency
# ----------------------------------------------------------------------------
def first_val(dct):
    ks = sorted(dct, key=lambda s: tuple(map(int, s.split("/")))[::-1])
    return ks[0], dct[ks[0]]


for nm in wellnames:
    r = wells[nm]["raw"]
    if "Transient_Pressure" in r:
        for zk, blk in r["Transient_Pressure"].items():
            d0, p0 = first_val(blk["Data"])
            zk = float(zk)
            ps = float(interp_profile(wells[nm]["z"], wells[nm]["P"], np.array([zk]))[0])
            ts = float(interp_profile(wells[nm]["z"], wells[nm]["T"], np.array([zk]))[0])
            eh = r.get("Transient_Enthalpy", {}).get("Data", {})
            dh, h0 = first_val(eh) if eh else ("", np.nan)
            mf = r.get("Production_Mass_Flow", {}).get("Data", {})
            dm, m0 = first_val(mf) if mf else ("", np.nan)
            pcons_rows.append(dict(well=nm, gauge_elev=zk, first_pressure_date=d0, first_transient_P_bar=p0,
                                   survey_P_at_elev_bar=round(ps, 1), diff_bar=round(p0 - ps, 1),
                                   survey_T_at_elev=round(ts, 1), first_enthalpy_date=dh, first_enthalpy_kJkg=h0,
                                   liquid_T_equiv_of_enthalpy=round(float(T_from_hf(h0)), 0) if h0 < 1700 else "steam/2-phase (h>hf max)",
                                   first_massflow_date=dm, first_massflow_kgs=m0,
                                   Tsat_at_gauge_P=round(float(tsat_c(p0)), 1)))
write_csv(os.path.join(OUT, "pressure_consistency_check.csv"), pcons_rows)

# ----------------------------------------------------------------------------
# Calibration targets at model-block (layer-centre) resolution
# ----------------------------------------------------------------------------
tgt = []
tp_layer = []
for nm in wellnames:
    w = wells[nm]
    for b in wj[nm]["Blocks"]:
        lo = layer_of_block(b["block"])
        zc_layer = layers[lo, 2]; zb = layers[lo, 1]; zt = ltop[lo]
        top_i, bot_i = b["interval"]  # well intersection with block
        zmid = 0.5 * (top_i + bot_i)
        # average over the well's interval inside the block (limited to surveyed range)
        zz = np.linspace(max(bot_i, w["z"][-1]), min(top_i, w["z"][0]), 9) if min(top_i, w["z"][0]) > max(bot_i, w["z"][-1]) else np.array([])
        if len(zz) == 0:
            Tm = Pm = np.nan
        else:
            Tm = float(np.nanmean(interp_profile(w["z"], w["T"], zz))); Pm = float(np.nanmean(interp_profile(w["z"], w["P"], zz)))
        # also value at the layer-centre elevation (if surveyed)
        Tc = float(interp_profile(w["z"], w["T"], np.array([zc_layer]))[0]); Pc = float(interp_profile(w["z"], w["P"], np.array([zc_layer]))[0])
        tgt.append(dict(well=nm, block=b["block"], column=col_of_block(b["block"]), layer=lo,
                        layer_top=zt, layer_centre=zc_layer, layer_bottom=zb,
                        well_in_block_top=top_i, well_in_block_bottom=bot_i,
                        T_mean_over_well_interval=round(Tm, 1) if not np.isnan(Tm) else "",
                        P_mean_over_well_interval_bar=round(Pm, 1) if not np.isnan(Pm) else "",
                        T_at_layer_centre=round(Tc, 1) if not np.isnan(Tc) else "",
                        P_at_layer_centre_bar=round(Pc, 2) if not np.isnan(Pc) else "",
                        feed_zone_block=any(b["block"] in str(wj[nm]["Feed_Zones"][k]) for k in wj[nm]["Feed_Zones"])))
write_csv(os.path.join(OUT, "NS_calibration_targets.csv"), tgt)

# T,P at common elevations for inter-well comparison
common = [200, 0, -250, -500, -750, -1000, -1500, -2000, -2500]
tp_rows = []
for nm in wellnames:
    w = wells[nm]
    for zq in common:
        t = float(interp_profile(w["z"], w["T"], np.array([zq]))[0]); p = float(interp_profile(w["z"], w["P"], np.array([zq]))[0])
        if not (np.isnan(t) or np.isnan(p)):
            tp_rows.append(dict(well=nm, x=w["x"], y=w["y"], elevation=zq, T=round(t, 1), P_bar=round(p, 1)))
write_csv(os.path.join(OUT, "T_P_at_common_elevations.csv"), tp_rows)

# ----------------------------------------------------------------------------
# Terrain / geometry diagnostics
# ----------------------------------------------------------------------------
geo_rows = []
cn = np.array([[cols[c][0], cols[c][1]] for c in cols]); cnames = list(cols)
for nm in wellnames:
    w = wells[nm]
    d = np.hypot(cn[:, 0] - w["x"], cn[:, 1] - w["y"]); c = cnames[int(np.argmin(d))]
    def dseg(pt, a, b):
        a = np.array(a, float); b = np.array(b, float); pt = np.array(pt, float)
        t = np.clip(((pt - a) @ (b - a)) / ((b - a) @ (b - a)), 0, 1); return float(np.linalg.norm(pt - (a + t * (b - a))))

    def dist_feature(key_filter):
        best = 1e12
        for k, v in mapping.items():
            if key_filter(k):
                dd = v["Data"]; segs = dd if isinstance(dd[0][0], list) else [dd]
                for sg in segs:
                    for i in range(len(sg) - 1):
                        best = min(best, dseg((w["x"], w["y"]), sg[i], sg[i + 1]))
        return round(best, 0)
    geo_rows.append(dict(well=nm, nearest_column=c, dist_to_col_centre_m=round(float(d.min()), 0),
                         dist_Fault1_m=dist_feature(lambda k: k == "Fault 1"), dist_Fault2_m=dist_feature(lambda k: k == "Fault 2"),
                         dist_Fault3_m=dist_feature(lambda k: k == "Fault 3"),
                         dist_rivers_m=dist_feature(lambda k: "iver" in k), dist_lake_m=dist_feature(lambda k: "Lake" in k),
                         column_surface_elev=surf.get(c, np.nan), well_collar_elev=w["ztop"],
                         collar_minus_surface=round(w["ztop"] - surf.get(c, np.nan), 1)))
write_csv(os.path.join(OUT, "well_geometry_check.csv"), geo_rows)

# ----------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------
cmap = plt.get_cmap("tab10")
col = {nm: cmap(i % 10) for i, nm in enumerate(wellnames)}

# 1. map
fig, ax = plt.subplots(figsize=(9, 8))
sx = [cols[c][0] for c in cnames]; sy = [cols[c][1] for c in cnames]; sz = [surf.get(c, np.nan) for c in cnames]
sc = ax.scatter(sx, sy, c=sz, s=60, cmap="terrain", marker="s", alpha=.8)
plt.colorbar(sc, label="Model surface elevation (mRL)", shrink=.8)
for k, v in mapping.items():
    d = v["Data"]
    segs = d if isinstance(d[0][0], list) else [d]
    c = "tab:blue" if "Lake" in k or "iver" in k else ("red" if "Fault" in k else "gray")
    for j, sg in enumerate(segs):
        p = np.array(sg)
        ax.plot(p[:, 0], p[:, 1], lw=1.4 if "Fault" in k else 1, ls="--" if "Fault" in k else "-", color=c,
                label=k if j == 0 else None)
for nm in wellnames:
    w = wells[nm]
    ax.plot(w["track"][:, 0], w["track"][:, 1], "-", color="k", lw=1)
    ax.plot(w["x"], w["y"], "o", color=col[nm], ms=8, mec="k")
    ax.annotate(nm.replace("WELL-", "W"), (w["x"], w["y"]), xytext=(5, 5), textcoords="offset points", fontweight="bold")
ax.set_xlim(2814000, 2828000); ax.set_ylim(6327000, 6340000)
ax.set_aspect("equal"); ax.set_xlabel("NZTM Easting (m)"); ax.set_ylabel("Northing (m)"); ax.legend(loc="lower right"); ax.set_title("Well locations, model columns, lake, faults")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "01_well_map.png"), dpi=130); plt.close(fig)

# 2. T and P profiles all wells
fig, axs = plt.subplots(1, 3, figsize=(16, 7), sharey=True)
for nm in wellnames:
    w = wells[nm]
    axs[0].plot(w["T"], w["z"], "-o", ms=3, color=col[nm], label=nm)
    axs[1].plot(w["P"], w["z"], "-o", ms=3, color=col[nm], label=nm)
    axs[2].plot(w["P"] - np.interp(w["z"], [-4000, 4000][::1], [0, 0]) , w["z"], alpha=0)  # placeholder (keeps shared y)
# pressure minus hydrostatic reference (cold, 20C column from 1 bar at 300 m)
zr = np.linspace(-3000, 350, 50)
for nm in wellnames:
    w = wells[nm]
    ref = 1.013 + 998.0 * G * (300 - w["z"]) / BAR
    axs[2].plot(w["P"] - ref, w["z"], "-o", ms=3, color=col[nm], label=nm)
axs[0].set_xlabel("Temperature (degC)"); axs[1].set_xlabel("Pressure (bar)"); axs[2].set_xlabel("P - cold hydrostatic ref. (bar)")
axs[0].set_ylabel("Elevation (mRL)")
for a in axs: a.grid(alpha=.3)
axs[0].legend(ncol=2, fontsize=8)
axs[2].set_title("Excess over 20C-water column from 300 mRL")
fig.suptitle("Downhole surveys, 29/07/2015"); fig.tight_layout(); fig.savefig(os.path.join(FIG, "02_T_P_profiles.png"), dpi=130); plt.close(fig)

# 3. per-well panels: T, BPD, feed zones, P vs hydrostatic
fig, axs = plt.subplots(2, 5, figsize=(22, 11), sharey=True)
for a, nm in zip(axs.ravel(), wellnames):
    w = wells[nm]; z, T, P = w["z"], w["T"], w["P"]
    a.plot(T, z, "r-o", ms=3, label="T survey")
    a.plot(tsat_c(np.maximum(P, 0.05)), z, "b--", label="Tsat(P_survey) ~ BPD")
    a.plot(P, z, "g-", label="P survey (bar)")
    for fzr in [r for r in fzrows if r["well"] == nm]:
        a.axhspan(fzr["bottom"], fzr["top"], color="orange", alpha=.4)
    for rr in [r for r in pcons_rows if r["well"] == nm]:
        a.plot(rr["first_transient_P_bar"], rr["gauge_elev"], "k*", ms=12, label="First transient P")
    a.set_title(nm); a.grid(alpha=.3)
axs[0, 0].legend(fontsize=7); [a.set_xlabel("degC / bar") for a in axs[1]]; [a.set_ylabel("Elevation (mRL)") for a in axs[:, 0]]
fig.suptitle("Per-well natural-state profiles (orange = feed zones)"); fig.tight_layout(); fig.savefig(os.path.join(FIG, "03_well_panels.png"), dpi=110); plt.close(fig)

# 4. T maps by elevation (well points) and P maps
for var, fn, lab in (("T", "04_T_maps_by_elevation.png", "Temperature (degC)"), ("P", "05_P_maps_by_elevation.png", "Pressure (bar)")):
    zs = [0, -500, -1000, -1500, -2000, -2500]
    fig, axs = plt.subplots(2, 3, figsize=(17, 11), sharex=True, sharey=True)
    allv = [r["T" if var == "T" else "P_bar"] for r in tp_rows if r["elevation"] in zs]
    vmin, vmax = min(allv), max(allv)
    for a, zq in zip(axs.ravel(), zs):
        pts = [r for r in tp_rows if r["elevation"] == zq]
        a.scatter(sx, sy, c="lightgray", s=25, marker="s")
        s = a.scatter([r["x"] for r in pts], [r["y"] for r in pts], c=[r["T" if var == "T" else "P_bar"] for r in pts], s=260, cmap="inferno" if var == "T" else "viridis", vmin=vmin, vmax=vmax, ec="k")
        for r in pts:
            a.annotate(f'{r["well"].replace("WELL-", "W")}\n{r["T" if var == "T" else "P_bar"]:.0f}', (r["x"], r["y"]), xytext=(9, 4), textcoords="offset points", fontsize=8)
        a.set_title(f"{zq} mRL"); a.set_aspect("equal")
    fig.colorbar(s, ax=axs, label=lab, shrink=.7); fig.savefig(os.path.join(FIG, fn), dpi=120); plt.close(fig)

# 5. cross sections: wells projected on two lines
def section(name, p0, p1, fn, width=1500):
    p0 = np.array(p0, float); p1 = np.array(p1, float); u = (p1 - p0) / np.linalg.norm(p1 - p0); L = np.linalg.norm(p1 - p0)
    fig, axs = plt.subplots(1, 2, figsize=(16, 7), sharey=True)
    used = []
    for nm in wellnames:
        w = wells[nm]; q = np.array([w["x"], w["y"]]) - p0; s = q @ u; off = abs(q[0] * u[1] - q[1] * u[0])
        if -500 <= s <= L + 500 and off <= width:
            used.append((nm, s, off))
            for a, v in zip(axs, ("T", "P")):
                sc = a.scatter(np.full_like(w["z"], s), w["z"], c=w[v], s=30, cmap="inferno" if v == "T" else "viridis", vmin=10 if v == "T" else 0, vmax=300 if v == "T" else 280)
            axs[0].annotate(nm.replace("WELL-", "W"), (s, w["ztop"] + 60), ha="center", fontweight="bold")
    for a, lab in zip(axs, ("T (degC)", "P (bar)")):
        a.set_xlabel("Distance along section (m)"); a.grid(alpha=.3)
    plt.colorbar(axs[0].collections[0], ax=axs[0], label="T (degC)"); plt.colorbar(axs[1].collections[0], ax=axs[1], label="P (bar)")
    axs[0].set_ylabel("Elevation (mRL)"); fig.suptitle(f"{name} (wells within {width} m)"); fig.tight_layout()
    fig.savefig(os.path.join(FIG, fn), dpi=120); plt.close(fig)
    return used


secA = section("Section A: WSW-ENE", (2817500, 6332000), (2823500, 6334500), "06_section_A_WSW-ENE.png", 1800)
secB = section("Section B: SSW-NNE", (2819200, 6332200), (2821500, 6337800), "06_section_B_SSW-NNE.png", 1800)

# 6. isothermal/convective layer & BPD cross-plot for all wells
fig, ax = plt.subplots(1, 2, figsize=(14, 6))
for nm in wellnames:
    w = wells[nm]
    ax[0].plot(w["T"], w["P"], "-o", ms=3, color=col[nm], label=nm)
    ax[1].plot(w["T"], w["z"], "-", color=col[nm], label=nm)
tt = np.linspace(10, 330, 100); ax[0].plot(tt, psat_bar(tt), "k--", label="Saturation")
ax[0].set_xlabel("T (degC)"); ax[0].set_ylabel("P (bar)"); ax[0].legend(fontsize=7, ncol=2); ax[0].grid(alpha=.3); ax[0].set_title("P-T state of surveyed fluid vs saturation curve")
ax[1].set_title("Temperature vs elevation"); ax[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "07_PT_state.png"), dpi=130); plt.close(fig)

# 7. production-test initial values vs survey (pressure)
print("DONE. Wells:", len(wellnames), "targets:", len(tgt))
print("sections", [u[0] for u in secA], [u[0] for u in secB])
