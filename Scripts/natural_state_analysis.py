"""
Natural-state data analysis for the Teaching_Model_3D (TE14515) geothermal field.

Reads ONLY from the Baseline folder (read-only) and writes all results to Output/.
Field data used:  Baseline/Teaching_Model_3D/Data/Data_json/WELL-*_data.json
                  (stable downhole T & P profiles, feed zones, well tracks, casing)
Geometry used:    Baseline/Teaching_Model_3D/Models/TE14515/gTE14515.dat (mulgrid)
                  - grid columns, layers and surface only, to map data onto blocks.
Map features:     Baseline/Teaching_Model_3D/Data/Mapping/Mapping.json
The reference model's rock types / boundary conditions / NS results are NOT used
as data.

Run from the repository root:   python Scripts/natural_state_analysis.py
"""
import ast
import csv
import json
import os
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib.path import Path as MplPath

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "Baseline" / "Teaching_Model_3D"
DATA = BASE / "Data"
GEO = BASE / "Models" / "TE14515" / "gTE14515.dat"
OUT = ROOT / "Output"
FIG = OUT / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

G = 9.81
SURVEY_DATE = "29/07/2015"

# ----------------------------------------------------------------------------
# Thermodynamic helpers
# ----------------------------------------------------------------------------
_N = [0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
      0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
      -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
      0.65017534844798e3]


def tsat_c(p_bar):
    """Saturation temperature (degC) from pressure (bar abs), IAPWS-IF97 region 4."""
    p = np.asarray(p_bar, float) / 10.0  # MPa
    b = p ** 0.25
    E = b ** 2 + _N[2] * b + _N[5]
    F = _N[0] * b ** 2 + _N[3] * b + _N[6]
    Gg = _N[1] * b ** 2 + _N[4] * b + _N[7]
    D = 2 * Gg / (-F - np.sqrt(F ** 2 - 4 * E * Gg))
    T = (_N[9] + D - np.sqrt((_N[9] + D) ** 2 - 4 * (_N[8] + _N[9] * D))) / 2
    return T - 273.15


_RHO_T = np.array([0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240,
                   260, 280, 300, 320])
_RHO = np.array([999.8, 998.2, 992.2, 983.2, 971.8, 958.4, 943.1, 926.1, 907.4,
                 887.0, 864.7, 840.2, 813.6, 783.6, 750.3, 712.1, 667.1])


def rho_liq(T):
    """Approximate liquid water density (kg/m3) at saturation, adequate for
    hydrostatic diagnostics."""
    return np.interp(T, _RHO_T, _RHO)


# ----------------------------------------------------------------------------
# Geometry (mulgrid) parsing
# ----------------------------------------------------------------------------
def read_geometry(path):
    lines = path.read_text().splitlines()
    sec = None
    verts, cols, layers, surf = {}, {}, [], {}
    i = 0
    while i < len(lines):
        ln = lines[i]
        key = ln.strip()
        if key in ("VERTICES", "GRID", "CONNECTIONS", "LAYERS", "SURFA", "WELLS"):
            sec = key
            i += 1
            continue
        if not key:
            i += 1
            continue
        if sec == "VERTICES":
            verts[ln[0:3]] = (float(ln[3:13]), float(ln[13:23]))
        elif sec == "GRID":
            name = ln[0:3]
            nv = int(ln[4:6])
            vn = [lines[i + k + 1][0:3] for k in range(nv)]
            cols[name] = np.array([verts[v] for v in vn])
            i += nv
        elif sec == "LAYERS":
            p = ln.split()
            layers.append((int(p[0]), float(p[1]), float(p[2])))
        elif sec == "SURFA":
            p = ln.split()
            surf[p[0]] = float(p[1])
        i += 1
    # layer tops: top of layer k = bottom of layer k-1
    lay = []
    for k, (idx, bot, cen) in enumerate(layers):
        top = layers[k - 1][1] if k > 0 else bot
        lay.append(dict(index=idx, top=top, bottom=bot, centre=cen))
    return verts, cols, lay[1:], surf  # drop layer 0 (atmosphere/top marker)


VERTS, COLS, LAYERS, SURF = read_geometry(GEO)
COL_PATHS = {c: MplPath(v) for c, v in COLS.items()}
COL_CENTRES = {c: v.mean(axis=0) for c, v in COLS.items()}


def find_column(x, y):
    for c, p in COL_PATHS.items():
        if p.contains_point((x, y)):
            return c
    return None


def find_layer(z):
    for L in LAYERS:
        if L["bottom"] <= z <= L["top"]:
            return L
    return None


def model_block(x, y, z):
    """Return (column, layer index, layer centre, block name, above_surface)."""
    c = find_column(x, y)
    L = find_layer(z)
    if c is None or L is None:
        return c, None, None, None, None
    s = SURF.get(c)
    above = s is not None and z > s
    # PyTOUGH style block name: column name + layer index as 2 chars
    return c, L["index"], L["centre"], f"{c}{L['index']:2d}", above


# ----------------------------------------------------------------------------
# Field data
# ----------------------------------------------------------------------------
def _lit(s):
    return ast.literal_eval(s) if isinstance(s, str) else s


def load_wells():
    wells = {}
    for f in sorted(DATA.glob("Data_json/WELL-*_data.json"),
                    key=lambda p: int(p.stem.split("-")[1].split("_")[0])):
        d = json.loads(f.read_text())
        name = d["Well_Name"]
        track = np.array([_lit(t) for t in d["Well_Track"]["Data"]], float)
        fz = [_lit(t) for t in d["Feed_Zones"]["Data"]]
        design = [_lit(t) for t in d["Well_Design"]["Data"]]
        T = d["Downhole_Temperature"][SURVEY_DATE]
        P = d["Downhole_Pressure"][SURVEY_DATE]
        z = np.array([float(k) for k in T["Data"]])
        t = np.array(list(T["Data"].values()), float)
        p = np.array([P["Data"][k] for k in T["Data"]], float)
        o = np.argsort(-z)
        wells[name] = dict(
            name=name, aut2=d["AUT2_Name"], track=track, feedzones=fz,
            design=design, z=z[o], T=t[o], P=p[o],
            T_stable=T.get("StableFlag"), P_stable=P.get("StableFlag"),
            role=("injector" if "Injection_Mass_Flow" in d else
                  "producer" if "Production_Mass_Flow" in d else "monitor"),
            transient_P={k: v["Data"] for k, v in d.get("Transient_Pressure", {}).items()},
            enthalpy=d.get("Transient_Enthalpy", {}).get("Data", {}),
            prod=d.get("Production_Mass_Flow", {}).get("Data", {}),
        )
    return wells


def xy_at_z(track, z):
    """Interpolate well track x,y at elevation z (track sorted top->bottom)."""
    zt = track[:, 2]
    o = np.argsort(zt)
    x = np.interp(z, zt[o], track[o, 0])
    y = np.interp(z, zt[o], track[o, 1])
    return x, y


def profile_at(w, z, var="T"):
    zz, vv = w["z"][::-1], w[var][::-1]
    if z < zz.min() or z > zz.max():
        return np.nan
    return float(np.interp(z, zz, vv))


# ----------------------------------------------------------------------------
# Analysis
# ----------------------------------------------------------------------------
def classify(w):
    """Simple thermal regime classification for a profile."""
    T, z, P = w["T"], w["z"], w["P"]
    imax = int(np.argmax(T))
    Tmax, zmax = T[imax], z[imax]
    inversion = Tmax - T[imax:].min() if imax < len(T) - 1 else 0.0
    margin = np.nanmin(tsat_c(np.maximum(P, 0.5)) - T)
    if Tmax < 40:
        return "Cold, near-isothermal - downflow / recharge"
    if margin < 5:
        return (f"At boiling point below {z[np.argmin(tsat_c(np.maximum(P, 0.5)) - T)]:.0f} mRL "
                "- two-phase / steam zone (shallow upflow)")
    if inversion > 20:
        return (f"Temperature reversal of {inversion:.0f} C below {zmax:.0f} mRL - "
                "lateral hot outflow over cooler fluid")
    if inversion > 5:
        return (f"Weak reversal ({inversion:.0f} C) below {zmax:.0f} mRL - "
                "distal outflow, conductive at depth")
    if Tmax < 150:
        return "Conductive / gently increasing - field margin"
    return "Convective - near-isothermal at depth"


def main():
    wells = load_wells()
    faults = json.loads((DATA / "Mapping" / "Mapping.json").read_text())

    # --------------------------------------------------------------- tables
    tgt_rows = []
    summ_rows = []
    fz_rows = []
    for w in wells.values():
        wh = w["track"][0]
        # per-point calibration targets
        for zi, Ti, Pi in zip(w["z"], w["T"], w["P"]):
            x, y = xy_at_z(w["track"], zi)
            c, li, lc, blk, above = model_block(x, y, zi)
            ts = float(tsat_c(Pi)) if Pi > 0 else np.nan
            in_fz = any(min(a, b) <= zi <= max(a, b) for a, b, _ in w["feedzones"])
            near_fz = any(min(a, b) - 30 <= zi <= max(a, b) + 30 for a, b, _ in w["feedzones"])
            cased = any(not s and min(a, b) <= zi <= max(a, b) for a, b, _, s, _ in w["design"])
            # Suggested weight for NS calibration
            if above:
                wt, note = 0.0, "above model surface in this column"
            elif zi > wh[2] - 60 and Ti < 60:
                wt, note = 0.3, "near-surface / wellhead - affected by casing & air"
            elif near_fz:
                wt, note = 1.0, "feed zone - most reliable formation temperature"
            elif cased:
                wt, note = 0.5, "behind casing - conductive interpretation"
            else:
                wt, note = 0.8, "open hole"
            tgt_rows.append(dict(
                well=w["name"], aut2=w["aut2"], x=round(x, 1), y=round(y, 1),
                elev_mRL=zi, T_C=Ti, P_bar=Pi, Tsat_C=round(ts, 1),
                boiling_margin_C=round(ts - Ti, 1), in_feedzone=in_fz,
                model_column=c, model_layer=li, layer_centre_mRL=lc, block=blk,
                column_surface_mRL=SURF.get(c), suggested_weight=wt, note=note))

        imax = int(np.argmax(w["T"]))
        x0, y0 = wh[0], wh[1]
        c0 = find_column(x0, y0)
        dz = w["z"][0] - w["z"][-1]
        # deep gradient over the bottom ~500 m
        sel = w["z"] <= w["z"][-1] + 500
        g_deep = (np.polyfit(w["z"][sel], w["T"][sel], 1)[0] * -1000
                  if sel.sum() > 1 and dz > 1000 else np.nan)
        # water level: extrapolate upper liquid pressures to 1 bar
        wl = np.nan
        upper = (w["P"] > 3) & (w["P"] < 30)
        if upper.sum() >= 2:
            a, b = np.polyfit(w["P"][upper], w["z"][upper], 1)
            wl = a * 1.013 + b
        summ_rows.append(dict(
            well=w["name"], aut2=w["aut2"], role=w["role"], x=x0, y=y0,
            wellhead_mRL=wh[2], TD_mRL=w["track"][-1, 2],
            deviated=len(w["track"]) > 2, column=c0,
            column_surface_mRL=SURF.get(c0),
            surface_misfit_m=round(wh[2] - SURF.get(c0, np.nan), 1),
            Tmax_C=w["T"][imax], Tmax_elev_mRL=w["z"][imax],
            T_bottom_C=w["T"][-1], deep_gradient_C_per_km=round(g_deep, 1),
            est_water_level_mRL=round(wl, 0),
            min_boiling_margin_C=round(float(np.nanmin(tsat_c(np.maximum(w["P"], 0.5)) - w["T"])), 1),
            regime=classify(w)))
        for a, b, q in w["feedzones"]:
            zm = 0.5 * (a + b)
            x, y = xy_at_z(w["track"], zm)
            c, li, lc, blk, _ = model_block(x, y, zm)
            fz_rows.append(dict(
                well=w["name"], top_mRL=max(a, b), bottom_mRL=min(a, b),
                mid_mRL=zm, quality=q, T_C=round(profile_at(w, zm, "T"), 1),
                P_bar=round(profile_at(w, zm, "P"), 1), column=c, layer=li,
                layer_centre_mRL=lc, block=blk))

    def write_csv(name, rows):
        with open(OUT / name, "w", newline="") as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            wr.writeheader()
            wr.writerows(rows)

    write_csv("NS_calibration_targets.csv", tgt_rows)
    write_csv("well_summary.csv", summ_rows)
    write_csv("feedzones.csv", fz_rows)

    # Temperature / pressure at model layer centres (for direct block comparison)
    lay_rows = []
    for L in LAYERS:
        if L["centre"] > 300:
            continue
        r = dict(layer=L["index"], centre_mRL=L["centre"])
        for w in wells.values():
            r[w["name"] + "_T"] = round(profile_at(w, L["centre"], "T"), 1)
        for w in wells.values():
            r[w["name"] + "_P"] = round(profile_at(w, L["centre"], "P"), 1)
        lay_rows.append(r)
    write_csv("T_P_at_layer_centres.csv", lay_rows)

    # Pre-production pressure cross-check (downhole profile vs first transient P)
    pchk = []
    for w in wells.values():
        for zs, series in w["transient_P"].items():
            first_date = list(series.keys())[0]
            pchk.append(dict(well=w["name"], elev_mRL=float(zs),
                             profile_P_bar=round(profile_at(w, float(zs), "P"), 1),
                             first_transient_date=first_date,
                             first_transient_P_bar=series[first_date]))
    write_csv("pressure_consistency_check.csv", pchk)

    # Hydrostatic diagnostic: measured dP/dz vs rho(T) g
    hyd = []
    for w in wells.values():
        z, P, T = w["z"], w["P"], w["T"]
        for k in range(len(z) - 1):
            if z[k] - z[k + 1] < 40:
                continue
            dpdz = (P[k + 1] - P[k]) / (z[k] - z[k + 1]) * 1e5  # Pa/m
            rho_eff = dpdz / G
            rho_T = rho_liq(0.5 * (T[k] + T[k + 1]))
            hyd.append(dict(well=w["name"], z_top=z[k], z_bot=z[k + 1],
                            rho_from_dPdz=round(rho_eff, 0),
                            rho_liquid_at_T=round(float(rho_T), 0),
                            ratio=round(rho_eff / rho_T, 3)))
    write_csv("hydrostatic_check.csv", hyd)

    # --------------------------------------------------------------- figures
    names = list(wells)
    cmap = plt.get_cmap("tab10")
    colour = {n: cmap(i % 10) for i, n in enumerate(names)}

    def draw_base(ax, labels=True):
        for cn, v in COLS.items():
            ax.fill(v[:, 0], v[:, 1], fill=False, lw=0.2, color="0.75")
        for k, feat in faults.items():
            parts = feat["Data"] if feat["type"] == "multiplot" else [feat["Data"]]
            isf = k.startswith("Fault")
            for part in parts:
                d = np.array(part)
                ax.plot(d[:, 0], d[:, 1], "--" if isf else "-",
                        color="k" if isf else "tab:blue", lw=1.2 if isf else 0.8)
            if isf:
                ax.text(np.mean(d[:, 0]), 6339500 if d[:, 1].max() > 6339500 else d[:, 1].max(),
                        k, fontsize=7)
        for n, w in wells.items():
            tr = w["track"]
            ax.plot(tr[:, 0], tr[:, 1], "-", color="k", lw=1)
            ax.plot(tr[0, 0], tr[0, 1], "o", color=colour[n], mec="k", ms=6)
            if labels:
                ax.annotate(n.replace("WELL-", "W"), tr[0, :2], xytext=(4, 4),
                            textcoords="offset points", fontsize=8)
        ax.set_xlim(2815000, 2825000)
        ax.set_ylim(6330000, 6340000)
        ax.set_aspect("equal")
        ax.ticklabel_format(useOffset=False, style="plain")
        ax.tick_params(labelsize=7)

    # 1. Location map
    fig, ax = plt.subplots(figsize=(8, 8))
    draw_base(ax)
    ax.set_title("Well locations, model grid (gTE14515), lake & faults")
    fig.savefig(FIG / "01_well_map.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # 2. Temperature profiles + BPD
    fig, axs = plt.subplots(1, 2, figsize=(13, 9), sharey=True)
    for n, w in wells.items():
        axs[0].plot(w["T"], w["z"], "-o", ms=3, color=colour[n], label=n)
        for a, b, q in w["feedzones"]:
            axs[0].plot([profile_at(w, 0.5 * (a + b))], [0.5 * (a + b)], "s",
                        ms=9, mfc="none", mec=colour[n], mew=2)
        axs[1].plot(w["P"], w["z"], "-o", ms=3, color=colour[n], label=n)
        axs[0].plot(tsat_c(np.maximum(w["P"], 1.0)), w["z"], ":", color=colour[n], lw=0.8)
    axs[0].set_xlabel("Temperature (°C)   [dotted = Tsat at measured P]")
    axs[1].set_xlabel("Pressure (bar abs)")
    axs[0].set_ylabel("Elevation (mRL)")
    axs[0].set_title("Stable downhole temperature, 29/07/2015 (squares = feed zones)")
    axs[1].set_title("Downhole pressure, 29/07/2015")
    for ax in axs:
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.savefig(FIG / "02_T_P_profiles.png", dpi=150, bbox_inches="tight")
    plt.close(fig)

    # 3. Individual well panels with layer boundaries
    fig, axs = plt.subplots(2, 5, figsize=(20, 11), sharey=True)
    for ax, (n, w) in zip(axs.flat, wells.items()):
        for L in LAYERS:
            if L["bottom"] > -3600:
                ax.axhline(L["bottom"], color="0.85", lw=0.5)
        for a, b, q, s, spec in w["design"]:
            if not s:
                ax.axhspan(min(a, b), max(a, b), xmin=0, xmax=0.03, color="0.3")
        for a, b, q in w["feedzones"]:
            ax.axhspan(min(a, b), max(a, b), color="orange", alpha=0.5)
            ax.text(5, 0.5 * (a + b), q, fontsize=7, va="center")
        ax.plot(w["T"], w["z"], "-o", ms=3, color="tab:red", label="T")
        ax.plot(tsat_c(np.maximum(w["P"], 1.0)), w["z"], "--", color="tab:blue",
                lw=1, label="Tsat(P)")
        ax.set_title(f"{n} ({w['role']})", fontsize=10)
        ax.set_xlim(0, 300)
        ax.set_ylim(-3000, 400)
        ax.grid(alpha=0.3)
    axs[0, 0].set_ylabel("Elevation (mRL)")
    axs[1, 0].set_ylabel("Elevation (mRL)")
    axs[0, 0].legend(fontsize=7)
    fig.suptitle("Temperature vs boiling curve; orange = feed zones; grey bar = blank casing; "
                 "grey lines = model layer bottoms")
    fig.savefig(FIG / "03_well_panels.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 4. Temperature maps at elevation slices
    elevs = [0, -250, -500, -750, -1000, -1500, -2000]
    fig, axs = plt.subplots(2, 4, figsize=(20, 10))
    for ax, ze in zip(axs.flat, elevs):
        xs, ys, ts, ls = [], [], [], []
        for n, w in wells.items():
            t = profile_at(w, ze)
            if np.isnan(t):
                continue
            x, y = xy_at_z(w["track"], ze)
            xs.append(x); ys.append(y); ts.append(t); ls.append(n)
        draw_base(ax, labels=False)
        if len(xs) >= 3:
            tri = mtri.Triangulation(xs, ys)
            cf = ax.tricontourf(tri, ts, levels=np.arange(0, 290, 20), cmap="inferno",
                                alpha=0.7)
        sc = ax.scatter(xs, ys, c=ts, vmin=0, vmax=280, cmap="inferno", s=80,
                        edgecolor="w", zorder=5)
        for x, y, t, n in zip(xs, ys, ts, ls):
            ax.annotate(f"{n.replace('WELL-', 'W')}\n{t:.0f}", (x, y), xytext=(5, -10),
                        textcoords="offset points", fontsize=7, color="k",
                        bbox=dict(fc="w", alpha=0.6, lw=0))
        ax.set_title(f"T at {ze} mRL (linear between wells)", fontsize=10)
    axs.flat[-1].axis("off")
    fig.colorbar(sc, ax=axs.flat[-1], fraction=0.5, label="Temperature (°C)")
    fig.savefig(FIG / "04_T_maps_by_elevation.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 5. Pressure at elevation, minus cold hydrostatic, and maps
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    for ax, ze in zip(axs, [-500, -1000, -1500]):
        draw_base(ax, labels=False)
        for n, w in wells.items():
            p = profile_at(w, ze, "P")
            if np.isnan(p):
                continue
            x, y = xy_at_z(w["track"], ze)
            ax.scatter(x, y, c=[p], cmap="viridis", vmin=60, vmax=170, s=90,
                       edgecolor="k", zorder=5)
            ax.annotate(f"{n.replace('WELL-', 'W')}\n{p:.1f} bar", (x, y),
                        xytext=(5, -10), textcoords="offset points", fontsize=7)
        ax.set_title(f"Pressure at {ze} mRL", fontsize=10)
    fig.savefig(FIG / "05_P_maps_by_elevation.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    # 6. Cross-sections: T contoured on wells projected onto a line
    sections = {
        "A_WSW-ENE": ((2817500, 6333000), (2823800, 6335200)),
        "B_SSW-NNE": ((2818500, 6331500), (2822500, 6338200)),
    }
    for sname, (p0, p1) in sections.items():
        p0, p1 = np.array(p0, float), np.array(p1, float)
        u = (p1 - p0) / np.linalg.norm(p1 - p0)
        nrm = np.array([-u[1], u[0]])
        fig, ax = plt.subplots(figsize=(12, 7))
        X, Z, T = [], [], []
        for n, w in wells.items():
            xs_, ys_ = xy_at_z(w["track"], w["z"])
            d = (np.c_[xs_, ys_] - p0) @ u
            off = np.abs((np.c_[xs_, ys_] - p0) @ nrm).mean()
            X += list(d); Z += list(w["z"]); T += list(w["T"])
            ax.plot(d, w["z"], "k-", lw=0.8)
            ax.text(d[-1], w["z"][-1] - 60, f"{n.replace('WELL-', 'W')}\n({off/1000:.1f} km off)",
                    ha="center", va="top", fontsize=7, rotation=90)
        tri = mtri.Triangulation(X, Z)
        # mask long skinny triangles
        xt, zt = np.array(X)[tri.triangles], np.array(Z)[tri.triangles]
        tri.set_mask((xt.max(1) - xt.min(1)) > 3500)
        cf = ax.tricontourf(tri, T, levels=np.arange(0, 290, 20), cmap="inferno")
        ax.tricontour(tri, T, levels=np.arange(0, 290, 20), colors="w", linewidths=0.4)
        fig.colorbar(cf, ax=ax, label="Temperature (°C)")
        ax.set_xlabel(f"Distance along section {sname} (m) — all wells projected")
        ax.set_ylabel("Elevation (mRL)")
        ax.set_title(f"Interpolated isotherms, section {sname} (indicative only)")
        fig.savefig(FIG / f"06_section_{sname}.png", dpi=130, bbox_inches="tight")
        plt.close(fig)
        # inset map of section line
    fig, ax = plt.subplots(figsize=(7, 7))
    draw_base(ax)
    for sname, (p0, p1) in sections.items():
        ax.plot([p0[0], p1[0]], [p0[1], p1[1]], "r-", lw=2)
        ax.text(p1[0], p1[1], sname, color="r")
    ax.set_title("Section lines")
    fig.savefig(FIG / "06_section_lines.png", dpi=130, bbox_inches="tight")
    plt.close(fig)

    print("Wrote:", *sorted(p.name for p in OUT.iterdir()), sep="\n  ")
    for r in summ_rows:
        print({k: r[k] for k in ("well", "column", "surface_misfit_m", "Tmax_C", "Tmax_elev_mRL",
                                 "T_bottom_C", "deep_gradient_C_per_km", "est_water_level_mRL",
                                 "min_boiling_margin_C", "regime")})
    for r in fz_rows:
        print(r)
    for r in pchk:
        print(r)


if __name__ == "__main__":
    main()
