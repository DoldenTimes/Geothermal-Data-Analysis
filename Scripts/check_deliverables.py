"""Acceptance check for a natural-state analysis run (CLAUDE.md v4).

Dataset-agnostic: the only input paths it uses come from
Output/Additional_Output/input_inventory.csv, which the analysis script writes (CLAUDE.md §2.1).

Run it from the repository root after the final run of Scripts/natural_state_analysis.py:

    python Scripts/check_deliverables.py            # static checks
    python Scripts/check_deliverables.py --rerun    # also re-run the analysis script in a temp copy

Options:
    --root DIR       repository to check (default: current directory)
    --baseline DIR   Baseline directory (default: <root>/Baseline). Inventory paths that start
                     with "Baseline/" are resolved against it.

The script only reads Baseline/. It writes Output/Additional_Output/acceptance_check.txt and
exits with 1 if any check FAILs, 0 otherwise. WARN items must be explained in the process log.

Checks (ID: what it catches):
  I1  input inventory: required roles present, paths exist, every Baseline file covered
  F1  required files (in their folders) and figures exist
  F2  figures_index.csv lists every standard figure, and every listed file exists
  F3  standard files are not misplaced; Calibration_Input holds only standard files
  S1  required columns, in the required order, in every standard CSV
  S2  one pressure-consistency row per well; hydrostatic band labels
  S3  UTF-8 with LF line endings; integer feed-zone numbers; ISO dates
  G1  survey-point blocks match the well-to-block file
  G2  feed-zone blocks match the well-to-block file
  G3  layer centres equal the grid file's layer centres, and lie inside the survey range
  G4  model_surface_mRL equals the grid file's surface section at the wellhead column
  T1  Tsat recomputed with IAPWS-IF97 Region 4 agrees with Tsat_degC to 0.1 K
  T2  boiling_margin_degC = Tsat_degC - T_degC
  D1  well_summary.csv values recomputed from the targets with the CLAUDE.md §4.4 definitions
  W1  weights in [0, 1]; zero above the model surface; weight_P only at feed zones
  W2  natural-state confidence differs between wells (WARN if uniform)
  W3  weights = location factor x confidence factor; verdict and confidence agree across files
  R1  every number in the report can be found in an output table (FAIL above 25 % unmatched)
  R2  the analysis script writes the report
  R3  (--rerun) re-running the script reproduces every CSV and the report
  R4  no unfilled template text ({name}, nan) in the report
  P1  the process log has the prompt, departures and the acceptance-check output
"""
import argparse
import bisect
import glob
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

# ----------------------------------------------------------------------------------------
# Specification (keep in step with CLAUDE.md v4, sections 2, 4.4, 5 and 8)
# ----------------------------------------------------------------------------------------
ROLES = {"well_data_json", "well_data_dat", "mapping", "well_groups", "geometry_dat",
         "well_blocks_json", "reference_rock", "reference_bc", "reference_results",
         "reference_other", "out_of_scope"}
SINGLE_ROLES = ["geometry_dat", "well_blocks_json"]          # exactly one path each
REQUIRED_ROLES = SINGLE_ROLES + ["well_data_json", "mapping"]

CAL = os.path.join("Output", "Calibration_Input")
ADD = os.path.join("Output", "Additional_Output")
FOLDER = {   # standard file -> the folder it must be in
    "NS_data_interpretation.md": CAL,
    "NS_calibration_targets.csv": CAL,
    "feedzones.csv": CAL,
    "T_P_at_layer_centres.csv": CAL,
    "well_summary.csv": CAL,
    "figures_index.csv": CAL,
    "input_inventory.csv": ADD,
    "pressure_consistency_check.csv": ADD,
    "hydrostatic_check.csv": ADD,
    "well_geometry_check.csv": ADD,
    "T_P_at_elevation_bands.csv": ADD,
    "run_console_log.txt": ADD,
    "process_log.md": ADD,
}
SCRIPT = os.path.join("Scripts", "natural_state_analysis.py")

REQUIRED_COLUMNS = {
    "input_inventory.csv": ["role", "path", "note"],
    "NS_calibration_targets.csv": [
        "well", "x", "y", "elevation_mRL", "T_degC", "P_bar", "Tsat_degC", "boiling_margin_degC",
        "column", "layer", "block", "in_feedzone", "behind_casing", "above_model_surface",
        "ns_confidence", "weight", "weight_reason", "weight_P"],
    "feedzones.csv": [
        "well", "feedzone", "top_mRL", "bottom_mRL", "mid_mRL", "quality", "T_degC", "P_bar",
        "Tsat_degC", "boiling_margin_degC", "block", "ns_confidence", "weight_T", "weight_P"],
    "T_P_at_layer_centres.csv": [
        "well", "layer", "layer_centre_mRL", "block", "T_degC", "P_bar", "Tsat_degC",
        "boiling_margin_degC"],
    "well_summary.csv": [
        "well", "Tmax_degC", "Tmax_elev_mRL", "T_at_TD_degC", "TD_mRL", "inversion_degC",
        "inversion_top_mRL", "inversion_base_mRL", "shallow_gradient_degC_per_km",
        "deep_gradient_degC_per_km", "z_100C_mRL", "z_200C_mRL", "min_boiling_margin_degC",
        "min_boiling_margin_elev_mRL", "water_level_mRL", "pivot_elev_mRL", "pivot_P_bar", "regime",
        "ns_verdict", "ns_confidence"],
    "pressure_consistency_check.csv": [
        "well", "gauge_elev_mRL", "survey_P_at_gauge_bar", "first_gauge_date", "first_gauge_P_bar",
        "diff_survey_minus_first_bar", "verdict", "confidence"],
    "hydrostatic_check.csv": [
        "well", "z_top_mRL", "z_bottom_mRL", "length_m", "band", "dP_measured_bar", "dP_IF97_bar",
        "ratio_measured_to_IF97", "rho_implied_kgm3"],
    "well_geometry_check.csv": [
        "well", "wellhead_mRL", "model_surface_mRL", "wellhead_minus_surface_m",
        "max_horizontal_deviation_m", "nearest_fault", "dist_nearest_fault_m",
        "nearest_water_body", "dist_nearest_water_body_m"],
    "T_P_at_elevation_bands.csv": [
        "well", "band_top_mRL", "band_bottom_mRL", "n_points", "T_mean_degC", "T_field_median_degC",
        "T_minus_field_median_degC", "P_mean_bar", "P_minus_trend_bar"],
    "figures_index.csv": [
        "figure_id", "file", "shows", "wells", "elevation_or_line", "data_sources", "claim_type"],
}
HYDRO_BANDS = {"cased", "open"}
# Standard figures in Calibration_Input/figures: id -> (globs, minimum count per glob).
# "per_well" = one per well in the targets; "gauged" = required only if any well has a gauge.
STANDARD_FIGURES = {
    "F1": [("well_map*.png", 1)],
    "F2": [("profile_*.png", "per_well"), ("profile_panels_*.png", 1)],
    "F3": [("maps_temperature*.png", 1), ("maps_pressure*.png", 1)],
    "F4": [("section_*.png", 1)],
    "F5": [("PT_state*.png", 1)],
    "F6": [("conceptual_model*.png", 1)],
    "F7": [("gauge_check*.png", "gauged")],
}
CLAIM_TYPES = {"obs", "interp", "mixed"}
CONF_FACTOR = {"high": 1.0, "medium": 0.75, "low": 0.5}
LOC_FACTORS = [0.0, 0.3, 0.5, 0.8, 1.0]
Q_FACTORS = [0.5, 0.8, 1.0]
REPORT_UNMATCHED_FAIL = 0.25
REPORT_UNMATCHED_WARN = 0.05
P_ATM = 1.01325
P_CRIT = 220.64

RESULTS = []


def record(cid, status, msg):
    RESULTS.append((cid, status, msg))


def is_true(v):
    return str(v).strip().lower() in ("true", "1", "yes")


def opath(root, name):
    """Absolute path of a standard output file in its required folder."""
    return os.path.join(root, FOLDER[name], name)


def read(root, name, **kw):
    return pd.read_csv(opath(root, name), **kw)


# ----------------------------------------------------------------------------------------
# IAPWS-IF97 Region 4 (saturation temperature), independent of the analysis script
# ----------------------------------------------------------------------------------------
_N = [None, 0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2, 0.12020824702470e5,
      -0.32325550322333e7, 0.14915108613530e2, -0.48232657361591e4, 0.40511340542057e6,
      -0.23855557567849, 0.65017534844798e3]


def tsat_c(p_bar):
    """Saturation temperature (degC) for pressure in bar, IF97 eq. 31. NaN outside its range."""
    p = p_bar / 10.0  # MPa
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


def validate_if97():
    # IF97 verification values for Region 4, plus the two values quoted in CLAUDE.md.
    for p_bar, t_k in [(1.0, 372.755919), (10.0, 453.035632), (100.0, 584.149488)]:
        assert abs(tsat_c(p_bar) + 273.15 - t_k) < 1e-5, "IF97 Region 4 implementation error"
    assert abs(tsat_c(1.01325) - 99.97) < 0.1 and abs(tsat_c(100.0) - 311.0) < 0.1


# ----------------------------------------------------------------------------------------
# Inventory and geometry helpers
# ----------------------------------------------------------------------------------------
def resolve(root, base, rel):
    rel = rel.replace("\\", "/").strip()
    p = os.path.join(root, rel)
    if not os.path.exists(p) and rel.startswith("Baseline/"):
        p = os.path.join(base, rel[len("Baseline/"):])
    return os.path.normpath(p)


def check_inventory(root, base):
    """I1. Returns {role: [absolute paths]} (empty dict if the inventory is unusable)."""
    path = opath(root, "input_inventory.csv")
    if not os.path.exists(path):
        record("I1", "FAIL", f"{FOLDER['input_inventory.csv']}/input_inventory.csv not found; "
                             "geometry checks cannot run")
        return {}
    inv = pd.read_csv(path).fillna("")
    probs, roles = [], {}
    for _, r in inv.iterrows():
        role = str(r.get("role", "")).strip()
        if role not in ROLES:
            probs.append(f"unknown role '{role}'")
            continue
        p = resolve(root, base, str(r.get("path", "")))
        if not os.path.exists(p):
            probs.append(f"{role}: path does not exist: {r.get('path')}")
            continue
        roles.setdefault(role, []).append(p)
    for role in REQUIRED_ROLES:
        if role not in roles:
            probs.append(f"required role '{role}' missing")
    for role in SINGLE_ROLES:
        if len(roles.get(role, [])) > 1:
            probs.append(f"role '{role}' has {len(roles[role])} paths; exactly one is required")
    record("I1", "FAIL" if probs else "PASS",
           "; ".join(probs[:6]) if probs else f"{len(inv)} inventory rows, all roles valid and paths present")

    covered = [p for ps in roles.values() for p in ps]
    uncovered = []
    for dirpath, _, files in os.walk(base):
        for f in files:
            fp = os.path.normpath(os.path.join(dirpath, f))
            if not any(fp == c or fp.startswith(c + os.sep) for c in covered):
                uncovered.append(os.path.relpath(fp, base))
    record("I1", "WARN" if uncovered else "PASS",
           f"{len(uncovered)} Baseline files not covered by the inventory; first: {uncovered[:5]}"
           if uncovered else "every Baseline file is covered by the inventory")
    return roles


def read_section(dat_path, header, ntok):
    """Lines of a MULgraph section with exactly ntok whitespace-separated tokens."""
    rows, on = [], False
    with open(dat_path) as f:
        for line in f:
            if line.startswith(header):
                on = True
                continue
            if on:
                p = line.split()
                if len(p) != ntok:
                    break
                rows.append(p)
    return rows


def read_layers(dat_path):
    """Layer section of a MULgraph grid file: 'index bottom_mRL centre_mRL'; layer 0 is the top."""
    layers = {int(p[0]): (float(p[1]), float(p[2])) for p in read_section(dat_path, "LAYERS", 3)}
    if not layers:
        raise ValueError("no LAYERS section found in the geometry file")
    return {i: dict(top=layers[i - 1][0], bottom=layers[i][0], centre=layers[i][1])
            for i in sorted(layers) if i > 0}


def read_surface(dat_path):
    """Surface section: 'column elevation_mRL'."""
    return {p[0]: float(p[1]) for p in read_section(dat_path, "SURFA", 2)}


def blocks_at(wells_json, well, z, tol=0.5):
    """Blocks whose interval contains z (both, if z is on a boundary)."""
    hits = set()
    for b in wells_json[well]["Blocks"]:
        top, bottom = max(b["interval"]), min(b["interval"])
        if bottom - tol <= z <= top + tol:
            hits.add(b["block"])
    return hits


def load_well_blocks(path):
    """The well-to-block file; some values may be stored as Python-literal strings."""
    import ast
    wj = json.load(open(path))
    for w, d in wj.items():
        for k, v in list(d.items()):
            if isinstance(v, str):
                try:
                    d[k] = ast.literal_eval(v)
                except (ValueError, SyntaxError):
                    pass
    return wj


# ----------------------------------------------------------------------------------------
# Files, figures and formats
# ----------------------------------------------------------------------------------------
def n_wells_and_gauged(root):
    n_wells, gauged = 0, False
    if os.path.exists(opath(root, "NS_calibration_targets.csv")):
        n_wells = read(root, "NS_calibration_targets.csv", usecols=["well"]).well.nunique()
    if os.path.exists(opath(root, "pressure_consistency_check.csv")):
        pc = read(root, "pressure_consistency_check.csv")
        gauged = "gauge_elev_mRL" in pc.columns and pc.gauge_elev_mRL.notna().any()
    return n_wells, gauged


def figure_files(figdir, pat):
    files = glob.glob(os.path.join(figdir, pat))
    if pat == "profile_*.png":   # group panels must not stand in for per-well profiles
        files = [f for f in files if not os.path.basename(f).startswith("profile_panels_")]
    return files


def check_files(root):
    req = [SCRIPT] + [os.path.join(FOLDER[n], n) for n in FOLDER]
    missing = [p.replace("\\", "/") for p in req if not os.path.exists(os.path.join(root, p))]
    record("F1", "FAIL" if missing else "PASS",
           f"missing files: {missing}" if missing else f"all {len(req)} required files present in their folders")
    n_wells, gauged = n_wells_and_gauged(root)
    figdir = os.path.join(root, CAL, "figures")
    short = []
    for fid, specs in STANDARD_FIGURES.items():
        for pat, n in specs:
            need = n_wells if n == "per_well" else (1 if gauged else 0) if n == "gauged" else n
            k = len(figure_files(figdir, pat))
            if k < need:
                short.append(f"{fid} {pat}: {k}/{need}")
    record("F1", "FAIL" if short else "PASS",
           f"standard figures missing from {CAL.replace(os.sep, '/')}/figures: {short}" if short
           else "all standard figure types present")


def check_figure_index(root):
    path = opath(root, "figures_index.csv")
    if not os.path.exists(path):
        record("F2", "FAIL", "figures_index.csv not found")
        return
    idx = pd.read_csv(path, dtype=str).fillna("")
    _, gauged = n_wells_and_gauged(root)
    probs = []
    need = [f for f in STANDARD_FIGURES if f != "F7" or gauged]
    ids = set(idx.get("figure_id", pd.Series(dtype=str)).astype(str).str.strip())
    miss = [f for f in need if f not in ids]
    if miss:
        probs.append(f"figure ids not indexed: {miss}")
    bad_ids = sorted(ids - set(STANDARD_FIGURES))
    if bad_ids:
        probs.append(f"unknown figure ids: {bad_ids}")
    listed = set()
    for f in idx.get("file", []):
        p = os.path.normpath(os.path.join(root, str(f).strip()))
        listed.add(p)
        if not os.path.exists(p):
            probs.append(f"indexed file does not exist: {f}")
    if "claim_type" in idx.columns:
        bad = sorted(set(idx.claim_type.astype(str).str.strip()) - CLAIM_TYPES)
        if bad:
            probs.append(f"claim_type values other than {sorted(CLAIM_TYPES)}: {bad}")
    record("F2", "FAIL" if probs else "PASS",
           "; ".join(probs[:6]) if probs else f"figure index lists {len(idx)} files covering {sorted(ids)}")
    unlisted = [os.path.basename(p) for p in glob.glob(os.path.join(root, CAL, "figures", "*.png"))
                if os.path.normpath(p) not in listed]
    record("F2", "WARN" if unlisted else "PASS",
           f"{len(unlisted)} figures in Calibration_Input/figures not in the index: {unlisted[:6]}"
           if unlisted else "every Calibration_Input figure is indexed")


def check_placement(root):
    out = os.path.join(root, "Output")
    misplaced = []
    for dirpath, _, files in os.walk(out):
        for f in files:
            if f in FOLDER and os.path.normpath(dirpath) != os.path.normpath(os.path.join(root, FOLDER[f])):
                misplaced.append(os.path.relpath(os.path.join(dirpath, f), root).replace("\\", "/"))
    record("F3", "FAIL" if misplaced else "PASS",
           f"standard files outside their required folder: {misplaced[:6]}" if misplaced
           else "no standard file is misplaced")
    cal = os.path.join(root, CAL)
    extra = []
    if os.path.isdir(cal):
        allowed = {n for n, d in FOLDER.items() if d == CAL} | {"figures"}
        extra = [f for f in os.listdir(cal) if f not in allowed]
        figdir = os.path.join(cal, "figures")
        if os.path.isdir(figdir):
            pats = [p for specs in STANDARD_FIGURES.values() for p, _ in specs]
            extra += ["figures/" + f for f in os.listdir(figdir)
                      if not any(glob.fnmatch.fnmatch(f, p) for p in pats)]
    record("F3", "WARN" if extra else "PASS",
           f"non-standard items in Calibration_Input (belong in Additional_Output): {extra[:6]}"
           if extra else "Calibration_Input holds only standard files")


def check_formats(root):
    probs, crlf = [], []
    for folder in (CAL, ADD):
        for path in glob.glob(os.path.join(root, folder, "**", "*"), recursive=True):
            if not path.endswith((".csv", ".md", ".txt")) or path.endswith("acceptance_check.txt"):
                continue
            raw = open(path, "rb").read()
            rel = os.path.relpath(path, root).replace("\\", "/")
            try:
                raw.decode("utf-8")
            except UnicodeDecodeError:
                probs.append(f"{rel}: not UTF-8")
            if b"\r\n" in raw:
                crlf.append(rel)
    if os.path.exists(opath(root, "feedzones.csv")):
        fz = read(root, "feedzones.csv")
        for well, g in fz.groupby("well", sort=False):
            vals = g.feedzone.tolist()
            try:
                ints = [int(v) for v in vals]
                ok = all(float(v) == i for v, i in zip(vals, ints)) and ints == list(range(1, len(ints) + 1))
            except (ValueError, TypeError):
                ok = False
            if not ok:
                probs.append(f"feedzones.csv: {well} feedzone values {vals} are not 1, 2, 3, ...")
    iso = re.compile(r"^\d{4}-\d{2}-\d{2}$")
    if os.path.exists(opath(root, "pressure_consistency_check.csv")):
        pc = read(root, "pressure_consistency_check.csv", dtype=str)
        for c in [c for c in pc.columns if c.endswith("_date")]:
            bad = [v for v in pc[c].dropna() if not iso.match(v.strip())]
            if bad:
                probs.append(f"pressure_consistency_check.csv {c}: non-ISO dates {bad[:3]}")
    record("S3", "FAIL" if probs else "PASS",
           "; ".join(probs[:6]) if probs else "UTF-8, integer feed zones, ISO dates")
    # WARN, not FAIL: git with core.autocrlf=true can add CRLF on checkout after a clean LF write
    record("S3", "WARN" if crlf else "PASS",
           f"{len(crlf)} text outputs have CRLF line endings (write with newline='\\n' / "
           f"lineterminator='\\n'): {crlf[:4]}" if crlf else "LF line endings")


def check_columns(root):
    for name, cols in REQUIRED_COLUMNS.items():
        path = opath(root, name)
        if not os.path.exists(path):
            continue
        have = list(pd.read_csv(path, nrows=0).columns)
        missing = [c for c in cols if c not in have]
        if missing:
            record("S1", "FAIL", f"{name}: missing columns {missing}")
        elif have[:len(cols)] != cols:
            record("S1", "FAIL", f"{name}: required columns present but not first and in order")
        else:
            record("S1", "PASS", f"{name}: columns OK")


def check_tables(root):
    wells = set(read(root, "NS_calibration_targets.csv").well)
    pc = read(root, "pressure_consistency_check.csv")
    diff = wells.symmetric_difference(set(pc.well))
    dup = pc.well.duplicated().any()
    record("S2", "FAIL" if (diff or dup) else "PASS",
           f"pressure_consistency_check.csv must have one row per well; mismatched: {sorted(diff)}"
           f"{'; duplicate wells' if dup else ''}" if (diff or dup) else "one pressure-consistency row per well")
    hy = read(root, "hydrostatic_check.csv")
    if "band" in hy.columns:
        bad = sorted(set(hy.band.dropna().astype(str)) - HYDRO_BANDS)
        record("S2", "WARN" if bad else "PASS",
               f"hydrostatic band labels other than {sorted(HYDRO_BANDS)}: {bad}" if bad
               else "hydrostatic band labels OK")


# ----------------------------------------------------------------------------------------
# Geometry
# ----------------------------------------------------------------------------------------
def check_geometry(root, roles):
    if not roles.get("geometry_dat") or not roles.get("well_blocks_json"):
        record("G*", "FAIL", "geometry checks skipped: inventory lacks geometry_dat or well_blocks_json")
        return
    wj = load_well_blocks(roles["well_blocks_json"][0])
    dat = roles["geometry_dat"][0]
    layers = read_layers(dat)
    tgt = read(root, "NS_calibration_targets.csv")

    # G1: survey points
    bad, unknown = [], sorted(set(tgt.well) - set(wj))
    for _, r in tgt.iterrows():
        if is_true(r.get("above_model_surface", False)) or r.well not in wj:
            continue
        ok = blocks_at(wj, r.well, r.elevation_mRL)
        if ok and r.block not in ok:
            bad.append(f"{r.well}@{r.elevation_mRL}: {r.block} (block file {'/'.join(sorted(ok))})")
    msg = (f"{len(bad)}/{len(tgt)} survey-point blocks differ from the well-to-block file; first: {bad[:3]}"
           if bad else f"all {len(tgt)} survey-point blocks match the well-to-block file")
    if unknown:
        msg += f"; wells not in the block file: {unknown}"
    record("G1", "FAIL" if (bad or unknown) else "PASS", msg)

    # G2: feed zones, matched by order within each well
    fz = read(root, "feedzones.csv")
    if "block" not in fz.columns:
        record("G2", "FAIL", "feedzones.csv has no 'block' column")
    else:
        bad = []
        for well, g in fz.groupby("well", sort=False):
            fzj = list(wj.get(well, {}).get("Feed_Zones", {}).values())
            if len(fzj) != len(g):
                bad.append(f"{well}: {len(g)} feed zones in feedzones.csv, {len(fzj)} in the block file")
                continue
            for i, (_, r) in enumerate(g.iterrows()):
                ok = set(fzj[i].keys())
                if r.block not in ok:
                    bad.append(f"{well} feed zone {i + 1}: {r.block} (block file {'/'.join(sorted(ok))})")
        record("G2", "FAIL" if bad else "PASS",
               f"feed-zone blocks differ: {bad[:6]}" if bad else f"all {len(fz)} feed-zone blocks match")

    # G3: layer centres and no extrapolation
    lc = read(root, "T_P_at_layer_centres.csv")
    wrong = []
    for lay, c in lc.groupby("layer").layer_centre_mRL.first().items():
        true_c = layers.get(int(lay), {}).get("centre")
        if true_c is None or abs(c - true_c) > 0.5:
            wrong.append((int(lay), c, true_c))
    if wrong:
        errs = [abs(a - b) for _, a, b in wrong if b is not None]
        record("G3", "FAIL", f"{len(wrong)} layer centres differ from the grid file "
                             f"(largest error {max(errs) if errs else 'n/a'} m); "
                             f"e.g. layer {wrong[0][0]}: {wrong[0][1]} vs {wrong[0][2]}")
    else:
        record("G3", "PASS", f"all {lc.layer.nunique()} layer centres match the grid file")
    rng = tgt.groupby("well").elevation_mRL.agg(["min", "max"])
    outside = [(w, z) for w, z in zip(lc.well, lc.layer_centre_mRL)
               if w in rng.index and not (rng.loc[w, "min"] - 0.5 <= z <= rng.loc[w, "max"] + 0.5)]
    record("G3", "FAIL" if outside else "PASS",
           f"{len(outside)} layer-centre rows lie outside the survey range (extrapolation)"
           if outside else "no layer-centre row lies outside its well's survey range")

    # G4: model surface at the wellhead column, from the surface section
    surf = read_surface(dat)
    grid_top = layers[min(layers)]["top"]
    widths = {len(c) for c in surf}
    geo = read(root, "well_geometry_check.csv")
    bad, n = [], 0
    for _, r in geo.iterrows():
        blocks = wj.get(r.well, {}).get("Blocks", [])
        if not blocks:
            bad.append(f"{r.well}: not in the well-to-block file")
            continue
        top_block = max(blocks, key=lambda b: max(b["interval"]))["block"]
        # column names in the surface section share one width; otherwise strip the layer number
        col = top_block[:next(iter(widths))] if len(widths) == 1 else re.sub(r"\d+$", "", top_block)
        truth = surf.get(col, grid_top)
        n += 1
        if pd.isna(r.model_surface_mRL) or abs(r.model_surface_mRL - truth) > 0.5:
            bad.append(f"{r.well} (column {col}): {r.model_surface_mRL} vs surface section {truth}")
    record("G4", "FAIL" if bad else "PASS",
           f"model_surface_mRL differs from the surface section for {len(bad)}/{len(geo)} wells: {bad[:4]}"
           if bad else f"model_surface_mRL matches the surface section for all {n} wells")


# ----------------------------------------------------------------------------------------
# Thermodynamics
# ----------------------------------------------------------------------------------------
def check_thermo(root):
    tgt = read(root, "NS_calibration_targets.csv")
    errs, bm_errs = [], []
    for _, r in tgt.iterrows():
        ts = tsat_c(r.P_bar)
        if math.isnan(ts) or pd.isna(r.Tsat_degC):
            continue
        if abs(ts - r.Tsat_degC) > 0.1:
            errs.append(f"{r.well}@{r.elevation_mRL}: {r.Tsat_degC} vs {ts:.2f}")
        if pd.notna(r.boiling_margin_degC) and abs(r.Tsat_degC - r.T_degC - r.boiling_margin_degC) > 0.06:
            bm_errs.append(f"{r.well}@{r.elevation_mRL}")
    record("T1", "FAIL" if errs else "PASS",
           f"Tsat differs from IF97 by > 0.1 K at {len(errs)} points: {errs[:3]}" if errs
           else "Tsat agrees with an independent IF97 Region 4 within 0.1 K at every point")
    record("T2", "FAIL" if bm_errs else "PASS",
           f"boiling margin != Tsat - T at {len(bm_errs)} points" if bm_errs else "boiling margins consistent")


# ----------------------------------------------------------------------------------------
# D1: §4.4 definitions recomputed from the targets
# ----------------------------------------------------------------------------------------
def _fit_slope(z, v):
    z, v = np.asarray(z, float), np.asarray(v, float)
    return np.polyfit(z, v, 1)


def _crossing(z, t, x):
    for i in range(len(t) - 1):
        if t[i] < x <= t[i + 1]:
            return z[i] + (x - t[i]) * (z[i + 1] - z[i]) / (t[i + 1] - t[i])
    return float("nan")


def expected_summary(g, fz_bottom):
    """Per-well characteristics from one well's target rows (CLAUDE.md §4.4)."""
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
    if len(zl):
        w = zl >= zl[0] - 300
        if w.sum() >= 3:
            e["shallow_gradient_degC_per_km"] = -_fit_slope(zl[w], tl[w])[0] * 1000

    bounds = [b for b in (inv[2], fz_bottom) if not np.isnan(b)]
    sel = z <= min(bounds) if bounds else np.ones_like(z, bool)
    e["deep_gradient_degC_per_km"] = np.nan
    if sel.sum() >= 3 and z[sel].max() - z[sel].min() >= 300:
        e["deep_gradient_degC_per_km"] = -_fit_slope(z[sel], t[sel])[0] * 1000

    e["z_100C_mRL"], e["z_200C_mRL"] = _crossing(z, t, 100.0), _crossing(z, t, 200.0)

    ok = ~np.isnan(bml)
    if ok.any():
        j = int(np.argmin(np.where(ok, bml, np.inf)))
        e["min_boiling_margin_degC"], e["min_boiling_margin_elev_mRL"] = bml[j], zl[j]
    else:
        e["min_boiling_margin_degC"] = e["min_boiling_margin_elev_mRL"] = np.nan

    e["water_level_mRL"] = np.nan
    if len(zl) >= 3:
        b, a = _fit_slope(zl[:3], pl[:3])
        if b != 0:
            e["water_level_mRL"] = (P_ATM - a) / b

    if ok.any() and np.nanmin(bml) <= 1:
        reg = "two-phase"
    elif e["Tmax_degC"] < 100:
        reg = "recharge"
    elif not np.isnan(inv[0]) and inv[0] >= 10:
        reg = "outflow"
    elif e["Tmax_degC"] >= 200 and np.isnan(inv[0]):
        reg = "upflow"
    else:
        reg = "conductive"
    e["regime"] = reg
    return e


TOL = {"Tmax_degC": 0.15, "Tmax_elev_mRL": 0.6, "T_at_TD_degC": 0.15, "TD_mRL": 0.6,
       "inversion_degC": 0.15, "inversion_top_mRL": 0.6, "inversion_base_mRL": 0.6,
       "shallow_gradient_degC_per_km": "grad", "deep_gradient_degC_per_km": "grad",
       "z_100C_mRL": 1.5, "z_200C_mRL": 1.5, "min_boiling_margin_degC": 0.15,
       "min_boiling_margin_elev_mRL": 0.6, "water_level_mRL": 2.0}


def _conf_factor(v):
    return CONF_FACTOR.get(str(v).strip().lower())


def check_definitions(root):
    tgt = read(root, "NS_calibration_targets.csv")
    ws = read(root, "well_summary.csv").set_index("well")
    fz = read(root, "feedzones.csv")
    probs, n = [], 0
    for well, g in tgt.groupby("well", sort=False):
        if well not in ws.index:
            probs.append(f"{well}: not in well_summary.csv")
            continue
        fzw = fz[fz.well == well]
        fzb = fzw.bottom_mRL.min() if len(fzw) else np.nan
        e = expected_summary(g, fzb)
        row = ws.loc[well]
        for col, tol in TOL.items():
            exp, got = e[col], pd.to_numeric(row.get(col), errors="coerce")
            n += 1
            if np.isnan(exp) and pd.isna(got):
                continue
            if np.isnan(exp) != pd.isna(got):
                probs.append(f"{well} {col}: {got} (expected {'blank' if np.isnan(exp) else round(exp, 2)})")
                continue
            t = max(1.0, 0.01 * abs(exp)) if tol == "grad" else tol
            if abs(got - exp) > t:
                probs.append(f"{well} {col}: {got} (expected {exp:.2f})")
        reg = str(row.get("regime", "")).lower()
        n += 1
        keys = ["recharge", "downflow"] if e["regime"] == "recharge" else [e["regime"]]
        if not any(k in reg for k in keys):
            probs.append(f"{well} regime: '{row.get('regime')}' (expected {e['regime']})")
        # pivot: highest weight_P quality factor, then thickest, then deepest
        if len(fzw):
            cf = _conf_factor(fzw.ns_confidence.iloc[0]) or 1.0
            q = (fzw.weight_P / cf).round(2)
            cand = fzw[q == q.max()].assign(th=lambda d: (d.top_mRL - d.bottom_mRL).abs())
            cand = cand.sort_values(["th", "mid_mRL"], ascending=[False, True])
            pv = cand.iloc[0]
            n += 1
            if abs(pd.to_numeric(row.get("pivot_elev_mRL"), errors="coerce") - pv.mid_mRL) > 0.6 or \
                    abs(pd.to_numeric(row.get("pivot_P_bar"), errors="coerce") - pv.P_bar) > 0.15:
                probs.append(f"{well} pivot: {row.get('pivot_elev_mRL')} mRL / {row.get('pivot_P_bar')} bar "
                             f"(expected feed zone {pv.feedzone}: {pv.mid_mRL} mRL / {pv.P_bar} bar)")
    record("D1", "FAIL" if probs else "PASS",
           f"{len(probs)}/{n} well_summary values differ from the §4.4 definitions; first: {probs[:5]}"
           if probs else f"all {n} well_summary values match the §4.4 definitions")


# ----------------------------------------------------------------------------------------
# Weights
# ----------------------------------------------------------------------------------------
def _near(v, allowed, tol=0.004):
    return any(abs(v - a) <= tol for a in allowed)


def check_weights(root):
    tgt = read(root, "NS_calibration_targets.csv")
    probs = []
    if not tgt.weight.between(0, 1).all():
        probs.append("weight outside [0, 1]")
    above = tgt.above_model_surface.map(is_true)
    if (tgt.loc[above, "weight"] > 0).any():
        probs.append("non-zero weight above the model surface")
    if "weight_P" in tgt.columns:
        infz = tgt.in_feedzone.map(is_true)
        n = int((tgt.loc[~infz, "weight_P"].fillna(0) > 0).sum())
        if n:
            probs.append(f"weight_P > 0 at {n} points outside feed zones")
    else:
        probs.append("no weight_P column")
    record("W1", "FAIL" if probs else "PASS", "; ".join(probs) if probs else "weights valid")
    conf = tgt.groupby("well").ns_confidence.first()
    if conf.nunique() == 1:
        record("W2", "WARN", f"ns_confidence is '{conf.iloc[0]}' for every well, so it changes no relative "
                             "weight; the report must justify this")
    else:
        record("W2", "PASS", f"ns_confidence takes {conf.nunique()} values across wells")


def check_weight_rules(root):
    tgt = read(root, "NS_calibration_targets.csv")
    fz = read(root, "feedzones.csv")
    ws = read(root, "well_summary.csv").set_index("well")
    pc = read(root, "pressure_consistency_check.csv").set_index("well")
    probs = []
    # location factor x confidence factor, by the ordered rule in §5.5
    bad = []
    for _, r in tgt.iterrows():
        cf = _conf_factor(r.ns_confidence)
        if cf is None:
            bad.append(f"{r.well}: ns_confidence '{r.ns_confidence}'")
            continue
        loc = r.weight / cf
        if is_true(r.above_model_surface):
            allowed = [0.0]
        elif is_true(r.in_feedzone):
            allowed = [1.0, 0.3]
        elif is_true(r.behind_casing):
            allowed = [0.5, 0.3]
        else:
            allowed = [0.8, 0.5, 0.3]
        if not _near(loc, [a for a in allowed]):
            bad.append(f"{r.well}@{r.elevation_mRL}: weight {r.weight} = {loc:.3f} x {cf} (allowed {allowed})")
        if is_true(r.in_feedzone) and not is_true(r.above_model_surface) and pd.notna(r.weight_P) \
                and not _near(r.weight_P / cf, Q_FACTORS):
            bad.append(f"{r.well}@{r.elevation_mRL}: weight_P {r.weight_P} is not a quality factor x {cf}")
    if bad:
        probs.append(f"{len(bad)} target weights break §5.5; first: {bad[:3]}")
    for _, r in fz.iterrows():
        cf = _conf_factor(r.ns_confidence)
        if cf is None or not _near(r.weight_P / cf, Q_FACTORS) or not _near(r.weight_T / cf, LOC_FACTORS):
            probs.append(f"feedzones.csv {r.well} feed zone {r.feedzone}: weight_T {r.weight_T}, "
                         f"weight_P {r.weight_P} inconsistent with confidence '{r.ns_confidence}'")
    # one confidence and verdict per well across files
    per_well = tgt.groupby("well").ns_confidence.agg(lambda s: set(s.astype(str).str.lower()))
    for well, s in per_well.items():
        vals = {"targets": s}
        if well in ws.index:
            vals["well_summary"] = {str(ws.loc[well, "ns_confidence"]).lower()}
        if well in pc.index:
            vals["pressure_consistency"] = {str(pc.loc[well, "confidence"]).lower()}
        f = fz[fz.well == well]
        if len(f):
            vals["feedzones"] = set(f.ns_confidence.astype(str).str.lower())
        allv = set().union(*vals.values())
        if len(allv) > 1:
            probs.append(f"{well}: ns_confidence differs across files {vals}")
        if well in ws.index and well in pc.index and \
                str(ws.loc[well, "ns_verdict"]).strip().lower() != str(pc.loc[well, "verdict"]).strip().lower():
            probs.append(f"{well}: ns_verdict '{ws.loc[well, 'ns_verdict']}' vs pressure_consistency "
                         f"verdict '{pc.loc[well, 'verdict']}'")
    record("W3", "FAIL" if probs else "PASS",
           "; ".join(probs[:5]) if probs else "weights follow §5.5; verdict and confidence agree across files")


# ----------------------------------------------------------------------------------------
# Report
# ----------------------------------------------------------------------------------------
NUM = re.compile(r"(?<![\w.\-−])[-−]?\d[\d,]*(?:\.\d+)?")


def _numbers(text):
    vals = []
    for tok in NUM.findall(text):
        t = tok.replace("−", "-").replace(",", "")
        try:
            vals.append((tok, float(t), len(t.split(".")[1]) if "." in t else 0))
        except ValueError:
            pass
    return vals


def output_text_files(root):
    files = []
    for folder in (CAL, ADD):
        files += glob.glob(os.path.join(root, folder, "**", "*"), recursive=True)
    return [f for f in files if f.endswith((".csv", ".json", ".txt")) and not f.endswith("acceptance_check.txt")]


def check_report(root):
    raw = open(opath(root, "NS_data_interpretation.md"), encoding="utf-8").read()
    rep = re.sub(r"§\s*\d+(\.\d+)*", " ", raw)                 # section references
    rep = re.sub(r"`[^`]*`", " ", rep)                          # file names and code
    rep = re.sub(r"\d{1,4}[/-]\d{1,2}[/-]\d{1,4}", " ", rep)    # dates
    rep = re.sub(r"\(.*?\.(png|csv|md|json|txt)\)", " ", rep)  # links
    pool = set()
    for path in output_text_files(root):
        for _, v, _ in _numbers(open(path, encoding="utf-8", errors="ignore").read()):
            pool.add(abs(v))
    pool = sorted(pool)

    def found(v, d):
        tol = 0.5 * 10 ** (-d) + 1e-9
        i = bisect.bisect_left(pool, abs(v) - tol)
        return i < len(pool) and pool[i] <= abs(v) + tol

    checked, missing = 0, []
    for tok, v, d in _numbers(rep):
        if d == 0 and abs(v) < 100:
            continue  # small integers: counts, list numbers, well numbers
        checked += 1
        if not found(v, d):
            missing.append(tok)
    frac = len(missing) / checked if checked else 0
    status = "FAIL" if frac > REPORT_UNMATCHED_FAIL else "WARN" if frac > REPORT_UNMATCHED_WARN else "PASS"
    record("R1", status, f"{len(missing)}/{checked} report numbers ({frac:.0%}) not found in any output table; "
                         f"first: {missing[:12]}")
    src = open(os.path.join(root, SCRIPT), encoding="utf-8").read()
    ok = "NS_data_interpretation.md" in src
    record("R2", "PASS" if ok else "FAIL",
           "analysis script writes the report" if ok
           else "analysis script never names NS_data_interpretation.md: the report is not generated")
    prose = re.sub(r"```.*?```", " ", raw, flags=re.S)
    prose = re.sub(r"`[^`]*`", " ", prose)
    hits = re.findall(r"\{[A-Za-z_][^{}\n]{0,60}\}", prose) + re.findall(r"\bnan\b", prose, flags=re.I)
    record("R4", "FAIL" if hits else "PASS",
           f"unfilled template text in the report: {hits[:6]}" if hits else "no unfilled template text in the report")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def check_rerun(root, base):
    tmp = tempfile.mkdtemp(prefix="nscheck_")
    try:
        for folder in (CAL, ADD):
            os.makedirs(os.path.join(tmp, folder, "figures"))
        shutil.copytree(base, os.path.join(tmp, "Baseline"))
        shutil.copytree(os.path.join(root, "Scripts"), os.path.join(tmp, "Scripts"))
        p = subprocess.run([sys.executable, SCRIPT], cwd=tmp, capture_output=True, text=True, timeout=1800)
        if p.returncode != 0:
            record("R3", "FAIL", f"re-run exited {p.returncode}: {p.stderr.strip()[-300:]}")
            return
        diff = []
        files = [opath(root, "NS_data_interpretation.md")]
        for folder in (CAL, ADD):
            files += sorted(glob.glob(os.path.join(root, folder, "*.csv")))
        for f in files:
            rel = os.path.relpath(f, root)
            g = os.path.join(tmp, rel)
            if not os.path.exists(g):
                diff.append(rel.replace("\\", "/") + " (not regenerated)")
            elif sha(f) != sha(g):
                a = open(f, encoding="utf-8").read().splitlines()
                b = open(g, encoding="utf-8").read().splitlines()
                n = sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))
                if not (f.endswith(".md") and n <= 1):   # allow one timestamp line in the report
                    diff.append(f"{os.path.basename(f)} ({n} lines differ)")
        record("R3", "FAIL" if diff else "PASS",
               f"re-run does not reproduce: {diff}" if diff else "re-run reproduces every CSV and the report")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def check_process_log(root):
    txt = open(opath(root, "process_log.md"), encoding="utf-8").read().lower()
    miss = [k for k in ("prompt", "verbatim", "departure", "acceptance") if k not in txt]
    record("P1", "WARN" if miss else "PASS",
           f"process log lacks: {miss}" if miss else "process log has prompt, departures and acceptance check")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".")
    ap.add_argument("--baseline", default=None)
    ap.add_argument("--rerun", action="store_true")
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    base = os.path.abspath(a.baseline or os.path.join(root, "Baseline"))
    validate_if97()

    roles = {}
    try:
        roles = check_inventory(root, base)
    except Exception as e:
        record("I1", "FAIL", f"inventory check crashed: {type(e).__name__}: {e}")
    for fn in (check_files, check_figure_index, check_placement, check_columns, check_tables,
               check_formats, check_thermo, check_definitions, check_weights, check_weight_rules,
               check_report, check_process_log):
        try:
            fn(root)
        except Exception as e:  # a crash in one check must not hide the others
            record(fn.__name__, "FAIL", f"check crashed: {type(e).__name__}: {e}")
    try:
        check_geometry(root, roles)
    except Exception as e:
        record("G*", "FAIL", f"geometry check crashed: {type(e).__name__}: {e}")
    if a.rerun:
        check_rerun(root, base)

    lines = [f"{s:4}  {c:3}  {m}" for c, s, m in RESULTS]
    n_fail = sum(s == "FAIL" for _, s, _ in RESULTS)
    n_warn = sum(s == "WARN" for _, s, _ in RESULTS)
    lines.append(f"\n{n_fail} FAIL, {n_warn} WARN, {len(RESULTS) - n_fail - n_warn} PASS")
    text = "\n".join(lines)
    print(text)
    out = os.path.join(root, ADD, "acceptance_check.txt")
    if os.path.isdir(os.path.dirname(out)):
        with open(out, "w", encoding="utf-8", newline="\n") as f:
            f.write(text + "\n")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
