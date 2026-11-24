"""Acceptance check for a natural-state analysis run (DRAFT for CLAUDE.md v3).

Dataset-agnostic: the only input paths it uses come from Output/input_inventory.csv, which
the analysis script writes (CLAUDE.md v3, section 2.1).

Intended location: Scripts/check_deliverables.py. Run it from the repository root after
the final run of Scripts/natural_state_analysis.py:

    python Scripts/check_deliverables.py            # static checks
    python Scripts/check_deliverables.py --rerun    # also re-run the analysis script in a temp copy

Options:
    --root DIR       repository to check (default: current directory)
    --baseline DIR   Baseline directory (default: <root>/Baseline). Inventory paths that start
                     with "Baseline/" are resolved against it.

The script only reads Baseline/. It writes Output/acceptance_check.txt and exits with 1 if
any check FAILs, 0 otherwise. WARN items must be explained in the process log.

Checks (ID: what it catches):
  I1  input inventory: required roles present, paths exist, every Baseline file covered
  F1  required files and figures exist
  S1  required columns, in the required order, in every standard CSV
  S2  one pressure-consistency row per well; hydrostatic band labels
  G1  survey-point blocks match the well-to-block file
  G2  feed-zone blocks match the well-to-block file
  G3  layer centres equal the grid file's layer centres, and lie inside the survey range
  T1  Tsat recomputed with IAPWS-IF97 Region 4 agrees with Tsat_degC to 0.1 K
  T2  boiling_margin_degC = Tsat_degC - T_degC
  W1  weights in [0, 1]; zero above the model surface; weight_P only at feed zones
  W2  natural-state confidence differs between wells (WARN if uniform)
  R1  every number in the report can be found in an output table (FAIL above 25 % unmatched)
  R2  the analysis script writes the report
  R3  (--rerun) re-running the script reproduces every CSV and the report
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

import pandas as pd

# ----------------------------------------------------------------------------------------
# Specification (keep in step with CLAUDE.md v3, sections 2 and 5)
# ----------------------------------------------------------------------------------------
ROLES = {"well_data_json", "well_data_dat", "mapping", "well_groups", "geometry_dat",
         "well_blocks_json", "reference_rock", "reference_bc", "reference_results",
         "reference_other", "out_of_scope"}
SINGLE_ROLES = ["geometry_dat", "well_blocks_json"]          # exactly one path each
REQUIRED_ROLES = SINGLE_ROLES + ["well_data_json", "mapping"]

REQUIRED_FILES = [
    "Scripts/natural_state_analysis.py",
    "Output/input_inventory.csv",
    "Output/NS_data_interpretation.md",
    "Output/NS_calibration_targets.csv",
    "Output/feedzones.csv",
    "Output/T_P_at_layer_centres.csv",
    "Output/well_summary.csv",
    "Output/pressure_consistency_check.csv",
    "Output/hydrostatic_check.csv",
    "Output/well_geometry_check.csv",
    "Output/run_console_log.txt",
    "Output/process_log.md",
]
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
}
HYDRO_BANDS = {"cased", "open"}
# Figure glob in Output/figures -> minimum count ("per_well" = one per well in the targets).
REQUIRED_FIGURES = {
    "well_map*.png": 1, "profile_*.png": "per_well", "PT_state*.png": 1,
    "maps_temperature*.png": 1, "maps_pressure*.png": 1, "section_*.png": 1,
    "conceptual_model*.png": 1,
}
REPORT_UNMATCHED_FAIL = 0.25
REPORT_UNMATCHED_WARN = 0.05

RESULTS = []


def record(cid, status, msg):
    RESULTS.append((cid, status, msg))


def is_true(v):
    return str(v).strip().lower() in ("true", "1", "yes")


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
    path = os.path.join(root, "Output", "input_inventory.csv")
    if not os.path.exists(path):
        record("I1", "FAIL", "Output/input_inventory.csv not found; geometry checks cannot run")
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


def read_layers(dat_path):
    """Layer section of a MULgraph grid file: 'index bottom_mRL centre_mRL'; layer 0 is the top."""
    layers, on = {}, False
    with open(dat_path) as f:
        for line in f:
            if line.startswith("LAYERS"):
                on = True
                continue
            if on:
                p = line.split()
                if len(p) != 3:
                    break
                layers[int(p[0])] = (float(p[1]), float(p[2]))
    if not layers:
        raise ValueError("no LAYERS section found in the geometry file")
    return {i: dict(top=layers[i - 1][0], bottom=layers[i][0], centre=layers[i][1])
            for i in sorted(layers) if i > 0}


def blocks_at(wells_json, well, z, tol=0.5):
    """Blocks whose interval contains z (both, if z is on a boundary)."""
    hits = set()
    for b in wells_json[well]["Blocks"]:
        top, bottom = max(b["interval"]), min(b["interval"])
        if bottom - tol <= z <= top + tol:
            hits.add(b["block"])
    return hits


# ----------------------------------------------------------------------------------------
# Checks
# ----------------------------------------------------------------------------------------
def check_files(root):
    missing = [p for p in REQUIRED_FILES if not os.path.exists(os.path.join(root, p))]
    record("F1", "FAIL" if missing else "PASS",
           f"missing files: {missing}" if missing else f"all {len(REQUIRED_FILES)} required files present")
    n_wells = 0
    tp = os.path.join(root, "Output", "NS_calibration_targets.csv")
    if os.path.exists(tp):
        n_wells = pd.read_csv(tp, usecols=["well"]).well.nunique()
    figdir = os.path.join(root, "Output", "figures")
    short = []
    for pat, n in REQUIRED_FIGURES.items():
        need = n_wells if n == "per_well" else n
        k = len(glob.glob(os.path.join(figdir, pat)))
        if k < need:
            short.append(f"{pat}: {k}/{need}")
    record("F1", "FAIL" if short else "PASS",
           f"figures missing: {short}" if short else "all required figure types present")


def check_columns(root):
    for name, cols in REQUIRED_COLUMNS.items():
        path = os.path.join(root, "Output", name)
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
    out = os.path.join(root, "Output")
    wells = set(pd.read_csv(os.path.join(out, "NS_calibration_targets.csv")).well)
    pc = pd.read_csv(os.path.join(out, "pressure_consistency_check.csv"))
    diff = wells.symmetric_difference(set(pc.well))
    dup = pc.well.duplicated().any()
    record("S2", "FAIL" if (diff or dup) else "PASS",
           f"pressure_consistency_check.csv must have one row per well; mismatched: {sorted(diff)}"
           f"{'; duplicate wells' if dup else ''}" if (diff or dup) else "one pressure-consistency row per well")
    hy = pd.read_csv(os.path.join(out, "hydrostatic_check.csv"))
    if "band" in hy.columns:
        bad = sorted(set(hy.band.dropna().astype(str)) - HYDRO_BANDS)
        record("S2", "WARN" if bad else "PASS",
               f"hydrostatic band labels other than {sorted(HYDRO_BANDS)}: {bad}" if bad
               else "hydrostatic band labels OK")


def check_geometry(root, roles):
    if not roles.get("geometry_dat") or not roles.get("well_blocks_json"):
        record("G*", "FAIL", "geometry checks skipped: inventory lacks geometry_dat or well_blocks_json")
        return
    wj = json.load(open(roles["well_blocks_json"][0]))
    layers = read_layers(roles["geometry_dat"][0])
    tgt = pd.read_csv(os.path.join(root, "Output", "NS_calibration_targets.csv"))

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
    fz = pd.read_csv(os.path.join(root, "Output", "feedzones.csv"))
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
    lc = pd.read_csv(os.path.join(root, "Output", "T_P_at_layer_centres.csv"))
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


def check_thermo(root):
    tgt = pd.read_csv(os.path.join(root, "Output", "NS_calibration_targets.csv"))
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


def check_weights(root):
    tgt = pd.read_csv(os.path.join(root, "Output", "NS_calibration_targets.csv"))
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


def check_report(root):
    out = os.path.join(root, "Output")
    rep = open(os.path.join(out, "NS_data_interpretation.md"), encoding="utf-8").read()
    rep = re.sub(r"§\s*\d+(\.\d+)*", " ", rep)                 # section references
    rep = re.sub(r"`[^`]*`", " ", rep)                          # file names and code
    rep = re.sub(r"\d{1,4}[/-]\d{1,2}[/-]\d{1,4}", " ", rep)    # dates
    rep = re.sub(r"\(.*?\.(png|csv|md|json|txt)\)", " ", rep)  # links
    pool = set()
    for path in glob.glob(os.path.join(out, "*")):
        if path.endswith((".csv", ".json", ".txt")) and not path.endswith("acceptance_check.txt"):
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
    src = open(os.path.join(root, "Scripts", "natural_state_analysis.py"), encoding="utf-8").read()
    ok = "NS_data_interpretation.md" in src
    record("R2", "PASS" if ok else "FAIL",
           "analysis script writes the report" if ok
           else "analysis script never names NS_data_interpretation.md: the report is not generated")


def sha(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def check_rerun(root, base):
    tmp = tempfile.mkdtemp(prefix="nscheck_")
    try:
        os.makedirs(os.path.join(tmp, "Output", "figures"))
        shutil.copytree(base, os.path.join(tmp, "Baseline"))
        shutil.copytree(os.path.join(root, "Scripts"), os.path.join(tmp, "Scripts"))
        p = subprocess.run([sys.executable, "Scripts/natural_state_analysis.py"], cwd=tmp,
                           capture_output=True, text=True, timeout=1800)
        if p.returncode != 0:
            record("R3", "FAIL", f"re-run exited {p.returncode}: {p.stderr.strip()[-300:]}")
            return
        diff = []
        for f in sorted(glob.glob(os.path.join(root, "Output", "*.csv"))) + \
                [os.path.join(root, "Output", "NS_data_interpretation.md")]:
            g = os.path.join(tmp, "Output", os.path.basename(f))
            if not os.path.exists(g):
                diff.append(os.path.basename(f) + " (not regenerated)")
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
    txt = open(os.path.join(root, "Output", "process_log.md"), encoding="utf-8").read().lower()
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
    check_files(root)
    for fn in (check_columns, check_tables, check_thermo, check_weights, check_report, check_process_log):
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
    out = os.path.join(root, "Output", "acceptance_check.txt")
    if os.path.isdir(os.path.dirname(out)):
        with open(out, "w", encoding="utf-8") as f:
            f.write(text + "\n")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
