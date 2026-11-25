"""Compare the natural-state analysis runs of Haiku, Sonnet and Opus across three iterations.

Run from the repository root:

    python Scripts/compare_model_runs.py            # extract runs, run checks, write tables + figures
    python Scripts/compare_model_runs.py --no-rerun # skip the (slower) re-run reproducibility check

Each run lives on its own git branch (origin/<Model><iteration>). The script extracts every branch
into a temporary folder, runs the iteration-3 acceptance check (Scripts/check_deliverables.py) on a
copy of each, pulls comparable numbers out of each run's own output tables, and writes:

    Output/model_comparison/*.csv          comparison tables
    Output/model_comparison/acceptance/    acceptance-check output for each run
    Output/figures/model_comparison/*.png  comparison figures

It only reads Baseline/ (for the model-surface ground truth) and never writes there.
The two written reports (Output/model_comparison_iteration3.md and
Output/model_comparison_across_iterations.md) quote numbers from these tables.
"""
import argparse
import importlib.util
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import matplotlib.image as mpimg  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "Output", "model_comparison")
FIG = os.path.join(ROOT, "Output", "figures", "model_comparison")
BASE = os.path.join(ROOT, "Baseline")
GEOM = os.path.join(BASE, "Teaching_Model_3D", "Models", "TE14515", "gTE14515.dat")
WELLS_JSON = os.path.join(BASE, "Teaching_Model_3D", "Models", "TE14515", "gTE14515_wells.json")

MODELS = ["Haiku", "Sonnet", "Opus"]
ITERS = [1, 2, 3]
RUNS = [(m, i) for i in ITERS for m in MODELS]
COLOURS = {"Haiku": "#2a9d8f", "Sonnet": "#e9a23b", "Opus": "#5b5bd6"}
WELLS = [f"WELL-{k}" for k in range(1, 11)]
CONF_SCORE = {"high": 3, "medium": 2, "low": 1}
CONF_COL = {"high": "#2e7d32", "medium": "#f2b134", "low": "#d1495b"}

plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.titleweight": "bold", "figure.dpi": 110})

# The iteration-3 acceptance check is the common yardstick.
_spec = importlib.util.spec_from_file_location("chk", os.path.join(ROOT, "Scripts", "check_deliverables.py"))
chk = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(chk)


def run_id(m, i):
    return f"{m}{i}"


# ---------------------------------------------------------------------------------------------
# 1. Extract every branch
# ---------------------------------------------------------------------------------------------
def git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True, check=True).stdout


def extract(tmp):
    meta = []
    for m, i in RUNS:
        rid, ref = run_id(m, i), f"origin/{run_id(m, i)}"
        dest = os.path.join(tmp, rid)
        os.makedirs(dest)
        data = git("archive", ref, "Output", "Scripts", "CLAUDE.md")
        tarfile.open(fileobj=io.BytesIO(data)).extractall(dest, filter="data")
        commit, date, subject = git("log", "-1", "--format=%h|%ad|%s", "--date=iso", ref).decode().strip().split("|", 2)
        status = git("diff-tree", "--no-commit-id", "-r", "--name-status", ref).decode().split("\n")
        touched = {l.split("\t")[-1] for l in status if l.strip() and l[0] in "AM"}
        present = [os.path.relpath(os.path.join(dp, f), dest).replace("\\", "/")
                   for dp, _, fs in os.walk(os.path.join(dest, "Output")) for f in fs]
        carried = [p for p in present if p not in touched]
        meta.append(dict(run=rid, model=m, iteration=i, commit=commit, commit_date=date[:16],
                         subject=subject, files_unchanged_from_previous_commit=";".join(carried)))
    return pd.DataFrame(meta)


# ---------------------------------------------------------------------------------------------
# 2. Acceptance check (iteration-3 version) on every run
# ---------------------------------------------------------------------------------------------
def acceptance(tmp, rerun):
    os.makedirs(os.path.join(OUT, "acceptance"), exist_ok=True)
    rows = []
    for m, i in RUNS:
        rid = run_id(m, i)
        work = os.path.join(tmp, "chk_" + rid)
        shutil.copytree(os.path.join(tmp, rid), work)
        cmd = [sys.executable, os.path.join(ROOT, "Scripts", "check_deliverables.py"), "--root", work,
               "--baseline", BASE] + (["--rerun"] if rerun else [])
        p = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        txt = p.stdout
        with open(os.path.join(OUT, "acceptance", f"acceptance_{rid}.txt"), "w", encoding="utf-8") as f:
            f.write(txt)
        for line in txt.splitlines():
            mm = re.match(r"^(PASS|WARN|FAIL)\s+(\S+)\s+(.*)$", line)
            if not mm:
                continue
            status, cid, msg = mm.groups()
            note = ""
            # git stores text with LF; a script that writes CRLF on Windows then differs only in line
            # endings ("0 lines differ"). That is an artefact of extracting from git, not of the run.
            if cid == "R3" and status == "FAIL" and "lines differ" in msg and \
                    all(x == "0" for x in re.findall(r"\((\d+) lines differ\)", msg)) and \
                    "not regenerated" not in msg:
                status, note = "PASS", "line-ending-only differences (git LF vs CRLF), treated as PASS"
            if cid.startswith("check_"):
                cid = {"check_tables": "S2", "check_thermo": "T1", "check_weights": "W1"}.get(cid, cid)
            rows.append(dict(run=rid, model=m, iteration=i, check=cid, status=status, message=msg[:400], note=note))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------
# 3. Pull comparable per-well numbers out of each run's own tables
# ---------------------------------------------------------------------------------------------
STD = ["Tmax_degC", "Tmax_elev_mRL", "T_at_TD_degC", "inversion_degC", "shallow_gradient_degC_per_km",
       "deep_gradient_degC_per_km", "z_100C_mRL", "z_200C_mRL", "min_boiling_margin_degC",
       "water_level_mRL", "pivot_P_bar", "regime", "ns_verdict", "ns_confidence"]
COLMAP = {
    "Sonnet1": dict(Tmax="Tmax_degC", z_Tmax="Tmax_elev_mRL", T_bottom="T_at_TD_degC",
                    inversion_drop_C="inversion_degC", z_T100="z_100C_mRL", z_T200="z_200C_mRL",
                    min_boiling_margin_C="min_boiling_margin_degC", water_level_est="water_level_mRL"),
    "Opus1": dict(Tmax_C="Tmax_degC", Tmax_elev_mRL="Tmax_elev_mRL", T_bottom_C="T_at_TD_degC",
                  deep_gradient_C_per_km="deep_gradient_degC_per_km", est_water_level_mRL="water_level_mRL",
                  min_boiling_margin_C="min_boiling_margin_degC", regime="regime"),
    "Haiku2": dict(max_T_degC="Tmax_degC", max_T_elev_mRL="Tmax_elev_mRL", T_at_TD_degC="T_at_TD_degC",
                   inversion_drop_K="inversion_degC", shallow_grad_K_per_km="shallow_gradient_degC_per_km",
                   deep_grad_K_per_km="deep_gradient_degC_per_km", T100_crossings_mRL="z_100C_mRL",
                   T200_crossings_mRL="z_200C_mRL", min_boiling_margin_K="min_boiling_margin_degC",
                   est_water_level_mRL="water_level_mRL", regime="regime", ns_verdict="ns_verdict",
                   ns_confidence="ns_confidence"),
    "Sonnet2": dict(T_max_degC="Tmax_degC", T_max_elev_mRL="Tmax_elev_mRL", T_at_TD_degC="T_at_TD_degC",
                    inversion_size_degC="inversion_degC", cap_gradient_K_km="shallow_gradient_degC_per_km",
                    deep_gradient_K_km="deep_gradient_degC_per_km", elev_100C_first_mRL="z_100C_mRL",
                    elev_200C_first_mRL="z_200C_mRL", min_boiling_margin_degC="min_boiling_margin_degC",
                    water_level_est_mRL="water_level_mRL", pivot_P_bar="pivot_P_bar", regime="regime"),
}


def first_number(v):
    if pd.isna(v):
        return np.nan
    mm = re.search(r"-?\d+(\.\d+)?", str(v))
    return float(mm.group(0)) if mm else np.nan


def regime_class(text):
    t = str(text).lower()
    if t in ("", "nan"):
        return ""
    for key, cls in [("two-phase", "two-phase"), ("upflow", "upflow"), ("outflow", "outflow"),
                     ("reversal", "outflow"), ("recharge", "recharge/downflow"), ("downflow", "recharge/downflow"),
                     ("cold", "recharge/downflow"), ("conductive", "conductive")]:
        if key in t:
            return cls
    return "other"


def well_metrics(tmp, rid):
    out = os.path.join(tmp, rid, "Output")
    if rid == "Haiku1":   # no well summary; per-well numbers live in a JSON file
        d = json.load(open(os.path.join(out, "NS_summary_stats.json")))
        df = pd.DataFrame([dict(well=w, Tmax_degC=d[w]["T_max"], min_boiling_margin_degC=d[w]["min_Tsat_minus_T"])
                           for w in WELLS])
    else:
        ws = pd.read_csv(os.path.join(out, "well_summary.csv"))
        if rid in COLMAP:
            ws = ws.rename(columns=COLMAP[rid])
        df = ws[["well"] + [c for c in STD if c in ws.columns]].copy()
        if rid == "Sonnet2":   # NS verdicts sit in the pressure-consistency table in this run
            pc = pd.read_csv(os.path.join(out, "pressure_consistency_check.csv"))
            df = df.merge(pc[["well", "verdict", "confidence_level"]].rename(
                columns={"verdict": "ns_verdict", "confidence_level": "ns_confidence"}), on="well", how="left")
    for c in STD:
        if c not in df.columns:
            df[c] = np.nan
    for c in STD[:11]:
        df[c] = df[c].map(first_number)
    df["regime_class"] = df["regime"].map(regime_class)
    df["ns_confidence"] = df["ns_confidence"].astype(str).str.strip().str.lower().replace({"nan": ""})
    df["ns_verdict"] = df["ns_verdict"].astype(str).str.strip().str.lower().replace({"nan": ""})
    df.insert(0, "run", rid)
    return df[["run", "well"] + STD + ["regime_class"]]


def model_surface_truth():
    lines = open(GEOM).read().splitlines()
    s, e = lines.index("SURFA"), lines.index("WELLS")
    surf = {l.split()[0]: float(l.split()[1]) for l in lines[s + 1:e] if l.strip()}
    wj = json.load(open(WELLS_JSON))
    rows = []
    for w in WELLS:
        b0 = wj[w]["Blocks"][0]
        col = re.match(r"[a-z]+", b0["block"]).group(0)
        rows.append(dict(well=w, column=col, wellhead_mRL=max(b0["interval"]), surface_mRL_SURFA=surf[col]))
    t = pd.DataFrame(rows)
    t["wellhead_minus_surface_m_truth"] = t.wellhead_mRL - t.surface_mRL_SURFA
    return t


def surface_by_run(tmp, rid):
    out = os.path.join(tmp, rid, "Output")
    f = os.path.join(out, "well_geometry_check.csv")
    g = pd.read_csv(f) if os.path.exists(f) else pd.DataFrame()
    if "model_surface_mRL" in g.columns:
        return g[["well", "model_surface_mRL"]]
    ws_f = os.path.join(out, "well_summary.csv")
    if os.path.exists(ws_f):
        ws = pd.read_csv(ws_f)
        for c in ["column_surface_mRL", "wellhead_surface_SURFA_mRL"]:
            if c in ws.columns:
                return ws[["well", c]].rename(columns={c: "model_surface_mRL"})
    return pd.DataFrame(columns=["well", "model_surface_mRL"])


def point_checks(tmp, rid, wj):
    """Survey-point block agreement with the well-to-block file and Tsat accuracy (independent IF97)."""
    out = os.path.join(tmp, rid, "Output")
    res = dict(run=rid, n_points=np.nan, block_checked=np.nan, block_mismatch=np.nan,
               tsat_checked=np.nan, tsat_max_err_K=np.nan)
    f = os.path.join(out, "NS_calibration_targets.csv")
    if not os.path.exists(f):
        return res
    t = pd.read_csv(f).rename(columns={"elev_mRL": "elevation_mRL", "Tsat_C": "Tsat_degC"})
    if "elevation_mRL" not in t.columns:   # Sonnet1: one row per block, not per survey point
        return res
    res["n_points"] = len(t)
    if "block" in t.columns:
        n = bad = 0
        for _, r in t.iterrows():
            if chk.is_true(r.get("above_model_surface", False)) or r.well not in wj:
                continue
            ok = chk.blocks_at(wj, r.well, r.elevation_mRL)
            if ok:
                n += 1
                bad += r.block not in ok
        res.update(block_checked=n, block_mismatch=bad)
    if {"P_bar", "Tsat_degC"} <= set(t.columns):
        e = [abs(chk.tsat_c(p) - ts) for p, ts in zip(t.P_bar, t.Tsat_degC)
             if pd.notna(ts) and not np.isnan(chk.tsat_c(p))]
        res.update(tsat_checked=len(e), tsat_max_err_K=max(e) if e else np.nan)
    return res


# ---------------------------------------------------------------------------------------------
# 4. Size / effort / report metrics
# ---------------------------------------------------------------------------------------------
def report_metrics(tmp, rid):
    root = os.path.join(tmp, rid)
    out = os.path.join(root, "Output")
    rep = open(os.path.join(out, "NS_data_interpretation.md"), encoding="utf-8").read()
    words = len(re.findall(r"\S+", rep))
    secs = re.split(r"\n(?=## )", rep)
    s1 = next((s for s in secs if re.match(r"## *1\b|## *Key|## *1\.", s.strip())), "")
    n_concl = len(re.findall(r"^\s*\d+\.\s", s1, flags=re.M))
    n_obs = len(re.findall(r"\[(obs|observation)\]", rep, flags=re.I))
    n_inf = len(re.findall(r"\[(inf|interp|interpretation)\]", rep, flags=re.I))
    n_action = len(re.findall(r"\*\*(Action|Calibration)[:\*]", rep))
    n_ref = len(re.findall(r"\[?reference model\]?", rep, flags=re.I))
    script = os.path.join(root, "Scripts", "natural_state_analysis.py")
    sl = open(script, encoding="utf-8").read().splitlines() if os.path.exists(script) else []
    plog = os.path.join(out, "process_log.md")
    pw = len(re.findall(r"\S+", open(plog, encoding="utf-8").read())) if os.path.exists(plog) else 0
    files = [f for f in os.listdir(out) if os.path.isfile(os.path.join(out, f))]
    figs = [f for f in os.listdir(os.path.join(out, "figures"))] if os.path.isdir(os.path.join(out, "figures")) else []
    tables = [f for f in files if f.endswith((".csv", ".json"))]
    kell = any("kell" in open(os.path.join(out, f), encoding="utf-8", errors="ignore").readline().lower()
               for f in tables if f.endswith(".csv")) or \
        ("kell" in open(os.path.join(out, "NS_summary_stats.json")).read().lower()
         if os.path.exists(os.path.join(out, "NS_summary_stats.json")) else False)
    return dict(run=rid, report_words=words, report_h2_sections=len(re.findall(r"^## ", rep, flags=re.M)),
                key_conclusions=n_concl, obs_labels=n_obs, interp_labels=n_inf, action_lines=n_action,
                reference_model_mentions=n_ref, script_lines=len(sl), process_log_words=pw,
                output_tables=len(tables), output_files=len(files), figures=len(figs),
                kell_density_in_outputs=kell)


def coverage(tmp, rid):
    """Share of the iteration-3 standard deliverables (files, required columns, figure types) present."""
    root = os.path.join(tmp, rid)
    f_ok = sum(os.path.exists(os.path.join(root, p)) for p in chk.REQUIRED_FILES)
    c_tot = c_ok = 0
    for name, cols in chk.REQUIRED_COLUMNS.items():
        p = os.path.join(root, "Output", name)
        have = list(pd.read_csv(p, nrows=0).columns) if os.path.exists(p) else []
        c_tot += len(cols)
        c_ok += sum(c in have for c in cols)
    import glob
    g_ok = sum(len(glob.glob(os.path.join(root, "Output", "figures", g))) > 0 for g in chk.REQUIRED_FIGURES)
    return dict(run=rid, files_present=f_ok, files_required=len(chk.REQUIRED_FILES),
                columns_present=c_ok, columns_required=c_tot,
                figure_types_present=g_ok, figure_types_required=len(chk.REQUIRED_FIGURES),
                pct_files=round(100 * f_ok / len(chk.REQUIRED_FILES), 1),
                pct_columns=round(100 * c_ok / c_tot, 1),
                pct_figure_types=round(100 * g_ok / len(chk.REQUIRED_FIGURES), 1))


# ---------------------------------------------------------------------------------------------
# 5. Iteration-3 specific tables
# ---------------------------------------------------------------------------------------------
def iter3_tables(tmp):
    w, h, g, fz = [], [], [], []
    for m in MODELS:
        out = os.path.join(tmp, run_id(m, 3), "Output")
        t = pd.read_csv(os.path.join(out, "NS_calibration_targets.csv"))
        t["model"] = m
        w.append(t[["model", "well", "elevation_mRL", "weight", "weight_reason", "weight_P", "ns_confidence",
                    "behind_casing", "in_feedzone"]])
        hh = pd.read_csv(os.path.join(out, "hydrostatic_check.csv"))
        hh["model"] = m
        h.append(hh[["model", "well", "z_top_mRL", "z_bottom_mRL", "length_m", "band", "ratio_measured_to_IF97",
                     "rho_implied_kgm3"]])
        pc = pd.read_csv(os.path.join(out, "pressure_consistency_check.csv"))
        pc["model"] = m
        g.append(pc[["model", "well", "gauge_elev_mRL", "survey_P_at_gauge_bar", "first_gauge_P_bar",
                     "diff_survey_minus_first_bar", "verdict", "confidence"]])
        f = pd.read_csv(os.path.join(out, "feedzones.csv"))
        f["model"] = m
        fz.append(f[["model", "well", "feedzone", "quality", "T_degC", "P_bar", "boiling_margin_degC",
                     "block", "weight_T", "weight_P"]])
    return pd.concat(w), pd.concat(h), pd.concat(g), pd.concat(fz)


# ---------------------------------------------------------------------------------------------
# 6. Figures
# ---------------------------------------------------------------------------------------------
def savefig(fig, name):
    fig.savefig(os.path.join(FIG, name), bbox_inches="tight", dpi=130)
    plt.close(fig)


def fig_acceptance(acc):
    checks = ["I1", "F1", "S1", "S2", "G*", "G1", "G2", "G3", "T1", "T2", "W1", "W2", "R1", "R2", "R3", "P1"]
    order = {"FAIL": 0, "WARN": 1, "PASS": 2}
    runs = [run_id(m, i) for m, i in RUNS]
    mat = np.full((len(runs), len(checks)), np.nan)
    for (r, c), grp in acc.groupby(["run", "check"]):
        if c in checks:
            mat[runs.index(r), checks.index(c)] = min(order[s] for s in grp.status)
    # G* (geometry skipped) stands in for G1-G3 when the inventory is missing
    for k, r in enumerate(runs):
        if not np.isnan(mat[k, checks.index("G*")]):
            for c in ["G1", "G2", "G3"]:
                mat[k, checks.index(c)] = 0
    mat = np.delete(mat, checks.index("G*"), axis=1)
    checks.remove("G*")
    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#d1495b", "#f2b134", "#2e7d32"])
    fig, ax = plt.subplots(figsize=(11.5, 4.8))
    ax.imshow(np.nan_to_num(mat, nan=-1), cmap=ListedColormap(["#dddddd", "#d1495b", "#f2b134", "#2e7d32"]),
              vmin=-1, vmax=2, aspect="auto")
    labels = {"I1": "Input\ninventory", "F1": "Files &\nfigures", "S1": "Column\nnames", "S2": "Table\nshape",
              "G1": "Survey\nblocks", "G2": "Feed-zone\nblocks", "G3": "Layer\ncentres", "T1": "Tsat\n(IF97)",
              "T2": "Boiling\nmargin", "W1": "Weight\nrules", "W2": "Confidence\nvaries", "R1": "Report\nnumbers\nin tables",
              "R2": "Script\nwrites\nreport", "R3": "Re-run\nrepro-\nduces", "P1": "Process\nlog"}
    ax.set_xticks(range(len(checks)), [labels[c] for c in checks], fontsize=7)
    ax.set_yticks(range(len(runs)), [f"{r[:-1]}  (iter {r[-1]})" for r in runs])
    for k in range(len(runs)):
        for j in range(len(checks)):
            v = mat[k, j]
            ax.text(j, k, {0: "FAIL", 1: "WARN", 2: "PASS"}.get(v, "–") if not np.isnan(v) else "n/a",
                    ha="center", va="center", fontsize=6.5, color="white" if v in (0, 2) else "black")
    for y in (2.5, 5.5):
        ax.axhline(y, color="white", lw=3)
    ax.set_title("Iteration-3 acceptance check applied to all nine runs")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    savefig(fig, "x_acceptance_heatmap.png")
    _ = cmap


def fig_scorecard_trends(rep, cov, acc):
    df = rep.merge(cov, on="run")
    df["model"] = df.run.str[:-1]
    df["iteration"] = df.run.str[-1].astype(int)
    passes = acc.groupby("run").status.apply(lambda s: (s == "PASS").sum() / len(s) * 100).rename("pct_pass")
    df = df.merge(passes, on="run")
    panels = [("pct_pass", "Acceptance checks passed (%)"), ("pct_columns", "Required columns present (%)"),
              ("script_lines", "Analysis script (lines)"), ("report_words", "Main report (words)"),
              ("output_tables", "Output tables (CSV/JSON)"), ("figures", "Figures produced")]
    fig, axs = plt.subplots(2, 3, figsize=(11, 6))
    for ax, (col, title) in zip(axs.flat, panels):
        for m in MODELS:
            d = df[df.model == m].sort_values("iteration")
            ax.plot(d.iteration, d[col], "-o", color=COLOURS[m], label=m, lw=2)
            for _, r in d.iterrows():
                ax.annotate(f"{r[col]:.0f}", (r.iteration, r[col]), textcoords="offset points", xytext=(5, 3),
                            fontsize=7, color=COLOURS[m])
        ax.set_xticks(ITERS, ["Iter 1", "Iter 2", "Iter 3"])
        ax.set_title(title, fontsize=9)
        ax.grid(axis="y", alpha=.3)
    axs[0, 0].legend(frameon=False)
    fig.suptitle("How each model's output changed across the three iterations", fontweight="bold")
    fig.tight_layout()
    savefig(fig, "x_trends.png")
    return df


def fig_agreement(wm):
    """How often the three models give the same number for the same well, per iteration."""
    metrics = [("Tmax_degC", "Max temperature", 0.15, "°C"), ("T_at_TD_degC", "T at total depth", 0.15, "°C"),
               ("inversion_degC", "Inversion size", 0.15, "°C"), ("z_100C_mRL", "100 °C isotherm", 1.0, "m"),
               ("z_200C_mRL", "200 °C isotherm", 1.0, "m"),
               ("min_boiling_margin_degC", "Min boiling margin", 0.15, "°C"),
               ("shallow_gradient_degC_per_km", "Shallow gradient", 1.0, "°C/km"),
               ("deep_gradient_degC_per_km", "Deep gradient", 1.0, "°C/km"),
               ("water_level_mRL", "Water level", 1.0, "m"), ("pivot_P_bar", "Pivot pressure", 0.15, "bar"),
               ("regime_class", "Regime", None, ""), ("ns_confidence", "NS confidence", None, "")]
    rows = []
    for it in ITERS:
        sub = wm[wm.run.str.endswith(str(it))]
        for col, lab, tol, unit in metrics:
            piv = sub.pivot(index="well", columns="run", values=col)
            if tol is None:
                piv = piv.replace("", np.nan)
            n_models = int(piv.notna().any().sum())
            piv = piv.loc[:, piv.notna().any()]
            full = piv.dropna(how="any") if n_models >= 2 else piv.iloc[0:0]
            if tol is None:
                agree = full.nunique(axis=1) == 1
                maxd = np.nan
            else:
                rng = full.max(axis=1) - full.min(axis=1)
                agree = rng <= tol
                maxd = rng.max() if len(rng) else np.nan
            rows.append(dict(iteration=it, metric=col, label=lab, unit=unit, tolerance=tol,
                             models_reporting=n_models, wells_compared=int(len(full)),
                             wells_all_agree=int(agree.sum()), max_difference=maxd))
    ag = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(7.4, 6.2))
    labs = [m[1] for m in metrics]
    for j, it in enumerate(ITERS):
        for k, (col, lab, tol, unit) in enumerate(metrics):
            r = ag[(ag.iteration == it) & (ag.metric == col)].iloc[0]
            if r.models_reporting < 2 or r.wells_compared == 0:
                colr, txt = "#e6e6e6", ("reported by 1 model" if r.models_reporting == 1 else "not reported")
            else:
                frac = r.wells_all_agree / r.wells_compared
                colr = plt.cm.RdYlGn(0.15 + 0.75 * frac)
                txt = f"{r.wells_all_agree}/{r.wells_compared} wells agree"
                if r.models_reporting < 3:
                    txt += f"\n({r.models_reporting} models)"
                elif pd.notna(r.max_difference) and r.wells_all_agree < r.wells_compared:
                    txt += f"\nmax diff {r.max_difference:.1f} {unit}"
            ax.add_patch(plt.Rectangle((j, k), 0.96, 0.92, color=colr))
            ax.text(j + 0.48, k + 0.46, txt, ha="center", va="center", fontsize=7)
    ax.set_xlim(0, 3)
    ax.set_ylim(len(metrics), 0)
    ax.set_xticks([0.48, 1.48, 2.48], ["Iteration 1", "Iteration 2", "Iteration 3"])
    ax.xaxis.tick_top()
    ax.set_yticks([k + 0.46 for k in range(len(metrics))], labs)
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title("Do the three models give the same answer for the same well?\n"
                 "(wells where every reporting model agrees within rounding)", fontsize=9.5, pad=24)
    savefig(fig, "x_agreement.png")
    return ag


def fig_confidence_grid(wm, runs, name, title):
    fig, ax = plt.subplots(figsize=(1.0 + 0.62 * len(runs), 4.4))
    for j, rid in enumerate(runs):
        d = wm[wm.run == rid].set_index("well")
        for k, w in enumerate(WELLS):
            conf = d.ns_confidence.get(w, "") if w in d.index else ""
            ver = d.ns_verdict.get(w, "") if w in d.index else ""
            col = CONF_COL.get(conf, "#e6e6e6")
            ax.add_patch(plt.Rectangle((j, k), 0.95, 0.9, color=col))
            short = {"natural state": "NS", "probably natural state": "prob. NS",
                     "possibly disturbed": "disturbed?"}.get(ver, "–" if not ver else ver[:9])
            ax.text(j + 0.475, k + 0.45, short, ha="center", va="center", fontsize=6.5,
                    color="white" if conf in ("high", "low") else "black")
    ax.set_xlim(0, len(runs))
    ax.set_ylim(len(WELLS), 0)
    ax.set_xticks([j + 0.475 for j in range(len(runs))], [f"{r[:-1]}\niter {r[-1]}" for r in runs], fontsize=7.5)
    ax.set_yticks([k + 0.45 for k in range(len(WELLS))], WELLS)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=CONF_COL[c], label=c) for c in ["high", "medium", "low"]] +
              [Patch(color="#e6e6e6", label="not given")], loc="upper center", bbox_to_anchor=(0.5, -0.1),
              ncol=4, frameon=False, fontsize=7)
    ax.set_title(title, fontsize=9.5)
    savefig(fig, name)


def fig_i3_metrics(wm, surf):
    i3 = wm[wm.run.str.endswith("3")].copy()
    i3["model"] = i3.run.str[:-1]
    metrics = [("shallow_gradient_degC_per_km", "Shallow gradient (°C/km)"),
               ("deep_gradient_degC_per_km", "Deep gradient (°C/km)"),
               ("water_level_mRL", "Water level (mRL)"),
               ("min_boiling_margin_degC", "Min boiling margin (°C)")]
    fig, axs = plt.subplots(1, 5, figsize=(14, 4.6), sharey=True)
    off = {"Haiku": -0.22, "Sonnet": 0, "Opus": 0.22}
    for ax, (col, lab) in zip(axs, metrics):
        for m in MODELS:
            d = i3[i3.model == m].set_index("well").reindex(WELLS)
            ax.scatter(d[col], np.arange(len(WELLS)) + off[m], color=COLOURS[m], s=28, label=m, zorder=3)
        ax.set_title(lab, fontsize=9)
        ax.grid(axis="x", alpha=.3)
        for k in range(len(WELLS)):
            ax.axhline(k + 0.5, color="#eee", lw=.6)
    ax = axs[4]
    for m in MODELS:
        d = surf[surf.run == run_id(m, 3)].set_index("well").reindex(WELLS)
        ax.scatter(d.wellhead_minus_surface_m, np.arange(len(WELLS)) + off[m], color=COLOURS[m], s=28, zorder=3)
    t = surf.drop_duplicates("well").set_index("well").reindex(WELLS)
    ax.scatter(t.wellhead_minus_surface_m_truth, np.arange(len(WELLS)), marker="|", s=260, color="k",
               label="grid file (SURFA)", zorder=4)
    ax.set_title("Wellhead minus model surface (m)", fontsize=9)
    ax.grid(axis="x", alpha=.3)
    axs[0].set_yticks(range(len(WELLS)), WELLS)
    axs[0].invert_yaxis()
    h1, l1 = axs[0].get_legend_handles_labels()
    h2, l2 = ax.get_legend_handles_labels()
    fig.legend(h1 + h2, l1 + l2, loc="upper center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.95))
    fig.suptitle("Iteration 3: the well-summary numbers on which the three models differ "
                 "(Tmax, T at TD, inversions, isotherms and pivot pressures are identical)", fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    savefig(fig, "i3_well_metric_differences.png")


def fig_i3_weights(w3, h3):
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.2))
    ax = axs[0]
    piv = w3.groupby(["well", "model"]).weight.mean().unstack().reindex(index=WELLS, columns=MODELS)
    x = np.arange(len(WELLS))
    for k, m in enumerate(MODELS):
        ax.bar(x + (k - 1) * 0.27, piv[m], width=0.27, color=COLOURS[m], label=m)
    ax.set_xticks(x, [w.replace("WELL-", "W") for w in WELLS])
    ax.set_title("Mean temperature weight per well", fontsize=9)
    ax.legend(frameon=False, fontsize=7)
    ax = axs[1]
    tot = w3.groupby("model").agg(T_weight_sum=("weight", "sum"), n=("weight", "size"),
                                  P_points=("weight_P", lambda s: (s > 0).sum()),
                                  zero_T=("weight", lambda s: (s == 0).sum())).reindex(MODELS)
    ax.bar(MODELS, tot.T_weight_sum, color=[COLOURS[m] for m in MODELS])
    for k, (m, r) in enumerate(tot.iterrows()):
        ax.text(k, r.T_weight_sum, f"{r.T_weight_sum:.1f}\n{int(r.P_points)} P-targets\n{int(r.zero_T)} zero-weight",
                ha="center", va="bottom", fontsize=7)
    ax.margins(y=0.3)
    ax.set_title(f"Total temperature weight ({int(tot.n.iloc[0])} survey points)", fontsize=9)
    ax = axs[2]
    for k, m in enumerate(MODELS):
        d = h3[h3.model == m]
        for band, mk in (("cased", "o"), ("open", "s")):
            dd = d[d.band == band]
            ax.scatter(np.full(len(dd), k) + (0.12 if band == "open" else -0.12), dd.ratio_measured_to_IF97,
                       marker=mk, color=COLOURS[m], alpha=.75, s=22,
                       label=f"{band}" if k == 0 else None, edgecolor="k" if band == "open" else "none")
    ax.axhline(1, color="k", lw=.8, ls="--")
    ax.set_xticks(range(3), [f"{m}\n({(h3.model == m).sum()} windows)" for m in MODELS])
    ax.set_ylabel("measured ΔP / IF97 ΔP")
    ax.set_title("Hydrostatic windows (circles cased, squares open)", fontsize=9)
    fig.suptitle("Iteration 3: calibration weights and hydrostatic tests", fontweight="bold")
    fig.tight_layout()
    savefig(fig, "i3_weights_hydrostatic.png")
    return tot


def fig_montage(tmp, runs, rel_candidates, name, title, ncols=3):
    n = len(runs)
    nrows = int(np.ceil(n / ncols))
    fig, axs = plt.subplots(nrows, ncols, figsize=(5.2 * ncols, 4.0 * nrows))
    for ax, rid in zip(np.atleast_1d(axs).flat, runs):
        ax.axis("off")
        figdir = os.path.join(tmp, rid, "Output", "figures")
        path = None
        for pat in rel_candidates:
            hits = sorted(f for f in os.listdir(figdir) if re.fullmatch(pat, f)) if os.path.isdir(figdir) else []
            if hits:
                path = os.path.join(figdir, hits[0])
                break
        if path:
            ax.imshow(mpimg.imread(path))
            ax.set_title(f"{rid[:-1]}, iteration {rid[-1]}  ({os.path.basename(path)})", fontsize=8)
        else:
            ax.text(0.5, 0.5, f"{rid[:-1]}, iteration {rid[-1]}\nno such figure", ha="center", va="center")
    for ax in list(np.atleast_1d(axs).flat)[n:]:
        ax.axis("off")
    fig.suptitle(title, fontweight="bold")
    fig.tight_layout()
    savefig(fig, name)


def fig_i3_scorecard(rep, extra):
    i3 = rep[rep.run.str.endswith("3")].set_index("run").reindex([run_id(m, 3) for m in MODELS])
    cols = [("report_words", "Report words"), ("key_conclusions", "Key conclusions"),
            ("obs_labels", "[obs] labels"), ("interp_labels", "[interp] labels"),
            ("script_lines", "Script lines"), ("output_tables", "Output tables"), ("figures", "Figures"),
            ("process_log_words", "Process-log words")]
    fig, axs = plt.subplots(1, len(cols), figsize=(14, 2.8))
    for ax, (c, lab) in zip(axs, cols):
        ax.bar(range(3), i3[c], color=[COLOURS[m] for m in MODELS])
        for k, v in enumerate(i3[c]):
            ax.text(k, v, f"{v:.0f}", ha="center", va="bottom", fontsize=7)
        ax.set_xticks(range(3), [m[0] for m in MODELS])
        ax.set_title(lab, fontsize=8.5)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.margins(y=0.2)
    fig.suptitle("Iteration 3 at a glance (H = Haiku, S = Sonnet, O = Opus)", fontweight="bold")
    fig.tight_layout()
    savefig(fig, "i3_scorecard.png")
    _ = extra


# ---------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-rerun", action="store_true")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    os.makedirs(FIG, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="cmp_runs_")
    try:
        meta = extract(tmp)
        meta.to_csv(os.path.join(OUT, "runs.csv"), index=False)
        print("extracted", len(meta), "runs")

        acc = acceptance(tmp, rerun=not a.no_rerun)
        acc.to_csv(os.path.join(OUT, "acceptance_results.csv"), index=False)
        accsum = acc.groupby("run").status.value_counts().unstack(fill_value=0).reindex(columns=["PASS", "WARN", "FAIL"],
                                                                                         fill_value=0)
        accsum = accsum.reindex([run_id(m, i) for m, i in RUNS])
        accsum.to_csv(os.path.join(OUT, "acceptance_summary.csv"))
        print(accsum)

        wm = pd.concat([well_metrics(tmp, run_id(m, i)) for m, i in RUNS], ignore_index=True)
        wm.to_csv(os.path.join(OUT, "well_metrics_all_runs.csv"), index=False)

        truth = model_surface_truth()
        surf = []
        for m, i in RUNS:
            s = surface_by_run(tmp, run_id(m, i)).merge(truth, on="well", how="right")
            s["run"] = run_id(m, i)
            s["wellhead_minus_surface_m"] = s.wellhead_mRL - s.model_surface_mRL
            s["surface_correct"] = np.where(s.model_surface_mRL.isna(), np.nan,
                                            (s.model_surface_mRL - s.surface_mRL_SURFA).abs() < 0.5)
            surf.append(s)
        surf = pd.concat(surf, ignore_index=True)[["run", "well", "column", "wellhead_mRL", "model_surface_mRL",
                                                   "surface_mRL_SURFA", "wellhead_minus_surface_m",
                                                   "wellhead_minus_surface_m_truth", "surface_correct"]]
        surf.to_csv(os.path.join(OUT, "model_surface_check.csv"), index=False)

        wj = json.load(open(WELLS_JSON))
        pts = pd.DataFrame([point_checks(tmp, run_id(m, i), wj) for m, i in RUNS])
        rep = pd.DataFrame([report_metrics(tmp, run_id(m, i)) for m, i in RUNS])
        cov = pd.DataFrame([coverage(tmp, run_id(m, i)) for m, i in RUNS])
        sc = surf.groupby("run").surface_correct.agg(lambda s: f"{int(np.nansum(s))}/{int(s.notna().sum())}"
                                                       if s.notna().any() else "n/a").rename("model_surface_correct")
        overview = meta[["run", "model", "iteration", "commit", "commit_date"]].merge(accsum.reset_index(), on="run") \
            .merge(cov, on="run").merge(rep, on="run").merge(pts, on="run").merge(sc, on="run")
        overview.to_csv(os.path.join(OUT, "run_overview.csv"), index=False)

        # NS verdicts side by side
        ns = wm.pivot(index="well", columns="run", values="ns_confidence").reindex(index=WELLS)
        ns = ns[[run_id(m, i) for m, i in RUNS]]
        ns.to_csv(os.path.join(OUT, "ns_confidence_by_run.csv"))
        rg = wm.pivot(index="well", columns="run", values="regime_class").reindex(index=WELLS)
        rg[[run_id(m, i) for m, i in RUNS]].to_csv(os.path.join(OUT, "regime_by_run.csv"))

        # iteration-3 detail
        w3, h3, g3, f3 = iter3_tables(tmp)
        i3 = wm[wm.run.str.endswith("3")]
        side = i3.melt(id_vars=["run", "well"], value_vars=STD + ["regime_class"], var_name="metric")
        side = side.pivot_table(index=["well", "metric"], columns="run", values="value", aggfunc="first").reset_index()
        side = side[["well", "metric"] + [run_id(m, 3) for m in MODELS]]

        def differ(r):
            v = [r[run_id(m, 3)] for m in MODELS]
            try:
                f = [float(x) for x in v if pd.notna(x) and str(x) != ""]
                return (max(f) - min(f) > 0.15) if len(f) > 1 else (len(f) != 3 and len(f) > 0)
            except (TypeError, ValueError):
                return len({str(x) for x in v}) > 1
        side["models_differ"] = side.apply(differ, axis=1)
        side.to_csv(os.path.join(OUT, "iter3_well_summary_side_by_side.csv"), index=False)
        side.groupby("metric").models_differ.sum().rename("wells_where_models_differ") \
            .to_csv(os.path.join(OUT, "iter3_metric_disagreement_counts.csv"))
        wsum = w3.groupby(["model", "well"]).agg(n_points=("weight", "size"), mean_T_weight=("weight", "mean"),
                                                 sum_T_weight=("weight", "sum"),
                                                 n_P_targets=("weight_P", lambda s: (s > 0).sum()),
                                                 ns_confidence=("ns_confidence", "first")).reset_index()
        wsum.to_csv(os.path.join(OUT, "iter3_weights_by_well.csv"), index=False)
        reasons = w3.assign(reason=w3.weight_reason.astype(str).str.split(r"[;,(]").str[0].str.strip().str.lower())
        reasons.groupby(["model", "reason"]).size().unstack(0, fill_value=0) \
            .to_csv(os.path.join(OUT, "iter3_weight_reasons.csv"))
        hsum = h3.groupby(["model", "band"]).ratio_measured_to_IF97.agg(["count", "median", "min", "max"]).round(4)
        hsum.to_csv(os.path.join(OUT, "iter3_hydrostatic_by_band.csv"))
        g3.to_csv(os.path.join(OUT, "iter3_pressure_gauge_test.csv"), index=False)
        f3w = f3.pivot_table(index=["well", "feedzone"], columns="model", values=["weight_P", "weight_T"]).round(3)
        f3w.to_csv(os.path.join(OUT, "iter3_feedzone_weights.csv"))
        # weight per point: how often do the models give the same weight?
        wp = w3.pivot_table(index=["well", "elevation_mRL"], columns="model", values="weight")
        same = pd.DataFrame(dict(
            pair=["Haiku-Sonnet", "Haiku-Opus", "Sonnet-Opus", "all three"],
            points_same_T_weight=[int((abs(wp.Haiku - wp.Sonnet) < 1e-6).sum()),
                                  int((abs(wp.Haiku - wp.Opus) < 1e-6).sum()),
                                  int((abs(wp.Sonnet - wp.Opus) < 1e-6).sum()),
                                  int(((abs(wp.Haiku - wp.Sonnet) < 1e-6) & (abs(wp.Sonnet - wp.Opus) < 1e-6)).sum())],
            points=len(wp)))
        same.to_csv(os.path.join(OUT, "iter3_point_weight_agreement.csv"), index=False)

        # figures
        fig_acceptance(acc)
        trends = fig_scorecard_trends(rep, cov, acc)
        trends.to_csv(os.path.join(OUT, "trend_metrics.csv"), index=False)
        ag = fig_agreement(wm)
        ag.round(3).to_csv(os.path.join(OUT, "model_agreement_by_iteration.csv"), index=False)
        fig_confidence_grid(wm, [run_id(m, 3) for m in MODELS], "i3_ns_confidence.png",
                            "Iteration 3: natural-state verdict and confidence")
        fig_confidence_grid(wm, [run_id(m, i) for i in (2, 3) for m in MODELS], "x_ns_confidence.png",
                            "Natural-state verdict and confidence (iteration 1 gave none)")
        fig_i3_metrics(wm, surf)
        tot = fig_i3_weights(w3, h3)
        tot.to_csv(os.path.join(OUT, "iter3_weight_totals.csv"))
        fig_i3_scorecard(rep, None)
        fig_montage(tmp, [run_id(m, 3) for m in MODELS], [r"conceptual_model.*\.png"], "i3_conceptual_models.png",
                    "Iteration 3: the three conceptual-model sketches", ncols=3)
        fig_montage(tmp, [run_id(m, 3) for m in MODELS], [r"profile_WELL-3\.png"], "i3_profile_WELL-3.png",
                    "Iteration 3: the same well (WELL-3) as drawn by each model", ncols=3)
        fig_montage(tmp, [run_id(m, 3) for m in MODELS], [r"section_T.*\.png", r"section_.*\.png"],
                    "i3_sections.png", "Iteration 3: temperature cross-sections through the hottest wells", ncols=3)
        fig_montage(tmp, [run_id(m, i) for m, i in RUNS],
                    [r"well_map.*\.png", r"01_well_map\.png"], "x_well_maps.png",
                    "The well map from every run (rows: iteration 1, 2, 3; columns: Haiku, Sonnet, Opus)")
        fig_montage(tmp, [run_id(m, i) for m, i in RUNS],
                    [r"conceptual_model.*\.png", r"section_hot.*\.png", r"06_section_A.*\.png", r"05_T_P_maps.*\.png"],
                    "x_conceptual_models.png",
                    "Conceptual-model sketch from every run (fallback: a section or map where no sketch exists)")
        print("tables in", os.path.relpath(OUT, ROOT), "; figures in", os.path.relpath(FIG, ROOT))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
