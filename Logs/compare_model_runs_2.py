"""Compare the second-iteration natural-state analyses made by three models.

The runs live on the git branches origin/Opus2, origin/Sonnet2 and origin/Haiku2.
Each branch tip holds that model's Scripts/natural_state_analysis.py and Output/.

This script:
  1. extracts each branch's Scripts/ and Output/ to a temporary directory (git archive);
  2. re-runs each model's script against a temporary copy of Baseline/ and checks that
     the committed outputs are reproduced (Baseline/ itself is only read);
  3. compares deliverables, per-well characteristics, verdicts, weights, block
     assignment and layer centres between the runs;
  4. writes the tables to Output/model_comparison/ and the console output to
     Output/model_comparison/compare_console_log.txt.

Run from the repository root:  python Scripts/compare_model_runs.py  [--no-rerun]
"""
import hashlib
import io
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time

import pandas as pd

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
OUT = os.path.join(ROOT, "Output", "model_comparison")
MODELS = ["Opus2", "Sonnet2", "Haiku2"]
REQUIRED_FILES = [
    "Scripts/natural_state_analysis.py",
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
REQUIRED_TARGET_COLS = ["well", "x", "y", "elevation_mRL", "T_degC", "P_bar", "Tsat_degC",
                        "boiling_margin_degC", "column", "layer", "block", "in_feedzone",
                        "behind_casing", "above_model_surface", "ns_confidence", "weight",
                        "weight_reason"]
# Column names differ between the runs; map them onto common names.
SUMMARY_COLS = {
    "Opus2": dict(Tmax="Tmax_degC", Tmax_elev="Tmax_elev_mRL", T_TD="T_at_TD_degC", inversion="inversion_degC",
                  shallow_grad="shallow_gradient_degC_per_km", deep_grad="deep_gradient_degC_per_km",
                  z100="z_100C_mRL", z200="z_200C_mRL", min_bm="min_boiling_margin_degC",
                  water_level="water_level_mRL", regime="regime"),
    "Sonnet2": dict(Tmax="T_max_degC", Tmax_elev="T_max_elev_mRL", T_TD="T_at_TD_degC", inversion="inversion_size_degC",
                    shallow_grad="cap_gradient_K_km", deep_grad="deep_gradient_K_km",
                    z100="elev_100C_first_mRL", z200="elev_200C_first_mRL", min_bm="min_boiling_margin_degC",
                    water_level="water_level_est_mRL", regime="regime"),
    "Haiku2": dict(Tmax="max_T_degC", Tmax_elev="max_T_elev_mRL", T_TD="T_at_TD_degC", inversion="inversion_drop_K",
                   shallow_grad="shallow_grad_K_per_km", deep_grad="deep_grad_K_per_km",
                   z100="T100_crossings_mRL", z200="T200_crossings_mRL", min_bm="min_boiling_margin_K",
                   water_level="est_water_level_mRL", regime="regime"),
}

LOG = io.StringIO()


def log(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def extract(branch, dest):
    data = subprocess.run(["git", "archive", f"origin/{branch}", "Scripts", "Output", "CLAUDE.md"],
                          cwd=ROOT, capture_output=True, check=True).stdout
    with tarfile.open(fileobj=io.BytesIO(data)) as t:
        t.extractall(dest, filter="data")
    head = subprocess.run(["git", "log", "-1", "--format=%h %s", f"origin/{branch}"], cwd=ROOT,
                          capture_output=True, text=True, check=True).stdout.strip()
    return head


def first_number(v):
    """Haiku writes isotherm crossings as '-14; -992'. Take the first (shallowest) value."""
    if pd.isna(v):
        return float("nan")
    try:
        return float(str(v).split(";")[0])
    except ValueError:
        return float("nan")


def main():
    rerun = "--no-rerun" not in sys.argv
    os.makedirs(OUT, exist_ok=True)
    tmp = tempfile.mkdtemp(prefix="model_cmp_")
    log(f"Temporary work directory: {tmp}")
    run_dir = {}

    # ---- 1. extract and check deliverables -------------------------------------------------
    deliv_rows, file_rows = [], []
    # Compare git blob ids (the working tree may have CRLF line endings).
    def blob(ref):
        return subprocess.run(["git", "rev-parse", ref], cwd=ROOT, capture_output=True, text=True).stdout.strip()
    claude_main = blob("main:CLAUDE.md")
    for m in MODELS:
        d = os.path.join(tmp, m)
        os.makedirs(d)
        head = extract(m, d)
        run_dir[m] = d
        log(f"{m}: branch tip {head}")
        figs = sorted(os.listdir(os.path.join(d, "Output", "figures")))
        outs = [f for f in os.listdir(os.path.join(d, "Output")) if os.path.isfile(os.path.join(d, "Output", f))]
        tgt = pd.read_csv(os.path.join(d, "Output", "NS_calibration_targets.csv"))
        missing_cols = [c for c in REQUIRED_TARGET_COLS if c not in tgt.columns]
        with open(os.path.join(d, "Scripts", "natural_state_analysis.py"), encoding="utf-8") as f:
            script_lines = sum(1 for _ in f)
        with open(os.path.join(d, "Output", "NS_data_interpretation.md"), encoding="utf-8") as f:
            report_lines = sum(1 for _ in f)
        row = dict(model=m, branch_tip=head,
                   claude_md_same_as_main=blob(f"origin/{m}:CLAUDE.md") == claude_main,
                   required_files_present=sum(os.path.exists(os.path.join(d, p)) for p in REQUIRED_FILES),
                   required_files_total=len(REQUIRED_FILES),
                   output_files=len(outs), figures=len(figs),
                   targets_rows=len(tgt), targets_missing_required_cols=";".join(missing_cols) or "none",
                   targets_extra_cols=";".join(c for c in tgt.columns if c not in REQUIRED_TARGET_COLS) or "none",
                   script_lines=script_lines, report_lines=report_lines,
                   has_conceptual_sketch=any("conceptual" in f for f in figs),
                   n_sections=sum(f.startswith("section") for f in figs))
        deliv_rows.append(row)
        for f in figs:
            file_rows.append(dict(model=m, figure=f))
    deliv = pd.DataFrame(deliv_rows)

    # ---- 2. reproducibility: re-run each script against a copy of Baseline ------------------
    repro_rows = []
    if rerun:
        for m in MODELS:
            r = os.path.join(tmp, "rerun_" + m)
            os.makedirs(os.path.join(r, "Scripts"))
            os.makedirs(os.path.join(r, "Output", "figures"))
            shutil.copytree(os.path.join(ROOT, "Baseline"), os.path.join(r, "Baseline"))
            shutil.copy(os.path.join(run_dir[m], "Scripts", "natural_state_analysis.py"), os.path.join(r, "Scripts"))
            t0 = time.time()
            p = subprocess.run([sys.executable, "Scripts/natural_state_analysis.py"], cwd=r,
                               capture_output=True, text=True)
            dt = time.time() - t0
            committed = os.path.join(run_dir[m], "Output")
            regenerated = os.path.join(r, "Output")
            same, diff, not_regen = [], [], []
            for f in sorted(os.listdir(committed)):
                if not os.path.isfile(os.path.join(committed, f)) or f.endswith(".txt"):
                    continue
                g = os.path.join(regenerated, f)
                if not os.path.exists(g):
                    not_regen.append(f)
                elif sha(g) == sha(os.path.join(committed, f)):
                    same.append(f)
                else:
                    diff.append(f)
            n_fig_new = len(os.listdir(os.path.join(regenerated, "figures")))
            repro_rows.append(dict(model=m, exit_code=p.returncode, runtime_s=round(dt, 1),
                                   files_identical=len(same), files_differ=";".join(diff) or "none",
                                   files_not_regenerated=";".join(not_regen) or "none",
                                   figures_regenerated=n_fig_new))
            log(f"{m}: rerun exit {p.returncode} in {dt:.1f} s; identical {len(same)}; differ {diff}; "
                f"not regenerated {not_regen}")
            if diff:
                for f in diff:
                    a = open(os.path.join(committed, f), encoding="utf-8").read().splitlines()
                    b = open(os.path.join(regenerated, f), encoding="utf-8").read().splitlines()
                    nd = [i for i, (x, y) in enumerate(zip(a, b)) if x != y]
                    log(f"   {f}: {len(nd)} differing line(s), first: {a[nd[0]][:120] if nd else ''}")
    repro = pd.DataFrame(repro_rows)

    # ---- 3. per-well characteristics side by side --------------------------------------------
    summ = {}
    for m in MODELS:
        s = pd.read_csv(os.path.join(run_dir[m], "Output", "well_summary.csv")).set_index("well")
        c = SUMMARY_COLS[m]
        df = pd.DataFrame({k: s[v] for k, v in c.items()})
        for k in ("z100", "z200"):
            df[k] = df[k].map(first_number)
        verdict_col = "ns_verdict" if "ns_verdict" in s.columns else None
        summ[m] = df
    wells = summ["Opus2"].index
    rows = []
    for w in wells:
        for k in SUMMARY_COLS["Opus2"]:
            vals = {m: summ[m].loc[w, k] if w in summ[m].index else None for m in MODELS}
            r = dict(well=w, quantity=k, **vals)
            if k != "regime":
                nums = [float(v) for v in vals.values() if v is not None and not pd.isna(v)]
                r["spread"] = round(max(nums) - min(nums), 2) if len(nums) > 1 else None
            rows.append(r)
    well_cmp = pd.DataFrame(rows)
    spread = (well_cmp[well_cmp.quantity != "regime"].groupby("quantity")["spread"]
              .agg(["max", "median"]).reset_index())

    # ---- 4. natural-state verdicts and weights -----------------------------------------------
    verd = []
    for m in MODELS:
        pc = pd.read_csv(os.path.join(run_dir[m], "Output", "pressure_consistency_check.csv"))
        ws = pd.read_csv(os.path.join(run_dir[m], "Output", "well_summary.csv"))
        tg = pd.read_csv(os.path.join(run_dir[m], "Output", "NS_calibration_targets.csv"))
        if "verdict" in pc.columns and "confidence_level" in pc.columns:      # Sonnet
            v = pc.set_index("well")[["verdict", "confidence_level"]]
        else:
            v = ws.set_index("well")[["ns_verdict", "ns_confidence"]]
        v.columns = ["verdict", "confidence"]
        mw = tg.groupby("well")["weight"].agg(["mean", "sum"])
        for w in v.index:
            verd.append(dict(model=m, well=w, verdict=v.loc[w, "verdict"], confidence=v.loc[w, "confidence"],
                             mean_weight=round(mw.loc[w, "mean"], 3), sum_weight=round(mw.loc[w, "sum"], 2)))
    verd = pd.DataFrame(verd)
    verd_wide = verd.pivot(index="well", columns="model", values=["verdict", "confidence", "mean_weight"])
    verd_wide.columns = [f"{b}_{a}" for a, b in verd_wide.columns]
    verd_wide = verd_wide.reset_index()

    tgt = {m: pd.read_csv(os.path.join(run_dir[m], "Output", "NS_calibration_targets.csv")) for m in MODELS}
    wsum = pd.DataFrame([dict(model=m, points=len(t), weight_sum=round(t.weight.sum(), 1),
                              weight_mean=round(t.weight.mean(), 3),
                              weight_min=round(t.weight.min(), 3), weight_max=round(t.weight.max(), 3),
                              n_distinct_weights=t.weight.round(4).nunique(),
                              in_feedzone=int(t.in_feedzone.sum()), behind_casing=int(t.behind_casing.sum()),
                              pressure_weight_column="weight_P" in t.columns)
                         for m, t in tgt.items()])

    # ---- 5. point-level agreement: blocks, flags, Tsat ----------------------------------------
    key = ["well", "elevation_mRL"]
    base = tgt["Opus2"][key + ["block", "Tsat_degC", "in_feedzone", "behind_casing"]]
    agree_rows, block_diff = [], []
    for m in ["Sonnet2", "Haiku2"]:
        mm = base.merge(tgt[m][key + ["block", "Tsat_degC", "in_feedzone", "behind_casing"]], on=key,
                        suffixes=("_Opus2", "_" + m))
        agree_rows.append(dict(pair=f"Opus2 vs {m}", matched_points=len(mm),
                               block_same=int((mm.block_Opus2 == mm["block_" + m]).sum()),
                               in_feedzone_same=int((mm.in_feedzone_Opus2 == mm["in_feedzone_" + m]).sum()),
                               behind_casing_same=int((mm.behind_casing_Opus2 == mm["behind_casing_" + m]).sum()),
                               max_abs_Tsat_diff_degC=round((mm.Tsat_degC_Opus2 - mm["Tsat_degC_" + m]).abs().max(), 3)))
        dd = mm[mm.block_Opus2 != mm["block_" + m]]
        for _, r in dd.iterrows():
            block_diff.append(dict(model=m, well=r.well, elevation_mRL=r.elevation_mRL,
                                   block_Opus2=r.block_Opus2, block_other=r["block_" + m]))
    agree = pd.DataFrame(agree_rows)
    block_diff = pd.DataFrame(block_diff)

    # The model's own well-to-block file is the arbiter for WELL-2.
    import json
    wj = json.load(open(os.path.join(ROOT, "Baseline", "Teaching_Model_3D", "Models", "TE14515",
                                     "gTE14515_wells.json")))
    w2_cols = sorted({b["block"][:3] for b in wj["WELL-2"]["Blocks"]})
    log(f"gTE14515_wells.json WELL-2 columns: {w2_cols}")

    # ---- 6. layer centres ---------------------------------------------------------------------
    lc = {m: pd.read_csv(os.path.join(run_dir[m], "Output", "T_P_at_layer_centres.csv")) for m in MODELS}
    cen = pd.DataFrame({m: lc[m].groupby("layer").layer_centre_mRL.first() for m in MODELS}).dropna()
    # Layer centre from the LAYERS section of the geometry file (columns: index, bottom, centre).
    lay = {}
    with open(os.path.join(ROOT, "Baseline", "Teaching_Model_3D", "Models", "TE14515", "gTE14515.dat")) as f:
        on = False
        for line in f:
            if line.startswith("LAYERS"):
                on = True
                continue
            if on:
                p = line.split()
                if len(p) != 3:
                    break
                lay[int(p[0])] = (float(p[1]), float(p[2]))
    cen["geometry_file_centre"] = [lay[i][1] for i in cen.index]
    for m in MODELS:
        cen[f"{m}_minus_file"] = cen[m] - cen["geometry_file_centre"]
    cen = cen.reset_index()
    lcm = lc["Opus2"].merge(lc["Haiku2"], on=["well", "layer"], suffixes=("_Opus2", "_Haiku2"))
    layer_stats = pd.DataFrame([dict(
        model=m, rows=len(lc[m]),
        max_abs_centre_error_m=round(cen[f"{m}_minus_file"].abs().max(), 1),
        layers_with_centre_error=int((cen[f"{m}_minus_file"].abs() > 0.01).sum()),
        layers_compared=len(cen)) for m in MODELS])
    tdiff = (lcm.T_degC_Opus2 - lcm.T_degC_Haiku2).abs()
    log(f"Layer-centre T, Opus2 vs Haiku2 at the same layer: median |dT| {tdiff.median():.2f} C, "
        f"max {tdiff.max():.2f} C over {len(lcm)} pairs")
    layer_stats["T_vs_Opus2_median_abs_dT"] = [0.0, round((lc['Opus2'].merge(lc['Sonnet2'], on=['well', 'layer']).pipe(
        lambda x: (x.T_degC_x - x.T_degC_y).abs().median())), 2), round(tdiff.median(), 2)]
    layer_stats["T_vs_Opus2_max_abs_dT"] = [0.0, round((lc['Opus2'].merge(lc['Sonnet2'], on=['well', 'layer']).pipe(
        lambda x: (x.T_degC_x - x.T_degC_y).abs().max())), 2), round(tdiff.max(), 2)]

    # ---- 7. hydrostatic headline numbers -------------------------------------------------------
    hyd = []
    for m in MODELS:
        h = pd.read_csv(os.path.join(run_dir[m], "Output", "hydrostatic_check.csv"))
        if m == "Opus2":
            n = len(h)
            lens = h.length_m
        elif m == "Sonnet2":
            n = len(h)
            lens = h.length_m
        else:
            h = h[h.status == "single-phase"]
            n = len(h)
            lens = h.dz_m
        hyd.append(dict(model=m, windows=n, min_window_m=lens.min(), max_window_m=lens.max()))
    hyd = pd.DataFrame(hyd)

    # ---- write -------------------------------------------------------------------------------
    out = {
        "deliverables.csv": deliv, "reproducibility.csv": repro, "well_metrics_side_by_side.csv": well_cmp,
        "well_metric_spread_summary.csv": spread, "ns_verdicts_and_weights.csv": verd_wide,
        "weight_summary.csv": wsum, "point_agreement.csv": agree, "block_disagreements.csv": block_diff,
        "layer_centre_check.csv": cen, "layer_centre_summary.csv": layer_stats, "hydrostatic_windows.csv": hyd,
        "figures_by_model.csv": pd.DataFrame(file_rows),
    }
    for name, df in out.items():
        df.to_csv(os.path.join(OUT, name), index=False)
        log(f"wrote Output/model_comparison/{name} ({len(df)} rows)")
    log("\n== deliverables\n" + deliv.T.to_string())
    log("\n== reproducibility\n" + repro.to_string())
    log("\n== weights\n" + wsum.to_string())
    log("\n== verdicts\n" + verd_wide.to_string())
    log("\n== point agreement\n" + agree.to_string())
    log("\n== layer centres\n" + layer_stats.to_string())
    log("\n== metric spread\n" + spread.to_string())
    log("\n== hydrostatic windows\n" + hyd.to_string())
    shutil.rmtree(tmp, ignore_errors=True)
    with open(os.path.join(OUT, "compare_console_log.txt"), "w", encoding="utf-8") as f:
        f.write(LOG.getvalue())


if __name__ == "__main__":
    main()
