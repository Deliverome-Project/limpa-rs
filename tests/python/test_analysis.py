# SPDX-License-Identifier: GPL-3.0-or-later
"""limpa wrapper — runs the real R driver on a small simulated Spectronaut Normal Report.

The simulation plants three things the wrapper must get right:

* censoring by a logistic detection curve (low-intensity precursors go missing),
* 40 truly differential proteins (logFC 2) among 300, so DE has a known answer,
* poison rows that must never reach the quantification: values flagged ``EG.IsImputed`` and
  values above the q-value cutoff, both set to an absurd intensity.

Tests that need R are skipped when Rscript / limpa are unavailable (e.g. CI without R).
"""

from __future__ import annotations

import os

import numpy as np
import pandas as pd
import pytest

import limpa_rs as lp
from limpa_rs.api import COMPARATORS

needs_limpa = pytest.mark.skipif(
    os.environ.get("LIMPA_RS_INTEGRATION") != "1",
    reason="set LIMPA_RS_INTEGRATION=1 with pinned R restored",
)

N_PROT, N_PREC, N_DE = 300, 3, 40
RUNS = [f"run_{g}{i}" for g in "AB" for i in (1, 2, 3)]
POISON = 2.0**40  # log2 = 40; any leak shows up as a protein far above the ~20 range


def _simulate(rng: np.random.Generator) -> pd.DataFrame:
    rows = []
    for p in range(N_PROT):
        prot = f"P{p:05d}"
        mu = rng.normal(18, 2)
        lfc = 2.0 if p < N_DE else 0.0
        for k in range(N_PREC):
            offset = rng.normal(0, 1)
            for run in RUNS:
                x = (
                    mu
                    + offset
                    + (lfc if run.startswith("run_B") else 0)
                    + rng.normal(0, 0.3)
                )
                detected = rng.random() < lp.detection_probability(x, -12.0, 0.75)
                rows.append(
                    {
                        "R.FileName": run,
                        "R.Condition": run[4],
                        "R.Replicate": int(run[-1]),
                        "PG.ProteinGroups": prot,
                        "PG.ProteinAccessions": prot,
                        "EG.ModifiedSequence": f"_PEP{p}K{k}_",
                        "FG.Charge": 2,
                        "EG.TotalQuantity (Settings)": 2.0**x if detected else np.nan,
                        "EG.Qvalue": 0.001 if detected else np.nan,
                        "PG.Qvalue": 0.001,
                        "EG.IsImputed": False,
                    }
                )
    df = pd.DataFrame(rows)
    # Poison: an imputed value and a failed-q value, both absurdly high, on two proteins.
    imp = (df["PG.ProteinGroups"] == "P00299") & (df["R.FileName"] == "run_A1")
    df.loc[imp, ["EG.TotalQuantity (Settings)", "EG.Qvalue", "EG.IsImputed"]] = [
        POISON,
        0.001,
        True,
    ]
    badq = (df["PG.ProteinGroups"] == "P00298") & (df["R.FileName"] == "run_B1")
    df.loc[badq, ["EG.TotalQuantity (Settings)", "EG.Qvalue"]] = [POISON, 0.2]
    return df.dropna(
        subset=["EG.TotalQuantity (Settings)"]
    )  # Spectronaut omits undetected rows


@pytest.fixture(scope="module")
def report(tmp_path_factory) -> str:
    path = tmp_path_factory.mktemp("limpa") / "Report.tsv"
    _simulate(np.random.default_rng(7)).to_csv(path, sep="\t", index=False)
    return str(path)


def test_detection_probability_is_logistic():
    b0, b1 = -10.0, 0.5
    assert lp.detection_probability(20.0, b0, b1) == pytest.approx(
        0.5
    )  # midpoint at -b0/b1
    p = lp.detection_probability([10, 20, 30], b0, b1)
    assert np.all(np.diff(p) > 0)


def test_missing_report_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        lp.run_limpa_spectronaut(tmp_path / "nope.tsv")


@needs_limpa
def test_quantification_complete_and_poison_excluded(report):
    res = lp.run_limpa_spectronaut(report)
    n_reported = pd.read_csv(report, sep="\t")["PG.ProteinGroups"].nunique()
    assert res.protein_log2.shape == (n_reported, len(RUNS))
    assert not res.protein_log2.isna().any().any()  # DPC-Quant gives complete data
    assert (res.protein_se > 0).all().all()
    # limpa orders runs by first appearance in the file; matrices and samples must agree.
    assert sorted(res.samples.index) == sorted(RUNS)
    assert list(res.protein_log2.columns) == list(res.samples.index)
    assert list(res.protein_se.columns) == list(res.samples.index)
    assert res.summary["imputed_values_filtered"] is True
    # Neither the imputed nor the failed-q poison value may leak into a protein estimate.
    assert res.protein_log2.loc["P00299", "run_A1"] < 30
    assert res.protein_log2.loc["P00298", "run_B1"] < 30
    assert 0 < res.summary["precursor_missing_fraction"] < 0.6
    assert res.summary["dpc_beta1"] > 0  # detection rises with intensity


@needs_limpa
def test_de_recovers_planted_proteins(report):
    res = lp.run_limpa_spectronaut(report, formula="~ R.Condition")
    tt = res.de["R.ConditionB"]
    hits = set(tt.index[tt["adj.P.Val"] < 0.05])
    truth = {f"P{p:05d}" for p in range(N_DE)}
    assert len(hits & truth) >= 0.8 * N_DE  # power
    assert len(hits - truth) <= 0.1 * max(len(hits), 1)  # FDR roughly controlled
    assert tt.loc[sorted(truth & hits)[0], "logFC"] == pytest.approx(2, abs=0.6)


@needs_limpa
def test_samples_table_and_contrasts(report):
    samples = pd.DataFrame({"run": RUNS, "group": ["ctrl"] * 3 + ["trt"] * 3})
    res = lp.run_limpa_spectronaut(
        report, samples=samples, formula="~ 0 + group", contrasts=["grouptrt-groupctrl"]
    )
    assert list(res.de) == ["grouptrt-groupctrl"]
    # Covariates must land on the right run regardless of run order in the report.
    expected = samples.set_index("run")["group"]
    assert res.samples["group"].to_dict() == expected.to_dict()
    assert (
        res.samples["group"]
        == res.samples["R.Condition"].map({"A": "ctrl", "B": "trt"})
    ).all()


@needs_limpa
def test_no_replicates_is_a_clear_error(tmp_path):
    df = _simulate(np.random.default_rng(1))
    df = df[df["R.FileName"].isin(["run_A1", "run_B1"])]
    path = tmp_path / "two_runs.tsv"
    df.to_csv(path, sep="\t", index=False)
    lp.run_limpa_spectronaut(path)  # quantification alone is fine with 2 runs
    with pytest.raises(RuntimeError, match="no replicates"):
        lp.run_limpa_spectronaut(path, formula="~ R.Condition")


@needs_limpa
def test_missing_required_column_is_a_clear_error(tmp_path):
    df = _simulate(np.random.default_rng(2)).drop(columns=["PG.Qvalue"])
    path = tmp_path / "noq.tsv"
    df.to_csv(path, sep="\t", index=False)
    with pytest.raises(RuntimeError, match=r"PG\.Qvalue"):
        lp.run_limpa_spectronaut(path)


def test_unknown_comparator_rejected(tmp_path):
    path = tmp_path / "r.tsv"
    path.write_text("x\n")
    with pytest.raises(ValueError, match="unknown compare"):
        lp.run_limpa_spectronaut(
            path, formula="~ R.Condition", compare=["maxlfq", "nope"]
        )


@needs_limpa
def test_comparators_run_on_same_input(report):
    res = lp.run_limpa_spectronaut(
        report, formula="~ R.Condition", compare=list(COMPARATORS)
    )
    assert set(res.compare) == set(lp.COMPARATORS)
    # MaxLFQ keeps NAs (proteins unobserved in a run); limpa's matrix is complete.
    assert res.maxlfq_log2.shape == res.protein_log2.shape
    assert set(res.maxlfq_log2.index) == set(res.protein_log2.index)
    truth = {f"P{p:05d}" for p in range(N_DE)}
    for method, tabs in res.compare.items():
        tt = tabs["R.ConditionB"]
        assert {"logFC", "P.Value", "adj.P.Val"} <= set(tt.columns), method
        hits = set(tt.index[tt["adj.P.Val"] < 0.05])
        assert len(hits & truth) >= 0.8 * N_DE, method
        # every pipeline must see the poison-free data: no absurd fold changes
        assert tt["logFC"].abs().max() < 15, method


@needs_limpa
def test_both_engines_agree_through_installed_api(report):
    reference = lp.run_limpa_spectronaut(
        report, engine="reference", seed=141, formula="~ R.Condition"
    )
    rust = lp.run_limpa_spectronaut(report, cores=4, seed=141, formula="~ R.Condition")
    pd.testing.assert_frame_equal(
        reference.protein_log2, rust.protein_log2, atol=1e-3, rtol=0
    )
    pd.testing.assert_frame_equal(
        reference.protein_se, rust.protein_se, atol=1e-3, rtol=0
    )
    pd.testing.assert_frame_equal(reference.proteins, rust.proteins)
    pd.testing.assert_frame_equal(reference.samples, rust.samples)
    a, b = (
        reference.de["R.ConditionB"].sort_index(),
        rust.de["R.ConditionB"].sort_index(),
    )
    np.testing.assert_allclose(
        a[["logFC", "P.Value", "adj.P.Val"]],
        b[["logFC", "P.Value", "adj.P.Val"]],
        atol=1e-3,
        rtol=0,
    )
    assert set(a.index[a["adj.P.Val"] < 0.05]) == set(b.index[b["adj.P.Val"] < 0.05])


@needs_limpa
def test_compact_matrix_input(report, tmp_path):
    # Generate the same file contract as PR141's streaming trimmer, using its
    # reference reader so this tests package transport rather than a second filter.
    import json
    from limpa_rs.runtime import _run

    script = tmp_path / "matrix.R"
    script.write_text(
        "\n".join(
            [
                "suppressPackageStartupMessages(library(limpa))",
                "args <- commandArgs(trailingOnly=TRUE)",
                'y <- readSpectronaut(args[1], annotation.columns=c("PG.ProteinGroups","PG.ProteinAccessions"), q.cutoffs=.01, filter.columns="EG.IsImputed", filter.values=TRUE)',
                "wide <- data.frame(precursor=rownames(y$E),y$E,check.names=FALSE)",
                "ann <- data.frame(precursor=rownames(y$E),y$genes,check.names=FALSE)",
                'names(ann)[names(ann)=="PG.ProteinGroups"] <- "protein_id"',
                'nanoparquet::write_parquet(wide,file.path(args[2],"precursor_log2.parquet"))',
                'nanoparquet::write_parquet(ann,file.path(args[2],"precursors.parquet"))',
                'data.table::fwrite(data.frame(run=colnames(y$E),y$targets),file.path(args[2],"samples.tsv"),sep="\\t")',
            ]
        )
    )
    _run([str(script), report, str(tmp_path)])
    (tmp_path / "manifest.json").write_text(
        json.dumps({"imputed_flag_present": True, "source_report": report})
    )
    a = lp.run_limpa_spectronaut(report)
    b = lp.run_limpa_spectronaut(matrix_dir=tmp_path, cores=4)
    pd.testing.assert_frame_equal(a.protein_log2, b.protein_log2, atol=1e-3, rtol=0)
    pd.testing.assert_frame_equal(a.protein_se, b.protein_se, atol=1e-3, rtol=0)
