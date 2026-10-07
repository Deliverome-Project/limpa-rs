# SPDX-License-Identifier: GPL-3.0-or-later
import json
import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import limpa_rs as lp
from limpa_rs import api, runtime


def test_wheel_contains_runtime_and_binary():
    assert runtime._binary().is_file()
    for name in (
        "R/limpa_rs.R",
        "scripts/limpa_spectronaut.R",
        "renv.lock",
        "renv/activate.R",
        "validation/reference.R",
        "validation/modules/compare.R",
        "validation/fixtures/limpa-1.4.2.R",
    ):
        assert (runtime.ASSETS / name).is_file(), name
    lock = json.loads((runtime.ASSETS / "renv.lock").read_text())
    assert lock["Packages"]["limpa"]["Version"] == "1.4.2"
    assert (
        lock["Packages"]["limpa"]["RemoteSha"]
        == "d9ad5218207434cfcbe080328db5e79d39a1623a"
    )
    assert subprocess.run([str(runtime._binary())], capture_output=True).returncode != 0


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(cores=0),
        dict(cores=1.5),
        dict(cores=True),
        dict(engine="newest"),
        dict(seed=-1),
        dict(q_cutoff=float("nan")),
        dict(reference_cores=0),
        dict(samples=pd.DataFrame({"run": ["a", "a"]})),
    ],
)
def test_validation_precedes_r(kwargs):
    with pytest.raises(ValueError):
        lp.run_limpa_spectronaut("absent.tsv", **kwargs)


def test_failure_cleans_temporary_outputs(tmp_path, monkeypatch):
    report = tmp_path / "input.tsv"
    report.touch()
    scratch = tmp_path / "output"
    scratch.mkdir()
    monkeypatch.setattr(api.tempfile, "mkdtemp", lambda **kwargs: str(scratch))

    def fail(*args, **kwargs):
        raise RuntimeError("simulated R failure")

    monkeypatch.setattr(api, "_run", fail)
    with pytest.raises(RuntimeError, match="simulated"):
        lp.run_limpa_spectronaut(report)
    assert not scratch.exists()


def test_existing_outputs_are_not_mixed(tmp_path):
    report = tmp_path / "input.tsv"
    report.touch()
    out = tmp_path / "out"
    out.mkdir()
    (out / "old.tsv").touch()
    with pytest.raises(FileExistsError):
        lp.run_limpa_spectronaut(report, outdir=out)


def test_mismatched_r_project_rejected(tmp_path):
    (tmp_path / "renv").mkdir()
    (tmp_path / "renv/activate.R").touch()
    (tmp_path / "renv.lock").write_text("{}")
    with pytest.raises(RuntimeError, match="lock differs"):
        runtime._check_project(tmp_path)


def test_identifier_strings_preserved(tmp_path):
    (tmp_path / "summary.json").write_text("{}")
    for name in ("protein_log2.tsv", "protein_se.tsv", "protein_annotation.tsv"):
        (tmp_path / name).write_text("protein_id\tx\n001\t1\nNA\tNA\n")
    (tmp_path / "samples.tsv").write_text("run\tgroup\n001\tA\nNA\tB\n")
    result = api._read_outputs(tmp_path, True)
    assert list(result.protein_log2.index) == ["001", "NA"]
    assert list(result.samples.index) == ["001", "NA"]
    assert np.isnan(result.protein_log2.loc["NA", "x"])


@pytest.mark.skipif(
    os.environ.get("LIMPA_RS_INTEGRATION") != "1", reason="requires restored R"
)
def test_installed_reference_equivalence():
    assert lp.limpa_available()
    lp.validate()


def _valid_result():
    e = pd.DataFrame(
        [[18.0, 19.0]], index=pd.Index(["P1"], name="protein_id"), columns=["a", "b"]
    )
    return lp.LimpaResult(
        e,
        e * 0 + 0.2,
        pd.DataFrame(index=e.index),
        pd.DataFrame(index=pd.Index(["a", "b"], name="run")),
        {"n_proteins": 1, "n_runs": 2},
    )


@pytest.mark.parametrize("mutation", ["nonfinite", "negative_se", "ordering", "counts"])
def test_invalid_outputs_rejected(mutation):
    result = _valid_result()
    if mutation == "nonfinite":
        result.protein_log2.iloc[0, 0] = np.nan
    elif mutation == "negative_se":
        result.protein_se.iloc[0, 0] = -0.1
    elif mutation == "ordering":
        result.samples = result.samples.iloc[::-1]
    else:
        result.summary["n_runs"] = 3
    with pytest.raises(RuntimeError):
        api._validate_result(result)


def test_failed_persistent_run_never_publishes_partial_files(tmp_path, monkeypatch):
    report = tmp_path / "input.tsv"
    report.touch()
    destination = tmp_path / "final"

    def fail(args, **kwargs):
        out = Path(
            next(arg.split("=", 1)[1] for arg in args if arg.startswith("--outdir="))
        )
        (out / "protein_log2.tsv").write_text("incomplete")
        raise RuntimeError("interrupted run")

    monkeypatch.setattr(api, "_run", fail)
    with pytest.raises(RuntimeError, match="interrupted"):
        lp.run_limpa_spectronaut(report, outdir=destination)
    assert not destination.exists()
    assert not list(tmp_path.glob(".limpa-rs-*"))
