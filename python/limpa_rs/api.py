# SPDX-License-Identifier: GPL-3.0-or-later
"""Pandas interface adapted from deliverome-analysis PR #141; see NOTICE."""

from __future__ import annotations

import json
import os
import shutil
import warnings
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from .runtime import ASSETS, _run

R_SCRIPT = ASSETS / "scripts" / "limpa_spectronaut.R"
COMPARATORS = ("maxlfq_proda",)


@dataclass
class LimpaResult:
    """Outputs of one limpa run. All matrices are indexed by ``protein_id``."""

    protein_log2: pd.DataFrame
    protein_se: pd.DataFrame
    proteins: pd.DataFrame
    samples: pd.DataFrame
    summary: dict
    de: dict[str, pd.DataFrame] = field(default_factory=dict)
    dpc_points: pd.DataFrame | None = (
        None  # precursor, mu, n_detected, nsamples (DPC fit set)
    )
    maxlfq_log2: pd.DataFrame | None = (
        None  # MaxLFQ protein x run with NAs (with ``compare``)
    )
    # comparator DE tables: method -> {coefficient/contrast -> topTable-shaped DataFrame}
    compare: dict[str, dict[str, pd.DataFrame]] = field(default_factory=dict)
    outdir: Path | None = None

    @property
    def dpc(self) -> tuple[float, float]:
        """(beta0, beta1) of the fitted detection probability curve."""
        return self.summary["dpc_beta0"], self.summary["dpc_beta1"]


def detection_probability(log2_intensity, beta0: float, beta1: float):
    """P(precursor detected | true log2 intensity) under the fitted DPC (a logistic curve)."""
    x = np.asarray(log2_intensity, dtype=float)
    return np.exp(-np.logaddexp(0, -(beta0 + beta1 * x)))


def run_limpa_spectronaut(
    report: str | os.PathLike | None = None,
    *,
    matrix_dir: str | os.PathLike | None = None,
    cores: int = 1,
    engine: str = "rust",
    seed: int = 141,
    reference_cores: int = 1,
    outdir: str | os.PathLike | None = None,
    protein_id: str = "PG.ProteinGroups",
    q_cutoff: float = 0.01,
    samples: pd.DataFrame | None = None,
    formula: str | None = None,
    contrasts: list[str] | None = None,
    compare: list[str] | None = None,
    verbose: bool = False,
) -> LimpaResult:
    """Quantify proteins with limpa from a Spectronaut Normal Report, optionally test DE.

    Args:
        report: local path to the Spectronaut Normal Report TSV (download S3 keys first with
            ``processed.download(key)``). Give this or ``matrix_dir``.
        matrix_dir: output of :func:`deliverome_analysis.spectronaut_trim.trim_spectronaut_report`
            -- the compact precursor matrix; much faster and lighter for large cohorts.
        cores: Rust threads (or R worker processes for engine="reference").
        engine: "rust" (default) or pinned "reference" for comparison.
        seed: shared DPC subsampling seed, default 141.
        reference_cores: initialization chunks to match PR141's parallel R fitter.
            Keep 1 for serial-reference starts regardless of Rust thread count.
        outdir: where the R driver writes its TSVs. Default: a temp dir (outputs are still
            returned as DataFrames); pass a path to keep them.
        protein_id: report column used to group precursors into proteins.
        q_cutoff: EG.Qvalue / PG.Qvalue threshold; values above it are treated as missing.
        samples: optional run-level table with a ``run`` column matching ``R.FileName`` plus
            any covariates for ``formula``. Spectronaut's own ``R.Condition`` /
            ``R.Replicate`` are always available without it.
        formula: R model formula over the sample columns, e.g. ``"~ group"`` or
            ``"~ 0 + group"``. Omit for quantification only.
        contrasts: limma contrasts over the design's column names, e.g. ``["grouptrt-groupctrl"]``
            (use with a ``~ 0 + ...`` formula). Without it, one table per non-intercept
            coefficient.
        compare: comparator pipelines run on the same filtered precursors. Currently
            ``["maxlfq_proda"]`` (MaxLFQ via iq, then proDA), the runner-up kept for ongoing
            evaluation; the matrix is returned as ``maxlfq_log2`` and tables under
            ``compare["maxlfq_proda"]``. Needs ``formula`` for DE; without it only MaxLFQ is
            returned. (The full six-method comparison is S5E7: analysis/sprint5/s5e7-dia-limpa-benchmark/.)
        verbose: echo R's progress messages.

    Returns:
        LimpaResult with protein_log2, protein_se, proteins, samples, summary and de tables.
    """
    # Validate arguments before touching the filesystem or R, so bad calls fail the same way
    # on machines without R (e.g. CI).
    if (report is None) == (matrix_dir is None):
        raise ValueError("give exactly one of report or matrix_dir")
    if contrasts and not formula:
        raise ValueError("contrasts need a formula")
    unknown = set(compare or ()) - set(COMPARATORS)
    if unknown:
        raise ValueError(
            f"unknown compare method(s) {sorted(unknown)}; choose from {COMPARATORS}"
        )
    if engine not in ("rust", "reference"):
        raise ValueError("engine must be rust or reference")
    for name, value in (("cores", cores), ("reference_cores", reference_cores)):
        if isinstance(value, bool) or not isinstance(value, int) or value < 1:
            raise ValueError(f"{name} must be a positive integer")
    if (
        isinstance(seed, bool)
        or not isinstance(seed, int)
        or not 0 <= seed <= 2147483647
    ):
        raise ValueError("seed must be an integer between 0 and 2147483647")
    if not 0 <= q_cutoff <= 1:
        raise ValueError("q_cutoff must be between 0 and 1")
    if samples is not None and (
        "run" not in samples
        or samples["run"].isna().any()
        or samples["run"].duplicated().any()
    ):
        raise ValueError("samples needs a unique, nonmissing 'run' column")
    # absolute paths: R runs from the repo root, not the caller's cwd
    src = Path(report if report is not None else matrix_dir).resolve()
    if not src.exists():
        raise FileNotFoundError(src)
    keep = outdir is not None
    destination = Path(outdir).resolve() if keep else None
    if destination is not None:
        if destination.exists() and (
            not destination.is_dir() or any(destination.iterdir())
        ):
            raise FileExistsError(
                "outdir must be empty to prevent mixing old and new results"
            )
        destination.parent.mkdir(parents=True, exist_ok=True)
    # Write alongside the destination, then atomically rename only a fully validated run.
    out = Path(
        tempfile.mkdtemp(prefix=".limpa-rs-", dir=destination.parent if keep else None)
    )

    cmd = [
        str(R_SCRIPT),
        f"--report={src}" if report is not None else f"--matrix-dir={src}",
        f"--cores={cores}",
        f"--engine={engine}",
        f"--seed={seed}",
        f"--reference-cores={reference_cores}",
        f"--outdir={out}",
        f"--protein-id={protein_id}",
        f"--q-cutoff={q_cutoff}",
    ]
    try:
        if samples is not None:
            if "run" not in samples.columns:
                raise ValueError(
                    "samples needs a 'run' column matching Spectronaut R.FileName"
                )
            samples_path = out / "samples_input.tsv"
            samples.to_csv(samples_path, sep="\t", index=False)
            cmd.append(f"--samples={samples_path}")
        if formula:
            cmd.append(f"--formula={formula}")
        if contrasts:
            cmd.append(f"--contrasts={','.join(contrasts)}")
        if compare:
            cmd.append(f"--compare={','.join(compare)}")
        proc = _run(cmd, verbose=verbose)
        for line in proc.stderr.splitlines():
            if line.startswith("WARNING:"):
                warnings.warn(line, RuntimeWarning, stacklevel=2)
        result = _read_outputs(out, keep)
        _validate_result(result)
        if destination is not None:
            os.replace(out, destination)  # cannot replace a nonempty directory
            result.outdir = destination
        return result
    finally:
        shutil.rmtree(out, ignore_errors=True)


def _read_outputs(out: Path, keep: bool) -> LimpaResult:
    def tsv(name: str, index: str | None = "protein_id") -> pd.DataFrame:
        columns = pd.read_csv(out / name, sep="\t", nrows=0).columns
        df = pd.read_csv(
            out / name,
            sep="\t",
            keep_default_na=False,
            na_values={col: ["NA"] for col in columns if col != index},
            dtype={index: str} if index else None,
        )
        return df.set_index(index) if index else df

    summary = json.loads((out / "summary.json").read_text())
    de = {name: tsv(fname) for name, fname in (summary.get("de_tables") or {}).items()}
    compare: dict[str, dict[str, pd.DataFrame]] = {}
    for key, fname in (summary.get("compare_tables") or {}).items():
        method, coef = key.split("::", 1)
        compare.setdefault(method, {})[coef] = tsv(fname)
    maxlfq = out / "maxlfq_log2.tsv"
    return LimpaResult(
        protein_log2=tsv("protein_log2.tsv"),
        protein_se=tsv("protein_se.tsv"),
        proteins=tsv("protein_annotation.tsv"),
        samples=tsv("samples.tsv", index="run"),
        summary=summary,
        de=de,
        dpc_points=tsv("dpc_points.tsv", index=None)
        if (out / "dpc_points.tsv").exists()
        else None,
        maxlfq_log2=tsv("maxlfq_log2.tsv") if maxlfq.exists() else None,
        compare=compare,
        outdir=out if keep else None,
    )


def _validate_result(result: LimpaResult) -> None:
    """Reject incomplete or misaligned scientific outputs before publishing them."""
    e, se = result.protein_log2, result.protein_se
    if (
        e.empty
        or not e.index.is_unique
        or not e.columns.is_unique
        or not e.index.equals(se.index)
        or not e.columns.equals(se.columns)
        or not e.index.equals(result.proteins.index)
        or list(e.columns) != list(result.samples.index)
    ):
        raise RuntimeError(
            "LIMPA output dimensions, identifiers or ordering are inconsistent"
        )
    try:
        values, errors = e.to_numpy(dtype=float), se.to_numpy(dtype=float)
    except (TypeError, ValueError) as exc:
        raise RuntimeError("LIMPA returned nonnumeric protein outputs") from exc
    if (
        not np.isfinite(values).all()
        or not np.isfinite(errors).all()
        or (errors <= 0).any()
    ):
        raise RuntimeError(
            "LIMPA returned nonfinite estimates or invalid standard errors"
        )
    if result.summary.get("n_proteins") != len(e) or result.summary.get(
        "n_runs"
    ) != len(e.columns):
        raise RuntimeError("LIMPA summary counts disagree with output matrices")
