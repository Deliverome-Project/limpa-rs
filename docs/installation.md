# Install and use from an analysis

`limpa-rs` is an experimental **Python package with a bundled Rust executable**.
Python >=3.11 and R 4.6.x are required. macOS and Linux are supported. Installing
from source also needs a Rust toolchain compatible with the locked dependencies
(the tested compiler is Rust 1.99); installing a matching wheel does not need Rust.
Release wheels target Apple Silicon macOS 11+ and Linux x86-64 (tested on the
Ubuntu GitHub runner); the Linux wheel is not advertised as manylinux-portable.
Build from source for other compatible machines, including Intel Macs.
R remains necessary for DPC estimation, global parameters and differential analysis.
There is no CRAN/Bioconductor `library(limpaRs)` package or public PyPI release.

## Install

Install a matching wheel from the
[GitHub release](https://github.com/Deliverome-Project/limpa-rs/releases/tag/v0.1.0-rc.1),
or build the tagged source:

```sh
uv pip install "git+https://github.com/Deliverome-Project/limpa-rs.git@v0.1.0-rc.1"
python -m limpa_rs setup-r
python -m limpa_rs doctor
python -m limpa_rs validate
```

Run these inside your analysis Python environment. `python -m pip install "git+https://github.com/Deliverome-Project/limpa-rs.git@v0.1.0-rc.1"`
works too. Installation compiles the Rust executable in release mode with
`Cargo.lock` enforced and includes the R driver, lockfile and validation fixtures.
The explicit `setup-r` step downloads/restores R dependencies; importing the package
never installs software. It does not install R itself.

For a reproducible dependency in a `uv` analysis project, use a reviewed commit:

```sh
uv add "limpa-rs @ git+https://github.com/Deliverome-Project/limpa-rs.git@<reviewed-commit-sha>"
uv run python -m limpa_rs setup-r
uv run python -m limpa_rs doctor
```

Replace the placeholder with the full reviewed commit SHA. The repository is public; no GitHub token is needed to read it. Commit your
analysis project's `pyproject.toml` and `uv.lock` through its normal review process.
No changes to `deliverome-analysis` are required just to use this package.

## R isolation and configuration

By default, R packages are restored beneath
`~/.cache/limpa-rs/<lockfile-hash>` (or `$XDG_CACHE_HOME/limpa-rs/...`). A new lock
gets a new directory. Analysis runs use that explicit project, independent of
where Python was launched. No R libraries are written inside Python site-packages.
The **entire** R lock must match; both the wrapper and the R bridge check the
assessed LIMPA version. R 4.6.x is enforced at runtime.

- `LIMPA_RS_RSCRIPT`: explicit Rscript path; `DELIVEROME_RSCRIPT` is also accepted
  for compatibility with deliverome-analysis.
- `LIMPA_RS_R_PROJECT`: optional location of the pinned R project. Set it before
  setup and analysis. An existing PR141 project can be reused if its lock matches
  exactly. A different lock is rejected, never overwritten.
- `LIMPA_RS_BIN`: developer-only override of the packaged executable; run validation
  after overriding it. Normal users do not need this setting.

`python -m limpa_rs doctor` reports the selected paths and validated LIMPA version.
If setup fails, fix the reported R/network/system-library issue and rerun setup.

## Same analysis interface as PR141

Change the import:

```python
from limpa_rs import run_limpa_spectronaut

result = run_limpa_spectronaut(
    "/path/to/NormalReport.tsv",
    cores=4,
    seed=141,
    formula="~ R.Condition",  # omit for quantification only
    outdir="/path/to/new-results",
)
protein_log2 = result.protein_log2
protein_se = result.protein_se
annotations = result.proteins
de_tables = result.de
```

The `LimpaResult` fields match `deliverome_analysis.limpa`: pandas protein matrices,
annotations, sample metadata, summary, DPC points and optional DE/comparator tables.
`detection_probability` and `limpa_available` are also available. This is interface
compatibility, not Python class identity with Deliverome's original `LimpaResult`.
S3 download remains the analysis project's responsibility; pass local paths here.

Normal Report columns and q-value/imputation filtering match PR141. The input must
be precursor-level, with `EG.IsImputed` included when Spectronaut has imputed values.
The output directory must be empty; use a new directory for each run. Without
`outdir`, outputs are returned as DataFrames and temporary files are cleaned up,
including on failure. Persistent results are staged beside the requested destination
and published only after estimates, uncertainty and IDs have been checked. Local
macOS/Linux filesystems supporting atomic directory renames and hard links are required.

For large 384-run exports, keep using the existing streaming trimmer:

```python
from deliverome_analysis.spectronaut_trim import trim_spectronaut_report
from limpa_rs import run_limpa_spectronaut

# Prepare the compact matrix using PR141's existing trimmer and its documented API.
result = run_limpa_spectronaut(matrix_dir="/path/to/trimmed", cores=4)
```

Give exactly one of `report` or `matrix_dir`. Existing `samples`, `formula`,
`contrasts`, `protein_id`, `q_cutoff` and `compare=["maxlfq_proda"]` options remain.
Additional options are `engine="rust"` (default) or `"reference"`, `seed=141`, and
`reference_cores=1`. For a matched four-worker R comparison, use `cores=4` in both
runs and `reference_cores=4` for Rust; starting-value chunks affect early stopping.

## Equivalence and release status

`python -m limpa_rs validate` checks the **installed** executable with one and four
threads against frozen LIMPA 1.4.2 outputs. See `validation/README.md` for candidate
LIMPA upgrade assessment in a separate environment. Passing these synthetic checks
does not establish full real 384-sample equivalence or >=20x full-workload speedup.
The package remains experimental until that assessment is complete.

Build source archives and wheels with `uv build` (or `python -m build`). Test wheels
from outside the checkout; developer source imports do not exercise bundled assets.
CI builds from the source archive, installs the resulting wheel, and runs package
and R integration tests. There is no automatic publication to PyPI. Public releases are explicit, tagged GitHub releases.
When distributing binaries, provide the matching corresponding-source bundle (including
vendored Rust dependencies) and retained notices; see `LICENSING.md`.
