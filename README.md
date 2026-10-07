# limpa-rs

Rust acceleration of LIMPA 1.4.2 protein quantification from precursor-level
proteomics data. The numerical core preserves LIMPA's likelihood, 16-node normal
quadrature, sum-to-zero peptide effects, protein priors, R-compatible BFGS stopping
rules, and observed-Hessian standard errors. Independent proteins run in parallel.

**This is a hybrid implementation, not a complete R-free LIMPA rewrite.** R still
reads/filters inputs, estimates the detection probability curve and empirical
Bayes hyperparameters, supplies reference starting values, and performs limma DE.
Keeping those stages shared isolates the equivalence test to the expensive core.

## Input formats

The numerical core is independent of the instrument and upstream search software.
It works with precursor-level log2 intensities, missing observations, and
precursor-to-protein assignments—not an already aggregated protein matrix.

- **Spectronaut Normal Reports:** the Python convenience function reads the report
  and applies the supported quality and imputation filters.
- **Prepared precursor matrices:** the same Python function accepts `matrix_dir`.
  See the [compact matrix format](docs/installation.md#compact-matrix-input) for
  the required files. Filtering and normalization must be done upstream.
- **Other data readers:** the source-level R bridge accepts a LIMPA-compatible
  `EList` with a precursor matrix and protein annotations (example below).

The Python function retains the name `run_limpa_spectronaut` for compatibility;
it is not a universal report reader. Other search-engine exports need an adapter
or preparation into the matrix format. Biological suitability and numerical
agreement must be assessed for each new input workflow.

Target: at least 20× faster quantification on approximately 11,000 proteins and
384 samples, with matrix-to-protein processing within 900 seconds on a local
machine. See [BENCHMARKS.md](BENCHMARKS.md) for measured results and the exact timing
boundaries. Do not interpret a subset speedup as a measured full-workload speedup.

**Release status:** v0.1.0-rc.1 is a public release candidate. Production use on
an actual ~384-sample cohort has not yet been validated. See
[release acceptance criteria](docs/release-readiness.md) before replacing a validated
analysis workflow.

## Install as a Python package

Python >=3.11, R 4.6.x, macOS/Linux. From this checkout, in your analysis environment:

```sh
uv pip install .
python -m limpa_rs setup-r
python -m limpa_rs doctor
python -m limpa_rs validate
```

Source installation uses the pinned Rust toolchain; a matching compiled wheel bundles the executable.
The setup step explicitly restores the pinned R environment into a user cache.

```python
from limpa_rs import run_limpa_spectronaut

result = run_limpa_spectronaut("NormalReport.tsv", cores=4, seed=141)
result.protein_log2
result.protein_se
```

The Python API returns pandas matrices, annotations, and optional differential-analysis
tables. See [installation and analysis integration](docs/installation.md) for
tagged Git installs, reproducible version pins, DE, compact matrices, setup paths
and troubleshooting. This is an **experimental installable package**, not a public
PyPI release or a completed real 384-sample validation.

## Build and run from source

Install Rust and R 4.6.x, then:

```sh
cargo build --release --locked
scripts/install-r-deps.sh
Rscript scripts/limpa_spectronaut.R \
  --report=/path/to/NormalReport.tsv --outdir=/path/to/output \
  --engine=rust --cores=4 --seed=141
```

For large reports, prepare a compact precursor matrix and supply
`--matrix-dir=/path/to/matrix` in place of `--report`. See the
[compact matrix format](docs/installation.md#compact-matrix-input). This avoids
loading the full long-format report into R. Input preparation is separate from
quantification and must be timed separately for report-to-results benchmarks.

Optional flags include `--samples`, `--formula`, `--contrasts`, `--protein-id`,
`--q-cutoff`, and `--compare=maxlfq_proda`.
Use `--engine=reference --cores=1 --seed=141` for the original serial fitter.
The seed makes DPC random subsampling reproducible across separate invocations.
Rust's `--cores` controls worker threads; its starting values match serial R by
default, independently of the thread count. For comparison with a four-worker R
reference, add `--reference-cores=4` to the Rust run and use
`--engine=reference --cores=4` for the reference run. The parallel reference
computes initial imputation separately in each worker chunk, so changing chunks
can change early-stopped fits. Current real-data validations use normalized
Spectronaut quantities without a second normalization step.

The R bridge can also be sourced directly from a checkout using the pinned R
environment. Here `precursor_log2` is a numeric precursor-by-sample matrix with
unique row/column names and `NA` for missing observations; `protein_ids` supplies
one nonmissing protein assignment per row. Inputs must already be filtered and
appropriately normalized; this call does not read or filter vendor reports:

```r
source("R/limpa_rs.R")
Sys.setenv(LIMPA_RS_BIN = normalizePath("target/release/limpa-rs"))
library(limma)
y <- new("EList", list(
  E = precursor_log2,
  genes = data.frame(protein_id = protein_ids, row.names = rownames(precursor_log2))
))
set.seed(141)
dpcfit <- limpa::dpc(y)
protein <- dpc_quant_rust(y, "protein_id", dpcfit, cores = 4L)
```

## Equivalence and tests

```sh
cargo test --locked
cargo clippy --all-targets --locked -- -D warnings
Rscript validation/test_modules.R
Rscript validation/reference.R check
export LIMPA_RS_ROOT="$PWD"
Rscript validation/equivalence.R
Rscript validation/end_to_end.R
Rscript validation/edge_cases.R
export LIMPA_REPORT=/path/to/astral_NormalReport.tsv
Rscript validation/real_data.R
Rscript validation/scale.R
Rscript validation/parallel_reference.R
```

The R validations require the pinned `renv.lock` environment and
`LIMPA_RS_ROOT` pointing to this checkout. Private data
and generated matrices stay outside git; only aggregate benchmark metrics are
committed. Rust CI runs independently, including an R-generated synthetic numerical
oracle. A separate R workflow restores the pinned environment and runs the frozen
reference suite and synthetic end-to-end checks. Neither workflow claims that
private-data or full-workload validation ran.

The production bridge enforces LIMPA **1.4.2**; the lock also pins its exact source
commit. See [the upgrade validation guide](validation/README.md) for reusable
modules that compare future Rust builds and candidate LIMPA versions against
frozen outputs without changing the production pin.

Acceptance is absolute error ≤0.001 for every protein log2 estimate and standard
error, with matching IDs, dimensions, ordering, observation counts and annotations.
The end-to-end synthetic check also compares log fold changes, p-values, adjusted
p-values and the exact set of calls at adjusted p<0.05, and poisons imputed and
failed-q inputs to verify filtering. These finite-data checks are evidence, not a
proof of equivalence for every possible dataset. Outputs are not bit-identical.

## Numerical design and limitations

The reference builds a dense design matrix for every protein. Rust evaluates the
same sums directly. For standard errors it exploits the sample Hessian's diagonal
plus rank-one structure and factors a peptide Schur complement, avoiding a dense
sample-by-sample inverse. Work is buffered in batches of 64 proteins; only these
fits run concurrently, and results retain input protein order.

Matching the statistical optimum alone was insufficient: a Newton solver changed
some weakly identified estimates beyond the fixed tolerance. The default therefore
ports R's BFGS iteration and termination rules. `fit_newton` remains an experimental
library function, is never selected by the CLI or bridge, and has no blanket
output-equivalence claim.

- Completely unobserved precursor rows use the original R fitter with a warning;
  their peptide effects are weakly identified and optimizer dependent. Raw Rust
  callers receive an error for those rows.
- LIMPA uses a separate estimator when *every* protein has one precursor. The R
  bridge delegates that case unchanged to LIMPA.
- The core currently supports unweighted observations.
  It is not a replacement for all of LIMPA's public APIs.
- Reaching BFGS's 100-iteration limit fails explicitly; the original R implementation
  can return `convergence=1` without the caller checking it.
- High precursor counts still make BFGS and the peptide factorization expensive.
- Differential-expression computation and comparator pipelines remain in R;
  their runtime is additional to the quantification benchmarks.
- The native CLI publishes its output only after success and refuses to overwrite
  an existing file. The Python API validates complete, finite, aligned outputs before
  atomically publishing the result directory; a failed run publishes no partial result.

## Provenance

Original LIMPA methods: Li, Cobbold & Smyth (2025), and Li & Smyth (2023).
See [CITATION.cff](CITATION.cff) for the full references.

Reference: LIMPA 1.4.2, limma 3.68.5, and R 4.6.1.
The optimizer is adapted from R Core's `vmmin`. Full source provenance and
contributor acknowledgments are retained in [NOTICE](NOTICE); see also [LICENSE](LICENSE).

## License and porting practice

Licensed under **GPL-3.0-or-later**; see [LICENSE](LICENSE),
[upstream notices and rationale](LICENSING.md), and [NOTICE](NOTICE).
The vendored renv bootstrap retains its separate MIT license.
[PORTING.md](PORTING.md) records the public Claude-assisted rewrites reviewed,
practices adopted here, and validation still outstanding.
