# limpa-rs

Rust acceleration of LIMPA 1.4.2 protein quantification for Deliverome Spectronaut
DIA data. The numerical core preserves LIMPA's likelihood, 16-node normal
quadrature, sum-to-zero peptide effects, protein priors, R-compatible BFGS stopping
rules, and observed-Hessian standard errors. Independent proteins run in parallel.

**This is a hybrid implementation, not a complete R-free LIMPA rewrite.** R still
reads/filters inputs, estimates the detection probability curve and empirical
Bayes hyperparameters, supplies reference starting values, and performs limma DE.
Keeping those stages shared isolates the equivalence test to the expensive core.
The Spectronaut entry point preserves the output files used by
[deliverome-analysis PR #141](https://github.com/Deliverome-Project/deliverome-analysis/pull/141).

Target: at least 20× faster quantification on approximately 11,000 proteins and
384 samples, with matrix-to-protein processing within 900 seconds on a local
machine. See [BENCHMARKS.md](BENCHMARKS.md) for measured results and the exact timing
boundaries. Do not interpret a subset speedup as a measured full-workload speedup.

## Build and run

Install Rust and R 4.6.x, then:

```sh
cargo build --release --locked
scripts/install-r-deps.sh
Rscript scripts/limpa_spectronaut.R \
  --report=/path/to/NormalReport.tsv --outdir=/path/to/output \
  --engine=rust --cores=4 --seed=141
```

For a large 384-run report, first use PR #141's streaming
`deliverome_analysis.spectronaut_trim.trim_spectronaut_report`, then supply
`--matrix-dir=/path/to/trimmed` in place of `--report`. This prevents loading a
potentially 100+ GB full Spectronaut export as an R long table. Trimming is separate
from quantification and must be timed separately when assessing report-to-results
latency. The compact matrix itself is roughly 0.9 GB for 303,817 × 384 doubles.

Optional flags from PR #141 remain available: `--samples`, `--formula`,
`--contrasts`, `--protein-id`, `--q-cutoff`, and `--compare=maxlfq_proda`.
Use `--engine=reference --cores=1 --seed=141` for the original serial fitter.
The seed makes DPC random subsampling reproducible across separate invocations.
Rust's `--cores` controls worker threads; its starting values match serial R by
default, independently of the thread count. To reproduce PR #141's R worker
chunks, add `--reference-cores=4` when comparing with `--engine=reference --cores=4`.
This is necessary because PR #141 computes its initial imputation separately in
each worker chunk. Changing initialization chunks can change early-stopped fits.
All current validations use normalized Spectronaut quantities without a second
normalization step.

The R bridge can also be sourced directly:

```r
source("R/limpa_rs.R")
Sys.setenv(LIMPA_RS_BIN = normalizePath("target/release/limpa-rs"))
set.seed(141)
dpcfit <- limpa::dpc(y)
protein <- dpc_quant_rust(y, "PG.ProteinGroups", dpcfit, cores = 4L)
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

The R validations require the pinned `renv.lock` environment. They can also be run
from the restored PR #141 checkout while `LIMPA_RS_ROOT` points here. Private data
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
- The core assumes unweighted observations, matching this Spectronaut pipeline.
  It is not a replacement for all of LIMPA's public APIs.
- Reaching BFGS's 100-iteration limit fails explicitly; the original R implementation
  can return `convergence=1` without the caller checking it.
- High precursor counts still make BFGS and the peptide factorization expensive.
- Differential-expression computation and comparator pipelines remain in R;
  their runtime is additional to the quantification benchmarks.
- A failed CLI run can leave an incomplete output file. The bridge checks process
  success before reading it and cleans its temporary files. The CLI refuses to
  overwrite an existing output file.

## Provenance

Reference: LIMPA 1.4.2, limma 3.68.5, R 4.6.1 and PR #141 commit
`32c563c95410637e0431faf27eab1e7a9786464b`. The driver and R package lock are copied
from that commit with a selectable Rust backend and deterministic seed added.
The optimizer is adapted from R Core's `vmmin`; see [NOTICE](NOTICE) and [LICENSE](LICENSE).

## License and porting practice

Licensed under **GPL-3.0-or-later**; see [LICENSE](LICENSE),
[upstream notices and rationale](LICENSING.md), and [NOTICE](NOTICE).
The vendored renv bootstrap retains its separate MIT license.
[PORTING.md](PORTING.md) records the public Claude-assisted rewrites reviewed,
practices adopted here, and validation still outstanding.
