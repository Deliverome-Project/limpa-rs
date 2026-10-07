# Validation record

Measured October 7, 2026 on an Apple M2, 16 GiB RAM, eight logical CPUs,
macOS/aarch64. Rust runs use four worker threads and a release build with thin LTO;
R 4.6.1, LIMPA 1.4.2, limma 3.68.5, rustc 1.99.0. Dependencies are pinned in
`renv.lock` and `Cargo.lock`. These are single-run measurements, not confidence
intervals or portable performance guarantees.

## Real Astral export: full output comparison

Input is PR #141's cached Spectronaut Astral Zoom HeLa export:

```
processed://scratch/ellaine/spectronaut_benchmarks/run_4/20261006_211223_astral_zoom_1/astral_zoom_1_Report_BGS Factory Report (Normal).tsv
ETag: 3fe65feea59289c24af657d073c7b950-112
```

All 303,817 precursors / 11,010 proteins / three samples were used. The inputs were
read and filtered once; both quantifiers used the same detection-curve fit.

| Measurement | Result |
|---|---:|
| Read/filter | 3.838 s |
| DPC fit | 0.552 s |
| Original serial R quantification | 22.273 s |
| Rust quantification including R bridge | 2.071 s |
| Quantification speedup | 10.75× |
| Maximum absolute protein log2 difference | 0.000000067484 |
| Maximum absolute SE difference | 0.000000000087859 |

Every protein/sample identifier, ordering, annotation, observation count and sample
table matched. Both numerical differences pass the **predeclared absolute 0.001**
tolerance. Three-sample speedup is below the 20× target, which is for 384 samples.

## Expanded 384-sample workload

This is an **expanded workload, not a real 384-sample experiment**. The three real
precursor columns are repeated to 384, then independent normal technical noise
(SD 0.1 log2) is added using seed 384141. The input preserves the original
303,817 precursors, 11,010 protein groups, precursor counts per protein and empirical
missingness patterns. It does not model new missingness patterns or biological
heterogeneity that could occur in 384 distinct real samples.

| Full expanded workload measurement | Result |
|---|---:|
| Global DPC and hyperparameter estimation | 10.155 s |
| Quantification including starting values, IPC and output readback | 54.644 s |
| Rust executable within that quantification | 45.308 s |
| Total compact-matrix-to-protein time | **64.799 s** |
| Time target | ≤900 s: pass on this workload |

The timing excludes generating the expanded matrix, reading/trimming the raw
Spectronaut long report, R startup/build/install, DE and comparator pipelines.
**Do not describe this as an end-to-end 100+ GB report benchmark.**

A deterministic 48-protein subset spanning the sorted IDs was also fit in serial R
and Rust with the same input, hyperparameters and starting-value procedure:

| Matched 48-protein × 384-sample comparison | Result |
|---|---:|
| Original serial R | 114.582 s |
| Rust including bridge | 0.243 s |
| Measured subset speedup | **471.53×** |
| Maximum log2 difference | 0.000492562 |
| Maximum SE difference | 0.000000020167 |

The subset passes the numerical and 20× speed thresholds. **The complete 11,010 ×
384 R baseline has not been run.** A full-workload speedup and all-protein numerical
agreement at 384 samples have therefore not been demonstrated. No extrapolated
reference time is reported as a measurement.

## Parallel reference initialization

PR #141's four-worker R path computes initial imputation separately in each
protein chunk. Comparing that path with global serial initialization gave a
maximum 0.001999834 log2 difference: **failed**, despite a large speedup. That
failed comparison is preserved in `benchmarks/parallel-initialization-mismatch.tsv`.

Rust now accepts the same initialization groups (`--reference-cores=4`) while
keeping computational thread count separate. This preserves either the serial or
parallel reference's starts without changing the scientific tolerance. The
four-worker matched-initialization comparison uses 48 proteins × 384 samples:

| Matched parallel comparison | Result |
|---|---:|
| PR #141 reference, four R workers | 47.755 s |
| Rust, four threads including bridge | 0.532 s |
| Measured subset speedup | **89.77×** |
| Maximum log2 difference | 0.000427150 |
| Maximum SE difference | 0.000000045994 |
| Fixed 0.001 equivalence threshold | Pass |

Matching the starts alone initially left a 0.001056595 log2 discrepancy. Preserving
R's column-major arithmetic reduction order and expression grouping in the
likelihood/gradient resolved that last failure; no tolerance was relaxed.
The 384-sample full-run timing above uses serial-compatible starts, while this
comparison deliberately uses the parallel reference's chunk-specific starts.

## Additional checks

- Synthetic oracle at 3, 6, 32 and 384 samples, with 1–50 precursors per protein,
  missing cells and entire missing protein/sample combinations. Fully unobserved
  precursor rows are removed for this standard input test; separate tests cover
  exact reference fallback for them.
- End-to-end Spectronaut simulation: 299 retained proteins, 870 precursors, six
  samples, 40 planted differential proteins; imputed and failed-q poison values
  are excluded. All protein estimates and SEs pass the same 0.001 tolerance.
- DE logFC, p-value and adjusted p-value errors pass 0.001, with the identical set
  of significant calls at adjusted p<0.05. Annotations, samples and DPC points
  match exactly. This does not imply identical threshold calls for every dataset.
- Five Rust tests cover an R-generated numerical oracle, analytic complete-data
  solution, finite-difference derivatives, structured versus dense Hessian
  inversion, and invalid/unsupported input rejection.
- R fallback tests verify exact results for all-missing precursor rows and the
  all-singleton estimator.

## Why the compatibility optimizer is the default

An initial structured Newton solver found the same model's optimum more accurately
but did not always reproduce the reference's early stopping behavior. On synthetic
inputs with completely unobserved precursor rows, discrepancies reached roughly
0.011 log2 and 0.071 SE. That implementation was **not accepted as equivalent**.
The default now reproduces R's BFGS rules; completely unobserved precursor inputs
use the original R fitter explicitly. The tolerance was not increased to accept
failed results. The Newton solver remains an experimental library entry point.

Raw validation tables are under `benchmarks/`. To reproduce, run the scripts listed
in README using the pinned R environment. The private input and expanded matrices
are never committed.

Follow-up practice review: three additional CLI regression tests now cover identical
outputs across worker counts/batches, malformed inputs, and overwrite protection.
These supplement the five numerical unit tests; no numerical code or tolerance
changed. See PORTING.md for the review sources and validation still outstanding.
