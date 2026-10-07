# Reference pin and upgrade checks

The production reference is **LIMPA 1.4.2**, source commit
`d9ad5218207434cfcbe080328db5e79d39a1623a` in `bioc/limpa`.
`renv.lock` pins this revision and the rest of the R environment (R 4.6.1,
limma 3.68.5). `scripts/install-r-deps.sh` restores and verifies it.
Both R bridge functions and both Spectronaut engines reject other versions.
When installed package metadata provides a source SHA, it must also match;
standard repository installations without SHA metadata are checked by version.
Use the lockfile restore to obtain the exact assessed source revision.

## Regression checks for future Rust changes

From the restored project root, after `cargo build --release --locked`:

```sh
Rscript validation/test_modules.R
Rscript validation/reference.R check
```

`check` compares the pinned R fitter and Rust with 1 and 4 threads to the
**committed, frozen** `fixtures/limpa-1.4.2.R` oracle. It includes synthetic
3/6/32/384-sample inputs, varying precursor counts, missing protein/sample
combinations, the all-missing-row fallback, singleton dispatch, and the bridge's
hyperparameter and annotation handling. The fixture contains the actual inputs,
reference outputs and package provenance, using hexadecimal doubles to preserve
exact values across text serialization. It contains no experimental data.
Do not regenerate it to make a failing test pass.

The reusable modules are:

- `modules/fixtures.R`: deterministic synthetic inputs, reference/Rust evaluation,
  and extraction of estimates, uncertainty and metadata.
- `modules/compare.R`: saved-output comparison, fixed absolute tolerance of 0.001
  for every log2 estimate and standard error, exact IDs/order/dimensions, missing
  value patterns, observation counts, annotations, targets and hyperparameters.
- `test_modules.R`: deliberately corrupts outputs to prove failures are detected,
  including numerical drift, nonfinite values, reordering and metadata changes.

The original `equivalence.R`, `end_to_end.R`, `edge_cases.R`, `real_data.R`,
`scale.R` and `parallel_reference.R` remain additional validation. The frozen
suite does not replace end-to-end filtering/DE tests, private Astral data checks,
or the missing full real 384-sample assessment.

## Assess a newer LIMPA without changing the production pin

Use a **separate R project/library** containing the candidate version; do not
update this repository's restored library or lockfile to run the assessment.
The script does not load this project's renv automatically when called by
absolute path from the candidate project's working directory. Ensure the intended
candidate library is active there. Replace `X.Y.Z` below with its exact version:

```sh
Rscript /path/to/limpa-rs/validation/reference.R capture X.Y.Z /tmp/candidate.rds
Rscript /path/to/limpa-rs/validation/reference.R compare \
  /path/to/limpa-rs/validation/fixtures/limpa-1.4.2.R /tmp/candidate.rds
```

`capture` checks the installed candidate version, evaluates it on the SAME stored
inputs, saves provenance and outputs, and immediately compares them with the old
oracle. Failure exits nonzero while retaining the candidate for inspection. It
refuses to overwrite an existing snapshot. `compare` needs only base R and can
compare any two saved snapshots with identical inputs. Only load trusted RDS or
R fixtures; R's `dget` evaluates R expressions.

A passing synthetic comparison is evidence, not automatic approval to upgrade.
Before promoting a version, run the end-to-end and real-data validations against
the candidate, review numerical and scientific differences, then update the lock,
bridge contract and a separately named oracle together in a PR. Retain the old
oracle as the record of the earlier assessment; do not loosen tolerances merely
to accommodate drift. Record reference/candidate versions and observed errors in
the PR. No production version bypass is provided.
