# Evidence-led Rust porting

Reviewed 2026-10-07. These are engineering lessons from public projects, not an
endorsement of their correctness or a claim that their throughput transfers here.

## Public examples reviewed

- [Bun's July 8, 2026 Rust rewrite](https://bun.com/blog/bun-in-rust) used a documented
  mapping from the old implementation, its existing language-independent test suite,
  and reviewers working separately from the implementation context. Apply that here
  by preserving LIMPA semantics before refactoring and requiring independent review
  before merging. We did not copy Bun code or import its agent instructions.
- [Anthropic's February 2026 C-compiler experiment](https://www.anthropic.com/engineering/building-c-compiler)
  emphasizes a strong verifier, regression checks, deterministic sampled tests and a
  known-good oracle. Its own evaluation reports significant limits despite impressive
  demos. Our R reference remains the oracle; a passing subset is not full equivalence.
- [FrankenSQLite's public repository](https://github.com/Dicklesworthstone/frankensqlite)
  describes comparison against original SQLite through a conformance harness. Its
  [testing guidance](https://github.com/Dicklesworthstone/frankensqlite/blob/main/AGENTS.md)
  includes error cases, generated tests and retaining rejected optimization evidence.
  Those are useful practices; repository claims alone do not establish production
  readiness. We retain failed numerical comparisons in BENCHMARKS.md.

## Contract for this port

1. Preserve filtering, missingness semantics, grouping, sample order, priors,
   optimization initialization/termination, and SE definitions. Compilation alone
   is not evidence of equivalence. List remaining R calls and fallback paths.
2. Freeze the reference version, data, seed, initialization groups and tolerances
   before comparing. Do not relax tolerance or drop a failing case to obtain a pass.
   Unsupported cases must fail or use an explicitly documented reference fallback.
3. Compare identical inputs under both implementations. Keep scalar likelihood /
   derivative / Hessian checks, an R-generated oracle and complete pipeline checks.
   Keep library tests separate from tests that exercise the real command-line process.
4. Check worker-count determinism with nonidentical proteins across a batch boundary.
   Exercise truncated input, invalid numeric parameters and existing-output protection.
   Require debug and release Rust tests on Linux and macOS; do not claim Windows
   support from these checks. Forbid unsafe code in our own two Rust crates. This
   does not assert that dependencies contain no unsafe code.
5. Record CPU, memory, versions, worker counts, input dimensions and timing boundary.
   Retain fixed seeds and raw aggregate metrics. Run repeated, uncontended trials
   before making release-grade speed claims; the current benchmark remains a single
   run. A subset ratio is not a full-workload speedup, and expanded data is not a
   real 384-sample experiment.
6. Have an independent reviewer challenge scientific semantics and the benchmark
   design before merging. This document does not claim that such a review has run.
   Fuzzing, a broader randomized R differential corpus and the complete real
   384-sample comparison remain future work; no passing status is assigned to them.

## Changes applied in this follow-up

Added CLI worker-count / batch-order and failure-path tests, Linux/macOS debug/release
CI coverage, explicit unsafe-code prohibition, SPDX headers, a documented GPLv3-or-later
choice and the missing renv MIT notice. Numerical algorithms and acceptance thresholds
are unchanged; the earlier timings remain measurements of that numerical implementation.
