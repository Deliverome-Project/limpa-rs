# Release acceptance

## v0.1.0-rc.1

This is a release candidate. Passing automated software checks is not the same as
validating a production scientific workflow on its intended cohort.

Completed checks include a pinned LIMPA/R dependency environment; frozen reference
outputs at 3, 6, 32 and 384 samples; matched estimates/uncertainty and metadata;
end-to-end filtering and differential-expression calls; installed-wheel tests;
error/failure cleanup; atomic result publication; Rust debug/release checks; and
full repository-history secret scanning. The original authors, licenses and
modified optimizer provenance are retained. Release assets include source and
third-party notices. See BENCHMARKS.md for measured data and limitations.

## Required before a stable production release for the intended workload

1. Obtain the actual ~384-sample Astral cohort and its sample/design metadata.
2. Run the pinned reference and Rust on the same complete cohort, with identical
   filtering, random seed and initialization chunks. Compare every estimate and SE
   at absolute tolerance 0.001, and require exact IDs, dimensions, order and metadata.
3. If DE is used, compare effect sizes, p-values, adjusted p-values and the exact
   significant calls at adjusted p<0.05 on the intended designs/contrasts.
4. Measure both complete runs on the same local machine and thread budget. Confirm
   >=20x speedup and <=900 seconds within a stated timing boundary. Record input
   provenance, hardware, versions and per-stage timings; do not extrapolate a
   subset speedup or label expanded data as an actual biological cohort.
5. Have the analysis owner review these results before replacing the existing
   analysis workflow. Retain its outputs for rollback and pin the adopted release.

The real three-run cohort is validated. The expanded 384-column workload is a
stress test, not evidence completing items 1–5. Until those items are satisfied,
do not advertise a production-ready full-cohort replacement or promote this tag
to a stable release. The reference engine remains available for side-by-side checks.

Releases are built from a tested commit. The release process checks CI success,
builds source and platform wheels, retains corresponding source and licenses,
records SHA-256 checksums, and only then publishes the tag and release assets.
