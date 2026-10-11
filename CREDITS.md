# Scientific credits and citation guidance

Cite the original LIMPA methods when publishing papers, preprints, or reports
using limpa-rs results. Also identify limpa-rs and its version or commit. The
statistical methods are the work of the original authors; limpa-rs accelerates
the protein-fitting core while retaining R stages.

- Mengbo Li, S. A. Cobbold, Gordon K. Smyth (2025). Quantification and differential
  analysis of mass spectrometry proteomics data with probabilistic recovery of
  information from missing values. https://doi.org/10.1101/2025.04.28.651125
- Mengbo Li, Gordon K. Smyth (2023). Neither random nor censored: estimating
  intensity-dependent probabilities for missing values in label-free proteomics.
  Bioinformatics 39(5). https://doi.org/10.1093/bioinformatics/btad200

Original LIMPA software authors: Mengbo Li and Gordon Smyth; contributor:
Pedro Baldoni. Reference version: LIMPA 1.4.2, pinned source commit
`d9ad5218207434cfcbe080328db5e79d39a1623a`. Reference environment: limma 3.68.5,
R 4.6.1. See NOTICE for R Core's vmmin optimizer and its J. C. Nash and
B. D. Ripley credits, source provenance, and dependency notices. No upstream
endorsement is implied.

## Suggested methods wording

Protein quantification used limpa-rs [version/commit], a Rust acceleration of
the LIMPA 1.4.2 protein-fitting core (Li, Cobbold & Smyth, 2025; Li & Smyth, 2023).
Input processing, detection-probability and empirical Bayes estimation, and
optional limma differential analysis remained in R.

Replace the version placeholder and describe only the stages actually used.
For additional methods used in a run, cite the corresponding upstream packages;
the run's CITATIONS.txt appends their R citation guidance and runtime versions.

## AI assistance disclosure

AI coding assistants have been used to help write this package's code and
documentation. OpenAI Codex assisted with the attribution documentation and
export changes in this update. This is a disclosure about this package's
development, not a claim about how the upstream software was written. A complete
historical inventory of AI tools, model versions, and per-file contributions is
not maintained. The human maintainer remains responsible for release decisions
and scientific claims; AI assistance does not establish correctness or imply
that every contribution has received independent human review.
