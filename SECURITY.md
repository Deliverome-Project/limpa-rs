# Reporting a security issue

Use GitHub's private vulnerability reporting for this repository when enabled.
Do not attach instrument reports, patient/sample metadata, access tokens or other
private inputs to a public issue. For ordinary bugs, provide a minimal synthetic
reproduction and the output of `python -m limpa_rs doctor` with local usernames
and private paths removed.

Only the current release candidate is maintained at this stage. Scientific
production acceptance is documented in `docs/release-readiness.md`.

## CI supply-chain defenses

CI installs Python build/test tools and runtime dependencies from `uv.lock`, with
hashes and a seven-day resolver cooldown. A standard-library check verifies every
locked Python version's age before installation, including the first lockfile;
missing metadata and network failures fail the check. Third-party Python source
builds are disabled in CI, and package builds use the locked build environment.
The final wheel is installed without resolving more dependencies.

Workflows use immutable action commits, read-only tokens, nonpersistent checkout
credentials, and timeouts. Security checks scan full Git history with Gitleaks,
audit the Python lock with pip-audit and Cargo.lock with RustSec, and lint workflow
security with zizmor. Downloaded scanner archives have fixed SHA-256 checksums.
Checks also run weekly after this workflow lands on main. Action and Cargo updates
are proposed by Dependabot and remain subject to review. To refresh Python
versions, run `uv lock --upgrade`, review the diff, and submit a pull request.

These measures reduce exposure to supply-chain campaigns such as Shai-Hulud; they
do not establish that every dependency is safe. The Python cooldown does not cover
Rust crates, R packages, or the dependencies inside third-party actions. Rust and
R remain pinned in their lockfiles; R dependencies have no automated advisory
scanner here. Changes to the R environment still require equivalence assessment.
The committed CI lock does not constrain dependencies resolved by downstream pip
users; analyses should keep their own reviewed dependency lock.

Current Rust audit warning: `paste` 1.0.15 is an unmaintained transitive dependency
(RUSTSEC-2024-0436). It is retained to avoid changing the assessed numerical stack
in a CI hardening change. Unmaintained warnings are visible but do not fail CI;
known vulnerabilities, unsoundness, and yanked releases do. Reassess the warning
when updating the numerical dependencies. GitHub secret scanning, push protection,
and Dependabot vulnerability alerts/security updates are enabled for this repo.
