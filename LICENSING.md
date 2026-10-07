# License and upstream provenance

limpa-rs author: **Rebecca Carlson**, with the Deliverome Project.
Copyright (c) 2026 Rebecca Carlson and limpa-rs contributors. This project
notice preserves the separate original LIMPA, R Core, and dependency credits.

The project license is **GPL-3.0-or-later**, consistently declared in `Cargo.toml`,
`LICENSE`, and SPDX headers on the Rust core and R bridge/driver. The existing GPLv3
license text remains unchanged. The public release retains all upstream notices.

This is a source-informed port: its BFGS implementation adapts R Core's `vmmin`,
whose notice permits GPL version 2 or later. LIMPA 1.4.2 declares `GPL (>=2)` in
its installed package metadata. We preserve those origins in `NOTICE`; generating
or translating code with an AI tool is not a reason to remove attribution or
label this work a clean-room implementation.

GPLv3-or-later is an appropriate choice for this combination: the R-derived code
allows a later GPL version, and our direct numerical dependency nalgebra is
Apache-2.0. GNU documents both the later-version route and Apache-2.0 compatibility
with GPLv3 in its [GPLv3 guide](https://www.gnu.org/licenses/quick-guide-gplv3.en.html).
We therefore do not offer the combined implementation under MIT or Apache alone.

| Component | Recorded license | Treatment |
|---|---|---|
| Project code and R-Core-derived optimizer | GPL-3.0-or-later; optimizer upstream GPL-2.0-or-later | Root LICENSE, SPDX headers and R Core copyright in NOTICE |
| LIMPA 1.4.2 reference package | GPL (>=2) | Installed through renv; provenance retained |
| Vendored renv 1.3.1 bootstrap | MIT | Preserve Posit's notice and complete text in LICENSES/renv-MIT.txt |
| nalgebra 0.34.2 | Apache-2.0 | Dependency's license remains in force |
| rayon 1.12.0 | MIT OR Apache-2.0 | Dependency's license remains in force |

`LICENSES/cargo-metadata.json` records license declarations for the locked Cargo
resolution, including optional packages. It is a metadata inventory, not a binary
SBOM or a substitute for dependency license texts. Dependencies include permissive
MIT, Apache-2.0, Zlib and Unicode-3.0 declarations. Refresh this inventory when the
lockfile changes. Binary releases include dependency license texts/notices and a corresponding-source
bundle with vendored Rust dependencies. These releases do not alter upstream grants.

The renv bootstrap remains unmodified. Its separate MIT text is reconstructed from
the installed 1.3.1 package's LICENSE fields and R's canonical MIT template, with
that provenance included alongside the text. The R optimizer's original notice is
available in [R's optim.c](https://github.com/wch/r-source/blob/trunk/src/appl/optim.c).

## Python distributions

Python wheels retain the root license, NOTICE, renv's MIT license and the
verbatim Rust dependency notices in `LICENSES/rust-dependencies.txt`. The latter
covers the locked resolution, including optional dependencies, and should be
refreshed when Cargo.lock changes. Wheels bundle the Rust executable and original
R/runtime sources; the matching source distribution includes Cargo.lock and all
Rust sources. The corresponding-source release bundle additionally vendors the
locked Rust dependency sources and Cargo offline configuration. Distribute that
matching bundle alongside binary wheels. Python
and R dependencies installed separately retain their own notices and licenses.
