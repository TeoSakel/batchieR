# AGENTS.md

Act as an expert R package developer and statistical software engineer. Work on `batchieR` with an emphasis on correctness, maintainability, reproducibility, stable interfaces, and efficient numerical code. Follow the standard R package toolchain built around `devtools`, `roxygen2`, `testthat`, `renv`, and `cli`.

## General guidelines

- Understand the affected code, tests, and public behavior before editing.
- Make small, idiomatic changes and preserve unrelated work.
- Add focused regression tests for behavior changes and bug fixes.
- Keep tests deterministic and numerical assertions tolerance-aware.
- Preserve model parameterization, constraints, priors, RNG behavior, and the dimensions, names, and ordering of posterior objects unless intentionally changing them.
- Test numerical routines against analytic cases, limiting cases, or small seeded simulations; report non-finite values and convergence problems clearly.
- Use roxygen2 comments for documentation and `cli` for user-facing conditions.
- Manage the environment with `renv`; avoid unnecessary dependencies and incidental `renv.lock` changes.
- Use the devtools pipeline rather than sourcing files or generating package artifacts with custom commands.
- Report what was validated and any remaining risk accurately.

## Tiered workflow

1. **Iterate:** Use `devtools::load_all()` and focused `devtools::test()` filters while developing.
2. **Integrate:** Use `devtools::test()` for the full testthat suite and `devtools::document()` when roxygen2 source changes. Use parameter recovery or simulation-based calibration for major statistical changes.
3. **Check:** Use `devtools::check()` for broad, metadata, dependency, or release changes, or when explicitly requested.

Escalate between tiers when the scope, risk, or test results justify it. At handoff, state the highest tier completed and its outcome. Do not run full package development checks on every iteration.
