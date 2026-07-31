# Development guide

`batchieR` is an R package for Bayesian latent-factor models of combination-response screens. Keep changes small, documented, reproducible, and compatible with standard Bayesian R workflows.

## Package workflow

- Use `renv` for dependency reproducibility and the devtools ecosystem for development: `usethis` for package setup, `roxygen2` for documentation, `testthat` for tests, and `devtools` for local checks.
- Never edit `NAMESPACE` or generated `.Rd` files by hand. Document exported objects with roxygen2, then run `devtools::document()`.
- Put user-facing functions in `R/`, Stan programs in `inst/stan/`, tests in `tests/testthat/`, and longer examples in vignettes.
- Declare runtime dependencies deliberately in `Imports`; keep development-only or optional integrations in `Suggests`. Avoid attaching packages in package code and use qualified calls (`pkg::fun`).
- Before finishing, run `devtools::test()`, `devtools::document()`, and `devtools::check()`. Add focused tests for every bug fix or new behavior.
- Use `cli::cli_abort`/`cli::cli_warn`/`cli::cli_inform` for user-facing messages.

## Bayesian and Stan conventions

- Keep the R interface independent of the Stan backend where practical. Validate data and priors in R, use stable parameter names, and make seeds, chains, warmup, sampling iterations, and parallelism explicit and reproducible.
- Return or provide conversion to `posterior` draw formats (prefer `draws_array`/`draws_df`) so results work naturally with `bayesplot`, `posterior`, and downstream tidy tooling.
- Expose pointwise log-likelihood draws with observation indexing suitable for `loo::loo()` and `loo::waic()`. Preserve observation identifiers through fitting, prediction, and diagnostics.
- Provide posterior predictive draws and generated quantities needed for model checking. Follow `bayesplot` naming conventions where useful and avoid custom plotting infrastructure when ecosystem tools suffice.
- Treat divergent transitions, low effective sample size, high R-hat, and problematic Pareto-k values as visible diagnostics, not silent warnings. Test dimensions, names, indexing, and numerical stability without requiring expensive sampling in routine unit tests.
- Keep Stan code vectorized and numerically stable, avoid unnecessary transformed-parameter storage, and document the statistical meaning and shape of data, parameters, priors, and generated quantities.

## Style and compatibility

- Follow tidyverse-style R code, clear snake_case names, informative errors, and no hidden global state.
- Preserve backward compatibility for exported functions and fitted-object structure; document intentional breaking changes in `NEWS.md`.
- Keep examples and tests fast. Mark tests requiring a compiler or full Stan sampling as optional, and use small deterministic fixtures for the default test suite.
