# Sampler benchmark

This development experiment compares sampling algorithms on fixed synthetic
combination screens. It is separate from the Merck active-learning replay:
every candidate receives the same generated responses and observation masks.
The directory is excluded from R package builds; the package's sampling API,
priors, initialization, and RNG behavior are unchanged.

## Run and inspect

Run from the repository root using its restored `renv` environment. The existing
development dependencies include `devtools`, `posterior`, `ggplot2`, `knitr`, and
`quarto`; the Quarto command-line tool must also be installed. Rendering
dependencies are checked before sampling starts.

```sh
Rscript experiments/samplers/run.R run \
  --profile smoke --output experiments/samplers/results/smoke

Rscript experiments/samplers/run.R run \
  --profile local --output experiments/samplers/results/local

# Resume an interrupted invocation without refitting committed jobs.
Rscript experiments/samplers/run.R run \
  --profile local --output experiments/samplers/results/local --resume

# Rebuild CSV tables and HTML from saved results, without fitting or loading
# the combination model. This also works for partial or failed runs.
Rscript experiments/samplers/run.R report \
  --output experiments/samplers/results/local
```

Open `report.html` in the selected output directory. Images, tables, styles, and
scripts are embedded, so the HTML can be shared as one file. A run is one
invocation over the configured scenarios, independently generated datasets,
observation masks, and candidate samplers. Each report describes the sampler,
effective settings, screen, mixing, recovery, SBC, analytic self-checks,
failures, and source provenance. One-candidate reports explicitly identify
that paired sampler comparisons are unavailable.

| Profile | Design: cells / drugs / doses | Rank | Datasets per scenario | Chains | Warmup / sampling |
| --- | --- | --- | --- | --- | --- |
| `smoke` | 2 / 3 / 2 | 2 | 1 | 2 | 2 / 8 |
| `local` | 4 / 8 / 3 | 4 | 2 | 4 | 250 / 500 |
| `calibration` | 4 / 8 / 3 | 4 | 100 | 4 | 1,000 / 4,000 |

Chains run sequentially, retaining every post-warmup draw (`thin = 1`). All
profiles use both requested plate reveal fractions, 25% and 75%. Smoke and
local include all four scenarios. Calibration includes the three matched-prior
scenarios and is a substantial, explicitly requested computation: 600 fits at
the default settings. Local targets tens of minutes, but elapsed time depends
on hardware and generated posterior geometry. Smoke verifies plumbing and
cannot assess convergence; local SBC results are exploratory.

The CLI accepts positive-integer overrides `--cells`, `--drugs` (at least 2),
`--doses`, `--rank`, `--datasets`, `--chains`, `--iter-sampling`, and `--seed`.
`--iter-warmup` may also be zero. Use `--scenarios` with a comma-separated subset.
No adaptive iteration extension is performed.

```sh
Rscript experiments/samplers/run.R run --profile calibration \
  --scenarios gamma,multiplicative_gamma,horseshoe \
  --datasets 100 --output experiments/samplers/results/calibration

# Reuse an already prepared Merck design. No downloads or upstream replay.
Rscript experiments/samplers/run.R run --profile local \
  --merck experiments/merck/work/merck_batchieR.rds \
  --output experiments/samplers/results/merck-design
```

`--merck` uses all training and holdout design rows from that RDS, so the local
profile no longer implies a small screen. All original responses, retrospective
truth, and original observation masks are discarded. No full-scale Merck run
is part of package tests.

## Synthetic data

The default local design has 1,104 rows: 96 singles plus 1,008 combinations.
Each unordered drug pair has a complete dose grid in every cell. A combination
plate is a cell-by-drug-pair grid; singles have their own plates. Singles are
initially observed, and nested 25%/75% masks reveal additional plates from a
shared random ordering. If imported plates mix singles with combinations,
observing a single reveals that entire plate. At least one candidate plate
remains held out, and actual row counts are reported.

All matched-prior scenarios use categorical doses, fixed intercept
`qlogis(0.8)`, and observation precision `Gamma(shape=25, rate=1)`:

- `gamma`: each Gaussian component uses `Gamma(shape=3, rate=0.3)` precision.
- `multiplicative_gamma`: same model, with package-default multiplicative-gamma
  shrinkage for cell factors.
- `horseshoe`: package-default component priors, including horseshoes and
  multiplicative-gamma cell factors. The fixed intercept and observation prior
  remain as specified above.

The generating model and fitted model are identical in these scenarios. An
internal prior snapshot supplies generating parameters and precisions; an
independently written predictor checks the model algebra against package
prediction. The current prior-prediction API is unchanged. Finite extreme
prior draws are kept, rather than rejecting them based on observed responses.
Non-finite draws are recorded as generation failures, not silently resampled.

`realism` generates rank-2 effects with standard-normal cell factors, cell
offset SD 0.2, main drug-factor SD 0.3, interaction drug-factor SD 0.4, negative
drug offsets `-0.8 * abs(rnorm(1))`, and observation noise SD 0.2. Each drug has
one vector scaled by ordered positive-dose weights `1/D, ..., 1`; compound
effects are therefore shared across doses. The fitted rank follows the profile
(4 locally), with package-default shrinkage. These are controlled synthetic
assumptions, not empirically estimated Merck parameters, and this scenario is
excluded from formal SBC.

Responses are generated on the Gaussian logit-response scale. Viability
summaries use `plogis()` without clipping generated observations. The report
records the fraction near viability boundaries to expose extreme priors.

## Metrics and calibration

Each fit retains unthinned draws for a design-selected panel: noise SD, eight
latent means, eight interaction contributions, and total observed-data
log-likelihood (smaller designs use fewer unique rows). Diagnostics preserve
chain and iteration boundaries. They include rank-normalized R-hat, bulk/tail
ESS, MCSE of the mean, and ESS divided by total fitting time, including
compilation, initialization, and warmup. Raw factor labels are not the primary
efficiency target: relabeling dimensions alone cannot improve these predictor
quantities. Predictor invariance under a column permutation does **not** imply
that the full prior permits that permutation.

Latent and interaction recovery uses all design rows. Held-out response and
viability RMSE use unrevealed plates. 50%, 80%, and 95% intervals report coverage
and average width. Observed viability predictions average `plogis()` of noisy
posterior predictive draws, not `plogis()` of a latent mean. Prediction and
likelihood evaluation run in row chunks; only compact diagnostic draws survive
in fit results. Generating truth is a recovery target, not an exact posterior
reference for estimating Monte Carlo bias.

For matched-prior SBC, each dataset contributes one randomized-tie rank per
quantity. The generating log-likelihood is independently evaluated using the
same observed rows and Gaussian normalization as `log_lik()`. The runner
selects M=100 draws, as evenly divided among chains as possible, using
deterministic temporal spacing. Each quantity needs full-chain R-hat <=1.01,
bulk and tail ESS >=100, and basic ESS >=90% of the equal-length spaced chain
segments. These checks only approximate independence. Short, nonconverged, or
still-correlated results are explicitly ineligible; there is no automatic
extension. Efficiency statistics continue to use unthinned draws.

ECDF difference plots compare eligible ranks to the discrete uniform CDF on
`0:100`. A 95% simultaneous envelope uses 10,000 independent null simulations
and the 95th percentile of the maximum absolute CDF difference over that
support. Its seed is recorded through the master seed and deterministic
`report-envelope` seed derivation. Envelopes with the same dataset count and
rank support are reused. Bands cover one quantity's rank support, not all
quantities simultaneously. Histograms explicitly bin the 101 possible ranks
into ten bins (11 ranks in the first, 10 in the others).

Quantities, scenarios, masks, and samplers are never pooled as independent
replications. Missing fits and ineligible ranks are counted. Any panel based on
an incomplete eligible subset is inconclusive because successful fits may be a
selected subset. Reports do not make automatic calibration-pass claims.

Every fitting invocation also runs an inexpensive analytic self-check: a
normal mean with a proper standard-normal prior and 20 observations of known
variance 1. Exact posterior draws are compared with deliberately incorrect
prior-only, shifted, and narrowed draws. The likelihood statistic detects the
prior-only failure that marginal parameter ranks can miss. These test adapters
are not candidate combination-model samplers. A correct reference can cross a
95% band by chance; that alone is not treated as a self-check failure.

Methods: [Stan SBC guidance](https://mc-stan.org/docs/stan-users-guide/simulation-based-calibration.html)
and [Modrák et al., choice of SBC test quantities](https://arxiv.org/abs/2211.02383).

## Add a candidate sampler

Add an entry to `sb_adapters()` in `config.R`, or supply an adapter list directly
to the experiment's `sb_run()` function. Each adapter has a unique lowercase
`id`, `name`, `description`, and `fit(model, data, settings)` function returning a
standard `combo_fit`. The description must state the update strategy,
initialization, and differences from baseline. The returned fit must honor the
budget and preserve the supplied model and masked data.

The existing adapter delegates directly to `fit_combo()`. Future adapters can
call an experimental implementation without changing the package's public
API. All candidates get the same generating data, masks, and fitting seeds;
generation, masking, fitting, evaluation, ranks, and report envelopes have
separate streams. Matching seeds do not imply identical random variates after
different algorithms consume different numbers of draws.

Candidate implementation code should live in fingerprinted experiment `.R`
files. For programmatic closures, captured configuration must be represented
in adapter metadata/description; the manifest records function source but
cannot identify changes in arbitrary external mutable state.

## Artifacts, resume, and validation

`manifest.rds` records the configuration, adapter descriptions and function
source, source-file MD5 fingerprints, R/package/Quarto versions, and optional
Merck input checksum. `jobs.rds` is the expected job table. `datasets/` retains
truth and seeds; `designs/` holds compact screen summaries. `fits/` contains
atomic per-job results, diagnostic draws, warnings, timing, and failures.
`self-checks.rds` retains the analytic checks. `metrics.csv`, `diagnostics.csv`,
`ranks.csv`, `coverage.csv`, and `status.csv` are regenerated from committed
results. `report-data.rds` is the render input; `report.qmd` and `report.R` are
staged alongside the standalone `report.html`.

Resume requires identical configuration, adapters, source fingerprints, input
checksum, and recorded environment. Committed failures remain failures;
interrupted jobs without committed results are rerun. Use a new output
directory after changing an algorithm or budget. Do not combine result files
from different runs. Report-only regeneration permits a newer renderer and
retains the original sampling provenance; it never refits or reruns analytic
checks. Rendering errors leave numerical results intact and write
`report-error.log`; fix the renderer/environment and run the report command.
User interrupts trigger a best-effort partial report; after a hard process
termination, use the report command manually.

Generated `results/` and `work/` are ignored by Git. Tests use the devtools
pipeline and skip experiment tests in installed/source package builds where
this directory is absent:

```r
devtools::load_all()
devtools::test(filter = "sampler-benchmark")
devtools::test()
```

Focused tests cover generation and RNG preservation, family matching, Merck
design import, complete masks, hidden-response isolation, independent truth,
baseline equivalence, chain ordering, ranks and envelopes, deliberate faults,
resume compatibility, and complete/partial/insufficient-draw HTML reports.

In environments where devtools exports a project-relative startup profile,
existing parallel-chain tests may need an explicit child-process library path.
This preserves the resolved renv libraries while avoiding relative startup
paths in worker directories:

```sh
Rscript -e 'devtools::load_all(); Sys.setenv(R_PROFILE_USER="/dev/null", R_LIBS=paste(.libPaths(), collapse=.Platform$path.sep)); devtools::test(stop_on_failure=TRUE)'
```
