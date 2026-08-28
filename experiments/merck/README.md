# Merck retrospective benchmark

This directory reproduces the Merck retrospective design from BatchIE v0.0.1
with the batchieR model and scorer. It is development infrastructure: the whole
`experiments/` tree is excluded from R package builds, and its dependencies are
not package dependencies.

The preparation is pinned to the Zenodo BatchIE source archive at commit
`6baa258cb5d77430dccab38b73b0908f91503c91`. Downloads are checked against the
published MD5 sums. Python 3.11 or newer is required because the archived
package declares that requirement. The upstream environment installs PyTorch
and Pyro and can therefore take substantial time and disk space.

## Prepare the data

From the package root, restore the R development environment and run:

```sh
Rscript -e 'renv::restore()'
python3.11 experiments/merck/prepare.py
```

The preparation command downloads the 3 MB screen and archived source, creates
`experiments/merck/work/.venv`, runs the original pairwise plate generator,
smoothers, sparse-cover initialization and 10% holdout, and converts the two
HDF5 screens to `work/merck_batchieR.rds`. A valid pinned result contains
279,832 training rows, 31,545 holdout rows, 702 plates, 133 initially observed
rows, and 168 non-control drug-dose entities.

To convert already prepared upstream screens directly:

```sh
Rscript --vanilla experiments/merck/convert.R \
  --training path/to/training.screen.h5 \
  --holdout path/to/test.screen.h5 \
  --output experiments/merck/work/merck_batchieR.rds
```

## Generated directory structure

All generated files live below `work/` or `results/`. Both directories are
ignored by Git, and the entire `experiments/` tree is excluded from package
builds. Nothing in these directories is required to install or load batchieR.

### Preparation workspace

A complete preparation has the following layout:

```text
work/
|-- .venv/                         # local Python environment for BatchIE v0.0.1
|-- batchie-v0.0.1-zenodo.zip      # checksum-verified upstream source archive
|-- merck_2016.screen.h5            # checksum-verified original Merck screen
|-- merck_training.screen.h5        # upstream prepared training screen
|-- merck_holdout.screen.h5         # upstream prepared 10% holdout screen
|-- preparation.json                # upstream revision, checksums, seed and date
|-- merck_batchieR.rds              # converted input used by the R benchmarks
`-- rlib/                            # local R library used by benchmark workers
```

`merck_batchieR.rds` contains the training and holdout tables plus their
provenance and dimensions. Masked training responses remain hidden in
`response`, while `truth` retains the retrospective value revealed when a plate
is selected. The controller installs the current batchieR source into `rlib/`
when its source fingerprint changes; `.batchieR-source-fingerprint` inside that
library records the installed version.

The downloads, Python environment, prepared HDF5 files and local R library can
all be regenerated. Keeping `merck_batchieR.rds` is sufficient to rerun a
benchmark without repeating the Python preparation. Removing `rlib/` only
causes the controller to reinstall the current package source on the next run.

## Run a replay

```sh
Rscript --vanilla experiments/merck/run.R run \
  --data experiments/merck/work/merck_batchieR.rds \
  --output experiments/merck/results/smoke \
  --profile smoke
```

Use `--resume` to continue an existing output directory. Completed round files
are authoritative; an interrupted round is rerun. The controller installs the
current package source into an ignored experiment library so multisession
workers run the exact source being benchmarked.

Profiles are fixed in `config.R`:

| Profile | Data and rounds | Sampler |
| --- | --- | --- |
| `smoke` | Initial plate plus 9 candidates; 1 round | Rank 2, 1 chain, 2 warmup, 3 retained draws |
| `scaled` | All plates; 3 rounds | Rank 12, 2 chains, 10 warmup, 20 sampling, thin 2 |
| `full` | All 702 plates; up to 78 rounds | Published rank 12, 2 chains, 2,000 warmup, 8,000 sampling, thin 40 |

All profiles select batches of nine, use scoring seed 12, compare posterior
responses on the viability scale, and cap PDBAL at 5,000 triplets. The full
profile is an HPC-scale experiment and is deliberately not run by package tests.

Each output contains a provenance manifest, resumable checkpoint, per-round
RDS and CSV summaries, stdout/stderr logs, and 250 ms process-tree RSS samples.
Fitted model objects are discarded after evaluation and scoring.

### Replay results

The value passed to `--output` names one independent run. For example,
`--output experiments/merck/results/scaled` produces:

```text
results/scaled/
|-- manifest.json
|-- checkpoint.rds
|-- round_metrics.csv
|-- round-001/
|   |-- result.rds
|   |-- resources.csv
|   |-- stdout.log
|   `-- stderr.log
|-- round-002/
|   `-- ...
`-- round-NNN/
    `-- ...
```

The files have the following roles:

- `manifest.json` identifies the profile, sampler and scoring configuration,
  input checksum, batchieR source revision and fingerprint, package and R
  versions, and prepared-data metadata. A completed run also records its
  completion time, status and number of completed rounds. Resume refuses to
  mix a different profile, input dataset or package source with an existing
  manifest.
- `checkpoint.rds` is the compact latest-round pointer and the accumulated set
  of revealed plates. It is written atomically after a successful round.
- `round_metrics.csv` is the human-readable, cumulative round summary. It
  includes holdout error, sampler diagnostics, timings, retained object sizes,
  heap measurements, selected plates, peak process-tree RSS, child exit status
  and process-cleanup status.
- `round-NNN/result.rds` is the authoritative completed-round record. It holds
  the selected and revealed plates, the candidate score table and that round's
  metrics. It deliberately does not retain the fitted model or posterior
  prediction matrices.
- `round-NNN/resources.csv` samples elapsed time, process count and total RSS
  for the controller, isolated round process and its worker descendants every
  250 ms.
- `round-NNN/stdout.log` and `stderr.log` capture output from the isolated
  round process. Empty logs are normal for a successful run with progress
  disabled.

The columns in `round_metrics.csv` are:

- `round`: one-based retrospective round number.
- `n_observed`: number of training rows whose responses were revealed before
  fitting that round.
- `n_candidates`: number of unrevealed candidate plates scored in that round.
- `n_selected`: number of plates selected for the next reveal, normally nine.
- `holdout_mse`: mean squared error of the posterior mean prediction against
  the held-out viability values. Posterior responses are transformed with
  `plogis()` before averaging.
- `holdout_rmse`: square root of `holdout_mse`.
- `sampler_rmse_median`: median training-response RMSE across retained Gibbs
  draws, on the model's logit-response scale.
- `fit_seconds`: elapsed wall time spent in `fit_combo()`.
- `evaluation_seconds`: elapsed wall time spent predicting and evaluating the holdout table.
- `scoring_seconds`: elapsed wall time spent scoring all candidate plates with PDBAL.
- `fit_bytes`: `object.size()` of the retained `combo_fit` object before it is discarded.
- `data_bytes`: `object.size()` of the round's prepared training table.
- `heap_before_mb`: R heap use reported by `gc()` immediately before fitting, in MiB.
- `heap_after_gc_mb`: R heap use after predictions and the fit are removed and
  garbage collection runs, in MiB.
- `workers_before`: worker count reported by the active `future` plan before fitting.
- `workers_after_fit`: worker count after `fit_combo()` returns; it should equal
  `workers_before`, confirming restoration of the caller's future plan.
- `selected_plates`: semicolon-separated plate identifiers, in descending selection priority.
- `peak_process_tree_rss_bytes`: maximum sampled combined RSS of the controller,
  isolated round process, and its worker descendants, in bytes.
- `exit_status`: isolated round-process exit status; zero indicates success.
- `cleanup_ok`: whether all observed descendant worker processes exited within
  the five-second cleanup grace period.

Round results, the checkpoint and cumulative metrics are written atomically.
On `--resume`, every `round-NNN/result.rds` is treated as committed; a round
directory without that file is treated as interrupted and rerun. Consequently,
do not combine or copy individual round directories between runs. To restart
from scratch, choose a new output name or remove the entire old run directory.

## Memory and cleanup probe

```sh
Rscript --vanilla experiments/merck/run.R leak-check \
  --data experiments/merck/work/merck_batchieR.rds \
  --output experiments/merck/results/leak \
  --repeats 3
```

This repeats an identical two-chain fit in one isolated process. Surviving
workers or remote errors fail the command. Post-GC heap growth is reported as
suspected retention only when it is monotonic and exceeds both 50 MiB and 10%;
that diagnostic does not fail the command because R may retain allocated pages.

A leak-probe output directory is flatter than a replay output:

```text
results/leak/
|-- leak_result.rds
|-- leak_metrics.csv
|-- leak_summary.json
|-- resources.csv
|-- stdout.log
`-- stderr.log
```

`leak_metrics.csv` records elapsed time, fit size, post-GC heap use and worker
counts for every repeat. `leak_summary.json` reports the retention flag, cleanup
status, peak process-tree RSS and the diagnostic rule; `leak_result.rds` retains
the same information in its native R representation. The resource samples and
logs have the same meaning as for replay rounds.

Run the converter self-test with:

```sh
Rscript --vanilla experiments/merck/test-harness.R
```
