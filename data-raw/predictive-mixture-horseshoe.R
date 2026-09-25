# Extract the ten-chain horseshoe example without rerunning MCMC.
# Run from the package root: Rscript data-raw/predictive-mixture-horseshoe.R
args <- commandArgs(trailingOnly = TRUE)
source_dir <- if (length(args)) args[1L] else
    "experiments/merck/results/predictive-mixtures-horseshoe-20260925"
destination <- if (length(args) > 1L) args[2L] else "inst/extdata"
devtools::load_all(quiet = TRUE)
stopifnot("status=complete" %in% readLines(file.path(source_dir, "experiment.status")))
job <- readRDS(file.path(source_dir, "job.rds"))
inputs <- readRDS(file.path(source_dir, "predictive-mixture-inputs-v1.rds"))
fit <- readRDS(file.path(source_dir, "predictive-mixture-fit-v1.rds"))
stopifnot(identical(job$model, job$generating_model),
          identical(fit$model, job$model), identical(fit$input$data, inputs$training),
          identical(sort(unique(fit$chain_id)), 1:10),
          all(table(fit$chain_id) == 800L))
source_sampling <- fit$sampling
keep <- fit$chain_id %in% 1:10
fit$draws <- fit$draws[keep]
fit$chain_id <- fit$chain_id[keep]
fit$draw_id <- fit$draw_id[keep]
fit$iteration <- fit$iteration[keep]
# Omit the large per-training-row cache; predictions are reconstructed from
# unchanged parameter draws. Keep RMSE, raw factors and shrinkage diagnostics.
fit$draws <- lapply(fit$draws, function(draw) {
    draw$fitted_mean <- NULL
    draw
})
fit$sampling$chains <- 10L
fit$sampling$parallel_chains <- 4L
fit$provenance <- list(source_sampling = source_sampling, selected_chains = 1:10,
                      omitted_draw_fields = "fitted_mean")
inputs$fixture_id <- "rank6-horseshoe-prior-ten-chains-v2"
inputs$model <- job$model
inputs$generating_model <- job$generating_model
inputs$settings <- list(generation_seed = job$generation_seed,
                        fixture_chains = 1:10, retained_per_chain = 800L)
inputs$provenance <- fit$provenance
inputs$provenance$source_md5 <- tools::md5sum(file.path(source_dir, c(
    "job.rds", "predictive-mixture-inputs-v1.rds", "predictive-mixture-fit-v1.rds"
)))
write_fixture <- function(value, name) {
    path <- file.path(destination, paste0("predictive-mixture-", name, "-v2.rds"))
    tmp <- paste0(path, ".tmp")
    saveRDS(value, tmp, compress = "xz", version = 3)
    if (!file.rename(tmp, path)) cli::cli_abort("Cannot save {.file {path}}.")
}
dir.create(destination, recursive = TRUE, showWarnings = FALSE)
write_fixture(inputs, "inputs")
write_fixture(fit, "fit")
cli::cli_inform("Saved ten chains with 800 draws each as v2 fixtures.")
