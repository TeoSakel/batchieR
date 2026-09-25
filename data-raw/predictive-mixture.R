# Build the fit/data fixtures used by the model-checking and predictive-mixture
# vignettes. This is a post-fit operation: it validates completed posterior
# shards and never runs the sampler or saves derived diagnostics.

args <- commandArgs(trailingOnly = TRUE)
source_job <- if (length(args)) args[[1L]] else
    "experiments/merck/results/full01-replication-ladder/jobs/r04/job.rds"
destination <- if (length(args) >= 2L) args[[2L]] else "inst/extdata"
if (!file.exists(source_job)) {
    stop("Required completed source job is unavailable: ", source_job)
}

write_atomic <- function(value, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- paste0(path, ".tmp")
    saveRDS(value, temporary, compress = "xz", version = 3)
    if (!file.rename(temporary, path)) {
        stop("Cannot finalize fixture: ", path, call. = FALSE)
    }
}
lean_draw <- function(draw) {
    draw <- draw[c("intercept", "beta", "precision", "components")]
    draw$components <- lapply(
        draw$components,
        function(component) if (is.null(component)) NULL else component["value"]
    )
    draw
}

job <- readRDS(source_job)
stopifnot(
    identical(job$id, "replicates-04"),
    identical(job$replicates, 4L),
    identical(job$model$rank, 6L),
    identical(job$config$chains, 10L),
    identical(job$config$retained_per_chain, 800L),
    nrow(job$data) == 26496L,
    identical(job$data, job$compiled$data)
)

draws <- list()
shard_hashes <- character()
for (chain in seq_len(job$config$chains)) {
    directory <- file.path(
        dirname(source_job), "chains", sprintf("chain-%02d", chain)
    )
    checkpoint_path <- file.path(directory, "checkpoint.rds")
    if (!file.exists(checkpoint_path)) {
        stop("Missing checkpoint: ", checkpoint_path, call. = FALSE)
    }
    checkpoint <- readRDS(checkpoint_path)
    stopifnot(
        isTRUE(checkpoint$complete),
        identical(checkpoint$fingerprint, job$fingerprint),
        checkpoint$retained_draws == 800L,
        checkpoint$chain == chain,
        !anyDuplicated(checkpoint$shards)
    )
    chain_draws <- list()
    for (shard in checkpoint$shards) {
        path <- file.path(directory, shard)
        value <- readRDS(path)
        stopifnot(
            identical(value$fingerprint, job$fingerprint),
            value$chain == chain
        )
        shard_hashes[paste0("chain-", chain, "/", shard)] <-
            unname(tools::md5sum(path))
        chain_draws <- c(chain_draws, lapply(value$draws, lean_draw))
    }
    stopifnot(length(chain_draws) == 800L)
    draws <- c(draws, chain_draws)
}

full_fit <- structure(list(
    model = job$model,
    compiled = job$compiled,
    draws = draws,
    chain_id = rep(seq_len(10L), each = 800L),
    draw_id = rep(seq_len(800L), times = 10L),
    iteration = rep(2000L + 40L * seq_len(800L), times = 10L),
    engine = "gibbs",
    sampling = list(
        chains = 10L, iter_warmup = 2000L, iter_sampling = 32000L,
        thin = 40L, seed = job$seed, parallel_chains = 4L
    ),
    control = list(),
    input = list(data = job$data, cell_data = NULL, compound_data = NULL),
    implementation = if (job$compiled$fast_iid) "gibbs_iid" else "gibbs_sparse"
), class = "combo_fit")
subset_fit <- function(object, keep) {
    result <- object
    result$draws <- result$draws[keep]
    result$chain_id <- result$chain_id[keep]
    result$draw_id <- result$draw_id[keep]
    result$iteration <- result$iteration[keep]
    result$sampling$chains <- length(unique(result$chain_id))
    result
}
selected_chains <- c(1L, 3L, 4L, 6L, 10L)
keep <- full_fit$chain_id %in% selected_chains & full_fit$draw_id %% 2L == 0L
fit <- subset_fit(full_fit, keep)
fit$draw_id <- ave(fit$draw_id, fit$chain_id, FUN = seq_along)
long_fit <- subset_fit(full_fit, full_fit$chain_id %in% c(1L, 10L))

design <- job$reference_data[setdiff(names(job$reference_data), "response")]
truth_mean <- job$truth_mean[seq_len(nrow(design))]
old_kind <- RNGkind()
on.exit(do.call(RNGkind, as.list(old_kind)), add = TRUE)
RNGkind("L'Ecuyer-CMRG", "Inversion", "Rejection")
set.seed(20260924L)
holdout <- design
holdout$response <- rnorm(
    nrow(design), truth_mean, 1 / sqrt(job$truth$precision)
)
inputs <- list(
    format_version = 1L,
    fixture_id = "rank6-horseshoe-bounded-truth-v1",
    design = design,
    training = job$data,
    holdout = holdout,
    truth = job$truth,
    truth_mean = truth_mean,
    model = job$model,
    settings = list(
        transform = "plogis", blocks = 4L,
        holdout_seed = 20260924L, predictive_seed = 20260925L,
        fixture_chains = selected_chains,
        retained_per_chain = 400L,
        source_retained_per_chain = 800L
    ),
    provenance = list(
        job_fingerprint = job$fingerprint,
        source_job_fingerprint = job$source_job_fingerprint,
        shard_hashes = shard_hashes,
        sampling = full_fit$sampling,
        chain_seeds = job$seeds,
        pedagogical_chain_selection = list(
            core = c(1L, 3L, 4L), departures = c(6L, 10L),
            selected_without_holdout = TRUE
        ),
        rng_kind = RNGkind(),
        empirical_intercept = mean(job$data$response),
        provisional = FALSE,
        previous_holdout_used = FALSE
    )
)

write_atomic(
    inputs,
    file.path(destination, "predictive-mixture-inputs-v1.rds")
)
write_atomic(
    fit,
    file.path(destination, "predictive-mixture-fit-v1.rds")
)
write_atomic(
    long_fit,
    file.path(destination, "predictive-mixture-long-fit-v1.rds")
)
message("Wrote synthetic-data and combo_fit fixtures; no derived analysis saved.")
