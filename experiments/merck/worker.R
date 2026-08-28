merck_run_round <- function(args) {
    source(file.path(args$script_dir, "common.R"))
    library(package = "batchieR", lib.loc = args$package_library, character.only = TRUE)
    options(parallelly.availableCores.methods = "system")
    dataset <- readRDS(args$data_path)
    dataset <- merck_subset_profile(dataset, args$profile)
    training <- dataset$training
    revealed <- unique(args$revealed_plates)
    training$response <- NA_real_
    training$response[training$plate %in% revealed] <-
        training$truth[training$plate %in% revealed]
    candidates <- sort(setdiff(unique(training$plate), revealed))
    if (!length(candidates)) stop("No unobserved candidate plates remain", call. = FALSE)

    model <- combo_model(rank = args$profile$rank)
    workers_before <- future::nbrOfWorkers()
    heap_before <- merck_heap_mb(gc())
    fit_start <- proc.time()[["elapsed"]]
    fit <- fit_combo(
        model,
        training,
        chains = args$profile$chains,
        iter_warmup = args$profile$iter_warmup,
        iter_sampling = args$profile$iter_sampling,
        thin = args$profile$thin,
        seed = args$fit_seed,
        parallel_chains = args$profile$parallel_chains,
        refresh = 0L
    )
    fit_seconds <- proc.time()[["elapsed"]] - fit_start
    workers_after_fit <- future::nbrOfWorkers()
    if (!identical(workers_after_fit, workers_before)) {
        stop("fit_combo() did not restore the previous future plan", call. = FALSE)
    }

    evaluation_start <- proc.time()[["elapsed"]]
    holdout_epred <- posterior_epred(fit, newdata = dataset$holdout)
    holdout_prediction <- colMeans(stats::plogis(holdout_epred))
    holdout_mse <- mean((holdout_prediction - dataset$holdout$viability)^2)
    evaluation_seconds <- proc.time()[["elapsed"]] - evaluation_start

    candidate_data <- lapply(candidates, function(plate) {
        training[training$plate == plate, , drop = FALSE]
    })
    names(candidate_data) <- candidates
    scoring_start <- proc.time()[["elapsed"]]
    set.seed(args$scoring_seed)
    scores <- score_plate_pdbal(
        fit,
        candidate_data = candidate_data,
        reference_grid = training,
        response_transform = stats::plogis,
        max_triplets = args$max_triplets
    )
    scoring_seconds <- proc.time()[["elapsed"]] - scoring_start
    selected <- head(scores$plate, min(args$batch_size, nrow(scores)))
    fit_bytes <- as.numeric(object.size(fit))
    data_bytes <- as.numeric(object.size(training))
    fit_summary <- summary(fit)
    rm(fit, holdout_epred, holdout_prediction, candidate_data)
    heap_after <- merck_heap_mb(gc())

    result <- list(
        round = args$round,
        revealed_before = revealed,
        selected_plates = selected,
        revealed_after = unique(c(revealed, selected)),
        scores = scores,
        metrics = list(
            round = args$round,
            n_observed = sum(!is.na(training$response)),
            n_candidates = length(candidates),
            n_selected = length(selected),
            holdout_mse = holdout_mse,
            holdout_rmse = sqrt(holdout_mse),
            sampler_rmse_median = unname(fit_summary$rmse[[2L]]),
            fit_seconds = fit_seconds,
            evaluation_seconds = evaluation_seconds,
            scoring_seconds = scoring_seconds,
            fit_bytes = fit_bytes,
            data_bytes = data_bytes,
            heap_before_mb = heap_before,
            heap_after_gc_mb = heap_after,
            workers_before = workers_before,
            workers_after_fit = workers_after_fit,
            selected_plates = paste(selected, collapse = ";")
        )
    )
    merck_atomic_save_rds(result, args$result_path)
    invisible(result)
}
