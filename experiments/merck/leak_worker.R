merck_run_leak_probe <- function(args) {
    source(file.path(args$script_dir, "common.R"))
    library(package = "batchieR", lib.loc = args$package_library, character.only = TRUE)
    options(parallelly.availableCores.methods = "system")
    dataset <- merck_subset_profile(readRDS(args$data_path), args$profile)
    training <- dataset$training
    training$response <- NA_real_
    initial <- dataset$metadata$initial_plate
    training$response[training$plate == initial] <- training$truth[training$plate == initial]
    workers_before <- future::nbrOfWorkers()
    rows <- vector("list", args$repeats)
    for (repeat_index in seq_len(args$repeats)) {
        heap_before <- merck_heap_mb(gc())
        started <- proc.time()[["elapsed"]]
        fit <- fit_combo(
            combo_model(rank = 2L),
            training,
            chains = 2L,
            iter_warmup = 0L,
            iter_sampling = 3L,
            thin = 1L,
            seed = 9000L + repeat_index,
            parallel_chains = 2L,
            refresh = 0L
        )
        elapsed <- proc.time()[["elapsed"]] - started
        fit_bytes <- as.numeric(object.size(fit))
        rm(fit)
        heap_after <- merck_heap_mb(gc())
        workers_after <- future::nbrOfWorkers()
        if (!identical(workers_after, workers_before)) {
            stop("A repeated fit did not restore the previous future plan", call. = FALSE)
        }
        rows[[repeat_index]] <- data.frame(
            repeat_id = repeat_index,
            elapsed_seconds = elapsed,
            fit_bytes = fit_bytes,
            heap_before_mb = heap_before,
            heap_after_gc_mb = heap_after,
            workers_before = workers_before,
            workers_after = workers_after
        )
    }
    metrics <- do.call(rbind, rows)
    retained <- metrics$heap_after_gc_mb
    suspected <- all(diff(retained) > 0) &&
        retained[[length(retained)]] - retained[[1L]] > 50 &&
        retained[[length(retained)]] > 1.1 * retained[[1L]]
    result <- list(
        metrics = metrics,
        suspected_retention = suspected,
        rule = "monotonic post-GC growth exceeding both 50 MiB and 10 percent"
    )
    merck_atomic_save_rds(result, args$result_path)
    invisible(result)
}
