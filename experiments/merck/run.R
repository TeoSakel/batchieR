#!/usr/bin/env Rscript

merck_process_tree <- function(root_handle) {
    handles <- list(root_handle)
    children <- tryCatch(ps::ps_children(root_handle, recursive = TRUE), error = function(e) list())
    c(handles, children)
}

merck_handle_key <- function(handle) {
    pid <- tryCatch(ps::ps_pid(handle), error = function(e) NA_integer_)
    created <- tryCatch(ps::ps_create_time(handle), error = function(e) NA_real_)
    paste(pid, created, sep = "@")
}

merck_process_rss <- function(handle) {
    tryCatch(as.numeric(ps::ps_memory_info(handle)[["rss"]]), error = function(e) 0)
}

merck_monitor_callr <- function(process, resource_path, grace_seconds = 5) {
    merck_require(c("ps", "callr"))
    root <- ps::ps_handle(process$get_pid())
    controller <- ps::ps_handle(Sys.getpid())
    known <- list()
    samples <- list()
    sample_index <- 0L
    repeat {
        handles <- merck_process_tree(root)
        for (handle in handles) known[[merck_handle_key(handle)]] <- handle
        all_handles <- c(list(controller), handles)
        sample_index <- sample_index + 1L
        samples[[sample_index]] <- data.frame(
            timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
            elapsed_seconds = NA_real_,
            n_processes = length(all_handles),
            total_rss_bytes = sum(vapply(all_handles, merck_process_rss, numeric(1)))
        )
        if (!process$is_alive()) break
        Sys.sleep(0.25)
    }
    started <- as.POSIXct(samples[[1L]]$timestamp, tz = "UTC")
    resources <- do.call(rbind, samples)
    resources$elapsed_seconds <- as.numeric(
        as.POSIXct(resources$timestamp, tz = "UTC") - started,
        units = "secs"
    )
    process$wait(timeout = 1000)
    deadline <- Sys.time() + grace_seconds
    still_running <- function() {
        Filter(function(handle) {
            tryCatch(ps::ps_is_running(handle), error = function(e) FALSE)
        }, known)
    }
    orphans <- still_running()
    while (length(orphans) && Sys.time() < deadline) {
        Sys.sleep(0.1)
        orphans <- still_running()
    }
    orphan_pids <- if (length(orphans)) {
        vapply(orphans, ps::ps_pid, integer(1))
    } else {
        integer()
    }
    if (length(orphans)) {
        invisible(lapply(orphans, function(handle) {
            tryCatch(ps::ps_kill(handle), error = function(e) NULL)
        }))
    }
    remote_error <- NULL
    tryCatch(
        invisible(process$get_result()),
        error = function(error) {
            remote_error <<- conditionMessage(error)
        }
    )
    merck_atomic_write_csv(resources, resource_path)
    list(
        exit_status = process$get_exit_status(),
        peak_rss_bytes = max(resources$total_rss_bytes),
        orphan_pids = orphan_pids,
        cleanup_ok = !length(orphan_pids),
        remote_error = remote_error
    )
}

merck_expected_shape <- function(dataset) {
    actual <- c(
        n_training = nrow(dataset$training),
        n_holdout = nrow(dataset$holdout),
        n_plates = length(unique(dataset$training$plate)),
        n_initial = sum(dataset$training$observed),
        n_treatments = dataset$metadata$n_treatments
    )
    expected <- c(
        n_training = 279832L,
        n_holdout = 31545L,
        n_plates = 702L,
        n_initial = 133L,
        n_treatments = 168L
    )
    if (!identical(as.integer(actual), as.integer(expected))) {
        stop(
            "Prepared Merck shape does not match the pinned workflow: ",
            paste(names(actual), actual, sep = "=", collapse = ", "),
            call. = FALSE
        )
    }
    invisible(actual)
}

merck_metrics_frame <- function(results) {
    rows <- lapply(results, function(result) {
        metric <- result$metrics
        as.data.frame(metric, stringsAsFactors = FALSE)
    })
    if (length(rows)) do.call(rbind, rows) else data.frame()
}

merck_package_fingerprint <- function(package_root) {
    files <- c(
        list.files(file.path(package_root, "R"), full.names = TRUE),
        file.path(package_root, c("DESCRIPTION", "NAMESPACE"))
    )
    hashes <- unname(tools::md5sum(files))
    temporary <- tempfile("batchieR-source-hashes-")
    on.exit(unlink(temporary), add = TRUE)
    writeLines(hashes, temporary)
    unname(tools::md5sum(temporary))
}

merck_ensure_package <- function(package_root, library) {
    dir.create(library, recursive = TRUE, showWarnings = FALSE)
    marker <- file.path(library, ".batchieR-source-fingerprint")
    fingerprint <- merck_package_fingerprint(package_root)
    installed <- file.exists(file.path(library, "batchieR", "DESCRIPTION"))
    current <- if (file.exists(marker)) readLines(marker, warn = FALSE) else character()
    if (!installed || !identical(current, fingerprint)) {
        message("Installing the current batchieR source into the experiment library")
        status <- system2(
            file.path(R.home("bin"), "R"),
            c(
                "CMD", "INSTALL", "--no-multiarch", "--with-keep.source",
                paste0("--library=", shQuote(library)),
                shQuote(package_root)
            )
        )
        if (!identical(status, 0L)) stop("Could not install batchieR for the experiment", call. = FALSE)
        writeLines(fingerprint, marker)
    }
    normalizePath(library)
}

merck_run_experiment <- function(args, script_dir, package_root) {
    merck_require(c("callr", "ps", "jsonlite"))
    profiles <- merck_profiles()
    profile_name <- args$profile %||% "smoke"
    if (!profile_name %in% names(profiles)) stop("Unknown profile: ", profile_name, call. = FALSE)
    if (is.null(args$data) || is.null(args$output)) {
        stop("run requires --data and --output", call. = FALSE)
    }
    dataset <- readRDS(args$data)
    merck_expected_shape(dataset)
    package_library <- merck_ensure_package(
        package_root,
        file.path(dirname(normalizePath(args$data)), "rlib")
    )
    profile <- profiles[[profile_name]]
    output <- normalizePath(args$output, mustWork = FALSE)
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    existing <- sort(list.files(output, pattern = "^round-[0-9]+$", full.names = TRUE))
    if (length(existing) && !isTRUE(args$resume)) {
        stop("Output already contains rounds; use --resume", call. = FALSE)
    }
    completed <- list()
    revealed <- dataset$metadata$initial_plate
    for (directory in existing) {
        path <- file.path(directory, "result.rds")
        if (!file.exists(path)) next
        result <- readRDS(path)
        completed[[length(completed) + 1L]] <- result
        revealed <- result$revealed_after
    }
    manifest_path <- file.path(output, "manifest.json")
    data_md5 <- unname(tools::md5sum(args$data))
    source_fingerprint <- merck_package_fingerprint(package_root)
    if (!file.exists(manifest_path)) {
        manifest <- list(
            started_at = format(Sys.time(), tz = "UTC", usetz = TRUE),
            profile = profile_name,
            profile_settings = profile,
            batch_size = MERCK_BATCH_SIZE,
            scoring_seed = MERCK_SCORING_SEED,
            max_triplets = MERCK_MAX_TRIPLETS,
            data_md5 = data_md5,
            package_source_fingerprint = source_fingerprint,
            package_revision = merck_git_revision(package_root),
            package_version = unname(read.dcf(
                file.path(package_root, "DESCRIPTION"), fields = "Version"
            )[[1L]]),
            R_version = R.version.string,
            data_metadata = dataset$metadata
        )
        merck_atomic_write_json(manifest, manifest_path)
    } else {
        manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
        if (!identical(manifest$profile, profile_name) ||
                !identical(manifest$data_md5, data_md5) ||
                !identical(manifest$package_source_fingerprint, source_fingerprint)) {
            stop("Resume manifest does not match profile, data, or package source", call. = FALSE)
        }
    }
    candidate_total <- length(setdiff(unique(
        merck_subset_profile(dataset, profile)$training$plate
    ), dataset$metadata$initial_plate))
    total_rounds <- min(profile$max_rounds, ceiling(candidate_total / MERCK_BATCH_SIZE))
    if (length(completed) >= total_rounds) {
        message("All ", total_rounds, " configured rounds are already complete")
        return(invisible(completed))
    }
    for (round in seq.int(length(completed) + 1L, total_rounds)) {
        round_dir <- file.path(output, sprintf("round-%03d", round))
        dir.create(round_dir, recursive = TRUE, showWarnings = FALSE)
        result_path <- file.path(round_dir, "result.rds")
        worker_args <- list(
            script_dir = script_dir,
            package_root = package_root,
            package_library = package_library,
            data_path = normalizePath(args$data),
            result_path = result_path,
            profile = profile,
            round = round,
            revealed_plates = revealed,
            fit_seed = 1000L + round,
            scoring_seed = MERCK_SCORING_SEED,
            max_triplets = MERCK_MAX_TRIPLETS,
            batch_size = MERCK_BATCH_SIZE
        )
        message("Starting round ", round, " of ", total_rounds)
        process <- callr::r_bg(
            function(worker_path, worker_args) {
                source(worker_path)
                merck_run_round(worker_args)
            },
            args = list(file.path(script_dir, "worker.R"), worker_args),
            stdout = file.path(round_dir, "stdout.log"),
            stderr = file.path(round_dir, "stderr.log"),
            libpath = c(package_library, .libPaths()),
            system_profile = FALSE,
            user_profile = FALSE,
            supervise = TRUE
        )
        monitoring <- merck_monitor_callr(
            process,
            file.path(round_dir, "resources.csv")
        )
        if (!identical(monitoring$exit_status, 0L) ||
                !monitoring$cleanup_ok || !is.null(monitoring$remote_error)) {
            stop(
                "Round ", round, " failed (exit ", monitoring$exit_status,
                "; orphan PIDs: ", paste(monitoring$orphan_pids, collapse = ", "),
                "; remote error: ", monitoring$remote_error %||% "none", ")",
                call. = FALSE
            )
        }
        if (!file.exists(result_path)) stop("Round worker produced no result", call. = FALSE)
        result <- readRDS(result_path)
        result$metrics$peak_process_tree_rss_bytes <- monitoring$peak_rss_bytes
        result$metrics$exit_status <- monitoring$exit_status
        result$metrics$cleanup_ok <- monitoring$cleanup_ok
        merck_atomic_save_rds(result, result_path)
        completed[[length(completed) + 1L]] <- result
        revealed <- result$revealed_after
        merck_atomic_save_rds(
            list(round = round, revealed_plates = revealed),
            file.path(output, "checkpoint.rds")
        )
        merck_atomic_write_csv(
            merck_metrics_frame(completed),
            file.path(output, "round_metrics.csv")
        )
    }
    manifest <- jsonlite::read_json(manifest_path, simplifyVector = TRUE)
    manifest$completed_at <- format(Sys.time(), tz = "UTC", usetz = TRUE)
    manifest$status <- "complete"
    manifest$completed_rounds <- length(completed)
    merck_atomic_write_json(manifest, manifest_path)
    message("Completed ", total_rounds, " round(s)")
    invisible(completed)
}

merck_run_leak_check <- function(args, script_dir, package_root) {
    merck_require(c("callr", "ps", "jsonlite"))
    if (is.null(args$data) || is.null(args$output)) {
        stop("leak-check requires --data and --output", call. = FALSE)
    }
    dataset <- readRDS(args$data)
    merck_expected_shape(dataset)
    package_library <- merck_ensure_package(
        package_root,
        file.path(dirname(normalizePath(args$data)), "rlib")
    )
    output <- normalizePath(args$output, mustWork = FALSE)
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    result_path <- file.path(output, "leak_result.rds")
    repeats <- as.integer(args$repeats %||% 3L)
    worker_args <- list(
        script_dir = script_dir,
        package_root = package_root,
        package_library = package_library,
        data_path = normalizePath(args$data),
        result_path = result_path,
        profile = merck_profiles()$smoke,
        repeats = repeats
    )
    process <- callr::r_bg(
        function(worker_path, worker_args) {
            source(worker_path)
            merck_run_leak_probe(worker_args)
        },
        args = list(file.path(script_dir, "leak_worker.R"), worker_args),
        stdout = file.path(output, "stdout.log"),
        stderr = file.path(output, "stderr.log"),
        libpath = c(package_library, .libPaths()),
        system_profile = FALSE,
        user_profile = FALSE,
        supervise = TRUE
    )
    monitoring <- merck_monitor_callr(process, file.path(output, "resources.csv"))
    if (!identical(monitoring$exit_status, 0L) ||
            !monitoring$cleanup_ok || !is.null(monitoring$remote_error)) {
        stop(
            "Leak probe failed or left worker processes: ",
            monitoring$remote_error %||% "no remote error",
            call. = FALSE
        )
    }
    result <- readRDS(result_path)
    result$peak_process_tree_rss_bytes <- monitoring$peak_rss_bytes
    result$cleanup_ok <- monitoring$cleanup_ok
    merck_atomic_save_rds(result, result_path)
    merck_atomic_write_csv(result$metrics, file.path(output, "leak_metrics.csv"))
    merck_atomic_write_json(
        list(
            suspected_retention = result$suspected_retention,
            cleanup_ok = result$cleanup_ok,
            peak_process_tree_rss_bytes = result$peak_process_tree_rss_bytes,
            rule = result$rule
        ),
        file.path(output, "leak_summary.json")
    )
    message(
        "Leak probe complete; suspected retention: ",
        result$suspected_retention,
        "; cleanup: ", result$cleanup_ok
    )
    invisible(result)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

if (sys.nframe() == 0L) {
    script_dir <- local({
        command <- commandArgs(trailingOnly = FALSE)
        file_arg <- grep("^--file=", command, value = TRUE)
        dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
    })
    source(file.path(script_dir, "common.R"))
    source(file.path(script_dir, "config.R"))
    values <- commandArgs(trailingOnly = TRUE)
    if (!length(values) || !values[[1L]] %in% c("run", "leak-check")) {
        stop(
            "Usage: run.R run|leak-check --data FILE --output DIR [--profile NAME] [--resume]",
            call. = FALSE
        )
    }
    command <- values[[1L]]
    args <- merck_parse_args(values[-1L], flags = "resume")
    package_root <- normalizePath(file.path(script_dir, "..", ".."))
    if (command == "run") {
        merck_run_experiment(args, script_dir, package_root)
    } else {
        merck_run_leak_check(args, script_dir, package_root)
    }
}
