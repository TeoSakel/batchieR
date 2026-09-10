#!/usr/bin/env Rscript
# Run in the repository development environment; load the current checkout.
sb_main <- function(arguments = commandArgs(trailingOnly = TRUE)) {
    script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
    root <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
    environment <- new.env(parent = globalenv())
    for (file in c("config.R", "data.R", "sbc.R", "metrics.R", "report.R", "runner.R")) {
        sys.source(file.path(root, "experiments/samplers", file), envir = environment)
    }
    if (!length(arguments) || !arguments[1L] %in% c("run", "report")) {
        cli::cli_abort("Usage: Rscript experiments/samplers/run.R run|report --output PATH [--profile smoke|local|calibration] [--resume]")
    }
    command <- arguments[1L]
    args <- list()
    i <- 2L
    while (i <= length(arguments)) {
        key <- sub("^--", "", arguments[i])
        if (!startsWith(arguments[i], "--") || key %in% names(args)) cli::cli_abort("Invalid or duplicate option: {arguments[i]}")
        if (key == "resume") {
            args[[key]] <- TRUE
            i <- i + 1L
        } else {
            if (i == length(arguments)) cli::cli_abort("Missing value for {key}.")
            args[[key]] <- arguments[i + 1L]
            i <- i + 2L
        }
    }
    allowed <- c("output", "profile", "seed", "scenarios", "merck", "resume",
                 "cells", "drugs", "doses", "rank", "datasets", "chains", "iter-warmup", "iter-sampling")
    unknown <- setdiff(names(args), allowed)
    if (length(unknown)) cli::cli_abort("Unknown option(s): {.and {unknown}}")
    if (is.null(args$output)) cli::cli_abort("--output is required.")
    if (command == "report") {
        environment$sb_finish_report(root, normalizePath(args$output, mustWork = TRUE))
    } else {
        devtools::load_all(root, quiet = TRUE)
        names(args) <- gsub("-", "_", names(args), fixed = TRUE)
        fields <- intersect(names(args), names(environment$sb_profiles()$local))
        overrides <- lapply(args[fields], as.numeric)
        config <- environment$sb_config(
            profile = if (is.null(args$profile)) "local" else args$profile,
            seed = if (is.null(args$seed)) 20260910L else as.numeric(args$seed),
            overrides = overrides,
            scenarios = if (!is.null(args$scenarios)) strsplit(args$scenarios, ",", fixed = TRUE)[[1L]] else NULL,
            merck = args$merck)
        environment$sb_run(root, args$output, config, resume = isTRUE(args$resume))
    }
}

if (sys.nframe() == 0L) sb_main()
