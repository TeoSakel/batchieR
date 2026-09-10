sb_config <- function(profile = "local", seed = 20260910L, overrides = list(), scenarios = NULL, merck = NULL) {
    config <- sb_profiles()[[profile]]
    if (is.null(config)) cli::cli_abort("Unknown profile: {profile}")
    unknown <- setdiff(names(overrides), names(config))
    if (length(unknown)) cli::cli_abort("Unknown settings: {.and {unknown}}")
    config[names(overrides)] <- overrides
    for (name in names(config)) {
        x <- config[[name]]
        lower <- if (name == "iter_warmup") 0 else 1
        if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x) || x < lower || x > .Machine$integer.max) {
            cli::cli_abort("{name} must be an integer >= {lower}.")
        }
        config[[name]] <- as.integer(x)
    }
    if (config$drugs < 2L) cli::cli_abort("At least two drugs are required.")
    if (length(seed) != 1L || !is.finite(seed) || seed < 1 || seed != floor(seed) || seed > .Machine$integer.max) {
        cli::cli_abort("Seed must be a positive integer.")
    }
    if (is.null(scenarios)) scenarios <- if (profile == "calibration") head(sb_scenarios(), 3L) else sb_scenarios()
    if (!length(scenarios) || anyNA(scenarios) || anyDuplicated(scenarios) || any(!scenarios %in% sb_scenarios())) {
        cli::cli_abort("Select unique known scenarios.")
    }
    c(config, list(profile = profile, seed = as.integer(seed), scenarios = scenarios,
                   merck = if (!is.null(merck)) normalizePath(merck, mustWork = TRUE) else NULL))
}

sb_preflight <- function() {
    for (package in c("devtools", "posterior", "ggplot2", "knitr", "quarto")) {
        if (!requireNamespace(package, quietly = TRUE)) cli::cli_abort("Missing development dependency: {package}")
    }
    if (is.null(quarto::quarto_path())) cli::cli_abort("Quarto CLI is required to render run reports.")
}

sb_fingerprint <- function(root) {
    files <- sort(c(list.files(file.path(root, "R"), full.names = TRUE, pattern = "\\.R$"),
                    file.path(root, c("DESCRIPTION", "NAMESPACE", "renv.lock")),
                    list.files(file.path(root, "experiments/samplers"), full.names = TRUE,
                               pattern = "\\.(R|qmd)$")))
    stats::setNames(unname(tools::md5sum(files)), substring(files, nchar(root) + 2L))
}

sb_manifest <- function(root, config, adapters) {
    adapter_metadata <- lapply(adapters, function(a) {
        list(id = a$id, name = a$name, description = a$description,
             function_source = deparse(a$fit))
    })
    versions <- vapply(c("batchieR", "posterior", "ggplot2", "knitr", "quarto"),
                       function(p) as.character(utils::packageVersion(p)), character(1))
    list(format_version = 1L, config = config, adapters = adapter_metadata,
         source = sb_fingerprint(root), versions = versions, R = R.version.string,
         quarto = as.character(quarto::quarto_version()),
         platform = R.version$platform, rng = c("Mersenne-Twister", "Inversion", "Rejection"),
         merck_md5 = if (!is.null(config$merck)) unname(tools::md5sum(config$merck)) else NULL)
}

sb_fit_job <- function(dataset, mask, adapter, settings, evaluation_seed, rank_seed) {
    warnings <- character()
    started <- proc.time()[["elapsed"]]
    result <- tryCatch(withCallingHandlers({
        if (!is.null(dataset$error)) cli::cli_abort("Generation failed: {dataset$error}")
        fit_started <- proc.time()[["elapsed"]]
        fit <- adapter$fit(dataset$model, sb_fit_data(dataset, mask), settings)
        fit_seconds <- proc.time()[["elapsed"]] - fit_started
        if (!inherits(fit, "combo_fit")) cli::cli_abort("Adapter must return a combo_fit.")
        if (!identical(fit$sampling$chains, settings$chains) ||
            !identical(fit$sampling$iter_sampling, settings$iter_sampling) ||
            !identical(fit$sampling$iter_warmup, settings$iter_warmup) ||
            !identical(fit$sampling$thin, settings$thin)) {
            cli::cli_abort("Adapter did not honor the requested sampling budget.")
        }
        evaluation_started <- proc.time()[["elapsed"]]
        value <- sb_evaluate(fit, dataset, mask, fit_seconds, evaluation_seed, rank_seed)
        value$metrics$evaluation_seconds <- proc.time()[["elapsed"]] - evaluation_started
        c(list(status = "ok", settings = fit$sampling), value)
    }, warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
    }), error = function(e) list(status = "failed", error = conditionMessage(e)))
    c(result, list(warnings = unique(warnings), seconds = proc.time()[["elapsed"]] - started,
                   requested_settings = settings, evaluation_seed = evaluation_seed, rank_seed = rank_seed))
}

sb_dataset_summary <- function(dataset) {
    if (!is.null(dataset$error)) return(list(error = dataset$error))
    list(n_rows = nrow(dataset$design), cells = length(unique(dataset$design$cell)),
         plates = length(unique(dataset$design$plate)),
         treatments = length(unique(c(dataset$indices$treatment_1, dataset$indices$treatment_2)[
             c(dataset$indices$treatment_1, dataset$indices$treatment_2) > 0L])),
         response_quantiles = stats::quantile(dataset$response, c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1)),
         viability_quantiles = stats::quantile(stats::plogis(dataset$response), c(0, 0.01, 0.25, 0.5, 0.75, 0.99, 1)),
         fraction_near_boundary = mean(stats::plogis(dataset$response) < 0.01 | stats::plogis(dataset$response) > 0.99),
         observed = vapply(dataset$masks, sum, integer(1)), seed = dataset$seed, mask_seed = dataset$mask_seed)
}

sb_run <- function(root, output, config = sb_config(), adapters = sb_adapters(), resume = FALSE) {
    root <- normalizePath(root, mustWork = TRUE)
    sb_validate_adapters(adapters)
    sb_preflight()
    manifest <- sb_manifest(root, config, adapters)
    dir.create(output, recursive = TRUE, showWarnings = FALSE)
    output <- normalizePath(output, mustWork = TRUE)
    manifest_path <- file.path(output, "manifest.rds")
    if (file.exists(manifest_path)) {
        if (!resume) cli::cli_abort("Output already contains a run; use --resume or a new directory.")
        previous <- readRDS(manifest_path)
        if (!identical(previous$specification, manifest)) {
            cli::cli_abort("Resume refused: configuration, adapters, source, input, or environment changed.")
        }
    } else {
        if (length(list.files(output, all.files = TRUE, no.. = TRUE))) cli::cli_abort("New run output must be empty.")
        previous <- list(specification = manifest, started = format(Sys.time(), tz = "UTC", usetz = TRUE))
        sb_save(previous, manifest_path)
    }
    jobs <- expand.grid(scenario = config$scenarios, dataset = seq_len(config$datasets),
                        mask = c("reveal_25", "reveal_75"), adapter = names(adapters),
                        stringsAsFactors = FALSE)
    jobs$id <- with(jobs, sprintf("%s-%03d-%s-%s", scenario, dataset, mask, adapter))
    sb_save(jobs, file.path(output, "jobs.rds"))
    # This exit handler also runs on user interrupts; committed results survive.
    on.exit({
        tryCatch(sb_finish_report(root, output), error = function(e) {
            writeLines(conditionMessage(e), file.path(output, "report-error.log"))
            cli::cli_warn("Results were saved, but report rendering failed: {conditionMessage(e)}")
        })
    }, add = TRUE)
    checks_path <- file.path(output, "self-checks.rds")
    if (!file.exists(checks_path)) {
        sb_save(sb_self_checks(sb_seed(config$seed, "self-checks")), checks_path)
    }
    design <- sb_design(config$cells, config$drugs, config$doses, config$merck)
    for (scenario in config$scenarios) for (replicate in seq_len(config$datasets)) {
        dataset_id <- sprintf("%s-%03d", scenario, replicate)
        dataset_path <- file.path(output, "datasets", paste0(dataset_id, ".rds"))
        generation_seed <- sb_seed(config$seed, "generation", scenario, replicate)
        mask_seed <- sb_seed(config$seed, "mask", replicate)
        if (!file.exists(dataset_path)) {
            dataset <- tryCatch(sb_generate(design, scenario, config$rank, generation_seed, mask_seed),
                                error = function(e) list(error = conditionMessage(e), seed = generation_seed, mask_seed = mask_seed))
            sb_save(dataset, dataset_path)
        } else dataset <- readRDS(dataset_path)
        summary_path <- file.path(output, "designs", paste0(dataset_id, ".rds"))
        sb_save(sb_dataset_summary(dataset), summary_path)
        rows <- which(jobs$scenario == scenario & jobs$dataset == replicate)
        for (i in rows) {
            job <- jobs[i, ]
            path <- file.path(output, "fits", paste0(job$id, ".rds"))
            if (file.exists(path)) next
            cli::cli_inform("[{i}/{nrow(jobs)}] {job$id}")
            settings <- list(chains = config$chains, iter_warmup = config$iter_warmup,
                             iter_sampling = config$iter_sampling, thin = 1L,
                             parallel_chains = 1L, refresh = 0L,
                             seed = sb_seed(config$seed, "fit", scenario, replicate, job$mask))
            result <- sb_fit_job(dataset, job$mask, adapters[[job$adapter]], settings,
                                 sb_seed(config$seed, "evaluation", scenario, replicate, job$mask),
                                 sb_seed(config$seed, "ranks", scenario, replicate, job$mask))
            result$job <- job
            sb_save(result, path)
            if (result$status == "failed") cli::cli_warn("{job$id}: {result$error}")
        }
    }
    invisible(output)
}

sb_finish_report <- function(root, output) {
    bundle <- sb_collect(output)
    sb_save(bundle, file.path(output, "report-data.rds"))
    for (name in c("metrics", "diagnostics", "ranks", "coverage", "status")) {
        sb_csv(bundle[[name]], file.path(output, paste0(name, ".csv")))
    }
    sb_render_report(root, output)
    invisible(bundle)
}

sb_render_report <- function(root, output) {
    sb_preflight()
    # Copy the template and renderer so report rendering never loads package code.
    for (name in c("report.qmd", "report.R")) {
        if (!file.copy(file.path(root, "experiments/samplers", name), file.path(output, name), overwrite = TRUE)) {
            cli::cli_abort("Could not stage report renderer.")
        }
    }
    # devtools tests may export a project-relative startup profile. The report
    # child only needs the already resolved library paths, never project startup.
    variables <- c("R_PROFILE_USER", "R_LIBS")
    previous <- Sys.getenv(variables, unset = NA_character_)
    profile <- tempfile("report-profile-")
    writeLines(character(), profile)
    on.exit({
        unlink(profile)
        for (name in variables) {
            if (is.na(previous[[name]])) Sys.unsetenv(name)
            else do.call(Sys.setenv, stats::setNames(list(previous[[name]]), name))
        }
    }, add = TRUE)
    Sys.setenv(R_PROFILE_USER = profile, R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
    quarto::quarto_render(input = file.path(output, "report.qmd"), output_file = "report.html", quiet = TRUE)
    unlink(file.path(output, "report-error.log"))
    invisible(file.path(output, "report.html"))
}
