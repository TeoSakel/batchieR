# Experiment configuration; not part of the installed package API.
sb_profiles <- function() {
    list(
        smoke = list(cells = 2L, drugs = 3L, doses = 2L, rank = 2L,
                     datasets = 1L, chains = 2L, iter_warmup = 2L,
                     iter_sampling = 8L),
        local = list(cells = 4L, drugs = 8L, doses = 3L, rank = 4L,
                     datasets = 2L, chains = 4L, iter_warmup = 250L,
                     iter_sampling = 500L),
        calibration = list(cells = 4L, drugs = 8L, doses = 3L, rank = 4L,
                           datasets = 100L, chains = 4L, iter_warmup = 1000L,
                           iter_sampling = 4000L)
    )
}

sb_scenarios <- function() c("gamma", "multiplicative_gamma", "horseshoe", "realism")

sb_model <- function(scenario, rank) {
    if (!scenario %in% sb_scenarios()) cli::cli_abort("Unknown scenario: {scenario}")
    args <- list(rank = rank, mean = batchieR::fixed_mean(stats::qlogis(0.8)),
                 family = batchieR::gaussian_response(
                     precision = batchieR::gamma_precision(shape = 25, rate = 1)))
    if (scenario %in% c("gamma", "multiplicative_gamma")) {
        component <- batchieR::combo_gaussian_component(
            shrinkage = batchieR::gamma_precision(shape = 3, rate = 0.3))
        for (name in c("cell_offset", "cell_factors", "treatment_offset",
                       "treatment_main_factors", "treatment_interaction_factors")) {
            args[[name]] <- component
        }
        if (scenario == "multiplicative_gamma") {
            args$cell_factors <- batchieR::combo_gaussian_component(
                shrinkage = batchieR::multiplicative_gamma())
        }
    }
    do.call(batchieR::combo_model, args)
}

sb_adapters <- function() {
    list(gibbs = list(
        id = "gibbs", name = "Current batchieR Gibbs sampler",
        description = paste(
            "Fixed block order: mean, cell offsets, treatment offsets, cell factors,",
            "interaction factors, main treatment factors, hyperparameters, observation precision.",
            "Each chain uses the package's zero-initialized component values and default",
            "precision initialization, with an independent RNG stream.",
            "This is the unmodified fit_combo() baseline."),
        fit = function(model, data, settings) {
            do.call(batchieR::fit_combo, c(list(model = model, data = data), settings))
        }
    ))
}

sb_validate_adapters <- function(adapters) {
    if (!length(adapters) || is.null(names(adapters)) || anyDuplicated(names(adapters))) {
        cli::cli_abort("Adapters must be a nonempty, uniquely named list.")
    }
    for (id in names(adapters)) {
        a <- adapters[[id]]
        if (!grepl("^[a-z][a-z0-9_]*$", id) || !identical(a$id, id) ||
            !is.function(a$fit) || !is.character(a$name) || length(a$name) != 1L ||
            is.na(a$name) || !nzchar(a$name) || !is.character(a$description) ||
            length(a$description) != 1L || is.na(a$description) || !nzchar(a$description)) {
            cli::cli_abort("Adapter {id} needs an ID, name, description, and fit function.")
        }
    }
    invisible(adapters)
}

# Stable seed derivation does not depend on candidate ordering or failures.
sb_seed <- function(seed, ...) {
    bytes <- utf8ToInt(paste(c(seed, ...), collapse = ":"))
    value <- 0
    for (byte in bytes) value <- (value * 131 + byte) %% 2147483646
    as.integer(value + 1)
}

sb_with_seed <- function(seed, code) {
    existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (existed) old <- get(".Random.seed", envir = .GlobalEnv)
    kind <- RNGkind()
    on.exit({
        do.call(RNGkind, as.list(kind))
        if (existed) assign(".Random.seed", old, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
            rm(".Random.seed", envir = .GlobalEnv)
    })
    RNGkind("Mersenne-Twister", "Inversion", "Rejection")
    set.seed(seed)
    force(code)
}

sb_internal <- function(name) getFromNamespace(name, "batchieR")

sb_save <- function(object, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- tempfile(".commit-", tmpdir = dirname(path))
    on.exit(unlink(temporary))
    saveRDS(object, temporary)
    if (!file.rename(temporary, path)) cli::cli_abort("Cannot commit {path}.")
    invisible(path)
}
