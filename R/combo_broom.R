# Broom-compatible interfaces for combination fits.

#' @importFrom generics tidy
#' @export
generics::tidy

#' @importFrom generics glance
#' @export
generics::glance

#' @importFrom generics augment
#' @export
generics::augment

#' Tidy posterior parameters of a combination-response fit
#'
#' Returns summaries of all public posterior variables by default, retaining
#' canonical parameter names and their order. Selection follows
#' [posterior_draws()]. Chain boundaries are preserved for diagnostics.
#'
#' @param x A `combo_fit` object.
#' @param robust Use posterior median and MAD instead of mean and SD.
#'   MAD uses the defaults of [stats::mad()]. This does not change intervals
#'   or convergence diagnostics.
#' @param conf.int Include equal-tailed credible interval bounds.
#' @param conf.level Credible interval probability, strictly between zero and one.
#' @inheritParams posterior_draws
#' @param ... Unused arguments are ignored with a warning.
#' @return A tibble with `term`, `estimate`, `std.error`, `rhat`, `ess_bulk`,
#'   and `ess_tail`, plus `conf.low` and `conf.high` when requested. Here
#'   `std.error` is posterior SD or MAD, not Monte Carlo standard error.
#'   Metadata columns from [parameter_map()] follow, with its `variable`
#'   identifier represented by `term`. Unavailable diagnostics remain `NA`.
#' @export
tidy.combo_fit <- function(
    x,
    robust = FALSE,
    conf.int = TRUE,
    conf.level = 0.95,
    components = NULL,
    variable = NULL,
    include = c("public", "all"),
    ...
) {
    if (...length()) cli::cli_warn("Unused arguments in {.code ...} are ignored.")
    if (!is_logical(robust)) {
        cli::cli_abort("{.arg robust} must be TRUE or FALSE.")
    }
    if (!is_logical(conf.int)) {
        cli::cli_abort("{.arg conf.int} must be TRUE or FALSE.")
    }
    if (!is_positive_number(conf.level) || conf.level >= 1) {
        cli::cli_abort("{.arg conf.level} must be a finite number strictly between zero and one.")
    }
    include <- match.arg(include)
    draws <- posterior::as_draws_array(
        x, components = components, variable = variable, include = include
    )
    measures <- c(
        if (robust) c("median", "mad") else c("mean", "sd"),
        posterior::default_convergence_measures()
    )
    result <- posterior::summarize_draws(draws, measures)
    names(result)[1:3] <- c("term", "estimate", "std.error")
    result <- tibble::as_tibble(result)
    if (conf.int) {
        bounds <- posterior::summarize_draws(
            draws,
            interval = function(z) combo_broom_interval(z, conf.level)
        )
        result$conf.low <- bounds$lower
        result$conf.high <- bounds$upper
    }
    map <- parameter_map(x, include)
    map <- map[match(result$term, map$variable), setdiff(names(map), "variable"), drop = FALSE]
    tibble::as_tibble(c(as.list(result), as.list(map)))
}

#' Glance at a combination-response fit
#'
#' Reports dimensions, observation scale, sampler RMSE, and convergence
#' extrema across public posterior variables. No LOO or information criteria
#' are calculated. Missing diagnostics are counted separately and excluded
#' from extrema; infinite values are retained. An entirely unavailable
#' diagnostic has an `NA` extremum.
#'
#' @inheritParams tidy.combo_fit
#' @return A one-row tibble with `nobs` (observed responses), `n_cells`,
#'   `n_treatments`, `n_chains`, `n_draws` (total retained draws),
#'   `n_variables` (public variables), `sigma` (posterior mean observation SD),
#'   `sampler_rmse` (mean retained sampler RMSE), `rhat_max`, `ess_bulk_min`,
#'   `ess_tail_min`, and unavailable counts `rhat_unavailable`,
#'   `ess_bulk_unavailable`, and `ess_tail_unavailable`.
#' @export
glance.combo_fit <- function(x, ...) {
    if (...length()) cli::cli_warn("Unused arguments in {.code ...} are ignored.")
    diagnostics <- posterior::summarize_draws(
        posterior::as_draws_array(x),
        posterior::default_convergence_measures()
    )
    tibble::tibble(
        nobs = length(x$compiled$observed_rows),
        n_cells = length(x$compiled$cells),
        n_treatments = nrow(x$compiled$treatments),
        n_chains = x$sampling$chains,
        n_draws = length(x$draws),
        n_variables = nrow(diagnostics),
        sigma = mean(vapply(x$draws, function(draw) 1 / sqrt(draw$precision), numeric(1))),
        sampler_rmse = mean(vapply(x$draws, function(draw) draw$last_rmse, numeric(1))),
        rhat_max = combo_broom_extreme(diagnostics$rhat, max),
        ess_bulk_min = combo_broom_extreme(diagnostics$ess_bulk, min),
        ess_tail_min = combo_broom_extreme(diagnostics$ess_tail, min),
        rhat_unavailable = sum(is.na(diagnostics$rhat)),
        ess_bulk_unavailable = sum(is.na(diagnostics$ess_bulk)),
        ess_tail_unavailable = sum(is.na(diagnostics$ess_tail))
    )
}

#' Augment observations with expected-response summaries
#'
#' Adds summaries from [posterior_epred()] without generating observation
#' noise or changing RNG state. Unlike `broom.mixed`'s augmentation of brms
#' fits using predictive draws, uncertainty here concerns the expected response.
#' All input columns and rows are retained, including missing-response rows.
#'
#' @inheritParams tidy.combo_fit
#' @param data Original fitting data, optionally with additional columns. The
#'   modeling columns and their row order must match the original data. Use
#'   `newdata` for other observations.
#' @param newdata Optional prediction data, with or without a response column.
#'   Cannot be supplied together with explicit `data`. Entities must already
#'   occur in the fitted mappings, as required by [posterior_epred()].
#' @param se.fit Include posterior SD of expected responses.
#' @param interval `"none"` or `"confidence"`. The latter adds equal-tailed
#'   credible intervals for expected responses, excluding observation noise.
#' @return A tibble containing the input columns and `.fitted` (posterior mean
#'   expected response), optionally `.se.fit` (posterior SD) and `.lower` and
#'   `.upper` (credible interval bounds). If `response` exists, `.resid` is
#'   `response - .fitted`, with missing residuals retained. Collisions with
#'   generated column names are errors. Empty input returns a zero-row tibble
#'   with the corresponding columns.
#' @export
augment.combo_fit <- function(
    x,
    data = x$input$data,
    newdata = NULL,
    se.fit = TRUE,
    interval = c("none", "confidence"),
    conf.level = 0.95,
    ...
) {
    if (...length()) cli::cli_warn("Unused arguments in {.code ...} are ignored.")
    validate_combo_fit(x)
    if (!is_logical(se.fit)) {
        cli::cli_abort("{.arg se.fit} must be TRUE or FALSE.")
    }
    if (!is_positive_number(conf.level) || conf.level >= 1) {
        cli::cli_abort("{.arg conf.level} must be a finite number strictly between zero and one.")
    }
    interval <- match.arg(interval)
    if (!missing(data) && !is.null(newdata)) {
        cli::cli_abort("Supply only one of {.arg data} and {.arg newdata}.")
    }
    input <- if (is.null(newdata)) data else newdata
    validation_input <- input
    if (!is.null(newdata) && inherits(input, "data.frame") && !"response" %in% names(input)) {
        validation_input$response <- rep(NA_real_, nrow(input))
    }
    parsed <- validate_experiments(validation_input, require_response = is.null(newdata))
    if (is.null(newdata) && !missing(data)) {
        original <- validate_experiments(x$input$data)
        columns <- c("cell", "drug_1", "dose_1", "drug_2", "dose_2", "response")
        encoder <- x$compiled$mean$sources$observation$encoder
        features <- if (is.null(encoder)) character() else all.vars(encoder$terms)
        same_features <- all(vapply(features, function(name) identical(input[[name]], x$input$data[[name]]), logical(1)))
        if (!identical(parsed[columns], original[columns]) || !same_features) {
            cli::cli_abort("{.arg data} must retain the original modeling columns and row order; use {.arg newdata} for other observations.")
        }
    }
    generated <- c(
        ".fitted",
        if (se.fit) ".se.fit",
        if (interval == "confidence") c(".lower", ".upper"),
        if ("response" %in% names(input)) ".resid"
    )
    collisions <- intersect(names(input), generated)
    if (length(collisions)) {
        cli::cli_abort("Input contains generated column name(s): {.and {collisions}}.")
    }
    result <- tibble::as_tibble(input)
    if (!nrow(result)) {
        for (name in generated) result[[name]] <- numeric()
        return(result)
    }
    expected <- posterior_epred(x, newdata = newdata)
    result[[".fitted"]] <- unname(colMeans(expected))
    if (se.fit) result[[".se.fit"]] <- unname(apply(expected, 2L, stats::sd))
    if (interval == "confidence") {
        bounds <- vapply(
            seq_len(ncol(expected)),
            function(i) combo_broom_interval(expected[, i], conf.level),
            c(lower = 0, upper = 0)
        )
        result[[".lower"]] <- unname(bounds["lower", ])
        result[[".upper"]] <- unname(bounds["upper", ])
    }
    if ("response" %in% names(input)) {
        result[[".resid"]] <- parsed$response - result[[".fitted"]]
    }
    result
}

combo_broom_interval <- function(x, level) {
    stats::setNames(
        stats::quantile(x, c((1 - level) / 2, (1 + level) / 2), names = FALSE),
        c("lower", "upper")
    )
}

combo_broom_extreme <- function(x, fun) {
    if (all(is.na(x))) NA_real_ else fun(x, na.rm = TRUE)
}
