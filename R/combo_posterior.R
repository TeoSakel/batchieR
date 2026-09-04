# Posterior draws and ecosystem integration for combo fits.

utils::globalVariables(".data")

# Prediction interface and workers ------------------------

#' Posterior predictions for combination-response fits
#'
#' These functions return posterior linear predictors, expected responses, or
#' replicated observations. For the supported Gaussian identity-link model,
#' `posterior_linpred()` and `posterior_epred()` are identical.
#'
#' @param object A `combo_fit` object.
#' @param newdata Optional combination-screen rows. When `NULL`, predictions
#'   are returned for every input row, including rows with missing responses.
#' @param type Either `"draws"` for the draws-by-rows matrix or `"mean"` for
#'   posterior means by row.
#' @param observation Whether [predict()] should include Gaussian observation noise.
#' @param ... Arguments passed between methods.
#' @return A numeric draws-by-rows matrix, except `predict(..., type = "mean")`,
#'   which returns one posterior mean per row.
#' @name posterior_predictions
NULL

#' @rdname posterior_predictions
#' @export
posterior_linpred <- function(object, ...) {
    UseMethod("posterior_linpred")
}

#' @rdname posterior_predictions
#' @export
posterior_linpred.combo_fit <- function(object, newdata = NULL, ...) {
    validate_combo_fit(object)
    indices <- combo_prediction_indices(object, newdata)
    predictions <- lapply(
        object$draws,
        predict_combo_draw,
        indices = indices,
        rank = object$model$rank
    )
    predictions <- do.call(rbind, predictions)
    colnames(predictions) <- as.character(indices$row_id)
    predictions
}

#' @rdname posterior_predictions
#' @export
posterior_epred <- function(object, ...) {
    UseMethod("posterior_epred")
}

#' @rdname posterior_predictions
#' @export
posterior_epred.combo_fit <- function(object, newdata = NULL, ...) {
    link <- object$model$family$link
    if (!identical(link, "identity")) {
        cli::cli_abort("No inverse-link implementation is available for {link}")
    }
    posterior_linpred(object, newdata = newdata, ...)
}

#' @rdname posterior_predictions
#' @export
posterior_predict <- function(object, ...) {
    UseMethod("posterior_predict")
}

#' @rdname posterior_predictions
#' @export
posterior_predict.combo_fit <- function(object, newdata = NULL, ...) {
    expected <- posterior_epred(object, newdata = newdata, ...)
    result <- expected
    draw_sd <- 1 / sqrt(vapply(object$draws, function(draw) draw[["precision"]], numeric(1)))
    for (response in seq_len(ncol(expected))) {
        result[, response] <- stats::rnorm(
            nrow(expected),
            expected[, response],
            draw_sd
        )
    }
    result
}

#' @rdname posterior_predictions
#' @export
predict.combo_fit <- function(
    object,
    newdata = NULL,
    type = c("draws", "mean"),
    observation = FALSE,
    ...
) {
    type <- match.arg(type)
    if (!is.logical(observation) || length(observation) != 1L || is.na(observation)) {
        cli::cli_abort("observation must be TRUE or FALSE")
    }
    predictions <- if (isTRUE(observation)) {
        posterior_predict(object, newdata = newdata, ...)
    } else {
        posterior_epred(object, newdata = newdata, ...)
    }
    if (type == "draws") {
        return(predictions)
    }
    colMeans(predictions)
}

predict_combo_draw <- function(draw, indices, rank) {
    n <- length(indices$cell)
    mu <- rep(draw$alpha, n)
    W0 <- draw$components$cell_offset
    if (!is.null(W0)) {
        mu <- mu + W0$value[indices$cell, 1L]
    }
    V0 <- draw$components$treatment_offset
    if (!is.null(V0)) {
        mu <- mu +
            get_snapshot_rows(V0, indices$treatment_1, 1L)[, 1L] +
            get_snapshot_rows(V0, indices$treatment_2, 1L)[, 1L]
    }
    W <- draw$components$cell_factors
    if (!is.null(W)) {
        w <- W$value[indices$cell, , drop = FALSE]
        V1 <- draw$components$treatment_main_factors
        if (!is.null(V1)) {
            main <- get_snapshot_rows(V1, indices$treatment_1, rank) +
                get_snapshot_rows(V1, indices$treatment_2, rank)
            mu <- mu + rowSums(w * main)
        }
        V2 <- draw$components$treatment_interaction_factors
        if (!is.null(V2)) {
            interaction <- get_snapshot_rows(
                V2,
                indices$treatment_1,
                rank
            ) * get_snapshot_rows(
                V2,
                indices$treatment_2,
                rank
            )
            mu <- mu + rowSums(w * interaction)
        }
    }
    mu
}

combo_prediction_indices <- function(fit, newdata) {
    if (is.null(newdata)) {
        return(list(
            cell = fit$compiled$all_cell,
            treatment_1 = fit$compiled$all_treatment_1,
            treatment_2 = fit$compiled$all_treatment_2,
            row_id = seq_len(nrow(fit$compiled$data))
        ))
    }
    parsed <- validate_experiments(newdata, require_response = FALSE)
    cell <- match(parsed$cell, fit$compiled$cells)
    treatment_1 <- match(parsed$key_1, fit$compiled$treatments$key)
    treatment_2 <- match(parsed$key_2, fit$compiled$treatments$key)
    treatment_1[is.na(parsed$key_1)] <- 0L
    treatment_2[is.na(parsed$key_2)] <- 0L
    unseen_cells <- unique(parsed$cell[is.na(cell)])
    unseen_treatments <- unique(c(
        parsed$key_1[is.na(treatment_1) & !is.na(parsed$key_1)],
        parsed$key_2[is.na(treatment_2) & !is.na(parsed$key_2)]
    ))
    if (length(unseen_cells) || length(unseen_treatments)) {
        cli::cli_abort(c(
            "newdata contains entities absent from the fitted mappings:",
            if (length(unseen_cells)) {
                "*" = "cells: {.and {unseen_cells}}"
            },
            if (length(unseen_treatments)) {
                "*" = "treatments: {.and {unseen_treatments}}"
            }
        ))
    }
    list(
        cell = cell,
        treatment_1 = treatment_1,
        treatment_2 = treatment_2,
        row_id = seq_len(nrow(newdata))
    )
}


# Posterior draw interface and workers --------------------

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws
as_draws.combo_fit <- function(x, ...) {
    posterior::as_draws(posterior_draws(x, ...))
}

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_array
as_draws_array.combo_fit <- function(x, ...) {
    posterior::as_draws_array(posterior_draws(x, ...))
}

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_df
as_draws_df.combo_fit <- function(x, ...) {
    posterior::as_draws_df(posterior_draws(x, ...))
}

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.combo_fit <- function(x, ...) {
    posterior::as_draws_matrix(posterior_draws(x, ...))
}

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_list
as_draws_list.combo_fit <- function(x, ...) {
    posterior::as_draws_list(posterior_draws(x, ...))
}

#' @rdname posterior_draws
#' @exportS3Method posterior::as_draws_rvars
as_draws_rvars.combo_fit <- function(x, ...) {
    posterior::as_draws_rvars(posterior_draws(x, ...))
}

#' Extract posterior draws
#'
#' Converts retained Gibbs snapshots to a numeric iterations-by-chains-by-variables
#' array. Public variables include the intercept, observation scale and precision,
#' enabled component values, mean coefficients, and shrinkage precisions. Use the
#' `posterior::as_draws*()` methods when a `posterior` draws format is required.
#'
#' @param fit A `combo_fit` object.
#' @param components Optional character vector of model components to retain.
#' @param variable Optional character vector of exact posterior variable names.
#'   It is mutually exclusive with `components`.
#' @param include `"public"` for the stable fitted-model interface or `"all"`
#'   to additionally include raw nodes, fitted means, and sampler diagnostics.
#' @param x A `combo_fit` object passed to a [posterior] conversion generic.
#' @param ... Selection arguments forwarded to `posterior_draws()`.
#' @return `posterior_draws()` returns a numeric three-dimensional array. The
#'   `as_draws_*()` methods return the corresponding `posterior` draws format.
#' @export
posterior_draws <- function(
    fit,
    components = NULL,
    variable = NULL,
    include = c("public", "all")
) {
    validate_combo_fit(fit)
    include <- match.arg(include)
    if (!is.null(components) && !is.null(variable)) {
        cli::cli_abort("components and variable are mutually exclusive")
    }
    draws <- combo_draws_array(fit, include)
    map <- parameter_map(fit, include)
    selected <- map$variable
    if (!is.null(components)) {
        valid <- names(fit$model$components)
        unknown <- setdiff(components, valid)
        if (length(unknown)) {
            cli::cli_abort("Unknown components: {.and {unknown}}")
        }
        selected <- map$variable[map$component %in% components]
        if (!length(selected)) {
            cli::cli_abort("Selected components have no posterior variables")
        }
    } else if (!is.null(variable)) {
        unknown <- setdiff(variable, map$variable)
        if (length(unknown)) {
            cli::cli_abort("Unknown variables: {.and {unknown}}")
        }
        selected <- variable
    }
    draws[, , selected, drop = FALSE]
}

combo_draws_array <- function(fit, include) {
    chain_sizes <- tabulate(
        fit$chain_id,
        nbins = max(fit$chain_id)
    )
    if (!length(chain_sizes) || any(chain_sizes != chain_sizes[1L])) {
        cli::cli_abort("combo_fit contains unbalanced posterior chains")
    }
    flattened <- lapply(
        seq_along(fit$draws),
        function(index) {
            combo_flatten_snapshot(
                fit,
                fit$draws[[index]],
                index,
                include
            )
        }
    )
    variable_names <- names(flattened[[1L]])
    consistent <- vapply(
        flattened,
        function(draw) identical(names(draw), variable_names),
        logical(1)
    )
    if (!all(consistent)) {
        cli::cli_abort("Posterior snapshots do not share a stable variable schema")
    }
    values <- array(
        NA_real_,
        dim = c(
            chain_sizes[1L],
            length(chain_sizes),
            length(variable_names)
        ),
        dimnames = list(
            iteration = as.character(seq_len(chain_sizes[1L])),
            chain = as.character(seq_along(chain_sizes)),
            variable = variable_names
        )
    )
    for (index in seq_along(flattened)) {
        values[fit$draw_id[index], fit$chain_id[index], ] <-
            flattened[[index]]
    }
    values
}

#' Map posterior variables to model entities
#'
#' Returns stable variable names and their component, parameter, entity, latent
#' dimension, and metadata-feature meanings. Rows are aligned with the variable
#' dimension returned by [posterior_draws()].
#'
#' The canonical names can be passed directly to [posterior::rename_variables()]
#' or used to construct a programmatic renaming call. Renaming changes only the
#' extracted draws; the fit and its parameter map remain unchanged.
#'
#' @inheritParams posterior_draws
#' @return A data frame with one row per posterior variable.
#' @examples
#' \dontrun{
#' draws <- posterior::as_draws(fit)
#' map <- parameter_map(fit)
#'
#' # Inspect canonical names before choosing a presentation-specific name.
#' map[map$component == "cell_offset",
#'     c("variable", "entity_label", "dimension")]
#'
#' # Rename one selected variable programmatically in the extracted draws.
#' target <- map[map$component == "cell_offset", ][1L, ]
#' rename_args <- stats::setNames(
#'     as.list(target$variable),
#'     paste0(target$entity_label, "_cell_offset")
#' )
#' renamed_draws <- do.call(
#'     posterior::rename_variables,
#'     c(list(.x = draws), rename_args)
#' )
#' }
#' @export
parameter_map <- function(fit, include = c("public", "all")) {
    validate_combo_fit(fit)
    include <- match.arg(include)
    snapshot <- fit$draws[[1L]]
    result <- list(
        combo_basic_map_row("alpha", "global", "value"),
        combo_basic_map_row("sigma", "observation", "scale"),
        combo_basic_map_row("observation_precision", "observation", "precision")
    )
    for (component_name in names(snapshot$components)) {
        component <- snapshot$components[[component_name]]
        if (is.null(component)) {
            next
        }
        result[[length(result) + 1L]] <- combo_matrix_map(
            fit,
            component_name,
            "value",
            component_name,
            component$value
        )
        if (length(component$beta)) {
            beta_map <- combo_basic_map_row(
                combo_vector_variable_names(
                    paste0(component_name, "_mean_coefficient"),
                    component$beta
                ),
                component_name,
                "mean_coefficient"
            )
            beta_map$feature <- names(component$beta)
            result[[length(result) + 1L]] <- beta_map
        }
        if (!is.null(component$beta_precision)) {
            result[[length(result) + 1L]] <- combo_basic_map_row(
                paste0(component_name, "_beta_precision"),
                component_name,
                "beta_precision"
            )
        }
        result[[length(result) + 1L]] <- combo_hyperparameter_map(
            component_name,
            "global_precision",
            paste0(component_name, "_global_precision"),
            component$global_precision
        )
        if (!is.null(component$local_precision)) {
            result[[length(result) + 1L]] <- combo_matrix_map(
                fit,
                component_name,
                "local_precision",
                paste0(component_name, "_local_precision"),
                component$local_precision
            )
        }
        if (include == "all") {
            raw_map <- combo_basic_map_row(
                combo_matrix_variable_names(
                    paste0(component_name, "_raw"),
                    component$raw
                ),
                component_name,
                "raw"
            )
            n_dims <- ncol(component$raw)
            n_nodes <- nrow(component$raw)
            node_names <- rownames(component$raw)
            raw_map$entity_type <- "structure_node"
            raw_map$entity_index <- rep(seq_len(n_nodes), times = n_dims)
            raw_map$entity_key <- rep(node_names, times = n_dims)
            raw_map$entity_label <- raw_map$entity_key
            raw_map$dimension <- rep(seq_len(n_dims), each = n_nodes)
            result[[length(result) + 1L]] <- raw_map
        }
    }
    if (include == "all") {
        fitted_map <- combo_basic_map_row(
            combo_vector_variable_names("fitted_mean", snapshot$Mu),
            "observation",
            "fitted_mean"
        )
        fitted_map$entity_type <- "observation"
        fitted_map$entity_index <- seq_along(snapshot$Mu)
        fitted_map$entity_key <- as.character(fit$compiled$observed_rows)
        fitted_map$entity_label <- fitted_map$entity_key
        result <- c(
            result,
            list(
                fitted_map,
                combo_basic_map_row(
                    c(
                        "sampler_rmse",
                        "sampler_step",
                        "sampler_iteration"
                    ),
                    "sampler",
                    c("rmse", "step", "iteration")
                )
            )
        )
    }
    result <- do.call(rbind, result)
    rownames(result) <- NULL
    expected <- names(combo_flatten_snapshot(fit, snapshot, 1L, include))
    if (!identical(result$variable, expected)) {
        cli::cli_abort("Internal posterior variable map is inconsistent")
    }
    result
}

combo_flatten_snapshot <- function(fit, snapshot, draw_index, include) {
    result <- c(
        alpha = snapshot$alpha,
        sigma = 1 / sqrt(snapshot$precision),
        observation_precision = snapshot$precision
    )
    for (component_name in names(snapshot$components)) {
        component <- snapshot$components[[component_name]]
        if (is.null(component)) {
            next
        }
        values <- as.numeric(component$value)
        names(values) <- combo_matrix_variable_names(
            component_name,
            component$value
        )
        result <- c(result, values)
        if (length(component$beta)) {
            beta <- as.numeric(component$beta)
            names(beta) <- combo_vector_variable_names(
                paste0(component_name, "_mean_coefficient"),
                component$beta
            )
            result <- c(result, beta)
        }
        if (!is.null(component$beta_precision)) {
            beta_precision <- component$beta_precision
            names(beta_precision) <-
                paste0(component_name, "_beta_precision")
            result <- c(result, beta_precision)
        }
        global <- as.numeric(component$global_precision)
        names(global) <- combo_vector_variable_names(
            paste0(component_name, "_global_precision"),
            component$global_precision
        )
        result <- c(result, global)
        if (!is.null(component$local_precision)) {
            local <- as.numeric(component$local_precision)
            names(local) <- combo_matrix_variable_names(
                paste0(component_name, "_local_precision"),
                component$local_precision
            )
            result <- c(result, local)
        }
        if (include == "all") {
            raw <- as.numeric(component$raw)
            names(raw) <- combo_matrix_variable_names(
                paste0(component_name, "_raw"),
                component$raw
            )
            result <- c(result, raw)
        }
    }
    if (include == "all") {
        fitted_mean <- as.numeric(snapshot$Mu)
        names(fitted_mean) <- combo_vector_variable_names(
            "fitted_mean",
            snapshot$Mu
        )
        diagnostics <- c(
            sampler_rmse = snapshot$last_rmse,
            sampler_step = snapshot$n_steps,
            sampler_iteration = fit$iteration[draw_index]
        )
        result <- c(result, fitted_mean, diagnostics)
    }
    result
}

combo_component_entities <- function(fit, component_name) {
    compiled <- fit$compiled$components[[component_name]]
    if (compiled$side == "cell") {
        keys <- fit$compiled$cells
        return(list(type = "cell", key = keys, label = keys))
    }
    treatments <- fit$compiled$treatments
    list(
        type = "treatment",
        key = treatments$key,
        label = paste0(treatments$drug, " @ ", treatments$dose)
    )
}

combo_matrix_map <- function(fit, component_name, parameter, prefix, x) {
    entity <- combo_component_entities(fit, component_name)
    entity_index <- rep(seq_len(nrow(x)), times = ncol(x))
    dimension <- rep(seq_len(ncol(x)), each = nrow(x))
    result <- combo_map_template(length(entity_index))
    result$variable <- combo_matrix_variable_names(prefix, x)
    result$component <- component_name
    result$parameter <- parameter
    result$entity_type <- entity$type
    result$entity_index <- entity_index
    result$entity_key <- entity$key[entity_index]
    result$entity_label <- entity$label[entity_index]
    result$dimension <- dimension
    result
}

combo_hyperparameter_map <- function(component_name, parameter, prefix, x) {
    result <- combo_basic_map_row(
        combo_vector_variable_names(prefix, x),
        component_name,
        parameter
    )
    result$dimension <- seq_along(x)
    result
}

# Model checking interface and workers --------------------

#' Pointwise log-likelihood for combination-response fits
#'
#' Computes the Gaussian log-likelihood for each retained draw and each row
#' with an observed response. Missing-response rows are omitted.
#'
#' @param object A `combo_fit` object.
#' @param newdata Optional combination-screen data containing responses. When
#'   `NULL`, the observed rows from the fitted input are used.
#' @param ... Reserved for future methods.
#' @return A numeric draws-by-observed-rows matrix. Column names and the
#'   `row_ids` attribute identify the corresponding rows in the input data.
#' @name log_lik
NULL

#' @rdname log_lik
#' @export
log_lik <- function(object, ...) {
    UseMethod("log_lik")
}

#' @rdname log_lik
#' @export
log_lik.combo_fit <- function(object, newdata = NULL, ...) {
    validate_combo_fit(object)
    observed <- combo_observed_data(object, newdata)
    expected <- posterior_epred(object, newdata = observed$data)
    draw_sd <- 1 / sqrt(vapply(object$draws, function(draw) draw[["precision"]], numeric(1)))
    result <- expected
    for (response in seq_len(ncol(expected))) {
        result[, response] <- stats::dnorm(
            observed$response[response],
            mean = expected[, response],
            sd = draw_sd,
            log = TRUE
        )
    }
    colnames(result) <- as.character(observed$row_id)
    attr(result, "row_ids") <- observed$row_id
    result
}

#' Posterior predictive checks for combination-response fits
#'
#' Generates Gaussian replicated responses and produces predictive-check plots
#' using `bayesplot` and `ggplot2`.
#' Only observed rows are included. Pooled checks can be conditioned with
#' grouping, predictive errors, test statistics, and PSIS-LOO calibration.
#'
#' @param object A `combo_fit` object.
#' @param type Check type: `"dens_overlay"`, `"ecdf_overlay"`, `"intervals"`,
#'   `"stat"`, `"stat_2d"`, `"error"`, or `"loo_pit"`.
#' @param ndraws `NULL`, or a positive number of posterior draws, capped at
#'   the retained draws. By default overlays use at most 50 draws and other
#'   non-LOO checks use all draws. LOO-PIT always uses all draws.
#' @param newdata Optional data containing responses. LOO-PIT requires `NULL`
#'   because it is only defined for observations used to fit the model.
#' @param stat Statistic for statistic or error checks. `NULL` uses `"mean"`,
#'   except `type = "stat_2d"` uses `c("mean", "sd")`.
#' @param group `NULL`, one source-data column name, or a vector aligned with
#'   all source rows. Supported by grouped overlays, intervals, statistics,
#'   and errors.
#' @param x `NULL`, one numeric source-data column name, or a finite numeric
#'   vector aligned with all source rows. Used by intervals and
#'   errors. Errors default to posterior mean fitted responses.
#' @param ... Arguments passed to the selected plotting function. For `type = "error"`, they are passed to [ggplot2::geom_point()].
#' @return A `ggplot` object.
#' @details Density overlays assume smooth continuous responses. Prefer an
#' ECDF for rounded, censored, bounded, or point-mass responses. LOO-PIT uses
#' simultaneous confidence bands; inspect [loo::pareto_k_table()] when warned.
#' @exportS3Method bayesplot::pp_check
pp_check.combo_fit <- function(
    object,
    type = c(
        "dens_overlay", "ecdf_overlay", "intervals", "stat",
        "stat_2d", "error", "loo_pit"
    ),
    ndraws = NULL,
    newdata = NULL,
    stat = NULL,
    group = NULL,
    x = NULL,
    ...
) {
    if (!requireNamespace("bayesplot", quietly = TRUE)) {
        cli::cli_abort("{.pkg bayesplot} is required for {.fn pp_check}.")
    }
    type <- match.arg(type)
    combo_validate_ppc_options(type, ndraws, newdata, group, x)
    observed <- combo_observed_data(object, newdata)
    group <- combo_ppc_variable(group, "group", observed)
    x <- combo_ppc_variable(x, "x", observed, numeric = TRUE)

    if (identical(type, "loo_pit")) {
        return(combo_ppc_loo_pit(object, observed, ...))
    }

    yrep <- posterior_predict(object, newdata = observed$data)
    yrep <- combo_ppc_draws(yrep, type, ndraws)
    y <- observed$response
    if (identical(type, "dens_overlay")) combo_warn_ppc_density(y)
    if (is.null(stat)) {
        stat <- if (identical(type, "stat_2d")) c("mean", "sd") else "mean"
    }
    if (identical(type, "stat") && length(stat) != 1L) {
        cli::cli_abort("{.arg stat} must contain one statistic for {.val stat} checks.")
    }
    if (identical(type, "stat_2d") && length(stat) != 2L) {
        cli::cli_abort("{.arg stat} must contain two statistics for {.val stat_2d} checks.")
    }

    if (identical(type, "error")) {
        return(combo_ppc_error_plot(object, observed, yrep, x, group, stat, ...))
    }

    if (!is.null(group)) {
        return(switch(
            type,
            dens_overlay = bayesplot::ppc_dens_overlay_grouped(y, yrep, group, ...),
            ecdf_overlay = bayesplot::ppc_ecdf_overlay_grouped(y, yrep, group, ...),
            intervals = bayesplot::ppc_intervals_grouped(y, yrep, x, group, ...),
            stat = bayesplot::ppc_stat_grouped(y, yrep, group, stat = stat, ...)
        ))
    }

    switch(
        type,
        dens_overlay = bayesplot::ppc_dens_overlay(y, yrep, ...),
        ecdf_overlay = bayesplot::ppc_ecdf_overlay(y, yrep, ...),
        intervals = bayesplot::ppc_intervals(y, yrep, x = x, ...),
        stat = bayesplot::ppc_stat(y, yrep, stat = stat, ...),
        stat_2d = bayesplot::ppc_stat_2d(y, yrep, stat = stat, ...)
    )
}

combo_validate_ppc_options <- function(type, ndraws, newdata, group, x) {
    if (!is.null(ndraws)) param_positive_integer(ndraws, "ndraws")
    if (identical(type, "loo_pit")) {
        if (!is.null(newdata)) {
            cli::cli_abort("LOO-PIT is only available for fitted observations; {.arg newdata} must be {.code NULL}.")
        }
        if (!is.null(ndraws)) {
            cli::cli_abort("LOO-PIT requires every retained draw; {.arg ndraws} must be {.code NULL}.")
        }
    }
    grouped_types <- c("dens_overlay", "ecdf_overlay", "intervals", "stat", "error")
    if (!is.null(group) && !type %in% grouped_types) {
        cli::cli_abort("{.arg group} is not supported for a {.val {type}} check.")
    }
    if (!is.null(x) && !type %in% c("intervals", "error")) {
        cli::cli_abort("{.arg x} is only supported for interval and error checks.")
    }
    invisible(NULL)
}

combo_ppc_variable <- function(value, name, observed, numeric = FALSE) {
    if (is.null(value)) return(NULL)
    if (is.character(value) && length(value) == 1L) {
        if (!value %in% names(observed$data)) {
            cli::cli_abort("Column {.field {value}} supplied to {.arg {name}} is not in the checking data.")
        }
        value <- observed$data[[value]]
    } else {
        if (length(value) != observed$n_rows) {
            cli::cli_abort("{.arg {name}} must name one column or have length {observed$n_rows}.")
        }
        value <- value[observed$row_id]
    }
    if (anyNA(value)) {
        cli::cli_abort("{.arg {name}} must not contain missing values on observed rows.")
    }
    if (numeric && (!is.numeric(value) || any(!is.finite(value)))) {
        cli::cli_abort("{.arg {name}} must resolve to finite numeric values.")
    }
    value
}

combo_ppc_error_plot <- function(object, observed, yrep, x, group, stat, ...) {
    stat_fn <- match.fun(stat)
    error <- vapply(seq_along(observed$response), function(index) {
        value <- stat_fn(observed$response[index] - yrep[, index])
        if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
            cli::cli_abort("{.arg stat} must return one finite numeric predictive-error summary.")
        }
        value
    }, numeric(1))
    default_x <- is.null(x)
    if (default_x) {
        x <- colMeans(posterior_epred(object, newdata = observed$data))
    }
    plot_data <- data.frame(x = x, error = error)
    if (!is.null(group)) plot_data$group <- group
    plot <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(x = .data[["x"]], y = .data[["error"]])
    ) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
        ggplot2::geom_point(...) +
        ggplot2::labs(
            x = if (default_x) "Posterior mean fitted response" else "Checking variable",
            y = "Posterior predictive error"
        )
    if (!is.null(group)) {
        plot <- plot + ggplot2::facet_wrap(stats::as.formula("~ group"))
    }
    plot
}

combo_ppc_draws <- function(yrep, type, ndraws) {
    target <- if (is.null(ndraws)) {
        if (type %in% c("dens_overlay", "ecdf_overlay")) 50L else nrow(yrep)
    } else {
        param_positive_integer(ndraws, "ndraws")
    }
    target <- min(target, nrow(yrep))
    if (target == nrow(yrep)) return(yrep)
    yrep[sample.int(nrow(yrep), target, replace = FALSE), , drop = FALSE]
}

combo_warn_ppc_density <- function(y) {
    counts <- table(y)
    if (length(counts) && max(counts) > 1L && max(counts) / length(y) > 0.02) {
        cli::cli_warn(c(
            "The response contains repeated values with more than 2% empirical mass.",
            "i" = "A density overlay can hide rounding, censoring, or point masses; also inspect an ECDF overlay."
        ))
    }
    invisible(NULL)
}

combo_loo_pit_values <- function(y, yrep, weights) {
    if (!identical(dim(yrep), dim(weights)) || ncol(yrep) != length(y)) {
        cli::cli_abort("Internal LOO-PIT inputs are not aligned")
    }
    indicators <- sweep(yrep, 2L, y, FUN = "<=")
    colSums(weights * indicators)
}

combo_ppc_loo_pit <- function(object, observed, ...) {
    if (!requireNamespace("loo", quietly = TRUE)) {
        cli::cli_abort("{.pkg loo} is required for a LOO-PIT check.")
    }
    pointwise <- log_lik(object)
    r_eff <- loo::relative_eff(exp(pointwise), chain_id = object$chain_id)
    loo_result <- loo::loo(pointwise, r_eff = r_eff, save_psis = TRUE)
    k_table <- loo::pareto_k_table(loo_result)
    threshold <- attr(k_table, "k_threshold")
    bad <- sum(loo::pareto_k_values(loo_result) > threshold, na.rm = TRUE)
    if (bad) {
        cli::cli_warn(c(
            "{bad} observation{?s} exceed{?s/} the Pareto-k reliability threshold ({format(threshold, digits = 2)}).",
            "i" = "Inspect {.code loo::pareto_k_table(loo::loo(fit))}; the LOO-PIT approximation may be unreliable."
        ))
    }
    yrep <- posterior_predict(object, newdata = observed$data)
    dots <- list(...)
    if (is.null(dots$plot_diff)) dots$plot_diff <- TRUE
    weights <- stats::weights(loo_result$psis_object, normalize = TRUE, log = FALSE)
    pit <- combo_loo_pit_values(observed$response, yrep, weights)
    do.call(
        bayesplot::ppc_loo_pit_ecdf,
        c(list(pit = pit), dots)
    )
}

#' Leave-one-out cross-validation for combination-response fits
#'
#' Estimates out-of-sample predictive performance using Pareto-smoothed
#' importance-sampling leave-one-out cross-validation (PSIS-LOO). Use the
#' result to assess a model or compare candidate models fitted to the same
#' observations. The suggested `loo` package is required.
#'
#' The printed result includes:
#'
#' - `elpd_loo`, the expected log predictive density; larger values indicate
#'   better expected predictive performance.
#' - `p_loo`, an estimate of the effective number of model parameters.
#' - `looic`, an information criterion equal to `-2 * elpd_loo`; smaller values
#'   indicate better expected predictive performance.
#' - Pareto `k` diagnostics. Large values identify influential observations for
#'   which the approximation may be unreliable and should be investigated.
#'
#' Compare models with [loo::loo_compare()]. Differences are meaningful only
#' when the models use the same response observations.
#'
#' @param x A `combo_fit` object.
#' @param newdata Optional combination-screen data containing responses. When
#'   `NULL`, the observed rows from the fitted input are evaluated.
#' @param ... Additional options passed to [loo::loo()].
#' @param r_eff Optional relative effective sample sizes used to account for
#'   autocorrelation. It is normally left as `NULL` and computed automatically.
#' @return A `loo` object containing predictive-performance estimates,
#'   uncertainty measures, pointwise contributions, and Pareto `k` diagnostics.
#'   Its `combo_row_ids` attribute maps pointwise results back to input rows.
#' @examples
#' \dontrun{
#' loo_result <- loo::loo(fit)
#' print(loo_result)
#' loo::pareto_k_table(loo_result)
#'
#' comparison <- loo::loo_compare(
#'     model_a = loo::loo(fit_a),
#'     model_b = loo::loo(fit_b)
#' )
#' print(comparison)
#' }
#' @exportS3Method loo::loo
loo.combo_fit <- function(x, newdata = NULL, ..., r_eff = NULL) {
    if (!requireNamespace("loo", quietly = TRUE)) {
        cli::cli_abort("loo is required for loo()")
    }
    pointwise <- log_lik(x, newdata = newdata)
    if (is.null(r_eff)) {
        r_eff <- loo::relative_eff(exp(pointwise), chain_id = x$chain_id)
    }
    result <- loo::loo(pointwise, r_eff = r_eff, ...)
    attr(result, "combo_row_ids") <- attr(pointwise, "row_ids")
    result
}

combo_observed_data <- function(object, newdata = NULL) {
    if (is.null(newdata)) {
        return(with(
            object$compiled,
            list(
                data = data[observed_rows, , drop = FALSE],
                response = response,
                row_id = observed_rows,
                n_rows = nrow(data)
            )
        ))
    }
    parsed <- validate_experiments(newdata)
    observed <- which(!is.na(parsed$response))
    if (!length(observed)) {
        cli::cli_abort("newdata contains no observed responses")
    }
    list(
        data = parsed$data[observed, , drop = FALSE],
        response = parsed$response[observed],
        row_id = observed,
        n_rows = nrow(parsed$data)
    )
}

# Shared utilities ---------------------------------------

validate_combo_fit <- function(fit) {
    if (!inherits(fit, "combo_fit")) {
        cli::cli_abort("fit must be a combo_fit")
    }
    invisible(fit)
}

get_snapshot_rows <- function(snapshot, index, n_dimensions) {
    result <- matrix(0, nrow = length(index), ncol = n_dimensions)
    selected <- index > 0L
    if (!is.null(snapshot) && any(selected)) {
        result[selected, ] <- snapshot$value[index[selected], , drop = FALSE]
    }
    result
}

combo_matrix_variable_names <- function(prefix, x) {
    rows <- rep(seq_len(nrow(x)), times = ncol(x))
    dimensions <- rep(seq_len(ncol(x)), each = nrow(x))
    if (ncol(x) == 1L) {
        paste0(prefix, "[", rows, "]")
    } else {
        paste0(prefix, "[", rows, ",", dimensions, "]")
    }
}

combo_vector_variable_names <- function(prefix, x) {
    paste0(prefix, "[", seq_along(x), "]")
}

combo_map_template <- function(n) {
    data.frame(
        variable = rep(NA_character_, n),
        component = rep(NA_character_, n),
        parameter = rep(NA_character_, n),
        entity_type = rep(NA_character_, n),
        entity_index = rep(NA_integer_, n),
        entity_key = rep(NA_character_, n),
        entity_label = rep(NA_character_, n),
        dimension = rep(NA_integer_, n),
        feature = rep(NA_character_, n),
        stringsAsFactors = FALSE
    )
}

combo_basic_map_row <- function(variable, component, parameter) {
    result <- combo_map_template(length(variable))
    result$variable <- variable
    result$component <- component
    result$parameter <- parameter
    result
}
