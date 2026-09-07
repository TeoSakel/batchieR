# Prior prediction for combination models.

#' Draw predictions from the prior
#'
#' Draws latent conditional means or replicated responses from the generative prior defined by a [combo_model()].
#'
#' @param model A combination-model specification created by [combo_model()].
#' @param data A data frame containing the `cell`, `drug_1`, `dose_1`, `drug_2`,
#'   and `dose_2` design columns. An optional numeric `response` column is used
#'   only as described for `intercept`.
#' @param cell_data An optional data frame of cell-level covariates.
#' @param compound_data An optional data frame of compound-level covariates.
#' @param draws Number of prior-predictive draws. Must be a positive integer.
#' @param seed `NULL`, or a positive integer used to initialize the R
#'   random-number generator.
#' @param intercept An optional finite intercept. It overrides the response mean
#'   for [empirical_mean()] and is required for an intercept-bearing
#'   [formula_mean()]. It must be `NULL` for [fixed_mean()] or an intercept-free
#'   formula.
#' @param type Whether to draw replicated `response` values including observation
#'   noise, or latent conditional `mean` values.
#' @return A numeric matrix with one row per draw and one column per input design row.
#'   With `type = "response"`, it matches the shape and column naming convention
#'   of [posterior_predict()]; with `type = "mean"`, it contains latent conditional means.
#'
#' @examples
#' design <- data.frame(
#'     cell = c("A", "A"),
#'     drug_1 = c("X", "Y"),
#'     dose_1 = c(1, 1),
#'     drug_2 = c(NA, "X"),
#'     dose_2 = c(NA, 1)
#' )
#' model <- combo_model(
#'     rank = 1L,
#'     mean = fixed_mean(0)
#' )
#' prior_predict(model, design, draws = 5L, seed = 123L)
#'
#' @export
prior_predict <- function(
    model,
    data,
    cell_data = NULL,
    compound_data = NULL,
    draws = 100L,
    seed = NULL,
    intercept = NULL,
    type = c("response", "mean")
) {
    if (!inherits(model, "combo_model")) {
        cli::cli_abort("model must be constructed by combo_model()")
    }
    type <- match.arg(type)
    if (!is_positive_integer(draws)) {
        cli::cli_abort("draws must be integer-valued and at least 1")
    }
    if (!is.null(seed) && !is_positive_integer(seed)) {
        cli::cli_abort("seed must be NULL or one finite integer")
    }
    draws <- as.integer(draws)
    if (!is.null(seed)) {
        seed <- as.integer(seed)
    }

    intercept <- prior_intercept(model, data, intercept)
    design <- data
    if (inherits(design, "data.frame")) {
        design$response <- NULL
    }
    compiled <- compile_combo_design(
        model,
        design,
        cell_data = cell_data,
        compound_data = compound_data
    )
    if (!is.null(seed)) {
        set.seed(seed)
    }
    indices <- list(
        cell = compiled$all_cell,
        treatment_1 = compiled$all_treatment_1,
        treatment_2 = compiled$all_treatment_2,
        mean_design = compiled$mean$design,
        row_id = seq_len(nrow(compiled$data))
    )
    result <- matrix(
        NA_real_,
        nrow = draws,
        ncol = nrow(compiled$data),
        dimnames = list(NULL, as.character(indices$row_id))
    )
    for (draw in seq_len(draws)) {
        snapshot <- prior_snapshot(compiled, intercept)
        expected <- predict_combo_draw(
            snapshot,
            indices = indices,
            rank = model$rank
        )
        result[draw, ] <- if (type == "mean") {
            expected
        } else {
            stats::rnorm(
                length(expected),
                mean = expected,
                sd = 1 / sqrt(snapshot$precision)
            )
        }
    }
    result
}

prior_intercept <- function(model, data, intercept) {
    specification <- model$mean
    if (specification$type == "fixed") {
        if (!is.null(intercept)) {
            cli::cli_abort("intercept must be NULL when the model uses fixed_mean()")
        }
        return(specification$value)
    }
    if (specification$type == "formula") {
        has_intercept <- attr(stats::terms(specification$formula), "intercept") == 1L
        if (!has_intercept) {
            if (!is.null(intercept)) {
                cli::cli_abort("intercept must be NULL for an intercept-free formula_mean()")
            }
            return(0)
        }
        if (is.null(intercept)) {
            cli::cli_abort(
                "formula_mean() with an intercept requires intercept for prior prediction"
            )
        }
        return(prior_validate_intercept(intercept))
    }
    if (!is.null(intercept)) {
        return(prior_validate_intercept(intercept))
    }
    if (!inherits(data, "data.frame") || !"response" %in% names(data)) {
        cli::cli_abort(
            "empirical_mean() requires intercept or at least one observed response"
        )
    }
    response <- data$response
    if (!is.numeric(response)) {
        cli::cli_abort("response must be numeric when used for the empirical intercept")
    }
    observed <- response[!is.na(response)]
    if (!length(observed)) {
        cli::cli_abort(
            "empirical_mean() requires intercept or at least one observed response"
        )
    }
    if (any(!is.finite(observed))) {
        cli::cli_abort("Nonmissing responses must be finite numeric values")
    }
    mean(observed)
}

prior_validate_intercept <- function(intercept) {
    if (!is_number(intercept)) {
        cli::cli_abort("intercept must be NULL or one finite numeric value")
    }
    as.numeric(intercept)
}

prior_snapshot <- function(compiled, intercept) {
    precision <- with(
        compiled$observation_prior,
        if (type == "fixed") precision else stats::rgamma(1L, shape = shape, rate = rate)
    )
    component_names <- names(compiled$components)
    components <- lapply(
        component_names,
        function(name) prior_draw_component(compiled$components[[name]])
    )
    names(components) <- component_names
    beta <- if (length(compiled$mean$beta_mean)) {
        stats::setNames(
            compiled$mean$beta_mean +
                batchieR_rmvnorm(compiled$mean$beta_precision),
            compiled$mean$coefficient_names
        )
    } else {
        numeric()
    }
    list(
        intercept = intercept,
        beta = beta,
        precision = precision,
        components = components
    )
}

prior_draw_component <- function(compiled) {
    if (is.null(compiled)) {
        return(NULL)
    }
    n_entities <- compiled$n_entities
    n_dimensions <- compiled$n_dimensions
    shrinkage <- prior_draw_shrinkage(
        compiled$shrinkage,
        n_entities,
        n_dimensions
    )
    draw <- prior_draw(
        compiled$structure,
        shrinkage,
        n_dimensions
    )
    list(
        value = draw$value,
        global_precision = shrinkage$global,
        local_precision = shrinkage$local,
        raw = draw$raw
    )
}

prior_draw_shrinkage <- function(
    specification,
    n_entities,
    n_dimensions
) {
    type <- specification$type
    result <- list(global = rep(1, n_dimensions), local = NULL)
    if (type == "fixed") {
        result$global <- rep(specification$precision, n_dimensions)
    } else if (type == "gamma") {
        result$global <- stats::rgamma(
            n_dimensions,
            shape = specification$shape,
            rate = specification$rate
        )
    } else if (type == "global_half_cauchy") {
        result$global <- 1 / rhcauchy(n_dimensions, specification$scale)^2
    } else if (type == "local_half_cauchy") {
        result$local <- matrix(
            1 / rhcauchy(n_entities * n_dimensions, specification$scale)^2,
            nrow = n_entities,
            ncol = n_dimensions
        )
    } else if (type == "horseshoe") {
        result$global <- 1 / rhcauchy(n_dimensions, specification$global_scale)^2
        result$local <- matrix(
            1 / rhcauchy(n_entities * n_dimensions, specification$local_scale)^2,
            nrow = n_entities,
            ncol = n_dimensions
        )
    } else if (type == "multiplicative_gamma") {
        delta <- stats::rgamma(
            n_dimensions,
            shape = specification$shape,
            rate = specification$rate
        )
        result$global <- cumprod(delta)
    } else {
        cli::cli_abort("Unsupported shrinkage prior: {type}")
    }
    result
}

prior_draw_structure <- function(
    structure,
    shrinkage,
    n_dimensions
) {
    n_nodes <- length(structure$nodes)
    if (structure$iid) {
        precision <- matrix(
            shrinkage$global,
            nrow = n_nodes,
            ncol = n_dimensions,
            byrow = TRUE
        )
        if (!is.null(shrinkage$local)) {
            precision <- precision * shrinkage$local
        }
        return(matrix(
            stats::rnorm(
                n_nodes * n_dimensions,
                sd = 1 / sqrt(as.numeric(precision))
            ),
            nrow = n_nodes,
            ncol = n_dimensions
        ))
    }
    raw <- matrix(0, nrow = n_nodes, ncol = n_dimensions)
    for (dim in seq_len(n_dimensions)) {
        raw[, dim] <- batchieR_rmvnorm(shrinkage$global[dim] * structure$Q)
    }
    dimnames(raw) <- list(structure$nodes, NULL)
    raw
}

prior_draw <- function(structure, shrinkage, n_dim) {
    raw <- prior_draw_structure(structure, shrinkage, n_dim)
    value <- raw[structure$modeled_index, , drop = FALSE] / structure$modeled_scale
    dimnames(value) <- list(
        structure$entity_names,
        if (n_dim > 1L) paste0("dim_", seq_len(n_dim)) else "value"
    )
    list(value = value, raw = raw)
}
