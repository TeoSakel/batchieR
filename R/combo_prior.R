# Prior prediction for combination models.

#' Draw predictions from the prior
#'
#' Draws latent conditional means or replicated responses from the generative prior defined by a [combo_model()].
#'
#' @details
#' Each `beta_offset` entry may be an unnamed numeric scalar, which is recycled,
#' or a named numeric vector matched to compiled coefficient names. Omitted
#' coefficients and list entries default to zero. A nonzero or named entry
#' requires the corresponding component to have coefficients.
#'
#' Formula-generated model-matrix columns are centered and scaled to unit sample
#' standard deviation before coefficients are applied. `beta_offset` therefore
#' operates on this standardized design scale, not on raw covariates. For an
#' ordinary numeric covariate, one coefficient unit is the effect of a one-standard-
#' deviation increase, holding the other design columns fixed.
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
#' @param intercept For a model using [empirical_intercept()], an optional
#'   finite numeric value for the global intercept. It takes precedence over
#'   the mean of nonmissing responses. If neither is available, an error is
#'   raised. This argument must be `NULL` when the model uses
#'   [fixed_intercept()].
#' @param type Whether to draw replicated `response` values including observation
#'   noise, or latent conditional `mean` values.
#' @param beta_offset A named list with optional `cell_offset` and
#'   `treatment_offset` entries giving prior means for coefficients in the
#'   cell-offset and treatment-offset mean formulas, respectively.
#'
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
#'     global_intercept = fixed_intercept(0)
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
    type = c("response", "mean"),
    beta_offset = list(cell_offset = 0, treatment_offset = 0)
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

    alpha <- prior_intercept(model, data, intercept)
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
    beta_offset <- prior_resolve_beta_offsets(compiled, beta_offset)

    if (!is.null(seed)) {
        set.seed(seed)
    }
    indices <- list(
        cell = compiled$all_cell,
        treatment_1 = compiled$all_treatment_1,
        treatment_2 = compiled$all_treatment_2,
        row_id = seq_len(nrow(compiled$data))
    )
    result <- matrix(
        NA_real_,
        nrow = draws,
        ncol = nrow(compiled$data),
        dimnames = list(NULL, as.character(indices$row_id))
    )
    for (draw in seq_len(draws)) {
        snapshot <- prior_snapshot(compiled, alpha, beta_offset)
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
    specification <- model$global_intercept
    if (specification$type == "fixed") {
        if (!is.null(intercept)) {
            cli::cli_abort("intercept must be NULL when the model uses fixed_intercept()")
        }
        return(specification$value)
    }
    if (!is.null(intercept)) {
        if (!is.numeric(intercept) || length(intercept) != 1L ||
                is.na(intercept) || !is.finite(intercept)) {
            cli::cli_abort("intercept must be NULL or one finite numeric value")
        }
        return(as.numeric(intercept))
    }
    if (!inherits(data, "data.frame") || !"response" %in% names(data)) {
        cli::cli_abort(
            "empirical_intercept() requires intercept or at least one observed response"
        )
    }
    response <- data$response
    if (!is.numeric(response)) {
        cli::cli_abort("response must be numeric when used for the empirical intercept")
    }
    observed <- response[!is.na(response)]
    if (!length(observed)) {
        cli::cli_abort(
            "empirical_intercept() requires intercept or at least one observed response"
        )
    }
    if (any(!is.finite(observed))) {
        cli::cli_abort("Nonmissing responses must be finite numeric values")
    }
    mean(observed)
}

prior_resolve_beta_offsets <- function(compiled, beta_offset) {
    offset_components <- c("cell_offset", "treatment_offset")
    if (!is.list(beta_offset)) {
        cli::cli_abort("beta_offset must be a named list")
    }
    offset_names <- names(beta_offset)
    if (length(beta_offset) &&
            (is.null(offset_names) || is_invalid_key(offset_names))) {
        cli::cli_abort("beta_offset must have unique, nonempty names")
    }
    unknown <- setdiff(offset_names, offset_components)
    if (length(unknown)) {
        cli::cli_abort("beta_offset contains unknown entries: {.and {unknown}}")
    }

    result <- lapply(compiled$components, function(component) {
        if (is.null(component)) {
            return(numeric())
        }
        coefficient_names <- colnames(component$mean$X)
        stats::setNames(numeric(length(coefficient_names)), coefficient_names)
    })
    for (component_name in offset_components) {
        value <- if (component_name %in% offset_names) {
            beta_offset[[component_name]]
        } else {
            0
        }
        coefficient_names <- names(result[[component_name]])
        result[[component_name]] <- prior_resolve_beta_offset(
            value,
            coefficient_names,
            component_name
        )
    }
    result
}

prior_resolve_beta_offset <- function(
    value,
    coef_names,
    offset_name
) {
    valid_value <- is.numeric(value) && is.null(dim(value)) &&
        length(value) > 0L && !anyNA(value) && all(is.finite(value))
    if (!valid_value) {
        cli::cli_abort(
            "beta_offset {offset_name} must be a finite numeric scalar or named vector"
        )
    }

    value_names <- names(value)
    if (is.null(value_names)) {
        if (length(value) != 1L) {
            cli::cli_abort("beta_offset {offset_name} must be scalar or have coefficient names")
        }
        if (!length(coef_names) && value != 0) {
            cli::cli_abort(
                "beta_offset {offset_name} is nonzero but its component has no coefficients"
            )
        }
        return(stats::setNames(
            rep(as.numeric(value), length(coef_names)),
            coef_names
        ))
    }
    if (is_invalid_key(value_names)) {
        cli::cli_abort(
            "beta_offset {offset_name} must have unique, nonempty coefficient names"
        )
    }
    unknown <- setdiff(value_names, coef_names)
    if (length(unknown)) {
        cli::cli_abort(
            "beta_offset {offset_name} has unknown coefficients: {.and {unknown}}"
        )
    }
    result <- stats::setNames(numeric(length(coef_names)), coef_names)
    result[value_names] <- as.numeric(value)
    result
}

prior_snapshot <- function(compiled, alpha, beta_offset) {
    precision <- with(
        compiled$observation_prior,
        if (type == "fixed") precision else stats::rgamma(1L, shape = shape, rate = rate)
    )
    component_names <- names(compiled$components)
    components <- lapply(
        component_names,
        function(name) prior_draw_component(compiled$components[[name]], beta_offset[[name]])
    )
    names(components) <- component_names
    list(
        alpha = alpha,
        precision = precision,
        components = components
    )
}

prior_draw_component <- function(compiled, beta_offset) {
    if (is.null(compiled)) {
        return(NULL)
    }
    n_entities <- compiled$n_entities
    n_dimensions <- compiled$n_dimensions
    mean_state <- prior_draw_mean(compiled$mean, beta_offset)
    shrinkage <- prior_draw_shrinkage(
        compiled$shrinkage,
        n_entities,
        n_dimensions
    )
    draw <- prior_draw(
        compiled$structure,
        shrinkage,
        n_dimensions,
        mean_state$value
    )
    list(
        value = draw$value,
        beta = mean_state$beta,
        beta_precision = mean_state$precision,
        global_precision = shrinkage$global,
        local_precision = shrinkage$local,
        raw = draw$raw
    )
}

prior_draw_mean <- function(mean, beta_offset) {
    n_coef <- ncol(mean$X)
    if (!n_coef) {
        return(list(
            beta = numeric(),
            precision = NULL,
            value = numeric(nrow(mean$X))
        ))
    }
    spec <- mean$beta_precision
    precision <- if (spec$type == "fixed") {
        spec$precision
    } else {
        stats::rgamma(1L, shape = spec$shape, rate = spec$rate)
    }

    beta <- stats::rnorm(n_coef, mean = beta_offset, sd = 1 / sqrt(precision))
    names(beta) <- colnames(mean$X)
    list(
        beta = beta,
        precision = precision,
        value = as.numeric(mean$X %*% beta)
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

prior_draw <- function(structure, shrinkage, n_dim, mean = 0) {
    raw <- prior_draw_structure(structure, shrinkage, n_dim)
    value <- raw[structure$modeled_index, , drop = FALSE] / structure$modeled_scale
    value <- value + mean
    dimnames(value) <- list(
        structure$entity_names,
        if (n_dim > 1L) paste0("dim_", seq_len(n_dim)) else "value"
    )
    list(value = value, raw = raw)
}
