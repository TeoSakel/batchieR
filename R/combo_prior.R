# Prior prediction for combination models.

#' Draw prior-predictive responses
#'
#' Draws replicated responses from the generative prior defined by a [combo_model()].
#' The design and metadata inputs follow [fit_combo()], but a `response` column is optional.
#' Responses do not otherwise condition the draws; for a model using [empirical_intercept()],
#' their only possible use is to supply the fallback intercept.
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
#'
#' @return A numeric matrix with one row per draw and one column per input design row,
#'         matching the shape and column naming convention of [posterior_predict()].
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
    intercept = NULL
) {
    if (!inherits(model, "combo_model")) {
        cli::cli_abort("model must be constructed by combo_model()")
    }
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

    alpha <- combo_prior_intercept(model, data, intercept)
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
        row_id = seq_len(nrow(compiled$data))
    )
    result <- matrix(
        NA_real_,
        nrow = draws,
        ncol = nrow(compiled$data),
        dimnames = list(NULL, as.character(indices$row_id))
    )
    for (draw in seq_len(draws)) {
        snapshot <- combo_prior_snapshot(compiled, alpha)
        expected <- predict_combo_draw(
            snapshot,
            indices = indices,
            rank = model$rank
        )
        result[draw, ] <- stats::rnorm(
            length(expected),
            mean = expected,
            sd = 1 / sqrt(snapshot$precision)
        )
    }
    result
}

combo_prior_intercept <- function(model, data, intercept) {
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

combo_prior_snapshot <- function(compiled, alpha) {
    precision_prior <- compiled$observation_prior
    components <- lapply(compiled$components, combo_draw_component_prior)
    list(
        alpha = alpha,
        precision = stats::rgamma(
            1L,
            shape = precision_prior$shape,
            rate = precision_prior$rate
        ),
        components = components
    )
}

combo_draw_component_prior <- function(compiled) {
    if (is.null(compiled)) {
        return(NULL)
    }
    n_entities <- compiled$n_entities
    n_dimensions <- compiled$n_dimensions
    mean_state <- combo_draw_mean_prior(compiled$mean)
    shrinkage <- combo_draw_shrinkage_prior(
        compiled$shrinkage,
        n_entities,
        n_dimensions
    )
    raw <- combo_draw_structured_prior(
        compiled$structure,
        shrinkage,
        n_dimensions
    )
    values <- sweep(
        raw[compiled$structure$modeled_index, , drop = FALSE],
        1L,
        compiled$structure$modeled_scale,
        "/"
    )
    values <- values + mean_state$value
    dimnames(values) <- list(
        compiled$structure$entity_names,
        if (n_dimensions > 1L) {
            paste0("dim_", seq_len(n_dimensions))
        } else {
            "value"
        }
    )
    list(
        value = values,
        beta = mean_state$beta,
        mean_precision = mean_state$precision,
        global_precision = shrinkage$global,
        local_precision = shrinkage$local,
        raw = raw
    )
}

combo_draw_mean_prior <- function(mean) {
    n_coef <- ncol(mean$X)
    if (!n_coef) {
        return(list(
            beta = numeric(),
            precision = NULL,
            value = numeric(nrow(mean$X))
        ))
    }
    spec <- mean$shrinkage
    precision <- if (spec$type == "fixed") {
        spec$precision
    } else {
        stats::rgamma(1L, shape = spec$shape, rate = spec$rate)
    }

    beta <- stats::rnorm(n_coef, sd = 1 / sqrt(precision))
    names(beta) <- colnames(mean$X)
    list(
        beta = beta,
        precision = precision,
        value = as.numeric(mean$X %*% beta)
    )
}

combo_draw_shrinkage_prior <- function(
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

combo_draw_structured_prior <- function(
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
