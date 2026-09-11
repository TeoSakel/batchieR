# Plate-level scoring for combination models.

#' Rank candidate plates using Gaussian PDBAL
#'
#' Uses posterior-distance-based active learning (PDBAL) to rank experimental
#' plates for a fitted combination-response model. Posterior draws act as
#' competing response models, while `reference_grid` determines the response
#' surface on which disagreements among those models matter.
#'
#' Let \eqn{\pi_n} be the current posterior and let \eqn{d} be the chosen
#' distance between two posterior hypotheses. Diameter-based active learning
#' aims to make the posterior average diameter small:
#'
#' \deqn{\operatorname{avg-diam}(\pi_n) =
#' \mathbb{E}_{\theta,\theta' \sim \pi_n}[d(\theta,\theta')].}
#'
#' For a candidate plate \eqn{P}, PDBAL uses the regularized objective
#'
#' \deqn{s_n(P) = \mathbb{E}_{\theta,\theta',\theta^\star \sim \pi_n}
#' \left[d(\theta,\theta')
#' L_{\theta^\star}(\theta,\theta';P)
#' \exp\{2H_{\theta^\star}(P)\}\right],}
#'
#' where the likelihood-overlap term is
#'
#' \deqn{L_{\theta^\star}(\theta,\theta';P) =
#' \mathbb{E}_{y_P \sim p_{\theta^\star}(\cdot\mid P)}
#' [p_\theta(y_P\mid P)p_{\theta'}(y_P\mid P)]}
#'
#' and \eqn{H_\theta(P)} is the predictive Shannon entropy of the plate under
#' \eqn{\theta}. Thus the selected plate \eqn{\arg\min_P s_n(P)} is a
#' regularized minimizer of posterior diameter: the distance prioritizes
#' posterior hypotheses that remain far apart, likelihood overlap emphasizes
#' hypotheses the plate would fail to distinguish, and the entropy factor
#' regularizes that overlap for the concentration of the plate predictive
#' distribution.
#'
#' For each candidate plate, the Gaussian approximation:
#'
#' 1. enumerates or samples triplets of posterior draws;
#' 2. measures their separation on transformed reference-grid responses;
#' 3. measures the overlap of their candidate-response distributions; and
#' 4. combines separation, overlap, and predictive entropy into the PDBAL objective.
#'
#' The implementation returns a ranking-equivalent, unnormalized log Monte
#' Carlo approximation to \eqn{s_n(P)}. Lower scores are preferred.
#' `max_triplets` bounds the approximation cost.
#' When subsampling is required it uses R's current random-number generator;
#' call [set.seed()] before this function for reproducible results.
#'
#' `response_transform` is applied only to the posterior expected responses on
#' `reference_grid`. Candidate likelihoods remain Gaussian on the response
#' scale used to fit the model. The built-in `"mean_squared_error"` distance is
#' the pairwise mean squared difference between posterior response vectors.
#' Alternatively, `distance` may be a function taking the transformed
#' draw-by-row prediction matrix and returning a finite, nonnegative, symmetric
#' draw-by-draw matrix with a zero diagonal.
#'
#' This function ports the Gaussian PDBAL scoring approach used by BATCHIE.
#'
#' @param fit A fitted `combo_fit` object with at least three posterior draws.
#' @param candidate_data A non-empty, uniquely named list of candidate-plate
#'   data frames. Each data frame follows the prediction data contract used by
#'   [posterior_epred()]. List names identify plates in the returned scores.
#' @param reference_grid A non-empty prediction data frame defining the response
#'   surface used to compare posterior draws.
#' @param response_transform A function applied to the draw-by-row reference
#'   predictions before computing distances. It must return a finite numeric
#'   matrix with the same dimensions. Use [plogis()] when a model fitted to
#'   logit responses should be compared on the viability scale.
#' @param distance Either `"mean_squared_error"` or a function accepting the
#'   transformed reference prediction matrix and returning a draw-by-draw
#'   distance matrix.
#' @param max_triplets Positive integer limiting the posterior-draw triplets.
#' @param distance_factor Positive finite exponent applied to posterior-model
#'   distances in the PDBAL objective.
#'
#' @return A data frame with one row per candidate plate, sorted from best to
#'   worst. `plate` is the corresponding name from `candidate_data`;
#'   `n_experiments` is the number of rows on that plate; `score` is the
#'   ranking-equivalent, unnormalized log PDBAL objective; `rank` is its
#'   ascending rank, with input order breaking ties.
#'
#' @references
#' Tosh, C. et al. (2025). "A Bayesian active learning platform for scalable
#' combination drug screens." *Nature Communications*, 16, 156.
#' \doi{10.1038/s41467-024-55287-7}.
#'
#' @examples
#' \dontrun{
#' candidates <- list(
#'     confirmatory = candidate_rows[c(1, 2, 3), ],
#'     exploratory = candidate_rows[c(4, 5, 6), ]
#' )
#'
#' set.seed(42)
#' result <- score_plate_pdbal(
#'     fit,
#'     candidate_data = candidates,
#'     reference_grid = reference_rows,
#'     response_transform = plogis
#' )
#'
#' mean_absolute_distance <- function(predictions) {
#'     as.matrix(stats::dist(predictions, method = "manhattan")) /
#'         ncol(predictions)
#' }
#' custom_result <- score_plate_pdbal(
#'     fit,
#'     candidates,
#'     reference_rows,
#'     distance = mean_absolute_distance
#' )
#' }
#' @export
score_plate_pdbal <- function(
    fit,
    candidate_data,
    reference_grid,
    response_transform = identity,
    distance = "mean_squared_error",
    max_triplets = 5000L,
    distance_factor = 1
) {
    validate_combo_fit(fit)
    if (!identical(fit$model$family$name, "gaussian") ||
            !identical(fit$model$family$link, "identity")) {
        cli::cli_abort(
            "score_plate_pdbal() requires a Gaussian identity-link fit"
        )
    }
    if (!is.list(candidate_data) || !length(candidate_data) ||
            is.null(names(candidate_data)) ||
            is_invalid_key(names(candidate_data))) {
        cli::cli_abort(
            "candidate_data must be a non-empty, uniquely named list"
        )
    }
    invisible(lapply(
        candidate_data,
        validate_scoring_data,
        label = "Each candidate_data element"
    ))
    validate_scoring_data(reference_grid, "reference_grid")
    if (!is.function(response_transform)) {
        cli::cli_abort("response_transform must be a function")
    }
    if (is.character(distance)) {
        if (length(distance) != 1L || is.na(distance) ||
                !identical(distance, "mean_squared_error")) {
            cli::cli_abort(
                'distance must be "mean_squared_error" or a function'
            )
        }
    } else if (!is.function(distance)) {
        cli::cli_abort(
            'distance must be "mean_squared_error" or a function'
        )
    }
    max_triplets <- param_positive_integer(max_triplets, "max_triplets")
    distance_factor <- param_scalar_positive(distance_factor, "distance_factor")

    reference_predictions <- posterior_epred(fit, newdata = reference_grid)
    transformed_predictions <- response_transform(reference_predictions)
    if (!is.matrix(transformed_predictions) ||
            !is.numeric(transformed_predictions) ||
            !identical(dim(transformed_predictions), dim(reference_predictions)) ||
            any(!is.finite(transformed_predictions))) {
        cli::cli_abort(
            "response_transform must return a finite numeric matrix with unchanged dimensions"
        )
    }
    distance_matrix <- if (is.character(distance)) {
        pairwise_mean_squared_error(transformed_predictions)
    } else {
        distance(transformed_predictions)
    }
    distance_matrix <- validate_distance_matrix(
        distance_matrix,
        nrow(reference_predictions)
    )

    means_by_plate <- lapply(
        candidate_data,
        function(candidate) posterior_epred(fit, newdata = candidate)
    )
    precision_draws <- posterior_draws(
        fit,
        select = "observation_precision"
    )
    precisions <- as.numeric(precision_draws[, , "observation_precision"])
    if (length(precisions) != nrow(reference_predictions) ||
            any(!is.finite(precisions)) || any(precisions <= 0)) {
        cli::cli_abort(
            "All posterior observation precisions must be finite and positive"
        )
    }
    variances_by_plate <- lapply(
        candidate_data,
        function(candidate) {
            matrix(
                1 / precisions,
                nrow = length(precisions),
                ncol = nrow(candidate)
            )
        }
    )

    scores <- pdbal_gaussian_scores(
        means_by_plate = means_by_plate,
        variances_by_plate = variances_by_plate,
        distance_matrix = distance_matrix,
        max_triplets = max_triplets,
        distance_factor = distance_factor
    )
    score_table <- data.frame(
        plate = names(scores),
        n_experiments = vapply(candidate_data, nrow, integer(1L)),
        score = unname(scores),
        stringsAsFactors = FALSE
    )
    score_table$rank <- rank(score_table$score, ties.method = "first")
    score_table <- score_table[order(score_table$rank), ]
    rownames(score_table) <- NULL
    score_table
}

pdbal_gaussian_scores <- function(means_by_plate, variances_by_plate,
                                  distance_matrix, max_triplets = 5000L,
                                  distance_factor = 1) {
    if (!is.list(means_by_plate) || !length(means_by_plate) ||
            is.null(names(means_by_plate)) ||
            is_invalid_key(names(means_by_plate))) {
        cli::cli_abort(
            "means_by_plate must be a non-empty, uniquely named list"
        )
    }
    if (!is.list(variances_by_plate) ||
            !identical(names(variances_by_plate), names(means_by_plate))) {
        cli::cli_abort(
            "variances_by_plate must have the same names as means_by_plate"
        )
    }
    max_triplets <- param_positive_integer(max_triplets, "max_triplets")
    distance_factor <- param_scalar_positive(distance_factor, "distance_factor")

    distance_matrix <- as.matrix(distance_matrix)
    if (!is.numeric(distance_matrix) ||
            nrow(distance_matrix) != ncol(distance_matrix)) {
        cli::cli_abort("distance_matrix must be a square numeric matrix")
    }
    if (nrow(distance_matrix) < 3L) {
        cli::cli_abort("PDBAL requires at least three posterior draws")
    }
    distance_matrix <- validate_distance_matrix(
        distance_matrix,
        nrow(distance_matrix)
    )
    n_draws <- nrow(distance_matrix)

    for (plate_name in names(means_by_plate)) {
        means <- as.matrix(means_by_plate[[plate_name]])
        variances <- as.matrix(variances_by_plate[[plate_name]])
        if (!is.numeric(means) || !is.numeric(variances) ||
                !identical(dim(means), dim(variances)) ||
                nrow(means) != n_draws || ncol(means) < 1L) {
            cli::cli_abort(
                "Each mean and variance matrix must have matching draw-by-experiment dimensions"
            )
        }
        if (any(!is.finite(means)) || any(!is.finite(variances)) ||
                any(variances <= 0)) {
            cli::cli_abort(
                "Plate means must be finite and variances must be finite and positive"
            )
        }
        means_by_plate[[plate_name]] <- means
        variances_by_plate[[plate_name]] <- variances
    }

    triplets <- draw_triplets(n_draws, max_triplets)
    idx1 <- triplets[, 1L]
    idx2 <- triplets[, 2L]
    idx3 <- triplets[, 3L]
    log_distances <- cbind(
        distance_factor * log(distance_matrix[cbind(idx1, idx2)]),
        distance_factor * log(distance_matrix[cbind(idx2, idx3)]),
        distance_factor * log(distance_matrix[cbind(idx1, idx3)])
    )

    scores <- vapply(names(means_by_plate), function(plate_name) {
        means <- means_by_plate[[plate_name]]
        variances <- variances_by_plate[[plate_name]]
        mean1 <- means[idx1, , drop = FALSE]
        mean2 <- means[idx2, , drop = FALSE]
        mean3 <- means[idx3, , drop = FALSE]
        variance1 <- variances[idx1, , drop = FALSE]
        variance2 <- variances[idx2, , drop = FALSE]
        variance3 <- variances[idx3, , drop = FALSE]

        alpha <- variance1 * variance2 +
            variance2 * variance3 +
            variance1 * variance3
        discrepancy <- variance1 * (mean2 - mean3)^2 +
            variance2 * (mean1 - mean3)^2 +
            variance3 * (mean1 - mean2)^2
        log_normalizer <- rowSums(-0.5 * log(alpha))
        exponential_factor <- 0.5 * variance1 * variance2 *
            variance3 / alpha^2
        log_likelihood <- rowSums(-exponential_factor * discrepancy)

        # Up to a ranking-invariant constant, twice the predictive entropy is
        # the sum of log variances. Each distance is weighted by the entropy
        # of the third draw under its corresponding role assignment.
        log_entropies <- cbind(
            rowSums(log(variance3)),
            rowSums(log(variance1)),
            rowSums(log(variance2))
        )
        log_weighted_distances <- matrixStats::rowLogSumExps(
            log_distances + log_entropies
        )

        matrixStats::logSumExp(
            log_normalizer + log_likelihood + log_weighted_distances
        )
    }, numeric(1L))

    attr(scores, "triplets") <- triplets
    attr(scores, "n_possible") <- attr(triplets, "n_possible")
    scores
}

draw_triplets <- function(n_draws, max_triplets = 5000L) {
    n_draws <- param_positive_integer(n_draws, "n_draws")
    max_triplets <- param_positive_integer(max_triplets, "max_triplets")
    if (n_draws < 3L) {
        cli::cli_abort("PDBAL requires at least three posterior draws")
    }

    n_possible <- choose(n_draws, 3L)
    if (n_possible <= max_triplets) {
        triplets <- t(utils::combn(n_draws, 3L))
    } else {
        ranks <- sample.int(n_possible, max_triplets, replace = FALSE) - 1
        triplets <- t(vapply(ranks, unrank_triplet, integer(3L)))
    }
    colnames(triplets) <- c("draw_1", "draw_2", "draw_3")
    attr(triplets, "n_possible") <- n_possible
    triplets
}

unrank_triplet <- function(rank) {
    result <- integer(3L)
    for (position in 3:1) {
        candidate <- position - 1L
        while (choose(candidate + 1L, position) <= rank) {
            candidate <- candidate + 1L
        }
        result[position] <- candidate + 1L
        rank <- rank - choose(candidate, position)
    }
    result
}

pairwise_mean_squared_error <- function(predictions) {
    predictions <- validate_prediction_matrix(predictions)
    squared_norm <- rowSums(predictions^2)
    distances <- outer(squared_norm, squared_norm, `+`) -
        2 * tcrossprod(predictions)
    # Remove floating-point noise and normalize over the reference grid.
    distances <- pmax(distances, 0)
    diag(distances) <- 0
    distances / ncol(predictions)
}

# Common Utilities for scoring functions ----------------------------------

validate_prediction_matrix <- function(predictions) {
    if (!is.matrix(predictions) || !is.numeric(predictions) ||
            nrow(predictions) < 1L || ncol(predictions) < 1L ||
            any(!is.finite(predictions))) {
        cli::cli_abort(
            "predictions must be a finite numeric matrix with at least one row and column"
        )
    }
    predictions
}

validate_distance_matrix <- function(distance_matrix, n_draws, tolerance = 1e-10) {
    distance_matrix <- as.matrix(distance_matrix)
    if (!is.numeric(distance_matrix) ||
            !identical(dim(distance_matrix), c(n_draws, n_draws)) ||
            any(!is.finite(distance_matrix))) {
        cli::cli_abort(
            "distance must return a finite numeric draw-by-draw matrix"
        )
    }
    if (any(distance_matrix < -tolerance)) {
        cli::cli_abort("distance matrix must be non-negative")
    }
    if (!isTRUE(all.equal(
        distance_matrix,
        t(distance_matrix),
        tolerance = tolerance,
        check.attributes = FALSE
    ))) {
        cli::cli_abort("distance matrix must be symmetric")
    }
    if (any(abs(diag(distance_matrix)) > tolerance)) {
        cli::cli_abort("distance matrix must have a zero diagonal")
    }
    distance_matrix <- pmax(distance_matrix, 0)
    diag(distance_matrix) <- 0
    distance_matrix
}

validate_scoring_data <- function(data, label) {
    if (!inherits(data, "data.frame") || nrow(data) < 1L) {
        cli::cli_abort("{label} must be a data frame with at least one row")
    }
    validate_experiments(data, require_response = FALSE)
    invisible(data)
}