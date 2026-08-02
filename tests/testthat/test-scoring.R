test_that("mean squared error compares posterior response vectors", {
    predictions <- rbind(
        c(0, 0),
        c(1, 1),
        c(2, 0)
    )

    expect_equal(
        pairwise_mean_squared_error(predictions),
        matrix(
            c(0, 1, 2, 1, 0, 1, 2, 1, 0),
            nrow = 3L,
            byrow = TRUE
        )
    )
})

test_that("Gaussian PDBAL scoring agrees with a three-draw known case", {
    distance_matrix <- matrix(1, nrow = 3L, ncol = 3L)
    diag(distance_matrix) <- 0
    scores <- pdbal_gaussian_scores(
        means_by_plate = list(plate = matrix(0, nrow = 3L, ncol = 1L)),
        variances_by_plate = list(plate = matrix(1, nrow = 3L, ncol = 1L)),
        distance_matrix = distance_matrix
    )

    expect_equal(as.numeric(scores), 0.5 * log(3))
    expect_identical(attr(scores, "n_possible"), 1)
    expect_equal(
        unname(attr(scores, "triplets")),
        matrix(c(1L, 2L, 3L), nrow = 1L),
        ignore_attr = TRUE
    )
})

test_that("triplet enumeration and sampling follow R RNG conventions", {
    set.seed(12)
    state <- .Random.seed
    exhaustive <- draw_triplets(4L, max_triplets = 4L)

    expect_identical(.Random.seed, state)
    expect_identical(nrow(exhaustive), 4L)
    expect_identical(attr(exhaustive, "n_possible"), 4)

    set.seed(42)
    first <- draw_triplets(6L, max_triplets = 2L)
    state_after <- .Random.seed
    set.seed(42)
    second <- draw_triplets(6L, max_triplets = 2L)

    expect_identical(first, second)
    expect_identical(.Random.seed, state_after)
    expect_identical(nrow(first), 2L)
    expect_identical(attr(first, "n_possible"), 20)
})

test_that("distance matrix validation enforces the callable contract", {
    expect_error(
        validate_distance_matrix(matrix(0, 2L, 2L), 3L),
        "draw-by-draw"
    )

    asymmetric <- matrix(c(0, 1, 2, 0), 2L, 2L)
    expect_error(
        validate_distance_matrix(asymmetric, 2L),
        "symmetric"
    )

    negative <- matrix(c(0, -1, -1, 0), 2L, 2L)
    expect_error(
        validate_distance_matrix(negative, 2L),
        "non-negative"
    )

    nonzero_diagonal <- matrix(c(1, 0, 0, 0), 2L, 2L)
    expect_error(
        validate_distance_matrix(nonzero_diagonal, 2L),
        "zero diagonal"
    )

    nonfinite <- matrix(c(0, Inf, Inf, 0), 2L, 2L)
    expect_error(
        validate_distance_matrix(nonfinite, 2L),
        "finite numeric"
    )
})

test_that("score_plate_pdbal ranks named candidate data frames", {
    fit <- combo_test_fit()
    columns <- c("cell", "drug_1", "dose_1", "drug_2", "dose_2")
    prediction_data <- combo_test_data()[, columns]
    candidates <- list(
        focused = prediction_data[c(1L, 3L), ],
        broad = prediction_data[c(4L, 6L), ]
    )

    result <- score_plate_pdbal(
        fit,
        candidate_data = candidates,
        reference_grid = prediction_data
    )

    expect_s3_class(result, "data.frame")
    expect_named(
        result,
        c("plate", "n_experiments", "score", "rank")
    )
    expect_setequal(result$plate, names(candidates))
    expect_equal(result$n_experiments, c(2L, 2L))
    expect_true(all(is.finite(result$score)))
    expect_identical(result$rank, c(1L, 2L))
})

test_that("custom distance receives transformed reference predictions", {
    fit <- combo_test_fit()
    columns <- c("cell", "drug_1", "dose_1", "drug_2", "dose_2")
    prediction_data <- combo_test_data()[, columns]
    candidates <- list(plate = prediction_data[1:2, ])
    received <- NULL
    custom_distance <- function(predictions) {
        received <<- predictions
        pairwise_mean_squared_error(predictions)
    }

    result <- score_plate_pdbal(
        fit,
        candidates,
        prediction_data,
        response_transform = plogis,
        distance = custom_distance
    )

    expect_equal(received, plogis(posterior_epred(fit, prediction_data)))
    expect_s3_class(result, "data.frame")
})

test_that("score_plate_pdbal validates its public inputs", {
    fit <- combo_test_fit()
    columns <- c("cell", "drug_1", "dose_1", "drug_2", "dose_2")
    prediction_data <- combo_test_data()[, columns]
    candidates <- list(plate = prediction_data[1:2, ])

    expect_error(
        score_plate_pdbal(posterior_draws(fit), candidates, prediction_data),
        "fit must be a combo_fit"
    )
    expect_error(
        score_plate_pdbal(fit, unname(candidates), prediction_data),
        "uniquely named"
    )
    expect_error(
        score_plate_pdbal(fit, candidates, prediction_data[0, ]),
        "at least one row"
    )
    expect_error(
        score_plate_pdbal(
            fit,
            candidates,
            prediction_data,
            response_transform = 1
        ),
        "response_transform must be a function"
    )
    expect_error(
        score_plate_pdbal(
            fit,
            candidates,
            prediction_data,
            response_transform = function(x) as.numeric(x)
        ),
        "unchanged dimensions"
    )
    expect_error(
        score_plate_pdbal(fit, candidates, prediction_data, distance = "other"),
        "mean_squared_error"
    )
    expect_error(
        score_plate_pdbal(
            fit,
            candidates,
            prediction_data,
            distance = function(x) matrix(0, 2L, 2L)
        ),
        "draw-by-draw"
    )
    expect_error(
        score_plate_pdbal(fit, candidates, prediction_data, max_triplets = 0),
        "positive integer"
    )
    expect_error(
        score_plate_pdbal(fit, candidates, prediction_data, distance_factor = 0),
        "positive number"
    )

    unseen <- prediction_data[1, ]
    unseen$cell <- "unseen"
    expect_error(
        score_plate_pdbal(fit, list(plate = unseen), prediction_data),
        "cells: unseen",
        fixed = TRUE
    )

    invalid_precision <- fit
    invalid_precision$draws[[1L]]$precision <- 0
    expect_error(
        score_plate_pdbal(invalid_precision, candidates, prediction_data),
        "precisions must be finite and positive"
    )

    too_few_draws <- combo_test_fit(
        chains = 1L,
        iter_sampling = 2L,
        thin = 1L
    )
    expect_error(
        score_plate_pdbal(too_few_draws, candidates, prediction_data),
        "at least three posterior draws"
    )

    incompatible <- fit
    incompatible$model$family$link <- "other"
    expect_error(
        score_plate_pdbal(incompatible, candidates, prediction_data),
        "Gaussian identity-link"
    )
})
