mean_test_data <- function() {
    data.frame(
        cell = c("A", "B", "C", "A", "B", "C"),
        drug_1 = c("X", "Y", "Z", "Y", "Z", "X"),
        dose_1 = 1,
        drug_2 = c(NA, "X", "Y", "Z", NA, "Y"),
        dose_2 = c(NA, 1, 1, 1, NA, 1),
        response = seq(0.1, 0.6, by = 0.1),
        obs_x = c(-2, -1, 0, 1, 2, 3),
        stringsAsFactors = FALSE
    )
}

mean_test_cell_data <- function() {
    data.frame(
        cell = c("A", "B", "C"),
        cell_x = c(-1, 0, 2),
        cell_group = factor(c("low", "high", "high"), levels = c("low", "high"))
    )
}

mean_test_compound_data <- function() {
    data.frame(
        drug = c("X", "Y", "Z"),
        compound_x = c(2, -1, 0.5),
        compound_group = factor(
            c("old", "new", "new"), levels = c("old", "new")
        )
    )
}

mean_only_fit_model <- function(mean) {
    combo_model(
        rank = 1L, mean = mean,
        cell_offset = NULL, cell_factors = NULL, treatment_offset = NULL,
        treatment_main_factors = NULL,
        treatment_interaction_factors = NULL
    )
}

test_that("global mean constructors validate their natural API", {
    expect_identical(combo_model()$mean$type, "empirical")
    expect_identical(fixed_mean(1.25)$value, 1.25)
    expect_identical(deparse(formula_mean()$formula), "~1")
    expect_error(fixed_mean("one"), "finite number")
    expect_identical(deparse(formula_mean(y ~ x)$formula), "y ~ x")
    expect_error(formula_mean(beta_mean = NA_real_), "finite numeric")
    expect_error(formula_mean(beta_precision = -1), "positive")
    expect_error(formula_mean(beta_precision = c(1, 2)), "must have")

    exports <- getNamespaceExports("batchieR")
    expect_true(all(c(
        "empirical_mean", "fixed_mean", "formula_mean"
    ) %in% exports))
    expect_false(any(c(
        "empirical_intercept", "fixed_intercept"
    ) %in% exports))
})

test_that("formula mean ignores the left-hand side and excludes it from dot", {
    data <- mean_test_data()
    data$excluded <- seq_len(nrow(data))
    model <- mean_only_fit_model(formula_mean(excluded ~ .))

    X <- model.matrix(model, data)

    expect_identical(colnames(X), c("(Intercept)", "obs_x"))
    expect_false("excluded" %in% colnames(X))
})

test_that("one design combines observation, cell, and additive compound terms", {
    data <- mean_test_data()
    cell_data <- mean_test_cell_data()
    compound_data <- mean_test_compound_data()
    model <- mean_only_fit_model(formula_mean(
        ~ obs_x + cell_x + cell_group + compound_x + compound_group
    ))
    X <- model.matrix(model, data, cell_data, compound_data)

    expect_identical(
        colnames(X),
        c(
            "(Intercept)", "obs_x", "cell_x", "cell_grouphigh",
            "compound_x", "compound_groupnew"
        )
    )
    expect_equal(X[, "obs_x"], data$obs_x)
    expect_equal(X[, "cell_x"], cell_data$cell_x[match(data$cell, cell_data$cell)])
    expect_equal(X[, "cell_grouphigh"], c(0, 1, 1, 0, 1, 1))
    compound_row <- match(compound_data$drug, compound_data$drug)
    names(compound_row) <- compound_data$drug
    first <- compound_row[data$drug_1]
    second <- compound_row[data$drug_2]
    compound_x <- compound_data$compound_x[first]
    compound_group <- as.numeric(compound_data$compound_group[first] == "new")
    present <- !is.na(second)
    compound_x[present] <- compound_x[present] +
        compound_data$compound_x[second[present]]
    compound_group[present] <- compound_group[present] +
        as.numeric(compound_data$compound_group[second[present]] == "new")
    expect_equal(X[, "compound_x"], unname(compound_x))
    expect_equal(X[, "compound_groupnew"], unname(compound_group))
})

test_that("mean feature ownership is unambiguous", {
    data <- mean_test_data()
    cell_data <- mean_test_cell_data()
    compound_data <- mean_test_compound_data()

    colliding <- transform(compound_data, cell_x = compound_x)
    expect_error(
        model.matrix(
            mean_only_fit_model(formula_mean(~ cell_x)),
            data, cell_data, colliding
        ),
        "unique across"
    )
    expect_error(
        model.matrix(
            mean_only_fit_model(formula_mean(~ cell_x * compound_x)),
            data, cell_data, compound_data
        ),
        "exactly one"
    )
    expect_error(
        model.matrix(
            mean_only_fit_model(formula_mean(~ unknown)),
            data, cell_data, compound_data
        ),
        "unknown variables"
    )
    expect_error(
        model.matrix(
            mean_only_fit_model(formula_mean(~ cell_x)),
            data, transform(cell_data, response = cell_x), compound_data
        ),
        "reserved observation columns"
    )
    expect_no_error(model.matrix(
        mean_only_fit_model(formula_mean(~ cell_x:cell_group)),
        data, cell_data, compound_data
    ))
})

test_that("factor encodings are retained for prediction", {
    data <- mean_test_data()
    data$batch <- factor(
        c("reference", "other", "reference", "other", "reference", "other"),
        levels = c("reference", "other")
    )
    fit <- fit_combo(
        mean_only_fit_model(formula_mean(~ obs_x + batch)),
        data,
        chains = 1L, iter_warmup = 5L, iter_sampling = 5L, seed = 12L
    )
    newdata <- data[2L, setdiff(names(data), "response"), drop = FALSE]

    expect_named(fit$draws[[1L]]$beta, c("obs_x", "batchother"))
    expect_no_error(prediction <- posterior_epred(fit, newdata))
    expect_identical(dim(prediction), c(5L, 1L))
    unseen <- newdata
    unseen$batch <- factor("new", levels = c("reference", "other", "new"))
    expect_error(posterior_epred(fit, unseen), "new level")
})

test_that("named coefficient priors align to formula columns", {
    data <- transform(mean_test_data(), z = obs_x^2)
    vector_model <- mean_only_fit_model(formula_mean(
        ~ 0 + z + obs_x,
        beta_mean = c(obs_x = 7, z = 3),
        beta_precision = c(obs_x = 2, z = 5)
    ))
    vector_compiled <- compile_combo_design(vector_model, data)
    expect_equal(vector_compiled$mean$beta_mean, c(z = 3, obs_x = 7))
    expect_equal(diag(vector_compiled$mean$beta_precision), c(5, 2))

    P <- matrix(c(4, 1, 1, 3), 2L, dimnames = list(
        c("obs_x", "z"), c("obs_x", "z")
    ))
    matrix_model <- mean_only_fit_model(formula_mean(
        ~ 0 + z + obs_x,
        beta_mean = 0,
        beta_precision = Matrix::Matrix(P, sparse = TRUE)
    ))
    matrix_compiled <- compile_combo_design(matrix_model, data)
    expect_equal(
        matrix_compiled$mean$beta_precision,
        P[c("z", "obs_x"), c("z", "obs_x")]
    )

    expect_error(
        model.matrix(mean_only_fit_model(formula_mean(
            ~ obs_x, beta_mean = c(other = 0)
        )), data),
        "exactly the mean coefficients"
    )
})

test_that("saturated fixed-effect blocks warn when residual offsets remain", {
    data <- mean_test_data()
    cell_data <- transform(
        mean_test_cell_data(),
        cell_group = factor(cell, levels = c("A", "B", "C"))
    )
    model <- combo_model(
        rank = 1L,
        mean = formula_mean(~ cell_group),
        cell_offset = combo_gaussian_component(shrinkage = fixed_scale()),
        cell_factors = NULL, treatment_offset = NULL,
        treatment_main_factors = NULL,
        treatment_interaction_factors = NULL
    )

    expect_warning(
        compile_combo_design(model, data, cell_data = cell_data),
        "cell_offset is active"
    )
    expect_no_warning(model.matrix(model, data, cell_data))
})

test_that("the mean Gibbs block matches analytic Gaussian moments", {
    x <- seq(-1.5, 1.5, length.out = 12L)
    data <- data.frame(
        cell = "A", drug_1 = NA, dose_1 = NA,
        drug_2 = NA, dose_2 = NA,
        response = 0.7 + 1.4 * x + sin(seq_along(x)) / 10,
        x = x
    )
    model <- mean_only_fit_model(formula_mean(
        ~ x, beta_mean = c(x = -0.2), beta_precision = c(x = 2.5)
    ))
    compiled <- compile_combo_model(model, data)
    state <- init_gibbs_state(compiled)
    state$precision <- 4
    X <- compiled$mean$observed_design
    prior_P <- diag(c(0, 2.5))
    posterior_P <- 4 * crossprod(X) + prior_P
    posterior_mean <- solve(
        posterior_P,
        4 * crossprod(X, compiled$response) + c(0, -0.5)
    )
    posterior_covariance <- solve(posterior_P)

    set.seed(410L)
    samples <- matrix(NA_real_, 6000L, 2L)
    for (iteration in seq_len(nrow(samples))) {
        state <- gibbs_update_mean(state)
        samples[iteration, ] <- c(state$intercept, state$beta)
    }
    expect_equal(colMeans(samples), as.numeric(posterior_mean), tolerance = 0.015)
    expect_lt(
        max(abs(stats::cov(samples) - posterior_covariance)),
        0.001
    )
})

test_that("seeded fits recover pure regression and formula-plus-offset signals", {
    set.seed(811L)
    x <- seq(-2, 2, length.out = 60L)
    pure_data <- data.frame(
        cell = "A", drug_1 = NA, dose_1 = NA,
        drug_2 = NA, dose_2 = NA,
        response = 1.2 + 1.8 * x + stats::rnorm(length(x), sd = 0.2),
        x = x
    )
    pure_fit <- fit_combo(
        mean_only_fit_model(formula_mean(~ x)), pure_data,
        chains = 1L, iter_warmup = 200L, iter_sampling = 500L, seed = 812L
    )
    expect_equal(
        mean(vapply(pure_fit$draws, function(draw) draw$intercept, numeric(1))),
        1.2, tolerance = 0.15
    )
    expect_equal(
        mean(vapply(pure_fit$draws, function(draw) draw$beta[["x"]], numeric(1))),
        1.8, tolerance = 0.15
    )

    cell_data <- data.frame(
        cell = paste0("C", seq_len(8L)),
        cell_x = seq(-1.5, 1.5, length.out = 8L)
    )
    residual <- c(-0.3, 0.1, 0.2, -0.1, 0.15, -0.2, 0.25, -0.1)
    indexed <- rep(seq_len(8L), each = 8L)
    offset_data <- data.frame(
        cell = cell_data$cell[indexed],
        drug_1 = NA, dose_1 = NA, drug_2 = NA, dose_2 = NA,
        response = 0.5 + 1.3 * cell_data$cell_x[indexed] +
            residual[indexed] + stats::rnorm(length(indexed), sd = 0.15)
    )
    offset_model <- combo_model(
        rank = 1L, mean = formula_mean(~ cell_x),
        cell_offset = combo_gaussian_component(
            shrinkage = fixed_scale(precision = 8)
        ),
        cell_factors = NULL, treatment_offset = NULL,
        treatment_main_factors = NULL,
        treatment_interaction_factors = NULL
    )
    offset_fit <- fit_combo(
        offset_model, offset_data, cell_data = cell_data,
        chains = 1L, iter_warmup = 300L, iter_sampling = 600L, seed = 813L
    )
    expect_equal(
        mean(vapply(offset_fit$draws, function(draw) draw$beta[["cell_x"]], numeric(1))),
        1.3, tolerance = 0.3
    )

    default_fit <- fit_combo(
        mean_only_fit_model(empirical_mean()), pure_data,
        chains = 1L, iter_warmup = 2L, iter_sampling = 3L, seed = 814L
    )
    empirical <- mean(pure_data$response)
    expect_equal(
        vapply(default_fit$draws, function(draw) draw$intercept, numeric(1)),
        rep(empirical, 3L)
    )
})

test_that("global beta variables are exposed under the mean component", {
    data <- transform(mean_test_data(), response = 1 + 0.5 * obs_x)
    fit <- fit_combo(
        mean_only_fit_model(formula_mean(~ obs_x)), data,
        chains = 1L, iter_warmup = 2L, iter_sampling = 3L, seed = 22L
    )
    map <- parameter_map(fit)
    beta_row <- map[map$variable == "beta[1]", ]

    expect_identical(beta_row$component, "mean")
    expect_identical(beta_row$parameter, "coefficient")
    expect_identical(beta_row$feature, "obs_x")
    expect_false(any(map$parameter == "beta_precision"))
    expect_true(all(c("intercept", "beta[1]") %in% map$variable))
})
