test_that("model and component specifications form the public fit contract", {
    model <- combo_test_model()

    expect_s3_class(model, "combo_model")
    expect_s3_class(model$components$cell_offset, "combo_gaussian_component")
    expect_s3_class(model$components$cell_offset, "combo_component")
    expect_s3_class(iid(), "structural_prior")
    expect_s3_class(fixed_scale(), "param_shrinkage")
    expect_s3_class(gamma_precision(), "param_shrinkage")
    expect_s3_class(global_half_cauchy(), "param_shrinkage")
    expect_s3_class(local_half_cauchy(), "param_shrinkage")
    expect_s3_class(horseshoe(), "param_shrinkage")
    expect_s3_class(multiplicative_gamma(), "param_shrinkage")
    expect_s3_class(empirical_intercept(), "param_intercept")
    expect_s3_class(fixed_intercept(0), "param_intercept")
    expect_s3_class(categorical(), "combo_dose")
    expect_s3_class(nested(), "combo_dose")
    expect_s3_class(gaussian(), "combo_family")
    expect_output(print(model), "<combo_model>", fixed = TRUE)
})

test_that("fit_combo retains deterministic chain and thinning metadata", {
    fit <- combo_test_fit()
    repeated <- combo_test_fit()

    expect_s3_class(fit, "combo_fit")
    expect_identical(fit$draws, repeated$draws)
    expect_identical(fit$chain_id, c(1L, 1L, 2L, 2L))
    expect_identical(fit$draw_id, c(1L, 2L, 1L, 2L))
    expect_identical(fit$iteration, c(4L, 6L, 4L, 6L))
    expect_identical(fit$sampling$seed, 812L)
    expect_identical(fit$implementation, "gibbs_iid")
    expect_length(fit$draws, 4L)
    expect_identical(fit$compiled$observed_rows, c(1L, 2L, 3L, 4L, 6L))
    expect_true(is.na(fit$input$data$response[5L]))
})

test_that("fit_combo validates engine and sampling controls", {
    model <- combo_test_model()
    data <- combo_test_data()

    expect_error(fit_combo(model, data, engine = "other"), 'engine must be "gibbs"', fixed = TRUE)
    expect_error(fit_combo(model, data, chains = 0), "integer-valued")
    expect_error(fit_combo(model, data, iter_warmup = -1), "integer-valued")
    expect_error(fit_combo(model, data, iter_sampling = 1, thin = 2), "integer-valued")
    expect_error(fit_combo(model, data, thin = 1.5), "integer-valued")
    expect_error(fit_combo(model, data, seed = 0), "seed must be NULL")
    expect_error(fit_combo(model, data, control = 1), "named list")
    expect_error(fit_combo(model, data, control = list(extra = 1)), "Unused Gibbs")
    expect_error(
        fit_combo(model, transform(data, response = NA_real_)),
        "at least one observed response"
    )
    expect_error(
        fit_combo(list(), data),
        "model must be constructed by combo_model"
    )
})

test_that("fit_combo selects sparse updates for structured components", {
    Q <- Matrix::Diagonal(2L)
    dimnames(Q) <- list(c("A", "B"), c("A", "B"))
    model <- combo_test_model(precision(Q))
    fit <- combo_test_fit(model = model, chains = 1L)

    expect_identical(fit$implementation, "gibbs_sparse")
})

test_that("combo fit print and summary methods report retained results", {
    fit <- combo_test_fit()
    summary_fit <- summary(fit)

    expect_s3_class(summary_fit, "summary.combo_fit")
    expect_identical(summary_fit$n_observed, 5L)
    expect_identical(summary_fit$n_draws, 4L)
    expect_length(summary_fit$rmse, 3L)
    expect_length(summary_fit$observation_precision, 3L)
    expect_output(print(fit), "retained draws: 4", fixed = TRUE)
    expect_output(print(summary_fit), "Combination-model Gibbs fit", fixed = TRUE)
    expect_invisible(print(fit))
    expect_invisible(print(summary_fit))
})
test_that("the default latent-factor model fits every component", {
    fit <- combo_test_fit(
        model = combo_model(rank = 1L),
        chains = 1L,
        iter_warmup = 1L,
        iter_sampling = 2L,
        thin = 1L
    )

    expect_true(all(vapply(fit$draws[[1L]]$components, Negate(is.null), logical(1))))
    expect_identical(dim(posterior_epred(fit)), c(2L, 6L))
    map <- parameter_map(fit)
    expect_true(any(map$component == "cell_factors"))
    expect_true(any(map$component == "treatment_main_factors"))
    expect_true(any(map$component == "treatment_interaction_factors"))
})

test_that("fit_combo compiles cell and compound metadata means", {
    mean_component <- function() {
        combo_gaussian_component(
            mean = ~ 0 + feature,
            mean_shrinkage = fixed_scale(),
            shrinkage = fixed_scale()
        )
    }
    model <- combo_model(
        rank = 1L,
        cell_offset = mean_component(),
        cell_factors = NULL,
        treatment_offset = mean_component(),
        treatment_main_factors = NULL,
        treatment_interaction_factors = NULL
    )
    cell_data <- data.frame(cell = c("A", "B"), feature = c(-1, 1))
    compound_data <- data.frame(drug = c("X", "Y"), feature = c(0, 1))
    fit <- fit_combo(
        model,
        combo_test_data(),
        cell_data = cell_data,
        compound_data = compound_data,
        chains = 1L,
        iter_warmup = 1L,
        iter_sampling = 2L,
        seed = 29L
    )

    expect_identical(fit$input$cell_data, cell_data)
    expect_identical(fit$input$compound_data, compound_data)
    expect_named(fit$draws[[1L]]$components$cell_offset$beta, "feature")
    expect_named(fit$draws[[1L]]$components$treatment_offset$beta, "feature")
})
