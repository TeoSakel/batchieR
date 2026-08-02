test_that("posterior draws and parameter maps share a stable schema", {
    fit <- combo_test_fit()
    draws <- posterior_draws(fit)
    map <- parameter_map(fit)

    expect_true(is.array(draws))
    expect_identical(dim(draws)[1:2], c(2L, 2L))
    expect_identical(dimnames(draws)$variable, map$variable)
    expect_named(
        map,
        c(
            "variable", "component", "parameter", "entity_type",
            "entity_index", "entity_key", "entity_label", "dimension",
            "feature"
        )
    )

    selected <- posterior_draws(fit, variable = c("alpha", "sigma"))
    expect_identical(dimnames(selected)$variable, c("alpha", "sigma"))

    component <- posterior_draws(fit, components = "cell_offset")
    expect_true(all(grepl("^cell_offset", dimnames(component)$variable)))

    all_map <- parameter_map(fit, include = "all")
    all_draws <- posterior_draws(fit, include = "all")
    expect_identical(dimnames(all_draws)$variable, all_map$variable)
    expect_true(all(c("sampler_rmse", "sampler_step", "sampler_iteration") %in% all_map$variable))
    expect_gt(nrow(all_map), nrow(map))
})

test_that("posterior draw selection rejects malformed requests and fits", {
    fit <- combo_test_fit()

    expect_error(
        posterior_draws(fit, components = "cell_offset", variable = "alpha"),
        "mutually exclusive"
    )
    expect_error(posterior_draws(fit, components = "unknown"), "Unknown components")
    expect_error(posterior_draws(fit, components = "cell_factors"), "no posterior variables")
    expect_error(posterior_draws(fit, variable = "unknown"), "Unknown variables")
    expect_error(posterior_draws(list()), "fit must be a combo_fit")
    expect_error(parameter_map(list()), "fit must be a combo_fit")

    unbalanced <- fit
    unbalanced$chain_id[4L] <- 1L
    expect_error(posterior_draws(unbalanced), "unbalanced posterior chains")

    inconsistent <- fit
    value <- inconsistent$draws[[2L]]$components$cell_offset$value
    inconsistent$draws[[2L]]$components$cell_offset$value <- rbind(value, value[1L, , drop = FALSE])
    expect_error(posterior_draws(inconsistent), "stable variable schema")
})

test_that("posterior conversion methods support posterior formats", {
    skip_if_not_installed("posterior")
    fit <- combo_test_fit()

    expect_s3_class(posterior::as_draws(fit), "draws_array")
    expect_s3_class(posterior::as_draws_array(fit), "draws_array")
    expect_s3_class(posterior::as_draws_df(fit), "draws_df")
    expect_s3_class(posterior::as_draws_matrix(fit), "draws_matrix")
    expect_s3_class(posterior::as_draws_list(fit), "draws_list")
    expect_s3_class(posterior::as_draws_rvars(fit), "draws_rvars")
})

test_that("prediction methods preserve row mappings and supported semantics", {
    fit <- combo_test_fit()
    linear <- posterior_linpred(fit)
    expected <- posterior_epred(fit)

    expect_identical(dim(linear), c(4L, 6L))
    expect_equal(expected, linear)
    expect_equal(predict(fit), expected)
    expect_equal(predict(fit, type = "mean"), colMeans(expected))

    set.seed(91)
    replicated <- posterior_predict(fit)
    expect_identical(dim(replicated), dim(expected))
    expect_false(isTRUE(all.equal(replicated, expected)))

    single <- combo_test_data()[1L, setdiff(names(combo_test_data()), "response"), drop = FALSE]
    expect_identical(dim(posterior_epred(fit, single)), c(4L, 1L))

    unseen_cell <- single
    unseen_cell$cell <- "C"
    expect_error(posterior_epred(fit, unseen_cell), "cells: C", fixed = TRUE)

    unseen_treatment <- single
    unseen_treatment$drug_1 <- "Z"
    expect_error(posterior_epred(fit, unseen_treatment), "treatments:", fixed = TRUE)

    expect_error(predict(fit, observation = 1), "observation must be TRUE or FALSE")
    expect_error(predict(fit, type = "invalid"), "'arg' should be one of")
})

test_that("log_lik matches the Gaussian density on observed rows", {
    fit <- combo_test_fit()
    observed_rows <- fit$compiled$observed_rows
    observed <- fit$input$data[observed_rows, , drop = FALSE]
    expected <- posterior_epred(fit, newdata = observed)
    draw_sd <- 1 / sqrt(vapply(fit$draws, function(draw) draw$precision, numeric(1)))

    manual <- expected
    for (column in seq_len(ncol(expected))) {
        manual[, column] <- stats::dnorm(
            observed$response[column],
            expected[, column],
            draw_sd,
            log = TRUE
        )
    }

    pointwise <- log_lik(fit)
    expect_equal(as.numeric(pointwise), as.numeric(manual))
    expect_identical(attr(pointwise, "row_ids"), observed_rows)
    expect_identical(colnames(pointwise), as.character(observed_rows))
    expect_identical(attr(log_lik(fit, fit$input$data), "row_ids"), observed_rows)

    missing <- fit$input$data
    missing$response <- NA_real_
    expect_error(log_lik(fit, missing), "no observed responses")
})

test_that("loo and posterior predictive checks integrate with suggested packages", {
    fit <- combo_test_fit()

    skip_if_not_installed("loo")
    pointwise <- log_lik(fit)
    loo_result <- suppressWarnings(
        loo::loo(fit, r_eff = rep(1, ncol(pointwise)))
    )
    expect_s3_class(loo_result, "loo")
    expect_identical(attr(loo_result, "combo_row_ids"), fit$compiled$observed_rows)

    skip_if_not_installed("bayesplot")
    for (type in c("dens_overlay", "ecdf_overlay", "intervals", "stat")) {
        plot <- suppressWarnings(
            bayesplot::pp_check(fit, type = type, ndraws = 4L)
        )
        expect_s3_class(plot, "ggplot")
    }
})
