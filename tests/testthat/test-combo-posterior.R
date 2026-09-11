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

    selected <- posterior_draws(fit, select = c("intercept", "sigma"))
    expect_identical(dimnames(selected)$variable, c("intercept", "sigma"))

    component <- posterior_draws(fit, select = "cell_offset")
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
        posterior_draws(fit, variable = "intercept"),
        "unused argument"
    )
    expect_error(posterior_draws(fit, select = "cell_factors"), "Unmatched names or prefixes")
    expect_error(posterior_draws(fit, select = "unknown"), "Unmatched names or prefixes")
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


test_that("posterior predictive draw selection uses type-specific defaults", {
    yrep <- matrix(seq_len(300), nrow = 100L)

    set.seed(10L)
    expect_identical(nrow(combo_ppc_draws(yrep, "dens_overlay", NULL)), 50L)
    expect_identical(nrow(combo_ppc_draws(yrep, "ecdf_overlay", NULL)), 50L)
    expect_identical(nrow(combo_ppc_draws(yrep, "intervals", NULL)), 100L)
    expect_identical(nrow(combo_ppc_draws(yrep, "stat", 7L)), 7L)
})

test_that("posterior predictive selectors align to observed source rows", {
    fit <- combo_test_fit()
    observed <- combo_observed_data(fit)
    group <- c("a", "a", "b", "b", NA, "b")

    expect_identical(
        combo_ppc_variable("cell", "group", observed),
        fit$input$data$cell[fit$compiled$observed_rows]
    )
    expect_identical(
        combo_ppc_variable(group, "group", observed),
        group[fit$compiled$observed_rows]
    )
    expect_error(combo_ppc_variable("unknown", "group", observed), "not in")
    expect_error(combo_ppc_variable(letters[1:2], "group", observed), "length 6")
    expect_error(combo_ppc_variable("cell", "x", observed, numeric = TRUE), "numeric")
})

test_that("posterior predictive check options reject incompatible inputs", {
    skip_if_not_installed("bayesplot", minimum_version = "1.13.0")
    fit <- combo_test_fit()

    expect_error(
        bayesplot::pp_check(fit, type = "stat_2d", stat = "mean"),
        "two statistics"
    )
    expect_error(
        bayesplot::pp_check(fit, type = "stat", stat = c("mean", "sd")),
        "one statistic"
    )
    expect_error(
        bayesplot::pp_check(fit, type = "stat_2d", group = "cell"),
        "group.*not supported"
    )
    expect_error(
        bayesplot::pp_check(fit, type = "ecdf_overlay", x = "dose_1"),
        "x.*only supported"
    )
    expect_error(
        bayesplot::pp_check(fit, type = "loo_pit", ndraws = 2L),
        "requires every retained draw"
    )
    expect_error(
        bayesplot::pp_check(fit, type = "loo_pit", newdata = fit$input$data),
        "only available for fitted observations"
    )
    expect_error(
        suppressWarnings(bayesplot::pp_check(
            fit, type = "error", stat = function(x) range(x)
        )),
        "one finite numeric"
    )
})

test_that("density checks warn when repeated response mass can mislead", {
    expect_warning(
        combo_warn_ppc_density(c(0, 0, seq_len(20))),
        "repeated values"
    )
    expect_no_warning(combo_warn_ppc_density(seq_len(20)))
})

test_that("extended posterior predictive checks integrate with bayesplot", {
    skip_if_not_installed("bayesplot", minimum_version = "1.13.0")
    fit <- combo_test_fit(iter_sampling = 10L, thin = 1L)

    for (type in c(
        "dens_overlay", "ecdf_overlay", "intervals", "stat",
        "stat_2d", "error"
    )) {
        plot <- expect_no_warning(suppressMessages(
            bayesplot::pp_check(fit, type = type)
        ))
        expect_s3_class(plot, "ggplot")
    }

    grouped <- list(
        ecdf_overlay = list(group = "cell"),
        intervals = list(group = "cell", x = "dose_1"),
        stat = list(group = "cell", stat = "sd"),
        error = list(group = "cell", x = "dose_1")
    )
    for (type in names(grouped)) {
        plot <- suppressMessages(suppressWarnings(do.call(
            bayesplot::pp_check,
            c(list(object = fit, type = type), grouped[[type]])
        )))
        expect_s3_class(plot, "ggplot")
    }
})

test_that("LOO-PIT values use aligned normalized importance weights", {
    y <- c(0, 1)
    yrep <- matrix(c(-1, 0, 1, 2, 0, 1, 2, 3), nrow = 4L)
    weights <- matrix(0.25, nrow = 4L, ncol = 2L)

    expect_equal(combo_loo_pit_values(y, yrep, weights), c(0.5, 0.5))
    expect_error(
        combo_loo_pit_values(y, yrep, weights[, 1L, drop = FALSE]),
        "not aligned"
    )
})

test_that("LOO-PIT returns a calibration plot", {
    skip_if_not_installed("bayesplot", minimum_version = "1.13.0")
    skip_if_not_installed("loo", minimum_version = "2.0.0")
    fit <- combo_test_fit(iter_sampling = 10L, thin = 1L)

    plot <- suppressMessages(suppressWarnings(
        bayesplot::pp_check(fit, type = "loo_pit")
    ))
    expect_s3_class(plot, "ggplot")
})

test_that("prediction observation flags use strict scalar logical validation", {
    fit <- combo_test_fit()
    for (value in list(NULL, logical(), c(TRUE, FALSE), NA, 0, 1L, "TRUE", list(TRUE))) {
        expect_error(predict(fit, observation = value),
            "observation must be TRUE or FALSE", class = "rlang_error")
    }
    expect_equal(predict(fit, observation = FALSE), posterior_epred(fit))
    set.seed(17)
    expected <- posterior_predict(fit)
    set.seed(17)
    expect_identical(predict(fit, observation = TRUE), expected)
})


test_that("selection expands literal prefixes in requested order without duplicates", {
    fit <- combo_test_fit()
    draws <- posterior_draws(fit)
    offsets <- c("cell_offset[1]", "cell_offset[2]", "cell_offset_global_precision[1]")
    expected <- c("sigma", "cell_offset[2]", offsets[c(1, 3)], "intercept")
    select <- c("sigma", "cell_offset[2]", "cell_offset", "intercept", "sigma")
    expect_equal(posterior_draws(fit, select = select), draws[, , expected, drop = FALSE])
    expect_equal(posterior_draws(fit, select = "cell_offset["),
                 draws[, , offsets[1:2], drop = FALSE])
    expect_equal(posterior_draws(fit, select = "sig"), draws[, , "sigma", drop = FALSE])
    expect_identical(summary(fit, select = select)$posterior_summary$variable, expected)
    expect_identical(tidy(fit, select = select)$term, expected)
    for (convert in list(posterior::as_draws, posterior::as_draws_array,
                         posterior::as_draws_matrix, posterior::as_draws_df,
                         posterior::as_draws_list, posterior::as_draws_rvars)) {
        expect_equal(convert(fit, select = select), convert(draws[, , expected, drop = FALSE]))
    }
    for (invalid in list(character(), "", NA_character_, c("sigma", NA_character_),
                         1, TRUE, list("sigma"))) {
        expect_error(posterior_draws(fit, select = invalid), "nonempty character vector")
    }
    expect_error(posterior_draws(fit, select = c("sigma", "missing")), "missing")
    expect_error(posterior_draws(fit, select = "cell_offset.*"), "Unmatched")
    expect_error(posterior_draws(fit, select = "mean"), "Unmatched")
    expect_error(posterior_draws(fit, components = "cell_offset"), "unused argument")
    expect_error(summary(fit, variable = "sigma"), "use.*select")
    expect_error(summary(fit, components = "cell_offset"), "use.*select")
    expect_error(tidy(fit, variable = "sigma"), "use.*select")
    expect_error(tidy(fit, components = "cell_offset"), "use.*select")
})

test_that("prefix selection respects public and expert visibility", {
    fit <- combo_test_fit(model = combo_model(rank = 2L))
    public <- posterior_draws(fit)
    all <- posterior_draws(fit, include = "all")
    expect_equal(posterior_draws(fit, select = "cell_"),
                 public[, , startsWith(dimnames(public)$variable, "cell_"), drop = FALSE])
    expert_names <- dimnames(all)$variable[startsWith(dimnames(all)$variable, "cell_factors")]
    expect_equal(posterior_draws(fit, select = "cell_factors", include = "all"),
                 all[, , expert_names, drop = FALSE])
    expect_error(posterior_draws(fit, select = "cell_factors"), 'include = "all"', fixed = TRUE)
    expect_error(posterior_draws(fit, select = c("missing", "cell_factors")), "missing")
    expect_error(posterior_draws(fit, select = c("missing", "cell_factors")),
                 'include = "all"', fixed = TRUE)
    selected <- c("sigma", "main_effect")
    expected <- c("sigma", dimnames(public)$variable[startsWith(dimnames(public)$variable, "main_effect")])
    expect_equal(posterior_draws(fit, select = selected), public[, , expected, drop = FALSE])
})
