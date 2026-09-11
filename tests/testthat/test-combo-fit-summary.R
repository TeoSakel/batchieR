test_that("fit summaries agree with posterior while preserving legacy fields", {
    fit <- combo_test_fit(iter_sampling = 40L, thin = 1L)
    rng <- .Random.seed
    result <- summary(fit)
    expected <- posterior::summarize_draws(posterior::as_draws_array(fit),
        c("mean", "sd", "quantile2", "rhat", "ess_bulk", "ess_tail"))

    expect_equal(result$posterior_summary, expected)
    expect_equal(result$rmse, stats::quantile(
        vapply(fit$draws, function(x) x$last_rmse, numeric(1)),
        c(0.05, 0.5, 0.95)
    ))
    expect_equal(result$observation_precision, stats::quantile(
        vapply(fit$draws, function(x) x$precision, numeric(1)),
        c(0.05, 0.5, 0.95)
    ))
    output <- paste(capture.output(print(fit)), collapse = "\n")
    expect_match(output, paste0("max R-hat: ", format(max(expected$rhat, na.rm = TRUE), digits = 4)), fixed = TRUE)
    expect_match(output, paste0("min bulk ESS: ", format(min(expected$ess_bulk, na.rm = TRUE), digits = 4)), fixed = TRUE)
    expect_match(output, paste0("min tail ESS: ", format(min(expected$ess_tail, na.rm = TRUE), digits = 4)), fixed = TRUE)
    output <- paste(capture.output(print(result)), collapse = "\n")
    expect_match(output, "ess_tail", fixed = TRUE)
    expect_match(output, tail(expected$variable, 1), fixed = TRUE)
    expect_identical(.Random.seed, rng)
})

test_that("summary supports posterior selections and custom measures", {
    fit <- combo_test_fit(iter_sampling = 20L, thin = 1L)
    result <- summary(fit, select = c("sigma", "intercept"))
    expect_identical(result$posterior_summary$variable, c("sigma", "intercept"))
    result <- summary(fit, select = "cell_offset")
    expect_equal(result$posterior_summary, posterior::summarize_draws(
        posterior::as_draws_array(fit, select = "cell_offset"),
        c("mean", "sd", "quantile2", "rhat", "ess_bulk", "ess_tail")
    ))
    result <- summary(fit, "mean", "mcse_mean", include = "all", select = "sampler_rmse")
    expect_identical(names(result$posterior_summary), c("variable", "mean", "mcse_mean"))
    expect_identical(result$posterior_summary$variable, "sampler_rmse")
    expect_invisible(print(result, n = 1, width = 80))
})

test_that("short, single, constant and separated chains retain diagnostic behavior", {
    for (chains in c(1L, 2L)) {
        fit <- combo_test_fit(chains = chains, iter_sampling = 1L, thin = 1L)
        expect_equal(summary(fit)$posterior_summary,
            posterior::summarize_draws(posterior::as_draws_array(fit),
        c("mean", "sd", "quantile2", "rhat", "ess_bulk", "ess_tail")))
        expect_output(print(fit), "unavailable", fixed = TRUE)
        expect_output(print(fit), "max R-hat: NA", fixed = TRUE)
    }
    fit <- combo_test_fit(iter_sampling = 40L, thin = 1L)
    for (i in seq_along(fit$draws)) {
        fit$draws[[i]]$intercept <- 1
        fit$draws[[i]]$precision <- fit$draws[[i]]$precision + 100 * fit$chain_id[i]
    }
    result <- summary(fit)$posterior_summary
    expect_true(is.na(result$rhat[result$variable == "intercept"]))
    expect_gt(result$rhat[result$variable == "observation_precision"], 1.1)
    expect_equal(result, posterior::summarize_draws(posterior::as_draws_array(fit),
        c("mean", "sd", "quantile2", "rhat", "ess_bulk", "ess_tail")))
    expect_output(print(fit), "1 unavailable", fixed = TRUE)
    expect_identical(combo_diagnostic_range(c(NA_real_, Inf), max), "Inf (1 unavailable)")
})


test_that("print and summary report the configured model rank", {
    for (rank in c(1L, 2L)) {
        fit <- combo_test_fit(model = combo_model(rank = rank))
        result <- summary(fit)
        expect_identical(result$rank, rank)
        expect_output(print(fit), paste0("rank: ", rank), fixed = TRUE)
        expect_output(print(result), paste0("rank: ", rank), fixed = TRUE)
    }
})

test_that("robust summaries select median and MAD while retaining quantiles and diagnostics", {
    fit <- combo_test_fit(model = combo_model(rank = 1L), iter_sampling = 20L, thin = 1L)
    original <- fit
    rng <- .Random.seed
    variables <- c("sigma", "main_effect[3]")
    regular <- summary(fit, select = variables)$posterior_summary
    robust <- summary(fit, robust = TRUE, select = variables)$posterior_summary
    expect_named(regular, c("variable", "mean", "sd", "q5", "q95", "rhat", "ess_bulk", "ess_tail"))
    expect_named(robust, c("variable", "median", "mad", "q5", "q95", "rhat", "ess_bulk", "ess_tail"))
    draws <- posterior_draws(fit, select = variables)
    expect_equal(regular$mean, unname(apply(draws, 3L, mean)))
    expect_equal(regular$sd, unname(apply(draws, 3L, stats::sd)))
    expect_equal(robust$median, unname(apply(draws, 3L, stats::median)))
    expect_equal(robust$mad, unname(apply(draws, 3L, stats::mad)))
    retained <- c("variable", "q5", "q95", "rhat", "ess_bulk", "ess_tail")
    expect_identical(regular[retained], robust[retained])
    tidy_robust <- tidy(fit, robust = TRUE, select = variables)
    expect_equal(robust$median, tidy_robust$estimate)
    expect_equal(robust$mad, tidy_robust$std.error)
    expect_identical(summary(fit, robust = TRUE, select = variables, .cores = 1L)$posterior_summary,
                     robust)
    expect_identical(fit, original)
    expect_identical(.Random.seed, rng)
})

test_that("custom summary measures override robust defaults and robust is strictly logical", {
    fit <- combo_test_fit(iter_sampling = 12L, thin = 1L)
    for (robust in c(FALSE, TRUE)) {
        result <- summary(fit, "mean", "median", "sd", "mad", robust = robust,
                          select = "sigma")$posterior_summary
        expect_named(result, c("variable", "mean", "median", "sd", "mad"))
        result <- summary(fit, spread = function(x) diff(range(x)), robust = robust,
                          select = "sigma")$posterior_summary
        expect_named(result, c("variable", "spread"))
    }
    for (invalid in list(NULL, logical(), NA, 1, "TRUE", c(TRUE, FALSE))) {
        expect_error(summary(fit, robust = invalid), "robust.*TRUE or FALSE")
    }
})
