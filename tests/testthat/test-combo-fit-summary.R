test_that("fit summaries agree with posterior while preserving legacy fields", {
    fit <- combo_test_fit(iter_sampling = 40L, thin = 1L)
    rng <- .Random.seed
    result <- summary(fit)
    expected <- posterior::summarize_draws(posterior::as_draws_array(fit))

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
    result <- summary(fit, variable = c("sigma", "intercept"))
    expect_identical(result$posterior_summary$variable, c("sigma", "intercept"))
    result <- summary(fit, components = "cell_offset")
    expect_equal(result$posterior_summary, posterior::summarize_draws(
        posterior::as_draws_array(fit, components = "cell_offset")
    ))
    result <- summary(fit, "mean", "mcse_mean", include = "all", variable = "sampler_rmse")
    expect_identical(names(result$posterior_summary), c("variable", "mean", "mcse_mean"))
    expect_identical(result$posterior_summary$variable, "sampler_rmse")
    expect_invisible(print(result, n = 1, width = 80))
})

test_that("short, single, constant and separated chains retain diagnostic behavior", {
    for (chains in c(1L, 2L)) {
        fit <- combo_test_fit(chains = chains, iter_sampling = 1L, thin = 1L)
        expect_equal(summary(fit)$posterior_summary,
            posterior::summarize_draws(posterior::as_draws_array(fit)))
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
    expect_equal(result, posterior::summarize_draws(posterior::as_draws_array(fit)))
    expect_output(print(fit), "1 unavailable", fixed = TRUE)
    expect_identical(combo_diagnostic_range(c(NA_real_, Inf), max), "Inf (1 unavailable)")
})
