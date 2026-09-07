test_that("tidy agrees with posterior and retains canonical parameter metadata", {
    fit <- combo_test_fit(iter_sampling = 40L, thin = 1L)
    draws <- posterior::as_draws_array(fit)
    expected <- posterior::summarize_draws(draws)
    result <- generics::tidy(fit)
    map <- parameter_map(fit)

    expect_s3_class(result, "tbl_df")
    expect_identical(names(result), c(
        "term", "estimate", "std.error", "rhat", "ess_bulk", "ess_tail",
        "conf.low", "conf.high", setdiff(names(map), "variable")
    ))
    expect_identical(result$term, map$variable)
    expect_equal(result$estimate, expected$mean, ignore_attr = TRUE)
    expect_equal(result$std.error, expected$sd, ignore_attr = TRUE)
    for (name in c("rhat", "ess_bulk", "ess_tail")) {
        expect_equal(result[[name]], expected[[name]], ignore_attr = TRUE)
    }
    expect_equal(result[setdiff(names(map), "variable")],
        tibble::as_tibble(map[setdiff(names(map), "variable")]))
    for (j in seq_along(result$term)) {
        bounds <- stats::quantile(draws[, , j], c(0.025, 0.975), names = FALSE)
        expect_equal(c(result$conf.low[j], result$conf.high[j]), bounds)
    }
    expect_equal(result$std.error[result$term == "intercept"], 0)
    expect_true(is.na(result$rhat[result$term == "intercept"]))
})

test_that("robust tidy summaries only change location and spread", {
    fit <- combo_test_fit(iter_sampling = 20L, thin = 1L)
    regular <- tidy(fit, conf.level = 0.8)
    robust <- tidy(fit, robust = TRUE, conf.level = 0.8)
    draws <- posterior_draws(fit)
    expect_equal(robust$estimate, unname(apply(draws, 3L, stats::median)))
    expect_equal(robust$std.error, unname(apply(draws, 3L, stats::mad)))
    unchanged <- setdiff(names(regular), c("estimate", "std.error"))
    expect_identical(robust[unchanged], regular[unchanged])
    expect_equal(robust$conf.low, unname(apply(draws, 3L, stats::quantile, probs = 0.1)))
    expect_equal(robust$conf.high, unname(apply(draws, 3L, stats::quantile, probs = 0.9)))
    expect_false(any(c("conf.low", "conf.high") %in% names(tidy(fit, conf.int = FALSE))))
})

test_that("tidy selections preserve requested order and map alignment", {
    fit <- combo_test_fit()
    selected <- tidy(fit, variable = c("sigma", "intercept"))
    expect_identical(selected$term, c("sigma", "intercept"))
    expect_identical(selected$component, c("observation", "mean"))
    map <- parameter_map(fit)
    selected <- tidy(fit, components = "cell_offset")
    expect_identical(selected$term, map$variable[map$component == "cell_offset"])
    all_map <- parameter_map(fit, include = "all")
    expect_identical(tidy(fit, include = "all")$term, all_map$variable)
    expect_identical(tidy(fit, include = "all", variable = "sampler_rmse")$parameter, "rmse")
    expect_error(tidy(fit, variable = "absent"), "Unknown variables")
    expect_error(tidy(fit, components = "absent"), "Unknown components")
    expect_error(tidy(fit, components = "cell_factors"), "no posterior variables")
    expect_error(tidy(fit, variable = "sigma", components = "cell_offset"), "mutually exclusive")
})

test_that("glance reports dimensions and diagnostics without removing unavailable counts", {
    fit <- combo_test_fit(iter_sampling = 40L, thin = 1L)
    result <- generics::glance(fit)
    diagnostics <- posterior::summarize_draws(posterior::as_draws_array(fit))
    expect_s3_class(result, "tbl_df")
    expect_identical(names(result), c(
        "nobs", "n_cells", "n_treatments", "n_chains", "n_draws", "n_variables",
        "sigma", "sampler_rmse", "rhat_max", "ess_bulk_min", "ess_tail_min",
        "rhat_unavailable", "ess_bulk_unavailable", "ess_tail_unavailable"
    ))
    expect_identical(nrow(result), 1L)
    expect_identical(unname(unlist(result[1:6])), c(5L, 2L, 2L, 2L, 80L, nrow(diagnostics)))
    expect_equal(result$sigma, mean(vapply(fit$draws, function(x) 1 / sqrt(x$precision), numeric(1))))
    expect_equal(result$sampler_rmse, mean(vapply(fit$draws, function(x) x$last_rmse, numeric(1))))
    expect_equal(result$rhat_max, max(diagnostics$rhat, na.rm = TRUE))
    expect_equal(result$ess_bulk_min, min(diagnostics$ess_bulk, na.rm = TRUE))
    expect_equal(result$ess_tail_min, min(diagnostics$ess_tail, na.rm = TRUE))
    expect_equal(result$rhat_unavailable, sum(is.na(diagnostics$rhat)))
    expect_equal(result$ess_bulk_unavailable, sum(is.na(diagnostics$ess_bulk)))
    expect_equal(result$ess_tail_unavailable, sum(is.na(diagnostics$ess_tail)))
    expect_identical(combo_broom_extreme(c(NA_real_, Inf), max), Inf)
    expect_identical(combo_broom_extreme(c(NA_real_, -Inf), min), -Inf)
})

test_that("short and single-chain fits preserve posterior behavior", {
    for (chains in c(1L, 2L)) {
        fit <- combo_test_fit(chains = chains, iter_sampling = 1L, thin = 1L)
        result <- tidy(fit)
        expected <- posterior::summarize_draws(posterior::as_draws_array(fit))
        expect_equal(result$rhat, expected$rhat, ignore_attr = TRUE)
        expect_equal(result$ess_bulk, expected$ess_bulk, ignore_attr = TRUE)
        expect_equal(result$ess_tail, expected$ess_tail, ignore_attr = TRUE)
        overall <- glance(fit)
        expect_identical(overall$rhat_max, NA_real_)
        expect_identical(overall$ess_bulk_min, NA_real_)
        expect_identical(overall$ess_tail_min, NA_real_)
        expect_identical(overall$rhat_unavailable, nrow(result))
        if (chains == 1L) expect_true(all(is.na(augment(fit)$.se.fit)))
    }
    fit <- combo_test_fit(chains = 1L, iter_sampling = 20L, thin = 1L)
    capture_warnings <- function(expr) {
        warnings <- character()
        value <- withCallingHandlers(expr, warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
        list(value = value, warnings = warnings)
    }
    result <- capture_warnings(tidy(fit))
    expected <- capture_warnings(posterior::summarize_draws(posterior::as_draws_array(fit)))
    expect_identical(result$warnings, expected$warnings)
    for (name in c("rhat", "ess_bulk", "ess_tail")) {
        expect_equal(result$value[[name]], expected$value[[name]], ignore_attr = TRUE)
    }
})

test_that("augment matches expected-response draws and preserves original data", {
    fit <- combo_test_fit(iter_sampling = 20L, thin = 1L)
    data <- fit$input$data
    data$label <- factor(letters[seq_len(nrow(data))])
    result <- generics::augment(fit, data = tibble::as_tibble(data), interval = "confidence", conf.level = 0.8)
    expected <- posterior_epred(fit)
    expect_s3_class(result, "tbl_df")
    expect_identical(result[names(data)], tibble::as_tibble(data))
    expect_equal(result$.fitted, unname(colMeans(expected)))
    expect_equal(result$.se.fit, unname(apply(expected, 2L, stats::sd)))
    expect_equal(result$.lower, unname(apply(expected, 2L, stats::quantile, probs = 0.1)))
    expect_equal(result$.upper, unname(apply(expected, 2L, stats::quantile, probs = 0.9)))
    expect_equal(result$.resid, data$response - result$.fitted)
    expect_true(is.na(result$.resid[5L]))
    expect_identical(names(augment(fit, se.fit = FALSE)), c(names(fit$input$data), ".fitted", ".resid"))
    expect_false(any(c(".lower", ".upper") %in% names(augment(fit))))
})

test_that("augment handles new, repeated, reordered, and empty prediction rows", {
    fit <- combo_test_fit()
    newdata <- fit$input$data[c(6, 2, 6), ]
    newdata$tag <- c("first", "second", "third")
    result <- augment(fit, newdata = newdata)
    expect_identical(result[names(newdata)], tibble::as_tibble(newdata))
    expect_equal(result$.fitted, unname(colMeans(posterior_epred(fit)))[c(6, 2, 6)])
    expect_equal(result$.resid, newdata$response - result$.fitted)
    newdata$response <- NULL
    result <- augment(fit, newdata = newdata, interval = "confidence")
    expect_false(".resid" %in% names(result))
    one <- augment(fit, newdata = newdata[1, ], interval = "confidence")
    expect_equal(one, result[1, ])
    empty <- augment(fit, newdata = newdata[FALSE, ], interval = "confidence")
    expect_identical(empty, result[FALSE, ])
    empty_response <- augment(fit, newdata = fit$input$data[FALSE, ])
    expect_identical(empty_response, augment(fit)[FALSE, ])
})

test_that("augment rejects altered training data, collisions, and unknown entities", {
    fit <- combo_test_fit()
    data <- fit$input$data
    expect_error(augment(fit, data = data, newdata = data), "only one")
    expect_error(augment(fit, data = data[6:1, ]), "original modeling columns")
    changed <- data
    changed$response[1] <- changed$response[1] + 1
    expect_error(augment(fit, data = changed), "original modeling columns")
    changed <- data
    changed$response <- NULL
    expect_error(augment(fit, data = changed), "missing required columns")
    for (name in c(".fitted", ".se.fit", ".lower", ".upper", ".resid")) {
        changed <- data
        changed[[name]] <- 1
        expect_error(augment(fit, newdata = changed, interval = "confidence"), "generated column")
    }
    changed <- data
    changed$cell[1] <- "absent"
    expect_error(augment(fit, newdata = changed), "entities absent")
    changed <- data
    changed$drug_1[1] <- "absent"
    expect_error(augment(fit, newdata = changed), "entities absent")
    changed <- data
    changed$dose_2[1] <- 1
    expect_error(augment(fit, newdata = changed), "jointly present")
})

test_that("augmentation validates observation covariates and tidy maps mean coefficients", {
    data <- combo_test_data()
    data$covariate <- seq_len(nrow(data))
    model <- combo_test_model()
    model$mean <- formula_mean(~ covariate)
    fit <- fit_combo(model, data, chains = 1L, iter_warmup = 2L,
        iter_sampling = 10L, seed = 812L, refresh = 0L)
    result <- tidy(fit)
    map <- parameter_map(fit)
    expect_identical(result$feature, map$feature)
    expect_true("covariate" %in% result$feature)
    expect_no_error(augment(fit, data = data))
    changed <- data
    changed$covariate[1] <- -10
    expect_error(augment(fit, data = changed), "original modeling columns")
    expect_equal(augment(fit, newdata = changed)$.fitted,
        unname(colMeans(posterior_epred(fit, newdata = changed))))
    changed$covariate <- NULL
    expect_error(augment(fit, data = changed), "original modeling columns")
})

test_that("broom verbs validate controls and preserve fits and RNG state", {
    fit <- combo_test_fit()
    original <- fit
    rng <- .Random.seed
    for (value in list(NA, NULL, 1, c(TRUE, FALSE), "TRUE")) {
        expect_error(tidy(fit, robust = value), "robust.*TRUE or FALSE")
        expect_error(tidy(fit, conf.int = value), "conf.int.*TRUE or FALSE")
        expect_error(augment(fit, se.fit = value), "se.fit.*TRUE or FALSE")
    }
    for (level in list(0, 1, -1, Inf, NA_real_, NULL, c(0.5, 0.9), "0.95", TRUE, list(0.95))) {
        expect_error(tidy(fit, conf.level = level), "conf.level")
        expect_error(augment(fit, conf.level = level), "conf.level")
    }
    for (value in list("prediction", NA_character_, 1, c("confidence", "none"))) {
        expect_error(augment(fit, interval = value), "arg")
    }
    expect_identical(augment(fit, interval = "conf"), augment(fit, interval = "confidence"))
    for (value in list(NULL, "n", c("none", "confidence"))) {
        expect_identical(augment(fit, interval = value), augment(fit))
    }
    for (verb in list(tidy, glance, augment)) {
        expect_warning(result <- verb(fit, typo = TRUE), "Unused arguments",
            class = "rlang_warning")
        expect_identical(result, verb(fit))
    }
    expect_no_error(tidy(fit, robust = TRUE))
    expect_no_error(augment(fit, interval = "confidence"))
    expect_identical(fit, original)
    expect_identical(.Random.seed, rng)
})

test_that("broom and generics dispatch to the package methods", {
    skip_if_not_installed("broom")
    fit <- combo_test_fit()
    expect_identical(broom::tidy(fit), generics::tidy(fit))
    expect_identical(broom::glance(fit), generics::glance(fit))
    expect_identical(broom::augment(fit), generics::augment(fit))
    expect_identical(tidy, generics::tidy)
    expect_identical(glance, generics::glance)
    expect_identical(augment, generics::augment)
})
