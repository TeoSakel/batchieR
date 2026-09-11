aggregate_test_data <- function() {
    data <- combo_test_data()
    data <- rbind(data, data[3L, ], data[1L, ], data[3L, ], data[1L, ])
    data[7L, c("drug_1", "dose_1", "drug_2", "dose_2")] <-
        data[3L, c("drug_2", "dose_2", "drug_1", "dose_1")]
    data$response[8:9] <- NA_real_
    data$dose_1[9L] <- 2
    data[10L, c("drug_1", "dose_1", "drug_2", "dose_2")] <-
        data[1L, c("drug_2", "dose_2", "drug_1", "dose_1")]
    rownames(data) <- paste0("observation-", rev(seq_len(nrow(data))))
    data
}

aggregate_test_fit <- function(model = combo_model(rank = 2L), data = aggregate_test_data(),
                               chains = 2L, iter_sampling = 12L) {
    fit_combo(model, data, chains = chains, iter_warmup = 2L,
              iter_sampling = iter_sampling, thin = 1L, seed = 914L, refresh = 0L)
}

# Independent scalar calculation, including absent treatments and disabled terms.
aggregate_expected <- function(fit) {
    combos <- combo_map(fit)
    result <- lapply(c("main", "interaction"), function(kind) {
        t(vapply(fit$draws, function(draw) {
            vapply(seq_len(nrow(combos)), function(id) {
                W <- draw$components$cell_factors
                V <- draw$components[[if (kind == "main") {
                    "treatment_main_factors"
                } else "treatment_interaction_factors"]]
                if (is.null(W) || is.null(V)) return(0)
                take <- function(index) if (is.na(index)) rep(0, fit$model$rank) else V$value[index, ]
                a <- take(combos$treatment_1_id[id])
                b <- take(combos$treatment_2_id[id])
                sum(W$value[combos$cell_id[id], ] * if (kind == "main") (a + b) else (a * b))
            }, numeric(1))
        }, numeric(nrow(combos))))
    })
    stats::setNames(result, c("main", "interaction"))
}

expect_aggregate_draws <- function(fit) {
    expected <- aggregate_expected(fit)
    map <- parameter_map(fit)
    for (kind in c("main", "interaction")) {
        selected <- map[map$component == paste0(kind, "_effect"), ]
        if (!nrow(selected)) next
        actual <- posterior_draws(fit, select = selected$variable)
        expect_equal(as.numeric(actual), as.numeric(expected[[kind]][, selected$entity_index, drop = FALSE]),
                     tolerance = 1e-12)
    }
    rows <- row_map(fit)
    baseline <- fit
    for (i in seq_along(baseline$draws)) {
        baseline$draws[[i]]$components$cell_factors <- NULL
        baseline$draws[[i]]$components$treatment_main_factors <- NULL
        baseline$draws[[i]]$components$treatment_interaction_factors <- NULL
    }
    reconstructed <- posterior_epred(baseline) +
        expected$main[, rows$combo_id, drop = FALSE] +
        expected$interaction[, rows$combo_id, drop = FALSE]
    expect_equal(posterior_epred(fit), reconstructed, tolerance = 1e-12)
}

test_that("index maps distinguish entities, combinations, and observations", {
    fit <- aggregate_test_fit()
    original <- fit
    rng <- .Random.seed
    cells <- cell_map(fit)
    treatments <- treatment_map(fit)
    combos <- combo_map(fit)
    rows <- row_map(fit)
    expect_named(cells, c("cell_id", "cell"))
    expect_named(treatments, c("treatment_id", "drug", "dose"))
    expect_named(combos, c("combo_id", "cell_id", "cell", "treatment_1_id", "drug_1", "dose_1",
                           "treatment_2_id", "drug_2", "dose_2"))
    expect_named(rows, c("row_id", "combo_id", "observed"))
    expect_identical(cells$cell, fit$compiled$cells)
    expect_identical(treatments$drug, fit$compiled$treatments$drug)
    expect_identical(treatments$dose, fit$compiled$treatments$dose)
    expect_identical(rows$row_id, 1:10)
    expect_identical(rows$combo_id, c(1:6, 3L, 1L, 7L, 1L))
    expect_identical(rows$observed, !is.na(fit$input$data$response))
    expect_identical(combos$combo_id, 1:7)
    expect_identical(combos$cell, cells$cell[combos$cell_id])
    expect_true(all(combos$treatment_1_id <= combos$treatment_2_id, na.rm = TRUE))
    expect_true(all(is.na(combos$drug_2[is.na(combos$treatment_2_id)])))
    x_dose <- if (combos$drug_1[7L] == "X") combos$dose_1[7L] else combos$dose_2[7L]
    expect_identical(x_dose, "2")
    expect_identical(fit, original)
    expect_identical(.Random.seed, rng)
    fit$draws <- NULL
    expect_identical(combo_map(fit), combos)
    expect_identical(row_map(fit), rows)
    for (fun in list(cell_map, treatment_map, combo_map, row_map)) {
        expect_error(fun(list()), "fit must be a combo_fit")
    }
})

test_that("public contributions replace factors while all preserves expert draws", {
    fit <- aggregate_test_fit()
    map <- parameter_map(fit)
    all_map <- parameter_map(fit, include = "all")
    draws <- posterior_draws(fit)
    all_draws <- posterior_draws(fit, include = "all")
    expect_identical(dimnames(draws)$variable, map$variable)
    expect_identical(dimnames(all_draws)$variable, all_map$variable)
    expect_false(any(grepl("_factors$", map$component)))
    main <- subset(map, component == "main_effect")
    interaction <- subset(map, component == "interaction_effect")
    expect_identical(main$variable, paste0("main_effect[", 1:7, "]"))
    expect_identical(interaction$variable, paste0("interaction_effect[", c(3L, 6L, 7L), "]"))
    expect_true(all(main$entity_type == "combo" & main$parameter == "contribution"))
    expect_true(all(is.na(main$dimension) & is.na(main$feature)))
    expect_identical(main$entity_index, 1:7)
    expect_equal(draws, all_draws[, , map$variable, drop = FALSE])
    for (name in names(fit$draws[[1L]]$components)) {
        for (parameter in c("value", "global_precision", "local_precision", "raw")) {
            selected <- all_map$variable[all_map$component == name & all_map$parameter == parameter]
            if (!length(selected)) next
            expected <- t(vapply(fit$draws, function(x) as.numeric(x$components[[name]][[parameter]]),
                                 numeric(length(selected))))
            expect_equal(as.numeric(posterior_draws(fit, select = selected, include = "all")),
                         as.numeric(expected))
        }
    }
    expect_error(posterior_draws(fit, select = "cell_factors"), 'include = "all"', fixed = TRUE)
    expect_error(summary(fit, select = "cell_factors[1,1]"), 'include = "all"', fixed = TRUE)
    expect_aggregate_draws(fit)
})

test_that("contributions support every shrinkage family and enabled factor subset", {
    specifications <- list(fixed_scale(), gamma_precision(), global_half_cauchy(),
                           local_half_cauchy(), horseshoe(), multiplicative_gamma())
    for (spec in specifications) {
        component <- combo_gaussian_component(shrinkage = spec)
        model <- combo_model(rank = 2L, cell_factors = component,
                             treatment_main_factors = component, treatment_interaction_factors = component)
        expect_aggregate_draws(aggregate_test_fit(model, chains = 1L, iter_sampling = 2L))
    }
    for (model in list(combo_model(rank = 1L),
                       combo_model(rank = 2L, treatment_interaction_factors = NULL),
                       combo_model(rank = 2L, treatment_main_factors = NULL), combo_test_model())) {
        fit <- aggregate_test_fit(model, chains = 1L, iter_sampling = 2L)
        expect_aggregate_draws(fit)
        map <- parameter_map(fit)
        expect_identical(any(map$component == "main_effect"), !is.null(model$components$treatment_main_factors))
        expect_identical(any(map$component == "interaction_effect"), !is.null(model$components$treatment_interaction_factors))
    }
    singles <- aggregate_test_data()[c(1L, 2L, 4L, 5L), ]
    fit <- aggregate_test_fit(data = singles, chains = 1L, iter_sampling = 2L)
    expect_false(any(parameter_map(fit)$component == "interaction_effect"))
    expect_aggregate_draws(fit)
})

test_that("structured and nested factors use modeled entity values", {
    hierarchy <- data.frame(node = c("root", "A", "B"), parent = c(NA, "root", "root"))
    Q <- Matrix::Matrix(matrix(c(2, -.5, -.5, 2), 2), sparse = TRUE,
                        dimnames = list(c("A", "B"), c("A", "B")))
    for (structure in list(tree(node, parent, hierarchy), precision(Q), gmrf(Q))) {
        component <- combo_gaussian_component(gamma_precision(), structure = structure)
        treatment <- combo_gaussian_component(gamma_precision())
        for (dose in list(categorical(), nested())) {
            model <- combo_model(rank = 2L, cell_factors = component, dose = dose,
                                 treatment_offset = treatment, treatment_main_factors = treatment,
                                 treatment_interaction_factors = treatment)
            fit <- aggregate_test_fit(model, chains = 1L, iter_sampling = 2L)
            expect_identical(fit$implementation, "gibbs_sparse")
            expect_aggregate_draws(fit)
        }
    }
})

test_that("observation covariates do not split shared factor contributions", {
    data <- aggregate_test_data()
    data$plate <- seq_len(nrow(data))
    fit <- aggregate_test_fit(combo_model(rank = 2L, mean = formula_mean(~ plate)), data)
    expect_identical(row_map(fit)$combo_id[c(3L, 7L)], c(3L, 3L))
    predictions <- posterior_epred(fit)
    beta <- vapply(fit$draws, function(x) x$beta[["plate"]], numeric(1))
    expect_equal(unname(predictions[, 7L] - predictions[, 3L]), beta * 4, tolerance = 1e-12)
    expect_aggregate_draws(fit)
})

test_that("aggregate summaries retain posterior diagnostics and symmetry invariance", {
    fit <- aggregate_test_fit(iter_sampling = 24L)
    rng <- .Random.seed
    original <- fit
    expected <- aggregate_expected(fit)
    ids <- parameter_map(fit)$entity_index[parameter_map(fit)$component == "interaction_effect"]
    reference <- array(expected$interaction[, ids, drop = FALSE],
                       dim = c(24L, 2L, length(ids)),
                       dimnames = list(NULL, NULL, paste0("interaction_effect[", ids, "]")))
    stats <- posterior::summarize_draws(posterior::as_draws_array(reference),
        c("mean", "sd", "quantile2", "rhat", "ess_bulk", "ess_tail"))
    actual <- summary(fit, select = "interaction_effect")$posterior_summary
    expect_equal(actual, stats, tolerance = 1e-10)
    tidy_result <- tidy(fit, select = "interaction_effect")
    expect_equal(tidy_result$std.error, stats$sd)
    expect_equal(tidy_result$rhat, stats$rhat)
    expect_equal(tidy_result$ess_bulk, stats$ess_bulk)
    for (operation in c("sign", "scale", "permutation")) {
        changed <- fit
        for (i in seq_along(changed$draws)) {
            if (changed$chain_id[i] != 2L) next
            for (name in c("cell_factors", "treatment_main_factors", "treatment_interaction_factors")) {
                multiplier <- if (operation == "scale") {
                    switch(name, cell_factors = 4, treatment_main_factors = .25,
                           treatment_interaction_factors = .5)
                } else if (operation == "sign" && name == "treatment_interaction_factors") -1 else 1
                for (field in c("value", "raw")) {
                    x <- changed$draws[[i]]$components[[name]][[field]]
                    if (operation == "permutation") x[,] <- x[, 2:1, drop = FALSE]
                    changed$draws[[i]]$components[[name]][[field]] <- x * multiplier
                }
            }
        }
        expect_equal(posterior_draws(changed), posterior_draws(fit), tolerance = 1e-12)
        expect_equal(summary(changed)$posterior_summary, summary(fit)$posterior_summary, tolerance = 1e-9)
        expect_equal(glance(changed), glance(fit), tolerance = 1e-9)
    }
    diagnostics <- posterior::summarize_draws(posterior::as_draws_array(fit),
        c("mean", "sd", "quantile2", "rhat", "ess_bulk", "ess_tail"))
    overall <- glance(fit)
    expect_identical(overall$n_variables, nrow(diagnostics))
    expect_equal(overall$rhat_max, max(diagnostics$rhat, na.rm = TRUE))
    expect_output(print(fit), paste0("max R-hat: ", format(overall$rhat_max, digits = 4)), fixed = TRUE)
    expect_identical(.Random.seed, rng)
    expect_identical(fit, original)
})

test_that("selection avoids evaluating unrelated aggregate combinations", {
    fit <- aggregate_test_fit()
    calculate <- combo_factor_contribution
    visited <- list()
    local_mocked_bindings(combo_factor_contribution = function(snapshot, indices, component) {
        visited[[length(visited) + 1L]] <<- indices
        calculate(snapshot, indices, component)
    })
    parameter_map(fit)
    combo_map(fit)
    row_map(fit)
    posterior_draws(fit, select = "sigma")
    posterior_draws(fit, include = "all", select = "cell_factors[1,1]")
    expect_length(visited, 0L)
    posterior_draws(fit, select = "main_effect[3]")
    expect_length(visited, length(fit$draws))
    expect_true(all(vapply(visited, nrow, integer(1)) == 1L))
})

test_that("aggregate diagnostics retain short, constant, and separated chain behavior", {
    short <- aggregate_test_fit(chains = 1L, iter_sampling = 1L)
    result <- summary(short, select = "main_effect")$posterior_summary
    expect_true(all(is.na(result$rhat)))
    fit <- aggregate_test_fit(iter_sampling = 16L)
    for (i in seq_along(fit$draws)) {
        fit$draws[[i]]$components$cell_factors$value[] <- 0
        fit$draws[[i]]$components$treatment_main_factors$value[] <- 1
        fit$draws[[i]]$components$treatment_interaction_factors$value[] <- 1
    }
    constant <- summary(fit, select = c("main_effect", "interaction_effect"))$posterior_summary
    expect_true(all(constant$sd == 0))
    expect_true(all(is.na(constant$rhat)))
    for (i in seq_along(fit$draws)) {
        fit$draws[[i]]$components$cell_factors$value[] <-
            10 * fit$chain_id[i] + sin(fit$draw_id[i])
    }
    expected <- aggregate_expected(fit)$main[, 1L]
    reference <- array(expected, c(16L, 2L, 1L),
                       dimnames = list(NULL, NULL, "main_effect[1]"))
    capture_warnings <- function(expr) {
        warnings <- character()
        value <- withCallingHandlers(expr, warning = function(w) {
            warnings <<- c(warnings, conditionMessage(w))
            invokeRestart("muffleWarning")
        })
        list(value = value, warnings = warnings)
    }
    diagnostics <- capture_warnings(posterior::summarize_draws(posterior::as_draws_array(reference),
        c("mean", "sd", "quantile2", "rhat", "ess_bulk", "ess_tail")))
    result <- capture_warnings(summary(fit, select = "main_effect[1]")$posterior_summary)
    expect_equal(result$value, diagnostics$value)
    expect_identical(result$warnings, diagnostics$warnings)
    expect_gt(result$value$rhat, 1.1)
})

test_that("maps also support treatment-free models and ignore factor-level ordering", {
    data <- data.frame(cell = factor(c("B", "A", "B"), levels = c("A", "B")),
                       drug_1 = NA, dose_1 = NA, drug_2 = NA, dose_2 = NA,
                       response = c(.1, .2, NA_real_))
    model <- combo_model(mean = fixed_mean(0), cell_offset = NULL, cell_factors = NULL,
                         treatment_offset = NULL, treatment_main_factors = NULL,
                         treatment_interaction_factors = NULL)
    fit <- aggregate_test_fit(model, data, chains = 1L, iter_sampling = 2L)
    expect_identical(cell_map(fit)$cell, c("B", "A"))
    expect_identical(treatment_map(fit)$treatment_id, integer())
    expect_identical(row_map(fit)$combo_id, c(1L, 2L, 1L))
    combos <- combo_map(fit)
    expect_true(all(is.na(combos$treatment_1_id) & is.na(combos$treatment_2_id)))
    expect_identical(parameter_map(fit)$variable, c("intercept", "sigma", "observation_precision"))
    expect_aggregate_draws(fit)
})
