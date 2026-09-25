mixture_test_fit <- function(chains = 3L, n = 40L, shifts = rep(0, chains), model = NULL, data = NULL) {
    if (is.null(data)) data <- combo_test_data()
    if (is.null(model)) model <- combo_model(rank = 1L, mean = fixed_mean(0),
        cell_offset = NULL, treatment_offset = NULL, cell_factors = NULL,
        treatment_main_factors = NULL, treatment_interaction_factors = NULL)
    compiled <- compile_combo_model(model, data)
    draws <- lapply(seq_len(chains * n), function(i) {
        chain <- (i - 1L) %/% n + 1L
        draw <- (i - 1L) %% n + 1L
        list(intercept = shifts[chain] + .01 * sin(draw * 1.71), beta = numeric(),
             precision = 1 + .01 * cos(draw * 2.31),
             components = lapply(compiled$components, function(x) NULL))
    })
    structure(list(model = model, compiled = compiled, draws = draws,
        chain_id = rep(seq_len(chains), each = n), draw_id = rep(seq_len(n), chains),
        iteration = rep(seq_len(n) * 2L + 10L, chains),
        sampling = list(chains = chains, thin = 2L, iter_warmup = 10L), input = list(data = data)),
        class = "combo_fit")
}

mixture_test_blocks <- function(means, noise = 0) {
    result <- array(0, c(nrow(means), 4L, ncol(means)))
    for (i in seq_len(nrow(means))) for (b in 1:4) result[i, b, ] <- means[i, ] + noise * (b - 2.5)
    result
}

test_that("public projections match manual chain and contiguous-block averages", {
    fit <- mixture_test_fit(chains = 2, n = 23, shifts = c(0, 2))
    fit$chain_id <- rep(c(8L, 3L), each = 23)
    p <- chain_projections(fit, response_transform = exp)
    mu <- posterior_epred(fit, p$reference_grid)
    expect_equal(p$chains, c(3, 8))
    expect_identical(rownames(p$mean), c("3", "8"))
    expect_equal(dim(p$block_mean), c(2, 4, nrow(p$reference_grid)))
    block <- rep(1:4, c(6, 6, 6, 5))
    manual_distances <- numeric()
    for (c in 1:2) {
        ids <- which(fit$chain_id == p$chains[c])
        expected <- exp(mu[ids, , drop = FALSE])
        expect_equal(unname(p$mean[c, ]), unname(colMeans(expected)))
        bm <- t(vapply(1:4, function(b) colMeans(expected[block == b, , drop = FALSE]),
                       numeric(ncol(expected))))
        expect_equal(p$block_mean[c, , ], unname(bm))
        pairs <- combn(1:4, 2)
        distances <- apply(pairs, 2, function(pair) sqrt(mean((bm[pair[1], ] - bm[pair[2], ])^2)))
        rows <- p$block_distances$chain == p$chains[c]
        expect_equal(p$block_distances$distance[rows], distances)
        expect_equal(p$block_distances$block_1[rows], pairs[1, ])
        expect_equal(p$block_distances$block_2[rows], pairs[2, ])
        manual_distances <- c(manual_distances, distances)
    }
    expect_equal(p$noise_floor, unname(quantile(manual_distances, .95)))
    expect_equal(p$source_draws$block, rep(block, 2))
    expect_equal(p$source_draws$iteration, rep(seq_len(23) * 2 + 10, 2))
    expect_equal(p$chain_draws$excluded, c(0, 0))
    expect_s3_class(hclust(dist(p$mean) / sqrt(ncol(p$mean)), method = "average"), "hclust")
})

test_that("public projection preparation is shared with mixtures and preserves order", {
    fit <- mixture_test_fit()
    reordered <- combo_mixture_subset(fit, c(120:81, 39:1, 80:41))
    p <- chain_projections(reordered)
    expected <- chain_projections(combo_mixture_subset(fit, c(1:39, 41:79, 81:119)))
    expect_equal(p$mean, expected$mean)
    expect_equal(p$block_mean, expected$block_mean)
    expect_equal(p$block_distances, expected$block_distances)
    expect_identical(p$source_draws, expected$source_draws)
    expect_equal(p$chain_draws$excluded, c(0, 1, 1))
    expect_equal(p$chain_draws$used, rep(39, 3))
    # QC must remain independent of PSIS, stacking and clustering.
    local_mocked_bindings(combo_mixture_psis = function(...) stop("PSIS called"),
                          combo_mixture_cluster = function(...) stop("clustering called"))
    expect_equal(chain_projections(fit)$chains, 1:3)
})

test_that("public projections respect reference weights, predictors and missing responses", {
    fit <- mixture_test_fit()
    p <- chain_projections(fit)
    expect_equal(nrow(p$reference_grid), 5)
    grid <- p$reference_grid[c(2, 1, 2), ]
    explicit <- chain_projections(fit, grid)
    expect_identical(explicit$reference_grid, grid)
    expect_equal(explicit$mean, p$mean[, c(2, 1, 2)])
    data <- combo_test_data()
    data <- rbind(data, data[1, ])
    data$x <- seq_len(nrow(data))
    model <- combo_model(rank = 1, mean = formula_mean(~ x),
        cell_offset = NULL, treatment_offset = NULL, cell_factors = NULL,
        treatment_main_factors = NULL, treatment_interaction_factors = NULL)
    fit <- mixture_test_fit(model = model, data = data)
    fit$draws <- lapply(fit$draws, function(draw) { draw$beta <- c(x = 2); draw })
    p <- chain_projections(fit)
    expect_equal(nrow(p$reference_grid), sum(!is.na(data$response)))
    mu <- posterior_epred(fit, p$reference_grid)
    expect_equal(unname(p$mean[1, ]), unname(colMeans(mu[1:40, , drop = FALSE])))
    expect_gt(diff(range(p$mean[1, ])), 0)
})

test_that("public projections validate inputs and preserve RNG on success and error", {
    fit <- mixture_test_fit()
    withr::local_seed(52)
    before <- .Random.seed
    chain_projections(fit)
    expect_identical(.Random.seed, before)
    expect_error(chain_projections(fit, response_transform = function(x) x + runif(1)), "deterministic")
    expect_identical(.Random.seed, before)
    rm(".Random.seed", envir = globalenv())
    chain_projections(fit)
    expect_false(exists(".Random.seed", globalenv(), inherits = FALSE))
    expect_error(chain_projections(fit, response_transform = function(x) x + runif(1)), "deterministic")
    expect_false(exists(".Random.seed", globalenv(), inherits = FALSE))
    expect_error(chain_projections(fit, response_transform = 1), "function")
    expect_error(chain_projections(fit, response_transform = as.vector), "shape-preserving")
    expect_error(chain_projections(fit, response_transform = function(x) x * Inf), "finite")
    expect_error(chain_projections(mixture_test_fit(n = 19)), "20 retained")
    expect_error(chain_projections(mixture_test_fit(chains = 1)), "two chains")
    fit$iteration[2] <- fit$iteration[1]
    expect_error(chain_projections(fit), "unique")
    fit <- mixture_test_fit()
    fit$draws[[1]]$precision <- 0
    expect_error(chain_projections(fit), "positive finite")
})

test_that("public block floors agree with mixture clustering including zero variation", {
    fit <- mixture_test_fit(shifts = c(0, 0, 2))
    p <- chain_projections(fit)
    clustering <- combo_mixture_cluster(p$mean, p$block_mean)
    expect_identical(p$noise_floor, clustering$noise_floor)
    fit$draws <- rep(fit$draws[1], length(fit$draws))
    p <- chain_projections(fit)
    expect_equal(p$block_distances$distance, rep(0, 18))
    expect_identical(p$noise_floor, .Machine$double.eps)
    skip_if_not_installed("loo")
    mixture <- suppressWarnings(predictive_mixture(fit))
    expect_equal(p[names(mixture$projection)], mixture$projection)
    expect_identical(p$noise_floor, mixture$clustering$noise_floor)
})

test_that("basins use deterministic noise-scaled mean-response separation", {
    means <- rbind(c(0, 0), c(.001, .001), c(2, 2))
    blocks <- mixture_test_blocks(means, .001)
    x <- combo_mixture_cluster(means, blocks)
    expect_s3_class(x$distance, "dist")
    expect_equal(as.numeric(x$distance), c(.001, 2, 1.999), tolerance = 1e-12)
    expect_equal(x$tree$height, c(.001, 1.9995), tolerance = 1e-12)
    expect_identical(x$basin, c("B1", "B1", "B2"))
    expect_true(all(c("threshold_1.5", "threshold_3", "threshold_6", "cut_0.5", "cut_2") %in% x$sensitivity$setting))
    zero <- combo_mixture_cluster(matrix(0, 3, 2), mixture_test_blocks(matrix(0, 3, 2)))
    expect_identical(zero$basin, rep("B1", 3))
    expect_identical(zero$noise_floor, .Machine$double.eps)
    rownames(means) <- c("10", "20", "30")
    permutation <- c(3, 1, 2)
    shuffled <- combo_mixture_cluster(means[permutation, ], blocks[permutation, , ])
    expect_identical(shuffled$basin, x$basin)
    expect_equal(shuffled$chains, c(10, 20, 30))
    tie <- combo_mixture_cluster(diag(3), mixture_test_blocks(diag(3)))
    expect_identical(tie$basin, rep("B1", 3))
    noisy <- combo_mixture_cluster(means, mixture_test_blocks(means, 10))
    expect_identical(noisy$basin, rep("B1", 3))
})

test_that("two-chain split includes the threshold boundary", {
    means <- matrix(c(0, 3 * .Machine$double.eps), 2, 1)
    expect_identical(combo_mixture_cluster(means, mixture_test_blocks(means))$basin, c("B1", "B2"))
    means[2, 1] <- 2.99 * .Machine$double.eps
    expect_identical(combo_mixture_cluster(means, mixture_test_blocks(means))$basin, c("B1", "B1"))
    expect_true(max(table(combo_mixture_blocks(23))) - min(table(combo_mixture_blocks(23))) <= 1)
    expect_true(all(diff(combo_mixture_blocks(23)) >= 0))
})

test_that("reference deduplication uses modeled predictors and unordered treatments", {
    data <- combo_test_data()
    data <- rbind(data, data[1, ], data[3, ])
    data$unused <- seq_len(nrow(data))
    data$drug_1[8] <- "Y"
    data$drug_2[8] <- "X"
    fit <- mixture_test_fit(data = data)
    expect_equal(nrow(combo_mixture_reference(fit)), 5L)
    data$x <- seq_len(nrow(data))
    model <- combo_model(rank = 1L, mean = formula_mean(~ x),
        cell_offset = NULL, treatment_offset = NULL, cell_factors = NULL,
        treatment_main_factors = NULL, treatment_interaction_factors = NULL)
    fit <- mixture_test_fit(model = model, data = data)
    expect_equal(nrow(combo_mixture_reference(fit)), 7L)
})

test_that("preparation preserves identities and chronological common prefixes", {
    fit <- mixture_test_fit()
    fit <- combo_mixture_subset(fit, c(81:120, 1:39, 41:80))
    result <- combo_mixture_prepare(fit)
    expect_equal(result$excluded$used, rep(39, 3))
    expect_equal(result$excluded$excluded, c(0, 1, 1))
    expect_equal(result$fit$chain_id, rep(1:3, each = 39))
    expect_equal(result$fit$draw_id, rep(1:39, 3))
    fit$iteration[2] <- fit$iteration[1]
    expect_error(combo_mixture_prepare(fit), "unique")
    expect_error(combo_mixture_prepare(mixture_test_fit(n = 19)), "20 retained")
    expect_error(combo_mixture_prepare(mixture_test_fit(chains = 1)), "two chains")
    fit <- mixture_test_fit()
    fit$draws[[1]]$precision <- 0
    expect_error(combo_mixture_prepare(fit), "positive finite")
})

test_that("transforms precede averages and validate elementwise repeatability", {
    fit <- mixture_test_fit()
    grid <- combo_mixture_reference(fit)
    p <- combo_mixture_projection(fit, grid, exp)
    mu <- posterior_epred(fit, grid)
    expect_equal(unname(p$mean[1, ]), unname(colMeans(exp(mu[1:40, , drop = FALSE]))))
    expect_gt(p$mean[1, 1], exp(mean(mu[1:40, 1])))
    expect_error(combo_mixture_projection(fit, grid, function(x) as.vector(x)), "shape-preserving")
    expect_error(combo_mixture_projection(fit, grid, function(x) x * Inf), "finite")
    expect_error(combo_mixture_projection(fit, grid, function(x) x - mean(x)), "elementwise")
    expect_error(combo_mixture_projection(fit, grid, function(x) x + runif(1)), "deterministic")
})

test_that("stable stacking shares indistinguishable columns and favors prediction", {
    lpd <- cbind(good = rep(-1000, 20), bad = rep(-1005, 20), copy = rep(-1000, 20))
    weights <- combo_mixture_weights(lpd)
    expect_equal(sum(weights), 1, tolerance = 1e-12)
    expect_equal(weights[1], weights[3], ignore_attr = TRUE)
    expect_gt(weights[1], .499)
    expect_lt(weights[2], .002)
    expect_equal(unname(combo_mixture_weights(matrix(1, 3, 1))), 1)
    expect_equal(combo_mixture_lme(matrix(c(-1000, -1001, -Inf, -Inf), 2)),
                 c(-1000 + log((1 + exp(-1)) / 2), -Inf))
    expect_error(combo_mixture_weights(matrix(NA_real_, 3, 2)), "finite")
    local_mocked_bindings(combo_mixture_optimize = function(...) cli::cli_abort("optimizer failed"))
    expect_error(combo_mixture_weights(lpd), "optimizer failed")
})

test_that("representative and pseudo-draw allocations are balanced and reproducible", {
    ids <- rep(1:7, each = 800)
    selected <- combo_mixture_balanced(ids, seq_along(ids), 800)
    expect_length(unique(selected), 800)
    expect_lte(diff(range(table(ids[selected]))), 1)
    final <- combo_mixture_balanced(ids, selected, 117)
    expect_lte(diff(range(table(ids[final]))), 1)
    expect_identical(combo_mixture_allocation(c(.5, .25, .25), 7), c(3L, 2L, 2L))
    expect_identical(combo_mixture_allocation(c(.9999, .0001), 40), c(40L, 0L))
})

test_that("PSIS preserves observed rows, likelihood centering, and chain efficiency", {
    skip_if_not_installed("loo")
    fit <- mixture_test_fit(shifts = c(0, 0, 2))
    membership <- data.frame(chain = 1:3, basin = c("B1", "B1", "B2"))
    seen <- list()
    original <- combo_mixture_loo_chunk
    local_mocked_bindings(combo_mixture_loo_chunk = function(ll, chain_id) {
        seen[[length(seen) + 1L]] <<- chain_id
        original(ll, chain_id)
    })
    result <- suppressWarnings(combo_mixture_psis(fit, membership, chunk_size = 2))
    expect_identical(as.integer(rownames(result$lpd)), fit$compiled$observed_rows)
    expect_equal(unique(result$pointwise$row_id), fit$compiled$observed_rows)
    expect_equal(seen[[1]], rep(1:2, each = 40))
    expect_equal(seen[[4]], rep(3, 40))
    ll <- log_lik(fit)
    first <- suppressWarnings(original(ll, fit$chain_id))
    second <- suppressWarnings(original(ll - 1000, fit$chain_id))
    expect_equal(first$r_eff, second$r_eff, tolerance = 1e-9)
    expect_equal(first$lpd - 1000, second$lpd, tolerance = 1e-8)
    constant <- suppressWarnings(original(matrix(-1, 40, 2), rep(1, 40)))
    expect_equal(constant$r_eff, c(1, 1))
    expect_error(original(matrix(Inf, 40, 2), rep(1, 40)), "Non-finite")
})

test_that("mixture predictions agree with their manually selected snapshots", {
    skip_if_not_installed("loo")
    fit <- mixture_test_fit(shifts = c(0, 0, 1))
    set.seed(184)
    seed <- .Random.seed
    x <- suppressWarnings(predictive_mixture(fit))
    expect_identical(.Random.seed, seed)
    y <- suppressWarnings(predictive_mixture(fit, ndraws = 25))
    expect_equal(x$weights$fitted, y$weights$fitted, tolerance = 1e-12)
    expect_false(inherits(x, "combo_fit"))
    expect_equal(sum(x$weights$count), 40)
    index <- match(paste(x$source_draws$chain, x$source_draws$draw), paste(fit$chain_id, fit$draw_id))
    manual <- combo_mixture_subset(fit, index)
    expect_equal(posterior_epred(x), posterior_epred(manual))
    expect_equal(posterior_linpred(x), posterior_linpred(manual))
    expect_equal(predict(x, type = "mean"), colMeans(posterior_epred(manual)))
    expect_equal(log_lik(x), log_lik(manual))
    expect_identical(attr(log_lik(x), "row_ids"), fit$compiled$observed_rows)
    set.seed(18)
    predicted <- posterior_predict(x)
    set.seed(18)
    expected <- posterior_predict(manual)
    expect_equal(predicted, expected)
    expect_equal(apply(predicted, 2, quantile, c(.05, .95)), apply(expected, 2, quantile, c(.05, .95)))
    expect_length(x$draws[[1]], 4)
    expect_named(summary(x), c("membership", "weights", "diagnostics", "reliable", "sensitivity", "provenance"))
    expect_output(print(x), "allocation_error")
    for (f in list(posterior_draws, parameter_map, posterior::as_draws,
                   posterior::as_draws_matrix, posterior::as_draws_array,
                   posterior::as_draws_df, posterior::as_draws_list,
                   posterior::as_draws_rvars, tidy, glance, stats::update, loo::loo)) {
        expect_error(f(x), "prediction-only")
    }
    expect_error(score_plate_pdbal(x, list()), "prediction-only")
    expect_error(predictive_mixture(fit, ndraws = 41), "ndraws")
    expect_error(predictive_mixture(fit, response_transform = function(z) z + runif(1)), "deterministic")
})

test_that("PSIS diagnostics flag adaptive and severe thresholds", {
    points <- data.frame(basin = "B1", row_id = 1:3, elpd_loo = -1,
                         pareto_k = c(.69, 1, NA), psis_n_eff = 20, r_eff = 1, draws = 100)
    expect_warning(x <- combo_mixture_psis_result(matrix(-1, 3, 1, dimnames = list(NULL, "B1")), points),
                   class = "combo_mixture_psis_warning")
    expect_false(x$reliable)
    expect_equal(x$diagnostics$count_ge_1, 1)
    expect_equal(x$diagnostics$flagged, 3)
})

test_that("prediction-preserving factor transformations leave projections unchanged", {
    set.seed(10)
    fit <- mixture_test_fit(model = combo_model(rank = 2, mean = fixed_mean(0)))
    snapshot <- gibbs_snapshot(init_gibbs_state(fit$compiled))
    fit$draws <- rep(list(snapshot), length(fit$draws))
    transformed <- fit
    transformed$draws <- lapply(transformed$draws, function(draw) {
        for (component in c("cell_factors", "treatment_main_factors", "treatment_interaction_factors")) {
            draw$components[[component]]$value <- draw$components[[component]]$value[, 2:1, drop = FALSE]
        }
        draw$components$cell_factors$value <- draw$components$cell_factors$value * 4
        draw$components$treatment_main_factors$value <- draw$components$treatment_main_factors$value / 4
        draw$components$treatment_interaction_factors$value <- -draw$components$treatment_interaction_factors$value / 2
        draw
    })
    expect_equal(posterior_epred(fit), posterior_epred(transformed), tolerance = 1e-12)
    grid <- combo_mixture_reference(fit)
    expect_equal(combo_mixture_projection(fit, grid, identity),
                 combo_mixture_projection(transformed, grid, identity), tolerance = 1e-12)
})

test_that("the public pipeline is invariant to chain storage order", {
    skip_if_not_installed("loo")
    fit <- mixture_test_fit(shifts = c(0, 0, 2))
    original <- suppressWarnings(predictive_mixture(fit, ndraws = 20))
    reordered <- combo_mixture_subset(fit, c(81:120, 40:1, 41:80))
    result <- suppressWarnings(predictive_mixture(reordered, ndraws = 20))
    expect_identical(original$membership, result$membership)
    expect_identical(original$source_draws, result$source_draws)
    expect_equal(original$weights, result$weights, tolerance = 1e-12)
    expect_equal(original$psis$pointwise, result$psis$pointwise, tolerance = 1e-12)
    expect_error(predictive_mixture(fit, ndraws = 1), "largest basin")
    expect_error(posterior::summarise_draws(original), "prediction-only")
})

test_that("mean basins retain separate predictive-spread diagnostics", {
    fit <- mixture_test_fit(chains = 2)
    fit$draws <- lapply(seq_along(fit$draws), function(i) {
        draw <- fit$draws[[i]]
        draw$precision <- if (fit$chain_id[i] == 1) 1 else .01
        draw
    })
    projection <- combo_mixture_projection(fit, combo_mixture_reference(fit), identity)
    expect_equal(projection$mean[1, ], projection$mean[2, ], ignore_attr = TRUE)
    expect_gt(min(projection$predictive_variance[2, ] - projection$predictive_variance[1, ]), 98)
    expect_identical(combo_mixture_cluster(projection$mean, projection$block_mean)$basin, c("B1", "B1"))
})

test_that("prediction mapping retains formula and entity metadata", {
    skip_if_not_installed("loo")
    data <- combo_test_data()
    data$x <- seq_len(nrow(data))
    fit <- mixture_test_fit(model = combo_model(rank = 1L, mean = formula_mean(~ x),
        cell_offset = NULL, treatment_offset = NULL, cell_factors = NULL,
        treatment_main_factors = NULL, treatment_interaction_factors = NULL), data = data)
    fit$draws <- lapply(fit$draws, function(draw) { draw$beta <- c(x = .2); draw })
    mixture <- suppressWarnings(predictive_mixture(fit))
    newdata <- data[c(1, 3, 5), ]
    newdata$x <- c(1.5, 7, 10)
    ids <- match(paste(mixture$source_draws$chain, mixture$source_draws$draw), paste(fit$chain_id, fit$draw_id))
    manual <- combo_mixture_subset(fit, ids)
    expect_equal(posterior_epred(mixture, newdata), posterior_epred(manual, newdata))
    expect_equal(log_lik(mixture, newdata), log_lik(manual, newdata))
    newdata$cell[1] <- "unseen"
    expect_error(posterior_epred(mixture, newdata), "absent")
})

test_that("all-constant likelihood and invalid efficiency have explicit behavior", {
    skip_if_not_installed("loo")
    ll <- matrix(-1, 40, 2)
    value <- combo_mixture_loo_chunk(ll, rep(1:2, each = 20))
    expect_equal(value$lpd, c(-1, -1), tolerance = 1e-12)
    expect_equal(value$r_eff, c(1, 1))
    expect_true(length(value$warnings) >= 1L)
    local_mocked_bindings(relative_eff = function(...) c(NA_real_, NA_real_), .package = "loo")
    ll <- ll + rep(sin(1:40), 2)
    expect_error(combo_mixture_loo_chunk(ll, rep(1:2, each = 20)), "Undefined PSIS relative")
})


test_that("stacking recovers small but useful mass in an imbalanced problem", {
    lpd <- cbind(core = c(rep(0, 999), -1000), alternative = c(rep(-1000, 999), 0))
    expect_equal(unname(combo_mixture_weights(lpd)), c(.999, .001), tolerance = 1e-6)
    expect_equal(unname(combo_mixture_weights(lpd - 10000)), c(.999, .001), tolerance = 1e-6)
    expect_equal(unname(combo_mixture_weights(rbind(c(0, -1000), c(-1000, 0)))), c(.5, .5), tolerance = 1e-6)
})
