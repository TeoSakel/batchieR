# The experiment is intentionally absent from installed/source package builds.
sampler_benchmark_environment <- function() {
    directory <- testthat::test_path("..", "..", "experiments", "samplers")
    skip_if_not(file.exists(file.path(directory, "config.R")), "Development experiment is not in this package build")
    env <- new.env(parent = globalenv())
    for (file in c("config.R", "data.R", "sbc.R", "metrics.R", "report.R", "runner.R")) {
        sys.source(file.path(directory, file), envir = env)
    }
    env$root <- normalizePath(file.path(directory, "../.."))
    env
}

test_that("synthetic generation preserves RNG, truth and complete nested plate masks", {
    b <- sampler_benchmark_environment()
    design <- b$sb_design(2, 3, 2)
    expect_equal(nrow(design), 36)
    set.seed(751)
    rng <- .Random.seed
    first <- b$sb_generate(design, "gamma", 2L, 12L, 13L)
    expect_identical(.Random.seed, rng)
    expect_identical(first, b$sb_generate(design, "gamma", 2L, 12L, 13L))
    expect_true(all(first$masks$reveal_25 <= first$masks$reveal_75))
    for (mask in first$masks) {
        expect_true(all(vapply(split(mask, design$plate), function(x) length(unique(x)) == 1L, logical(1))))
        expect_true(all(mask[is.na(design$drug_2)]))
        expect_true(any(!mask))
    }
    observed <- first$masks$reveal_25
    fit_data <- b$sb_fit_data(first, "reveal_25")
    expect_setequal(names(fit_data), c(names(design), "response"))
    expect_true(all(is.na(fit_data$response[!observed])))
    altered <- first
    altered$response[!observed] <- 1e8
    altered$mean <- 1e8
    expect_identical(fit_data, b$sb_fit_data(altered, "reveal_25"))
    # Same seed reproduces the established generative API exactly.
    expect_equal(first$response, as.numeric(prior_predict(first$model, design, draws = 1L, seed = 12L)))
    terms <- b$sb_truth_terms(first$snapshot, first$indices)
    expect_equal(terms$mean, first$mean)
    expect_equal(unname(terms$interaction[is.na(design$drug_2)]), rep(0, sum(is.na(design$drug_2))))
    permuted <- first$snapshot
    for (name in c("cell_factors", "treatment_main_factors", "treatment_interaction_factors")) {
        permuted$components[[name]]$value <- permuted$components[[name]]$value[, 2:1, drop = FALSE]
    }
    expect_equal(b$sb_truth_terms(permuted, first$indices), terms, tolerance = 1e-12)
})

test_that("matched priors cover shrinkage families and realism uses shared dose vectors", {
    b <- sampler_benchmark_environment()
    design <- b$sb_design(2, 3, 3)
    for (scenario in head(b$sb_scenarios(), 3)) {
        d <- b$sb_generate(design, scenario, 4L, 77L, 3L)
        expect_identical(d$model, d$generating_model)
        expect_true(d$matched_prior)
        expect_true(all(is.finite(d$response)))
    }
    gamma <- b$sb_model("gamma", 4L)
    mg <- b$sb_model("multiplicative_gamma", 4L)
    hs <- b$sb_model("horseshoe", 4L)
    expect_identical(gamma$components$cell_factors$shrinkage$type, "gamma")
    expect_identical(mg$components$cell_factors$shrinkage$type, "multiplicative_gamma")
    expect_identical(hs$components$treatment_main_factors$shrinkage$type, "horseshoe")
    d <- b$sb_generate(design, "realism", 4L, 77L, 3L)
    expect_false(d$matched_prior)
    expect_identical(d$generating_model$rank, 2L)
    expect_identical(d$model$rank, 4L)
    compiled <- batchieR:::compile_combo_design(d$generating_model, design)
    treatments <- compiled$treatments
    for (drug in unique(treatments$drug)) {
        index <- which(treatments$drug == drug)
        index <- index[order(as.numeric(treatments$dose[index]))]
        factors <- d$snapshot$components$treatment_main_factors$value[index, , drop = FALSE]
        expect_equal(factors[3L, ], 3 * factors[1L, ], tolerance = 1e-12)
    }
})

test_that("Merck import strips retrospective values and keeps mixed plates complete", {
    b <- sampler_benchmark_environment()
    design <- b$sb_design(2, 3, 2)
    design$truth <- design$response <- 99
    path <- tempfile(fileext = ".rds")
    on.exit(unlink(path))
    saveRDS(list(training = head(design, 20), holdout = tail(design, 16)), path)
    imported <- b$sb_design(1, 2, 1, merck = path)
    expect_false(any(c("truth", "response") %in% names(imported)))
    expect_identical(imported, design[setdiff(names(design), c("truth", "response"))])
    imported$plate[13L] <- imported$plate[1L]
    masks <- b$sb_masks(imported, 7L)
    expect_true(masks$reveal_25[13L])
})

test_that("baseline adapter reproduces direct fits and likelihood uses observed rows", {
    b <- sampler_benchmark_environment()
    d <- b$sb_generate(b$sb_design(2, 3, 2), "gamma", 2L, 8L, 4L)
    data <- b$sb_fit_data(d, "reveal_25")
    settings <- list(chains = 2L, iter_warmup = 2L, iter_sampling = 8L, thin = 1L,
                     parallel_chains = 1L, refresh = 0L, seed = 45L)
    adapter <- b$sb_adapters()$gibbs
    fit <- adapter$fit(d$model, data, settings)
    direct <- do.call(fit_combo, c(list(model = d$model, data = data), settings))
    expect_identical(fit, direct)
    value <- b$sb_evaluate(fit, d, "reveal_25", 1, 15L, 16L)
    expected <- rowSums(log_lik(fit))
    expect_equal(as.numeric(value$draws[, , "observed_log_likelihood"]), expected, tolerance = 1e-10)
    observed <- d$masks$reveal_25
    expect_equal(unname(value$truth["observed_log_likelihood"]),
                 sum(dnorm(d$response[observed], d$mean[observed], 1 / sqrt(d$snapshot$precision), log = TRUE)))
    expect_identical(dim(value$draws)[1:2], c(8L, 2L))
    expect_true(all(value$coverage$coverage >= 0 & value$coverage$coverage <= 1))
    expect_true(all(value$coverage$mean_width >= 0))
    expect_true(all(value$ranks$status == "insufficient_draws"))
    changed <- d
    changed$response[!observed] <- changed$response[!observed] + 100
    second <- b$sb_evaluate(fit, changed, "reveal_25", 1, 15L, 16L)
    expect_identical(value$draws, second$draws)
    expect_identical(value$truth["observed_log_likelihood"], second$truth["observed_log_likelihood"])
    expect_gt(second$metrics$response_holdout_rmse, value$metrics$response_holdout_rmse)
    # Exercise one-column interval bounds and reordered chain bookkeeping.
    expect_equal(b$sb_interval(matrix(c(0, 1, 2), ncol = 1), 1, 0.5), c(covered = 1, width = 1, n = 1))
    shuffled <- seq(length(fit$draws), 1L)
    perm <- fit
    for (field in c("draws", "chain_id", "draw_id", "iteration")) perm[[field]] <- fit[[field]][shuffled]
    expect_equal(b$sb_draw_array(matrix(expected[shuffled], ncol = 1), perm),
                 b$sb_draw_array(matrix(expected, ncol = 1), fit))
})

test_that("a one-row final evaluation chunk preserves matrix dimensions", {
    b <- sampler_benchmark_environment()
    base <- b$sb_design(2, 3, 2)
    design <- base[rep(seq_len(nrow(base)), length.out = 257L), ]
    rownames(design) <- NULL
    d <- b$sb_generate(design, "gamma", 2L, 19L, 14L)
    settings <- list(chains = 1L, iter_warmup = 2L, iter_sampling = 4L, thin = 1L,
                     parallel_chains = 1L, refresh = 0L, seed = 34L)
    result <- b$sb_fit_job(d, "reveal_25", b$sb_adapters()$gibbs, settings, 35L, 36L)
    expect_identical(result$status, "ok")
    expect_identical(dim(result$draws)[1:2], c(4L, 1L))
    expect_true(all(result$coverage$n[result$coverage$target == "interaction"] == 257L))
    expect_true(all(result$coverage$coverage >= 0 & result$coverage$coverage <= 1))
})

test_that("SBC rank, envelope, eligibility and negative controls are sensitive", {
    b <- sampler_benchmark_environment()
    expect_identical(b$sb_rank(-1, 0:3), 0L)
    expect_identical(b$sb_rank(4, 0:3), 4L)
    expect_true(is.na(b$sb_rank(NA_real_, 0:3)))
    tied <- b$sb_with_seed(12, replicate(500, b$sb_rank(1, c(0, 1, 1, 2))))
    expect_setequal(unique(tied), 1:3)
    expect_equal(b$sb_ecdf(0:100), rep(0, 101))
    envelope <- b$sb_envelope(100, seed = 4L)
    expect_identical(envelope, b$sb_envelope(100, seed = 4L))
    expect_true(all(envelope$lower <= 0 & envelope$upper >= 0))
    expect_equal(tail(envelope$upper, 1), 0)
    checks <- b$sb_self_checks(232L)
    expect_true(checks$passed)
    expect_false(checks$summary$parameter_flag[checks$summary$adapter == "prior_only"])
    expect_true(checks$summary$likelihood_flag[checks$summary$adapter == "prior_only"])
    draws <- b$sb_with_seed(19, array(rnorm(4000), c(1000, 4, 1), dimnames = list(NULL, NULL, "x")))
    diagnostics <- b$sb_summary_draws(draws)
    ranks <- b$sb_ranks(draws, c(x = 0), diagnostics, 42L, TRUE)
    expect_true(ranks$status %in% c("eligible", "residual_autocorrelation"))
    diagnostics$rhat <- 1.5
    expect_identical(b$sb_ranks(draws, c(x = 0), diagnostics, 42L, TRUE)$status, "inadequate_convergence_or_ess")
    expect_identical(b$sb_ranks(draws, c(x = 0), diagnostics, 42L, FALSE)$status, "not_matched_prior")
})

test_that("runner resumes committed fits, refuses changed specs and renders saved results", {
    b <- sampler_benchmark_environment()
    skip_if_not_installed("quarto")
    skip_if(is.null(quarto::quarto_path()), "Quarto CLI unavailable")
    output <- tempfile("sampler-benchmark-")
    on.exit(unlink(output, recursive = TRUE), add = TRUE)
    config <- b$sb_config("smoke", scenarios = "gamma")
    expect_error(b$sb_config(overrides = list(iter_sampling = 0)), "iter_sampling")
    expect_error(b$sb_config(scenarios = "missing"), "scenarios")
    adapters <- b$sb_adapters()
    calls <- 0L
    adapters$gibbs$fit <- function(model, data, settings) {
        calls <<- calls + 1L
        do.call(fit_combo, c(list(model = model, data = data), settings))
    }
    b$sb_run(b$root, output, config, adapters)
    expect_identical(calls, 2L)
    html_path <- file.path(output, "report.html")
    expect_true(file.exists(html_path))
    html <- paste(readLines(html_path, warn = FALSE), collapse = "\n")
    for (title in c("Run and samplers", "Screen and generating scenarios", "Mixing and efficiency",
                    "Recovery and prediction", "Simulation-based calibration", "Self-checks, failures, and reproducibility")) {
        expect_match(html, title, fixed = TRUE)
    }
    expect_match(html, "data:image/png;base64", fixed = TRUE)
    expect_match(html, "insufficient_draws", fixed = TRUE)
    b$sb_run(b$root, output, config, adapters, resume = TRUE)
    expect_identical(calls, 2L)
    changed <- config
    changed$iter_sampling <- 9L
    expect_error(b$sb_run(b$root, output, changed, adapters, resume = TRUE), "Resume refused")
    # Delete one committed fit to represent interruption, then report without fitting.
    fits <- list.files(file.path(output, "fits"), full.names = TRUE)
    unlink(fits[1L])
    b$sb_finish_report(b$root, output)
    expect_identical(calls, 2L)
    expect_identical(readRDS(file.path(output, "report-data.rds"))$completion, "PARTIAL / INCOMPLETE")
    expect_match(paste(readLines(html_path, warn = FALSE), collapse = "\n"), "PARTIAL / INCOMPLETE", fixed = TRUE)
    b$sb_run(b$root, output, config, adapters, resume = TRUE)
    expect_identical(calls, 3L)
    # Exercise eligible-rank report rendering with a clearly synthetic fixture.
    record <- readRDS(fits[1L])
    record$ranks$status <- "eligible"
    record$ranks$rank <- 50L
    saveRDS(record, fits[1L])
    b$sb_finish_report(b$root, output)
    expect_identical(calls, 3L)
    expect_gt(nrow(readRDS(file.path(output, "report-data.rds"))$ecdfs), 0)
    for (path in fits) {
        record <- readRDS(path)
        record$status <- "failed"
        record$error <- "Deliberate report-only failure fixture"
        saveRDS(record, path)
    }
    b$sb_finish_report(b$root, output)
    failed <- readRDS(file.path(output, "report-data.rds"))
    expect_identical(failed$completion, "COMPLETED WITH FAILURES")
    expect_equal(nrow(failed$metrics), 0)
    expect_identical(calls, 3L)
})
