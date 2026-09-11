test_that("model and component specifications form the public fit contract", {
    model <- combo_test_model()

    expect_s3_class(model, "combo_model")
    expect_identical(combo_model()$rank, 12L)
    expect_identical(combo_model(rank = 2L)$rank, 2L)
    expect_s3_class(model$components$cell_offset, "combo_gaussian_component")
    expect_s3_class(model$components$cell_offset, "combo_component")
    expect_s3_class(iid(), "structural_prior")
    expect_s3_class(fixed_scale(), "param_shrinkage")
    expect_s3_class(gamma_precision(), "param_shrinkage")
    expect_s3_class(global_half_cauchy(), "param_shrinkage")
    expect_s3_class(local_half_cauchy(), "param_shrinkage")
    expect_s3_class(horseshoe(), "param_shrinkage")
    expect_s3_class(multiplicative_gamma(), "param_shrinkage")
    expect_s3_class(empirical_mean(), "combo_mean")
    expect_s3_class(fixed_mean(0), "combo_mean")
    expect_s3_class(formula_mean(), "combo_mean")
    expect_s3_class(categorical(), "combo_dose")
    expect_s3_class(nested(), "combo_dose")
    expect_s3_class(gaussian_response(), "combo_family")
    expect_true("gaussian_response" %in% getNamespaceExports("batchieR"))
    expect_false("gaussian" %in% getNamespaceExports("batchieR"))
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
    map <- parameter_map(fit, include = "all")
    expect_true(any(map$component == "cell_factors"))
    expect_true(any(map$component == "treatment_main_factors"))
    expect_true(any(map$component == "treatment_interaction_factors"))
})

test_that("component configuration contains only structure and shrinkage", {
    expect_error(
        combo_gaussian_component(
            mean = ~ 0 + feature,
            mean_shrinkage = fixed_scale(),
            shrinkage = fixed_scale()
        ),
        "unused argument"
    )

    component <- combo_gaussian_component(shrinkage = fixed_scale())
    expect_named(component, c("structure", "shrinkage"))
})

test_that("fit_combo validates parallel and progress controls", {
    model <- combo_test_model()
    data <- combo_test_data()

    expect_error(
        fit_combo(model, data, parallel_chains = 0),
        "parallel_chains must be"
    )
    expect_error(
        fit_combo(model, data, parallel_chains = 1.5),
        "parallel_chains must be"
    )
    expect_error(fit_combo(model, data, refresh = -1), "refresh must be")
    expect_error(fit_combo(model, data, refresh = 1.5), "refresh must be")
})

test_that("chain RNG is deterministic across worker counts", {
    serial <- combo_test_fit(parallel_chains = 1L)
    parallel <- combo_test_fit(parallel_chains = 2L)
    capped <- combo_test_fit(parallel_chains = 10L)

    expect_identical(parallel$draws, serial$draws)
    expect_identical(capped$draws, serial$draws)
    expect_identical(parallel$chain_id, serial$chain_id)
    expect_identical(parallel$draw_id, serial$draw_id)
    expect_identical(parallel$iteration, serial$iteration)
    expect_identical(parallel$sampling$parallel_chains, 2L)
    expect_identical(capped$sampling$parallel_chains, 2L)
    expect_false(identical(
        serial$draws[serial$chain_id == 1L],
        serial$draws[serial$chain_id == 2L]
    ))
})

test_that("chain progress reports phases at the requested cadence", {
    events <- list()
    report <- batchieR:::combo_chain_reporter(
        refresh = 2L,
        iter_warmup = 2L,
        iter_sampling = 4L,
        chain = 1L,
        update = function(info) events[[length(events) + 1L]] <<- info
    )

    report(0L, "warmup")
    report(1L, "warmup")
    report(2L, "warmup")
    report(2L, "sampling")
    report(3L, "sampling")
    report(4L, "sampling")
    report(5L, "sampling")
    report(6L, "sampling")

    expect_identical(
        vapply(events, function(event) event[["current"]], integer(1)),
        c(0L, 2L, 2L, 4L, 6L)
    )
    expect_identical(
        vapply(events, function(event) event[["phase"]], character(1)),
        c("warmup", "warmup", "sampling", "sampling", "sampling")
    )
    expect_identical(
        vapply(events, function(event) event[["chain"]], integer(1)),
        rep(1L, 5L)
    )
    expect_identical(
        vapply(events, function(event) event[["amount"]], integer(1)),
        c(0L, 2L, 0L, 2L, 2L)
    )
})

test_that("fit_combo progress can be shown or disabled", {
    expect_silent(combo_test_fit(chains = 1L, refresh = 0L))
    expect_no_error(invisible(combo_test_fit(chains = 1L, refresh = 1L)))
    expect_no_error(invisible(combo_test_fit(parallel_chains = 2L, refresh = 1L)))
    expect_identical(future::nbrOfWorkers(), 1L)
})

test_that("progress keeps separate chain updates without changing draws", {
    local_mocked_bindings(combo_progress_terminal = function() FALSE)
    quiet <- combo_test_fit(refresh = 0L)
    for (workers in 1:2) {
        output <- capture.output(
            visible <- combo_test_fit(parallel_chains = workers, refresh = 1L),
            type = "message"
        )
        for (chain in 1:2) {
            lines <- output[startsWith(output, paste0("Chain ", chain, " "))]
            expected <- sprintf(
                "Chain %d Iteration: %d / 6 [%3.0f%%] (%s)",
                chain, c(0:2, 2:6), floor(100 * c(0:2, 2:6) / 6),
                c(rep("warmup", 3), rep("sampling", 5))
            )
            expect_identical(lines, expected)
        }
        expect_false(any(grepl("\r", output)))
        expect_match(tail(output, 3L)[[1]], "All 2 chains finished successfully.", fixed = TRUE)
        expect_match(tail(output, 2L)[[1]], "^Mean chain execution time: [0-9]+\\.[0-9] seconds\\.$")
        expect_match(tail(output, 1L), "^Total execution time: [0-9]+\\.[0-9] seconds\\.$")
        expect_identical(visible$draws, quiet$draws)
        expect_identical(visible$chain_id, quiet$chain_id)
        expect_identical(visible$iteration, quiet$iteration)
    }
    expect_identical(future::nbrOfWorkers(), 1L)
})

test_that("terminal progress keeps each live chain on its own row", {
    local_mocked_bindings(combo_progress_terminal = function() TRUE)
    quiet <- combo_test_fit(refresh = 0L)
    for (workers in 1:2) {
        output <- capture.output(
            visible <- combo_test_fit(parallel_chains = workers, refresh = 1L),
            type = "message"
        )
        expect_match(tail(output, 3L)[[1]], "All 2 chains finished successfully.", fixed = TRUE)
        # Each repaint moves up exactly two rows and replaces those same rows.
        bar_output <- head(output, -3L)
        frames <- strsplit(paste(bar_output, collapse = "\n"), "\033[2A", fixed = TRUE)[[1]]
        frames <- lapply(frames, function(frame) {
            rows <- strsplit(frame, "\n", fixed = TRUE)[[1]]
            rows <- rows[nzchar(rows)]
            sub("\r\033[2K", "", rows, fixed = TRUE)
        })
        expect_gt(length(frames), 2L)
        for (rows in frames) {
            expect_length(rows, 2L)
            expect_true(startsWith(rows[[1]], "Chain 1:"))
            expect_true(startsWith(rows[[2]], "Chain 2:"))
        }
        expect_true(all(grepl("100%$", tail(frames, 1L)[[1]])))
        # Completion of one chain must not make the other chain appear complete.
        expect_true(any(vapply(frames, function(rows) {
            sum(grepl("100%$", rows)) == 1L
        }, logical(1))))
        expect_identical(visible$draws, quiet$draws)
    }
    expect_silent(combo_test_fit(refresh = 0L))
})

test_that("interrupted terminal progress retains actual chain counts", {
    local_mocked_bindings(combo_progress_terminal = function() TRUE)
    output <- capture.output(
        expect_error(
            progressr::with_progress({
                progress <- progressr::progressor(steps = 8L)
                progress(amount = 1L, chain = 1L, current = 1L, phase = "sampling")
                stop("interrupted test fit")
            }, enable = TRUE, delay_conditions = c("message", "warning"),
            handlers = combo_progress_handler(2L, 4L, TRUE)),
            "interrupted test fit"
        ),
        type = "message"
    )
    expect_match(paste(output, collapse = "\n"), "25%")
    expect_false(any(grepl("100%", output)))
})

test_that("parallel chain failures identify the chain", {
    expect_error(
        batchieR:::combo_run_chains(
            compiled = list(),
            seed = 42L,
            chains = 2L,
            iter_warmup = 1L,
            iter_sampling = 1L,
            thin = 1L,
            parallel_chains = 2L,
            refresh = 0L
        ),
        "Chain 1 failed"
    )
})
