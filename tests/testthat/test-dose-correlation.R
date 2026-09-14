dose_test_treatments <- function() {
    data.frame(
        key = c("X:10", "Y:2", "X:0.1", "Z:1", "Y:20", "X:1"),
        drug = c("X", "Y", "X", "Z", "Y", "X"),
        dose = c("10", "2", "0.1", "1", "20", "1")
    )
}

dose_test_structures <- function() {
    drugs <- c("X", "Y", "Z")
    Q <- Matrix::Matrix(
        matrix(c(2, -0.5, 0, -0.5, 2, -0.5, 0, -0.5, 2), 3),
        sparse = TRUE, dimnames = list(drugs, drugs)
    )
    hierarchy <- data.frame(
        node = c("root", drugs), parent = c(NA, "root", "root", "root"),
        edge = c(1, 2, 3, 4)
    )
    list(iid(), precision(Q), gmrf(Q, ridge = 0.2), tree(node, parent, hierarchy, edge))
}

dose_test_covariance <- function(treatments, kernel, length_scale = 1) {
    distance <- abs(outer(log10(as.numeric(treatments$dose)),
                          log10(as.numeric(treatments$dose)), "-")) / length_scale
    K <- if (kernel == "exponential") exp(-distance) else
        (1 + sqrt(3) * distance) * exp(-sqrt(3) * distance)
    K[outer(treatments$drug, treatments$drug, "!=")] <- 0
    K
}

dose_test_model <- function(correlation, structure = NULL) {
    component <- combo_gaussian_component(shrinkage = gamma_precision(), structure = structure)
    combo_model(
        rank = 2L, mean = fixed_mean(0),
        treatment_offset = component, treatment_main_factors = component,
        treatment_interaction_factors = component,
        dose = nested(precision = 2, correlation = correlation)
    )
}

test_that("dose kernel constructors and model methods preserve their contract", {
    expect_s3_class(dose_kernel(), "combo_dose_kernel")
    expect_identical(dose_kernel()$kernel, "exponential")
    expect_identical(dose_kernel()$transform, "log10")
    expect_identical(dose_kernel()$length_scale, 1)
    expect_identical(nested(), nested(correlation = NULL))
    expect_error(nested(correlation = iid()), "correlation must be")
    expect_error(dose_kernel("unknown"), "kernel must be")
    expect_error(dose_kernel(transform = "sqrt"), "transform must be")
    for (value in list(0, -1, Inf, NA_real_, c(1, 2), "one")) {
        expect_error(dose_kernel(length_scale = value), "length_scale")
    }
    correlation <- dose_kernel("matern32", length_scale = 0.5)
    original <- dose_test_model(NULL)
    revised <- update(original, dose = nested(2, correlation))
    expect_identical(revised, dose_test_model(correlation))
    expect_null(original$dose$correlation)
    printed <- paste(capture.output(print(revised)), collapse = "\n")
    expect_match(printed, "matern32")
    expect_match(printed, "log10, length scale 0.5", fixed = TRUE)
    expect_false(any(grepl("correlation:", capture.output(print(original)))))
    expect_error(combo_model(dose = nested(correlation = correlation)), "entity-local")
})

test_that("kernel precisions match analytic correlations on irregular grids", {
    doses <- c(100, 0.1, 2, 0.4)
    for (kernel in c("exponential", "matern32")) {
        correlation <- dose_kernel(kernel, length_scale = 0.7)
        Q <- dose_kernel_precision(doses, correlation, "X")
        distance <- abs(outer(log10(doses), log10(doses), "-")) / 0.7
        K <- if (kernel == "exponential") exp(-distance) else
            (1 + sqrt(3) * distance) * exp(-sqrt(3) * distance)
        expect_equal(as.matrix(Q), solve(K), tolerance = 1e-10)
        expect_equal(as.matrix(solve(Q)), K, tolerance = 1e-10)
        expect_equal(as.matrix(dose_kernel_precision("4", correlation, "X")), matrix(1))
        identity <- dose_kernel(kernel, "identity", 0.7)
        expect_equal(Q, dose_kernel_precision(log10(doses), identity, "X"), tolerance = 1e-12)
    }
    Q <- dose_kernel_precision(sort(doses), dose_kernel(), "X")
    expect_true(all(as.matrix(Q)[abs(row(Q) - col(Q)) > 1] == 0))
})

test_that("correlated nesting preserves parent covariance and tree scaling", {
    treatments <- dose_test_treatments()
    compounds <- unique(treatments$drug)
    for (spec in dose_test_structures()) {
        base <- construct_structure(spec, compounds, "Compound")
        parent <- base$modeled_index[match(treatments$drug, compounds)]
        parent_covariance <- as.matrix(solve(base$Q))[parent, parent]
        legacy <- compile_treatment_structure(spec, nested(2), treatments)
        for (kernel in c("exponential", "matern32")) {
            compiled <- compile_treatment_structure(spec, nested(2, dose_kernel(kernel)), treatments)
            expected <- parent_covariance + dose_test_covariance(treatments, kernel) / 2
            expected <- expected / outer(legacy$modeled_scale, legacy$modeled_scale)
            index <- compiled$modeled_index
            actual <- as.matrix(solve(compiled$Q))[index, index] /
                outer(compiled$modeled_scale, compiled$modeled_scale)
            expect_equal(unname(actual), unname(expected), tolerance = 1e-10)
            expect_identical(compiled$nodes, legacy$nodes)
            expect_identical(compiled$modeled_index, legacy$modeled_index)
            expect_equal(compiled$modeled_scale, legacy$modeled_scale)
            expect_identical(compiled$entity_names, treatments$key)
            expect_false(compiled$iid)
            independent <- compile_treatment_structure(
                spec, nested(2, dose_kernel(kernel, length_scale = 1e-10)), treatments
            )
            expect_equal(independent$Q, legacy$Q, tolerance = 1e-12)
            expect_equal(independent$modeled_scale, legacy$modeled_scale)
        }
    }
})

test_that("correlated nesting still validates supplied compound precisions", {
    treatments <- dose_test_treatments()
    names <- unique(treatments$drug)
    Q <- Matrix::Diagonal(3)
    dimnames(Q) <- list(names, names)
    dose <- nested(correlation = dose_kernel())
    expect_error(compile_treatment_structure(precision(as.matrix(Q)), dose, treatments), "sparse")
    asymmetric <- methods::as(Q, "generalMatrix")
    asymmetric[1, 2] <- 0.1
    expect_error(compile_treatment_structure(precision(asymmetric), dose, treatments), "symmetric")
    Q[1, 1] <- -1
    expect_error(compile_treatment_structure(precision(Q), dose, treatments), "positive definite")
})

test_that("internal sorting does not change treatment mapping or covariance", {
    treatments <- dose_test_treatments()
    perm <- c(6, 4, 5, 3, 1, 2)
    for (kernel in c("exponential", "matern32")) {
        dose <- nested(3, dose_kernel(kernel))
        a <- compile_treatment_structure(dose_test_structures()[[4]], dose, treatments)
        b <- compile_treatment_structure(dose_test_structures()[[4]], dose, treatments[perm, ])
        covariance <- function(x) {
            as.matrix(solve(x$Q))[x$modeled_index, x$modeled_index] /
                outer(x$modeled_scale, x$modeled_scale)
        }
        expect_identical(b$entity_names, treatments$key[perm])
        expect_equal(unname(covariance(b)), unname(covariance(a)[perm, perm]), tolerance = 1e-10)
    }
})

test_that("dose coordinate errors identify the compound and kernel settings", {
    correlation <- dose_kernel()
    for (doses in list(c("low", "high"), c(NA, 1), c(Inf, 1))) {
        expect_error(dose_kernel_precision(doses, correlation, "X"), "compound X.*exponential.*finite numeric")
    }
    for (doses in list(c(0, 1), c(-1, 1))) {
        expect_error(dose_kernel_precision(doses, correlation, "X"), "strictly positive")
    }
    expect_error(dose_kernel_precision(c("1", "1.0"), correlation, "X"), "duplicate coordinates")
    # Distinct doubles can round to the same log10 coordinate.
    expect_error(dose_kernel_precision(c(1e100, 1e100 * (1 + 2e-16)), correlation, "X"),
                 "duplicate coordinates")
    for (kernel in c("exponential", "matern32")) {
        expect_no_error(dose_kernel_precision(c(-2, 0, 1), dose_kernel(kernel, "identity"), "X"))
        tiny <- dose_kernel(kernel, "identity", length_scale = 1e-300)
        expect_equal(as.matrix(dose_kernel_precision(c(-1e100, 1e100), tiny, "X")), diag(2))
        huge <- dose_kernel(kernel, length_scale = 1e300)
        expect_error(dose_kernel_precision(c(1, 10, 100), huge, "X"), "compound X.*length scale")
    }
    close <- dose_kernel_precision(c(0, 1e-7, 2e-7), dose_kernel(transform = "identity"), "X")
    expect_true(all(is.finite(sparse_values(close))))
})

test_that("old nested specifications retain compiled priors and seeded fits", {
    legacy <- structure(list(type = "nested", precision = 2), class = "combo_dose")
    treatments <- dose_test_treatments()
    treatments$dose <- letters[seq_len(nrow(treatments))]
    for (spec in dose_test_structures()) {
        expect_identical(compile_treatment_structure(spec, legacy, treatments),
                         compile_treatment_structure(spec, nested(2), treatments))
    }
    modern_model <- dose_test_model(NULL)
    legacy_model <- modern_model
    legacy_model$dose <- legacy
    modern <- combo_test_fit(model = modern_model, chains = 1L)
    old <- combo_test_fit(model = legacy_model, chains = 1L)
    expect_identical(modern$draws, old$draws)
})

test_that("correlated raw prior draws reproduce analytic normalized covariance", {
    treatments <- dose_test_treatments()[c(1, 3, 6), ]
    for (kernel in c("exponential", "matern32")) {
        spec <- dose_test_structures()[[4]]
        compiled <- compile_treatment_structure(spec, nested(2, dose_kernel(kernel)), treatments)
        set.seed(516)
        draws <- prior_draw(compiled, list(global = rep(2, 4000), local = NULL), 4000)$value
        # X has parent path variance 1 + 2; each terminal increment has variance 1/2.
        expected <- (matrix(3, 3, 3) + dose_test_covariance(treatments, kernel) / 2) / (3.5 * 2)
        expect_equal(unname(rowMeans(draws)), rep(0, 3), tolerance = 0.04)
        expect_equal(unname(stats::cov(t(draws))), expected, tolerance = 0.04)
    }
})

test_that("a Gaussian treatment update matches its analytic conditional posterior", {
    data <- data.frame(
        cell = "A", drug_1 = "X", dose_1 = c(1, 10),
        drug_2 = NA_character_, dose_2 = NA_real_, response = c(0.7, -0.1)
    )
    for (kernel in c("exponential", "matern32")) {
        model <- combo_model(
            mean = fixed_mean(0), cell_offset = NULL, cell_factors = NULL,
            treatment_main_factors = NULL, treatment_interaction_factors = NULL,
            treatment_offset = combo_gaussian_component(shrinkage = fixed_scale(2)),
            dose = nested(2, dose_kernel(kernel))
        )
        state <- init_gibbs_state(compile_combo_model(model, data))
        component <- state$components$treatment_offset
        component$raw[, 1] <- c(0.2, 0.3, -0.4)
        index <- component$compiled$structure$modeled_index
        component$values[, 1] <- component$raw[index, 1]
        # Hold parents fixed to test a scalar conditional update independently.
        component$compiled$structure$latent_index <- integer()
        state$components$treatment_offset <- component
        state$precision <- 3
        state$fitted_mean <- gibbs_reconstruct_mean(state)
        Q <- as.matrix(component$compiled$structure$Q) * 2
        node <- index[1]
        precision <- Q[node, node] + 3
        linear <- -sum(Q[node, -node] * component$raw[-node, 1]) + 3 * data$response[1]
        set.seed(799)
        values <- replicate(2500, gibbs_update_treatment_offset(state)$components$treatment_offset$values[1, 1])
        expect_equal(mean(values), linear / precision, tolerance = 0.025)
        expect_equal(stats::var(values), 1 / precision, tolerance = 0.02)
    }
})

test_that("both kernels fit all tensor components and predict registered missing doses", {
    data <- combo_test_data()
    extra <- data[c(1, 3), ]
    extra$dose_1 <- 10
    extra$response <- NA_real_
    data <- rbind(data, extra, data[1, ])
    # A control has neither treatment; it must not become a dose coordinate.
    control <- data[1, ]
    control$drug_1 <- NA_character_
    control$dose_1 <- NA_real_
    data <- rbind(data, control)
    swapped <- data
    swapped[c("drug_1", "dose_1", "drug_2", "dose_2")] <-
        data[c("drug_2", "dose_2", "drug_1", "dose_1")]
    for (kernel in c("exponential", "matern32")) {
        model <- dose_test_model(dose_kernel(kernel))
        fit <- function(input, specification = model) fit_combo(
            specification, input, chains = 1, iter_warmup = 3, iter_sampling = 4,
            seed = 63, refresh = 0
        )
        fitted <- fit(data)
        expect_identical(fitted$draws, fit(data)$draws)
        prediction <- posterior_epred(fitted)
        expect_identical(dim(prediction), c(4L, nrow(data)))
        expect_identical(colnames(prediction), as.character(seq_len(nrow(data))))
        expect_true(all(is.finite(prediction)))
        expect_equal(posterior_epred(fitted, swapped), prediction)
        independent <- fit(data, dose_test_model(NULL))
        for (name in names(fitted$draws[[1]]$components)) {
            a <- fitted$draws[[1]]$components[[name]]
            b <- independent$draws[[1]]$components[[name]]
            expect_identical(dim(a$value), dim(b$value))
            expect_identical(dimnames(a$value), dimnames(b$value))
            expect_identical(dimnames(a$raw), dimnames(b$raw))
            expect_true(all(is.finite(a$value)))
        }
        expect_identical(dimnames(prior_predict(model, data, draws = 4, seed = 3)), dimnames(prediction))
        unseen <- extra[1, ]
        unseen$dose_1 <- 100
        expect_error(posterior_epred(fitted, unseen), "absent from the fitted mappings")
        malformed <- data
        malformed$dose_1[nrow(data) - 2] <- "high"
        expect_error(fit(malformed), "finite numeric")
    }
})
