prior_test_design <- function() {
    data.frame(
        cell = c("A", "A", "B", "B"),
        drug_1 = c("X", "X", "Y", "X"), dose_1 = c(1, 2, 1, 1),
        drug_2 = c(NA, "X", NA, "Y"), dose_2 = c(NA, 2, NA, 1),
        stringsAsFactors = FALSE
    )
}

prior_mean_only_model <- function(mean = empirical_mean()) {
    combo_model(
        rank = 1L, mean = mean,
        cell_offset = NULL, cell_factors = NULL, treatment_offset = NULL,
        treatment_main_factors = NULL, treatment_interaction_factors = NULL,
        family = gaussian_response(
            precision = gamma_precision(shape = 2, rate = 3)
        )
    )
}

test_that("prior_predict returns reproducible draws-by-rows matrices", {
    design <- prior_test_design()
    model <- combo_test_model()
    first <- prior_predict(
        model, design, draws = 5L, seed = 104L, intercept = 0.25
    )
    second <- prior_predict(
        model, design, draws = 5L, seed = 104L, intercept = 0.25
    )

    expect_identical(first, second)
    expect_identical(dim(first), c(5L, 4L))
    expect_identical(colnames(first), as.character(seq_len(nrow(design))))
    expect_true(all(is.finite(first)))
})

test_that("prior_predict resolves empirical and fixed means", {
    design <- prior_test_design()
    design$response <- c(1, NA, 3, NA)
    empirical <- prior_mean_only_model()
    from_response <- prior_predict(
        empirical, design, draws = 3L, seed = 9L, type = "mean"
    )
    explicit <- prior_predict(
        empirical, transform(design, response = "ignored"),
        draws = 3L, seed = 9L, intercept = 2, type = "mean"
    )
    expect_identical(from_response, explicit)
    expect_equal(
        from_response,
        matrix(2, 3L, 4L,
            dimnames = list(NULL, as.character(seq_len(4L))))
    )

    no_response <- prior_test_design()
    expect_error(
        prior_predict(empirical, no_response),
        "requires intercept or at least one observed response"
    )
    expect_error(
        prior_predict(empirical, design, intercept = Inf),
        "one finite numeric value"
    )

    fixed <- prior_mean_only_model(fixed_mean(0.25))
    expect_equal(
        prior_predict(fixed, design, draws = 2L, seed = 3L, type = "mean"),
        matrix(0.25, 2L, 4L,
            dimnames = list(NULL, as.character(seq_len(4L))))
    )
    expect_error(
        prior_predict(fixed, design, intercept = 0.25),
        "must be NULL"
    )
})

test_that("formula prior prediction handles the flat intercept explicitly", {
    design <- transform(prior_test_design(), x = c(-1, 0, 1, 2))
    with_intercept <- prior_mean_only_model(
        formula_mean(~ x, beta_mean = c(x = 0.5), beta_precision = c(x = 4))
    )
    without_intercept <- prior_mean_only_model(
        formula_mean(~ 0 + x, beta_mean = c(x = 0.5), beta_precision = c(x = 4))
    )

    expect_error(prior_predict(with_intercept, design), "requires intercept")
    expect_no_error(
        prior_predict(with_intercept, design, intercept = 1, draws = 2L)
    )
    expect_error(
        prior_predict(without_intercept, design, intercept = 0),
        "intercept-free"
    )
    expect_no_error(prior_predict(without_intercept, design, draws = 2L))
})

test_that("formula prior slopes follow the configured Gaussian distribution", {
    design <- transform(prior_test_design()[1L, , drop = FALSE], x = 1)
    model <- prior_mean_only_model(formula_mean(
        ~ 0 + x,
        beta_mean = c(x = 1.5),
        beta_precision = c(x = 4)
    ))
    draws <- as.numeric(prior_predict(
        model, design, draws = 5000L, seed = 71L, type = "mean"
    ))

    expect_equal(mean(draws), 1.5, tolerance = 0.03)
    expect_equal(stats::var(draws), 0.25, tolerance = 0.025)
})

test_that("prior_predict samples observation precision and noise generatively", {
    design <- prior_test_design()[1:2, , drop = FALSE]
    model <- prior_mean_only_model(fixed_mean(1.5))
    set.seed(301L)
    expected <- matrix(NA_real_, nrow = 4L, ncol = nrow(design))
    for (draw in seq_len(4L)) {
        precision <- stats::rgamma(1L, shape = 2, rate = 3)
        expected[draw, ] <- stats::rnorm(
            nrow(design), mean = 1.5, sd = 1 / sqrt(precision)
        )
    }
    colnames(expected) <- c("1", "2")
    expect_identical(
        prior_predict(model, design, draws = 4L, seed = 301L), expected
    )
})

test_that("prior shrinkage samplers cover every supported family", {
    specifications <- list(
        fixed_scale(precision = 2), gamma_precision(shape = 2, rate = 3),
        global_half_cauchy(scale = 2), local_half_cauchy(scale = 2),
        horseshoe(global_scale = 2, local_scale = 3),
        multiplicative_gamma(shape = 2, rate = 3)
    )
    set.seed(18L)
    states <- lapply(
        specifications, prior_draw_shrinkage,
        n_entities = 4L, n_dimensions = 3L
    )
    for (state in states) {
        expect_length(state$global, 3L)
        expect_true(all(is.finite(state$global)))
        expect_true(all(state$global > 0))
        if (!is.null(state$local)) {
            expect_identical(dim(state$local), c(4L, 3L))
            expect_true(all(is.finite(state$local)))
            expect_true(all(state$local > 0))
        }
    }
    expect_equal(states[[1L]]$global, rep(2, 3L))
    expect_null(states[[3L]]$local)
})

test_that("prior_draw maps raw nodes to named entity values", {
    structure <- list(
        iid = TRUE, nodes = c("latent", "A", "B", "C"),
        modeled_index = c(2L, 3L, 4L), modeled_scale = c(2, 4, 5),
        entity_names = c("A", "B", "C")
    )
    shrinkage <- list(global = c(1, 1), local = NULL)
    set.seed(42L)
    draw <- prior_draw(structure, shrinkage, n_dim = 2L)

    expect_identical(dim(draw$value), c(3L, 2L))
    expect_identical(
        dimnames(draw$value),
        list(c("A", "B", "C"), c("dim_1", "dim_2"))
    )
    expect_equal(
        unname(draw$value),
        draw$raw[structure$modeled_index, , drop = FALSE] /
            structure$modeled_scale
    )
})

test_that("prior_predict supports structured component priors", {
    design <- prior_test_design()
    Q <- Matrix::Diagonal(2L, x = c(1, 2))
    dimnames(Q) <- list(c("A", "B"), c("A", "B"))
    operator <- Matrix::sparseMatrix(
        i = c(1, 1, 2, 2), j = c(1, 2, 1, 2),
        x = c(1, -1, -1, 1),
        dimnames = list(c("A", "B"), c("A", "B"))
    )
    hierarchy <- data.frame(
        node = c("root", "A", "B"), parent = c(NA, "root", "root"),
        stringsAsFactors = FALSE
    )
    structures <- list(
        iid(), precision(Q), gmrf(operator, ridge = 1),
        tree(node, parent, hierarchy)
    )
    for (component_structure in structures) {
        model <- combo_model(
            rank = 1L, mean = fixed_mean(0),
            cell_offset = combo_gaussian_component(
                shrinkage = fixed_scale(), structure = component_structure
            ),
            cell_factors = NULL, treatment_offset = NULL,
            treatment_main_factors = NULL,
            treatment_interaction_factors = NULL
        )
        predictions <- prior_predict(model, design, draws = 2L, seed = 7L)
        expect_identical(dim(predictions), c(2L, 4L))
        expect_true(all(is.finite(predictions)))
    }
})

test_that("prior prediction permits interaction self-combinations", {
    model <- combo_model(rank = 2L, mean = fixed_mean(0))
    predictions <- prior_predict(
        model, prior_test_design(), draws = 2L, seed = 8L
    )
    expect_identical(dim(predictions), c(2L, 4L))
})

test_that("prior_predict validates draw controls", {
    design <- prior_test_design()
    model <- prior_mean_only_model(fixed_mean(0))
    expect_error(prior_predict(list(), design), "combo_model")
    expect_error(prior_predict(model, design, draws = 0), "at least 1")
    expect_error(prior_predict(model, design, draws = 1.5), "at least 1")
    expect_error(prior_predict(model, design, seed = 0), "seed must be NULL")
    expect_error(prior_predict(model, design, type = "unknown"), "arg")
})
