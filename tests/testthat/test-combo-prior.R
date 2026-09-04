prior_test_design <- function() {
    data.frame(
        cell = c("A", "A", "B", "B"),
        drug_1 = c("X", "X", "Y", "X"),
        dose_1 = c(1, 2, 1, 1),
        drug_2 = c(NA, "X", NA, "Y"),
        dose_2 = c(NA, 2, NA, 1),
        stringsAsFactors = FALSE
    )
}

prior_intercept_only_model <- function(intercept = empirical_intercept()) {
    combo_model(
        rank = 1L,
        global_intercept = intercept,
        cell_offset = NULL,
        cell_factors = NULL,
        treatment_offset = NULL,
        treatment_main_factors = NULL,
        treatment_interaction_factors = NULL,
        family = gaussian_response(
            precision = gamma_precision(shape = 2, rate = 3)
        )
    )
}

test_that("prior_predict returns reproducible draws-by-rows matrices", {
    design <- prior_test_design()
    model <- combo_test_model()

    first <- prior_predict(
        model,
        design,
        draws = 5L,
        seed = 104L,
        intercept = 0.25
    )
    second <- prior_predict(
        model,
        design,
        draws = 5L,
        seed = 104L,
        intercept = 0.25
    )

    explicit_zero <- prior_predict(
        model,
        design,
        draws = 5L,
        seed = 104L,
        intercept = 0.25,
        beta_offset = list(cell_offset = 0, treatment_offset = 0)
    )

    omitted <- prior_predict(
        model,
        design,
        draws = 5L,
        seed = 104L,
        intercept = 0.25,
        beta_offset = list()
    )

    expect_identical(first, second)
    expect_identical(first, explicit_zero)
    expect_identical(first, omitted)
    expect_identical(dim(first), c(5L, 4L))
    expect_identical(colnames(first), as.character(seq_len(nrow(design))))
    expect_true(is.numeric(first))
    expect_true(all(is.finite(first)))
})

test_that("prior_predict can return latent conditional means", {
    design <- prior_test_design()
    model <- prior_intercept_only_model(fixed_intercept(0.25))

    means <- prior_predict(
        model, design, draws = 3L, seed = 9L, type = "mean"
    )

    expect_equal(
        means,
        matrix(
            0.25,
            nrow = 3L,
            ncol = nrow(design),
            dimnames = list(NULL, as.character(seq_len(nrow(design))))
        )
    )
    expect_error(prior_predict(model, design, type = "unknown"), "arg")
})

test_that("prior_predict resolves empirical and fixed intercepts", {
    design <- prior_test_design()
    empirical <- prior_intercept_only_model()
    design$response <- c(1, NA, 3, NA)

    from_response <- prior_predict(empirical, design, draws = 3L, seed = 9L)
    explicit <- prior_predict(
        empirical,
        transform(design, response = "ignored"),
        draws = 3L,
        seed = 9L,
        intercept = 2
    )
    expect_identical(from_response, explicit)

    no_response <- prior_test_design()
    expect_error(
        prior_predict(empirical, no_response),
        "requires intercept or at least one observed response"
    )
    no_response$response <- NA_real_
    expect_error(
        prior_predict(empirical, no_response),
        "requires intercept or at least one observed response"
    )
    expect_error(
        prior_predict(empirical, design, intercept = Inf),
        "one finite numeric value"
    )

    fixed <- prior_intercept_only_model(fixed_intercept(2))
    ignored <- transform(design, response = "not numeric")
    expect_identical(
        prior_predict(fixed, ignored, draws = 3L, seed = 9L),
        explicit
    )
    expect_error(
        prior_predict(fixed, design, intercept = 2),
        "must be NULL"
    )
})

test_that("prior_predict samples observation precision and noise generatively", {
    design <- prior_test_design()[1:2, , drop = FALSE]
    model <- prior_intercept_only_model(fixed_intercept(1.5))
    draws <- 4L

    set.seed(301L)
    expected <- matrix(NA_real_, nrow = draws, ncol = nrow(design))
    for (draw in seq_len(draws)) {
        precision <- stats::rgamma(1L, shape = 2, rate = 3)
        expected[draw, ] <- stats::rnorm(
            nrow(design),
            mean = 1.5,
            sd = 1 / sqrt(precision)
        )
    }
    colnames(expected) <- c("1", "2")

    expect_identical(
        prior_predict(model, design, draws = draws, seed = 301L),
        expected
    )
})

test_that("prior shrinkage samplers cover every supported family", {
    specifications <- list(
        fixed_scale(precision = 2),
        gamma_precision(shape = 2, rate = 3),
        global_half_cauchy(scale = 2),
        local_half_cauchy(scale = 2),
        horseshoe(global_scale = 2, local_scale = 3),
        multiplicative_gamma(shape = 2, rate = 3)
    )

    set.seed(18L)
    states <- lapply(
        specifications,
        prior_draw_shrinkage,
        n_entities = 4L,
        n_dimensions = 3L
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
    expect_identical(dim(states[[4L]]$local), c(4L, 3L))
    expect_identical(dim(states[[5L]]$local), c(4L, 3L))
})

test_that("prior_draw maps raw nodes to named entity values", {
    structure <- list(
        iid = TRUE,
        nodes = c("latent", "A", "B", "C"),
        modeled_index = c(2L, 3L, 4L),
        modeled_scale = c(2, 4, 5),
        entity_names = c("A", "B", "C")
    )
    shrinkage <- list(global = c(1, 1), local = NULL)

    set.seed(42L)
    draw <- prior_draw(structure, shrinkage, n_dim = 2L, mean = c(1, 2, 3))

    expect_identical(dim(draw$value), c(3L, 2L))
    expect_identical(
        dimnames(draw$value),
        list(c("A", "B", "C"), c("dim_1", "dim_2"))
    )
    expect_equal(
        unname(draw$value),
        draw$raw[structure$modeled_index, , drop = FALSE] /
            structure$modeled_scale + c(1, 2, 3)
    )
})

test_that("prior_predict supports structured component priors", {
    design <- prior_test_design()
    Q <- Matrix::Diagonal(2L, x = c(1, 2))
    dimnames(Q) <- list(c("A", "B"), c("A", "B"))
    operator <- Matrix::sparseMatrix(
        i = c(1, 1, 2, 2),
        j = c(1, 2, 1, 2),
        x = c(1, -1, -1, 1),
        dimnames = list(c("A", "B"), c("A", "B"))
    )
    hierarchy <- data.frame(
        node = c("root", "A", "B"),
        parent = c(NA, "root", "root"),
        stringsAsFactors = FALSE
    )
    structures <- list(
        iid(),
        precision(Q),
        gmrf(operator, ridge = 1),
        tree(node, parent, hierarchy)
    )

    for (component_structure in structures) {
        model <- combo_model(
            rank = 1L,
            global_intercept = fixed_intercept(0),
            cell_offset = combo_gaussian_component(
                shrinkage = fixed_scale(),
                structure = component_structure
            ),
            cell_factors = NULL,
            treatment_offset = NULL,
            treatment_main_factors = NULL,
            treatment_interaction_factors = NULL
        )
        predictions <- prior_predict(model, design, draws = 2L, seed = 7L)
        expect_identical(dim(predictions), c(2L, 4L))
        expect_true(all(is.finite(predictions)))
    }
})

test_that("prior_predict supports metadata means and nested doses", {
    design <- prior_test_design()
    cell_data <- data.frame(cell = c("A", "B"), feature = c(-1, 1))
    compound_data <- data.frame(drug = c("X", "Y"), feature = c(1, -1))
    model <- combo_model(
        rank = 1L,
        global_intercept = fixed_intercept(0),
        cell_offset = combo_gaussian_component(
            shrinkage = fixed_scale(),
            mean = ~ 0 + feature,
            beta_precision = gamma_precision()
        ),
        cell_factors = NULL,
        treatment_offset = combo_gaussian_component(
            shrinkage = gamma_precision(),
            mean = ~ 0 + feature,
            beta_precision = fixed_scale()
        ),
        treatment_main_factors = NULL,
        treatment_interaction_factors = NULL,
        dose = nested(precision = 2)
    )

    predictions <- prior_predict(
        model,
        design,
        cell_data = cell_data,
        compound_data = compound_data,
        draws = 3L,
        seed = 44L
    )
    expect_identical(dim(predictions), c(3L, 4L))
    expect_true(all(is.finite(predictions)))
})

test_that("prior_predict applies scalar and named beta offsets", {
    design <- data.frame(
        cell = c("A", "B", "C"),
        drug_1 = c("X", "Y", "Z"),
        dose_1 = 1,
        drug_2 = c(NA, "X", "Y"),
        dose_2 = c(NA, 1, 1),
        stringsAsFactors = FALSE
    )
    cell_data <- data.frame(
        cell = c("A", "B", "C"),
        cell_a = c(-1, 0, 1),
        cell_b = c(1, 1, -2)
    )
    compound_data <- data.frame(
        drug = c("X", "Y", "Z"),
        compound_a = c(2, -1, 0),
        compound_b = c(0, 1, 3)
    )
    mean_component <- function(formula) {
        combo_gaussian_component(
            shrinkage = fixed_scale(),
            mean = formula,
            beta_precision = fixed_scale()
        )
    }
    model <- combo_model(
        rank = 1L,
        global_intercept = fixed_intercept(0),
        cell_offset = mean_component(~ 0 + cell_a + cell_b),
        cell_factors = NULL,
        treatment_offset = mean_component(~ 0 + compound_a + compound_b),
        treatment_main_factors = NULL,
        treatment_interaction_factors = NULL
    )
    compiled <- compile_combo_design(
        model,
        design,
        cell_data = cell_data,
        compound_data = compound_data
    )

    scalar <- prior_resolve_beta_offsets(
        compiled,
        list(cell_offset = 0.4, treatment_offset = -0.2)
    )
    expect_equal(scalar$cell_offset, c(cell_a = 0.4, cell_b = 0.4))
    expect_equal(
        scalar$treatment_offset,
        c(compound_a = -0.2, compound_b = -0.2)
    )
    cell_only <- prior_resolve_beta_offsets(
        compiled,
        list(cell_offset = 0.4)
    )
    expect_equal(
        cell_only$treatment_offset,
        c(compound_a = 0, compound_b = 0)
    )

    offset <- list(
        cell_offset = c(cell_a = 0.5),
        treatment_offset = c(compound_b = -0.25)
    )
    resolved <- prior_resolve_beta_offsets(compiled, offset)
    expect_equal(resolved$cell_offset, c(cell_a = 0.5, cell_b = 0))
    expect_equal(
        resolved$treatment_offset,
        c(compound_a = 0, compound_b = -0.25)
    )

    baseline <- prior_predict(
        model,
        design,
        cell_data = cell_data,
        compound_data = compound_data,
        draws = 2L,
        seed = 72L
    )
    shifted <- prior_predict(
        model,
        design,
        cell_data = cell_data,
        compound_data = compound_data,
        draws = 2L,
        seed = 72L,
        beta_offset = offset
    )
    cell_effect <- as.numeric(
        compiled$components$cell_offset$mean$X %*% resolved$cell_offset
    )
    compound_effect <- c(0, as.numeric(
        compiled$components$treatment_offset$mean$X %*%
            resolved$treatment_offset
    ))
    expected_shift <- cell_effect[compiled$all_cell] +
        compound_effect[compiled$all_treatment_1 + 1L] +
        compound_effect[compiled$all_treatment_2 + 1L]
    expected <- matrix(expected_shift, nrow = 2L, ncol = 3L, byrow = TRUE)
    colnames(expected) <- colnames(shifted)
    expect_equal(shifted - baseline, expected)
})

test_that("prior_predict validates beta offsets", {
    design <- prior_test_design()
    model <- prior_intercept_only_model(fixed_intercept(0))

    expect_error(prior_predict(model, design, beta_offset = 1), "named list")
    expect_error(prior_predict(model, design, beta_offset = list(0)), "unique, nonempty")
    expect_error(
        prior_predict(model, design, beta_offset = list(other = 0)),
        "unknown entries"
    )
    expect_error(
        prior_predict(
            model,
            design,
            beta_offset = structure(
                list(0, 0),
                names = c("cell_offset", "cell_offset")
            )
        ),
        "unique, nonempty"
    )
    expect_error(
        prior_predict(model, design, beta_offset = list(cell_offset = c(0, 1))),
        "scalar or have coefficient names"
    )
    expect_error(
        prior_predict(model, design, beta_offset = list(cell_offset = NULL)),
        "finite numeric"
    )
    expect_error(
        prior_predict(model, design, beta_offset = list(cell_offset = NA_real_)),
        "finite numeric"
    )
    expect_error(
        prior_predict(model, design, beta_offset = list(cell_offset = Inf)),
        "finite numeric"
    )
    expect_error(
        prior_predict(model, design, beta_offset = list(cell_offset = "zero")),
        "finite numeric"
    )
    expect_error(
        prior_predict(
            model,
            design,
            beta_offset = list(
                cell_offset = c(feature = 0, feature = 0)
            )
        ),
        "unique, nonempty coefficient names"
    )
    expect_error(
        prior_predict(model, design, beta_offset = list(cell_offset = 1)),
        "nonzero but its component has no coefficients"
    )
    expect_error(
        prior_predict(
            model,
            design,
            beta_offset = list(cell_offset = c(feature = 0))
        ),
        "unknown coefficients"
    )
})

test_that("prior prediction permits interaction self-combinations", {
    design <- prior_test_design()
    model <- combo_model(rank = 2L, global_intercept = fixed_intercept(0))

    expect_no_error(
        predictions <- prior_predict(model, design, draws = 2L, seed = 8L)
    )
    expect_identical(dim(predictions), c(2L, 4L))
})

test_that("prior_predict validates draw controls", {
    design <- prior_test_design()
    model <- prior_intercept_only_model(fixed_intercept(0))

    expect_error(prior_predict(list(), design), "combo_model")
    expect_error(prior_predict(model, design, draws = 0), "at least 1")
    expect_error(prior_predict(model, design, draws = 1.5), "at least 1")
    expect_error(prior_predict(model, design, seed = 0), "seed must be NULL")
})
