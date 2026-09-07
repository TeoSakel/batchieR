test_that("experiment validation allows repeated cell identifiers", {
    experiments <- data.frame(
        cell = c("A", "A"),
        drug_1 = c("X", "Y"),
        dose_1 = c(1, 1),
        drug_2 = c(NA, NA),
        dose_2 = c(NA, NA),
        response = c(0.2, 0.3)
    )

    validated <- validate_experiments(experiments)

    expect_identical(validated$cell, c("A", "A"))
    experiments$cell[2L] <- NA
    expect_error(
        validate_experiments(experiments),
        "nonmissing and nonempty"
    )
})

test_that("numeric identifiers are stable across vector contexts", {
    expect_identical(
        combo_id_text(c(0.35, 10))[1L],
        combo_id_text(c(0.35, 1.08))[1L]
    )
    expect_identical(combo_id_text(-0), combo_id_text(0))
})

test_that("identical treatments in different positions compile once", {
    experiments <- data.frame(
        cell = rep("A", 3L),
        drug_1 = c("X", "Y", "Z"),
        dose_1 = c(0.35, 10, 0.01),
        drug_2 = c("Y", "X", NA),
        dose_2 = c(1.08, 0.35, NA),
        response = c(0.2, 0.3, 0.4)
    )

    compiled <- compile_combo_model(
        combo_model(rank = 1L),
        experiments
    )

    expect_equal(nrow(compiled$treatments), 4L)
    expect_equal(sum(compiled$treatments$drug == "X"), 1L)
})

test_that("clip supports elementwise lower bounds", {
    expect_equal(
        clip(c(0, 2), lower = c(1, 3), upper = 4),
        c(1, 3)
    )
    expect_error(
        clip(1, lower = c(0, 2), upper = 1),
        "lower must be less than or equal to upper",
        class = "rlang_error"
    )
})

test_that("rhcauchy draws nonnegative half-Cauchy values", {
    set.seed(92L)
    observed <- rhcauchy(4L, scale = 2)
    set.seed(92L)
    expected <- abs(stats::rcauchy(4L, location = 0, scale = 2))

    expect_identical(observed, expected)
    expect_true(all(observed >= 0))
})

test_that("numerical fallback reports a cli warning", {
    expect_warning(
        value <- rmvnorm_safe(
            matrix(-1, nrow = 1L),
            mu_part = 0,
            fallback = 7,
            label = "test precision"
        ),
        "Numerical instability in test precision",
        class = "rlang_warning"
    )
    expect_identical(value, 7)
})

test_that("number predicates separate finite scalars from positivity", {
    for (value in list(-2L, -0.5, 0, 0.5, 2L, .Machine$double.xmax)) {
        expect_identical(is_number(value), TRUE)
        expect_identical(is_positive_number(value), value > 0)
    }
    invalid <- list(NULL, numeric(), c(1, 2), NA_real_, NaN, Inf, -Inf,
        TRUE, FALSE, "1", list(1), factor("1"), 1 + 1i)
    for (value in invalid) {
        for (predicate in list(is_number, is_positive_number, is_integer,
            is_positive_integer, is_nonnegative_integer)) {
            expect_identical(predicate(value), FALSE)
        }
    }
})

test_that("integer predicates require representable scalar integers", {
    for (value in c(-.Machine$integer.max, -1, 0, 1, .Machine$integer.max)) {
        expect_identical(is_integer(value), TRUE)
        expect_identical(is_positive_integer(value), value > 0)
        expect_identical(is_nonnegative_integer(value), value >= 0)
    }
    for (value in c(-0.5, 0.5, .Machine$integer.max + 1, -.Machine$integer.max - 1)) {
        expect_identical(is_integer(value), FALSE)
        expect_identical(is_positive_integer(value), FALSE)
        expect_identical(is_nonnegative_integer(value), FALSE)
    }
})

test_that("numeric parameter helpers retain conversion and targeted conditions", {
    expect_identical(param_scalar_positive(2L, "scale"), 2)
    expect_identical(param_positive_integer(2, "count"), 2L)
    for (value in list(0, -1, NA_real_, Inf, NULL, c(1, 2), TRUE, "1", list(1), 1 + 1i)) {
        expect_error(param_scalar_positive(value, "scale"),
            "scale must be one finite positive number", class = "rlang_error")
    }
    expect_error(param_positive_integer(.Machine$integer.max + 1, "count"),
        "count must be one positive integer", class = "rlang_error")
    for (value in list(NULL, c(1, 2), TRUE, "1", list(1), 1 + 1i)) {
        expect_error(fixed_mean(value), "value must be one finite number",
            class = "rlang_error")
        expect_error(prior_validate_intercept(value), "intercept must be NULL or one finite numeric value",
            class = "rlang_error")
    }
    expect_identical(fixed_mean(-2L)$value, -2)
    expect_identical(prior_validate_intercept(-2L), -2)
})

test_that("logical predicates accept only nonmissing scalar flags", {
    for (value in list(TRUE, FALSE, c(flag = TRUE))) {
        expect_identical(is_logical(value), TRUE)
    }
    for (value in list(NULL, logical(), c(TRUE, FALSE), NA, NA_real_,
        0, 1L, "TRUE", list(TRUE), factor("TRUE"))) {
        expect_identical(is_logical(value), FALSE)
    }
})
