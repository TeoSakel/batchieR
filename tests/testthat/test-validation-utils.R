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
