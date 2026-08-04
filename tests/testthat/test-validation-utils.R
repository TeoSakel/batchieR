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
