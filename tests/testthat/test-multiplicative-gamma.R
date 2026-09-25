test_that("multiplicative-gamma increments accept rank-specific hyperparameters", {
    shrinkage <- multiplicative_gamma(shape = c(2, 3, 4), rate = c(1, 2, 3))
    expect_identical(shrinkage$shape, c(2, 3, 4))
    expect_identical(shrinkage$rate, c(1, 2, 3))

    model <- combo_model(
        rank = 3L,
        cell_factors = combo_gaussian_component(shrinkage = shrinkage)
    )
    compiled <- compile_combo_model(model, combo_test_data())
    expect_identical(compiled$components$cell_factors$shrinkage$shape, c(2, 3, 4))
    expect_identical(compiled$components$cell_factors$shrinkage$rate, c(1, 2, 3))

    expect_error(
        compile_combo_model(update(model, rank = 2L), combo_test_data()),
        "must have length 1 or match rank (2)", fixed = TRUE
    )
})
