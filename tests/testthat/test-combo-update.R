test_that("stats::update dispatches and preserves specifications without changes", {
    model <- combo_model()
    expect_s3_class(stats::update(model), "combo_model")
    expect_identical(stats::update(model), model)
    expect_identical(update(combo_test_model()), combo_test_model())
})

test_that("updates replace every constructor argument and preserve other values", {
    component <- combo_gaussian_component(shrinkage = gamma_precision(2, 3))
    model <- combo_model(
        treatment_offset = component,
        treatment_main_factors = component,
        treatment_interaction_factors = component
    )
    replacements <- list(
        rank = 4L,
        mean = fixed_mean(1.5),
        cell_offset = component,
        cell_factors = component,
        treatment_offset = combo_gaussian_component(shrinkage = fixed_scale()),
        treatment_main_factors = combo_gaussian_component(shrinkage = fixed_scale()),
        treatment_interaction_factors = combo_gaussian_component(shrinkage = fixed_scale()),
        dose = nested(precision = 2),
        family = gaussian_response(precision = gamma_precision(3, 4))
    )
    original <- unserialize(serialize(model, NULL))
    expect_setequal(names(replacements), names(formals(combo_model)))
    for (name in names(replacements)) {
        expected <- model
        if (name %in% names(model$components)) {
            expected$components[name] <- replacements[name]
        } else {
            expected[name] <- replacements[name]
        }
        expect_identical(do.call(update, c(list(model), replacements[name])), expected)
    }
    expect_identical(model, original)
})

test_that("updates preserve disabled components and formula environments", {
    model <- local({
        transform_covariate <- function(x) x * 2
        combo_model(
            mean = formula_mean(~ transform_covariate(feature), beta_precision = 3),
            treatment_interaction_factors = NULL
        )
    })
    revised <- update(model, rank = 4)
    expect_identical(revised$rank, 4L)
    expect_identical(revised$mean, model$mean)
    expect_identical(environment(revised$mean$formula), environment(model$mean$formula))
    expect_identical(revised$components, model$components)
    expect_identical(update(model), model)
    expect_identical(update(model, mean = empirical_mean())$mean, empirical_mean())
})

test_that("components are replaced in full and can be disabled and re-enabled", {
    Q <- Matrix::Diagonal(2L)
    dimnames(Q) <- list(c("A", "B"), c("A", "B"))
    model <- combo_test_model(precision(Q))
    replacement <- combo_gaussian_component(shrinkage = fixed_scale())
    revised <- update(model, cell_offset = replacement)
    expect_identical(revised$components$cell_offset, replacement)
    expect_identical(revised$components$cell_offset$structure$type, "iid")

    disabled <- update(model, cell_offset = NULL)
    expect_named(disabled$components, names(model$components))
    expect_null(disabled$components$cell_offset)
    expect_identical(
        update(disabled, cell_offset = model$components$cell_offset),
        model
    )
})

test_that("related changes are validated together by the constructor", {
    model <- combo_model()
    revised <- update(
        model,
        cell_factors = NULL,
        treatment_main_factors = NULL,
        treatment_interaction_factors = NULL
    )
    expect_identical(revised, combo_model(
        cell_factors = NULL,
        treatment_main_factors = NULL,
        treatment_interaction_factors = NULL
    ))
    component <- combo_gaussian_component(shrinkage = gamma_precision())
    expect_identical(
        update(
            model, dose = nested(), treatment_offset = component,
            treatment_main_factors = component, treatment_interaction_factors = component
        ),
        combo_model(
            dose = nested(), treatment_offset = component,
            treatment_main_factors = component, treatment_interaction_factors = component
        )
    )
    expect_error(update(model, rank = 0), "rank")
    expect_error(update(model, mean = NULL), "mean must be")
    expect_error(update(model, dose = NULL), "dose must be")
    expect_error(update(model, family = NULL), "family must be")
    expect_error(update(model, cell_offset = 1), "cell_offset must be")
    expect_error(update(model, cell_factors = NULL), "require cell_factors")
    expect_error(
        update(model, treatment_main_factors = NULL, treatment_interaction_factors = NULL),
        "contributes nothing"
    )
    expect_error(update(model, dose = nested()), "cannot use entity-local")
    expect_identical(model, combo_model())
})

test_that("update rejects unnamed, duplicate, unknown, and abbreviated arguments", {
    model <- combo_model()
    expect_error(update(model, 4), "must be named")
    expect_error(update(model, rank = 4, NULL), "must be named")
    expect_error(update(model, rank = 4, rank = 5), "Duplicate.*rank")
    expect_error(update(model, extra = 1), "Unknown.*extra")
    expect_error(update(model, ran = 4), "Unknown.*ran")
    expect_error(update(model, formula. = ~ . + feature), "Unknown.*formula")
    expect_error(update(model, evaluate = FALSE), "Unknown.*evaluate")
    expect_error(
        do.call(update, c(list(model), setNames(list(4), NA_character_))),
        "Unknown.*NA"
    )
})

test_that("updated and directly constructed models produce identical seeded fits", {
    revised <- update(combo_model(), rank = 2, treatment_interaction_factors = NULL)
    direct <- combo_model(rank = 2, treatment_interaction_factors = NULL)
    expect_identical(revised, direct)
    revised_fit <- combo_test_fit(model = revised, chains = 1L)
    direct_fit <- combo_test_fit(model = direct, chains = 1L)
    expect_identical(revised_fit$draws, direct_fit$draws)
    expect_identical(revised_fit$chain_id, direct_fit$chain_id)
    expect_identical(revised_fit$iteration, direct_fit$iteration)
})
