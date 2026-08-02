combo_test_data <- function() {
    data.frame(
        cell = c("A", "A", "A", "B", "B", "B"),
        drug_1 = c("X", "Y", "X", "X", "Y", "Y"),
        dose_1 = 1,
        drug_2 = c(NA, NA, "Y", NA, NA, "X"),
        dose_2 = c(NA, NA, 1, NA, NA, 1),
        response = c(0.20, 0.40, 0.10, 0.30, NA, 0.15),
        stringsAsFactors = FALSE
    )
}

combo_test_model <- function(structure = NULL) {
    combo_model(
        rank = 1L,
        cell_offset = combo_gaussian_component(
            shrinkage = gamma_precision(),
            structure = structure
        ),
        cell_factors = NULL,
        treatment_offset = combo_gaussian_component(
            shrinkage = gamma_precision()
        ),
        treatment_main_factors = NULL,
        treatment_interaction_factors = NULL
    )
}

combo_test_fit <- function(
    seed = 812L,
    chains = 2L,
    iter_warmup = 2L,
    iter_sampling = 4L,
    thin = 2L,
    model = combo_test_model()
) {
    fit_combo(
        model,
        combo_test_data(),
        chains = chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        thin = thin,
        seed = seed
    )
}
