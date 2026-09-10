# Regenerate inst/extdata/model_checking.csv from the repository root with:
# Rscript data-raw/model_checking.R
#
# This is the original toy screen for vignettes/model_checking.qmd, generated
# with batchieR 0.1.0 and seed 20260804. The vignette reads the committed CSV
# so its observations do not depend on future simulator or RNG changes.
# Regeneration is a deliberate maintenance action: rerender the vignette and
# review its interpretations before replacing the committed data.
devtools::load_all()
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
set.seed(20260804)

iid_component <- combo_gaussian_component(
    structure = iid(),
    shrinkage = gamma_precision(shape = 3, rate = 3)
)
full_model <- combo_model(
    rank = 1,
    mean = fixed_mean(0.7),
    cell_offset = iid_component,
    cell_factors = iid_component,
    treatment_offset = iid_component,
    treatment_main_factors = NULL,
    treatment_interaction_factors = iid_component
)

screen <- expand.grid(
    cell = c("A", "B", "C"),
    drug_1 = c("X", "Y", "Z"),
    dose_1 = 1:2,
    drug_2 = c(NA_character_, "X", "Y", "Z"),
    dose_2 = c(NA_real_, 1:2),
    replicate = 1:3,
    stringsAsFactors = FALSE
)
screen <- subset(
    screen,
    (is.na(drug_2) & is.na(dose_2)) |
        (sprintf("%s:%s", drug_1, drug_2) %in% c("X:Y", "X:Z", "Y:Z") &
            !is.na(dose_2))
)
screen$type <- ifelse(is.na(screen$drug_2), "single", "combination")
screen$plate <- paste0(screen$type, "_", screen$replicate)
screen$response <- prior_predict(full_model, screen, draws = 1)[1, ]

dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)
write.csv(screen, "inst/extdata/model_checking.csv", row.names = FALSE, na = "NA")
