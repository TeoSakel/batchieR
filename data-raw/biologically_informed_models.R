# Regenerate the triplicate screen from the repository root with:
# Rscript data-raw/biologically_informed_models.R
# Replicates share one prior parameter draw and have independent observation noise.
devtools::load_all()
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")

cells <- sprintf("C%02d", 1:8)
drugs <- sprintf("D%02d", 1:8)

cell_data <- data.frame(
    cell = cells,
    lineage = rep(c("lineage_A", "lineage_B"), each = 4),
    pathway_activity = c(-1.4, -1.0, -0.6, -0.2, 0.2, 0.6, 1.0, 1.4),
    growth_rate = c(-0.7, 0.8, -0.2, 0.4, 0.7, -0.8, 0.2, -0.4),
    stringsAsFactors = FALSE
)
compound_data <- data.frame(
    drug = drugs,
    mechanism = rep(sprintf("MOA_%d", 1:4), each = 2),
    target_signature = c(-1.5, -1.2, -0.6, -0.3, 0.3, 0.6, 1.2, 1.5),
    physchem_pc = c(0.6, -0.7, -0.4, 0.8, -0.8, 0.4, 0.7, -0.6),
    stringsAsFactors = FALSE
)

feature_mean <- formula_mean(
    ~ pathway_activity + growth_rate + target_signature + physchem_pc,
    beta_precision = 1
)

cell_tree <- data.frame(
    node = c(".cell_root", "lineage_A", "lineage_B", cells),
    parent = c(NA, ".cell_root", ".cell_root", rep("lineage_A", 4), rep("lineage_B", 4)),
    edge_length = c(1, 1, 1, rep(0.05, 8)),
    stringsAsFactors = FALSE
)
cell_tree <- tree(node, parent, cell_tree, edge_length = edge_length)

compound_tree <- data.frame(
    node = c(".drug_root", sprintf("MOA_%d", 1:4), drugs),
    parent = c(NA, rep(".drug_root", 4), rep(sprintf("MOA_%d", 1:4), each = 2)),
    edge_length = c(1, rep(1, 4), rep(0.05, 8)),
    stringsAsFactors = FALSE
)
compound_tree <- tree(node, parent, compound_tree, edge_length = edge_length)

tree_model <- combo_model(
    rank = 2,
    mean = feature_mean,
    cell_offset = combo_gaussian_component(
        structure = iid(),
        shrinkage = gamma_precision()
    ),
    cell_factors = combo_gaussian_component(
        structure = cell_tree,
        shrinkage = multiplicative_gamma()
    ),
    treatment_offset = combo_gaussian_component(
        structure = iid(),
        shrinkage = global_half_cauchy()
    ),
    treatment_main_factors = combo_gaussian_component(
        structure = compound_tree,
        shrinkage = global_half_cauchy()
    ),
    treatment_interaction_factors = combo_gaussian_component(
        structure = compound_tree,
        shrinkage = global_half_cauchy()
    ),
    dose = nested(precision = 4)
)


modeled_pairs <- c(
    # within pairs
    "D01:D02", "D03:D04",
    "D05:D06", "D07:D08",
    # cross pairs
    "D02:D03", "D04:D05",
    "D06:D07", "D01:D08",
    "D01:D03", "D03:D05",
    "D05:D07", "D01:D07"
)

screen <- expand.grid(
    cell = cells,
    drug_1 = drugs,
    dose_1 = 1:2,
    drug_2 = c(NA_character_, drugs),
    dose_2 = c(NA_real_, 1:2),
    replicate = 1:3,
    stringsAsFactors = FALSE
)
is_single_drug <- with(screen, is.na(drug_2) & is.na(dose_2))
is_modeled_pair <- with(
    screen,
    !is.na(dose_2) & paste(drug_1, drug_2, sep = ":") %in% modeled_pairs
)
screen <- screen[is_single_drug | is_modeled_pair, , drop = FALSE]
rownames(screen) <- NULL
screen$configuration <- rep(seq_len(sum(screen$replicate == 1)), times = 3)

# Use the tree model with the generating mean and observation precision.
synthetic_model <- tree_model
synthetic_model$mean <- formula_mean(
    ~ pathway_activity + growth_rate + target_signature + physchem_pc,
    beta_mean = c(
        pathway_activity = 0.42,
        growth_rate = 0,
        target_signature = -0.28,
        physchem_pc = 0
    ),
    beta_precision = 1000
)
synthetic_model$family <- gaussian_response(
    precision = gamma_precision(shape = 100, rate = 1.44)
)

screen$modeled_response <- as.numeric(
    prior_predict(
        synthetic_model,
        screen,
        cell_data = cell_data,
        compound_data = compound_data,
        draws = 1,
        type = "response",
        seed = 20260828,
        intercept = 1.15
    )
)
screen$viability <- plogis(screen$modeled_response)

# Retain the generating means for reproducibility, not for held-out RMSE.
screen$true_mean <- as.numeric(prior_predict(
    synthetic_model, screen, cell_data = cell_data, compound_data = compound_data,
    draws = 1, type = "mean", seed = 20260828, intercept = 1.15
))
write.csv(screen, "inst/extdata/biologically_informed_models.csv", row.names = FALSE)
