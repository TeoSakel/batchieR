sb_design <- function(cells, drugs, doses, merck = NULL) {
    columns <- c("cell", "drug_1", "dose_1", "drug_2", "dose_2", "plate")
    if (!is.null(merck)) {
        input <- readRDS(merck)
        if (!all(columns %in% names(input$training)) ||
            !all(columns %in% names(input$holdout))) {
            cli::cli_abort("Merck input must contain training and holdout design tables.")
        }
        design <- rbind(input$training[columns], input$holdout[columns])
    } else {
        cell <- sprintf("cell_%02d", seq_len(cells))
        drug <- sprintf("drug_%02d", seq_len(drugs))
        singles <- expand.grid(cell = cell, drug_1 = drug, dose_1 = seq_len(doses),
                               KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
        singles$drug_2 <- NA_character_
        singles$dose_2 <- NA_real_
        singles$plate <- paste0("single_", singles$cell, "_", singles$drug_1)
        pairs <- utils::combn(drug, 2L)
        combinations <- lapply(seq_len(ncol(pairs)), function(i) {
            rows <- expand.grid(cell = cell, dose_1 = seq_len(doses), dose_2 = seq_len(doses),
                                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
            rows$drug_1 <- pairs[1L, i]
            rows$drug_2 <- pairs[2L, i]
            rows$plate <- paste(rows$cell, rows$drug_1, rows$drug_2, sep = "__")
            rows[columns]
        })
        design <- rbind(singles[columns], do.call(rbind, combinations))
    }
    rownames(design) <- NULL
    if (anyNA(design$plate) || any(!nzchar(as.character(design$plate)))) {
        cli::cli_abort("Every design row needs a plate identifier.")
    }
    sb_internal("validate_experiments")(design, require_response = FALSE)
    design
}

sb_masks <- function(design, seed, fractions = c(0.25, 0.75)) {
    # With imported mixed plates, observing singles reveals the whole plate.
    singles <- is.na(design$drug_1) | is.na(design$drug_2)
    initial <- unique(design$plate[singles])
    candidates <- sort(setdiff(unique(design$plate), initial))
    if (length(candidates) < 2L) cli::cli_abort("Design needs at least two unobserved combination plates.")
    ordered <- sb_with_seed(seed, candidates[sample.int(length(candidates))])
    masks <- lapply(fractions, function(fraction) {
        n <- min(length(ordered) - 1L, max(1L, floor(fraction * length(ordered))))
        design$plate %in% c(initial, head(ordered, n))
    })
    names(masks) <- paste0("reveal_", as.integer(100 * fractions))
    masks
}

# Independent predictor algebra, deliberately avoiding package prediction helpers.
sb_truth_terms <- function(snapshot, indices) {
    n <- length(indices$cell)
    cell <- indices$cell
    a <- indices$treatment_1
    b <- indices$treatment_2
    components <- snapshot$components
    mean <- rep(snapshot$intercept, n)
    if (length(snapshot$beta)) {
        mean <- mean + as.numeric(indices$mean_design[, names(snapshot$beta), drop = FALSE] %*% snapshot$beta)
    }
    rows <- function(component, index) {
        rbind(rep(0, ncol(component$value)), component$value)[index + 1L, , drop = FALSE]
    }
    if (!is.null(components$cell_offset)) mean <- mean + components$cell_offset$value[cell, 1L]
    if (!is.null(components$treatment_offset)) {
        mean <- mean + as.numeric(rows(components$treatment_offset, a) + rows(components$treatment_offset, b))
    }
    interaction <- numeric(n)
    if (!is.null(components$cell_factors)) {
        w <- components$cell_factors$value[cell, , drop = FALSE]
        if (!is.null(components$treatment_main_factors)) {
            mean <- mean + rowSums(w * (rows(components$treatment_main_factors, a) +
                                            rows(components$treatment_main_factors, b)))
        }
        if (!is.null(components$treatment_interaction_factors)) {
            interaction <- rowSums(w * rows(components$treatment_interaction_factors, a) *
                                       rows(components$treatment_interaction_factors, b))
        }
    }
    list(mean = mean + interaction, interaction = interaction)
}

sb_indices <- function(compiled) {
    list(cell = compiled$all_cell, treatment_1 = compiled$all_treatment_1,
         treatment_2 = compiled$all_treatment_2, mean_design = compiled$mean$design,
         row_id = seq_along(compiled$all_cell))
}

sb_realism_snapshot <- function(compiled) {
    # Draw one drug vector, then multiply it by each drug's ordered dose weight.
    treatments <- compiled$treatments
    drugs <- unique(treatments$drug)
    drug_index <- match(treatments$drug, drugs)
    weight <- ave(as.numeric(treatments$dose), treatments$drug, FUN = function(x) {
        match(x, sort(unique(x))) / length(unique(x))
    })
    component <- function(value) list(value = value, raw = value,
                                      global_precision = rep(NA_real_, ncol(value)),
                                      local_precision = NULL)
    drug_value <- function(sd) {
        matrix(stats::rnorm(length(drugs) * 2L, sd = sd), ncol = 2L)[drug_index, , drop = FALSE] * weight
    }
    list(intercept = stats::qlogis(0.8), beta = numeric(), precision = 25,
         components = list(
             cell_offset = component(matrix(stats::rnorm(length(compiled$cells), sd = 0.2), ncol = 1L)),
             cell_factors = component(matrix(stats::rnorm(length(compiled$cells) * 2L), ncol = 2L)),
             treatment_offset = component(matrix(-0.8 * abs(stats::rnorm(length(drugs)))[drug_index] * weight, ncol = 1L)),
             treatment_main_factors = component(drug_value(0.3)),
             treatment_interaction_factors = component(drug_value(0.4))))
}

sb_generate <- function(design, scenario, rank, seed, mask_seed) {
    model <- sb_model(scenario, rank)
    generating_model <- if (scenario == "realism") sb_model(scenario, 2L) else model
    compiled <- sb_internal("compile_combo_design")(generating_model, design)
    indices <- sb_indices(compiled)
    generated <- sb_with_seed(seed, {
        snapshot <- if (scenario == "realism") sb_realism_snapshot(compiled) else
            sb_internal("prior_snapshot")(compiled, stats::qlogis(0.8))
        truth <- sb_truth_terms(snapshot, indices)
        package_mean <- sb_internal("predict_combo_draw")(snapshot, indices, generating_model$rank)
        if (!isTRUE(all.equal(truth$mean, package_mean, tolerance = 1e-10, check.attributes = FALSE))) {
            cli::cli_abort("Independent predictor disagrees with package predictor.")
        }
        response <- stats::rnorm(nrow(design), truth$mean, 1 / sqrt(snapshot$precision))
        if (any(!is.finite(c(response, truth$mean, truth$interaction, snapshot$precision)))) {
            cli::cli_abort("Non-finite generative draw; seed {seed} retained as a failure.")
        }
        list(snapshot = snapshot, mean = truth$mean, interaction = truth$interaction, response = response)
    })
    c(list(design = design, scenario = scenario, matched_prior = scenario != "realism",
           model = model, generating_model = generating_model, indices = indices,
           masks = sb_masks(design, mask_seed), seed = seed, mask_seed = mask_seed), generated)
}

sb_fit_data <- function(dataset, mask) {
    data <- dataset$design
    data$response <- ifelse(dataset$masks[[mask]], dataset$response, NA_real_)
    data
}

sb_panel <- function(design, n = 8L) {
    select <- function(rows) rows[unique(round(seq(1L, length(rows), length.out = min(n, length(rows)))))]
    combinations <- which(!is.na(design$drug_1) & !is.na(design$drug_2))
    list(mean = select(seq_len(nrow(design))), interaction = select(combinations))
}
