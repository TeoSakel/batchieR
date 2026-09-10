sb_draw_array <- function(values, fit) {
    counts <- table(fit$chain_id)
    if (!length(counts) || length(unique(counts)) != 1L || nrow(values) != length(fit$draws)) {
        cli::cli_abort("Adapter must return balanced chains with matching draws.")
    }
    if (anyDuplicated(paste(fit$chain_id, fit$iteration))) {
        cli::cli_abort("Adapter returned duplicate chain/iteration pairs.")
    }
    ordering <- order(fit$chain_id, fit$iteration)
    array(values[ordering, , drop = FALSE],
          dim = c(as.integer(counts[1L]), length(counts), ncol(values)),
          dimnames = list(NULL, names(counts), colnames(values)))
}

sb_interval <- function(draws, truth, level) {
    bounds <- apply(draws, 2L, stats::quantile,
                    probs = c((1 - level) / 2, (1 + level) / 2), names = FALSE)
    if (is.null(dim(bounds))) bounds <- matrix(bounds, nrow = 2L)
    c(covered = sum(truth >= bounds[1L, ] & truth <= bounds[2L, ]),
      width = sum(bounds[2L, ] - bounds[1L, ]), n = length(truth))
}

sb_evaluate <- function(fit, dataset, mask, fit_seconds, evaluation_seed, rank_seed) {
    if (!inherits(fit, "combo_fit")) cli::cli_abort("Adapter must return a combo_fit.")
    sb_internal("validate_combo_fit")(fit)
    expected_input <- sb_fit_data(dataset, mask)
    if (!identical(fit$input$data, expected_input) || !identical(fit$model, dataset$model)) {
        cli::cli_abort("Adapter changed the supplied model or masked data.")
    }
    panel <- sb_panel(dataset$design)
    quantities <- c("sigma", paste0("mean_row_", panel$mean),
                    paste0("interaction_row_", panel$interaction), "observed_log_likelihood")
    values <- matrix(NA_real_, length(fit$draws), length(quantities), dimnames = list(NULL, quantities))
    values[, "sigma"] <- vapply(fit$draws, function(draw) 1 / sqrt(draw$precision), numeric(1))
    values[, "observed_log_likelihood"] <- 0
    truth <- c(sigma = 1 / sqrt(dataset$snapshot$precision),
               stats::setNames(dataset$mean[panel$mean], paste0("mean_row_", panel$mean)),
               stats::setNames(dataset$interaction[panel$interaction], paste0("interaction_row_", panel$interaction)))
    observed <- dataset$masks[[mask]]
    truth <- c(truth, observed_log_likelihood = sum(stats::dnorm(
        dataset$response[observed], dataset$mean[observed], truth[["sigma"]], log = TRUE)))
    sums <- c(latent = 0, interaction = 0, response = 0, viability = 0)
    levels <- c(0.5, 0.8, 0.95)
    coverage <- expand.grid(target = c("latent_mean", "interaction", "response_holdout", "viability_holdout"),
                            level = levels, stringsAsFactors = FALSE)
    coverage$covered <- coverage$width <- coverage$n <- 0
    sb_with_seed(evaluation_seed, {
        for (start in seq.int(1L, nrow(dataset$design), by = 256L)) {
            ids <- seq.int(start, min(start + 255L, nrow(dataset$design)))
            data <- dataset$design[ids, , drop = FALSE]
            means <- batchieR::posterior_epred(fit, newdata = data)
            indices <- sb_internal("combo_prediction_indices")(fit, data)
            interactions <- do.call(rbind, lapply(fit$draws, function(draw) sb_truth_terms(draw, indices)$interaction))
            # Include independent observation noise when evaluating observed viability.
            predictions <- means + matrix(stats::rnorm(length(means)), nrow(means)) * values[, "sigma"]
            for (kind in c("mean", "interaction")) {
                selected <- intersect(ids, panel[[kind]])
                if (length(selected)) values[, paste0(kind, "_row_", selected)] <-
                    (if (kind == "mean") means else interactions)[, match(selected, ids), drop = FALSE]
            }
            is_observed <- observed[ids]
            if (any(is_observed)) {
                likelihood_data <- data[is_observed, , drop = FALSE]
                likelihood_data$response <- dataset$response[ids[is_observed]]
                values[, "observed_log_likelihood"] <- values[, "observed_log_likelihood"] +
                    rowSums(batchieR::log_lik(fit, newdata = likelihood_data))
            }
            holdout <- !is_observed
            sums["latent"] <- sums["latent"] + sum((colMeans(means) - dataset$mean[ids])^2)
            sums["interaction"] <- sums["interaction"] + sum((colMeans(interactions) - dataset$interaction[ids])^2)
            sums["response"] <- sums["response"] + sum((colMeans(means)[holdout] - dataset$response[ids[holdout]])^2)
            sums["viability"] <- sums["viability"] + sum((colMeans(stats::plogis(predictions))[holdout] -
                                                             stats::plogis(dataset$response[ids[holdout]]))^2)
            matrices <- list(latent_mean = means, interaction = interactions,
                             response_holdout = predictions[, holdout, drop = FALSE],
                             viability_holdout = stats::plogis(predictions[, holdout, drop = FALSE]))
            truths <- list(latent_mean = dataset$mean[ids], interaction = dataset$interaction[ids],
                           response_holdout = dataset$response[ids[holdout]],
                           viability_holdout = stats::plogis(dataset$response[ids[holdout]]))
            for (j in seq_len(nrow(coverage))) {
                target <- coverage$target[j]
                if (length(truths[[target]])) {
                    value <- sb_interval(matrices[[target]], truths[[target]], coverage$level[j])
                    coverage[j, c("covered", "width", "n")] <-
                        coverage[j, c("covered", "width", "n")] + value
                }
            }
        }
    })
    if (any(!is.finite(values))) cli::cli_abort("Non-finite posterior test quantities.")
    draws <- sb_draw_array(values, fit)
    diagnostics <- sb_summary_draws(draws)
    diagnostics$ess_bulk_per_second <- diagnostics$ess_bulk / fit_seconds
    diagnostics$ess_tail_per_second <- diagnostics$ess_tail / fit_seconds
    ranks <- sb_ranks(draws, truth, diagnostics, rank_seed, dataset$matched_prior)
    coverage$coverage <- coverage$covered / coverage$n
    coverage$mean_width <- coverage$width / coverage$n
    # A compact panel is retained at every iteration, without thinning.
    list(draws = draws, truth = truth, diagnostics = diagnostics, ranks = ranks, coverage = coverage,
         metrics = data.frame(
             n_rows = nrow(dataset$design), n_observed = sum(observed), n_holdout = sum(!observed),
             latent_rmse = sqrt(sums[["latent"]] / nrow(dataset$design)),
             interaction_rmse = sqrt(sums[["interaction"]] / nrow(dataset$design)),
             response_holdout_rmse = sqrt(sums[["response"]] / sum(!observed)),
             viability_holdout_rmse = sqrt(sums[["viability"]] / sum(!observed)),
             fit_seconds = fit_seconds, fit_bytes = as.numeric(object.size(fit))))
}
