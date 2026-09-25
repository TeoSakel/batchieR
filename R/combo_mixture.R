# Prediction-only combinations of mean-response basins.

#' Compare chains through expected-response projections
#'
#' Summarizes predictions by chain and by contiguous blocks, without clustering
#' or fitting predictive weights. Does not require the optional package `loo`.
#'
#' @inheritParams predictive_mixture
#' @details
#' Chains are sorted by their original IDs and draws by iteration. Each chain
#' uses its first `m` retained draws, where `m` is the shortest chain length.
#' Predictions are processed in chunks; the full draws-by-reference matrix is
#' not retained. The transformation is applied before averaging. This function
#' does not advance the caller's RNG, including when a transform is rejected.
#'
#' Each chain is divided into four contiguous blocks whose sizes differ by at
#' most one. `block_distances` contains all six RMS distances between block
#' projections per chain. `noise_floor` is their pooled 95th percentile (R's
#' default quantile type 7), bounded below by machine epsilon. This empirical
#' within-chain variation floor can reflect both Monte Carlo variation and
#' drift. It is not a calibrated Monte Carlo standard error, a convergence test,
#' or a dendrogram cut height. Blocks need not be independent.
#'
#' Use `stats::dist(mean) / sqrt(ncol(mean))` to obtain RMS distances for
#' [stats::hclust()]. Equal mean projections can conceal different predictive
#' distributions or unvisited regions. Native-scale spread summaries are
#' therefore returned separately; predictive variance uses the empirical
#' population variance of expected responses plus mean observation variance.
#' @return A list with:
#' \describe{
#'   \item{chains}{Sorted original chain IDs.}
#'   \item{mean}{Chain-by-reference matrix of transformed expected responses.}
#'   \item{block_mean}{Chain-by-block-by-reference array of transformed means.}
#'   \item{block_distances}{Data frame with `chain`, `block_1`, `block_2`, and
#'     RMS `distance`.}
#'   \item{noise_floor}{Empirical within-chain variation floor.}
#'   \item{native_mean, predictive_variance, noise_variance}{Native-scale chain
#'     summaries, with corresponding `block_native_mean`,
#'     `block_predictive_variance`, and `block_noise_variance`.}
#'   \item{reference_grid}{The prediction rows used, in matrix column order.
#'     Explicit duplicate rows retain their weight.}
#'   \item{chain_draws}{Per-chain counts: `retained`, `used`, and `excluded`.}
#'   \item{source_draws}{Retained `chain`, `draw`, `iteration`, and `block`
#'     mappings, in analysis order.}
#' }
#' Matrix row names identify original chains; reference dimensions follow
#' `reference_grid` row order. No posterior snapshots are returned.
#' @seealso [predictive_mixture()], [posterior_epred()]
#' @examples
#' \dontrun{
#' p <- chain_projections(fit, response_transform = plogis)
#' tree <- stats::hclust(stats::dist(p$mean) / sqrt(ncol(p$mean)),
#'                       method = "average")
#' groups <- stats::cutree(tree, k = 2)
#' plot(tree)
#' }
#' @export
chain_projections <- function(fit, reference_grid = NULL, response_transform = identity) {
    validate_combo_fit(fit)
    rng <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit({
        if (is.null(rng)) {
            if (exists(".Random.seed", globalenv(), inherits = FALSE)) {
                rm(".Random.seed", envir = globalenv())
            }
        } else assign(".Random.seed", rng, envir = globalenv())
    }, add = TRUE)
    prepared <- combo_mixture_projection_inputs(fit, reference_grid, response_transform)
    fit <- prepared$fit
    projection <- combo_mixture_projection(fit, prepared$reference_grid, response_transform)
    noise <- combo_mixture_block_noise(projection$block_mean, projection$chains)
    m <- prepared$excluded$used[1L]
    source_draws <- data.frame(
        chain = fit$chain_id,
        draw = fit$draw_id,
        iteration = fit$iteration,
        block = rep(combo_mixture_blocks(m), length(projection$chains))
    )
    c(
        projection,
        noise,
        list(
            reference_grid = prepared$reference_grid,
            chain_draws = prepared$excluded,
            source_draws = source_draws
        )
    )
}

#' Stack predictive basins from a multi-chain fit
#'
#' Groups chains by mean-response predictions and combines their predictive
#' distributions using PSIS-LOO stacking. This is a prediction-optimized
#' approximation, not a posterior over one common parameterization.
#'
#' @param fit A multi-chain `combo_fit`, with at least 20 retained draws in
#'   every chain. Unequal chains use their first common number of draws.
#' @param reference_grid Optional prediction rows. The default deduplicates
#'   observed training configurations, including modeled observation predictors.
#'   Explicit rows have equal weight, including duplicates.
#' @param response_transform A deterministic elementwise function, applied to
#'   each expected-response draw before averaging for clustering only. It must
#'   preserve dimensions and return finite numeric values.
#' @param ndraws Number of representatives per basin and final pseudo-draws.
#'   Defaults to the smallest chain size; must be between the largest basin's
#'   chain count and the smallest chain size.
#' @details
#' Four contiguous blocks per chain provide a 95th-percentile empirical noise
#' floor for RMS distances between transformed mean-response projections.
#' Average linkage is cut at the largest adjacent merge-height ratio when it
#' is at least three, with the lower height bounded by the noise floor. Two
#' chains split at three times that floor. This is a separation heuristic, not
#' a calibrated hypothesis test. Equal projections can hide different noise,
#' uncertainty, distribution shape, or unvisited basins.
#'
#' Stacking uses all observed training rows and retained analysis draws, with
#' chain-aware relative efficiencies. Partitions and compiled preprocessing
#' remain fixed. In particular, an empirical intercept is not re-estimated
#' after leaving out each observation. Weights are predictive optimization
#' weights, not posterior basin probabilities. Pareto-k warnings indicate that
#' finite returned weights may be unreliable; no automatic refitting occurs.
#'
#' Representatives and final allocations are deterministic. Largest-remainder
#' allocation approximates fitted weights with equally weighted pseudo-draws;
#' small positive weights can receive no draws. All prediction methods refer
#' to this realized empirical mixture, on the original model response scale.
#' Construction and noiseless prediction do not advance the caller's RNG.
#' User transforms are checked for repeatability and elementwise behavior on
#' the evaluated inputs; arbitrary function purity cannot be proven.
#' @return A `combo_predictive_mixture` with prediction snapshots, `membership`,
#'   `weights` (fitted and realized), `source_draws`, `representatives`,
#'   `clustering`, `projection`, `psis`, and `provenance`. `psis$pointwise`
#'   retains original observed-row IDs and diagnostics. It does not inherit
#'   from `combo_fit` and does not support parameter inference or active learning.
#' @seealso [chain_projections()], [posterior_predictions], [log_lik()]
#' @export
predictive_mixture <- function(fit, reference_grid = NULL,
                               response_transform = identity, ndraws = NULL) {
    validate_combo_fit(fit)
    if (!requireNamespace("loo", quietly = TRUE)) {
        cli::cli_abort("{.pkg loo} is required for {.fn predictive_mixture}.")
    }
    # Restore RNG even if a user-supplied transform is random or raises an error.
    rng <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
    on.exit({
        if (is.null(rng)) {
            if (exists(".Random.seed", globalenv(), inherits = FALSE)) {
                rm(".Random.seed", envir = globalenv())
            }
        } else assign(".Random.seed", rng, envir = globalenv())
    }, add = TRUE)
    combo_mixture_compute(fit, reference_grid, response_transform, ndraws,
        paste(deparse(substitute(response_transform)), collapse = " "))
}

# Internal providers permit the offline artifact builder to reuse predictions
# and identical basin likelihood calculations. Public calls use the workers
# below; both routes share preparation, clustering, weighting and allocation.
combo_mixture_compute <- function(fit, reference_grid, response_transform, ndraws,
    transform_label, project = combo_mixture_projection, compute_psis = combo_mixture_psis) {
    prepared <- combo_mixture_projection_inputs(fit, reference_grid, response_transform)
    fit <- prepared$fit
    reference_grid <- prepared$reference_grid
    m <- min(table(fit$chain_id))
    if (!is.null(ndraws)) {
        ndraws <- param_positive_integer(ndraws, "ndraws")
        if (ndraws > m) cli::cli_abort("{.arg ndraws} cannot exceed the smallest retained chain size ({m}).")
    }
    projection <- project(fit, reference_grid, response_transform)
    clustering <- combo_mixture_cluster(projection$mean, projection$block_mean)
    membership <- data.frame(chain = projection$chains, basin = clustering$basin)
    if (!is.null(ndraws) && ndraws < max(table(membership$basin))) {
        cli::cli_abort("{.arg ndraws} cannot be smaller than the largest basin's chain count.")
    }
    psis <- compute_psis(fit, membership)
    weights <- combo_mixture_weights(psis$lpd)
    combo_mixture_finish(fit, reference_grid, projection, clustering, psis,
                         weights, ndraws, prepared$excluded,
                         transform_label)
}

combo_mixture_projection_inputs <- function(fit, reference_grid, response_transform) {
    prepared <- combo_mixture_prepare(fit)
    if (!is.function(response_transform)) {
        cli::cli_abort("{.arg response_transform} must be a function.")
    }
    if (is.null(reference_grid)) reference_grid <- combo_mixture_reference(prepared$fit)
    combo_prediction_indices(prepared$fit, reference_grid)
    if (!nrow(reference_grid)) cli::cli_abort("{.arg reference_grid} must contain rows.")
    prepared$reference_grid <- reference_grid
    prepared
}

combo_mixture_prepare <- function(fit) {
    n <- length(fit$draws)
    for (name in c("chain_id", "draw_id", "iteration")) {
        x <- fit[[name]]
        if (!is.numeric(x) || length(x) != n || any(!is.finite(x)) ||
            any(x < 1 | x != floor(x))) {
            cli::cli_abort("Invalid fit {.field {name}} metadata.")
        }
    }
    chains <- sort(unique(fit$chain_id))
    if (length(chains) < 2L) cli::cli_abort("At least two chains are required.")
    groups <- lapply(chains, function(chain) {
        ids <- which(fit$chain_id == chain)
        if (anyDuplicated(fit$iteration[ids]) || anyDuplicated(fit$draw_id[ids])) {
            cli::cli_abort("Draw and iteration IDs must be unique within each chain.")
        }
        ids[order(fit$iteration[ids])]
    })
    m <- min(lengths(groups))
    if (m < 20L) cli::cli_abort("At least 20 retained draws per chain are required.")
    positions <- unlist(lapply(groups, utils::head, m), use.names = FALSE)
    excluded <- data.frame(chain = chains, retained = lengths(groups), used = m,
                           excluded = lengths(groups) - m)
    fit <- combo_mixture_subset(fit, positions)
    precision <- vapply(fit$draws, function(x) x$precision, numeric(1))
    if (any(!is.finite(precision) | precision <= 0)) {
        cli::cli_abort("Every draw must have positive finite observation precision.")
    }
    list(fit = fit, excluded = excluded)
}

combo_mixture_subset <- function(fit, positions) {
    fit$draws <- fit$draws[positions]
    for (name in c("chain_id", "draw_id", "iteration")) fit[[name]] <- fit[[name]][positions]
    fit$sampling$chains <- length(unique(fit$chain_id))
    fit
}

combo_mixture_reference <- function(fit) {
    rows <- fit$compiled$observed_rows
    index <- combo_prediction_indices(fit, NULL)
    key <- data.frame(cell = index$cell[rows],
                      first = pmin(index$treatment_1[rows], index$treatment_2[rows]),
                      second = pmax(index$treatment_1[rows], index$treatment_2[rows]))
    key <- cbind(key, as.data.frame(index$mean_design[rows, , drop = FALSE]))
    fit$compiled$data[rows[!duplicated(key)], , drop = FALSE]
}

combo_mixture_transform <- function(x, transform) {
    y <- transform(x)
    if (!is.numeric(y) || !identical(dim(y), dim(x)) || any(!is.finite(y)) ||
        !identical(y, transform(x))) {
        cli::cli_abort("{.arg response_transform} must be deterministic, finite, and shape-preserving.")
    }
    flat <- transform(as.numeric(x))
    if (!is.numeric(flat) || !is.null(dim(flat)) || length(flat) != length(x) ||
        !identical(as.numeric(y), as.numeric(flat))) {
        cli::cli_abort("{.arg response_transform} must act elementwise.")
    }
    # A whole-vector operation can preserve shape: check separated inputs too.
    probe <- unique(c(1L, length(x) %/% 2L + 1L, length(x)))
    scalar <- vapply(as.numeric(x)[probe], function(value) {
        result <- transform(value)
        if (!is.numeric(result) || length(result) != 1L) {
            cli::cli_abort("{.arg response_transform} must act elementwise.")
        }
        result
    }, numeric(1))
    if (!identical(as.numeric(y)[probe], unname(scalar))) {
        cli::cli_abort("{.arg response_transform} must act elementwise.")
    }
    y
}

combo_mixture_blocks <- function(n) {
    rep(seq_len(4L), length.out = n)[order(rep(seq_len(4L), length.out = n))]
}

combo_mixture_projection <- function(fit, reference, transform, chunk_size = 256L,
                                      predict_expected = posterior_epred) {
    chains <- sort(unique(fit$chain_id))
    positions <- lapply(chains, function(id) which(fit$chain_id == id))
    blocks <- lapply(positions, function(ids) split(ids, combo_mixture_blocks(length(ids))))
    nc <- length(chains)
    nr <- nrow(reference)
    means <- native_mean <- variance <- matrix(NA_real_, nc, nr,
                                               dimnames = list(chains, NULL))
    block_mean <- block_native_mean <- block_variance <- array(NA_real_, c(nc, 4L, nr))
    noise <- vapply(fit$draws, function(x) 1 / x$precision, numeric(1))
    for (start in seq.int(1L, nr, by = chunk_size)) {
        rows <- seq.int(start, min(nr, start + chunk_size - 1L))
        mu <- predict_expected(fit, newdata = reference[rows, , drop = FALSE])
        if (any(!is.finite(mu))) cli::cli_abort("Non-finite expected-response predictions.")
        z <- combo_mixture_transform(mu, transform)
        for (c in seq_len(nc)) {
            ids <- positions[[c]]
            means[c, rows] <- colMeans(z[ids, , drop = FALSE])
            native_mean[c, rows] <- colMeans(mu[ids, , drop = FALSE])
            variance[c, rows] <- colMeans(sweep(mu[ids, , drop = FALSE], 2L,
                                                native_mean[c, rows], "-")^2) + mean(noise[ids])
            for (b in seq_len(4L)) {
                selected <- blocks[[c]][[b]]
                block_mean[c, b, rows] <- colMeans(z[selected, , drop = FALSE])
                block_native_mean[c, b, rows] <- colMeans(mu[selected, , drop = FALSE])
                block_variance[c, b, rows] <- colMeans(sweep(mu[selected, , drop = FALSE], 2L,
                    block_native_mean[c, b, rows], "-")^2) + mean(noise[selected])
            }
        }
    }
    list(
        chains = chains,
        mean = means,
        block_mean = block_mean,
        native_mean = native_mean,
        predictive_variance = variance,
        block_native_mean = block_native_mean,
        block_predictive_variance = block_variance,
        noise_variance = vapply(positions, function(ids) mean(noise[ids]), numeric(1)),
        block_noise_variance = t(vapply(blocks, function(bs) { vapply(bs, function(ids) mean(noise[ids]), numeric(1))}, numeric(4)))
    )
}

combo_mixture_labels <- function(group, chains) {
    groups <- unique(group)
    counts <- vapply(groups, function(g) sum(group == g), integer(1))
    first <- vapply(groups, function(g) min(chains[group == g]), numeric(1))
    groups <- groups[order(-counts, first)]
    paste0("B", match(group, groups))
}

combo_mixture_block_noise <- function(blocks, chains) {
    pairs <- utils::combn(seq_len(4L), 2L)
    distances <- do.call(rbind, lapply(seq_along(chains), function(i) {
        block <- matrix(blocks[i, , ], nrow = 4L)
        data.frame(
            chain = chains[i],
            block_1 = pairs[1L, ],
            block_2 = pairs[2L, ],
            distance = as.numeric(stats::dist(block)) / sqrt(ncol(block))
        )
    }))
    if (any(!is.finite(distances$distance))) {
        cli::cli_abort("Non-finite within-chain block distances.")
    }
    noise_floor <- unname(stats::quantile(distances$distance, .95))
    list(block_distances = distances, noise_floor = max(noise_floor, .Machine$double.eps))
}

combo_mixture_cluster <- function(means, blocks) {
    chains <- as.numeric(rownames(means))
    if (is.null(rownames(means))) chains <- seq_len(nrow(means))
    ordering <- order(chains)
    chains <- chains[ordering]
    means <- means[ordering, , drop = FALSE]
    blocks <- blocks[ordering, , , drop = FALSE]
    distance <- stats::dist(means) / sqrt(ncol(means))
    block_noise <- combo_mixture_block_noise(blocks, chains)
    floor <- block_noise$noise_floor
    tree <- stats::hclust(distance, method = "average")
    partition <- function(threshold) {
        heights <- tree$height
        cut <- Inf
        ratio <- if (length(heights) == 1L) heights / floor else {
            heights[-1L] / pmax(heights[-length(heights)], floor)
        }
        j <- which.max(ratio)
        if (ratio[j] >= threshold) {
            lower <- if (length(heights) == 1L) floor else max(heights[j], floor)
            upper <- if (length(heights) == 1L) heights else heights[j + 1L]
            cut <- sqrt(lower) * sqrt(upper)
        }
        group <- if (is.finite(cut)) stats::cutree(tree, h = cut) else rep(1L, length(chains))
        list(basin = combo_mixture_labels(group, chains), cut = cut,
             gap_ratio = ratio[j], gap_index = j)
    }
    primary <- partition(3)
    sensitivity <- lapply(c(1.5, 3, 6), function(value) {
        data.frame(setting = paste0("threshold_", value), chain = chains,
                   basin = partition(value)$basin)
    })
    if (is.finite(primary$cut)) {
        for (multiplier in c(.5, 2)) {
            sensitivity[[length(sensitivity) + 1L]] <- data.frame(
                setting = paste0("cut_", multiplier), chain = chains,
                basin = combo_mixture_labels(stats::cutree(tree, h = primary$cut * multiplier), chains))
        }
    }
    c(primary, list(chains = chains, distance = distance, tree = tree,
                    noise_floor = floor, block_distances = block_noise$block_distances$distance,
                    sensitivity = do.call(rbind, sensitivity)))
}

combo_mixture_lme <- function(x) {
    maximum <- apply(x, 2L, max)
    result <- maximum
    finite <- is.finite(maximum)
    if (any(finite)) result[finite] <- maximum[finite] +
        log(colMeans(exp(sweep(x[, finite, drop = FALSE], 2L, maximum[finite], "-"))))
    result
}

combo_mixture_loo_chunk <- function(ll, chain_id) {
    if (any(!is.finite(ll))) cli::cli_abort("Non-finite training log likelihoods.")
    likelihood <- exp(sweep(ll, 2L, apply(ll, 2L, max), "-"))
    constant <- apply(likelihood, 2L, function(x) all(x == x[1L]))
    r_eff <- rep(1, ncol(ll))
    warnings <- character()
    if (any(!constant)) {
        r_eff[!constant] <- withCallingHandlers(
            loo::relative_eff(likelihood[, !constant, drop = FALSE],
                chain_id = match(chain_id, unique(chain_id)), cores = 1),
            warning = function(w) {
                warnings <<- unique(c(warnings, conditionMessage(w)))
                invokeRestart("muffleWarning")
            })
    }
    if (any(!is.finite(r_eff) | r_eff <= 0)) cli::cli_abort("Undefined PSIS relative efficiency.")
    # Aggregate numerical PSIS diagnostics ourselves, rather than repeating one
    # warning for each chunk. Unexpected warnings are preserved in the result.
    value <- withCallingHandlers(loo::loo(ll, r_eff = r_eff, cores = 1),
        warning = function(w) {
            warnings <<- unique(c(warnings, conditionMessage(w)))
            invokeRestart("muffleWarning")
        })
    lpd <- value$pointwise[, "elpd_loo"]
    if (any(!is.finite(lpd))) cli::cli_abort("Non-finite PSIS predictive densities.")
    list(lpd = lpd, k = loo::pareto_k_values(value),
         n_eff = loo::psis_n_eff_values(value), r_eff = r_eff, warnings = warnings)
}

combo_mixture_psis <- function(fit, membership, chunk_size = 256L,
                               log_likelihood = log_lik) {
    observed <- combo_observed_data(fit)
    basins <- unique(membership$basin)
    basins <- basins[order(as.integer(sub("B", "", basins)))]
    lpd <- matrix(NA_real_, length(observed$response), length(basins),
                  dimnames = list(as.character(observed$row_id), basins))
    points <- vector("list", length(basins))
    warnings <- character()
    for (b in seq_along(basins)) {
        ids <- which(fit$chain_id %in% membership$chain[membership$basin == basins[b]])
        subset <- combo_mixture_subset(fit, ids)
        current <- vector("list", ceiling(nrow(lpd) / chunk_size))
        for (chunk in seq_along(current)) {
            rows <- seq.int((chunk - 1L) * chunk_size + 1L, min(nrow(lpd), chunk * chunk_size))
            ll <- log_likelihood(subset, newdata = observed$data[rows, , drop = FALSE])
            result <- combo_mixture_loo_chunk(ll, subset$chain_id)
            lpd[rows, b] <- result$lpd
            current[[chunk]] <- data.frame(basin = basins[b], row_id = observed$row_id[rows],
                elpd_loo = result$lpd, pareto_k = result$k, psis_n_eff = result$n_eff,
                r_eff = result$r_eff, draws = length(ids))
            warnings <- unique(c(warnings, result$warnings))
        }
        points[[b]] <- do.call(rbind, current)
    }
    combo_mixture_psis_result(lpd, do.call(rbind, points), warnings)
}

combo_mixture_psis_result <- function(lpd, pointwise, warnings = character()) {
    pointwise$threshold <- pmin(.7, 1 - 1 / log10(pointwise$draws))
    invalid <- is.na(pointwise$pareto_k) | is.na(pointwise$psis_n_eff)
    # Infinite k is a diagnostic failure (including constant-tail fits), not a
    # reason to replace finite LOO densities or to hide the observation.
    pointwise$flagged <- invalid | pointwise$pareto_k > pointwise$threshold
    diagnostics <- do.call(rbind, lapply(colnames(lpd), function(b) {
        p <- pointwise[pointwise$basin == b, ]
        data.frame(basin = b, observations = nrow(p), draws = p$draws[1L],
            threshold = p$threshold[1L], maximum_k = if (anyNA(p$pareto_k)) NA_real_ else max(p$pareto_k),
            flagged = sum(p$flagged), count_ge_1 = sum(p$pareto_k >= 1, na.rm = TRUE),
            fraction_over_0_7 = mean(p$pareto_k > .7, na.rm = TRUE),
            elpd_loo = sum(p$elpd_loo))
    }))
    reliable <- !any(pointwise$flagged)
    if (!reliable) cli::cli_warn(c(
        "PSIS reliability threshold exceeded or undefined for {sum(pointwise$flagged)} basin-observation entries.",
        "i" = "Inspect {.field psis$diagnostics} and {.field psis$pointwise}; finite stacking weights are returned with unreliable status.",
        "i" = "{sum(diagnostics$count_ge_1)} entries have Pareto k >= 1; {sum(diagnostics$fraction_over_0_7 > .01, na.rm = TRUE)} basins have more than 1% above 0.7."
    ), class = "combo_mixture_psis_warning")
    list(lpd = lpd, pointwise = pointwise, diagnostics = diagnostics,
         reliable = reliable, warnings = warnings)
}

combo_mixture_weights <- function(lpd) {
    if (!is.matrix(lpd) || !nrow(lpd) || !ncol(lpd) || any(!is.finite(lpd))) {
        cli::cli_abort("Stacking requires finite pointwise predictive densities.")
    }
    columns <- lapply(seq_len(ncol(lpd)), function(i) unname(lpd[, i]))
    unique_columns <- !duplicated(columns)
    group <- match(columns, columns[unique_columns])
    # Optimize on unique columns so indistinguishable basins share mass equally.
    x <- lpd[, unique_columns, drop = FALSE]
    weight <- if (ncol(x) == 1L) 1 else combo_mixture_optimize(x)
    result <- weight[group] / tabulate(group)[group]
    stats::setNames(result, colnames(lpd))
}

combo_mixture_optimize <- function(lpd) {
    # Pairwise coordinate ascent on the concave simplex objective. Exact
    # one-dimensional line searches can reach a boundary without a barrier or
    # a vanishing logit gradient. Row centering avoids likelihood underflow.
    x <- exp(sweep(lpd, 1L, apply(lpd, 1L, max), "-"))
    weights <- rep(1 / ncol(x), ncol(x))
    for (iteration in seq_len(10000L)) {
        density <- as.numeric(x %*% weights)
        derivative <- colMeans(x / density)
        if (any(!is.finite(derivative))) break
        # Weighted average derivative equals one. This gap bounds remaining
        # improvement in average log score by concavity, including on faces.
        if (max(derivative) <= 1 + 1e-7) return(weights / sum(weights))
        to <- which.max(derivative)
        active <- which(weights > 0)
        from <- active[which.min(derivative[active])]
        if (from == to) break
        direction <- x[, to] - x[, from]
        upper <- weights[from]
        slope <- function(amount) {
            candidate <- density + amount * direction
            if (any(candidate <= 0)) return(-Inf)
            mean(direction / candidate)
        }
        endpoint <- slope(upper)
        amount <- if (endpoint >= 0) upper else {
            stats::uniroot(slope, c(0, upper), f.lower = mean(direction / density),
                f.upper = endpoint, tol = min(1e-12, upper * 1e-6))$root
        }
        if (!is.finite(amount) || amount <= 0) break
        weights[from] <- weights[from] - amount
        weights[to] <- weights[to] + amount
    }
    cli::cli_abort("Stacking optimizer failed to meet its simplex optimality tolerance after {iteration} iterations.")
}

combo_mixture_even <- function(ids, n) {
    if (!n) return(integer())
    ids[as.integer(round(seq.int(1L, length(ids), length.out = n)))]
}

combo_mixture_balanced <- function(chain_id, positions, n) {
    groups <- split(positions, factor(chain_id[positions], levels = sort(unique(chain_id[positions]))))
    allocation <- rep(n %/% length(groups), length(groups))
    if (n %% length(groups)) allocation[seq_len(n %% length(groups))] <-
        allocation[seq_len(n %% length(groups))] + 1L
    if (any(allocation > lengths(groups))) cli::cli_abort("Insufficient balanced representatives.")
    unlist(Map(combo_mixture_even, groups, allocation), use.names = FALSE)
}

combo_mixture_allocation <- function(weights, n) {
    exact <- n * weights
    count <- floor(exact)
    remainder <- as.integer(n - sum(count))
    if (remainder) {
        ids <- order(-(exact - count), seq_along(weights))[seq_len(remainder)]
        count[ids] <- count[ids] + 1L
    }
    as.integer(count)
}

combo_mixture_finish <- function(fit, reference, projection, clustering, psis,
                                 weights, ndraws, excluded, transform_label) {
    m <- min(table(fit$chain_id))
    n <- if (is.null(ndraws)) as.integer(m) else param_positive_integer(ndraws, "ndraws")
    membership <- data.frame(chain = projection$chains, basin = clustering$basin)
    if (n < max(table(membership$basin)) || n > m) {
        cli::cli_abort("{.arg ndraws} must be between the largest basin's chain count and the smallest retained chain size ({m}).")
    }
    basins <- names(weights)
    pools <- lapply(basins, function(b) combo_mixture_balanced(fit$chain_id,
        which(fit$chain_id %in% membership$chain[membership$basin == b]), n))
    names(pools) <- basins
    counts <- combo_mixture_allocation(weights, n)
    ids <- unlist(Map(function(pool, count) combo_mixture_balanced(fit$chain_id, pool, count),
                     pools, counts), use.names = FALSE)
    source <- data.frame(basin = rep(basins, counts), chain = fit$chain_id[ids],
                         draw = fit$draw_id[ids], iteration = fit$iteration[ids])
    representatives <- do.call(rbind, lapply(basins, function(b) {
        i <- pools[[b]]
        data.frame(basin = b, chain = fit$chain_id[i], draw = fit$draw_id[i], iteration = fit$iteration[i])
    }))
    snapshots <- lapply(fit$draws[ids], function(draw) {
        draw <- draw[c("intercept", "beta", "precision", "components")]
        draw$components <- lapply(draw$components, function(x) if (is.null(x)) NULL else x["value"])
        draw
    })
    # Prediction needs entity/mean mappings and original input rows, not prior
    # precision systems, sampler state, or parameter-extraction metadata.
    compiled <- fit$compiled[c("data", "observed_rows", "response", "cells", "treatments",
        "all_cell", "all_treatment_1", "all_treatment_2", "mean")]
    structure(list(model = fit$model[c("rank", "family")], compiled = compiled,
        draws = snapshots, membership = membership,
        weights = data.frame(basin = basins, fitted = unname(weights), count = counts,
            realized = counts / n, allocation_error = counts / n - weights,
            omitted_positive = weights > 0 & counts == 0),
        source_draws = source, representatives = representatives,
        reference_grid = reference, clustering = clustering, projection = projection,
        psis = psis, provenance = list(format_version = 1L, algorithm_version = "1.0.0",
            package_version = as.character(utils::packageVersion("batchieR")),
            loo_version = as.character(utils::packageVersion("loo")),
            source_sampling = fit$sampling, chain_draws = excluded,
            response_transform = transform_label, observed_row_ids = fit$compiled$observed_rows,
            preprocessing = "Fixed partitions and compiled preprocessing",
            empirical_intercept = identical(fit$model$mean$type, "empirical"))),
        class = "combo_predictive_mixture")
}

#' Methods for prediction-only mixtures
#' @param x,object A `combo_predictive_mixture`.
#' @param ... Reserved for future methods.
#' @return Printing returns its input invisibly. Summary returns mixture
#'   membership, fitted/realized weights, PSIS diagnostics and provenance; it
#'   does not calculate parameter convergence diagnostics.
#' @name combo_predictive_mixture_methods
NULL

#' @rdname combo_predictive_mixture_methods
#' @export
print.combo_predictive_mixture <- function(x, ...) {
    cli::cli_text("<combo_predictive_mixture> {length(x$draws)} prediction pseudo-draws, {nrow(x$weights)} mean-response basins")
    print(x$weights, row.names = FALSE)
    cli::cli_text("PSIS status: {if (x$psis$reliable) 'passed reported thresholds' else 'unreliable; inspect diagnostics'}.")
    cli::cli_text("Prediction-optimized approximation; weights are not posterior basin probabilities.")
    if (isTRUE(x$provenance$empirical_intercept)) cli::cli_text("PSIS conditions on the fitted empirical intercept.")
    invisible(x)
}

#' @rdname combo_predictive_mixture_methods
#' @export
summary.combo_predictive_mixture <- function(object, ...) {
    list(membership = object$membership, weights = object$weights,
         diagnostics = object$psis$diagnostics, reliable = object$psis$reliable,
         sensitivity = object$clustering$sensitivity, provenance = object$provenance)
}

# A private adapter shares existing prediction validation/workers only. It is
# never exposed as a posterior fit and contains no invented chain metadata.
combo_mixture_prediction_fit <- function(object) {
    structure(object[c("model", "compiled", "draws")], class = "combo_fit")
}

#' @rdname posterior_predictions
#' @export
posterior_linpred.combo_predictive_mixture <- function(object, newdata = NULL, ...) {
    posterior_linpred(combo_mixture_prediction_fit(object), newdata = newdata, ...)
}
#' @rdname posterior_predictions
#' @export
posterior_epred.combo_predictive_mixture <- function(object, newdata = NULL, ...) {
    posterior_epred(combo_mixture_prediction_fit(object), newdata = newdata, ...)
}
#' @rdname posterior_predictions
#' @export
posterior_predict.combo_predictive_mixture <- function(object, newdata = NULL, ...) {
    posterior_predict(combo_mixture_prediction_fit(object), newdata = newdata, ...)
}
#' @rdname posterior_predictions
#' @export
predict.combo_predictive_mixture <- function(object, newdata = NULL,
    type = c("draws", "mean"), observation = FALSE, ...) {
    stats::predict(combo_mixture_prediction_fit(object), newdata = newdata,
            type = type, observation = observation, ...)
}
#' @rdname log_lik
#' @export
log_lik.combo_predictive_mixture <- function(object, newdata = NULL, ...) {
    log_lik(combo_mixture_prediction_fit(object), newdata = newdata, ...)
}

combo_mixture_reject <- function(...) {
    cli::cli_abort(c("This operation is unavailable for a prediction-only mixture.",
        "i" = "It is a prediction-optimized approximation, not a posterior over one common parameterization. Use the original combo_fit for parameter inference or active learning."))
}

#' @exportS3Method posterior::as_draws
as_draws.combo_predictive_mixture <- function(x, ...) combo_mixture_reject()
#' @exportS3Method posterior::as_draws_array
as_draws_array.combo_predictive_mixture <- function(x, ...) combo_mixture_reject()
#' @exportS3Method posterior::as_draws_matrix
as_draws_matrix.combo_predictive_mixture <- function(x, ...) combo_mixture_reject()
#' @exportS3Method posterior::as_draws_df
as_draws_df.combo_predictive_mixture <- function(x, ...) combo_mixture_reject()
#' @exportS3Method posterior::as_draws_list
as_draws_list.combo_predictive_mixture <- function(x, ...) combo_mixture_reject()
#' @exportS3Method posterior::as_draws_rvars
as_draws_rvars.combo_predictive_mixture <- function(x, ...) combo_mixture_reject()
#' @exportS3Method loo::loo
loo.combo_predictive_mixture <- function(x, ...) combo_mixture_reject()
#' @exportS3Method stats::update
update.combo_predictive_mixture <- function(object, ...) combo_mixture_reject()
#' @export
tidy.combo_predictive_mixture <- function(x, ...) combo_mixture_reject()
#' @export
glance.combo_predictive_mixture <- function(x, ...) combo_mixture_reject()
