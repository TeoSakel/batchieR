# Ranks are discrete on 0:M. Random tie breaking also handles discrete quantities.
sb_rank <- function(truth, draws) {
    if (!is.finite(truth) || any(!is.finite(draws))) return(NA_integer_)
    less <- sum(draws < truth)
    equal <- sum(draws == truth)
    as.integer(less + if (equal) sample.int(equal + 1L, 1L) - 1L else 0L)
}

sb_envelope <- function(n, m = 100L, simulations = 10000L, seed = 1L, prob = 0.95) {
    if (n < 1L || m < 1L || simulations < 1L) cli::cli_abort("Envelope sizes must be positive.")
    null_cdf <- seq_len(m + 1L) / (m + 1L)
    maximum <- sb_with_seed(seed, replicate(simulations, {
        counts <- tabulate(sample.int(m + 1L, n, replace = TRUE), nbins = m + 1L)
        max(abs(cumsum(counts) / n - null_cdf))
    }))
    radius <- unname(stats::quantile(maximum, prob, type = 1L))
    data.frame(rank = 0:m, null_cdf = null_cdf,
               lower = pmax(0, null_cdf - radius) - null_cdf,
               upper = pmin(1, null_cdf + radius) - null_cdf,
               radius = radius)
}

sb_ecdf <- function(ranks, m = 100L) {
    counts <- tabulate(as.integer(ranks) + 1L, nbins = m + 1L)
    cumsum(counts) / length(ranks) - seq_len(m + 1L) / (m + 1L)
}

sb_summary_draws <- function(draws) {
    suppressWarnings(as.data.frame(posterior::summarise_draws(
        posterior::as_draws_array(draws), "mean", "sd", "rhat", "ess_bulk", "ess_tail", "mcse_mean")))
}

sb_ranks <- function(draws, truth, diagnostics, seed, matched_prior, m = 100L) {
    quantities <- dimnames(draws)[[3L]]
    result <- data.frame(quantity = quantities, truth = unname(truth[quantities]),
                         rank = NA_integer_, m = m, status = "insufficient_draws",
                         spaced_ess = NA_real_, stringsAsFactors = FALSE)
    if (!matched_prior) {
        result$status <- "not_matched_prior"
        return(result)
    }
    chains <- dim(draws)[2L]
    per_chain <- m %/% chains
    counts <- rep(per_chain, chains)
    if (m %% chains) counts[seq_len(m %% chains)] <- per_chain + 1L
    if (per_chain < 2L || any(counts > dim(draws)[1L])) return(result)
    # Balanced temporal spacing without using parameter values to select draws.
    positions <- lapply(counts, function(n) unique(round(seq(1L, dim(draws)[1L], length.out = n))))
    sb_with_seed(seed, {
        for (q in seq_along(quantities)) {
            full <- diagnostics[match(quantities[q], diagnostics$variable), ]
            values <- draws[, , q]
            if (!is.finite(result$truth[q]) || any(!is.finite(values))) {
                result$status[q] <- "nonfinite"
                next
            }
            if (!is.finite(full$rhat) || full$rhat > 1.01 ||
                !is.finite(full$ess_bulk) || !is.finite(full$ess_tail) ||
                min(full$ess_bulk, full$ess_tail) < m) {
                result$status[q] <- "inadequate_convergence_or_ess"
                next
            }
            spaced <- lapply(seq_len(chains), function(ch) values[positions[[ch]], ch])
            # Estimate independence with equal-length chain segments, even for
            # overrides where M is not exactly divisible by the chain count.
            equal <- vapply(spaced, head, numeric(per_chain), n = per_chain)
            ess <- suppressWarnings(posterior::ess_basic(equal))
            result$spaced_ess[q] <- ess
            if (!is.finite(ess) || ess < 0.9 * length(equal)) {
                result$status[q] <- "residual_autocorrelation"
                next
            }
            result$rank[q] <- sb_rank(result$truth[q], unlist(spaced, use.names = FALSE))
            result$status[q] <- "eligible"
        }
        result
    })
}

# Test-only reference adapters have a common scalar normal posterior interface.
sb_reference_adapters <- function() {
    list(
        exact = function(mu, sd, m) stats::rnorm(m, mu, sd),
        prior_only = function(mu, sd, m) stats::rnorm(m),
        shifted = function(mu, sd, m) stats::rnorm(m, mu + 3 * sd, sd),
        narrow = function(mu, sd, m) stats::rnorm(m, mu, 0.2 * sd)
    )
}

sb_self_checks <- function(seed, datasets = 300L, m = 100L) {
    adapters <- sb_reference_adapters()
    ranks <- sb_with_seed(seed, {
        result <- vector("list", datasets * length(adapters))
        index <- 0L
        for (i in seq_len(datasets)) {
            theta <- stats::rnorm(1L)
            y <- stats::rnorm(20L, theta, 1)
            posterior_sd <- 1 / sqrt(21)
            posterior_mean <- sum(y) / 21
            truth_ll <- sum(stats::dnorm(y, theta, 1, log = TRUE))
            for (id in names(adapters)) {
                draws <- adapters[[id]](posterior_mean, posterior_sd, m)
                ll <- vapply(draws, function(value) sum(stats::dnorm(y, value, 1, log = TRUE)), numeric(1))
                index <- index + 1L
                result[[index]] <- data.frame(adapter = id, dataset = i,
                    parameter_rank = sb_rank(theta, draws), likelihood_rank = sb_rank(truth_ll, ll))
            }
        }
        do.call(rbind, result)
    })
    envelope <- sb_envelope(datasets, m, seed = sb_seed(seed, "envelope"))
    summary <- do.call(rbind, lapply(names(adapters), function(id) {
        r <- ranks[ranks$adapter == id, ]
        parameter <- max(abs(sb_ecdf(r$parameter_rank, m)))
        likelihood <- max(abs(sb_ecdf(r$likelihood_rank, m)))
        data.frame(adapter = id, parameter_deviation = parameter,
                   likelihood_deviation = likelihood, radius = envelope$radius[1L],
                   parameter_flag = parameter > envelope$radius[1L],
                   likelihood_flag = likelihood > envelope$radius[1L])
    }))
    # An exact sampler can cross a 95% band by chance. Never call that a software
    # test failure; detect gross reference failures and strong deliberate faults.
    exact <- summary[summary$adapter == "exact", ]
    sensitivity <- with(summary, all(likelihood_flag[adapter == "prior_only"]) &&
                            all(parameter_flag[adapter %in% c("shifted", "narrow")]))
    list(seed = seed, datasets = datasets, m = m, ranks = ranks, summary = summary,
         passed = sensitivity && max(exact$parameter_deviation, exact$likelihood_deviation) < 0.2,
         note = paste("Exact normal posterior: prior N(0,1), 20 observations with variance 1.",
                      "A reference crossing a 95% band is possible by chance; these are sensitivity checks,",
                      "not evidence that the combination-model sampler is calibrated."))
}
