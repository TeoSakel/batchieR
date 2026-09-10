sb_csv <- function(data, path) {
    temporary <- tempfile(".csv-", tmpdir = dirname(path))
    on.exit(unlink(temporary))
    utils::write.csv(data, temporary, row.names = FALSE, na = "")
    if (!file.rename(temporary, path)) cli::cli_abort("Cannot commit {path}.")
}

sb_bind <- function(rows) {
    rows <- Filter(function(x) !is.null(x) && nrow(x) > 0L, rows)
    if (!length(rows)) return(data.frame())
    result <- do.call(rbind, rows)
    rownames(result) <- NULL
    result
}

sb_collect <- function(output) {
    manifest <- readRDS(file.path(output, "manifest.rds"))
    jobs <- readRDS(file.path(output, "jobs.rds"))
    records <- lapply(jobs$id, function(id) {
        path <- file.path(output, "fits", paste0(id, ".rds"))
        if (file.exists(path)) readRDS(path) else NULL
    })
    status <- jobs
    status$status <- vapply(records, function(r) if (is.null(r)) "pending" else r$status, character(1))
    status$error <- vapply(records, function(r) if (is.null(r$error)) "" else r$error, character(1))
    status$warnings <- vapply(records, function(r) paste(r$warnings, collapse = " | "), character(1))
    collect <- function(field) sb_bind(lapply(seq_along(records), function(i) {
        r <- records[[i]]
        if (is.null(r) || r$status != "ok" || is.null(r[[field]])) return(NULL)
        cbind(jobs[rep(i, nrow(r[[field]])), , drop = FALSE], r[[field]])
    }))
    tables <- lapply(c("metrics", "diagnostics", "ranks", "coverage"), collect)
    names(tables) <- c("metrics", "diagnostics", "ranks", "coverage")
    designs <- lapply(unique(sprintf("%s-%03d", jobs$scenario, jobs$dataset)), function(id) {
        path <- file.path(output, "designs", paste0(id, ".rds"))
        if (!file.exists(path)) return(NULL)
        d <- readRDS(path)
        if (!is.null(d$error)) return(NULL)
        data.frame(dataset = id, rows = d$n_rows, cells = d$cells, treatments = d$treatments,
                   plates = d$plates, observed_25 = d$observed[["reveal_25"]], observed_75 = d$observed[["reveal_75"]],
                   response_min = unname(d$response_quantiles[1L]), response_median = unname(d$response_quantiles[4L]),
                   response_max = unname(tail(d$response_quantiles, 1L)), boundary_fraction = d$fraction_near_boundary)
    })
    traces <- sb_bind(lapply(seq_along(records), function(i) {
        r <- records[[i]]
        if (is.null(r) || r$status != "ok" || jobs$dataset[i] != 1L) return(NULL)
        quantities <- dimnames(r$draws)[[3L]]
        select <- unique(c("sigma", quantities[startsWith(quantities, "mean_row_")][1L], "observed_log_likelihood"))
        sb_bind(lapply(select, function(quantity) {
            draw <- r$draws[, , quantity, drop = FALSE]
            data.frame(adapter = jobs$adapter[i], scenario = jobs$scenario[i], mask = jobs$mask[i],
                       quantity = quantity, iteration = rep(seq_len(dim(draw)[1L]), dim(draw)[2L]),
                       chain = rep(seq_len(dim(draw)[2L]), each = dim(draw)[1L]), value = as.numeric(draw))
        }))
    }))
    rank_summary <- ecdfs <- envelopes <- list()
    if (nrow(tables$ranks)) {
        groups <- split(tables$ranks, interaction(tables$ranks$scenario, tables$ranks$mask,
                        tables$ranks$adapter, tables$ranks$quantity, drop = TRUE))
        for (r in groups) {
            good <- r[r$status == "eligible", , drop = FALSE]
            rank_summary[[length(rank_summary) + 1L]] <- data.frame(
                scenario = r$scenario[1L], mask = r$mask[1L], adapter = r$adapter[1L], quantity = r$quantity[1L],
                eligible = nrow(good), unavailable = nrow(r) - nrow(good),
                expected = manifest$specification$config$datasets,
                missing_fits = manifest$specification$config$datasets - nrow(r))
            if (!nrow(good)) next
            key <- paste(nrow(good), good$m[1L], sep = "-")
            if (is.null(envelopes[[key]])) envelopes[[key]] <- sb_envelope(nrow(good), good$m[1L],
                seed = sb_seed(manifest$specification$config$seed, "report-envelope", key))
            envelope <- envelopes[[key]]
            ecdfs[[length(ecdfs) + 1L]] <- cbind(
                good[rep(1L, nrow(envelope)), c("scenario", "mask", "adapter", "quantity"), drop = FALSE],
                envelope, difference = sb_ecdf(good$rank, good$m[1L]), eligible = nrow(good),
                incomplete = nrow(good) < manifest$specification$config$datasets)
        }
    }
    checks_path <- file.path(output, "self-checks.rds")
    c(list(manifest = manifest, status = status, designs = sb_bind(designs), traces = traces,
           rank_summary = sb_bind(rank_summary), ecdfs = sb_bind(ecdfs),
           self_checks = if (file.exists(checks_path)) readRDS(checks_path) else NULL,
           completion = if (any(status$status == "pending")) "PARTIAL / INCOMPLETE" else
               if (any(status$status == "failed")) "COMPLETED WITH FAILURES" else "COMPLETED",
           generated = format(Sys.time(), tz = "UTC", usetz = TRUE)), tables)
}

sb_table <- function(data, digits = 3L) {
    if (is.null(data) || !nrow(data)) cat("\nUnavailable: no applicable completed results.\n\n") else
        cat("\n", as.character(knitr::kable(data, format = "html", digits = digits, escape = TRUE)), "\n\n", sep = "")
}

sb_report_identity <- function(b) {
    spec <- b$manifest$specification
    cat("**Status: ", b$completion, "**\n\nStarted: ", b$manifest$started,
        ". Report rendered: ", b$generated, ".\n\n", sep = "")
    sb_table(sb_bind(lapply(spec$adapters, function(a) data.frame(id = a$id, sampler = a$name, description = a$description))))
    settings <- spec$config
    sb_table(data.frame(setting = names(settings), value = vapply(settings, function(x) paste(x, collapse = ", "), character(1))))
    cat("Chains run sequentially; thin = 1. ESS/second includes compilation, initialization, and warmup.\n\n")
    sb_table(as.data.frame(table(b$status$status), stringsAsFactors = FALSE))
}

sb_report_design <- function(b) {
    cat("Matched-prior scenarios use fixed intercept logit(0.8) and Gamma(25, rate=1) observation precision. ",
        "`gamma` uses Gamma(3, rate=0.3) component precisions; `multiplicative_gamma` replaces cell-factor shrinkage; ",
        "`horseshoe` uses package-default component shrinkage. Generating and fitted priors match within each.\n\n", sep = "")
    cat("`realism` uses controlled rank-2 dose-weighted signals fitted at the configured rank. These are synthetic assumptions, not empirical Merck estimates. ",
        "Main factor SD=0.3, interaction factor SD=0.4, cell factor SD=1, cell offset SD=0.2, negative treatment offset scale=0.8, noise SD=0.2.\n\n", sep = "")
    cat("Singles reveal their entire plates; remaining plates follow a shared random ordering. Mask names denote requested reveal fractions; actual observed row counts appear below. ",
        "Boundary fraction counts synthetic viabilities below 0.01 or above 0.99. Extreme finite prior draws are retained.\n\n", sep = "")
    sb_table(b$designs)
}

sb_report_efficiency <- function(b) {
    d <- b$diagnostics
    if (!nrow(d)) return(sb_table(d))
    summary <- sb_bind(lapply(split(d, d$id), function(x) {
        finite <- function(v, fun) if (any(is.finite(v))) fun(v[is.finite(v)]) else NA_real_
        data.frame(id = x$id[1L], rhat_max = if (any(!is.na(x$rhat))) max(x$rhat, na.rm = TRUE) else NA_real_,
                   bulk_ess_min = finite(x$ess_bulk, min), tail_ess_min = finite(x$ess_tail, min),
                   median_bulk_ess_per_second = finite(x$ess_bulk_per_second, stats::median),
                   mcse_mean_max = finite(x$mcse_mean, max),
                   unavailable = sum(!is.finite(x$rhat) | !is.finite(x$ess_bulk) | !is.finite(x$ess_tail)))
    }))
    sb_table(summary)
    if (length(unique(d$adapter)) > 1L) {
        base <- if ("gibbs" %in% d$adapter) "gibbs" else sort(unique(d$adapter))[1L]
        keys <- c("scenario", "dataset", "mask", "variable")
        reference <- d[d$adapter == base, c(keys, "ess_bulk_per_second")]
        names(reference)[ncol(reference)] <- "reference_ess_per_second"
        paired <- merge(d[d$adapter != base, ], reference, by = keys)
        paired$ratio <- paired$ess_bulk_per_second / paired$reference_ess_per_second
        cat("Paired bulk ESS/second ratios relative to ", base, " (same dataset, mask, quantity).\n\n", sep = "")
        sb_table(paired[c("scenario", "dataset", "mask", "adapter", "variable", "ratio")])
    } else cat("One candidate sampler: paired sampler comparisons are unavailable.\n\n")
    if (nrow(b$traces)) {
        cat("Representative traces: first dataset per scenario, post-warmup iteration indices.\n\n")
        for (scenario in unique(b$traces$scenario)) {
            x <- b$traces[b$traces$scenario == scenario, ]
            print(ggplot2::ggplot(x, ggplot2::aes(iteration, value, colour = factor(chain))) +
                ggplot2::geom_line(linewidth = 0.3) +
                ggplot2::facet_wrap(~ adapter + mask + quantity, scales = "free_y", ncol = 3) +
                ggplot2::labs(title = scenario, colour = "Chain", x = "Retained iteration") + ggplot2::theme_bw())
        }
    }
}

sb_report_recovery <- function(b) {
    cat("Latent and interaction recovery covers the complete design; response prediction covers unrevealed plates. ",
        "Viability prediction averages transformed noisy predictive draws. Generating truth measures recovery; it is not an exact posterior reference.\n\n", sep = "")
    sb_table(b$metrics)
    if (nrow(b$coverage)) sb_table(b$coverage[c("id", "target", "level", "coverage", "mean_width", "n")])
}

sb_report_sbc <- function(b) {
    cat("Ranks compare generating truth with M=100 temporally spaced posterior draws. Eligibility requires full-chain R-hat <=1.01, ",
        "bulk and tail ESS >=100, and basic ESS >=90% of the balanced spaced sample size. These checks approximate independence; they do not prove it. ",
        "There is no adaptive extension. Fixed constants are excluded and ties randomized reproducibly.\n\n", sep = "")
    cat("Each panel uses one rank per independent dataset, separately by quantity, scenario, mask, and sampler. ",
        "95% simultaneous bands use 10,000 discrete-uniform simulations and cover rank support within a panel; they do not adjust across panels. ",
        "Observed-data log-likelihood includes Gaussian normalization and exactly the observed rows.\n\n", sep = "")
    if (b$manifest$specification$config$profile != "calibration") cat("**Exploratory run: this profile cannot establish calibration.**\n\n")
    cat("**Panels with failed, pending, or ineligible datasets are inconclusive: the successful subset may be selected.** ",
        "The realism track is not matched-prior SBC. No automatic calibration-pass claim is made.\n\n", sep = "")
    sb_table(b$rank_summary)
    if (nrow(b$ranks)) sb_table(as.data.frame(table(b$ranks$status), stringsAsFactors = FALSE))
    if (!nrow(b$ecdfs)) return(cat("No eligible SBC ranks; ECDF and histogram plots are unavailable.\n\n"))
    groups <- split(b$ecdfs, interaction(b$ecdfs$scenario, b$ecdfs$mask, b$ecdfs$adapter, drop = TRUE))
    for (name in names(groups)) {
        x <- groups[[name]]
        title <- paste(name, if (any(x$incomplete)) "— inconclusive subset" else "— diagnostic")
        print(ggplot2::ggplot(x, ggplot2::aes(rank, difference)) +
            ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper), fill = "grey85") +
            ggplot2::geom_hline(yintercept = 0, colour = "grey50") + ggplot2::geom_step() +
            ggplot2::facet_wrap(~ quantity, ncol = 3) +
            ggplot2::labs(title = title, y = "ECDF − discrete-uniform CDF", x = "SBC rank (0–100)") + ggplot2::theme_bw())
        r <- b$ranks[b$ranks$status == "eligible" & b$ranks$scenario == x$scenario[1L] &
                        b$ranks$mask == x$mask[1L] & b$ranks$adapter == x$adapter[1L], ]
        r$bin <- pmin(9L, floor(r$rank * 10 / 101))
        print(ggplot2::ggplot(r, ggplot2::aes(bin)) + ggplot2::geom_bar() +
            ggplot2::facet_wrap(~ quantity, ncol = 3) + ggplot2::scale_x_continuous(breaks = 0:9, limits = c(-0.5, 9.5)) +
            ggplot2::labs(title = title, x = "Rank bin (first: 11 ranks; others: 10)", y = "Datasets") + ggplot2::theme_bw())
    }
}

sb_report_checks <- function(b) {
    if (is.null(b$self_checks)) cat("Analytic self-checks unavailable: run interrupted before completion.\n\n") else {
        cat("Analytic harness self-check: **", if (b$self_checks$passed) "passed" else "FAILED", "**.\n\n", b$self_checks$note, "\n\n", sep = "")
        sb_table(b$self_checks$summary)
    }
    sb_table(b$status[b$status$status != "ok" | nzchar(b$status$warnings), ])
    cat("Fast mixing, recovery, and SBC alone do not prove correctness. Numerical warnings and invalid quantities remain visible. ",
        "CSV tables and compact draws are retained alongside this report.\n\n", sep = "")
    spec <- b$manifest$specification
    cat("R: ", spec$R, ". Platform: ", spec$platform, ". Quarto: ", spec$quarto, ".\n\n", sep = "")
    sb_table(data.frame(package = names(spec$versions), version = unname(spec$versions)))
    cat("<details><summary>Source fingerprints</summary>\n\n")
    sb_table(data.frame(file = names(spec$source), md5 = unname(spec$source)))
    cat("</details>\n\nMethods: [Stan SBC guidance](https://mc-stan.org/docs/stan-users-guide/simulation-based-calibration.html); ",
        "[Modrák et al., choice of SBC test quantities](https://arxiv.org/abs/2211.02383).\n\n", sep = "")
}
