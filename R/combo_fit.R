# Fitting and basic methods for combination models.

#' Fit a combination-response model
#'
#' Compiles a combination-model specification and draws posterior samples with
#' the Gibbs sampler.
#'
#' @details
#' `data` must contain the columns `cell`, `drug_1`, `dose_1`, `drug_2`,
#' `dose_2`, and `response`. Drug and dose must either both be present or both
#' be missing for each treatment position. Missing responses are allowed, but
#' at least one response must be observed.
#'
#' Each chain is initialized independently. After warmup, the sampler performs
#' `iter_sampling` transitions and retains every `thin`th state, for
#' `floor(iter_sampling / thin)` draws per chain.
#'
#' With progress enabled, successful fits end with a timing summary. Mean chain
#' time includes initialization, warmup, and sampling; total elapsed time also
#' includes worker setup and result collection, but excludes model compilation.
#'
#' @param model A combination-model specification created by `combo_model()`.
#' @param data A data frame of combination-screen observations. See Details for
#'   the required columns.
#' @param cell_data An optional data frame of cell-level covariates, containing a
#'   unique, nonmissing `cell` key and any covariates used by the global mean
#'   formula.
#' @param compound_data An optional data frame of compound-level covariates,
#'   containing a unique, nonmissing `drug` key and any covariates used by
#'   the global mean formula.
#' @param engine Sampling engine. Currently only `"gibbs"` is supported.
#' @param chains Number of independent Markov chains. Must be a positive integer.
#' @param iter_warmup Number of warmup transitions per chain. Must be a
#'   nonnegative integer.
#' @param iter_sampling Number of post-warmup transitions per chain. Must be an
#'   integer greater than or equal to `thin`.
#' @param thin Positive integer interval between retained post-warmup states.
#' @param seed `NULL`, or a positive integer used to initialize the R random-
#'   number generator.
#' @param control Named list of engine-specific control parameters. The Gibbs
#'   engine currently supports no additional control parameters, so this must
#'   be empty.
#' @param parallel_chains Maximum number of chains to run concurrently. Must be
#'   a positive integer and is capped at `chains`. Values greater than one
#'   use independent multisession workers. The default follows the `mc.cores`
#'   option.
#' @param refresh Nonnegative integer controlling how often each chain reports
#'   sampling progress through `progressr`. Terminals with multiline cursor
#'   support and RStudio 2026.06.0 or later show one live bar per chain. Other
#'   consoles print separate per-chain iteration updates. Set to `0` to
#'   disable progress output and the completion summary.
#'
#' @return A `combo_fit` object containing the original and compiled model,
#'   retained Gibbs-state snapshots in `draws`, chain and iteration indices,
#'   sampling settings, inputs, and the Gibbs implementation used.
#'
#' @export
fit_combo <- function(
    model,
    data,
    cell_data = NULL,
    compound_data = NULL,
    engine = "gibbs",
    chains = 4L,
    iter_warmup = 1000L,
    iter_sampling = 1000L,
    thin = 1L,
    seed = NULL,
    control = list(),
    parallel_chains = getOption("mc.cores", 1L),
    refresh = max((iter_warmup + iter_sampling) %/% 10L, 1L)
) {
    if (length(engine) != 1L || !identical(engine, "gibbs")) {
        cli::cli_abort("engine must be \"gibbs\"")
    }
    if (!is_positive_integer(chains)) {
        cli::cli_abort("chains must be integer-valued and at least 1")
    }
    if (!is_nonnegative_integer(iter_warmup)) {
        cli::cli_abort("iter_warmup must be integer-valued and at least 0")
    }
    if (!is_positive_integer(thin)) {
        cli::cli_abort("thin must be integer-valued and at least 1")
    }
    if (!is_positive_integer(iter_sampling) || iter_sampling < thin) {
        cli::cli_abort("iter_sampling must be integer-valued and at least thin")
    }
    if (!is.null(seed) && (!is_positive_integer(seed))) {
        cli::cli_abort("seed must be NULL or one finite integer")
    }
    if (!is.list(control)) {
        cli::cli_abort("control must be a named list")
    }
    if (length(control)) {
        labels <- names(control)
        if (is.null(labels) || any(!nzchar(labels))) {
            cli::cli_abort("control must be a named list")
        }
        cli::cli_abort("Unused Gibbs control parameter(s): {.and {labels}}")
    }
    if (!is_positive_integer(parallel_chains)) {
        cli::cli_abort("parallel_chains must be one finite positive integer")
    }
    if (!is_nonnegative_integer(refresh)) {
        cli::cli_abort("refresh must be one finite nonnegative integer")
    }
    chains <- as.integer(chains)
    iter_warmup <- as.integer(iter_warmup)
    iter_sampling <- as.integer(iter_sampling)
    thin <- as.integer(thin)
    if (!is.null(seed)) seed <- as.integer(seed)
    parallel_chains <- min(as.integer(parallel_chains), chains)
    refresh <- as.integer(refresh)
    compiled <- compile_combo_model(
        model,
        data,
        cell_data = cell_data,
        compound_data = compound_data
    )
    retained_per_chain <- iter_sampling %/% thin
    chain_draws <- combo_run_chains(
        compiled = compiled,
        seed = seed,
        chains = chains,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        thin = thin,
        parallel_chains = parallel_chains,
        refresh = refresh
    )
    snapshots <- unlist(chain_draws, recursive = FALSE)
    chain_id <- rep(seq_len(chains), each = retained_per_chain)
    draw_id <- rep(seq_len(retained_per_chain), times = chains)
    iteration <- rep(
        iter_warmup + seq.int(thin, iter_sampling, by = thin),
        times = chains
    )
    structure(
        list(
            model = model,
            compiled = compiled,
            draws = snapshots,
            chain_id = chain_id,
            draw_id = draw_id,
            iteration = iteration,
            engine = engine,
            sampling = list(
                chains = chains,
                iter_warmup = iter_warmup,
                iter_sampling = iter_sampling,
                thin = thin,
                seed = seed,
                parallel_chains = parallel_chains
            ),
            control = control,
            input = list(
                data = data,
                cell_data = cell_data,
                compound_data = compound_data
            ),
            implementation = if (compiled$fast_iid) "gibbs_iid" else "gibbs_sparse"
        ),
        class = "combo_fit"
    )
}


# Basic fit methods -------------------------------------------------------

#' Inspect a fitted combination-response model
#'
#' `print()` reports sampling dimensions, model rank, modeled entities, maximum R-hat,
#' minimum bulk and tail effective sample sizes (ESS), and counts of unavailable
#' diagnostics across public posterior variables.
#' `summary()` adds a per-variable [posterior::summarize_draws()] table with
#' means, medians, standard deviations, MADs, 5% and 95% quantiles, R-hat,
#' bulk ESS, and tail ESS by default. Chain boundaries are preserved.
#'
#' Diagnostics may be unavailable for short or constant chains. Unavailable
#' values are retained in the summary and counted separately in `print()`.
#' Diagnostics for latent factors can also reflect their non-identifiability.
#' The summary table is printed with all rows and columns by default.
#' Use `print(summary(fit), n = 20)` to limit the displayed rows.
#'
#' @param x,object A `combo_fit` object.
#' @param components,variable,include Selection arguments passed to
#'   [posterior_draws()] by `summary()`.
#' @param ... For `summary()`, arguments passed to
#'   [posterior::summarize_draws()]. Supplying summary functions replaces its
#'   default measures. For `print.summary.combo_fit()`, table printing options
#'   passed to `print()` (such as `n` and `width`). Unused by `print.combo_fit()`.
#' @return `print()` returns its input invisibly. `summary()` returns a
#'   `summary.combo_fit` object retaining the fit metadata (including `rank`,
#'   the configured latent dimension) and RMSE and
#'   observation-precision quantiles, with a `posterior_summary` element
#'   containing the draws summary table.
#' @name combo_fit_methods
NULL

#' @rdname combo_fit_methods
#' @export
print.combo_fit <- function(x, ...) {
    diagnostics <- posterior::summarize_draws(
        posterior::as_draws_array(x),
        posterior::default_convergence_measures()
    )
    chain_steps <- x[["sampling"]][["iter_warmup"]] + x[["sampling"]][["iter_sampling"]]
    output <- cli::cli_format_method({
        cli::cli_text("<combo_fit>")
        cli::cli_text("engine: {x[['engine']]}")
        cli::cli_text("rank: {x[['model']][['rank']]}")
        cli::cli_text("chains: {x[['sampling']][['chains']]}")
        cli::cli_text("iterations: {chain_steps} per chain")
        cli::cli_text("warmup: {x[['sampling']][['iter_warmup']]}")
        cli::cli_text("sampling: {x[['sampling']][['iter_sampling']]}")
        cli::cli_text("thin: {x[['sampling']][['thin']]}")
        cli::cli_text("retained draws: {length(x[['draws']])}")
        cli::cli_text("observations: {length(x[['compiled']][['observed_rows']])}")
        cli::cli_text("cells: {length(x[['compiled']][['cells']])}")
        cli::cli_text("treatments: {nrow(x[['compiled']][['treatments']])}")
        cli::cli_text("Convergence diagnostics ({nrow(diagnostics)} public variables):")
        cli::cli_text("max R-hat: {combo_diagnostic_range(diagnostics$rhat, max)}")
        cli::cli_text("min bulk ESS: {combo_diagnostic_range(diagnostics$ess_bulk, min)}")
        cli::cli_text("min tail ESS: {combo_diagnostic_range(diagnostics$ess_tail, min)}")
    })
    writeLines(output)
    invisible(x)
}

#' @rdname combo_fit_methods
#' @export
summary.combo_fit <- function(
    object,
    ...,
    components = NULL,
    variable = NULL,
    include = c("public", "all")
) {
    posterior_summary <- posterior::summarize_draws(
        posterior::as_draws_array(
            object, components = components, variable = variable, include = include
        ),
        ...
    )
    rmse <- vapply(object$draws, "[[", numeric(1), "last_rmse")
    precision <- vapply(object$draws, "[[", numeric(1), "precision")
    probs <- c(0.05, 0.5, 0.95)
    result <- list(
        engine = object$engine,
        rank = object$model$rank,
        n_observed = length(object$compiled$observed_rows),
        n_cells = length(object$compiled$cells),
        n_treatments = nrow(object$compiled$treatments),
        n_draws = length(object$draws),
        posterior_summary = posterior_summary,
        rmse = stats::quantile(rmse, probs),
        observation_precision = stats::quantile(precision, probs)
    )
    class(result) <- "summary.combo_fit"
    result
}

#' @rdname combo_fit_methods
#' @export
print.summary.combo_fit <- function(x, ...) {
    output <- cli::cli_format_method({
        cli::cli_text("Combination-model Gibbs fit")
        cli::cli_text("engine: {x[['engine']]}")
        cli::cli_text("rank: {x[['rank']]}")
        cli::cli_text("observations: {x[['n_observed']]}")
        cli::cli_text("cells: {x[['n_cells']]}")
        cli::cli_text("treatments: {x[['n_treatments']]}")
        cli::cli_text("draws: {x[['n_draws']]}")
        cli::cli_text(
            "RMSE (5%, 50%, 95%): {format(x[['rmse']], digits = 4)}"
        )
    })
    writeLines(output)
    options <- list(...)
    if (is.null(options$n)) options$n <- Inf
    if (is.null(options$width)) options$width <- Inf
    do.call(print, c(list(x = x$posterior_summary), options))
    invisible(x)
}

# Keep infinite diagnostics visible; omit only unavailable values from extrema.
combo_diagnostic_range <- function(x, fun) {
    available <- !is.na(x)
    value <- if (any(available)) format(fun(x[available]), digits = 4) else "NA"
    paste0(value, " (", sum(!available), " unavailable)")
}
