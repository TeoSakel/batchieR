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
#' @param model A combination-model specification created by `combo_model()`.
#' @param data A data frame of combination-screen observations. See Details for
#'   the required columns.
#' @param cell_data An optional data frame of cell-level covariates, containing a
#'   unique, nonmissing `cell` key and any covariates used by cell-component
#'   mean formulas.
#' @param compound_data An optional data frame of compound-level covariates,
#'   containing a unique, nonmissing `drug` key and any covariates used by
#'   treatment-component mean formulas.
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
    control = list()
) {
    if (length(engine) != 1L || !identical(engine, "gibbs")) {
        cli::cli_abort("engine must be \"gibbs\"")
    }
    run_values <- c(chains, iter_warmup, iter_sampling, thin)
    if (anyNA(run_values) || any(!is.finite(run_values)) ||
            chains < 1 || iter_warmup < 0 || iter_sampling < thin || thin < 1 ||
            any(run_values != as.integer(run_values))) {
        cli::cli_abort(
            "chains, iter_warmup, iter_sampling, and thin must be integer-valued and within their allowed ranges"
        )
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
    chains <- as.integer(chains)
    iter_warmup <- as.integer(iter_warmup)
    iter_sampling <- as.integer(iter_sampling)
    thin <- as.integer(thin)
    if (!is.null(seed)) seed <- as.integer(seed)
    compiled <- compile_combo_model(
        model,
        data,
        cell_data = cell_data,
        compound_data = compound_data
    )
    set.seed(seed)
    retained_per_chain <- iter_sampling %/% thin
    total <- chains * retained_per_chain
    snapshots <- vector("list", total)
    chain_id <- integer(total)
    draw_id <- integer(total)
    iteration <- integer(total)
    output <- 0L
    # TODO: Consider parallelizing chains in the future.
    for (chain in seq_len(chains)) {
        state <- init_gibbs_state(compiled)
        for (step in seq_len(iter_warmup)) state <- gibbs_step(state)
        for (step in seq_len(iter_sampling)) {
            state <- gibbs_step(state)
            if (step %% thin == 0L) {
                output <- output + 1L
                snapshots[[output]] <- gibbs_snapshot(state)
                chain_id[output] <- chain
                draw_id[output] <- step %/% thin
                iteration[output] <- iter_warmup + step
            }
        }
    }
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
                seed = seed
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
#' `print()` reports sampling dimensions and modeled entities. `summary()`
#' returns posterior quantiles of sampler RMSE and observation precision.
#'
#' @param x,object A `combo_fit` object.
#' @param ... Reserved for future methods.
#' @return `print()` returns its input invisibly. `summary()` returns a
#'   `summary.combo_fit` object.
#' @name combo_fit_methods
NULL

#' @rdname combo_fit_methods
#' @export
print.combo_fit <- function(x, ...) {
    output <- cli::cli_format_method({
    cli::cli_text("<combo_fit>")
    cli::cli_text("engine: {x[['engine']]}")
    cli::cli_text("chains: {x[['sampling']][['chains']]}")
    cli::cli_text(
        "iterations: {x[['sampling']][['iter_warmup']] + x[['sampling']][['iter_sampling']]} per chain"
    )
    cli::cli_text("warmup: {x[['sampling']][['iter_warmup']]}")
    cli::cli_text("sampling: {x[['sampling']][['iter_sampling']]}")
    cli::cli_text("thin: {x[['sampling']][['thin']]}")
    cli::cli_text("retained draws: {length(x[['draws']])}")
    cli::cli_text("observations: {length(x[['compiled']][['observed_rows']])}")
    cli::cli_text("cells: {length(x[['compiled']][['cells']])}")
    cli::cli_text("treatments: {nrow(x[['compiled']][['treatments']])}")
    })
    writeLines(output)
    invisible(x)
}

#' @rdname combo_fit_methods
#' @export
summary.combo_fit <- function(object, ...) {
    rmse <- vapply(
        object$draws,
        function(draw) draw$last_rmse,
        numeric(1)
    )
    precision <- vapply(
        object$draws,
        function(draw) draw$precision,
        numeric(1)
    )
    result <- list(
        engine = object$engine,
        n_observed = length(object$compiled$observed_rows),
        n_cells = length(object$compiled$cells),
        n_treatments = nrow(object$compiled$treatments),
        n_draws = length(object$draws),
        rmse = stats::quantile(rmse, c(0.05, 0.5, 0.95)),
        observation_precision = stats::quantile(
            precision,
            c(0.05, 0.5, 0.95)
        )
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
    cli::cli_text("observations: {x[['n_observed']]}")
    cli::cli_text("cells: {x[['n_cells']]}")
    cli::cli_text("treatments: {x[['n_treatments']]}")
    cli::cli_text("draws: {x[['n_draws']]}")
    cli::cli_text(
        "RMSE (5%, 50%, 95%): {format(x[['rmse']], digits = 4)}"
    )
    })
    writeLines(output)
    invisible(x)
}
