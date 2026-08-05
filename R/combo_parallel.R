# Parallel chain execution for combination models.

combo_run_chain <- function(
    compiled,
    iter_warmup,
    iter_sampling,
    thin,
    report = NULL
) {
    state <- init_gibbs_state(compiled)
    if (!is.null(report)) report(0L, "warmup")
    for (step in seq_len(iter_warmup)) {
        state <- gibbs_step(state)
        if (!is.null(report)) report(step, "warmup")
    }
    if (!is.null(report)) report(iter_warmup, "sampling")
    snapshots <- vector("list", iter_sampling %/% thin)
    for (step in seq_len(iter_sampling)) {
        state <- gibbs_step(state)
        if (step %% thin == 0L) {
            snapshots[[step %/% thin]] <- gibbs_snapshot(state)
        }
        if (!is.null(report)) report(iter_warmup + step, "sampling")
    }
    snapshots
}

combo_chain_task <- function(
    chain,
    compiled,
    iter_warmup,
    iter_sampling,
    thin,
    refresh,
    progress
) {
    report <- combo_chain_reporter(
        refresh = refresh,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        chain = chain,
        progress = progress
    )
    tryCatch(
        list(ok = TRUE, value = combo_run_chain(
            compiled,
            iter_warmup,
            iter_sampling,
            thin,
            report
        )),
        error = function(error) {
            list(ok = FALSE, message = conditionMessage(error))
        }
    )
}

combo_chain_reporter <- function(
    refresh,
    iter_warmup,
    iter_sampling,
    chain,
    progress = NULL,
    update = NULL
) {
    if (refresh == 0L || (is.null(progress) && is.null(update))) {
        return(NULL)
    }
    total <- iter_warmup + iter_sampling
    last <- -1L
    last_phase <- NULL
    function(current, phase) {
        boundary <- current == 0L || current == iter_warmup || current == total
        phase_change <- !identical(phase, last_phase)
        if (!boundary && !phase_change && current - last < refresh) {
            return(invisible(NULL))
        }
        amount <- max(current - max(last, 0L), 0L)
        last <<- current
        last_phase <<- phase
        info <- list(
            chain = chain,
            current = current,
            phase = phase,
            amount = amount
        )
        if (!is.null(update)) update(info)
        if (!is.null(progress)) {
            progress(
                amount = amount,
                message = paste0("Chain ", chain, ": ", phase)
            )
        }
        invisible(NULL)
    }
}

combo_run_chains <- function(
    compiled,
    seed,
    chains,
    iter_warmup,
    iter_sampling,
    thin,
    parallel_chains,
    refresh
) {
    strategy <- if (parallel_chains == 1L) {
        future::sequential
    } else {
        future::tweak(
            future::multisession,
            workers = parallel_chains
        )
    }
    previous_plan <- future::plan(strategy)
    on.exit(future::plan(previous_plan), add = TRUE)

    steps_per_chain <- iter_warmup + iter_sampling
    chain_task <- combo_chain_task
    environment(chain_task) <- list2env(
        list(
            combo_chain_reporter = combo_chain_reporter,
            combo_run_chain = combo_run_chain
        ),
        parent = environment(chain_task)
    )
    results <- progressr::with_progress(
        {
            progress <- lapply(seq_len(chains), function(chain) {
                progressr::progressor(
                    steps = steps_per_chain,
                    message = paste0("Chain ", chain, ": queued"),
                    label = paste0("Chain ", chain),
                    enable = refresh > 0L
                )
            })
            future.apply::future_lapply(
                seq_len(chains),
                function(chain) {
                    chain_task(
                        chain = chain,
                        compiled = compiled,
                        iter_warmup = iter_warmup,
                        iter_sampling = iter_sampling,
                        thin = thin,
                        refresh = refresh,
                        progress = progress[[chain]]
                    )
                },
                future.seed = if (is.null(seed)) TRUE else seed,
                future.scheduling = Inf,
                future.label = "batchieR-chain-%d"
            )
        },
        enable = refresh > 0L,
        handlers = progressr::handler_cli(
            format = "{message} {cli::pb_bar} {cli::pb_percent}"
        )
    )
    for (chain in seq_along(results)) {
        result <- results[[chain]]
        if (!isTRUE(result[["ok"]])) {
            cli::cli_abort(c("Chain {chain} failed.", "x" = result[["message"]]))
        }
        results[[chain]] <- result[["value"]]
    }
    results
}
