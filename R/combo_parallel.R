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
    started <- combo_clock()
    report <- combo_chain_reporter(
        refresh = refresh,
        iter_warmup = iter_warmup,
        iter_sampling = iter_sampling,
        chain = chain,
        progress = progress
    )
    tryCatch(
        {
            draws <- combo_run_chain(
                compiled, iter_warmup, iter_sampling, thin, report
            )
            list(ok = TRUE, value = draws, elapsed = combo_clock() - started)
        },
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
                message = paste0("Chain ", chain, ": ", phase),
                chain = chain,
                current = current,
                phase = phase
            )
        }
        invisible(NULL)
    }
}

combo_clock <- function() proc.time()[["elapsed"]]

combo_rstudio_version <- function() {
    version <- get0("RStudio.Version", envir = globalenv(), mode = "function")
    if (is.null(version)) return(NULL)
    tryCatch(version()$version, error = function(error) NULL)
}

combo_progress_terminal <- function(
    tty = isatty(stderr()),
    ansi = cli::is_ansi_tty(stderr()),
    dynamic = cli::is_dynamic_tty(stderr()),
    rstudio_version = combo_rstudio_version()
) {
    if (tty && ansi) return(TRUE)
    if (!dynamic || is.null(rstudio_version)) return(FALSE)
    # RStudio 2026.06 supports the multiline cursor movement used below.
    tryCatch(
        isTRUE(package_version(sub("\\+.*$", "", as.character(rstudio_version))) >= "2026.6.0"),
        error = function(error) FALSE
    )
}

combo_sampling_summary <- function(chain_elapsed, total_elapsed) {
    chains <- length(chain_elapsed)
    cli::cli_verbatim(sprintf(
        "All %d %s finished successfully.",
        chains, if (chains == 1L) "chain" else "chains"
    ))
    cli::cli_verbatim(sprintf(
        "Mean chain execution time: %.1f seconds.", mean(chain_elapsed)
    ))
    cli::cli_verbatim(sprintf(
        "Total execution time: %.1f seconds.", total_elapsed
    ))
}

combo_progress_handler <- function(chains, total, enable) {
    terminal <- combo_progress_terminal()
    current <- integer(chains)
    phase <- rep("queued", chains)
    visible <- FALSE

    render <- function() {
        if (!terminal) return()
        width <- max(1L, cli::console_width() - 1L)
        labels <- sprintf("Chain %d: %-8s ", seq_len(chains), phase)
        bar_width <- max(1L, min(30L, width - max(nchar(labels)) - 7L))
        filled <- floor(bar_width * current / total)
        lines <- paste0(
            labels, "[", strrep("=", filled),
            strrep("-", bar_width - filled), "] ",
            sprintf("%3.0f%%", floor(100 * current / total))
        )
        # Leave a spare column so a terminal cannot wrap into another bar's row.
        lines <- cli::ansi_strtrim(lines, width = width)
        if (visible) cat(sprintf("\033[%dA", chains), file = stderr())
        cat(paste0("\r\033[2K", lines, "\n"), sep = "", file = stderr())
        flush.console()
        visible <<- TRUE
    }
    hide <- function(...) {
        if (!visible) return()
        cat(sprintf("\033[%dA", chains),
            strrep("\r\033[2K\n", chains),
            sprintf("\033[%dA\r", chains), sep = "", file = stderr())
        visible <<- FALSE
    }
    progressr::make_progression_handler(
        "combo",
        enable = enable,
        # Each chain already applies the refresh cadence before sending events.
        interval = 0,
        times = Inf,
        reporter = list(
            reset = function(...) {
                current <<- integer(chains)
                phase <<- rep("queued", chains)
                visible <<- FALSE
            },
            initiate = function(...) render(),
            update = function(progression, ...) {
                chain <- progression$chain
                if (is.null(chain)) return()
                if (current[[chain]] == progression$current &&
                    phase[[chain]] == progression$phase) return()
                current[[chain]] <<- progression$current
                phase[[chain]] <<- progression$phase
                if (terminal) {
                    render()
                } else {
                    cli::cli_verbatim(sprintf(
                        "Chain %d Iteration: %d / %d [%3.0f%%] (%s)",
                        chain, current[[chain]], total,
                        floor(100 * current[[chain]] / total), phase[[chain]]
                    ))
                }
            },
            hide = hide,
            unhide = function(...) render(),
            finish = function(...) {
                # Keep actual counts, including when sampling exits early.
                if (!visible) render()
            }
        )
    )
}

combo_chain_worker_environment <- function() {
    namespace <- environment(combo_run_chains)
    if (!exists(".__DEVTOOLS__", envir = namespace, inherits = FALSE)) {
        return(namespace)
    }
    worker <- new.env(parent = baseenv())
    names <- ls(namespace, all.names = TRUE)
    names <- names[vapply(names, function(name) {
        value <- get(name, envir = namespace, inherits = FALSE)
        is.function(value) && identical(environment(value), namespace)
    }, logical(1))]
    # Clone package-owned closures into an ordinary environment. This allows
    # future to serialize development versions loaded with pkgload; namespace
    # closures would otherwise be treated as coming from an installed package.
    for (name in names) {
        assign(name, get(name, envir = namespace, inherits = FALSE), envir = worker)
    }
    for (name in names) {
        value <- get(name, envir = worker, inherits = FALSE)
        environment(value) <- worker
        assign(name, value, envir = worker)
    }
    worker
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
    started <- combo_clock()
    strategy <- if (parallel_chains == 1L) {
        future::sequential
    } else {
        future::tweak(
            future::multisession,
            workers = parallel_chains
        )
    }
    # Register restoration before starting a backend. Backend construction can
    # itself fail (for example, when local sockets are unavailable), and must
    # not leave a partially selected multisession plan behind.
    previous_plan <- future::plan()
    on.exit(future::plan(previous_plan), add = TRUE)
    future::plan(strategy)

    steps_per_chain <- iter_warmup + iter_sampling
    worker_environment <- combo_chain_worker_environment()
    chain_task <- worker_environment$combo_chain_task
    results <- progressr::with_progress(
        {
            # Keep one progressor alive for the entire fit. Progressors created
            # inside lapply() finish when each callback exits, before sampling.
            progress <- progressr::progressor(
                steps = chains * steps_per_chain,
                enable = refresh > 0L
            )
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
                        progress = progress
                    )
                },
                future.seed = if (is.null(seed)) TRUE else seed,
                future.scheduling = Inf,
                future.label = "batchieR-chain-%d"
            )
        },
        enable = refresh > 0L,
        # Do not buffer progression events themselves: they would cause a
        # hide/redraw cycle on every update.
        delay_conditions = c("message", "warning"),
        handlers = combo_progress_handler(
            chains = chains, total = steps_per_chain, enable = refresh > 0L
        )
    )
    chain_elapsed <- numeric(chains)
    for (chain in seq_along(results)) {
        result <- results[[chain]]
        if (!isTRUE(result[["ok"]])) {
            cli::cli_abort(c("Chain {chain} failed.", "x" = result[["message"]]))
        }
        chain_elapsed[[chain]] <- result[["elapsed"]]
        results[[chain]] <- result[["value"]]
    }
    if (refresh > 0L) {
        combo_sampling_summary(chain_elapsed, combo_clock() - started)
    }
    results
}
