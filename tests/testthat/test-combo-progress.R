test_that("multiline detection recognizes supported RStudio consoles", {
    supports <- function(version, dynamic = TRUE) {
        combo_progress_terminal(
            tty = FALSE, ansi = FALSE, dynamic = dynamic,
            rstudio_version = version
        )
    }
    expect_false(supports(NULL))
    expect_false(supports("2026.05.9"))
    expect_false(supports("unknown"))
    expect_true(supports("2026.06.0"))
    expect_true(supports("2026.06.0+548"))
    expect_true(supports(package_version("2026.08.2")))
    expect_false(supports("2026.08.2", dynamic = FALSE))
    expect_true(combo_progress_terminal(
        tty = TRUE, ansi = TRUE, dynamic = TRUE, rstudio_version = NULL
    ))
    expect_false(combo_progress_terminal(
        tty = TRUE, ansi = FALSE, dynamic = TRUE, rstudio_version = NULL
    ))
})

test_that("completion summaries use the measured mean and total", {
    output <- capture.output(combo_sampling_summary(c(2, 4), 10), type = "message")
    expect_identical(output, c(
        "All 2 chains finished successfully.",
        "Mean chain execution time: 3.0 seconds.",
        "Total execution time: 10.0 seconds."
    ))
    output <- capture.output(combo_sampling_summary(2, 3), type = "message")
    expect_identical(output[[1]], "All 1 chain finished successfully.")
})

test_that("chain timing measures elapsed execution without changing draws", {
    ticks <- c(10, 12.5)
    local_mocked_bindings(
        combo_clock = function() {
            value <- ticks[[1]]
            ticks <<- ticks[-1]
            value
        },
        combo_run_chain = function(...) list(draw = 42)
    )
    result <- combo_chain_task(1L, NULL, 2L, 4L, 1L, 0L, NULL)
    expect_true(result$ok)
    expect_identical(result$value, list(draw = 42))
    expect_identical(result$elapsed, 2.5)
})

test_that("failed fits never print a successful completion summary", {
    for (terminal in c(FALSE, TRUE)) {
        local_mocked_bindings(combo_progress_terminal = function() terminal)
        output <- capture.output(
            expect_error(
                combo_run_chains(
                    list(), 42L, 1L, 1L, 1L, 1L, 1L, 1L
                ),
                "Chain 1 failed"
            ),
            type = "message"
        )
        expect_false(any(grepl("finished successfully", output)))
        expect_false(any(grepl("execution time", output)))
    }
})

test_that("total timing spans setup and result collection", {
    ticks <- c(5, 15)
    local_mocked_bindings(
        combo_progress_terminal = function() FALSE,
        combo_clock = function() {
            value <- ticks[[1]]
            ticks <<- ticks[-1]
            value
        },
        combo_chain_worker_environment = function() {
            list(combo_chain_task = function(chain, ...) {
                list(ok = TRUE, value = list(chain = chain), elapsed = chain * 2)
            })
        }
    )
    output <- capture.output(
        draws <- combo_run_chains(NULL, 42L, 2L, 2L, 4L, 1L, 1L, 1L),
        type = "message"
    )
    expect_identical(tail(output, 3L), c(
        "All 2 chains finished successfully.",
        "Mean chain execution time: 3.0 seconds.",
        "Total execution time: 10.0 seconds."
    ))
    expect_identical(draws, list(list(chain = 1L), list(chain = 2L)))
})
