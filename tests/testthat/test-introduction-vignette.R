test_that("introduction vignette has the portable workflow contract", {
    source_path <- test_path("..", "..", "vignettes", "introduction.Rmd")
    if (!file.exists(source_path)) {
        source_path <- system.file(
            "doc",
            "introduction.Rmd",
            package = "batchieR",
            mustWork = TRUE
        )
    }
    vignette <- paste(readLines(source_path, warn = FALSE), collapse = "\n")

    expect_match(vignette, "rmarkdown::html_vignette", fixed = TRUE)
    expect_match(vignette, "%\\VignetteEngine{knitr::rmarkdown}", fixed = TRUE)
    expect_match(vignette, "only a Gaussian response", fixed = TRUE)
    expect_match(vignette, "bayesplot::pp_check", fixed = TRUE)
    expect_match(vignette, "vignette(\"model_checking\"", fixed = TRUE)
    expect_match(vignette, "future::availableCores()", fixed = TRUE)
    expect_match(vignette, "parallel_chains = fit_parallel_chains", fixed = TRUE)
    expect_match(vignette, "refresh = if \\(interactive\\(\\)\\) 10L? else 0L?")
    expect_match(vignette, "score_plate_pdbal", fixed = TRUE)
    expect_match(vignette, '"new_cell"', fixed = TRUE)
    expect_match(vignette, '"new_drug_pairs"', fixed = TRUE)
})
