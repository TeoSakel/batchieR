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
    expect_match(vignette, "vignette(\"model-checking\"", fixed = TRUE)
    expect_match(vignette, "score_plate_pdbal", fixed = TRUE)
    expect_match(vignette, "confirmatory_plate_6", fixed = TRUE)
    expect_match(vignette, "discovery_plate_7", fixed = TRUE)
})
