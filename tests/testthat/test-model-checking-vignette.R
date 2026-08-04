test_that("model-checking vignette covers the diagnostic workflow", {
    source_path <- test_path("..", "..", "vignettes", "model-checking.Rmd")
    if (!file.exists(source_path)) {
        source_path <- system.file(
            "doc",
            "model-checking.Rmd",
            package = "batchieR",
            mustWork = TRUE
        )
    }
    vignette <- paste(readLines(source_path, warn = FALSE), collapse = "\n")

    expect_match(vignette, "%\\VignetteEngine{knitr::rmarkdown}", fixed = TRUE)
    expect_match(vignette, "full_fit", fixed = TRUE)
    expect_match(vignette, "additive_fit", fixed = TRUE)
    expect_match(vignette, "type = \"loo_pit\"", fixed = TRUE)
    expect_match(vignette, "type = \"error\"", fixed = TRUE)
    expect_match(vignette, "type = \"stat_2d\"", fixed = TRUE)
    expect_match(vignette, "group = \"experiment_type\"", fixed = TRUE)
    expect_match(vignette, "From patterns to action", fixed = TRUE)
})
