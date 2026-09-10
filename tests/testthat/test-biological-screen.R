test_that("the biological screen contains distinct triplicate measurements", {
    screen <- read.csv(system.file(
        "extdata", "biologically_informed_models.csv",
        package = "batchieR", mustWork = TRUE
    ))
    groups <- split(screen, screen$configuration)
    expect_true(all(vapply(groups, nrow, integer(1)) == 3L))
    expect_true(all(vapply(groups, function(group) {
        identical(sort(group$replicate), 1:3)
    }, logical(1))))
    expect_true(all(vapply(groups, function(group) {
        length(unique(group$modeled_response)) == 3L &&
            length(unique(group$true_mean)) == 1L
    }, logical(1))))
    expect_true(all(is.finite(screen$modeled_response)))
    expect_equal(screen$viability, plogis(screen$modeled_response), tolerance = 1e-12)
})

test_that("biological vignette masking holds out all replicates together", {
    vignette <- test_path("..", "..", "vignettes", "biologically_informed_models.qmd")
    skip_if_not(file.exists(vignette), "Vignette source is not installed")
    lines <- readLines(vignette)
    start <- match("#| label: response-mask", lines) + 1L
    end <- start + which(lines[start:length(lines)] == "```")[1L] - 2L
    context <- new.env(parent = environment())
    context$screen <- read.csv(system.file(
        "extdata", "biologically_informed_models.csv",
        package = "batchieR", mustWork = TRUE
    ))
    eval(parse(text = lines[start:end]), envir = context)
    screen <- context$screen
    observed <- !is.na(screen$response)
    expect_true(all(table(screen$configuration[observed]) == 3L))
    expect_true(all(is.na(screen$response[screen$role != "training"])))
    expect_true(all(vapply(split(observed, screen$configuration), function(group) {
        length(unique(group)) == 1L
    }, logical(1))))
    for (plate in names(context$prediction_indices)) {
        selected <- context$prediction_indices[[plate]]
        expect_length(selected, 4L)
        expect_true(all(screen$replicate[selected] == 1L))
        expect_equal(sum(screen$role == plate), 12L)
        expect_equal(nrow(context$candidate_data[[plate]]), 8L)
    }
})
