horseshoe_vignette_fixture <- function(name) {
    readRDS(system.file("extdata", paste0("predictive-mixture-", name, "-v2.rds"),
                        package = "batchieR", mustWork = TRUE))
}

test_that("ten-chain fixture preserves draws and original run provenance", {
    inputs <- horseshoe_vignette_fixture("inputs")
    fit <- horseshoe_vignette_fixture("fit")
    expect_identical(inputs$model, inputs$generating_model)
    expect_identical(fit$model, inputs$model)
    expect_identical(fit$model$mean$type, "fixed")
    expect_equal(fit$model$mean$value, 0)
    expect_identical(sort(unique(fit$chain_id)), 1:10)
    expect_equal(as.integer(table(fit$chain_id)), rep(800L, 10))
    expect_identical(fit$draw_id, rep(1:800, 10))
    expect_identical(fit$iteration, rep(2000L + 40L * (1:800), 10))
    expect_identical(fit$input$data, inputs$training)
    expect_identical(fit$sampling$chains, 10L)
    expect_identical(fit$provenance$source_sampling$chains, 10L)
    expect_identical(fit$provenance$selected_chains, 1:10)
    expect_true(all(vapply(c("treatment_offset", "treatment_main_factors",
                            "treatment_interaction_factors"), function(component) {
        fit$model$components[[component]]$shrinkage$type == "horseshoe"
    }, logical(1))))
    diagnostics <- summary(fit, select = "sigma")$posterior_summary
    expect_true(all(is.finite(unlist(diagnostics[c("rhat", "ess_bulk", "ess_tail")]))))
})

test_that("new fixture responses reproduce the single unscreened prior draw", {
    inputs <- horseshoe_vignette_fixture("inputs")
    replicated <- do.call(rbind, replicate(5, inputs$design, simplify = FALSE))
    replicated$replicate <- rep(1:5, each = nrow(inputs$design))
    response <- prior_predict(inputs$generating_model, replicated, draws = 1,
                              type = "response", seed = inputs$settings$generation_seed)
    expect_equal(as.numeric(response[1, ]),
                 c(inputs$training$response, inputs$holdout$response), tolerance = 1e-12)
})

test_that("ten-chain clustering reflects its recalculated noise floor", {
    inputs <- horseshoe_vignette_fixture("inputs")
    fit <- horseshoe_vignette_fixture("fit")
    projection <- chain_projections(fit, reference_grid = inputs$design,
                                    response_transform = plogis)
    clustering <- batchieR:::combo_mixture_cluster(projection$mean, projection$block_mean)
    expect_identical(unname(clustering$basin), c("B1", "B1", "B2", rep("B1", 7)))
    expect_gte(clustering$gap_ratio, 3)
})
