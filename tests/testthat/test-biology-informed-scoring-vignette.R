test_that("prior-information vignette has a concept-first workflow", {
    source_path <- test_path(
        "..", "..", "vignettes", "biology-informed-scoring.Rmd"
    )
    if (!file.exists(source_path)) {
        source_path <- system.file(
            "doc",
            "biology-informed-scoring.Rmd",
            package = "batchieR",
            mustWork = TRUE
        )
    }
    vignette <- paste(readLines(source_path, warn = FALSE), collapse = "\n")

    expect_match(vignette, "Incorporating biological prior information", fixed = TRUE)
    expect_match(vignette, "## From biological knowledge to a prior", fixed = TRUE)
    expect_match(vignette, "\\theta = X\\beta + u", fixed = TRUE)
    expect_match(vignette, "The five configurable variables", fixed = TRUE)
    expect_match(vignette, "`combo_gaussian_component()` turns these choices", fixed = TRUE)
    expect_match(vignette, "limits all five configurable terms", fixed = TRUE)
    expect_match(vignette, "latent coordinates an intrinsic interpretation", fixed = TRUE)
    expect_match(vignette, "example-hierarchies", fixed = TRUE)
    expect_match(vignette, "## Four ways to inform the model", fixed = TRUE)
    expect_match(vignette, "### 1. Shrink unused latent dimensions", fixed = TRUE)
    expect_match(vignette, "### 2. Put measured features in the prior mean", fixed = TRUE)
    expect_match(vignette, "None is a uniquely implied", fixed = TRUE)
    expect_match(vignette, "isolates the effect of replacing", fixed = TRUE)
    expect_match(vignette, "Offsets are the natural place", fixed = TRUE)
    expect_match(vignette, "two fields relevant to feature", fixed = TRUE)
    expect_match(vignette, "matched to model", fixed = TRUE)
    expect_match(vignette, "Compound annotations are then attached", fixed = TRUE)
    expect_match(vignette, "modifies only `cell_offset` and", fixed = TRUE)
    expect_match(vignette, "specifications remain unchanged", fixed = TRUE)
    expect_match(vignette, "### 3. Partially pool doses within compounds", fixed = TRUE)
    expect_match(vignette, "### 4. Relate entities through a graph or tree", fixed = TRUE)
    expect_match(vignette, "dependence structure, not", fixed = TRUE)
    expect_match(vignette, "`cell_factors`, `treatment_main_factors`", fixed = TRUE)
    expect_match(vignette, "feature-informed offsets", fixed = TRUE)
    expect_match(vignette, "offsets do not use either kernel", fixed = TRUE)
    expect_match(vignette, "#### Choosing a representation", fixed = TRUE)
    expect_match(vignette, "#### A hierarchy with `tree()`", fixed = TRUE)
    expect_match(vignette, "#### Pairwise relationships with a GMRF", fixed = TRUE)
    expect_match(vignette, "#### Applying a structure to the factor components", fixed = TRUE)
    expect_match(vignette, "#### Technical requirements and scaling", fixed = TRUE)
    expect_match(vignette, "does not currently construct this operator", fixed = TRUE)
    expect_match(vignette, "exactly symmetric", fixed = TRUE)
    expect_match(vignette, "finite positive ridge", fixed = TRUE)
    expect_match(vignette, "`kernlab`", fixed = TRUE)
    expect_match(vignette, "Thresholding similarities to construct a graph", fixed = TRUE)
    expect_match(vignette, "x^\\top Lx", fixed = TRUE)
    expect_match(vignette, "nonzero pattern, not merely", fixed = TRUE)
    expect_match(vignette, "k-nearest-", fixed = TRUE)
    expect_match(vignette, "shared-nearest-neighbor (SNN)", fixed = TRUE)
    expect_match(vignette, "should be symmetrized", fixed = TRUE)

    expect_match(vignette, "screen-layout, include=FALSE", fixed = TRUE)
    expect_match(vignette, "generate-responses, include=FALSE", fixed = TRUE)
    expect_match(vignette, "response-mask, include=FALSE", fixed = TRUE)
    expect_match(vignette, "fit-models, include=FALSE", fixed = TRUE)
    expect_match(vignette, "structure-kernels, echo=FALSE", fixed = TRUE)

    expect_match(vignette, "default_model <- combo_model", fixed = TRUE)
    expect_match(vignette, "feature_model <- combo_model", fixed = TRUE)
    expect_match(vignette, "nested_dose_model <- combo_model", fixed = TRUE)
    expect_match(vignette, "gmrf_model <- combo_model", fixed = TRUE)
    expect_match(vignette, "tree_model <- combo_model", fixed = TRUE)
    expect_match(vignette, "~ 0 + pathway_activity + growth_rate", fixed = TRUE)
    expect_match(vignette, "~ 0 + target_signature + physchem_pc", fixed = TRUE)
    expect_match(vignette, "nested(precision = 4)", fixed = TRUE)
    expect_match(vignette, "cell_gmrf <- gmrf", fixed = TRUE)
    expect_match(vignette, "cell_tree <- tree", fixed = TRUE)
    expect_match(vignette, "posterior_predict", fixed = TRUE)
    expect_match(vignette, "score_plate_pdbal", fixed = TRUE)
    expect_match(vignette, "max_triplets = 1000", fixed = TRUE)

    concept_position <- regexpr("## From biological knowledge", vignette, fixed = TRUE)
    score_position <- regexpr("### Plate selection", vignette, fixed = TRUE)
    expect_gt(score_position, concept_position)
    expect_match(vignette, "set.seed(20260827)", fixed = TRUE)
    expect_match(vignette, "set.seed(20260900)", fixed = TRUE)
    seed_locations <- gregexpr("set.seed(", vignette, fixed = TRUE)[[1L]]
    expect_length(seed_locations[seed_locations > 0L], 2L)
})

test_that("introduction links the biology-informed workflow", {
    source_path <- test_path("..", "..", "vignettes", "introduction.Rmd")
    if (!file.exists(source_path)) {
        source_path <- system.file(
            "doc",
            "introduction.Rmd",
            package = "batchieR",
            mustWork = TRUE
        )
    }
    introduction <- paste(readLines(source_path, warn = FALSE), collapse = "\n")

    expect_match(
        introduction,
        "vignette(\"biology-informed-scoring\"",
        fixed = TRUE
    )
})
