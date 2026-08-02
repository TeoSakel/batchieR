named_sparse <- function(values, node_names) {
    Matrix::Matrix(
        values,
        sparse = TRUE,
        dimnames = list(node_names, node_names)
    )
}

compile_structure_spec <- function(spec, entity_names, label = "Cell") {
    construction <- construct_structure(spec, entity_names, label)
    finalize_structure(construction, label)
}

test_that("structure constructors use the generic structural prior class", {
    operator <- named_sparse(diag(2), c("a", "b"))
    hierarchy <- data.frame(node = "a", parent = NA_character_)

    expect_s3_class(iid(), "structural_prior")
    expect_s3_class(precision(operator), "structural_prior")
    expect_s3_class(gmrf(operator), "structural_prior")
    expect_s3_class(tree(node, parent, hierarchy), "structural_prior")
})

test_that("precision preserves a supplied positive-definite matrix", {
    Q <- named_sparse(matrix(c(2, 0.5, 0.5, 1), 2), c("b", "a"))

    compiled <- compile_structure_spec(precision(Q), c("a", "b"))

    expect_identical(compiled$type, "precision")
    expect_s4_class(compiled$Q, "dgCMatrix")
    expect_equal(as.matrix(compiled$Q), as.matrix(Q[c("a", "b"), c("a", "b")]))
    expect_gt(compiled$Q[1, 2], 0)

    component <- combo_gaussian_component(
        structure = precision(Q),
        shrinkage = gamma_precision()
    )
    expect_no_error(validate_component(component, "cell_offset", "offset"))
})

test_that("precision rejects invalid finished matrices", {
    node_names <- c("a", "b")
    valid <- named_sparse(diag(2), node_names)
    dense <- as.matrix(valid)
    asymmetric <- valid
    asymmetric[1, 2] <- 1e-12
    singular <- named_sparse(matrix(c(1, -1, -1, 1), 2), node_names)
    indefinite <- named_sparse(matrix(c(1, 2, 2, 1), 2), node_names)
    nonfinite <- valid
    nonfinite[1, 1] <- Inf
    unnamed <- valid
    dimnames(unnamed) <- NULL

    expect_error(compile_structure_spec(precision(dense), node_names), "sparse")
    expect_error(compile_structure_spec(precision(asymmetric), node_names), "symmetric")
    expect_error(compile_structure_spec(precision(singular), node_names), "positive definite")
    expect_error(compile_structure_spec(precision(indefinite), node_names), "positive definite")
    expect_error(compile_structure_spec(precision(nonfinite), node_names), "finite")
    expect_error(compile_structure_spec(precision(unnamed), node_names), "row and column")
    expect_error(compile_structure_spec(precision(valid), c("a", "c")), "exactly match")
})

test_that("gmrf adds only its diagonal ridge", {
    node_names <- c("a", "b")
    laplacian <- named_sparse(matrix(c(1, -1, -1, 1), 2), node_names)

    default <- compile_structure_spec(gmrf(laplacian), node_names)
    custom <- compile_structure_spec(gmrf(laplacian, ridge = 0.25), node_names)

    expect_identical(default$type, "gmrf")
    expect_equal(as.matrix(default$Q), as.matrix(laplacian + Matrix::Diagonal(2)))
    expect_equal(
        as.matrix(custom$Q),
        as.matrix(laplacian + Matrix::Diagonal(2, x = 0.25))
    )
    expect_error(compile_structure_spec(precision(laplacian), node_names), "positive definite")
    expect_error(gmrf(laplacian, ridge = 0), "positive")
    expect_error(gmrf(laplacian, ridge = Inf), "positive")
})

test_that("regularized GMRF operators must produce a valid precision", {
    operator <- named_sparse(matrix(c(-2, 0, 0, 1), 2), c("a", "b"))

    expect_error(
        compile_structure_spec(gmrf(operator, ridge = 1), c("a", "b")),
        "positive definite"
    )
})

test_that("generic precision validation checks modeled mappings", {
    Q <- named_sparse(diag(2), c("a", "b"))
    construction <- list(
        Q = Q,
        type = "precision",
        entity_names = c("a", "b"),
        modeled_index = c(1, 1),
        modeled_scale = c(1, 1)
    )

    expect_error(validate_precision(construction, "Cell"), "modeled indices")
    construction$modeled_index <- c(1, 2)
    construction$modeled_scale <- c(1, 0)
    expect_error(validate_precision(construction, "Cell"), "modeled scales")
})

test_that("tree and nested structures use the common finalizer", {
    hierarchy <- data.frame(
        node = c("root", "a", "b"),
        parent = c(NA, "root", "root"),
        edge = c(1, 2, 3)
    )
    tree_spec <- tree(node, parent, hierarchy, edge)
    compiled_tree <- compile_structure_spec(tree_spec, c("a", "b"))

    expect_identical(compiled_tree$modeled_index, c(2L, 3L))
    expect_equal(compiled_tree$modeled_scale, sqrt(c(3, 4)))

    treatments <- data.frame(
        key = c("a:1", "a:2", "b:1"),
        drug = c("a", "a", "b"),
        dose = c("1", "2", "1")
    )
    laplacian <- named_sparse(matrix(c(1, -1, -1, 1), 2), c("a", "b"))
    nested_gmrf <- compile_treatment_structure(gmrf(laplacian), nested(2), treatments)
    nested_tree <- compile_treatment_structure(tree_spec, nested(2), treatments)

    expect_identical(nested_gmrf$modeled_index, 3:5)
    expect_identical(nested_tree$modeled_index, 4:6)
    expect_identical(nested_gmrf$entity_names, treatments$key)
    expect_identical(nested_tree$entity_names, treatments$key)
    expect_true(inherits(nested_gmrf, "compiled_combo_structure"))
    expect_true(inherits(nested_tree, "compiled_combo_structure"))
})
