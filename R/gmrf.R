#' @rdname context_specifications
#' @export
gmrf <- function(operator, ridge = 1) {
    structure(
        list(
            type = "gmrf",
            operator = operator,
            ridge = param_scalar_positive(ridge, "ridge")
        ),
        class = "structural_prior"
    )
}

construct_gmrf_precision <- function(spec, entity_names, label) {
    operator <- align_structure_matrix(
        spec$operator,
        entity_names,
        label,
        "GMRF operator"
    )
    ridge <- Matrix::Diagonal(nrow(operator), x = spec$ridge)
    dimnames(ridge) <- dimnames(operator)
    list(
        Q = operator + ridge,
        type = "gmrf",
        entity_names = entity_names,
        modeled_index = seq_along(entity_names),
        modeled_scale = rep(1, length(entity_names))
    )
}

# TODO: implement graph_laplacian() to construct a GMRF operator from a graph adjacency matrix.
#       then either gmrf_from_graph() of gmrf.graph() for igraph objects can be implemented to construct a GMRF prior from a graph.