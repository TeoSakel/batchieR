# Sampler-agnostic initialization of combination-model components.

init_shrinkage <- function(spec, n_entities, n_dimensions) {
    type <- spec$type
    result <- list(type = type, spec = spec)
    if (type == "fixed") {
        result[["global"]] <- rep(spec$precision, n_dimensions)
    } else if (type == "gamma") {
        result[["global"]] <- rep(100, n_dimensions)
    } else if (type == "global_half_cauchy") {
        result[["global"]] <- rep(1, n_dimensions)
        result[["global_aux"]] <- rep(1, n_dimensions)
    } else if (type == "local_half_cauchy") {
        result[["global"]] <- rep(1, n_dimensions)
        result[["local"]] <- matrix(100, nrow = n_entities, ncol = n_dimensions)
        result[["local_aux"]] <- matrix(1, nrow = n_entities, ncol = n_dimensions)
    } else if (type == "horseshoe") {
        result[["global"]] <- rep(1, n_dimensions)
        result[["global_aux"]] <- rep(1, n_dimensions)
        result[["local"]] <- matrix(100, nrow = n_entities, ncol = n_dimensions)
        result[["local_aux"]] <- matrix(1, nrow = n_entities, ncol = n_dimensions)
    } else if (type == "multiplicative_gamma") {
        # Preserve the legacy diffuse initialization; subsequent updates make
        # global exactly cumprod(delta).
        result[["global"]] <- rep(100, n_dimensions)
        result[["delta"]] <- rep(1, n_dimensions)
    } else {
        cli::cli_abort("Unsupported shrinkage state: {type}")
    }
    result
}

#' Initialize all sampler state for one compiled model component
#'
#' Converts a static compiled component into mutable Gibbs state:
#'
#' - `compiled`: static metadata, dimensions, and structural prior.
#' - `values`: entity-level effects used to construct predictions.
#' - `raw`: node-level deviations governed by the structural prior; hierarchical
#'   structures may include latent nodes that have no corresponding entity.
#' - `shrinkage`: global, local, auxiliary, or dimension-specific precision state
#'   for the component deviations, depending on its shrinkage specification.
#'
#' Entity/node semantics:
#'
#' - An entity is a modeled cell or treatment whose `values` row enters the predictor.
#' - A node is a position in the structural prior graph. Every entity maps to one
#'   modeled node, but hierarchical graphs may also contain latent nodes used
#'   only to relate entity effects and share information.
#'
#' Offsets have one state column and factors have one per rank dimension.
#' Entities map to nodes through `modeled_index`; on those nodes,
#' `raw = values * modeled_scale`. A `NULL` compiled component is
#' disabled and remains `NULL` throughout the sampler.
#'
#' @param compiled A component produced by `combo_compile_component()`, or
#'   `NULL` for a disabled component.
#' @return An internal component-state list, or `NULL`.
#' @noRd
init_component_state <- function(compiled) {
    if (is.null(compiled)) {
        return(NULL)
    }
    values <- matrix(
        0,
        nrow = compiled$n_entities,
        ncol = compiled$n_dimensions,
        dimnames = list(
            compiled$structure$entity_names,
            if (compiled$n_dimensions > 1L) {
                paste0("dim_", seq_len(compiled$n_dimensions))
            } else {
                "value"
            }
        )
    )
    raw <- matrix(
        0,
        nrow = length(compiled$structure$nodes),
        ncol = compiled$n_dimensions,
        dimnames = list(compiled$structure$nodes, colnames(values))
    )
    list(
        compiled = compiled,
        values = values,
        raw = raw,
        shrinkage = init_shrinkage(
            compiled$shrinkage,
            compiled$n_entities,
            compiled$n_dimensions
        )
    )
}
