#' @rdname context_specifications
#' @export
tree <- function(node, parent, data, edge_length = 1) {
    data <- eval(substitute(data), envir = parent.frame())
    if (!inherits(data, "data.frame")) {
        stop("tree data must be a data frame", call. = FALSE)
    }
    nodes <- tree_column(substitute(node), data, "node")
    parents <- tree_column(substitute(parent), data, "parent")
    edge_lengths <- if (missing(edge_length)) {
        1
    } else {
        expression <- substitute(edge_length)
        if (is.symbol(expression) && as.character(expression) %in% names(data)) {
            data[[as.character(expression)]]
        } else {
            eval(expression, envir = parent.frame())
        }
    }
    structure(
        list(
            type = "tree",
            node = nodes,
            parent = parents,
            edge_length = edge_lengths
        ),
        class = "structural_prior"
    )
}

tree_column <- function(expression, data, label) {
    if (!is.symbol(expression)) {
        stop(label, " must name a column in data", call. = FALSE)
    }
    column_name <- as.character(expression)
    if (!column_name %in% names(data)) {
        stop(label, " column is absent from data: ", column_name, call. = FALSE)
    }
    data[[column_name]]
}

tree_precision_matrix <- function(tree) {
    n <- length(tree$nodes)
    diagonal <- numeric(n)
    ii <- integer()
    jj <- integer()
    xx <- numeric()
    for (node in seq_len(n)) {
        weight <- 1 / tree$edge_length[node]
        diagonal[node] <- diagonal[node] + weight
        parent <- tree$parent_index[node]
        if (parent != 0L) {
            diagonal[parent] <- diagonal[parent] + weight
            ii <- c(ii, node, parent)
            jj <- c(jj, parent, node)
            xx <- c(xx, -weight, -weight)
        }
    }
    Matrix::sparseMatrix(
        i = c(seq_len(n), ii),
        j = c(seq_len(n), jj),
        x = c(diagonal, xx),
        dims = c(n, n),
        dimnames = list(tree$nodes, tree$nodes),
        repr = "C"
    )
}

derive_tree <- function(tree_spec, modeled_index) {
    n_nodes <- length(tree_spec$nodes)
    parent_index <- tree_spec$parent_index
    children <- vector("list", n_nodes)
    for (node_index in seq_len(n_nodes)) {
        parent <- parent_index[node_index]
        if (parent != 0L) {
            children[[parent]] <- c(children[[parent]], node_index)
        }
    }

    depth <- integer(n_nodes)
    path_length <- numeric(n_nodes)
    resolved <- logical(n_nodes)
    resolve_node <- function(node_index) {
        if (resolved[node_index]) {
            return(invisible(NULL))
        }
        parent <- parent_index[node_index]
        if (parent != 0L) {
            resolve_node(parent)
            depth[node_index] <<- depth[parent] + 1L
            path_length[node_index] <<- path_length[parent] + tree_spec$edge_length[node_index]
        } else {
            depth[node_index] <<- 1L
            path_length[node_index] <<- tree_spec$edge_length[node_index]
        }
        resolved[node_index] <<- TRUE
        invisible(NULL)
    }
    for (node_index in seq_len(n_nodes)) {
        resolve_node(node_index)
    }

    if (anyDuplicated(modeled_index)) {
        stop("Each modeled parameter must map to a distinct tree leaf", call. = FALSE)
    }
    structure(
        c(
            tree_spec,
            list(
                children = children,
                depth = depth,
                path_length = path_length,
                modeled_index = as.integer(modeled_index),
                modeled_scale = sqrt(path_length[modeled_index]),
                latent_index = setdiff(seq_len(n_nodes), modeled_index),
                reverse_order = order(depth, decreasing = TRUE)
            )
        ),
        class = "batchie_gaussian_tree"
    )
}

validate_tree <- function(prior, entity_names, label) {
    nodes <- as.character(prior$node)
    parents <- as.character(prior$parent)
    parents[is.na(prior$parent) | !nzchar(parents)] <- NA_character_
    edge_lengths <- as.numeric(prior$edge_length)
    if (length(edge_lengths) == 1L) {
        edge_lengths <- rep(edge_lengths, length(nodes))
    }

    if (is_invalid_key(nodes)) {
        stop(label, " tree nodes must be unique non-empty strings", call. = FALSE)
    }
    if (length(parents) != length(nodes)) {
        stop(label, " tree parent column must have one value per node", call. = FALSE)
    }
    if (length(edge_lengths) != length(nodes) ||
            any(!is.finite(edge_lengths)) || any(edge_lengths <= 0)) {
        stop(label, " tree edge lengths must be finite and positive", call. = FALSE)
    }

    parent_index <- match(parents, nodes)
    invalid_parent <- !is.na(parents) & is.na(parent_index)
    if (any(invalid_parent)) {
        stop(label, " tree contains parents absent from its node column: ",
            paste(unique(parents[invalid_parent]), collapse = ", "),
            call. = FALSE
        )
    }
    parent_index[is.na(parent_index)] <- 0L
    if (any(parent_index == seq_along(nodes))) {
        stop(label, " tree nodes may not be their own parent", call. = FALSE)
    }

    seen <- logical(length(nodes))
    for (start in seq_along(nodes)) {
        seen[] <- FALSE
        current <- start
        while (current != 0L) {
            if (seen[current]) {
                stop(label, " tree contains a cycle", call. = FALSE)
            }
            seen[current] <- TRUE
            current <- parent_index[current]
        }
    }

    children <- split(
        seq_along(nodes)[parent_index != 0L],
        parent_index[parent_index != 0L]
    )
    has_children <- logical(length(nodes))
    has_children[as.integer(names(children))] <- TRUE

    missing_entities <- setdiff(entity_names, nodes)
    if (length(missing_entities)) {
        stop(label, " tree is missing modeled entities: ",
            paste(missing_entities, collapse = ", "),
            call. = FALSE
        )
    }
    entity_index <- match(entity_names, nodes)
    nonterminal_entities <- entity_names[has_children[entity_index]]
    if (length(nonterminal_entities)) {
        stop(label, " modeled entities must be terminal tree nodes: ",
            paste(nonterminal_entities, collapse = ", "),
            call. = FALSE
        )
    }

    list(
        nodes = nodes,
        parent_index = as.integer(parent_index),
        edge_length = edge_lengths,
        entity_index = as.integer(entity_index)
    )
}

construct_tree_precision <- function(spec, entity_names, label) {
    tree_spec <- validate_tree(spec, entity_names, label)
    tree <- derive_tree(tree_spec, tree_spec$entity_index)
    list(
        Q = tree_precision_matrix(tree),
        type = "tree",
        entity_names = entity_names,
        modeled_index = tree$modeled_index,
        modeled_scale = tree$modeled_scale
    )
}

append_nested_tree <- function(spec, compounds, treatments, precision) {
    tree_spec <- validate_tree(spec, compounds, "Compound")
    parent <- match(treatments$drug, tree_spec$nodes)
    n_treatments <- nrow(treatments)
    tree_spec$nodes <- c(
        tree_spec$nodes,
        paste0(".dose_", seq_len(n_treatments))
    )
    tree_spec$parent_index <- c(tree_spec$parent_index, parent)
    tree_spec$edge_length <- c(
        tree_spec$edge_length,
        rep(1 / precision, n_treatments)
    )
    modeled <- length(tree_spec$nodes) - n_treatments + seq_len(n_treatments)
    tree <- derive_tree(tree_spec, modeled)
    list(
        Q = tree_precision_matrix(tree),
        type = "tree",
        entity_names = treatments$key,
        modeled_index = tree$modeled_index,
        modeled_scale = tree$modeled_scale
    )
}
