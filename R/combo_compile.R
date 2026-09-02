# Compilation of model specifications and invariant data frames.

combo_id_text <- function(x) {
    missing <- is.na(x)
    if (inherits(x, "Date") || inherits(x, "POSIXt")) {
        result <- format(x, digits = 17, trim = TRUE)
        result[missing] <- NA_character_
        return(result)
    }
    if (is.factor(x)) {
        x <- as.character(x)
    }
    if (is.numeric(x)) {
        # `format()` chooses a common display precision from the full vector,
        # so the same number can receive different text when it appears in
        # different treatment-position columns. Seventeen significant digits
        # round-trip a double while remaining independent of vector context.
        x[x == 0 & !is.na(x)] <- 0
        result <- sprintf("%.17g", as.numeric(x))
        result[missing] <- NA_character_
        return(result)
    }
    result <- as.character(x)
    result[missing] <- NA_character_
    result
}

treatment_key <- function(drug, dose) {
    if (length(drug) != length(dose)) {
        cli::cli_abort("drug and dose must have the same length")
    }
    drug <- combo_id_text(drug)
    dose <- combo_id_text(dose)
    missing_pair <- is.na(drug) & is.na(dose)
    partial <- xor(is.na(drug), is.na(dose))
    if (any(partial)) {
        cli::cli_abort("Drug and dose must be jointly present or jointly missing")
    }
    result <- paste0(nchar(drug), ":", drug, "|", nchar(dose), ":", dose)
    result[missing_pair] <- NA_character_
    result
}

validate_experiments <- function(data, require_response = TRUE) {
    if (!inherits(data, "data.frame")) {
        cli::cli_abort("data must be a data frame")
    }
    # Required columns
    required <- c("cell", "drug_1", "dose_1", "drug_2", "dose_2")
    if (require_response) {
        required <- c(required, "response")
    }
    missing <- setdiff(required, names(data))
    if (length(missing)) {
        cli::cli_abort("data is missing required columns: {.and {missing}}")
    }
    # Response
    if (!"response" %in% names(data)) {
        data$response <- NA_real_
    }
    if (!is.numeric(data$response) && all(is.na(data$response))) {
        response <- rep(NA_real_, length(data$response))
    } else if (!is.numeric(data$response)) {
        cli::cli_abort("response must be numeric")
    } else {
        response <- as.numeric(data$response)
    }
    if (any(!is.finite(response[!is.na(response)]))) {
        cli::cli_abort("Nonmissing responses must be finite numeric values")
    }

    cell <- combo_id_text(data$cell)
    if (anyNA(cell) || any(!nzchar(cell))) {
        cli::cli_abort("cell values must be nonmissing and nonempty")
    }

    list(
        data = data,
        cell = cell,
        drug_1 = combo_id_text(data$drug_1),
        dose_1 = combo_id_text(data$dose_1),
        drug_2 = combo_id_text(data$drug_2),
        dose_2 = combo_id_text(data$dose_2),
        key_1 = treatment_key(data$drug_1, data$dose_1),
        key_2 = treatment_key(data$drug_2, data$dose_2),
        response = response
    )
}

combo_metadata_by_key <- function(data, key, values, label) {
    base <- data.frame(key_value = values, stringsAsFactors = FALSE, check.names = FALSE)
    names(base) <- key
    if (is.null(data)) {
        return(base)
    }
    if (!inherits(data, "data.frame")) {
        cli::cli_abort("{label} must be a data frame")
    }
    if (!key %in% names(data)) {
        cli::cli_abort("{label} must contain key column {key}")
    }
    metadata_key <- combo_id_text(data[[key]])
    if (is_invalid_key(metadata_key)) {
        cli::cli_abort("{label} key {key} must be nonmissing and unique")
    }
    data[[key]] <- metadata_key
    result <- merge(base, data, by = key, all.x = TRUE, sort = FALSE)
    # merge sends umatched rows to the end
    result <- result[match(values, result[[key]]), , drop = FALSE]
    rownames(result) <- NULL
    result
}

combo_treatment_metadata <- function(treatments, compound_data) {
    base <- data.frame(
        drug = treatments$drug,
        dose = treatments$dose,
        treatment = treatments$key,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )
    if (is.null(compound_data)) {
        return(base)
    }
    if (!inherits(compound_data, "data.frame") || !"drug" %in% names(compound_data)) {
        cli::cli_abort("compound_data must be a data frame containing drug")
    }
    compound_key <- combo_id_text(compound_data$drug)
    if (is_invalid_key(compound_key)) {
        cli::cli_abort("compound_data drug keys must be nonmissing and unique")
    }
    reserved <- intersect(setdiff(names(compound_data), "drug"), names(base))
    if (length(reserved)) {
        cli::cli_abort("compound_data uses reserved columns: {.and {reserved}}")
    }
    selected <- match(treatments$drug, compound_key)
    for (name in setdiff(names(compound_data), "drug")) {
        base[[name]] <- compound_data[[name]][selected]
    }
    base
}

combo_compile_mean <- function(component, metadata, label) {
    if (formula_is_zero(component$mean)) {
        values <- matrix(
            numeric(),
            nrow = nrow(metadata),
            ncol = 0L,
            dimnames = list(NULL, character())
        )
        return(list(
            X = values,
            center = numeric(),
            scale = numeric(),
            formula = component$mean,
            shrinkage = NULL
        ))
    }

    # Construct the design matrix and validate it.
    frame <- tryCatch(
        stats::model.frame(
            component$mean,
            data = metadata,
            na.action = stats::na.pass
        ),
        error = function(error) {
            cli::cli_abort(
                "{label} mean formula could not be evaluated: {conditionMessage(error)}"
            )
        }
    )
    X <- stats::model.matrix(component$mean, data = frame)
    if ("(Intercept)" %in% colnames(X)) {
        cli::cli_abort("{label} mean formula may not contain an intercept")
    }
    storage.mode(X) <- "double"
    if (any(!is.finite(X))) {
        cli::cli_abort("{label} mean design matrix contains missing or nonfinite values")
    }
    # Center and scale the design matrix
    center <- colMeans(X)
    centered <- sweep(X, 2L, center, `-`)
    scale <- apply(centered, 2L, stats::sd)
    if (any(!is.finite(scale)) || any(scale <= 0)) {
        bad <- colnames(X)[!is.finite(scale) | scale <= 0] # Columns with zero variance
        cli::cli_abort(
            "{label} mean design contains constant columns: {.and {bad}}"
        )
    }
    X <- sweep(centered, 2L, scale, `/`)
    names(center) <- colnames(X)
    names(scale) <- colnames(X)
    list(
        X = X,
        center = center,
        scale = scale,
        formula = component$mean,
        shrinkage = component$mean_shrinkage
    )
}


align_structure_matrix <- function(Q, entity_names, label, source) {
    if (!methods::is(Q, "sparseMatrix")) {
        cli::cli_abort("{label} {source} must be a Matrix sparse matrix")
    }
    row_names <- rownames(Q)
    column_names <- colnames(Q)
    valid_names <- !is.null(row_names) && !is.null(column_names) &&
        !is_invalid_key(row_names) && identical(row_names, column_names)
    if (!valid_names) {
        cli::cli_abort(
            "{label} {source} must have identical unique row and column names"
        )
    }
    if (length(row_names) != length(entity_names) ||
            !setequal(row_names, entity_names)) {
        cli::cli_abort("{label} {source} names must exactly match modeled entities")
    }
    Q[entity_names, entity_names, drop = FALSE]
}

construct_iid_precision <- function(entity_names) {
    Q <- Matrix::Diagonal(length(entity_names), x = 1)
    dimnames(Q) <- list(entity_names, entity_names)
    list(
        Q = Q,
        type = "iid",
        entity_names = entity_names,
        modeled_index = seq_along(entity_names),
        modeled_scale = rep(1, length(entity_names))
    )
}

construct_supplied_precision <- function(spec, entity_names, label) {
    list(
        Q = align_structure_matrix(spec$Q, entity_names, label, "precision"),
        type = "precision",
        entity_names = entity_names,
        modeled_index = seq_along(entity_names),
        modeled_scale = rep(1, length(entity_names))
    )
}

construct_structure <- function(spec, entity_names, label) {
    type <- structure_type(spec)
    if (type == "iid") {
        return(construct_iid_precision(entity_names))
    }
    if (type == "precision") {
        return(construct_supplied_precision(spec, entity_names, label))
    }
    if (type == "gmrf") {
        return(construct_gmrf_precision(spec, entity_names, label))
    }
    if (type == "tree") {
        return(construct_tree_precision(spec, entity_names, label))
    }
    cli::cli_abort("Unsupported structure type: {type}")
}

validate_precision <- function(construction, label) {
    Q <- construction$Q
    if (!methods::is(Q, "sparseMatrix")) {
        cli::cli_abort("{label} precision must be a Matrix sparse matrix")
    }
    row_names <- rownames(Q)
    column_names <- colnames(Q)
    valid_names <- !is.null(row_names) && !is.null(column_names) &&
        !is_invalid_key(row_names) && identical(row_names, column_names)
    if (nrow(Q) < 1L || !valid_names || !Matrix::isSymmetric(Q, checkDN = TRUE, tol = 0)) {
        cli::cli_abort("{label} precision must be a non-empty exactly symmetric matrix with valid names")
    }
    if (!all(is.finite(sparse_values(Q)))) {
        cli::cli_abort("{label} precision entries must be finite")
    }

    entity_names <- construction$entity_names
    modeled_index <- construction$modeled_index
    modeled_scale <- construction$modeled_scale
    if (is_invalid_key(entity_names)) {
        cli::cli_abort("{label} modeled entity names must be unique and non-empty")
    }
    valid_index <- length(modeled_index) == length(entity_names) &&
        !anyNA(modeled_index) &&
        all(modeled_index == as.integer(modeled_index)) &&
        all(modeled_index >= 1L & modeled_index <= nrow(Q)) &&
        !anyDuplicated(modeled_index)
    if (!valid_index) {
        cli::cli_abort("{label} modeled indices must be unique valid matrix rows")
    }
    valid_scale <- length(modeled_scale) == length(entity_names) &&
        is.numeric(modeled_scale) && all(is.finite(modeled_scale)) &&
        all(modeled_scale > 0)
    if (!valid_scale) {
        cli::cli_abort("{label} modeled scales must be finite positive numbers")
    }

    Q <- Matrix::drop0(Q)
    Q <- methods::as(Q, "generalMatrix")
    Q <- methods::as(Q, "CsparseMatrix") * 1
    tryCatch(
        Matrix::Cholesky(
            Matrix::forceSymmetric(Q),
            LDL = FALSE,
            perm = TRUE
        ),
        warning = function(warning) {
            cli::cli_abort(
                "{label} precision must be positive definite: {conditionMessage(warning)}"
            )
        },
        error = function(error) {
            cli::cli_abort(
                "{label} precision must be positive definite: {conditionMessage(error)}"
            )
        }
    )
    construction$Q <- Q
    construction$modeled_index <- as.integer(modeled_index)
    construction$modeled_scale <- as.numeric(modeled_scale)
    construction
}

finalize_structure <- function(construction, label) {
    construction <- validate_precision(construction, label)
    Q <- construction$Q
    entries <- Matrix::summary(Q)
    neighbors <- vector("list", nrow(Q))
    coefficients <- vector("list", nrow(Q))
    off_diagonal <- entries$i != entries$j
    for (node in seq_len(nrow(Q))) {
        selected <- off_diagonal & entries$i == node
        neighbors[[node]] <- as.integer(entries$j[selected])
        coefficients[[node]] <- as.numeric(entries$x[selected])
    }
    structure(
        list(
            type = construction$type,
            Q = Q,
            nodes = rownames(Q),
            diagonal = as.numeric(Matrix::diag(Q)),
            neighbors = neighbors,
            coefficients = coefficients,
            modeled_index = construction$modeled_index,
            modeled_scale = construction$modeled_scale,
            latent_index = setdiff(seq_len(nrow(Q)), construction$modeled_index),
            entity_names = construction$entity_names,
            iid = identical(construction$type, "iid")
        ),
        class = "compiled_combo_structure"
    )
}

compile_flat_structure <- function(spec, entity_names, label) {
    finalize_structure(construct_structure(spec, entity_names, label), label)
}

append_nested_precision <- function(base, treatments, relative_precision) {
    Q <- methods::as(base$Q, "generalMatrix")
    Q <- methods::as(Q, "CsparseMatrix") * 1
    n_base <- nrow(Q)
    n_treatments <- nrow(treatments)
    entity_position <- match(treatments$drug, base$entity_names)
    parent <- base$modeled_index[entity_position]
    leaf <- n_base + seq_len(n_treatments)
    node_names <- c(rownames(Q), paste0(".dose_", seq_len(n_treatments)))
    entries <- Matrix::summary(Q)
    edge_i <- c(parent, leaf, parent, leaf)
    edge_j <- c(parent, leaf, leaf, parent)
    edge_x <- c(
        rep(relative_precision, 2L * n_treatments),
        rep(-relative_precision, 2L * n_treatments)
    )
    augmented <- Matrix::sparseMatrix(
        i = c(entries$i, edge_i),
        j = c(entries$j, edge_j),
        x = c(entries$x, edge_x),
        dims = rep(length(node_names), 2L),
        dimnames = list(node_names, node_names),
        repr = "C"
    )
    list(
        Q = Matrix::drop0(augmented),
        type = paste0("nested_", base$type),
        entity_names = treatments$key,
        modeled_index = leaf,
        modeled_scale = rep(1, n_treatments)
    )
}

compile_treatment_structure <- function(spec, dose, treatments) {
    if (dose$type == "categorical") {
        return(compile_flat_structure(
            spec,
            treatments$key,
            "Treatment"
        ))
    }
    compounds <- unique(treatments$drug)
    type <- structure_type(spec)
    if (type == "tree") {
        construction <- append_nested_tree(
            spec,
            compounds,
            treatments,
            dose$precision
        )
        return(finalize_structure(construction, "Treatment"))
    }
    base <- construct_structure(spec, compounds, "Compound")
    augmented <- append_nested_precision(
        base,
        treatments,
        dose$precision
    )
    finalize_structure(augmented, "Treatment")
}

combo_compile_component <- function(
    component,
    name,
    kind,
    side,
    n_dimensions,
    metadata,
    entity_names,
    dose,
    treatments = NULL
) {
    if (is.null(component)) {
        return(NULL)
    }
    mean <- combo_compile_mean(component, metadata, name)
    structure <- if (side == "cell") {
        compile_flat_structure(
            component$structure,
            entity_names,
            "Cell"
        )
    } else {
        compile_treatment_structure(
            component$structure,
            dose,
            treatments
        )
    }
    shrinkage <- component$shrinkage
    if (shrinkage$type == "multiplicative_gamma") {
        for (parameter in c("shape", "rate")) {
            value <- shrinkage[[parameter]]
            if (length(value) == 1L) {
                shrinkage[[parameter]] <- rep(value, n_dimensions)
            } else if (length(value) != n_dimensions) {
                cli::cli_abort(
                    "{name} multiplicative_gamma() {parameter} must have length 1 or match rank ({n_dimensions})"
                )
            }
        }
    }
    list(
        name = name,
        kind = kind,
        side = side,
        n_entities = length(entity_names),
        n_dimensions = n_dimensions,
        mean = mean,
        structure = structure,
        shrinkage = shrinkage
    )
}

compile_combo_design <- function(
    model,
    data,
    cell_data = NULL,
    compound_data = NULL,
    require_response = FALSE
) {
    if (!inherits(model, "combo_model")) {
        cli::cli_abort("model must be constructed by combo_model()")
    }
    parsed <- validate_experiments(data, require_response = require_response)
    cells <- unique(parsed$cell)
    all_treatment_keys <- c(parsed$key_1, parsed$key_2)
    treatment_keys <- unique(all_treatment_keys)
    treatment_keys <- treatment_keys[!is.na(treatment_keys)]
    first_drug <- c(parsed$drug_1, parsed$drug_2)
    first_dose <- c(parsed$dose_1, parsed$dose_2)
    first_key <- c(parsed$key_1, parsed$key_2)
    treatment_position <- match(treatment_keys, first_key)
    treatments <- data.frame(
        key = treatment_keys,
        drug = first_drug[treatment_position],
        dose = first_dose[treatment_position],
        stringsAsFactors = FALSE
    )
    cell_metadata <- combo_metadata_by_key(
        cell_data,
        "cell",
        cells,
        "cell_data"
    )
    treatment_metadata <- combo_treatment_metadata(
        treatments,
        compound_data
    )
    cell_index <- match(parsed$cell, cells)
    treatment_index_1 <- match(parsed$key_1, treatment_keys)
    treatment_index_2 <- match(parsed$key_2, treatment_keys)
    treatment_index_1[is.na(parsed$key_1)] <- 0L
    treatment_index_2[is.na(parsed$key_2)] <- 0L

    specs <- list(
        cell_offset = list(kind = "offset", side = "cell", dims = 1L),
        cell_factors = list(kind = "factor", side = "cell", dims = model$rank),
        treatment_offset = list(kind = "offset", side = "treatment", dims = 1L),
        treatment_main_factors = list(
            kind = "factor", side = "treatment", dims = model$rank
        ),
        treatment_interaction_factors = list(
            kind = "factor", side = "treatment", dims = model$rank
        )
    )
    components <- vector("list", length(specs))
    names(components) <- names(specs)
    for (name in names(specs)) {
        definition <- specs[[name]]
        metadata <- if (definition$side == "cell") {
            cell_metadata
        } else {
            treatment_metadata
        }
        entities <- if (definition$side == "cell") {
            cells
        } else {
            treatment_keys
        }
        components[[name]] <- combo_compile_component(
            model$components[[name]],
            name,
            definition$kind,
            definition$side,
            definition$dims,
            metadata,
            entities,
            model$dose,
            treatments
        )
    }
    active_components <- Filter(Negate(is.null), components)
    fast_iid <- model$dose$type == "categorical" &&
        all(vapply(
            active_components,
            function(component) component$structure$iid,
            logical(1)
        ))
    structure(
        list(
            model = model,
            data = parsed$data,
            validated_response = parsed$response,
            all_cell = cell_index,
            all_treatment_1 = treatment_index_1,
            all_treatment_2 = treatment_index_2,
            cells = cells,
            treatments = treatments,
            cell_metadata = cell_metadata,
            treatment_metadata = treatment_metadata,
            components = components,
            observation_prior = model$family$precision,
            fast_iid = fast_iid
        ),
        class = "compiled_combo_design"
    )
}

compile_combo_model <- function(
    model,
    data,
    cell_data = NULL,
    compound_data = NULL
) {
    compiled <- compile_combo_design(
        model,
        data,
        cell_data = cell_data,
        compound_data = compound_data,
        require_response = TRUE
    )
    response <- compiled$validated_response
    observed <- which(!is.na(response))
    if (!length(observed)) {
        cli::cli_abort("fit_combo requires at least one observed response")
    }

    active_v2 <- !is.null(
        model$components$treatment_interaction_factors
    )
    self <- compiled$all_treatment_1[observed] > 0L &
        compiled$all_treatment_1[observed] ==
            compiled$all_treatment_2[observed]
    if (active_v2 && any(self)) {
        cli::cli_abort("Observed self-combinations are not conditionally Gaussian for treatment_interaction_factors")
    }

    cell_exposure <- tabulate(
        compiled$all_cell[observed],
        nbins = length(compiled$cells)
    )
    treatment_exposure <- tabulate(
        c(
            compiled$all_treatment_1[observed],
            compiled$all_treatment_2[observed]
        ),
        nbins = nrow(compiled$treatments)
    )
    for (name in names(compiled$components)) {
        component <- compiled$components[[name]]
        if (is.null(component)) {
            next
        }
        component$exposure_count <- if (component$side == "cell") {
            cell_exposure
        } else {
            treatment_exposure
        }
        compiled$components[[name]] <- component
    }

    alpha <- if (model$global_intercept$type == "empirical") {
        mean(response[observed])
    } else {
        model$global_intercept$value
    }
    structure(
        list(
            model = compiled$model,
            data = compiled$data,
            response = response[observed],
            observed_rows = observed,
            cell = compiled$all_cell[observed],
            treatment_1 = compiled$all_treatment_1[observed],
            treatment_2 = compiled$all_treatment_2[observed],
            all_cell = compiled$all_cell,
            all_treatment_1 = compiled$all_treatment_1,
            all_treatment_2 = compiled$all_treatment_2,
            cells = compiled$cells,
            treatments = compiled$treatments,
            cell_metadata = compiled$cell_metadata,
            treatment_metadata = compiled$treatment_metadata,
            components = compiled$components,
            alpha = alpha,
            observation_prior = compiled$observation_prior,
            fast_iid = compiled$fast_iid
        ),
        class = "compiled_combo_model"
    )
}
