#' Correlate dose deviations within a compound
#'
#' Builds a correlation specification for [nested()] from the concentrations
#' already present in the experiment design; no additional metadata are needed.
#' For each treatment component and factor dimension, raw dose-to-parent
#' deviations for compound `a` have covariance `K_a / (rho * precision)`, where
#' `rho` is the component's shrinkage precision and `precision` is supplied to
#' [nested()]. Different compounds have independent deviations conditional on
#' their parents. Tree-backed structures apply their existing path-variance
#' normalization after constructing these raw deviations.
#'
#' With `d = abs(t_i - t_j) / length_scale` for transformed doses `t`, the
#' exponential correlation is `exp(-d)` and the Matern 3/2 correlation is
#' `(1 + sqrt(3) * d) * exp(-sqrt(3) * d)`. Both have unit diagonal.
#' The exponential kernel provides local pooling with a sparse precision;
#' Matern 3/2 assumes smoother dose dependence. Neither imposes monotonicity.
#'
#' Doses must be finite numeric values or numeric text labels. Log10 requires
#' positive concentrations; identity also accepts already-transformed coordinates.
#' Within a compound, distinct treatment labels must not collapse to the same
#' numeric or transformed dose. Missing treatments retain their usual control
#' representation. A compound with one dose has correlation matrix `[1]`.
#'
#' The fixed length scale is shared across compounds and enabled treatment
#' components. A log10 length scale of one is one concentration decade; no
#' per-compound rescaling is performed. Component shrinkage and the relative
#' precision in [nested()] control amplitude. Entity-local shrinkage remains
#' incompatible with nested structures. Kernels do not enable prediction at
#' treatments absent from the fitted mappings; include unobserved treatments in
#' the fitting design with missing responses. No jitter is added to singular
#' or numerically unstable kernel constructions.
#'
#' @param kernel `"exponential"` (default) or `"matern32"`.
#' @param transform `"log10"` (default) or `"identity"`.
#' @param length_scale A finite positive scalar in transformed dose units.
#' @return A `combo_dose_kernel` specification for the `correlation` argument
#'   of [nested()].
#' @examples
#' component <- combo_gaussian_component(shrinkage = gamma_precision())
#' model <- combo_model(
#'     rank = 2,
#'     treatment_offset = component,
#'     treatment_main_factors = component,
#'     treatment_interaction_factors = component,
#'     dose = nested(precision = 4, correlation = dose_kernel("matern32"))
#' )
#' design <- data.frame(
#'     cell = "A", drug_1 = "X", dose_1 = c(0.1, 1, 10),
#'     drug_2 = NA_character_, dose_2 = NA_real_
#' )
#' prior_predict(model, design, draws = 3, seed = 42, intercept = 0)
#' @export
dose_kernel <- function(
    kernel = c("exponential", "matern32"),
    transform = c("log10", "identity"),
    length_scale = 1
) {
    kernel <- tryCatch(match.arg(kernel), error = function(error) {
        cli::cli_abort("kernel must be exponential or matern32")
    })
    transform <- tryCatch(match.arg(transform), error = function(error) {
        cli::cli_abort("transform must be log10 or identity")
    })
    structure(
        list(
            kernel = kernel,
            transform = transform,
            length_scale = param_scalar_positive(length_scale, "length_scale")
        ),
        class = "combo_dose_kernel"
    )
}

dose_kernel_label <- function(correlation, compound) {
    paste0(
        "Dose correlation for compound ", paste(compound, collapse = ", "),
        " (", correlation$kernel, ", ", correlation$transform,
        ", length scale ", format(correlation$length_scale), ")"
    )
}

# Coordinates are sorted here only; callers map the precision back to the
# original treatment order. Use each distinct treatment once, not each row.
dose_kernel_precision <- function(doses, correlation, compound) {
    label <- dose_kernel_label(correlation, compound)
    coordinates <- suppressWarnings(if (is.numeric(doses)) {
        as.numeric(doses)
    } else {
        as.numeric(as.character(doses))
    })
    if (any(!is.finite(coordinates))) {
        cli::cli_abort("{label}: doses must be finite numeric values")
    }
    if (correlation$transform == "log10") {
        if (any(coordinates <= 0)) {
            cli::cli_abort("{label}: log10 requires strictly positive concentrations")
        }
        coordinates <- log10(coordinates)
    }
    if (anyDuplicated(coordinates)) {
        cli::cli_abort("{label}: distinct dose labels map to duplicate coordinates")
    }
    position <- order(coordinates)
    coordinates <- coordinates[position]
    n <- length(coordinates)
    Q <- tryCatch({
        if (n == 1L) {
            Matrix::Diagonal(1L)
        } else if (correlation$kernel == "exponential") {
            # Stationary AR(1) transitions on an irregular grid:
            # z[1]^2 + sum((z[i+1] - r[i]*z[i])^2 / (1-r[i]^2)).
            distance <- diff(coordinates) / correlation$length_scale
            r <- exp(-distance)
            innovation <- -expm1(-2 * distance)
            weight <- 1 / innovation
            diagonal <- c(1, rep(0, n - 1L)) + c(r^2 * weight, 0) + c(0, weight)
            off <- -r * weight
            Matrix::sparseMatrix(
                i = c(seq_len(n), seq_len(n - 1L), 2:n),
                j = c(seq_len(n), 2:n, seq_len(n - 1L)),
                x = c(diagonal, off, off), dims = c(n, n)
            )
        } else {
            distance <- abs(outer(coordinates, coordinates, "-"))
            scaled <- sqrt(3) * distance / correlation$length_scale
            K <- exp(log1p(scaled) - scaled)
            K[is.infinite(scaled)] <- 0
            Matrix::Matrix(chol2inv(chol(K)), sparse = TRUE)
        }
    }, error = function(error) {
        cli::cli_abort("{label}: kernel construction failed: {conditionMessage(error)}")
    })
    if (any(!is.finite(sparse_values(Q)))) {
        cli::cli_abort("{label}: kernel precision contains non-finite values")
    }
    tryCatch(
        Matrix::Cholesky(Matrix::forceSymmetric(Q), LDL = FALSE),
        warning = function(error) {
            cli::cli_abort("{label}: kernel precision is not numerically positive definite: {conditionMessage(error)}")
        },
        error = function(error) {
            cli::cli_abort("{label}: kernel precision is not numerically positive definite: {conditionMessage(error)}")
        }
    )
    Q[order(position), order(position), drop = FALSE]
}

append_correlated_doses <- function(spec, dose, treatments) {
    compounds <- unique(treatments$drug)
    base <- validate_precision(construct_structure(spec, compounds, "Compound"), "Compound")
    n_base <- nrow(base$Q)
    n_treatments <- nrow(treatments)
    parent <- base$modeled_index[match(treatments$drug, compounds)]
    leaf <- n_base + seq_len(n_treatments)
    node_names <- c(rownames(base$Q), paste0(".dose_", seq_len(n_treatments)))
    groups <- lapply(compounds, function(compound) which(treatments$drug == compound))
    blocks <- lapply(
        seq_along(compounds),
        function(i) dose_kernel_precision(treatments$dose[groups[[i]]], dose$correlation, compounds[i])
    )
    R <- Matrix::bdiag(blocks)
    original_order <- order(unlist(groups, use.names = FALSE))
    R <- dose$precision * R[original_order, original_order, drop = FALSE]
    D <- Matrix::sparseMatrix(
        i = rep(seq_len(n_treatments), 2),
        j = c(parent, leaf),
        x = rep(c(-1, 1), each = n_treatments),
        dims = c(n_treatments, length(node_names))
    )
    Q <- Matrix::bdiag(base$Q, Matrix::Diagonal(n_treatments, x = 0))
    Q <- Q + Matrix::crossprod(D, R %*% D)
    # Store both triangles identically for the common precision finalizer.
    Q <- Matrix::forceSymmetric(Q, uplo = "U")
    dimnames(Q) <- list(node_names, node_names)
    scale <- if (base$type == "tree") {
        sqrt(base$modeled_scale[match(treatments$drug, compounds)]^2 + 1 / dose$precision)
    } else {
        rep(1, n_treatments)
    }
    finalize_structure(
        list(
            Q = Q,
            type = paste0("nested_correlated_", base$type),
            entity_names = treatments$key,
            modeled_index = leaf,
            modeled_scale = scale
        ),
        dose_kernel_label(dose$correlation, compounds)
    )
}
