# Public model specifications for the generic combination model.

#' Specify a combination-response model
#'
#' Defines the terms of the Bayesian latent-factor model fitted by [fit_combo()].
#' For cell `k` and treatments `i` and `j`, the linear predictor is
#'
#' ```
#' fitted_mean[k, i, j] = intercept + W0[k] + V0[i] + V0[j] +
#'     dot(W[k, ], V1[i, ] + V1[j, ] + V2[i, ] * V2[j, ])
#' ```
#'
#' where `*` denotes element-wise multiplication. The terms are:
#'
#' - `mean`: the global mean or metadata regression; see [mean_specifications].
#' - `cell_offset` (`W0`): baseline differences among cellular contexts.
#' - `treatment_offset` (`V0`): additive treatment effects shared across contexts.
#' - `cell_factors` (`W`): latent characteristics of each cellular context.
#' - `treatment_main_factors` (`V1`): latent treatment main effects. Their
#'   inner product with `W` allows treatment effects to vary by context.
#' - `treatment_interaction_factors` (`V2`): latent nonadditive combination
#'   effects. Two treatments are multiplied in latent space and projected onto
#'   `W`, allowing their interaction to vary by context.
#'
#' A single-treatment row has no second-treatment or interaction contribution.
#' Treatment order does not change the predictor. Setting a component to `NULL`
#' removes its term; treatment factor terms require `cell_factors`.
#' When `mean` uses a formula, enabled offsets remain zero-centered residual
#' entity deviations around that regression.
#'
#' See [combo_gaussian_component()] for configuring each component. The remaining model
#' specifications are described in [mean_specifications],
#' [dose_specifications], and [gaussian_response()].
#'
#' @param rank Positive integer shared latent dimension of `cell_factors`,
#'   `treatment_main_factors`, and `treatment_interaction_factors`. Defaults to 12.
#' @param mean Overall response mean or fixed-effect regression; see
#'   [mean_specifications].
#' @param cell_offset Component for cell-specific baselines (`W0`), or `NULL`.
#' @param cell_factors Component for latent cell characteristics (`W`), or `NULL`.
#' @param treatment_offset Component for shared additive treatment effects
#'   (`V0`), or `NULL`.
#' @param treatment_main_factors Component for context-specific treatment main
#'   effects (`V1`), or `NULL`.
#' @param treatment_interaction_factors Component for context-specific
#'   nonadditive effects (`V2`), or `NULL`.
#' @param dose Treatment-dose relationship; see [dose_specifications].
#' @param family Response distribution and link; see [gaussian_response()].
#' @return A `combo_model` specification consumed by [fit_combo()].
#' @examples
#' model <- combo_model(rank = 2L)
#' model
#' @export
combo_model <- function(
    rank = 12L,
    mean = empirical_mean(),
    cell_offset = combo_gaussian_component(shrinkage = gamma_precision(shape = 1.1, rate = 1.1)),
    cell_factors = combo_gaussian_component(
        shrinkage = multiplicative_gamma(
            shape = 2,
            rate = 1
        )
    ),
    treatment_offset = combo_gaussian_component(
        shrinkage = horseshoe(global_scale = 1, local_scale = 1)
    ),
    treatment_main_factors = combo_gaussian_component(
        shrinkage = horseshoe(global_scale = 1, local_scale = 1)
    ),
    treatment_interaction_factors = combo_gaussian_component(
        shrinkage = horseshoe(global_scale = 1, local_scale = 1)
    ),
    dose = categorical(),
    family = gaussian_response(
        link = "identity",
        precision = gamma_precision(shape = 1.1, rate = 1.1)
    )
) {
    rank <- param_positive_integer(rank, "rank")
    if (!inherits(mean, "combo_mean")) {
        cli::cli_abort(
            "mean must be empirical_mean(), fixed_mean(), or formula_mean()"
        )
    }
    if (!inherits(dose, "combo_dose")) {
        cli::cli_abort("dose must be categorical() or nested()")
    }
    if (!inherits(family, "combo_family")) {
        cli::cli_abort("family must be a supported family specification")
    }
    components <- list(
        cell_offset = cell_offset,
        cell_factors = cell_factors,
        treatment_offset = treatment_offset,
        treatment_main_factors = treatment_main_factors,
        treatment_interaction_factors = treatment_interaction_factors
    )
    kinds <- c(
        cell_offset = "offset",
        cell_factors = "factor",
        treatment_offset = "offset",
        treatment_main_factors = "factor",
        treatment_interaction_factors = "factor"
    )
    for (name in names(components)) {
        validate_component(components[[name]], name, kinds[[name]])
    }
    if (dose$type == "nested") {
        treatment_names <- c(
            "treatment_offset",
            "treatment_main_factors",
            "treatment_interaction_factors"
        )
        incompatible <- treatment_names[vapply(
            components[treatment_names],
            function(component) {
                !is.null(component) && component$shrinkage$type %in%
                    c("local_half_cauchy", "horseshoe")
            },
            logical(1)
        )]
        if (length(incompatible)) {
            cli::cli_abort(
                "nested() creates a non-IID treatment structure and cannot use entity-local half-Cauchy shrinkage in: {.and {incompatible}}"
            )
        }
    }
    factors_active <- !is.null(treatment_main_factors) ||
        !is.null(treatment_interaction_factors)
    if (factors_active && is.null(cell_factors)) {
        cli::cli_abort("treatment factor components require cell_factors")
    }
    if (!is.null(cell_factors) && !factors_active) {
        cli::cli_abort("cell_factors contributes nothing unless a treatment factor component is active")
    }
    structure(
        list(
            rank = rank,
            mean = mean,
            components = components,
            dose = dose,
            family = family
        ),
        class = "combo_model"
    )
}

#' Update a combination-model specification
#'
#' Creates a new specification by replacing selected arguments of [combo_model()].
#' Omitted arguments retain their stored values, including disabled components
#' and formula environments. Each supplied specification replaces its previous
#' value in full; nested fields are not merged. Set a component to `NULL` to
#' disable it. The complete result is validated by [combo_model()].
#'
#' This method updates specifications only and does not fit a model. Use
#' [fit_combo()] to fit the returned specification. Formula-update shorthand
#' and `evaluate = FALSE` are not supported. Calling `update(object)` without
#' replacements returns an identical specification. The original is unchanged.
#'
#' @param object A specification created by [combo_model()].
#' @param ... Replacements for any [combo_model()] arguments: `rank`, `mean`,
#'   `cell_offset`, `cell_factors`, `treatment_offset`, `treatment_main_factors`,
#'   `treatment_interaction_factors`, `dose`, or `family`. Each argument must
#'   have its full, unique name.
#' @return A validated `combo_model` specification.
#' @seealso [combo_model()], [fit_combo()]
#' @examples
#' model <- combo_model(rank = 2L)
#' update(model, rank = 4L)
#' update(model, treatment_interaction_factors = NULL)
#' update(model, treatment_offset = combo_gaussian_component(
#'     shrinkage = gamma_precision(shape = 2, rate = 1)
#' ))
#' @exportS3Method stats::update
update.combo_model <- function(object, ...) {
    replacements <- list(...)
    if (length(replacements)) {
        labels <- names(replacements)
        if (is.null(labels) || anyNA(labels) || any(!nzchar(labels))) {
            cli::cli_abort("All update() arguments must be named")
        }
        if (anyDuplicated(labels)) {
            duplicates <- unique(labels[duplicated(labels)])
            cli::cli_abort("Duplicate update() argument(s): {.and {duplicates}}")
        }
        unknown <- setdiff(labels, names(formals(combo_model)))
        if (length(unknown)) {
            cli::cli_abort("Unknown update() argument(s): {.and {unknown}}")
        }
    }
    arguments <- c(
        list(rank = object$rank, mean = object$mean),
        object$components,
        list(dose = object$dose, family = object$family)
    )
    # Single-bracket assignment preserves explicit NULL component replacements.
    arguments[names(replacements)] <- replacements
    do.call(combo_model, arguments)
}

#' @param x A `combo_model` object.
#' @param ... Reserved for future methods.
#' @return `x`, invisibly.
#' @rdname combo_model
#' @export
print.combo_model <- function(x, ...) {
    output <- cli::cli_format_method({
        cli::cli_text("<combo_model>")
        cli::cli_text("family: Gaussian(identity)")
        cli::cli_text("rank: {x[['rank']]}")
        mean <- x[["mean"]] # nolint: object_usage_linter.
        description <- switch(
            mean$type,
            empirical = "empirical observed-response mean",
            fixed = paste0("fixed at ", format(mean$value)),
            formula = paste(deparse(mean$formula), collapse = " ")
        )
        cli::cli_text("mean: {description}")
        if (x$dose$type == "categorical") {
            cli::cli_text("dose: categorical")
        } else {
            cli::cli_text("dose: nested (relative precision {x[['dose']][['precision']]})")
        }
        cli::cli_text("components:")
        for (name in names(x$components)) {
            component <- x$components[[name]]
            if (is.null(component)) {
                cli::cli_text("- {name}: disabled")
                next
            }
            structure <- structure_type(component[["structure"]]) # nolint: object_usage_linter.
            shrinkage <- component[["shrinkage"]][["type"]] # nolint: object_usage_linter.
            cli::cli_text("- {name}: structure {structure} | shrinkage {shrinkage}")
        }
    })
    writeLines(output)
    invisible(x)
}

#' Configure a Gaussian combination-model component
#'
#' Defines a conditionally Gaussian prior for one term in [combo_model()].
#' The model slot determines what the term means in the response predictor;
#' this function determines how its entity-level values are modeled.
#'
#' Every enabled component is a zero-centered residual deviation:
#'
#' ```
#' component value ~ Gaussian(0, precision_structure)
#' ```
#'
#' - `structure` describes which cell or treatment deviations are related and their
#'   relative precision. See [context_specifications].
#' - `shrinkage` controls the magnitude of deviations around zero;
#'   see [shrinkage_specifications].
#'
#' @param shrinkage Deviation-shrinkage specification; see [shrinkage_specifications].
#' @param structure Relationship structure among component entities. `NULL`
#'   selects IID; see [context_specifications].
#' @return A `combo_gaussian_component` specification, inheriting from the
#'   internal `combo_component` contract, for assignment in [combo_model()].
#' @export
combo_gaussian_component <- function(
    shrinkage,
    structure = NULL
) {
    if (missing(shrinkage) || !inherits(shrinkage, "param_shrinkage")) {
        cli::cli_abort("shrinkage must be an explicit shrinkage specification")
    }
    if (is.null(structure)) {
        # iid() is defined with the structural-prior code. A small independent
        # specification keeps this constructor source-order independent.
        structure <- structure(list(type = "iid"), class = "structural_prior")
    }
    if (!is_structure_spec(structure)) {
        cli::cli_abort("structure must be iid(), tree(), precision(), or gmrf()")
    }
    structure(
        list(
            structure = structure,
            shrinkage = shrinkage
        ),
        class = c("combo_gaussian_component", "combo_component")
    )
}

validate_component <- function(component, name, kind) {
    if (is.null(component)) {
        return(invisible(NULL))
    }
    if (!inherits(component, "combo_component")) {
        cli::cli_abort("{name} must be NULL or a combo_gaussian_component()")
    }
    q_type <- structure_type(component$structure)
    if (!q_type %in% c("iid", "tree", "precision", "gmrf")) {
        cli::cli_abort("{name} uses an unsupported structure: {q_type}")
    }
    shrinkage_type <- component$shrinkage$type
    allowed <- c(
        "fixed", "gamma", "global_half_cauchy", "local_half_cauchy",
        "horseshoe", "multiplicative_gamma"
    )
    if (!shrinkage_type %in% allowed) {
        cli::cli_abort("{name} uses unsupported shrinkage: {shrinkage_type}")
    }
    if (kind == "offset" && shrinkage_type == "multiplicative_gamma") {
        cli::cli_abort("{name} is scalar and cannot use multiplicative_gamma()")
    }
    if (q_type != "iid" && shrinkage_type %in% c("local_half_cauchy", "horseshoe")) {
        cli::cli_abort(
            "{name} uses a non-IID structure and cannot use entity-local shrinkage"
        )
    }
    invisible(component)
}

is_structure_spec <- function(x) {
    structure_specs <- c("structural_prior", "kernel_prior", "param_structure")
    for (spec in structure_specs) if (inherits(x, spec)) return(TRUE)
    FALSE
}

structure_type <- function(x) {
    if (!is_structure_spec(x) || is.null(x$type)) {
        cli::cli_abort("Invalid component structure")
    }
    as.character(x$type)
}

#' Context-structure specifications
#'
#' Specify how effects are related across cells or treatments:
#'
#' - `iid()` uses independent unit precision.
#' - `precision()` uses a supplied sparse positive-definite precision matrix.
#' - `gmrf()` adds a positive ridge to a sparse symmetric graph operator.
#' - `tree()` constructs an anchored Gaussian-tree precision from parent-child relationships.
#'
#' @param Q A named sparse positive-definite precision matrix.
#' @param operator A named sparse symmetric graph operator.
#' @param ridge A finite positive diagonal ridge added to `operator`.
#' @param node,parent Bare column names in `data` defining tree nodes and their
#'   parents. Missing parents denote anchored roots.
#' @param data A data frame containing the tree definition.
#' @param edge_length A finite positive scalar or a column in `data`.
#' @return A structural-prior specification for [combo_gaussian_component()].
#' @name context_specifications
NULL

#' @rdname context_specifications
#' @export
iid <- function() {
    structure(list(type = "iid"), class = "structural_prior")
}

#' @rdname context_specifications
#' @export
precision <- function(Q) {
    structure(
        list(type = "precision", Q = Q),
        class = "structural_prior"
    )
}

# Precision specifications for the Gaussian observation model

param_shrinkage <- function(type, ...) {
    structure(c(list(type = type), list(...)), class = "param_shrinkage")
}

#' Shrinkage specifications
#'
#' Configure the prior scale of zero-centered Gaussian component deviations.
#' All precision parameters use the shape-rate convention where a gamma
#' distribution is used.
#'
#' The general horseshoe prior applies to an IID deviation `u[n, d]` for entity
#' `n` and latent dimension `d`:
#'
#' ```
#' u[n, d] | tau[d], lambda[n, d] ~ Normal(0, tau[d]^2 * lambda[n, d]^2)
#' tau[d]       ~ half-Cauchy(0, global_scale)
#' lambda[n, d] ~ half-Cauchy(0, local_scale)
#' ```
#'
#' `tau[d]` controls the overall magnitude of dimension `d`, while
#' `lambda[n, d]` allows individual entity effects to escape that global
#' shrinkage. The constructors specialize this hierarchy as follows:
#'
#' - `horseshoe()` estimates both `tau[d]` and `lambda[n, d]`.
#' - `global_half_cauchy()` fixes every `lambda[n, d] = 1` and estimates only `tau[d]`.
#'   Because the scale is global, it can multiply any context precision matrix `Q`.
#' - `local_half_cauchy()` fixes every `tau[d] = 1` and estimates only `lambda[n, d]`.
#'   Entity-local scales are available only for IID contexts.
#'
#' Internally, the corresponding deviation precision is `1 / (tau[d]^2 * lambda[n, d]^2)`.
#' Thus smaller scales imply stronger shrinkage toward zero.
#'
#' `fixed_scale(precision = p)` fixes the precision multiplier:
#'
#' ```
#' u[, d] ~ Normal(0, inverse(p * Q))
#' ```
#'
#' `gamma_precision(shape = a, rate = b)` instead learns one precision
#' multiplier per component dimension:
#'
#' ```
#' rho[d] ~ Gamma(shape = a, rate = b)
#' u[, d] | rho[d] ~ Normal(0, inverse(rho[d] * Q))
#' ```
#'
#' For metadata-mean coefficients, `Q` is the identity matrix and the same
#' gamma prior regularizes the coefficient vector.
#'
#' `multiplicative_gamma(shape = a, rate = b)` is restricted to factor
#' components. It constructs dimension-specific precisions from cumulative
#' gamma increments:
#'
#' ```
#' delta[h] ~ Gamma(shape = a, rate = b)
#' rho[d] = product(delta[h], h = 1, ..., d)
#' u[, d] | rho[d] ~ Normal(0, inverse(rho[d] * Q))
#' ```
#'
#' Later dimensions inherit every earlier increment, encouraging progressively
#' stronger shrinkage and allowing unnecessary latent dimensions to collapse.
#'
#' @param precision Fixed positive precision multiplier `p`.
#' @param shape,rate Positive gamma shape `a` and rate `b`. For
#'   `multiplicative_gamma()`, the same pair is applied to every increment.
#' @param scale Positive half-Cauchy scale for the estimated global or local
#'   standard-deviation multiplier.
#' @param global_scale,local_scale Positive half-Cauchy scales for the
#'   horseshoe global and entity-local standard-deviation multipliers.
#' @return A shrinkage specification for [combo_gaussian_component()] or
#'   [gaussian_response()].
#' @name shrinkage_specifications
NULL

#' @rdname shrinkage_specifications
#' @export
fixed_scale <- function(precision = 1) {
    param_shrinkage(
        "fixed",
        precision = param_scalar_positive(precision, "precision")
    )
}

#' @rdname shrinkage_specifications
#' @export
gamma_precision <- function(shape = 1.1, rate = 1.1) {
    param_shrinkage(
        "gamma",
        shape = param_scalar_positive(shape, "shape"),
        rate = param_scalar_positive(rate, "rate")
    )
}

#' @rdname shrinkage_specifications
#' @export
global_half_cauchy <- function(scale = 1) {
    param_shrinkage(
        "global_half_cauchy",
        scale = param_scalar_positive(scale, "scale")
    )
}

#' @rdname shrinkage_specifications
#' @export
local_half_cauchy <- function(scale = 1) {
    param_shrinkage(
        "local_half_cauchy",
        scale = param_scalar_positive(scale, "scale")
    )
}

#' @rdname shrinkage_specifications
#' @export
horseshoe <- function(global_scale = 1, local_scale = 1) {
    param_shrinkage(
        "horseshoe",
        global_scale = param_scalar_positive(global_scale, "global_scale"),
        local_scale = param_scalar_positive(local_scale, "local_scale")
    )
}

#' @rdname shrinkage_specifications
#' @export
multiplicative_gamma <- function(shape = 2, rate = 1) {
    param_shrinkage(
        "multiplicative_gamma",
        shape = param_scalar_positive(shape, "shape"),
        rate = param_scalar_positive(rate, "rate")
    )
}

#' Gaussian response-family specification
#'
#' Defines the Gaussian identity-link observation model used by the Gibbs engine.
#'
#' @param link Response link. Only `"identity"` is currently supported.
#' @param precision A [gamma_precision()] specification for observation
#'   precision.
#' @return A `combo_family` specification for [combo_model()].
#' @export
gaussian_response <- function(
    link = "identity",
    precision = gamma_precision(shape = 1.1, rate = 1.1)
) {
    if (length(link) != 1L || !identical(as.character(link), "identity")) {
        cli::cli_abort("The Gibbs engine currently supports only the identity link")
    }
    if (!inherits(precision, "param_shrinkage") || !identical(precision$type, "gamma")) {
        cli::cli_abort("Gaussian observation precision must use gamma_precision()")
    }
    structure(
        list(name = "gaussian", link = "identity", precision = precision),
        class = "combo_family"
    )
}

# Mean specifications for the global regression

#' Global-mean specifications
#'
#' `empirical_mean()` fixes the intercept at the mean of the observed responses.
#' `fixed_mean()` fixes it at a supplied value. `formula_mean()` defines a
#' fixed-effect regression using observation, cell, and compound metadata. It
#' accepts one- or two-sided formulas; a left-hand side is ignored as an outcome
#' but remains excluded from `.` expansion.
#'
#' Formula terms retain the units and contrasts produced by [stats::model.matrix()];
#' no automatic centering or scaling is applied. The formula intercept, when
#' present, is learned with a flat prior. `beta_mean` and `beta_precision`
#' describe a proper Gaussian prior for the remaining coefficients.
#'
#' Non-reserved observation columns and covariates in the keyed cell and
#' compound metadata tables share one formula namespace. Their names must not
#' collide. Each term must use variables from only one source, so transformations
#' and within-source interactions are supported but cross-source interactions
#' are not. Cell terms enter once per observation; compound terms are evaluated
#' for each present drug and added. Use [stats::model.matrix()] on a
#' `combo_model` to inspect the exact coefficient names and raw-scale design.
#'
#' If the cell or compound terms, together with the formula intercept when
#' present, span every corresponding modeled entity while its offset component
#' remains active, [fit_combo()] warns that the fixed and residual effects will
#' be separated primarily by their priors. This check applies when there are
#' terms from that source and uses the rank of its entity-level design matrix.
#'
#' @param value A finite numeric scalar.
#' @param formula A one- or two-sided formula. Any left-hand side is ignored as
#'   an outcome but excluded from `.` expansion.
#' @param beta_mean A finite scalar or named numeric vector giving the prior
#'   mean of non-intercept coefficients.
#' @param beta_precision A positive scalar, named positive vector, or named
#'   symmetric positive-definite dense or sparse matrix giving the coefficient
#'   prior precision. Vector and matrix names refer to the non-intercept columns
#'   returned by [stats::model.matrix()].
#' @return A global-mean specification for [combo_model()].
#' @name mean_specifications
NULL

#' @rdname mean_specifications
#' @export
empirical_mean <- function() {
    structure(list(type = "empirical"), class = "combo_mean")
}

#' @rdname mean_specifications
#' @export
fixed_mean <- function(value) {
    if (!is_number(value)) {
        cli::cli_abort("value must be one finite number")
    }
    structure(
        list(type = "fixed", value = as.numeric(value)),
        class = "combo_mean"
    )
}

#' @rdname mean_specifications
#' @export
formula_mean <- function(
    formula = ~1,
    beta_mean = 0,
    beta_precision = 1
) {
    if (!inherits(formula, "formula") || !length(formula) %in% 2:3) {
        cli::cli_abort("formula must be a one- or two-sided formula")
    }
    valid_beta_mean <- is.numeric(beta_mean) && is.null(dim(beta_mean)) &&
        length(beta_mean) > 0L && !anyNA(beta_mean) && all(is.finite(beta_mean))
    if (!valid_beta_mean) {
        cli::cli_abort("beta_mean must be a finite numeric scalar or named vector")
    }
    if (!is.null(names(beta_mean)) && is_invalid_key(names(beta_mean))) {
        cli::cli_abort("beta_mean names must be unique and nonempty")
    }
    matrix_precision <- is.matrix(beta_precision) ||
        methods::is(beta_precision, "Matrix")
    vector_precision <- is.numeric(beta_precision) &&
        is.null(dim(beta_precision)) && length(beta_precision) > 0L
    if (!matrix_precision && !vector_precision) {
        cli::cli_abort(
            "beta_precision must be a positive scalar, named vector, or named matrix"
        )
    }
    if (vector_precision) {
        if (anyNA(beta_precision) || any(!is.finite(beta_precision)) ||
                any(beta_precision <= 0)) {
            cli::cli_abort("beta_precision values must be finite and positive")
        }
        if (length(beta_precision) > 1L &&
                (is.null(names(beta_precision)) ||
                    is_invalid_key(names(beta_precision)))) {
            cli::cli_abort(
                "A beta_precision vector must have unique, nonempty coefficient names"
            )
        }
    } else {
        precision <- as.matrix(beta_precision)
        row_names <- rownames(precision)
        column_names <- colnames(precision)
        if (!is.numeric(precision) || nrow(precision) != ncol(precision) ||
                anyNA(precision) || any(!is.finite(precision))) {
            cli::cli_abort("beta_precision matrix must be square and finite")
        }
        if (is.null(row_names) || is.null(column_names) ||
                is_invalid_key(row_names) || !identical(row_names, column_names)) {
            cli::cli_abort(
                "beta_precision matrix must have identical unique row and column names"
            )
        }
        if (!isTRUE(all.equal(precision, t(precision), tolerance = 1e-12))) {
            cli::cli_abort("beta_precision matrix must be symmetric")
        }
        tryCatch(
            chol(precision),
            error = function(error) {
                cli::cli_abort("beta_precision matrix must be positive definite")
            }
        )
    }
    structure(
        list(
            type = "formula",
            formula = formula,
            beta_mean = beta_mean,
            beta_precision = beta_precision
        ),
        class = "combo_mean"
    )
}

# Dose specifications for the combination model.

#' Dose specifications
#'
#' `categorical()` models each observed compound-dose treatment independently.
#' `nested()` connects doses to a latent compound parent with the supplied
#' relative precision.
#'
#' @param precision A finite positive relative precision for dose-to-compound
#'   edges.
#' @return A dose specification for [combo_model()].
#' @name dose_specifications
NULL

#' @rdname dose_specifications
#' @export
categorical <- function() {
    structure(list(type = "categorical"), class = "combo_dose")
}

#' @rdname dose_specifications
#' @export
nested <- function(precision = 1) {
    structure(
        list(
            type = "nested",
            precision = param_scalar_positive(precision, "precision")
        ),
        class = "combo_dose"
    )
}
