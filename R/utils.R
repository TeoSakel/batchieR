# useful functions used in multiple places

`%||%` <- function(a, b) if (is.null(a)) b else a

sparse_values <- function(A) {
    if ("x" %in% methods::slotNames(A)) {
        return(as.numeric(methods::slot(A, "x")))
    }
    rep(1, Matrix::nnzero(A))
}

# checkers
is_invalid_key <- function(x) {
    anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)
}

param_scalar_positive <- function(x, label) {
    if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0) {
        cli::cli_abort("{label} must be one finite positive number")
    }
    as.numeric(x)
}

is_integer <- function(x) {
    is.numeric(x) && !is.na(x) && is.finite(x) && x == as.integer(x)
}

is_positive_integer <- function(x) {
    is_integer(x) && x > 0
}

is_nonnegative_integer <- function(x) {
    is_integer(x) && x >= 0
}

param_positive_integer <- function(x, label) {
    if (!is_positive_integer(x)) {
        cli::cli_abort("{label} must be one positive integer")
    }
    as.integer(x)
}

formula_has_terms <- function(formula) {
    terms <- stats::terms(formula)
    length(attr(terms, "term.labels")) > 0L || attr(terms, "intercept") != 0L
}

formula_is_zero <- function(formula) {
    !formula_has_terms(formula)
}


clip <- function(x, lower = -Inf, upper = Inf) {
    if (any(lower > upper)) {
        cli::cli_abort("lower must be less than or equal to upper")
    }
    pmax(pmin(x, upper), lower)  # keep shape of x
}

rhcauchy <- function(n, scale = 1) {
    abs(stats::rcauchy(n, location = 0, scale = scale))
}

batchieR_rmvnorm <- function(precision, mu_part = NULL) {
    upper <- chol(precision)
    value <- backsolve(upper, stats::rnorm(nrow(precision)))
    if (!is.null(mu_part)) {
        value <- value + as.numeric(chol2inv(upper) %*% mu_part)
    }
    as.numeric(value)
}

rmvnorm_safe <- function(precision, mu_part, fallback, label) {
    tryCatch(
        batchieR_rmvnorm(precision, mu_part),
        error = function(error) {
            cli::cli_warn("Numerical instability in {label}; retaining previous value: {conditionMessage(error)}")
            fallback
        }
    )
}
