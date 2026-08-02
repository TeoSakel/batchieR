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
        stop(label, " must be one finite positive number", call. = FALSE)
    }
    as.numeric(x)
}

is_positive_integer <- function(x) {
    length(x) == 1L && !is.na(x) && is.finite(x) && x >= 1 && x == as.integer(x)
}

param_positive_integer <- function(x, label) {
    if (!is_positive_integer(x)) {
        stop(label, " must be one positive integer", call. = FALSE)
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
        stop("lower must be less than or equal to upper")
    }
    pmax(pmin(x, upper), lower)  # keep shape of x
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
            warning("Numerical instability in ", label, "; retaining previous value: ",
                conditionMessage(error),
                call. = FALSE
            )
            fallback
        }
    )
}
