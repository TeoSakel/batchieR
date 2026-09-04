# Generic Gibbs engine for compiled combination models.

init_gibbs_state <- function(compiled) {
    if (!inherits(compiled, "compiled_combo_model")) {
        cli::cli_abort("compiled must be a compiled_combo_model")
    }
    components <- lapply(compiled$components, init_component_state)
    structure(
        list(
            compiled = compiled,
            y = compiled$response,
            cell = compiled$cell,
            treatment_1 = compiled$treatment_1,
            treatment_2 = compiled$treatment_2,
            alpha = compiled$alpha,
            precision = 100,
            components = components,
            Mu = rep(compiled$alpha, length(compiled$response)),
            n_steps = 0L,
            last_rmse = NA_real_
        ),
        class = "combo_gibbs_state"
    )
}

gibbs_step <- function(state) {
    state$n_steps <- state$n_steps + 1L
    state$Mu <- gibbs_reconstruct_mean(state)
    state <- gibbs_update_cell_offset(state)
    state <- gibbs_update_treatment_offset(state)
    state <- gibbs_update_cell_factors(state)
    state <- gibbs_update_interactions(state)
    state <- gibbs_update_treatment_factors(state)
    state <- gibbs_update_hyperparameters(state)
    gibbs_update_obs_precision(state)
}

gibbs_reconstruct_mean <- function(state) {
    n <- length(state$y)
    mu <- rep(state$alpha, n)
    W0 <- state$components$cell_offset
    if (!is.null(W0)) mu <- mu + W0$values[state$cell, 1L]
    V0 <- state$components$treatment_offset
    if (!is.null(V0)) {
        mu <- mu +
            combo_component_rows(V0, state$treatment_1, 1L)[, 1L] +
            combo_component_rows(V0, state$treatment_2, 1L)[, 1L]
    }
    W <- state$components$cell_factors
    if (!is.null(W)) {
        w <- W$values[state$cell, , drop = FALSE]
        V1 <- state$components$treatment_main_factors
        if (!is.null(V1)) {
            v1 <- combo_component_rows(V1, state$treatment_1) +
                combo_component_rows(V1, state$treatment_2)
            mu <- mu + rowSums(w * v1)
        }
        V2 <- state$components$treatment_interaction_factors
        if (!is.null(V2)) {
            v2 <- combo_component_rows(V2, state$treatment_1) *
                combo_component_rows(V2, state$treatment_2)
            mu <- mu + rowSums(w * v2)
        }
    }
    mu
}

gibbs_update_cell_offset <- function(state) {
    component <- state$components$cell_offset
    if (is.null(component)) {
        return(state)
    }
    component <- combo_update_latent(component)
    for (entity in seq_len(nrow(component$values))) {
        idx <- which(state$cell == entity)
        old <- component$values[entity, 1L]
        prior <- combo_entity_prior(component, entity)
        precision <- prior$precision
        linear <- prior$linear
        if (length(idx)) {
            residual <- state$y[idx] - state$Mu[idx] + old
            precision <- precision + state$precision * length(idx)
            linear <- linear + state$precision * sum(residual)
        }
        value <- stats::rnorm(1L, linear / precision, 1 / sqrt(precision))
        component <- combo_store_entity(component, entity, value, prior)
        if (length(idx)) {
            state$Mu[idx] <- state$Mu[idx] + value - old
        }
    }
    state$components$cell_offset <- component
    state
}

gibbs_update_treatment_offset <- function(state) {
    component <- state$components$treatment_offset
    if (is.null(component)) {
        return(state)
    }
    component <- combo_update_latent(component)
    for (entity in seq_len(nrow(component$values))) {
        coefficient <- combo_treatment_coefficient(state, entity)
        idx <- which(coefficient > 0L)
        old <- component$values[entity, 1L]
        prior <- combo_entity_prior(component, entity)
        precision <- prior$precision
        linear <- prior$linear
        if (length(idx)) {
            x <- coefficient[idx]
            residual <- state$y[idx] - state$Mu[idx] + x * old
            precision <- precision + state$precision * sum(x^2)
            linear <- linear + state$precision * sum(x * residual)
        }
        value <- stats::rnorm(1L, linear / precision, 1 / sqrt(precision))
        component <- combo_store_entity(component, entity, value, prior)
        if (length(idx)) {
            state$Mu[idx] <- state$Mu[idx] + x * (value - old)
        }
    }
    state$components$treatment_offset <- component
    state
}

gibbs_update_cell_factors <- function(state) {
    component <- state$components$cell_factors
    if (is.null(component)) {
        return(state)
    }
    component <- combo_update_latent(component)
    rank <- ncol(component$values)
    for (entity in seq_len(nrow(component$values))) {
        idx <- which(state$cell == entity)
        old <- component$values[entity, ]
        prior <- combo_entity_prior(component, entity)
        precision <- diag(prior$precision, rank)
        linear <- prior$linear
        if (length(idx)) {
            x <- combo_factor_design_for_cell(state, idx)
            old_contribution <- as.numeric(x %*% old)
            residual <- state$y[idx] - state$Mu[idx] + old_contribution
            precision <- precision + state$precision * crossprod(x)
            linear <- linear +
                state$precision * as.numeric(crossprod(x, residual))
        }
        value <- rmvnorm_safe(
            precision,
            linear,
            old,
            "cell_factors update"
        )
        component <- combo_store_entity(component, entity, value, prior)
        if (length(idx)) {
            state$Mu[idx] <- state$Mu[idx] +
                as.numeric(x %*% value) - old_contribution
        }
    }
    state$components$cell_factors <- component
    state
}

gibbs_update_treatment_factors <- function(state) {
    component <- state$components$treatment_main_factors
    if (is.null(component)) {
        return(state)
    }
    component <- combo_update_latent(component)
    W <- state$components$cell_factors
    rank <- ncol(component$values)
    for (entity in seq_len(nrow(component$values))) {
        coefficient <- combo_treatment_coefficient(state, entity)
        idx <- which(coefficient > 0L)
        old <- component$values[entity, ]
        prior <- combo_entity_prior(component, entity)
        precision <- diag(prior$precision, rank)
        linear <- prior$linear
        if (length(idx)) {
            x <- W$values[state$cell[idx], , drop = FALSE] * coefficient[idx]
            old_contribution <- as.numeric(x %*% old)
            residual <- state$y[idx] - state$Mu[idx] + old_contribution
            precision <- precision + state$precision * crossprod(x)
            linear <- linear +
                state$precision * as.numeric(crossprod(x, residual))
        }
        value <- rmvnorm_safe(
            precision,
            linear,
            old,
            "treatment_main_factors update"
        )
        component <- combo_store_entity(component, entity, value, prior)
        if (length(idx)) {
            state$Mu[idx] <- state$Mu[idx] +
                as.numeric(x %*% value) - old_contribution
        }
    }
    state$components$treatment_main_factors <- component
    state
}

gibbs_update_interactions <- function(state) {
    component <- state$components$treatment_interaction_factors
    if (is.null(component)) {
        return(state)
    }
    component <- combo_update_latent(component)
    W <- state$components$cell_factors
    rank <- ncol(component$values)
    for (entity in seq_len(nrow(component$values))) {
        idx1 <- which(state$treatment_1 == entity)
        idx2 <- which(state$treatment_2 == entity)
        idx <- c(idx1, idx2)
        old <- component$values[entity, ]
        prior <- combo_entity_prior(component, entity)
        precision <- diag(prior$precision, rank)
        linear <- prior$linear
        if (length(idx)) {
            x1 <- if (length(idx1)) {
                W$values[state$cell[idx1], , drop = FALSE] *
                    combo_component_rows(
                        component,
                        state$treatment_2[idx1],
                        rank
                    )
            } else {
                matrix(numeric(), nrow = 0L, ncol = rank)
            }
            x2 <- if (length(idx2)) {
                W$values[state$cell[idx2], , drop = FALSE] *
                    combo_component_rows(
                        component,
                        state$treatment_1[idx2],
                        rank
                    )
            } else {
                matrix(numeric(), nrow = 0L, ncol = rank)
            }
            x <- rbind(x1, x2)
            old_contribution <- as.numeric(x %*% old)
            residual <- c(
                state$y[idx1] - state$Mu[idx1],
                state$y[idx2] - state$Mu[idx2]
            ) + old_contribution
            precision <- precision + state$precision * crossprod(x)
            linear <- linear +
                state$precision * as.numeric(crossprod(x, residual))
        }
        value <- rmvnorm_safe(
            precision,
            linear,
            old,
            "treatment_interaction_factors update"
        )
        component <- combo_store_entity(component, entity, value, prior)
        if (length(idx)) {
            state$Mu[idx] <- state$Mu[idx] +
                as.numeric(x %*% value) - old_contribution
        }
    }
    state$components$treatment_interaction_factors <- component
    state
}

gibbs_update_hyperparameters <- function(state) {
    for (name in names(state$components)) {
        component <- state$components[[name]]
        component <- gibbs_update_beta(component)
        component <- gibbs_update_shrinkage(component, length(state$y))
        component <- gibbs_update_beta_precision(component)
        state$components[[name]] <- component
    }
    state
}

gibbs_update_beta <- function(component) {
    if (is.null(component) || !length(component$beta)) {
        return(component)
    }
    compiled <- component$compiled
    structure <- compiled$structure
    X <- compiled$mean$X
    scale <- structure$modeled_scale
    A <- matrix(
        0,
        nrow = length(structure$nodes),
        ncol = ncol(X)
    )
    A[structure$modeled_index, ] <- X * scale
    z <- component$raw
    z[structure$modeled_index, 1L] <-
        z[structure$modeled_index, 1L] +
        scale * as.numeric(X %*% component$beta)
    prior_precision <- component$beta_precision$precision
    if (!is.null(component$shrinkage$local)) {
        weights <- component$shrinkage$global[1L] *
            component$shrinkage$local[, 1L]
        A_leaf <- A[structure$modeled_index, , drop = FALSE]
        z_leaf <- z[structure$modeled_index, 1L]
        precision <- crossprod(A_leaf, A_leaf * weights) +
            diag(prior_precision, ncol(A_leaf))
        linear <- as.numeric(crossprod(A_leaf, weights * z_leaf))
    } else {
        lambda <- component$shrinkage$global[1L]
        precision <- lambda *
            crossprod(A, as.matrix(structure$Q %*% A)) +
            diag(prior_precision, ncol(A))
        linear <- lambda * as.numeric(
            crossprod(A, as.numeric(structure$Q %*% z[, 1L]))
        )
    }
    component$beta <- stats::setNames(
        rmvnorm_safe(
            precision,
            linear,
            component$beta,
            paste0(compiled$name, " mean coefficient update")
        ),
        names(component$beta)
    )
    component$mean_value <- as.numeric(
        component$compiled$mean$X %*% component$beta
    )
    component <- combo_sync_component_raw(component)
    component
}

gibbs_update_shrinkage <- function(component, n_observations) {
    if (is.null(component)) {
        return(component)
    }
    state <- component$shrinkage
    spec <- state$spec
    type <- state$type
    n_nodes <- nrow(component$raw)
    quadratic <- combo_component_quadratic(component)
    if (type == "fixed") {
        return(component)
    }
    if (type == "gamma") {
        state$global <- stats::rgamma(
            length(state$global),
            shape = spec$shape + 0.5 * n_nodes,
            rate = spec$rate + 0.5 * quadratic + 1e-3
        )
    } else if (type == "global_half_cauchy") {
        state$global_aux <- stats::rgamma(
            length(state$global),
            shape = 1,
            rate = 1 / spec$scale^2 + state$global
        )
        state$global <- stats::rgamma(
            length(state$global),
            shape = 0.5 * (1 + n_nodes),
            rate = state$global_aux + 0.5 * quadratic + 1e-3
        )
    } else if (type %in% c("local_half_cauchy", "horseshoe")) {
        residual <- component$raw[
            component$compiled$structure$modeled_index,
            ,
            drop = FALSE
        ]
        local_scale <- if (type == "horseshoe") spec$local_scale else spec$scale
        state$local_aux <- matrix(
            stats::rgamma(
                length(state$local),
                shape = 1,
                rate = 1 / local_scale^2 + as.numeric(state$local)
            ),
            nrow = nrow(state$local),
            ncol = ncol(state$local)
        )
        global_matrix <- matrix(
            state$global,
            nrow = nrow(state$local),
            ncol = ncol(state$local),
            byrow = TRUE
        )
        state$local <- matrix(
            stats::rgamma(
                length(state$local),
                shape = 1,
                rate = as.numeric(
                    state$local_aux + 0.5 * global_matrix * residual^2 + 1e-3
                )
            ),
            nrow = nrow(state$local),
            ncol = ncol(state$local)
        )
        if (type == "horseshoe") {
            state$global_aux <- stats::rgamma(
                length(state$global),
                shape = 1,
                rate = 1 / spec$global_scale^2 + state$global
            )
            state$global <- stats::rgamma(
                length(state$global),
                shape = 0.5 * (1 + nrow(state$local)),
                rate = state$global_aux +
                    0.5 * colSums(state$local * residual^2) + 1e-3
            )
        }
    } else if (type == "multiplicative_gamma") {
        D <- length(state$delta)
        for (dimension in seq_len(D)) {
            tau_without <- cumprod(state$delta)[dimension:D] / state$delta[dimension]
            shape <- spec$shape[dimension] + 0.5 * n_nodes * (D - dimension + 1L)
            rate <- spec$rate[dimension] + 0.5 * sum(quadratic[dimension:D] * tau_without) + 1e-3
            state$delta[dimension] <- stats::rgamma(
                1L,
                shape = shape,
                rate = rate
            )
        }
        state$global <- cumprod(state$delta)
    }
    state$global <- combo_clip_precision(
        state$global,
        n_observations
    )
    if (!is.null(state$local)) {
        local_lower <- matrix(
            1 / sqrt(1 + component$compiled$exposure_count),
            nrow = nrow(state$local),
            ncol = ncol(state$local)
        )
        state$local <- matrix(
            combo_clip_precision(
                state$local,
                n_observations,
                lower = as.numeric(local_lower)
            ),
            nrow = nrow(state$local),
            ncol = ncol(state$local)
        )
    }
    component$shrinkage <- state
    component
}

gibbs_update_beta_precision <- function(component) {
    if (is.null(component) || !length(component$beta) ||
            component$beta_precision$spec$type == "fixed") {
        return(component)
    }
    spec <- component$beta_precision$spec
    component$beta_precision$precision <- stats::rgamma(
        1L,
        shape = spec$shape + 0.5 * length(component$beta),
        rate = spec$rate + 0.5 * sum(component$beta^2)
    )
    component
}

gibbs_update_obs_precision <- function(state, eps = 1e-3) {
    n_obs <- length(state$y)
    rss <- as.numeric(crossprod(state$y - state$Mu))
    prior <- state$compiled$observation_prior
    shape_post <- prior$shape + 0.5 * n_obs
    rate_post <- prior$rate + 0.5 * rss
    state$precision <- combo_clip_precision(
        stats::rgamma(1L, shape = shape_post, rate = rate_post),
        n_obs,
        lower = eps
    )
    state$last_rmse <- sqrt(rss / n_obs)
    state
}

gibbs_snapshot <- function(state) {
    values <- lapply(state$components, function(component) {
        if (is.null(component)) {
            return(NULL)
        }
        list(
            value = component[["values"]],
            beta = component[["beta"]],
            beta_precision = if (is.null(component[["beta_precision"]])) {
                NULL
            } else {
                component[["beta_precision"]][["precision"]]
            },
            global_precision = component[["shrinkage"]][["global"]],
            local_precision = component[["shrinkage"]][["local"]],
            raw = component[["raw"]]
        )
    })
    list(
        alpha = state$alpha,
        precision = state$precision,
        components = values,
        Mu = state$Mu,
        last_rmse = state$last_rmse,
        n_steps = state$n_steps
    )
}


# Auxiary functions ---------------------------------------

combo_entity_prior <- function(component, entity) {
    structure <- component$compiled$structure
    if (structure$iid) {
        lambda <- component$shrinkage$global
        if (!is.null(component$shrinkage$local)) {
            lambda <- lambda * component$shrinkage$local[entity, ]
        }
        mean <- component$mean_value[entity]
        return(list(
            precision = lambda,
            linear = lambda * mean,
            node = entity,
            scale = 1,
            mean = mean
        ))
    }
    node <- structure$modeled_index[entity]
    scale <- structure$modeled_scale[entity]
    lambda <- component$shrinkage$global
    if (!is.null(component$shrinkage$local)) {
        lambda <- lambda * component$shrinkage$local[entity, ]
    }
    precision <- lambda * structure$diagonal[node] * scale^2
    linear_delta <- -lambda * combo_neighbor_sum(component, node) * scale
    mean <- component[["mean_value"]][entity]
    list(
        precision = precision,
        linear = linear_delta + precision * mean,
        node = node,
        scale = scale,
        mean = mean
    )
}

combo_update_latent <- function(component) {
    if (is.null(component) || !length(component$compiled$structure$latent_index)) {
        return(component)
    }
    structure <- component$compiled$structure
    lambda <- component$shrinkage$global
    for (node in structure$latent_index) {
        diagonal <- structure$diagonal[node]
        component$raw[node, ] <- stats::rnorm(
            n = length(lambda),
            mean = -combo_neighbor_sum(component, node) / diagonal,
            sd = 1 / sqrt(lambda * diagonal)
        )
    }
    component
}

combo_sync_component_raw <- function(component) {
    idx <- component[["compiled"]][["structure"]][["modeled_index"]]
    scale <- component[["compiled"]][["structure"]][["modeled_scale"]]
    residual <- component[["values"]] - component[["mean_value"]]
    component[["raw"]][idx, ] <- residual * scale
    component
}

combo_factor_design_for_cell <- function(state, idx) {
    rank <- state$compiled$model$rank
    result <- matrix(0, nrow = length(idx), ncol = rank)
    V1 <- state$components$treatment_main_factors
    if (!is.null(V1)) {
        result <- result +
            combo_component_rows(V1, state$treatment_1[idx], rank) +
            combo_component_rows(V1, state$treatment_2[idx], rank)
    }
    V2 <- state$components$treatment_interaction_factors
    if (!is.null(V2)) {
        result <- result +
            combo_component_rows(V2, state$treatment_1[idx], rank) *
            combo_component_rows(V2, state$treatment_2[idx], rank)
    }
    result
}

combo_component_quadratic <- function(component) {
    if (component$compiled$structure$iid) {
        return(colSums(component$raw^2))
    }
    Q_raw <- component$compiled$structure$Q %*% component$raw
    pmax(colSums(component$raw * as.matrix(Q_raw)), 0)
}

combo_component_rows <- function(component, index, n_dimensions = NULL) {
    if (is.null(n_dimensions)) {
        n_dimensions <- if (is.null(component)) 1L else ncol(component[["values"]])
    }
    result <- matrix(0, nrow = length(index), ncol = n_dimensions)
    if (is.null(component)) return(result)
    selected <- index > 0L
    if (any(selected)) {
        result[selected, ] <- component$values[index[selected], , drop = FALSE]
    }
    result
}

combo_neighbor_sum <- function(component, node) {
    structure <- component$compiled$structure
    neighbors <- structure$neighbors[[node]]
    if (!length(neighbors)) return(numeric(ncol(component$raw)))
    colSums(
        component$raw[neighbors, , drop = FALSE] * structure$coefficients[[node]]
    )
}

combo_treatment_coefficient <- function(state, entity) {
    as.integer(state$treatment_1 == entity) + as.integer(state$treatment_2 == entity)
}

combo_store_entity <- function(component, entity, value, prior) {
    component[["values"]][entity, ] <- value
    component[["raw"]][prior$node, ] <- prior$scale * (value - prior$mean)
    component
}

combo_clip_precision <- function(x, n_obs, upper = 1e6, lower = NULL) {
    lower <- lower %||% 1 / sqrt(1 + n_obs)
    clip(as.numeric(x), lower = lower, upper = upper)
}
