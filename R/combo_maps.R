#' Look up fitted cell, treatment, combination, and observation indices
#'
#' These maps connect posterior variable indices to the fitted design without
#' extracting draws. Cells and treatments retain their compiled ordering; a
#' treatment is a drug-dose pair. Combinations identify a cell and an unordered
#' pair of treatments, including singles. Combo IDs follow first occurrence in
#' the original input, including rows with missing responses. Reversed treatment
#' order and replicate observations share a combo ID.
#'
#' Treatment IDs within a combination are in ascending order, with absent
#' treatments last and represented by `NA`. Observation row IDs are original
#' input positions, independent of data-frame row names. Multiple observations
#' can share a combination while having different observation-level covariates
#' and expected responses. Their factor contributions are shared.
#'
#' The `entity_index` column of [parameter_map()] refers to these IDs for
#' `entity_type` equal to `"cell"`, `"treatment"`, or `"combo"`.
#'
#' @param fit A `combo_fit` object.
#' @return A data frame:
#'   * `cell_map()`: `cell_id`, `cell`.
#'   * `treatment_map()`: `treatment_id`, `drug`, `dose`.
#'   * `combo_map()`: `combo_id`, `cell_id`, `cell`, `treatment_1_id`,
#'     `drug_1`, `dose_1`, `treatment_2_id`, `drug_2`, `dose_2`.
#'   * `row_map()`: `row_id`, `combo_id`, `observed`. `observed` indicates
#'     whether that row contributed a response to the fitted likelihood.
#' @examples
#' \dontrun{
#' combos <- combo_map(fit)
#' rows <- row_map(fit)
#' selected <- subset(combos, cell == "A" & !is.na(treatment_2_id))
#' map <- parameter_map(fit)
#' variables <- subset(map, component == "interaction_effect" &
#'     entity_index %in% selected$combo_id)$variable
#' summary(fit, select = variables)
#' }
#' @name combo_index_maps
NULL

#' @rdname combo_index_maps
#' @export
cell_map <- function(fit) {
    validate_combo_fit(fit)
    cells <- fit$compiled$cells
    data.frame(cell_id = seq_along(cells), cell = cells)
}

#' @rdname combo_index_maps
#' @export
treatment_map <- function(fit) {
    validate_combo_fit(fit)
    treatments <- fit$compiled$treatments
    data.frame(treatment_id = seq_len(nrow(treatments)),
               drug = treatments$drug,
               dose = treatments$dose)
}

#' @rdname combo_index_maps
#' @export
combo_map <- function(fit) {
    validate_combo_fit(fit)
    combos <- combo_index_map(fit)[["combos"]]
    treatments <- fit$compiled$treatments
    a <- combos$treatment_1
    b <- combos$treatment_2
    a[a == 0L] <- NA_integer_
    b[b == 0L] <- NA_integer_
    data.frame(
        combo_id = seq_len(nrow(combos)),
        cell_id = combos$cell,
        cell = fit$compiled$cells[combos$cell],
        treatment_1_id = a,
        drug_1 = treatments$drug[a],
        dose_1 = treatments$dose[a],
        treatment_2_id = b,
        drug_2 = treatments$drug[b],
        dose_2 = treatments$dose[b]
    )
}

#' @rdname combo_index_maps
#' @export
row_map <- function(fit) {
    validate_combo_fit(fit)
    ids <- combo_index_map(fit)[["row_combo"]]
    rows <- seq_along(ids)
    data.frame(
        row_id = rows,
        combo_id = ids,
        observed = rows %in% fit$compiled$observed_rows
    )
}

# Keep the sampler's zero sentinel internally; expose NA only in public maps.
combo_index_map <- function(fit) {
    a <- fit$compiled$all_treatment_1
    b <- fit$compiled$all_treatment_2
    first <- ifelse(a == 0L | b == 0L, pmax(a, b), pmin(a, b))
    second <- ifelse(a == 0L | b == 0L, 0L, pmax(a, b))
    keys <- paste(fit$compiled$all_cell, first, second, sep = ":")
    keep <- !duplicated(keys)
    combos <- data.frame(
        cell = fit$compiled$all_cell[keep],
        treatment_1 = first[keep],
        treatment_2 = second[keep]
    )
    list(combos = combos, row_combo = match(keys, keys[keep])
    )
}
