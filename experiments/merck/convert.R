#!/usr/bin/env Rscript

merck_convert_screen <- function(path, require_observed = FALSE) {
    merck_require("hdf5r")
    if (!file.exists(path)) stop("HDF5 screen does not exist: ", path, call. = FALSE)
    file <- hdf5r::H5File$new(path, mode = "r")
    on.exit(file$close_all())
    required <- c(
        "treatment_names", "treatment_doses", "observations",
        "observation_mask", "sample_names", "plate_names"
    )
    missing <- setdiff(required, names(file))
    if (length(missing)) {
        stop("HDF5 screen is missing: ", paste(missing, collapse = ", "), call. = FALSE)
    }
    observations <- as.numeric(file[["observations"]][])
    n <- length(observations)
    normalize_treatments <- function(value, label) {
        dimensions <- dim(value)
        if (identical(dimensions, c(2L, n))) return(t(value))
        if (identical(dimensions, c(n, 2L))) return(value)
        stop(label, " must have dimensions n_experiments x 2", call. = FALSE)
    }
    treatment_names <- normalize_treatments(
        file[["treatment_names"]][,], "treatment_names"
    )
    treatment_doses <- normalize_treatments(
        file[["treatment_doses"]][,], "treatment_doses"
    )
    sample_names <- as.character(file[["sample_names"]][])
    plate_names <- as.character(file[["plate_names"]][])
    mask <- as.logical(file[["observation_mask"]][])
    if (any(lengths(list(sample_names, plate_names, mask)) != n)) {
        stop("HDF5 screen arrays have inconsistent lengths", call. = FALSE)
    }
    if (any(!is.finite(observations))) {
        stop("HDF5 observations must be finite", call. = FALSE)
    }
    if (any(observations < 0 | observations > 1)) {
        stop("HDF5 observations must lie in [0, 1]", call. = FALSE)
    }
    if (require_observed && !all(mask)) {
        stop("Holdout screen must be fully observed", call. = FALSE)
    }
    plate_states <- split(mask, plate_names)
    if (any(vapply(plate_states, function(value) length(unique(value)) != 1L, logical(1)))) {
        stop("Every plate must be entirely observed or entirely masked", call. = FALSE)
    }
    control <- as.character(file$attr_open("control_treatment_name")$read())
    viability <- pmin(pmax(observations, 0.01), 0.99)
    truth <- stats::qlogis(viability)
    result <- data.frame(
        cell = sample_names,
        drug_1 = as.character(treatment_names[, 1L]),
        dose_1 = as.numeric(treatment_doses[, 1L]),
        drug_2 = as.character(treatment_names[, 2L]),
        dose_2 = as.numeric(treatment_doses[, 2L]),
        response = ifelse(mask, truth, NA_real_),
        plate = plate_names,
        truth = truth,
        viability = viability,
        observed = mask,
        stringsAsFactors = FALSE
    )
    for (position in 1:2) {
        drug <- paste0("drug_", position)
        dose <- paste0("dose_", position)
        is_control <- result[[drug]] == control | result[[dose]] <= 0
        if (anyNA(is_control)) stop("Treatment names and doses must be nonmissing", call. = FALSE)
        result[[drug]][is_control] <- NA_character_
        result[[dose]][is_control] <- NA_real_
    }
    list(data = result, control = control)
}

merck_convert <- function(training_path, holdout_path, output_path, manifest_path = NULL) {
    training <- merck_convert_screen(training_path)
    holdout <- merck_convert_screen(holdout_path, require_observed = TRUE)
    if (!identical(training$control, holdout$control)) {
        stop("Training and holdout controls differ", call. = FALSE)
    }
    initial <- unique(training$data$plate[training$data$observed])
    if (length(initial) != 1L) stop("Expected exactly one initially observed plate", call. = FALSE)
    preparation <- if (!is.null(manifest_path)) {
        merck_require("jsonlite")
        jsonlite::read_json(manifest_path, simplifyVector = TRUE)
    } else {
        list()
    }
    treatment_pairs <- unique(rbind(
        training$data[c("drug_1", "dose_1")],
        setNames(training$data[c("drug_2", "dose_2")], c("drug_1", "dose_1"))
    ))
    treatment_pairs <- treatment_pairs[!is.na(treatment_pairs$drug_1), , drop = FALSE]
    dataset <- list(
        training = training$data,
        holdout = holdout$data,
        metadata = list(
            format_version = 1L,
            control_treatment_name = training$control,
            initial_plate = initial,
            n_training = nrow(training$data),
            n_holdout = nrow(holdout$data),
            n_plates = length(unique(training$data$plate)),
            n_initial = sum(training$data$observed),
            n_cells = length(unique(training$data$cell)),
            n_treatments = nrow(treatment_pairs),
            preparation = preparation
        )
    )
    if (length(preparation) &&
            exists("UPSTREAM_COMMIT", inherits = TRUE) &&
            identical(preparation$upstream_commit, UPSTREAM_COMMIT)) {
        actual <- unlist(dataset$metadata[c(
            "n_training", "n_holdout", "n_plates", "n_initial", "n_treatments"
        )])
        expected <- c(279832L, 31545L, 702L, 133L, 168L)
        if (!identical(as.integer(actual), expected)) {
            stop("Pinned Merck preparation produced unexpected dimensions", call. = FALSE)
        }
    }
    merck_atomic_save_rds(dataset, output_path)
    message("Wrote ", normalizePath(output_path))
    invisible(dataset)
}

merck_convert_main <- function() {
    script_dir <- merck_script_dir()
    args <- merck_parse_args(commandArgs(trailingOnly = TRUE))
    required <- c("training", "holdout", "output")
    missing <- required[!required %in% names(args)]
    if (length(missing)) stop("Missing --", paste(missing, collapse = ", --"), call. = FALSE)
    merck_convert(args$training, args$holdout, args$output, args$manifest)
}

if (sys.nframe() == 0L) {
    command <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", command, value = TRUE)
    script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
    source(file.path(script_dir, "common.R"))
    source(file.path(script_dir, "config.R"))
    merck_convert_main()
}
