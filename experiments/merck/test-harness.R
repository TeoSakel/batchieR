#!/usr/bin/env Rscript

command <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command, value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", file_arg[[1L]])))
source(file.path(script_dir, "common.R"))
source(file.path(script_dir, "convert.R"))
merck_require("hdf5r")

write_screen <- function(path, observations, mask, plates) {
    file <- hdf5r::H5File$new(path, mode = "w")
    on.exit(file$close_all())
    n <- length(observations)
    file[["treatment_names"]] <- matrix(
        rep(c("drug", "control"), n), nrow = 2L
    )
    file[["treatment_doses"]] <- matrix(
        rep(c(0.35, 0), n), nrow = 2L
    )
    file[["observations"]] <- observations
    file[["observation_mask"]] <- mask
    file[["sample_names"]] <- rep("cell", n)
    file[["plate_names"]] <- plates
    invisible(file$create_attr("control_treatment_name", robj = "control"))
}

directory <- tempfile("batchieR-merck-test-")
dir.create(directory)
training_path <- file.path(directory, "training.h5")
holdout_path <- file.path(directory, "holdout.h5")
output_path <- file.path(directory, "converted.rds")
write_screen(
    training_path,
    observations = c(0, 0.2, 0.8, 1),
    mask = c(TRUE, TRUE, FALSE, FALSE),
    plates = c("initial", "initial", "candidate", "candidate")
)
write_screen(
    holdout_path,
    observations = c(0.3, 0.7),
    mask = c(TRUE, TRUE),
    plates = c("candidate", "candidate")
)
converted <- merck_convert(training_path, holdout_path, output_path)
stopifnot(
    file.exists(output_path),
    identical(dim(converted$training), c(4L, 10L)),
    identical(converted$metadata$initial_plate, "initial"),
    all(is.na(converted$training$drug_2)),
    all(is.na(converted$training$dose_2)),
    isTRUE(all.equal(converted$training$truth[[1L]], stats::qlogis(0.01))),
    isTRUE(all.equal(converted$training$truth[[4L]], stats::qlogis(0.99))),
    all(is.na(converted$training$response[3:4])),
    all(is.finite(converted$training$truth)),
    all(converted$holdout$observed)
)

bad_path <- file.path(directory, "bad.h5")
write_screen(
    bad_path,
    observations = c(0.2, 0.3),
    mask = c(TRUE, FALSE),
    plates = c("mixed", "mixed")
)
bad <- tryCatch(
    merck_convert_screen(bad_path),
    error = function(error) conditionMessage(error)
)
stopifnot(grepl("entirely observed or entirely masked", bad, fixed = TRUE))
unlink(directory, recursive = TRUE)
message("Merck converter self-test passed")
