merck_script_dir <- function() {
    command <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", command, value = TRUE)
    if (length(file_arg)) {
        return(dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))))
    }
    normalizePath(getwd())
}

merck_parse_args <- function(args, flags = character()) {
    result <- list()
    i <- 1L
    while (i <= length(args)) {
        key <- args[[i]]
        if (!startsWith(key, "--")) {
            stop("Unexpected argument: ", key, call. = FALSE)
        }
        name <- gsub("-", "_", substring(key, 3L), fixed = TRUE)
        if (name %in% flags) {
            result[[name]] <- TRUE
            i <- i + 1L
        } else {
            if (i == length(args)) stop("Missing value for ", key, call. = FALSE)
            result[[name]] <- args[[i + 1L]]
            i <- i + 2L
        }
    }
    result
}

merck_require <- function(packages) {
    missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
        stop(
            "Missing experiment-only R package(s): ",
            paste(missing, collapse = ", "),
            call. = FALSE
        )
    }
}

merck_atomic_save_rds <- function(object, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- paste0(path, ".tmp-", Sys.getpid())
    on.exit(unlink(temporary), add = TRUE)
    saveRDS(object, temporary)
    if (!file.rename(temporary, path)) {
        stop("Could not atomically commit ", path, call. = FALSE)
    }
    invisible(path)
}

merck_atomic_write_json <- function(object, path, pretty = TRUE) {
    merck_require("jsonlite")
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- paste0(path, ".tmp-", Sys.getpid())
    on.exit(unlink(temporary), add = TRUE)
    jsonlite::write_json(
        object,
        temporary,
        auto_unbox = TRUE,
        pretty = pretty,
        null = "null",
        na = "null"
    )
    if (!file.rename(temporary, path)) {
        stop("Could not atomically commit ", path, call. = FALSE)
    }
    invisible(path)
}

merck_atomic_write_csv <- function(object, path) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    temporary <- paste0(path, ".tmp-", Sys.getpid())
    on.exit(unlink(temporary), add = TRUE)
    utils::write.csv(object, temporary, row.names = FALSE)
    if (!file.rename(temporary, path)) {
        stop("Could not atomically commit ", path, call. = FALSE)
    }
    invisible(path)
}

merck_git_revision <- function(root) {
    result <- suppressWarnings(system2(
        "git", c("-C", shQuote(root), "rev-parse", "HEAD"),
        stdout = TRUE, stderr = FALSE
    ))
    if (!is.null(attr(result, "status"))) return(NA_character_)
    if (length(result)) result[[1L]] else NA_character_
}

merck_subset_profile <- function(dataset, profile) {
    candidate_plates <- sort(setdiff(
        unique(dataset$training$plate),
        dataset$metadata$initial_plate
    ))
    if (is.finite(profile$candidate_limit)) {
        candidate_plates <- head(candidate_plates, profile$candidate_limit)
    }
    keep <- c(dataset$metadata$initial_plate, candidate_plates)
    dataset$training <- dataset$training[dataset$training$plate %in% keep, , drop = FALSE]
    dataset$holdout <- dataset$holdout[dataset$holdout$plate %in% candidate_plates, , drop = FALSE]
    dataset$metadata$profile_candidate_plates <- candidate_plates
    dataset
}

merck_heap_mb <- function(collection = gc()) {
    sum(collection[, "(Mb)"])
}
