MERCK_DATA_URL <- paste0(
    "https://zenodo.org/records/12764821/files/",
    "merck_2016.screen.h5?download=1"
)
MERCK_DATA_MD5 <- "5d7cb090002a01ebcc0c06f0d5f72f67"
UPSTREAM_ARCHIVE_URL <- paste0(
    "https://zenodo.org/api/records/12765294/files/",
    "tansey-lab/batchie-v0.0.1-zenodo.zip/content"
)
UPSTREAM_ARCHIVE_MD5 <- "7deaf0ab4d9fb4c935e857d1937762bb"
UPSTREAM_COMMIT <- "6baa258cb5d77430dccab38b73b0908f91503c91"

merck_profiles <- function() {
    list(
        smoke = list(
            rank = 2L,
            chains = 1L,
            parallel_chains = 1L,
            iter_warmup = 2L,
            iter_sampling = 3L,
            thin = 1L,
            max_rounds = 1L,
            candidate_limit = 9L
        ),
        scaled = list(
            rank = 12L,
            chains = 2L,
            parallel_chains = 2L,
            iter_warmup = 10L,
            iter_sampling = 20L,
            thin = 2L,
            max_rounds = 3L,
            candidate_limit = Inf
        ),
        full = list(
            rank = 12L,
            chains = 2L,
            parallel_chains = 2L,
            iter_warmup = 2000L,
            iter_sampling = 8000L,
            thin = 40L,
            max_rounds = 78L,
            candidate_limit = Inf
        )
    )
}

MERCK_BATCH_SIZE <- 9L
MERCK_PREPARATION_SEED <- 0L
MERCK_SCORING_SEED <- 12L
MERCK_MAX_TRIPLETS <- 5000L
MERCK_HOLDOUT_FRACTION <- 0.1
