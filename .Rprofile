# callr children already inherit the parent's library paths. Re-activating
# renv can exceed the language server's background-session startup timeout.
if (is.na(Sys.getenv("CALLR_CHILD_R_LIBS", unset = NA))) {
  source("renv/activate.R")
}
