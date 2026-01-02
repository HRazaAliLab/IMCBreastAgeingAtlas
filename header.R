suppressPackageStartupMessages({
  library(here)
})

# ---- Required packages (core) ----
pkgs <- c(
  "data.table", "parallel",
  "ggplot2", "patchwork",
  "arrow", "fst", "latex2exp", "effsize", "fs",
  "grid"
)

checkPackages <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop(
      "Missing required R packages:\n  ",
      paste(missing, collapse = ", "),
      "\n\nInstall them (CRAN/Bioc) and re-run.\n",
      "Tip: run `install.packages(<pkg>)` for CRAN packages.\n",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

checkPackages(pkgs)

# Attach (quietly)
invisible(lapply(pkgs, function(p) {
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}))

# ---- Platform note ----
if (.Platform$OS.type == "windows") {
  warning(
    "This repository uses fork-based parallelism (e.g., mclapply) in some scripts, ",
    "which is not supported on Windows. Use Linux/macOS or replace mclapply with ",
    "a Windows-compatible backend where needed.",
    call. = FALSE
  )
}

# ---- data.table threads ----
data.table::setDTthreads(min(parallel::detectCores(), 64))

# ---- Load project utilities ----
source(here::here("code/generalUtilities.R"))
