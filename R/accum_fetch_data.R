#' Fetch optional accumulation datasets (on-demand)
#'
#' Downloads and caches optional datasets used for watershed accumulation presets
#' (e.g., MR NHDPlus CONUS). This avoids bundling large network files inside the
#' package. By default, data are stored in a user cache directory.
#'
#' This version downloads a ZIP directly from a GitHub repository (raw URL),
#' rather than using GitHub Releases. (You can switch to Releases later.)
#'
#' @param dataset Character name of the dataset to fetch. Currently supports
#'   `"mr_nhdplus_conus"`.
#' @param dest_dir Destination directory to store the cached dataset. If NULL,
#'   uses a package cache directory.
#' @param force Logical; if TRUE, re-download even if files already exist.
#' @param quiet Logical; if TRUE, suppress messages.
#' @param repo GitHub repository in `"owner/repo"` form that hosts the ZIP file
#'   in the repository (e.g., `"ryan-hill/StreamCatR-data"`).
#' @param branch Git branch or tag name to download from (default `"main"`).
#' @param asset Filename of the ZIP in the repo (default `"nhdplus_mr_conus_inputs.zip"`).
#'
#' @return A named list with `data_dir` and file paths for the dataset.
#' @export
sr_fetch_accum_data <- function(
    dataset = c("mr_nhdplus_conus"),
    dest_dir = NULL,
    force = FALSE,
    quiet = FALSE,
    repo = "ryan-hill/StreamCatR-data",
    branch = "main",
    asset = "nhdplus_mr_conus_inputs.zip"
) {
  dataset <- match.arg(dataset)

  if (is.null(dest_dir)) {
    dest_dir <- .sr_accum_cache_dir(dataset)
  }
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  # Expected files
  expected <- .sr_accum_expected_paths(dataset, dest_dir)

  if (!force && all(file.exists(unlist(expected$paths)))) {
    if (!quiet) message("Accumulation dataset already present at: ", dest_dir)
    return(expected)
  }

  # Build a raw GitHub URL for the ZIP
  # NOTE: this works when the file is committed to the repo at repo root (or provide path in `asset`)
  zip_url <- sprintf(
    "https://raw.githubusercontent.com/%s/%s/%s",
    repo, branch, asset
  )

  tmp <- tempfile(fileext = ".zip")
  if (!quiet) message("Downloading ", dataset, " from: ", zip_url)

  # Download
  ok <- try(
    utils::download.file(zip_url, destfile = tmp, mode = "wb", quiet = quiet),
    silent = TRUE
  )
  if (inherits(ok, "try-error") || !file.exists(tmp) || file.info(tmp)$size <= 0) {
    stop(
      "Download failed.\n",
      "Tried URL: ", zip_url, "\n",
      "Check that `repo`, `branch`, and `asset` are correct and that the file is accessible.",
      call. = FALSE
    )
  }

  if (!quiet) message("Unzipping to: ", dest_dir)
  utils::unzip(tmp, exdir = dest_dir)

  # Verify
  missing <- unlist(expected$paths)[!file.exists(unlist(expected$paths))]
  if (length(missing)) {
    stop(
      "Dataset was downloaded but expected files are missing:\n",
      paste(missing, collapse = "\n"),
      call. = FALSE
    )
  }

  expected
}
