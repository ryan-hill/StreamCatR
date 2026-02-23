#' Fetch optional accumulation datasets (on-demand)
#'
#' Downloads and caches optional datasets used for watershed accumulation presets
#' (e.g., MR NHDPlus CONUS). This avoids bundling large network files inside the
#' package. By default, data are stored in a user cache directory.
#'
#' @param dataset Character name of the dataset to fetch. Currently supports
#'   `"mr_nhdplus_conus"`.
#' @param dest_dir Destination directory to store the cached dataset. If NULL,
#'   uses a package cache directory.
#' @param force Logical; if TRUE, re-download even if files already exist.
#' @param quiet Logical; if TRUE, suppress messages.
#' @param repo GitHub repository in `"owner/repo"` form that hosts the Release asset.
#'   You can keep this the same as the package repo or use a separate data repo.
#' @param tag GitHub Release tag to download from (e.g., `"accum-data-v1"`).
#' @param asset Filename of the Release asset (ZIP) to download.
#'
#' @return A named list with `data_dir` and file paths for the dataset.
#' @export
sr_fetch_accum_data <- function(
    dataset = c("mr_nhdplus_conus"),
    dest_dir = NULL,
    force = FALSE,
    quiet = FALSE,
    repo = "YOUR_ORG/StreamCatR",
    tag = "accum-data-v1",
    asset = "mr_nhdplus_conus_accum_inputs.zip"
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

  # Download from GitHub Releases
  if (!requireNamespace("piggyback", quietly = TRUE)) {
    stop(
      "Package 'piggyback' is required to fetch data from GitHub Releases.\n",
      "Install it with: install.packages('piggyback')",
      call. = FALSE
    )
  }

  tmp <- tempfile(fileext = ".zip")
  if (!quiet) message("Downloading ", dataset, " from GitHub Releases...")

  piggyback::pb_download(
    file = asset,
    repo = repo,
    tag = tag,
    dest = tmp
  )

  if (!file.exists(tmp)) stop("Download failed for asset: ", asset, call. = FALSE)

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
