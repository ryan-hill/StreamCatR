#' Run zonal statistics for a single region
#'
#' Runs fast zonal summaries for one optimized zone set (`sr_zones`) and one predictor raster.
#' This is the packaged analogue of the validated 05 workflow, using cached zone artifacts
#' produced by [sr_optimize_zones()] and the true-blocking zonal engine modeled after 04c.
#'
#' @param predictor_path Path to predictor raster.
#' @param zones An `sr_zones` object produced by [sr_optimize_zones()] or constructed with [sr_zones()].
#' @param method Resampling method. Currently only `"near"` is supported.
#'   Zonal alignment uses nearest-cell center matching (no interpolation).
#'   If projection is required, this method is also used in `terra::project()`.
#' @param stats Character vector of statistics. Supported: `"sum"`, `"mean"`, `"count"` (alias: `"n"`).
#' @param pad_cells Integer predictor-cell padding around each zone window. Default `1L`, matching 04c.
#' @param exact_crop Logical; if `TRUE`, use `terra::crop()` for predictor windows, matching 04c most closely.
#' @param progress_every Print progress every N windows (0 disables).
#' @return A `data.frame` with one row per GRIDCODE.
#' @export
sr_zonal_region <- function(
    predictor_path,
    zones,                 # sr_zones object
    method = "near",
    stats = c("sum", "mean", "count"),
    pad_cells = 1L,
    exact_crop = TRUE,
    progress_every = 0L
) {
  stopifnot(file.exists(predictor_path))
  stopifnot(file.exists(zones$zidx_path))
  stopifnot(file.exists(zones$grid_index_path))
  stopifnot(!is.null(zones$wins_path) && file.exists(zones$wins_path))

  P <- terra::rast(predictor_path)

  # Project predictor once if needed (05 behavior)
  if (!.sr_same_crs(P, terra::crs(terra::rast(zones$zidx_path)))) {
    tmp_pred <- tempfile(fileext = ".tif")
    on.exit(if (file.exists(tmp_pred)) suppressWarnings(file.remove(tmp_pred)), add = TRUE)
    P <- terra::project(P, terra::crs(terra::rast(zones$zidx_path)), method = method, filename = tmp_pred, overwrite = TRUE)
  }

  Zidx <- terra::rast(zones$zidx_path)

  # ids in idx order (05 behavior)
  ids_tbl <- arrow::read_parquet(zones$grid_index_path, as_data_frame = TRUE)
  ids_tbl <- ids_tbl[order(ids_tbl$idx), , drop = FALSE]
  ids <- ids_tbl$GRIDCODE

  # windows
  wins <- arrow::read_parquet(zones$wins_path, as_data_frame = TRUE)

  # run
  res <- .sr_zonal_engine(
    predictor = P,
    indexed_zones = Zidx,
    ids = ids,
    windows = wins,
    method = method,
    stats = stats,
    pad_cells = pad_cells,
    exact_crop = exact_crop,
    progress_every = progress_every
  )

  # attach GRIDCODE labeling (05 style output)
  out <- data.frame(GRIDCODE = ids, stringsAsFactors = FALSE)
  if (!is.null(res$sum))  out$sum <- res$sum
  if (!is.null(res$n)) out$count <- res$n
  if (!is.null(res$mean)) out$mean <- res$mean

  out
}
