#' Run zonal statistics for a single region
#'
#' Loads cached zone artifacts from `zone_dir` for `region_id` and `blocksize`,
#' aligns the predictor (optionally projecting), and runs fast zonal summaries.
#'
#' @param region_id Character region identifier used in filenames.
#' @param zone_dir Directory containing zone artifacts produced by [sr_optimize_zones()].
#' @param predictor_path Path to predictor raster.
#' @param blocksize Integer block size used when building zone artifacts.
#' @param method Resampling method. Currently only `"near"` is supported.
#'   Zonal alignment uses nearest-cell center matching (no interpolation).
#'   If projection is required, this method is also used in `terra::project()`.
#' @param stats Character vector of statistics. Supported: `"sum"`, `"mean"`, `"count"` (alias: `"n"`).
#' @param progress_every Print progress every N windows (0 disables).
#' @param project_predictor_if_needed Logical; if TRUE, project predictor to match Zidx CRS.
#' @return A `data.table` with one row per GRIDCODE.
#' @export
sr_zonal_region <- function(
    predictor_path,
    zones,                 # sr_zones object
    method = "near",
    stats = c("sum", "mean", "n"),
    pad_cells = 1L,
    exact_crop = FALSE,
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
  if (!is.null(res$n))    out$n <- res$n
  if (!is.null(res$mean)) out$mean <- res$mean

  out
}
