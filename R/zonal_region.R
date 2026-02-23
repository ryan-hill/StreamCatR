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
    region_id,
    zone_dir,
    predictor_path,
    blocksize,
    method = "near",
    stats = c("sum", "mean", "count"),
    progress_every = 0L,
    project_predictor_if_needed = FALSE
) {
  blocksize <- as.integer(blocksize)

  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    stop("`method` must be a single, non-NA character string.", call. = FALSE)
  }
  if (!identical(method, "near")) {
    stop("Only method = 'near' is currently supported.", call. = FALSE)
  }

  zidx_path <- file.path(zone_dir, sprintf("%s_zone_index_block%d.tif", region_id, blocksize))
  ids_path  <- file.path(zone_dir, sprintf("%s_gridcode_index.parquet", region_id))
  wins_path <- file.path(zone_dir, sprintf("%s_nonempty_windows_%d.parquet", region_id, blocksize))

  if (!file.exists(zidx_path)) stop("Missing Zidx raster: ", zidx_path, call. = FALSE)
  if (!file.exists(ids_path))  stop("Missing ID index parquet: ", ids_path, call. = FALSE)
  if (!file.exists(wins_path)) stop("Missing nonempty windows parquet: ", wins_path, call. = FALSE)
  if (!file.exists(predictor_path)) stop("Missing predictor raster: ", predictor_path, call. = FALSE)

  Zidx <- terra::rast(zidx_path)
  R <- terra::rast(predictor_path)

  if (!terra::same.crs(R, Zidx)) {
    if (!project_predictor_if_needed) {
      stop(
        "CRS mismatch between predictor and Zidx. ",
        "Supply a predictor in the same CRS as Zidx (recommended), ",
        "or set project_predictor_if_needed = TRUE.",
        call. = FALSE
      )
    }
    message("Projecting predictor to match Zidx CRS (one-time)...")
    tmp_pred <- tempfile(fileext = ".tif")
    R <- terra::project(R, Zidx, method = method, filename = tmp_pred, overwrite = TRUE)
  }

  zones <- sr_zones(
    region_id       = region_id,
    blocksize       = blocksize,
    zone_dir        = zone_dir,
    grid_index_path = ids_path,
    zidx_path       = zidx_path,
    wins_path       = wins_path
  )

  out <- sr_zonal(
    predictor      = R,
    zones          = zones,
    stats          = stats,
    progress_every = progress_every
  )

  attr(out, "region_id")      <- region_id
  attr(out, "predictor_path") <- predictor_path
  attr(out, "blocksize")      <- blocksize
  attr(out, "paths") <- list(zidx = zidx_path, ids = ids_path, wins = wins_path, predictor = predictor_path)

  out
}
