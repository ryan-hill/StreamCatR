#' Run zonal statistics for a single region
#'
#' Loads cached zone artifacts from `zone_dir` for `region_id` and `blocksize`,
#' aligns the predictor (optionally projecting), and runs fast zonal summaries.
#'
#' @param region_id Character region identifier used in filenames.
#' @param zone_dir Directory containing zone artifacts produced by [sr_optimize_zones()].
#' @param predictor_path Path to predictor raster.
#' @param blocksize Integer block size used when building zone artifacts.
#' @param method Resampling method (default "near"). Used only if projecting predictor.
#' @param stats Character vector of statistics. Supported: `"sum"`, `"mean"`, `"n"`.
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
    stats = c("sum", "n"),
    progress_every = 0L,
    project_predictor_if_needed = FALSE
) {
  blocksize <- as.integer(blocksize)

  zidx_path <- file.path(zone_dir, sprintf("%s_zone_index_block%d.tif", region_id, blocksize))
  ids_path  <- file.path(zone_dir, sprintf("%s_gridcode_index.parquet", region_id))
  wins_path <- file.path(zone_dir, sprintf("%s_nonempty_windows_%d.parquet", region_id, blocksize))

  if (!file.exists(zidx_path)) stop("Missing Zidx raster: ", zidx_path)
  if (!file.exists(ids_path))  stop("Missing ID index parquet: ", ids_path)
  if (!file.exists(wins_path)) stop("Missing nonempty windows parquet: ", wins_path)
  if (!file.exists(predictor_path)) stop("Missing predictor raster: ", predictor_path)

  Zidx <- terra::rast(zidx_path)

  ids_tbl <- arrow::read_parquet(ids_path, as_data_frame = TRUE)
  if (!"GRIDCODE" %in% names(ids_tbl)) stop("Index parquet lacks GRIDCODE column: ", ids_path)
  if ("idx" %in% names(ids_tbl)) ids_tbl <- ids_tbl[order(ids_tbl$idx), , drop = FALSE]

  ids_raw <- ids_tbl$GRIDCODE
  ids <- if (max(ids_raw, na.rm = TRUE) <= .Machine$integer.max) {
    as.integer(ids_raw)
  } else {
    as.numeric(ids_raw)
  }

  wins_dt <- arrow::read_parquet(wins_path, as_data_frame = TRUE)
  R <- terra::rast(predictor_path)

  if (!terra::same.crs(R, Zidx)) {
    if (!project_predictor_if_needed) {
      stop(
        "CRS mismatch between predictor and Zidx. ",
        "Supply a predictor in the same CRS as Zidx (recommended), ",
        "or set project_predictor_if_needed = TRUE."
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
    wins_path       = if (file.exists(wins_path)) wins_path else NULL
  )

  out <- sr_zonal(
    predictor = R,
    zones = zones,
    method = method,
    stats = stats,
    progress_every = progress_every
  )

  attr(out, "region_id")      <- region_id
  attr(out, "predictor_path") <- predictor_path
  attr(out, "blocksize")      <- blocksize
  attr(out, "paths") <- list(zidx = zidx_path, ids = ids_path, wins = wins_path, predictor = predictor_path)

  out
}
