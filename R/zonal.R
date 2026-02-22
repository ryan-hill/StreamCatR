# small internal helper (avoids importing rlang)
#' @keywords internal
#' @noRd
`%||%` <- function(x, y) if (is.null(x)) y else x

#' Run fast zonal statistics for continuous predictors
#'
#' Computes summary statistics of a (continuous) predictor raster within catchment zones
#' prepared by [sr_optimize_zones()].
#'
#' @param predictor A `terra::SpatRaster` or a path to a raster readable by terra.
#' @param zones Either:
#'   (1) an `sr_zones` object returned by [sr_optimize_zones()], or
#'   (2) a list returned by [sr_optimize_zones()] (expects `zidx_path`, `grid_index_path`, and `wins_path`), or
#'   (3) a list with `indexed_zones_path` (alias of `zidx_path`), `grid_index_path`, and `wins_path`, and optionally `ids`.
#' @param stats Character vector of statistics to compute. Currently supported: `"sum"`, `"mean"`, `"n"`.
#' @param method Resampling method when aligning predictor to zones (default `"near"`).
#' @param progress_every Print progress every N windows (0 disables).
#' @return A `data.table` with one row per zone and columns for requested statistics.
#' @export
sr_zonal <- function(
    predictor,
    zones,
    method = "near",
    stats = c("sum", "mean"),
    progress_every = 0L
) {
  if (is.character(predictor) && length(predictor) == 1L) {
    predictor <- terra::rast(predictor)
  }
  if (!inherits(predictor, "SpatRaster")) {
    stop("`predictor` must be a terra::SpatRaster or a single file path.", call. = FALSE)
  }

  if (!inherits(zones, "sr_zones")) {
    stop("`zones` must be an sr_zones object returned by sr_optimize_zones().", call. = FALSE)
  }

  # Ensure windows exist; build + cache if missing
  zones <- .sr_ensure_windows(zones)

  indexed_zones <- terra::rast(zones$zidx_path)
  wins_path <- zones$wins_path
  grid_index_path <- zones$grid_index_path

  # ids: read parquet, order by idx, coerce away from integer64
  ids_tbl <- arrow::read_parquet(grid_index_path, as_data_frame = TRUE)
  if (!"GRIDCODE" %in% names(ids_tbl)) {
    stop("Index parquet lacks GRIDCODE column: ", grid_index_path, call. = FALSE)
  }
  if ("idx" %in% names(ids_tbl)) ids_tbl <- ids_tbl[order(ids_tbl$idx), , drop = FALSE]

  ids_raw <- ids_tbl$GRIDCODE
  ids <- if (max(ids_raw, na.rm = TRUE) <= .Machine$integer.max) as.integer(ids_raw) else as.numeric(ids_raw)

  .sr_zonal_engine(
    predictor      = predictor,
    indexed_zones  = indexed_zones,
    ids            = ids,
    windows        = wins_path,
    method         = method,
    stats          = stats,
    progress_every = progress_every
  )
}

#' @keywords internal
#' @noRd
.sr_zonal_engine <- function(
    predictor,
    indexed_zones,
    ids,
    windows,
    method = "near",
    stats = c("sum", "mean"),
    progress_every = 0L
) {
  stats <- unique(stats)
  ok <- c("sum", "mean", "n")
  bad <- setdiff(stats, ok)
  if (length(bad)) {
    stop("Unsupported `stats`: ", paste(bad, collapse = ", "),
         ". Supported: ", paste(ok, collapse = ", "), call. = FALSE)
  }

  wins_dt <- windows
  if (is.character(wins_dt) && length(wins_dt) == 1L) {
    if (!file.exists(wins_dt)) stop("Missing windows parquet: ", wins_dt, call. = FALSE)
    wins_dt <- arrow::read_parquet(wins_dt, as_data_frame = TRUE)
  }
  if (!is.data.frame(wins_dt)) {
    stop("`windows` must be a parquet path or a data.frame with columns row/col/nrows/ncols.", call. = FALSE)
  }
  req <- c("row", "col", "nrows", "ncols")
  if (!all(req %in% names(wins_dt))) {
    stop("`windows` must contain columns: ", paste(req, collapse = ", "), call. = FALSE)
  }

  if (!terra::same.crs(predictor, indexed_zones)) {
    stop("CRS mismatch between `predictor` and `indexed_zones`. Project predictor first.", call. = FALSE)
  }

  exZ <- terra::ext(indexed_zones)
  rZ  <- terra::res(indexed_zones)
  rxZ <- rZ[1]; ryZ <- rZ[2]

  exP <- terra::ext(predictor)
  rP  <- terra::res(predictor)
  rxP <- rP[1]; ryP <- rP[2]
  ncPred <- terra::ncol(predictor)
  nrPred <- terra::nrow(predictor)

  K <- length(ids)
  sumv <- numeric(K)
  n    <- integer(K)

  terra::readStart(indexed_zones)
  terra::readStart(predictor)
  on.exit({
    terra::readStop(indexed_zones)
    terra::readStop(predictor)
  }, add = TRUE)

  nwin <- nrow(wins_dt)

  for (i in seq_len(nwin)) {
    r0 <- as.integer(wins_dt$row[i])
    c0 <- as.integer(wins_dt$col[i])
    nr <- as.integer(wins_dt$nrows[i])
    nc <- as.integer(wins_dt$ncols[i])

    if (nr <= 0L || nc <= 0L) next

    idx <- terra::readValues(indexed_zones, row = r0, nrows = nr, col = c0, ncols = nc, mat = FALSE)
    if (!length(idx) || all(is.na(idx))) next

    xminZ <- exZ$xmin + (c0 - 1L) * rxZ
    xmaxZ <- xminZ + nc * rxZ
    ymaxZ <- exZ$ymax - (r0 - 1L) * ryZ
    yminZ <- ymaxZ - nr * ryZ

    col_min <- floor((xminZ - exP$xmin) / rxP) + 1L
    col_max <- ceiling((xmaxZ - exP$xmin) / rxP)
    row_min <- floor((exP$ymax - ymaxZ) / ryP) + 1L
    row_max <- ceiling((exP$ymax - yminZ) / ryP)

    col_min <- max(1L, min(ncPred, col_min))
    col_max <- max(1L, min(ncPred, col_max))
    row_min <- max(1L, min(nrPred, row_min))
    row_max <- max(1L, min(nrPred, row_max))

    if (col_max < col_min || row_max < row_min) next

    ncP <- as.integer(col_max - col_min + 1L)
    nrP <- as.integer(row_max - row_min + 1L)

    valsP <- terra::readValues(predictor, row = row_min, nrows = nrP, col = col_min, ncols = ncP, mat = FALSE)
    if (!length(valsP)) next

    xminP_sub <- exP$xmin + (col_min - 1L) * rxP
    ymaxP_sub <- exP$ymax - (row_min - 1L) * ryP

    acc_sum_n_map_centers1K(
      idx      = as.integer(idx),
      xminZ    = xminZ,
      ymaxZ    = ymaxZ,
      rxZ      = rxZ,
      ryZ      = ryZ,
      nx       = as.integer(nc),
      ny       = as.integer(nr),
      valsP    = as.numeric(valsP),
      xminP    = xminP_sub,
      ymaxP    = ymaxP_sub,
      rxP      = rxP,
      ryP      = ryP,
      ncP      = as.integer(ncP),
      nrP      = as.integer(nrP),
      sumv     = sumv,
      n        = n
    )

    if (progress_every > 0L && (i %% progress_every == 0L)) {
      message(sprintf("processed window %d / %d", i, nwin))
    }
  }

  out <- data.table::data.table(GRIDCODE = ids)

  if ("sum" %in% stats) out[["sum"]] <- sumv
  if ("n" %in% stats)   out[["n"]]   <- n

  if ("mean" %in% stats) {
    m <- rep(NA_real_, length(sumv))
    okn <- n > 0L
    m[okn] <- sumv[okn] / n[okn]
    out[["mean"]] <- m
  }

  out
}
