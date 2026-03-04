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
#' @param stats Character vector of statistics to compute. Currently supported: `"sum"`, `"mean"`, `"count"` (alias: `"n"`).
#' @param method Resampling method. Currently only `"near"` is supported.
#'   Zonal alignment uses nearest-cell center matching (no interpolation).
#'   If projection is required, this method is also used in `terra::project()`.
#' @param progress_every Print progress every N windows (0 disables).
#' @return A `data.table` with one row per zone and columns for requested statistics.
#' @export
sr_zonal <- function(
    predictor,
    zones,
    method = "near",
    stats = c("sum", "mean", "count"),
    pad_cells = 1L,
    exact_crop = TRUE,
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
}

#' @keywords internal
#' @noRd
.sr_zonal_engine <- function(
    predictor,
    indexed_zones,
    ids,
    windows,
    method = c("near"),
    stats = c("sum", "mean", "n"),
    pad_cells = 1L,
    exact_crop = TRUE,
    progress_every = 0L
) {

  # allow windows to be either a parquet path (character) or an in-memory data.frame
  if (is.character(windows) && length(windows) == 1L) {
    windows <- arrow::read_parquet(windows, as_data_frame = TRUE)
  }

  # keep using the variable name 'Zidx' in the existing body
  Zidx <- indexed_zones

  method <- match.arg(method)
  stats <- match.arg(stats, several.ok = TRUE)
  pad_cells <- as.integer(pad_cells)
  if (is.na(pad_cells) || pad_cells < 0L) pad_cells <- 0L

  acc <- .sr_init_accumulator(prefer_map = TRUE)
  acc_idx <- acc$acc_idx
  acc_map <- acc$acc_map
  use_cpp_map <- is.function(acc_map)

  # ids: ensure numeric/integer vector used by match
  # (04c used integer when possible, numeric otherwise)
  if (length(ids) == 0L) stop("ids is empty", call. = FALSE)
  ids <- if (max(ids, na.rm = TRUE) <= .Machine$integer.max) as.integer(ids) else as.numeric(ids)

  # Allocate accumulators: idx is 1..K (K = length(ids))
  K <- length(ids)
  sumv <- numeric(K)
  n    <- integer(K)

  exZ <- terra::ext(Zidx)
  resZ <- terra::res(Zidx)
  rxZ <- resZ[1]
  ryZ <- resZ[2]

  exP <- terra::ext(predictor)
  resP <- terra::res(predictor)
  rxP <- resP[1]
  ryP <- resP[2]

  # ReadStart for speed
  terra::readStart(Zidx)
  terra::readStart(predictor)
  on.exit({
    terra::readStop(Zidx)
    terra::readStop(predictor)
  }, add = TRUE)

  # Iterate windows (must include row/col/nrows/ncols or r0/c0/nr/nc depending on your windows parquet)
  # Here I assume windows has columns: row, col, nrows, ncols (1-based), matching 04c/05 scripts.
  stopifnot(all(c("row", "col", "nrows", "ncols") %in% names(windows)))

  nwin <- nrow(windows)
  if (nwin == 0L) {
    return(.sr_finish_stats(sumv, n, stats))
  }

  for (w in seq_len(nwin)) {
    r0 <- as.integer(windows$row[w])
    c0 <- as.integer(windows$col[w])
    nr <- as.integer(windows$nrows[w])
    nc <- as.integer(windows$ncols[w])
    if (nr <= 0L || nc <= 0L) next

    # ---- read idx window from Zidx ----
    idx <- terra::readValues(Zidx, row = r0, col = c0, nrows = nr, ncols = nc, mat = FALSE)
    if (!any(is.finite(idx))) next

    # ---- compute zone window extent (04c style) ----
    xminZ <- exZ$xmin + (c0 - 1L) * rxZ
    xmaxZ <- xminZ + nc * rxZ
    ymaxZ <- exZ$ymax - (r0 - 1L) * ryZ
    yminZ <- ymaxZ - nr * ryZ
    ewin  <- terra::ext(xminZ, xmaxZ, yminZ, ymaxZ)

    # ---- padding (critical 04c fix) ----
    padx <- rxP * pad_cells
    pady <- ryP * pad_cells
    ewin_pad <- terra::ext(ewin$xmin - padx, ewin$xmax + padx, ewin$ymin - pady, ewin$ymax + pady)

    # skip if no overlap (04c)
    if (!terra::relate(ewin_pad, exP, "intersects")) next
    ewin2 <- terra::intersect(ewin_pad, exP)
    if (is.null(ewin2)) next

    # ---- read predictor subwindow corresponding to ewin2 ----
    if (isTRUE(exact_crop)) {
      # Closest to 04c: crop semantics
      Pwin <- terra::crop(predictor, ewin2)
      if (is.null(Pwin)) next
      valsP <- terra::values(Pwin, mat = FALSE)
      if (!length(valsP)) next
      exPsub <- terra::ext(Pwin)
      resPsub <- terra::res(Pwin)
      xminP_sub <- exPsub$xmin
      ymaxP_sub <- exPsub$ymax
      rxP_sub <- resPsub[1]
      ryP_sub <- resPsub[2]
      ncP <- terra::ncol(Pwin)
      nrP <- terra::nrow(Pwin)
    } else {
      # Faster: compute row/col bounds from ewin2 then readValues directly
      col_min <- floor((ewin2$xmin - exP$xmin) / rxP) + 1L
      col_max <- ceiling((ewin2$xmax - exP$xmin) / rxP)
      row_min <- floor((exP$ymax - ewin2$ymax) / ryP) + 1L
      row_max <- ceiling((exP$ymax - ewin2$ymin) / ryP)

      # clamp
      col_min <- max(1L, col_min); row_min <- max(1L, row_min)
      col_max <- min(terra::ncol(predictor), col_max)
      row_max <- min(terra::nrow(predictor), row_max)

      ncP <- col_max - col_min + 1L
      nrP <- row_max - row_min + 1L
      if (ncP <= 0L || nrP <= 0L) next

      valsP <- terra::readValues(predictor, row = row_min, col = col_min, nrows = nrP, ncols = ncP, mat = FALSE)
      if (!length(valsP)) next

      xminP_sub <- exP$xmin + (col_min - 1L) * rxP
      ymaxP_sub <- exP$ymax - (row_min - 1L) * ryP
      rxP_sub <- rxP
      ryP_sub <- ryP
    }

    # ---- accumulate: prefer C++ map-centers, else R fallback + acc_idx ----
    if (use_cpp_map) {
      acc_map(
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
        rxP      = rxP_sub,
        ryP      = ryP_sub,
        ncP      = as.integer(ncP),
        nrP      = as.integer(nrP),
        sumv     = sumv,
        n        = n
      )
    } else {
      # 04c fallback: map zone cell centers to predictor indices, then acc_idx
      # Build x/y centers for zone window
      xs <- seq(xminZ + rxZ/2, xmaxZ - rxZ/2, by = rxZ)
      ys <- seq(ymaxZ - ryZ/2, yminZ + ryZ/2, by = -ryZ)

      # Map centers into predictor window (need a raster object for colFromX/rowFromY)
      Pmap <- if (isTRUE(exact_crop)) Pwin else terra::crop(predictor, ewin2)
      if (is.null(Pmap)) next

      col_vec <- terra::colFromX(Pmap, xs)
      row_vec <- terra::rowFromY(Pmap, ys)

      # Ensure valsP/ncP consistent with Pmap
      if (!isTRUE(exact_crop)) {
        valsP <- terra::values(Pmap, mat = FALSE)
        ncP <- terra::ncol(Pmap)
        nrP <- terra::nrow(Pmap)
      }

      v <- rep(NA_real_, length(idx))
      for (iy in seq_len(nr)) {
        rr <- row_vec[iy]
        if (is.na(rr)) next
        base <- (rr - 1L) * ncP
        out_idx <- (iy - 1L) * nc + seq_len(nc)
        for (ix in seq_len(nc)) {
          cc <- col_vec[ix]
          if (is.na(cc)) next
          v[out_idx[ix]] <- valsP[base + cc]
        }
      }

      acc_idx(as.integer(idx), as.numeric(v), sumv, n)
    }

    if (progress_every > 0L && (w %% progress_every == 0L)) {
      message(sprintf("zonal: processed window %d / %d", w, nwin))
    }
  }

  .sr_finish_stats(sumv, n, stats)
}

#' @keywords internal
#' @noRd
.sr_finish_stats <- function(sumv, n, stats) {
  out <- list()
  if ("sum" %in% stats)  out$sum <- sumv
  if ("n" %in% stats || "count" %in% stats) out$n <- n
  if ("mean" %in% stats) {
    m <- rep(NA_real_, length(sumv))
    ok <- n > 0L
    m[ok] <- sumv[ok] / n[ok]
    out$mean <- m
  }
  out
}
