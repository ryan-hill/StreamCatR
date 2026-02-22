#' Optimize catchments for fast zonal summaries
#'
#' Creates three cached artifacts used by fast zonal operations:
#' 1) GRIDCODE -> idx mapping parquet
#' 2) Zidx raster (idx 1..K)
#' 3) non-empty windows parquet (optional)
#'
#' @param catchment_path Path to a catchment raster (GRIDCODE values).
#' @param region_id Character region identifier used in output filenames.
#' @param blocksize Integer block size (pixels) for tiling and windows.
#' @param outdir Output directory for cached artifacts.
#' @param target_crs Target CRS (default "EPSG:5070").
#' @param overwrite_rasters If TRUE, overwrite existing cached outputs.
#' @param build_windows If TRUE, write non-empty windows parquet.
#' @param progress_every Progress frequency for streaming scan (0 disables).
#' @return An `sr_zones` object with paths to cached artifacts.
#' @export
sr_optimize_zones <- function(
    catchment_path,
    region_id,
    blocksize,
    outdir,
    target_crs        = "EPSG:5070",
    overwrite_rasters = TRUE,
    build_windows     = TRUE,
    progress_every    = 0L
) {
  stopifnot(is.character(catchment_path), length(catchment_path) == 1L)
  stopifnot(is.character(region_id), length(region_id) == 1L)
  stopifnot(is.character(outdir), length(outdir) == 1L)

  blocksize <- as.integer(blocksize)
  if (is.na(blocksize) || blocksize <= 0L) {
    stop("`blocksize` must be a positive integer.", call. = FALSE)
  }

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  grid_index_path <- file.path(outdir, sprintf("%s_gridcode_index.parquet", region_id))
  zidx_path       <- file.path(outdir, sprintf("%s_zone_index_block%d.tif", region_id, blocksize))
  wins_path       <- file.path(outdir, sprintf("%s_nonempty_windows_%d.parquet", region_id, blocksize))

  safe_unlink <- function(path, do = TRUE) {
    if (!isTRUE(do)) return(invisible(FALSE))
    sidecars <- c("", ".aux.xml", ".ovr", ".msk", ".msk.aux.xml")
    files <- file.path(dirname(path), paste0(basename(path), sidecars))
    suppressWarnings(file.remove(files[file.exists(files)]))
    invisible(TRUE)
  }

  src <- terra::rast(catchment_path)

  tmp_proj <- NULL
  if (!terra::same.crs(src, target_crs)) {
    tmp_proj <- tempfile(fileext = ".tif")
    src <- terra::project(src, target_crs, method = "near", filename = tmp_proj, overwrite = TRUE)
  }

  tmp_work <- tempfile(fileext = ".tif")
  terra::writeRaster(
    src, filename = tmp_work, overwrite = TRUE, datatype = "INT4S",
    gdal = c(
      "TILED=YES", "COMPRESS=ZSTD", "ZSTD_LEVEL=9",
      "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS",
      sprintf("BLOCKXSIZE=%d", blocksize),
      sprintf("BLOCKYSIZE=%d", blocksize)
    )
  )
  work <- terra::rast(tmp_work)

  # (Optional) cleanup of temp files on exit; safe because we only remove after function ends
  on.exit({
    if (!is.null(tmp_proj) && file.exists(tmp_proj)) suppressWarnings(file.remove(tmp_proj))
    if (!is.null(tmp_work) && file.exists(tmp_work)) suppressWarnings(file.remove(tmp_work))
  }, add = TRUE)

  work_trim <- try(terra::trim(work), silent = TRUE)
  if (inherits(work_trim, "try-error")) work_trim <- work

  if (terra::is.factor(work_trim)) {
    try({ levels(work_trim) <- NULL }, silent = TRUE)
  }

  if (file.exists(grid_index_path) && !overwrite_rasters) {
    message("Index exists, skipping: ", grid_index_path)
  } else {
    if (file.exists(grid_index_path)) safe_unlink(grid_index_path, TRUE)
    .sr_build_gridcode_index_parquet_stream(work_trim, grid_index_path, progress_every = progress_every)
  }

  ids_tbl <- arrow::read_parquet(grid_index_path, as_data_frame = TRUE)
  if (!"GRIDCODE" %in% names(ids_tbl)) {
    stop("Index parquet lacks GRIDCODE column: ", grid_index_path, call. = FALSE)
  }

  # Avoid integer64 surprises
  ids <- if (max(ids_tbl$GRIDCODE, na.rm = TRUE) <= .Machine$integer.max) {
    as.integer(ids_tbl$GRIDCODE)
  } else {
    as.numeric(ids_tbl$GRIDCODE)
  }

  Zidx <- terra::app(work_trim, fun = function(x) {
    out <- match(x, ids)
    out[is.na(x)] <- NA_integer_
    as.integer(out)
  })

  if (file.exists(zidx_path) && !overwrite_rasters) {
    message("Zidx exists, skipping: ", zidx_path)
  } else {
    if (file.exists(zidx_path)) safe_unlink(zidx_path, TRUE)
    terra::writeRaster(
      Zidx, zidx_path, overwrite = TRUE, datatype = "INT4S",
      gdal = c(
        "TILED=YES", "COMPRESS=ZSTD", "ZSTD_LEVEL=9",
        "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS",
        sprintf("BLOCKXSIZE=%d", blocksize),
        sprintf("BLOCKYSIZE=%d", blocksize)
      )
    )
  }

  wins_out <- NULL
  if (isTRUE(build_windows)) {
    if (file.exists(wins_path) && !overwrite_rasters) {
      message("Windows parquet exists, skipping: ", wins_path)
      wins_out <- wins_path
    } else {
      if (file.exists(wins_path)) safe_unlink(wins_path, TRUE)

      Zidx_disk <- terra::rast(zidx_path)
      nr <- terra::nrow(Zidx_disk)
      nc <- terra::ncol(Zidx_disk)

      occ <- terra::aggregate(
        !is.na(Zidx_disk),
        fact = c(blocksize, blocksize),
        fun  = function(x) as.integer(any(as.logical(x)))
      )

      vals <- terra::values(occ, mat = FALSE)
      if (any(vals == 1L, na.rm = TRUE)) {
        cells <- which(vals == 1L)
        rc <- terra::rowColFromCell(occ, cells)
        r0 <- (rc[, 1L] - 1L) * blocksize + 1L
        c0 <- (rc[, 2L] - 1L) * blocksize + 1L
        nrw <- pmin(blocksize, nr - r0 + 1L)
        ncw <- pmin(blocksize, nc - c0 + 1L)
        wins_dt <- data.frame(row = r0, col = c0, nrows = nrw, ncols = ncw)
      } else {
        wins_dt <- data.frame(row = integer(), col = integer(), nrows = integer(), ncols = integer())
      }

      arrow::write_parquet(wins_dt, wins_path, compression = "zstd")
      wins_out <- wins_path
    }
  }

  invisible(sr_zones(
    region_id       = region_id,
    blocksize       = blocksize,
    zone_dir        = outdir,
    grid_index_path = grid_index_path,
    zidx_path       = zidx_path,
    wins_path       = if (isTRUE(build_windows)) wins_path else NULL
  ))
}



