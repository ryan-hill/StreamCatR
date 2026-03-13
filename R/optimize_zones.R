#' Optimize zones raster for fast zonal operations
#'
#' Replicates the workflow of 01.build_optimized_catchments.R, but packaged and generic.
#'
#' @keywords internal
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

  # Ensure threading/memory settings are applied for this call
  try({
    terra::terraOptions(
      threads = max(1L, parallel::detectCores(logical = TRUE) - 1L),
      memfrac = 0.9
    )
    Sys.setenv(GDAL_NUM_THREADS = "ALL_CPUS")
    Sys.setenv(OMP_NUM_THREADS  = as.character(parallel::detectCores(logical = TRUE)))
  }, silent = TRUE)

  blocksize <- as.integer(blocksize)
  if (is.na(blocksize) || blocksize <= 0L) {
    stop("`blocksize` must be a positive integer.", call. = FALSE)
  }

  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

  # Timing ------------------------------------------------------------------
  t0 <- Sys.time()
  tick <- function(label) {
    now <- Sys.time()
    message(sprintf("[%.1fs] %s", as.numeric(difftime(now, t0, units = "secs")), label))
    t0 <<- now
  }
  tick("start")
  # Timing ------------------------------------------------------------------

  grid_index_path <- file.path(outdir, sprintf("%s_rasterid_index.parquet", region_id))
  zidx_path       <- file.path(outdir, sprintf("%s_zone_index.tif", region_id, blocksize))
  wins_path       <- file.path(outdir, sprintf("%s_nonempty_windows.parquet", region_id, blocksize))

  # ---- 1) Read source; project only if needed; ALWAYS write a tiled temp work raster ----
  src <- terra::rast(catchment_path)

  tmp_proj <- NULL
  if (!.sr_same_crs(src, target_crs)) {
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

  on.exit({
    if (!is.null(tmp_proj) && file.exists(tmp_proj)) suppressWarnings(file.remove(tmp_proj))
    if (!is.null(tmp_work) && file.exists(tmp_work)) suppressWarnings(file.remove(tmp_work))
  }, add = TRUE)

  tick("wrote temp work raster")

  # ---- 2) Trim NA edges; drop factor levels ----
  work_trim <- try(terra::trim(work), silent = TRUE)
  if (inherits(work_trim, "try-error")) work_trim <- work

  if (terra::is.factor(work_trim)) {
    try({ levels(work_trim) <- NULL }, silent = TRUE)
  }
  tick("trim finished")

  # ---- 2b) Short-circuit if raster has no cells (after trim) ----
  if (terra::nrow(work_trim) == 0L || terra::ncol(work_trim) == 0L || terra::ncell(work_trim) == 0L) {
    # Write empty index and empty outputs and return
    arrow::write_parquet(
      data.frame(GRIDCODE = integer(0), idx = integer(0)),
      grid_index_path, compression = "zstd"
    )

    if (file.exists(zidx_path)) .sr_safe_unlink(zidx_path, TRUE)
    terra::writeRaster(
      work_trim * NA, zidx_path, overwrite = TRUE, datatype = "INT4S",
      gdal = c(
        "TILED=YES", "COMPRESS=ZSTD", "ZSTD_LEVEL=9",
        "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS",
        sprintf("BLOCKXSIZE=%d", blocksize),
        sprintf("BLOCKYSIZE=%d", blocksize)
      )
    )

    if (isTRUE(build_windows)) {
      if (file.exists(wins_path)) .sr_safe_unlink(wins_path, TRUE)
      arrow::write_parquet(
        data.frame(row = integer(), col = integer(), nrows = integer(), ncols = integer()),
        wins_path, compression = "zstd"
      )
    }

    return(invisible(sr_zones(
      region_id       = region_id,
      blocksize       = blocksize,
      zone_dir        = outdir,
      grid_index_path = grid_index_path,
      zidx_path       = zidx_path,
      wins_path       = if (isTRUE(build_windows)) wins_path else NULL
    )))
  }

  # ---- 3) Build GRIDCODE -> idx mapping (streamed, like 01...R) ----
  built_now <- !(file.exists(grid_index_path) && !isTRUE(overwrite_rasters))
  if (isTRUE(built_now)) {
    if (file.exists(grid_index_path)) .sr_safe_unlink(grid_index_path, TRUE)
    ids_tbl <- .sr_build_gridcode_index_parquet_stream(
      Z = work_trim,
      out_file = grid_index_path,
      progress_every = progress_every
    )
  } else {
    message("Index exists, skipping: ", grid_index_path)
    ids_tbl <- arrow::read_parquet(grid_index_path, as_data_frame = TRUE)
  }
  if (!all(c("GRIDCODE", "idx") %in% names(ids_tbl))) {
    stop("Index parquet lacks GRIDCODE/idx columns: ", grid_index_path, call. = FALSE)
  }

  tick(sprintf("built/loaded index (%d codes)", nrow(ids_tbl)))

  # Handle empty mapping (all-NA raster)
  # Handle empty mapping (all-NA raster)
  if (nrow(ids_tbl) == 0L) {
    if (!file.exists(zidx_path) || isTRUE(overwrite_rasters)) {
      if (file.exists(zidx_path)) .sr_safe_unlink(zidx_path, TRUE)
      terra::writeRaster(
        work_trim * NA, zidx_path, overwrite = TRUE, datatype = "INT4S",
        gdal = c(
          "TILED=YES", "COMPRESS=ZSTD", "ZSTD_LEVEL=9",
          "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS",
          sprintf("BLOCKXSIZE=%d", blocksize),
          sprintf("BLOCKYSIZE=%d", blocksize)
        )
      )
    }

    if (isTRUE(build_windows) && (!file.exists(wins_path) || isTRUE(overwrite_rasters))) {
      if (file.exists(wins_path)) .sr_safe_unlink(wins_path, TRUE)
      arrow::write_parquet(
        data.frame(row = integer(), col = integer(), nrows = integer(), ncols = integer()),
        wins_path, compression = "zstd"
      )
    }

    return(invisible(sr_zones(
      region_id       = region_id,
      blocksize       = blocksize,
      zone_dir        = outdir,
      grid_index_path = grid_index_path,
      zidx_path       = zidx_path,
      wins_path       = if (isTRUE(build_windows)) wins_path else NULL
    )))
  }

  # ---- 4) Build Zidx by app + match (compute first, then write once) ----
  if (!all(c("GRIDCODE", "idx") %in% names(ids_tbl))) {
    stop("Index parquet lacks GRIDCODE/idx columns: ", grid_index_path, call. = FALSE)
  }

  # Numeric vector avoids integer64 issues and matches terra::app input type
  ids <- as.numeric(ids_tbl$GRIDCODE)

  # Compute Zidx without writing during app(); faster than streaming to disk
  Zidx <- terra::app(
    work_trim,
    fun = function(x) {
      out <- match(x, ids)
      out[is.na(x)] <- NA_integer_
      as.integer(out)
    }
  )

  # Now write once with compression and tiling
  if (file.exists(zidx_path)) .sr_safe_unlink(zidx_path, TRUE)
  gdal_opts <- c(
    "TILED=YES", "COMPRESS=ZSTD", "ZSTD_LEVEL=9",
    "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS",
    sprintf("BLOCKXSIZE=%d", blocksize),
    sprintf("BLOCKYSIZE=%d", blocksize)
  )
  terra::writeRaster(Zidx, zidx_path, overwrite = TRUE, datatype = "INT4S", gdal = gdal_opts)

  tick("Zidx computed and written")

  # ---- 5) Build non-empty windows from Zidx-on-disk (as in 01...R) ----
  if (isTRUE(build_windows)) {
    if (file.exists(wins_path) && !isTRUE(overwrite_rasters)) {
      message("Windows parquet exists, skipping: ", wins_path)
    } else {
      if (file.exists(wins_path)) .sr_safe_unlink(wins_path, TRUE)
      .sr_build_nonempty_windows_from_zidx_disk(
        zidx_path = zidx_path,
        blocksize = blocksize,
        out_parquet = wins_path
      )
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
  if (isTRUE(build_windows)) tick("windows built")
}
