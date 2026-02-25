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

  # Reproject if needed; write a single temp file with tiling/compression
  tmp_proj <- NULL
  if (!terra::same.crs(src, target_crs)) {
    tmp_proj <- tempfile(fileext = ".tif")
    src <- terra::project(
      src, target_crs, method = "near",
      filename = tmp_proj, overwrite = TRUE,
      gdal = c(
        "TILED=YES", "COMPRESS=ZSTD", "ZSTD_LEVEL=9",
        "BIGTIFF=YES", "NUM_THREADS=ALL_CPUS"
      )
    )
  }

  # Cleanup of temp projection on exit
  on.exit({
    if (!is.null(tmp_proj) && file.exists(tmp_proj)) suppressWarnings(file.remove(tmp_proj))
  }, add = TRUE)

  # Trim to remove exterior NA rows/cols; drop factor levels if present
  work_trim <- try(terra::trim(src), silent = TRUE)
  if (inherits(work_trim, "try-error")) work_trim <- src
  if (terra::is.factor(work_trim)) {
    try({ levels(work_trim) <- NULL }, silent = TRUE)
  }

  # Build or read GRIDCODE -> idx mapping
  ids_tbl <- NULL
  if (file.exists(grid_index_path) && !overwrite_rasters) {
    message("Index exists, skipping rebuild: ", grid_index_path)
    ids_tbl <- arrow::read_parquet(grid_index_path, as_data_frame = TRUE)
    if (!all(c("GRIDCODE", "idx") %in% names(ids_tbl))) {
      stop("Index parquet lacks GRIDCODE/idx columns: ", grid_index_path, call. = FALSE)
    }
  } else {
    # Get unique GRIDCODEs; older terra versions do not accept useNA argument
    fq <- suppressWarnings(terra::freq(work_trim))
    fq_df <- as.data.frame(fq)

    # Expect a column named "value" with unique values
    if (is.null(fq_df) || nrow(fq_df) == 0L || !("value" %in% names(fq_df))) {
      # Empty or all-NA raster
      ids_tbl <- data.frame(GRIDCODE = integer(), idx = integer())
      if (file.exists(grid_index_path)) safe_unlink(grid_index_path, TRUE)
      arrow::write_parquet(ids_tbl, grid_index_path, compression = "zstd")

      # Write empty Zidx with same geometry
      if (!file.exists(zidx_path) || overwrite_rasters) {
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

      if (isTRUE(build_windows) && (!file.exists(wins_path) || overwrite_rasters)) {
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

    # Drop NA values if present
    vals <- fq_df[["value"]]
    vals <- vals[!is.na(vals)]

    if (!length(vals)) {
      # All NA case
      ids_tbl <- data.frame(GRIDCODE = integer(), idx = integer())
      if (file.exists(grid_index_path)) safe_unlink(grid_index_path, TRUE)
      arrow::write_parquet(ids_tbl, grid_index_path, compression = "zstd")

      if (!file.exists(zidx_path) || overwrite_rasters) {
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

      if (isTRUE(build_windows) && (!file.exists(wins_path) || overwrite_rasters)) {
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

    # Avoid integer64 surprises: use integer if within range, else numeric
    if (max(abs(vals), na.rm = TRUE) <= .Machine$integer.max) {
      gridcodes <- as.integer(vals)
    } else {
      gridcodes <- as.numeric(vals)
    }

    gridcodes <- sort(unique(gridcodes))
    ids_tbl <- data.frame(GRIDCODE = gridcodes, idx = seq_along(gridcodes))

    # Write parquet (zstd)
    if (file.exists(grid_index_path)) safe_unlink(grid_index_path, TRUE)
    arrow::write_parquet(ids_tbl, grid_index_path, compression = "zstd")
  }

  # Build Zidx: prefer subs if available; else fall back to app+match (or fastmatch)
  has_subs <- isTRUE("subs" %in% getNamespaceExports("terra"))
  if (file.exists(zidx_path) && !overwrite_rasters) {
    message("Zidx exists, skipping: ", zidx_path)
  } else {
    if (file.exists(zidx_path)) safe_unlink(zidx_path, TRUE)
    if (has_subs) {
      # Use subs when available (some terra versions)
      lut <- data.frame(from = ids_tbl$GRIDCODE, to = ids_tbl$idx)
      Zidx <- terra::subs(work_trim, lut, by = "from", which = "to")
    } else {
      # Fallback: app with vectorized R function; optionally use fastmatch for speed
      use_fastmatch <- requireNamespace("fastmatch", quietly = TRUE)
      ids_vec <- ids_tbl$GRIDCODE
      if (use_fastmatch) {
        ffun <- function(x) {
          out <- fastmatch::fmatch(x, ids_vec)
          out[is.na(x)] <- NA_integer_
          as.integer(out)
        }
      } else {
        ffun <- function(x) {
          out <- match(x, ids_vec)
          out[is.na(x)] <- NA_integer_
          as.integer(out)
        }
      }
      Zidx <- terra::app(work_trim, fun = ffun)
    }

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

  # Compute non-empty windows directly from work_trim with "sum" reducer
  wins_out <- NULL
  if (isTRUE(build_windows)) {
    if (file.exists(wins_path) && !overwrite_rasters) {
      message("Windows parquet exists, skipping: ", wins_path)
      wins_out <- wins_path
    } else {
      if (file.exists(wins_path)) safe_unlink(wins_path, TRUE)

      nr <- terra::nrow(work_trim)
      nc <- terra::ncol(work_trim)

      occ <- terra::aggregate(!is.na(work_trim), fact = c(blocksize, blocksize), fun = "sum")

      vals <- terra::values(occ, mat = FALSE)
      if (any(vals > 0L, na.rm = TRUE)) {
        cells <- which(vals > 0L)
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

# Helper function (no longer used by sr_optimize_zones, kept for backward compatibility)
.sr_build_gridcode_index_parquet_stream <- function(Z, out_file, progress_every = 0L) {
  stopifnot(inherits(Z, "SpatRaster"))
  bs <- terra::blocks(Z)
  n_blocks <- if (!is.null(nrow(bs))) nrow(bs) else length(bs$row)
  if (is.null(n_blocks) || n_blocks <= 0L) stop("Unexpected blocks() structure")

  # Use an environment as a hash set to avoid repeated unions
  seen <- new.env(parent = emptyenv())

  terra::readStart(Z)
  on.exit(terra::readStop(Z), add = TRUE)

  for (i in seq_len(n_blocks)) {
    b_row   <- bs$row[i]
    b_nrows <- bs$nrows[i]
    z <- terra::readValues(Z, row = b_row, nrows = b_nrows, mat = FALSE)
    z <- z[!is.na(z)]
    if (!length(z)) next
    uz <- unique(z)
    for (v in uz) seen[[as.character(v)]] <- TRUE

    if (progress_every > 0L && (i %% progress_every == 0L)) {
      message(sprintf("  scanned block %d / %d", i, n_blocks))
    }
  }

  ids <- sort(as.integer(ls(seen)))
  dt  <- data.frame(GRIDCODE = ids, idx = seq_along(ids))
  arrow::write_parquet(dt, out_file, compression = "zstd")
  invisible(out_file)
}
