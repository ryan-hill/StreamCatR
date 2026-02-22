#' Optimize a predictor raster for fast zonal processing
#'
#' Reads a predictor raster, optionally projects/resamples it to a target CRS and
#' resolution, optionally trims NA-only edges, and writes a tiled/compressed GeoTIFF
#' suitable for fast windowed IO during zonal summaries.
#'
#' @param infile Path to an input raster readable by `terra`.
#' @param outdir Output directory for the optimized raster.
#' @param outname Base filename (without directory) for the optimized raster.
#' @param blocksize Integer tile size in pixels (used for GeoTIFF BLOCKXSIZE/BLOCKYSIZE).
#' @param target_crs Target CRS for output (e.g., `"EPSG:5070"`). If NULL, keep original.
#' @param resample Resampling method used when projecting/resampling (e.g., `"near"`, `"bilinear"`).
#' @param trim_edges Logical; if TRUE, trim NA-only edges using `terra::trim()` when possible.
#' @param datatype GDAL datatype string for output (e.g., `"FLT4S"`, `"INT2S"`).
#' @param overwrite Logical; if TRUE, overwrite existing output and common sidecar files.
#' @param gdal_compress Compression method (e.g., `"ZSTD"`, `"LZW"`).
#' @param zstd_level Integer ZSTD compression level (only used when `gdal_compress="ZSTD"`).
#' @param bigtiff Logical; if TRUE, write BIGTIFF=YES.
#' @param threads Character/Integer passed to GDAL `NUM_THREADS` (e.g., `"ALL_CPUS"`).
#' @param gdal_extra Character vector of additional GDAL creation options.
#'
#' @return Invisibly returns the output file path.
#' @export
sr_optimize_predictor <- function(
    infile,
    outdir        = dirname(infile),
    outname       = NULL,
    blocksize     = 512L,
    target_crs    = "EPSG:5070",  # ensure output is EPSG:5070
    resample      = "bilinear",   # bilinear for continuous rasters
    trim_edges    = TRUE,         # trim NA-only edges if present
    datatype      = "FLT4S",      # 32-bit float (use "FLT8S" for doubles)
    overwrite     = TRUE,
    gdal_compress = "ZSTD",       # "ZSTD", "LZW", or "DEFLATE"
    zstd_level    = 9,
    bigtiff       = TRUE,
    threads       = "ALL_CPUS",
    gdal_extra    = NULL          # e.g., c("PREDICTOR=3") for LZW/DEFLATE
) {
  stopifnot(file.exists(infile))
  blocksize <- as.integer(blocksize)

  # Load
  R <- terra::rast(infile)

  # Ensure EPSG:5070 (skip identity reprojection)
  if (!is.null(target_crs) && !terra::same.crs(R, target_crs)) {
    R <- terra::project(R, target_crs, method = resample)
  } else {
    # CRS already equivalent; optionally standardize label without resampling
    # (wrapped in try to avoid warnings on some builds)
    try(suppressWarnings(terra::crs(R) <- target_crs), silent = TRUE)
  }

  # Trim NA-only edges (skip gracefully if not applicable)
  if (isTRUE(trim_edges)) {
    R_trim <- try(terra::trim(R), silent = TRUE)
    if (!inherits(R_trim, "try-error")) {
      R <- R_trim
    } else {
      message("trim() skipped (likely all NA detected by trim or other issue).")
    }
  }

  # Output name and folder
  if (is.null(outname)) {
    base <- tools::file_path_sans_ext(basename(infile))
    outname <- sprintf("%s_%d.tif", base, blocksize)
  }
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  outfile <- file.path(outdir, outname)

  # GDAL options
  gdal_opts <- c(
    "TILED=YES",
    paste0("COMPRESS=", toupper(gdal_compress)),
    if (toupper(gdal_compress) == "ZSTD") paste0("ZSTD_LEVEL=", zstd_level) else NULL,
    paste0("BIGTIFF=", if (bigtiff) "YES" else "IF_NEEDED"),
    paste0("NUM_THREADS=", threads),
    paste0("BLOCKXSIZE=", blocksize),
    paste0("BLOCKYSIZE=", blocksize),
    gdal_extra
  )

  # Write
  terra::writeRaster(
    R,
    filename  = outfile,
    overwrite = overwrite,
    datatype  = datatype,
    gdal      = gdal_opts
  )

  invisible(outfile)
}
