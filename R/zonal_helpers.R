#' Ensure non-empty windows parquet exists for a zones object
#'
#' Builds windows if missing, writes to zones$wins_path (or default path), and returns updated zones.
#'
#' @keywords internal
#' @noRd
.sr_ensure_windows <- function(zones, progress_every = 500L, overwrite = FALSE) {
  if (!inherits(zones, "sr_zones")) {
    stop("`zones` must be an sr_zones object.", call. = FALSE)
  }

  # Decide where windows parquet should live
  wins_path <- zones$wins_path
  if (is.null(wins_path) || !nzchar(wins_path)) {
    wins_path <- file.path(
      zones$zone_dir,
      sprintf("%s_nonempty_windows_%d.parquet", zones$region_id, zones$blocksize)
    )
  }

  # If already there and not overwriting, just attach and return
  if (file.exists(wins_path) && !isTRUE(overwrite)) {
    zones$wins_path <- wins_path
    return(zones)
  }

  # Build windows and cache
  Zidx <- terra::rast(zones$zidx_path)

  # Use the fast occupancy aggregate method (your existing approach)
  nr <- terra::nrow(Zidx); nc <- terra::ncol(Zidx)

  occ <- terra::aggregate(
    !is.na(Zidx),
    fact = c(zones$blocksize, zones$blocksize),
    fun  = function(x) as.integer(any(as.logical(x)))
  )

  vals <- terra::values(occ, mat = FALSE)

  if (any(vals == 1L, na.rm = TRUE)) {
    cells <- which(vals == 1L)
    rc <- terra::rowColFromCell(occ, cells)
    r0 <- (rc[, 1L] - 1L) * zones$blocksize + 1L
    c0 <- (rc[, 2L] - 1L) * zones$blocksize + 1L
    nrw <- pmin(zones$blocksize, nr - r0 + 1L)
    ncw <- pmin(zones$blocksize, nc - c0 + 1L)
    wins_dt <- data.frame(row = r0, col = c0, nrows = nrw, ncols = ncw)
  } else {
    wins_dt <- data.frame(row = integer(), col = integer(), nrows = integer(), ncols = integer())
  }

  arrow::write_parquet(wins_dt, wins_path, compression = "zstd")

  zones$wins_path <- wins_path
  zones
}

# ------------- helper: unique GRIDCODE -> idx parquet (streamed) -------------
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

# ------------- helper: build non-empty windows parquet -----------------------
# Scans Zidx in tiles of size win_nrows x win_ncols; records windows with any non-NA.
.sr_build_nonempty_windows <- function(Zidx,
                                       win_nrows = 512L,
                                       win_ncols = 512L,
                                       out_parquet = NULL,
                                       progress_every = 500L) {
  stopifnot(inherits(Zidx, "SpatRaster"))
  stopifnot(terra::nlyr(Zidx) == 1L)

  nr <- terra::nrow(Zidx)
  nc <- terra::ncol(Zidx)

  # window grid
  row_starts <- seq.int(1L, nr, by = win_nrows)
  col_starts <- seq.int(1L, nc, by = win_ncols)

  # all windows as (row, col, nrows, ncols)
  win <- data.table::CJ(row = row_starts, col = col_starts)
  win[, nrows := pmin(win_nrows, nr - row + 1L)]
  win[, ncols := pmin(win_ncols, nc - col + 1L)]
  win[, win_id := .I]

  keep <- logical(nrow(win))

  terra::readStart(Zidx)
  on.exit(terra::readStop(Zidx), add = TRUE)

  for (i in seq_len(nrow(win))) {
    v <- terra::readValues(
      Zidx,
      row = win$row[i], nrows = win$nrows[i],
      col = win$col[i], ncols = win$ncols[i],
      mat = FALSE
    )
    keep[i] <- any(!is.na(v))

    if (progress_every > 0L && (i %% progress_every == 0L)) {
      message(sprintf("  scanned windows %d / %d", i, nrow(win)))
    }
  }

  out <- win[keep]

  if (!is.null(out_parquet)) {
    arrow::write_parquet(out, out_parquet, compression = "zstd")
  }

  out
}

.sr_auto_parallel_layout <- function(n_tasks,
                                     reserve_cores = 2L,
                                     prefer_threads_per_worker = 2L,
                                     max_workers = NULL,
                                     memfrac_total = 0.8,
                                     approx_mem_per_worker_gb = NULL) {
  total_cores  <- parallel::detectCores(logical = TRUE)
  usable_cores <- max(1L, total_cores - as.integer(reserve_cores))

  threads_per_worker <- max(1L, as.integer(prefer_threads_per_worker))
  workers_by_cpu <- floor(usable_cores / threads_per_worker)
  if (!is.null(max_workers)) workers_by_cpu <- min(workers_by_cpu, as.integer(max_workers))
  workers_by_cpu <- max(1L, workers_by_cpu)

  # Optional RAM cap (guarded; works even if ps doesn't export virtual_memory)
  workers_by_ram <- Inf
  available_gb <- NA_real_

  if (requireNamespace("ps", quietly = TRUE)) {
    ex <- getNamespaceExports("ps")
    mem_fun <- NULL
    if ("ps_system_memory" %in% ex) {
      mem_fun <- ps::ps_system_memory
    } else if ("virtual_memory" %in% ex) {
      mem_fun <- getExportedValue("ps", "virtual_memory")
    }
    if (!is.null(mem_fun)) {
      vm <- tryCatch(mem_fun(), error = function(e) NULL)
      if (is.list(vm)) {
        bytes <- vm$available
        if (is.null(bytes)) bytes <- vm$free
        if (is.numeric(bytes) && is.finite(bytes)) {
          available_gb <- as.numeric(bytes) / 1024^3
        }
      }
    }
  }

  if (!is.null(approx_mem_per_worker_gb) && is.finite(available_gb)) {
    max_workers_from_ram <- floor((available_gb * memfrac_total) / approx_mem_per_worker_gb)
    workers_by_ram <- max(1L, max_workers_from_ram)
  }

  workers <- min(n_tasks, workers_by_cpu, workers_by_ram)
  workers <- max(1L, as.integer(workers))

  terra_memfrac <- max(0.05, min(0.9, memfrac_total / workers))
  list(
    workers = workers,
    threads_per_worker = as.integer(threads_per_worker),
    gdal_threads  = as.integer(threads_per_worker),
    omp_threads   = as.integer(threads_per_worker),
    terra_threads = as.integer(threads_per_worker),
    terra_memfrac = terra_memfrac
  )
}
