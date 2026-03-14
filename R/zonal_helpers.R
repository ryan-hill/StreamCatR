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

  # Build windows (single source of truth)
  .sr_build_nonempty_windows_from_zidx_disk(
    zidx_path = zones$zidx_path,
    blocksize = zones$blocksize,
    out_parquet = wins_path
  )

  zones$wins_path <- wins_path
  zones
}


# ------------- helper: unique GRIDCODE -> idx parquet (streamed) -------------
# Helper function (no longer used by sr_optimize_zones, kept for backward compatibility)
.sr_build_gridcode_index_parquet_stream <- function(Z, out_file, progress_every = 0L) {
  stopifnot(inherits(Z, "SpatRaster"))

  # Early guard: empty raster (0 rows/cols) -> write empty parquet and return
  if (terra::nrow(Z) == 0L || terra::ncol(Z) == 0L || terra::ncell(Z) == 0L) {
    dt <- data.frame(GRIDCODE = integer(0), idx = integer(0))
    arrow::write_parquet(dt, out_file, compression = "zstd")
    return(invisible(dt))
  }

  bs <- terra::blocks(Z)

  # Normalize blocks() return: support both data.frame and list forms
  if (is.null(bs)) {
    # fallback: write empty and return
    dt <- data.frame(GRIDCODE = integer(0), idx = integer(0))
    arrow::write_parquet(dt, out_file, compression = "zstd")
    return(invisible(dt))
  }
  if (is.list(bs) && !is.data.frame(bs)) {
    bs <- as.data.frame(bs, stringsAsFactors = FALSE)
  }

  # Compute number of blocks robustly
  n_blocks <- if (!is.null(nrow(bs))) nrow(bs) else {
    if (!is.null(bs$row)) length(bs$row) else NA_integer_
  }

  # If still unknown or zero, treat as empty and return an empty parquet
  if (is.na(n_blocks) || n_blocks <= 0L) {
    dt <- data.frame(GRIDCODE = integer(0), idx = integer(0))
    arrow::write_parquet(dt, out_file, compression = "zstd")
    return(invisible(dt))
  }

  # Hash set of seen GRIDCODEs
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

  ids_chr <- ls(seen, all.names = TRUE)
  if (length(ids_chr) == 0L) {
    dt <- data.frame(GRIDCODE = integer(0), idx = integer(0))
    arrow::write_parquet(dt, out_file, compression = "zstd")
    return(invisible(dt))
  }

  # Keep the original behavior: integer GRIDCODEs in ascending order
  ids_int <- sort(as.integer(ids_chr))
  dt <- data.frame(GRIDCODE = ids_int, idx = seq_along(ids_int))
  arrow::write_parquet(dt, out_file, compression = "zstd")
  invisible(dt)
}

# ------------- helper: build non-empty windows parquet -----------------------
# Scans Zidx in tiles of size win_nrows x win_ncols; records windows with any non-NA.
.sr_build_nonempty_windows_from_zidx_disk <- function(zidx_path, blocksize, out_parquet) {
  Zidx_disk <- terra::rast(zidx_path)
  nr <- terra::nrow(Zidx_disk)
  nc <- terra::ncol(Zidx_disk)

  # Reduce by blocksize using a NA-safe presence test
  # sum(x, na.rm=TRUE) > 0 is robust even if x contains NA
  occ <- terra::aggregate(
    !is.na(Zidx_disk),
    fact = c(as.integer(blocksize), as.integer(blocksize)),
    fun  = function(x) as.integer(any(x != 0))   # fast and NA-safe
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

  arrow::write_parquet(wins_dt, out_parquet, compression = "zstd")
  invisible(out_parquet)
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
      mem_fun <- getExportedValue("ps", "ps_system_memory")
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

# Remove target and common GDAL sidecars when overwriting
.sr_safe_unlink <- function(path, do = TRUE) {
  if (!isTRUE(do)) return(invisible(FALSE))
  sidecars <- c("", ".aux.xml", ".ovr", ".msk", ".msk.aux.xml")
  files <- file.path(dirname(path), paste0(basename(path), sidecars))
  suppressWarnings(file.remove(files[file.exists(files)]))
  invisible(TRUE)
}


#' @keywords internal
#' @noRd
.sr_same_crs <- function(x, crs) {
  if (exists("same.crs", where = asNamespace("terra"), inherits = FALSE)) {
    return(terra::same.crs(x, crs))
  }
  identical(terra::crs(x), crs)
}

#' @keywords internal
#' @noRd
.sr_init_accumulator <- function(prefer_map = TRUE) {
  ns <- asNamespace("StreamCatR")

  has_idx <- exists("acc_sum_n_idx1K", envir = ns, mode = "function")
  has_map <- exists("acc_sum_n_map_centers1K", envir = ns, mode = "function")

  if (!has_idx) {
    stop(
      "Internal accumulator not available. Expected compiled acc_sum_n_idx1K in StreamCatR.",
      call. = FALSE
    )
  }

  invisible(list(
    acc_idx = get("acc_sum_n_idx1K", envir = ns, mode = "function"),
    acc_map = if (isTRUE(prefer_map) && has_map) {
      get("acc_sum_n_map_centers1K", envir = ns, mode = "function")
    } else {
      NULL
    }
  ))
}

#' @keywords internal
#' @noRd
.sr_parse_blocksize_from_filename <- function(path) {
  # For legacy names like <rid>_zone_index_blockNNNN.tif ->  NNNN
  m <- regexec("_block([0-9]+)\\.tif$", basename(path))
  reg <- regmatches(basename(path), m)[[1]]
  if (length(reg) >= 2) {
    as.integer(reg[2])
  } else {
    NA_integer_
  }
}

#' Prefer "new" naming, with clean fallback to "legacy"
#' Returns list(paths...) and inferred blocksize if from legacy names.
#' @keywords internal
#' @noRd
.sr_guess_zone_paths <- function(region_id, zone_dir) {
  # "new" convention
  new_zidx <- file.path(zone_dir, sprintf("%s_zone_index.tif", region_id))
  new_idx  <- file.path(zone_dir, sprintf("%s_rasterid_index.parquet", region_id))
  new_win  <- file.path(zone_dir, sprintf("%s_nonempty_windows.parquet", region_id))

  # "legacy" convention (blocksize suffix)
  # We don't know blocksize, so search for matching patterns and pick the most recent or any one
  legacy_zidx <- Sys.glob(file.path(zone_dir, sprintf("%s_zone_index_block*.tif", region_id)))
  legacy_idx  <- file.path(zone_dir, sprintf("%s_gridcode_index.parquet", region_id))
  legacy_win  <- Sys.glob(file.path(zone_dir, sprintf("%s_nonempty_windows_*.parquet", region_id)))

  if (file.exists(new_zidx) && file.exists(new_idx)) {
    list(
      zidx_path       = new_zidx,
      grid_index_path = new_idx,
      wins_path       = if (file.exists(new_win)) new_win else NULL,
      blocksize       = NA_integer_,  # not encoded in new names
      naming          = "new"
    )
  } else if (length(legacy_zidx) > 0L && file.exists(legacy_idx)) {
    # choose the first .tif (or the one with largest blocksize); match wins by blocksize if possible
    # pick the tif with the largest blocksize (usually the intentional one)
    bs_all <- vapply(legacy_zidx, .sr_parse_blocksize_from_filename, integer(1))
    pick <- if (any(!is.na(bs_all))) which.max(bs_all) else 1L
    zidx <- legacy_zidx[pick]
    bs   <- .sr_parse_blocksize_from_filename(zidx)

    win <- NULL
    if (length(legacy_win) > 0L) {
      if (!is.na(bs)) {
        # try to match wins file with same blocksize
        pat <- sprintf("_nonempty_windows_%d\\.parquet$", bs)
        sel <- grep(pat, legacy_win)
        if (length(sel) >= 1L) win <- legacy_win[sel[1]]
      }
      # fallback: just take the first windows file
      if (is.null(win)) win <- legacy_win[1]
    }

    list(
      zidx_path       = zidx,
      grid_index_path = legacy_idx,
      wins_path       = if (!is.null(win) && file.exists(win)) win else NULL,
      blocksize       = bs,
      naming          = "legacy"
    )
  } else {
    list(
      zidx_path       = new_zidx,
      grid_index_path = new_idx,
      wins_path       = if (file.exists(new_win)) new_win else NULL,
      blocksize       = NA_integer_,
      naming          = "none"
    )
  }
}


#' Infer blocksize from a windows parquet by using the most frequent window size
#' Interior tiles are full-size; edges are smaller.
#' @keywords internal
#' @noRd
.sr_blocksize_from_windows <- function(wins_path) {
  if (is.null(wins_path) || !file.exists(wins_path)) return(NA_integer_)
  wins <- tryCatch(arrow::read_parquet(wins_path, as_data_frame = TRUE),
                   error = function(e) NULL)
  if (is.null(wins) || nrow(wins) == 0L) return(NA_integer_)

  nr <- as.integer(wins$nrows)
  nc <- as.integer(wins$ncols)
  nr <- nr[is.finite(nr) & nr > 0L]
  nc <- nc[is.finite(nc) & nc > 0L]
  if (!length(nr) && !length(nc)) return(NA_integer_)

  # Mode helper
  mode_int <- function(x) {
    tx <- table(x)
    as.integer(names(tx)[which.max(tx)])
  }

  mr <- if (length(nr)) mode_int(nr) else NA_integer_
  mc <- if (length(nc)) mode_int(nc) else NA_integer_

  if (is.finite(mr) && is.finite(mc)) {
    if (mr == mc) mr else max(mr, mc)
  } else if (is.finite(mr)) {
    mr
  } else if (is.finite(mc)) {
    mc
  } else {
    NA_integer_
  }
}
