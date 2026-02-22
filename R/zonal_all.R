#' Run zonal statistics for all regions
#'
#' Applies [sr_zonal_region()] across multiple regions and (optionally) combines results.
#' Designed for large rasters and many regions; can parallelize across regions.
#'
#' @param regions Character vector of region identifiers.
#' @param zone_dir Directory containing zone artifacts produced by [sr_optimize_zones()].
#' @param predictor_path Path to the predictor raster.
#' @param blocksize Integer block size used when building the zone artifacts.
#' @param method Resampling method (default `"near"`). Used only if projecting predictor.
#' @param stats Character vector of statistics. Supported: `"sum"`, `"mean"`, `"n"`.
#'
#' @param plan_strategy Future strategy name (e.g., `"sequential"`, `"multisession"`).
#' @param reserve_cores Integer number of CPU cores to leave unused.
#' @param prefer_threads_per_worker Logical; if TRUE, prefer fewer workers with more threads each.
#' @param max_workers Maximum number of parallel workers to use.
#' @param memfrac_total Fraction of total RAM to budget for workers (0-1).
#' @param approx_mem_per_worker_gb Approximate memory (GB) needed per worker (used for worker count heuristics).
#'
#' @param out_terra_temp_base Optional directory for terra temp files (useful on HPC / scratch).
#' @param r_libs_user Optional library path to set `R_LIBS_USER` for workers.
#' @param combine Logical; if TRUE, row-bind all regional results into one table.
#' @param verbose Logical; if TRUE, print progress and configuration details.
#'
#' @return If `combine=TRUE`, a single `data.table`. Otherwise a named list of `data.table`s by region.
#' @export
sr_zonal_all <- function(regions,
                                 zone_dir,
                                 predictor_path,
                                 blocksize = 6144L,
                                 method = "near",
                                 stats = "sum",
                                 plan_strategy = c("multisession", "multicore", "cluster"),
                                 reserve_cores = 2L,
                                 prefer_threads_per_worker = 2L,
                                 max_workers = NULL,
                                 memfrac_total = 0.8,
                                 approx_mem_per_worker_gb = NULL,
                                 out_terra_temp_base = file.path(tempdir(), "terra_parallel"),
                                 #zonal_file = file.path("StreamCatR","zonal_stats_r","04.zonal_functions_true_blocking.R"),
                                 r_libs_user = NULL,  # optional: e.g., "C:/Users/ME/AppData/Local/R/libraries"
                                 combine = TRUE,
                                 verbose = TRUE) {
  stopifnot(length(regions) > 0)
  plan_strategy <- match.arg(plan_strategy)

  # Ensure future/future.apply available
  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install.packages(c('future','future.apply'))")
  }

  # Resolve zonal file path now to fail fast
  # zf <- tryCatch(normalizePath(zonal_file, mustWork = TRUE),
  #                error = function(e) stop("zonal_file not found: ", zonal_file))

  layout <- .sr_auto_parallel_layout(
    n_tasks = length(regions),
    reserve_cores = reserve_cores,
    prefer_threads_per_worker = prefer_threads_per_worker,
    max_workers = max_workers,
    memfrac_total = memfrac_total,
    approx_mem_per_worker_gb = approx_mem_per_worker_gb
  )
  if (verbose) {
    message(sprintf("Parallel layout: workers=%d, threads/worker=%d, terra_memfrac=%.3f",
                    layout$workers, layout$threads_per_worker, layout$terra_memfrac))
  }

  # Set future plan
  if (identical(plan_strategy, "multicore") && !future::supportsMulticore()) {
    warning("multicore not supported; falling back to multisession")
    plan_strategy <- "multisession"
  }
  if (identical(plan_strategy, "multisession")) {
    future::plan(future::multisession, workers = layout$workers)
  } else if (identical(plan_strategy, "multicore")) {
    future::plan(future::multicore, workers = layout$workers)
  } else {
    future::plan(future::cluster, workers = layout$workers)
  }
  on.exit(future::plan(future::sequential), add = TRUE)

  # Temp dirs for workers
  dir.create(out_terra_temp_base, recursive = TRUE, showWarnings = FALSE)

  # Optional: pre-check required window/parquet files
  win_path <- function(rid) file.path(zone_dir, sprintf("%s_nonempty_windows_%s.parquet", rid, blocksize))
  has_win <- vapply(regions, function(r) file.exists(win_path(r)), logical(1))
  if (!all(has_win)) {
    warning("Skipping regions without windows/parquet: ", paste(regions[!has_win], collapse = ", "))
    regions <- regions[has_win]
  }
  if (!length(regions)) stop("No regions to process after pre-checks.")

  # res_list <- future.apply::future_lapply(
  #   regions,
  #   FUN = function(rid) {
  #     # Ensure workers use your preferred library (optional)
  #     if (!is.null(r_libs_user)) .libPaths(c(r_libs_user, .libPaths()))
  #
  #     # Per-worker thread/env settings
  #     Sys.setenv(GDAL_NUM_THREADS = as.character(layout$gdal_threads),
  #                OMP_NUM_THREADS  = as.character(layout$omp_threads))
  #     wd <- file.path(out_terra_temp_base, paste0("worker_", Sys.getpid()))
  #     dir.create(wd, recursive = TRUE, showWarnings = FALSE)
  #     terra::terraOptions(threads = layout$terra_threads,
  #                         memfrac = layout$terra_memfrac,
  #                         tempdir = wd)
  #
  #     # Load required packages (don’t hard-require scaccum; 04 handles fallback)
  #     suppressPackageStartupMessages({
  #       library(terra); library(arrow); library(data.table)
  #       if (!requireNamespace("collapse", quietly = TRUE)) {
  #         stop("Need collapse. install.packages('collapse')")
  #       }
  #     })
  #     # Optional: attach scaccum if present (init in 04 will also handle it)
  #       # library(scaccum) # not strictly necessary if 04 binds the symbol
  #     }
  #
  #     # Source zonal functions; this should init acc_sum_n_idx1K (package or fallback)
  #     sys.source(zf, envir = environment())
  #
  #     if (isTRUE(verbose)) message(sprintf("[Worker %d] %s", Sys.getpid(), rid))
  #
  #     out <- sr_zonal_region(
  #       region_id       = rid,
  #       zone_dir        = zone_dir,
  #       predictor_path  = predictor_path,
  #       blocksize       = blocksize,
  #       method          = method,
  #       stats           = stats,
  #       progress_every  = 0L
  #     )
  #     list(region = rid, result = out)
  #   },
  #   future.packages = c("terra","arrow","data.table")  # core deps only
  # )
  res_list <- future.apply::future_lapply(
    regions,
    FUN = function(rid) {
      if (!is.null(r_libs_user)) .libPaths(c(r_libs_user, .libPaths()))

      Sys.setenv(GDAL_NUM_THREADS = as.character(layout$gdal_threads),
                 OMP_NUM_THREADS  = as.character(layout$omp_threads))
      Sys.setenv(OPENBLAS_NUM_THREADS   = as.character(layout$threads_per_worker),
                 MKL_NUM_THREADS        = as.character(layout$threads_per_worker),
                 BLIS_NUM_THREADS       = as.character(layout$threads_per_worker),
                 VECLIB_MAXIMUM_THREADS = as.character(layout$threads_per_worker))
      Sys.setenv(GDAL_CACHEMAX = "2048")

      wd <- file.path(out_terra_temp_base, paste0("worker_", Sys.getpid()))
      dir.create(wd, recursive = TRUE, showWarnings = FALSE)
      terra::terraOptions(
        threads = layout$terra_threads,
        memfrac = layout$terra_memfrac,
        tempdir = wd
      )

      # Optional: scaccum
      # if (requireNamespace("scaccum", quietly = TRUE)) {
      #   # optional attach:
      #   # suppressPackageStartupMessages(library(scaccum))
      # }

      out <- sr_zonal_region(
        region_id       = rid,
        zone_dir        = zone_dir,
        predictor_path  = predictor_path,
        blocksize       = blocksize,
        method          = method,
        stats           = stats,
        progress_every  = 0L
      )
      list(region = rid, result = out)
    },
    future.packages = c("StreamCatR", "terra", "arrow", "data.table")
  )

  names(res_list) <- vapply(res_list, `[[`, "", "region")

  if (!combine) return(res_list)

  all_df <- all(vapply(res_list, function(x) inherits(x$result, "data.frame"), logical(1)))
  if (all_df && requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::rbindlist(lapply(res_list, `[[`, "result"), fill = TRUE, use.names = TRUE))
  } else if (all_df) {
    out <- do.call(rbind, lapply(res_list, `[[`, "result"))
    rownames(out) <- NULL
    return(out)
  } else {
    return(res_list)
  }
}
