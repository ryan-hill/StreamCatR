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
#' @param stats Character vector of statistics. Supported: `"sum"`, `"mean"`, `"count"` (alias: `"n"`).
#'
#' @param plan_strategy Future strategy name (e.g., `"sequential"`, `"multisession"`).
#' @param reserve_cores Integer number of CPU cores to leave unused.
#' @param prefer_threads_per_worker Logical or integer tuning value used by the internal
#'   parallel-layout heuristic.
#' @param max_workers Maximum number of parallel workers to use.
#' @param memfrac_total Fraction of total RAM to budget for workers (0-1).
#' @param approx_mem_per_worker_gb Approximate memory (GB) needed per worker (used for worker count heuristics).
#'
#' @param out_terra_temp_base Optional directory for terra temp files (useful on HPC / scratch).
#' @param r_libs_user Optional library path to set `R_LIBS_USER` for workers.
#' @param combine Logical; if `TRUE`, row-bind all regional results into one table.
#' @param verbose Logical; if `TRUE`, print progress and configuration details.
#'
#' @return If `combine=TRUE`, a single `data.table` when `data.table` is available, otherwise a `data.frame`.
#'   If `combine=FALSE`, a named list of per-region results.
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
                         r_libs_user = NULL,
                         combine = TRUE,
                         verbose = TRUE) {
  stopifnot(length(regions) > 0)
  plan_strategy <- match.arg(plan_strategy)
  blocksize <- as.integer(blocksize)

  if (!requireNamespace("future", quietly = TRUE) ||
      !requireNamespace("future.apply", quietly = TRUE)) {
    stop("Please install.packages(c('future','future.apply'))")
  }

  layout <- .sr_auto_parallel_layout(
    n_tasks = length(regions),
    reserve_cores = reserve_cores,
    prefer_threads_per_worker = prefer_threads_per_worker,
    max_workers = max_workers,
    memfrac_total = memfrac_total,
    approx_mem_per_worker_gb = approx_mem_per_worker_gb
  )

  if (verbose) {
    message(sprintf(
      "Parallel layout: workers=%d, threads/worker=%d, terra_memfrac=%.3f",
      layout$workers, layout$threads_per_worker, layout$terra_memfrac
    ))
  }

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

  dir.create(out_terra_temp_base, recursive = TRUE, showWarnings = FALSE)

  win_path <- function(rid) {
    file.path(zone_dir, sprintf("%s_nonempty_windows_%d.parquet", rid, blocksize))
  }

  has_win <- vapply(regions, function(r) file.exists(win_path(r)), logical(1))
  if (!all(has_win)) {
    warning("Skipping regions without windows/parquet: ", paste(regions[!has_win], collapse = ", "))
    regions <- regions[has_win]
  }
  if (!length(regions)) {
    stop("No regions to process after pre-checks.")
  }

  res_list <- future.apply::future_lapply(
    regions,
    FUN = function(rid) {
      if (!is.null(r_libs_user)) {
        .libPaths(c(r_libs_user, .libPaths()))
      }

      Sys.setenv(
        GDAL_NUM_THREADS = as.character(layout$gdal_threads),
        OMP_NUM_THREADS = as.character(layout$omp_threads),
        OPENBLAS_NUM_THREADS = as.character(layout$threads_per_worker),
        MKL_NUM_THREADS = as.character(layout$threads_per_worker),
        BLIS_NUM_THREADS = as.character(layout$threads_per_worker),
        VECLIB_MAXIMUM_THREADS = as.character(layout$threads_per_worker),
        GDAL_CACHEMAX = "2048"
      )

      wd <- file.path(out_terra_temp_base, paste0("worker_", Sys.getpid()))
      dir.create(wd, recursive = TRUE, showWarnings = FALSE)

      terra::terraOptions(
        threads = layout$terra_threads,
        memfrac = layout$terra_memfrac,
        tempdir = wd
      )

      zones_obj <- sr_zones(
        region_id = rid,
        blocksize = blocksize,
        zone_dir = zone_dir,
        grid_index_path = file.path(zone_dir, sprintf("%s_gridcode_index.parquet", rid)),
        zidx_path = file.path(zone_dir, sprintf("%s_zone_index_block%d.tif", rid, blocksize)),
        wins_path = file.path(zone_dir, sprintf("%s_nonempty_windows_%d.parquet", rid, blocksize))
      )

      out <- sr_zonal_region(
        predictor_path = predictor_path,
        zones = zones_obj,
        method = method,
        stats = stats,
        pad_cells = 1L,
        exact_crop = TRUE,
        progress_every = 0L
      )

      list(region = rid, result = out)
    },
    future.packages = c("StreamCatR", "terra", "arrow", "data.table")
  )

  names(res_list) <- vapply(res_list, `[[`, "", "region")

  if (!combine) {
    return(res_list)
  }

  all_df <- all(vapply(res_list, function(x) inherits(x$result, "data.frame"), logical(1)))

  if (all_df && requireNamespace("data.table", quietly = TRUE)) {
    return(data.table::rbindlist(
      lapply(res_list, `[[`, "result"),
      fill = TRUE,
      use.names = TRUE
    ))
  }

  if (all_df) {
    out <- do.call(rbind, lapply(res_list, `[[`, "result"))
    rownames(out) <- NULL
    return(out)
  }

  res_list
}
