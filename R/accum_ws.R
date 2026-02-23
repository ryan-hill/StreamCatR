#' Accumulate zonal statistics upstream across a river network
#'
#' Given zonal outputs (e.g., from [sr_zonal_region()] or [sr_zonal_all()]) and a staged
#' upstream-pairs dataset (built by [sr_build_accum_framework()]), compute upstream (watershed)
#' accumulated statistics.
#'
#' Supports optional translation from raster GRIDCODE IDs to network IDs (e.g., COMID/FEATUREID).
#'
#' @param zonal Zonal results as a data.frame/data.table, Arrow table/dataset, or a file path to parquet.
#' @param framework Either a `pairs_dir` path or an `"sr_accum_framework"` object returned by
#'   [sr_build_accum_framework()].
#' @param translation Optional translation table (data.frame/Arrow table) or parquet path.
#'   Used to map `zonal_id_col` to `net_id_col` (e.g., GRIDCODE -> FEATUREID/COMID).
#'   If NULL, assumes zonal IDs already match network IDs.
#' @param zonal_id_col Column name in `zonal` identifying zones (default `"GRIDCODE"`).
#' @param translation_from Column name in `translation` matching `zonal_id_col` (default `"GRIDCODE"`).
#' @param translation_to Column name in `translation` providing network IDs (default `"FEATUREID"`).
#' @param net_id_col Output network ID column name (default `"COMID"`).
#' @param sum_col Name of the catchment-level sum column (default `"sum"`).
#' @param count_col Name of the catchment-level count column (default `"count"`).
#' @param mean Logical; if TRUE, compute mean as sum/count at both catchment and watershed scale.
#' @param ... Reserved for future options.
#'
#' @return A data.table with one row per network ID and columns for catchment and watershed metrics.
#' @export
sr_accumulate_ws <- function(
    zonal,
    framework,
    translation = NULL,
    zonal_id_col = "GRIDCODE",
    translation_from = "GRIDCODE",
    translation_to = "FEATUREID",
    net_id_col = "COMID",
    sum_col = "sum",
    count_col = "count",
    mean = TRUE,
    ...
) {
  pairs_dir <- if (inherits(framework, "sr_accum_framework")) framework$pairs_dir else framework
  if (!is.character(pairs_dir) || length(pairs_dir) != 1L) {
    stop("`framework` must be a pairs_dir path or an sr_accum_framework object.", call. = FALSE)
  }
  if (!dir.exists(pairs_dir)) stop("pairs_dir does not exist: ", pairs_dir, call. = FALSE)

  z <- .sr_read_table(zonal)

  # If translation provided (or embedded in framework schema), map to network IDs
  if (is.null(translation) && inherits(framework, "sr_accum_framework")) {
    tp <- framework$schema$translation_path %||% NULL
    if (!is.null(tp) && file.exists(tp)) translation <- tp
  }

  if (!is.null(translation)) {
    tr <- .sr_read_table(translation)

    # Join zonal -> translation -> network ID
    z <- dplyr::left_join(
      z,
      tr[, c(translation_from, translation_to)],
      by = setNames(translation_from, zonal_id_col)
    )
    if (!translation_to %in% names(z)) {
      stop("Translation join failed: missing column ", translation_to, call. = FALSE)
    }
    z[[net_id_col]] <- z[[translation_to]]

    # Aggregate if multiple zonal IDs map to same net ID
    z <- z |>
      dplyr::filter(!is.na(.data[[net_id_col]])) |>
      dplyr::group_by(.data[[net_id_col]]) |>
      dplyr::summarise(
        !!sum_col   := sum(.data[[sum_col]], na.rm = TRUE),
        !!count_col := sum(.data[[count_col]], na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    if (!net_id_col %in% names(z)) {
      # If no translation, assume zonal_id_col is the network id
      z[[net_id_col]] <- z[[zonal_id_col]]
    }
    z <- z |>
      dplyr::group_by(.data[[net_id_col]]) |>
      dplyr::summarise(
        !!sum_col   := sum(.data[[sum_col]], na.rm = TRUE),
        !!count_col := sum(.data[[count_col]], na.rm = TRUE),
        .groups = "drop"
      )
  }

  # Read staged upstream pairs dataset
  # Expect a dataset with columns: ID, UPCOMID (names configurable later)
  pairs <- arrow::open_dataset(pairs_dir, format = "parquet")
  # Convention (recommended): columns `COMID` and `UPCOMID`
  # If you choose different names, store them in framework$schema and use here.

  # Do the accumulation (placeholder: implement in helpers)
  out <- .sr_accumulate_with_pairs(
    cat_tbl = z,
    pairs = pairs,
    id_col = net_id_col,
    sum_col = sum_col,
    count_col = count_col
  )

  if (isTRUE(mean)) {
    out[["cat_mean"]] <- with(out, ifelse(cat_count > 0, cat_sum / cat_count, NA_real_))
    out[["ws_mean"]]  <- with(out, ifelse(ws_count  > 0, ws_sum  / ws_count,  NA_real_))
  }

  out
}
