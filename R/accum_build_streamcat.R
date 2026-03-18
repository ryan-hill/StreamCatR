#' Build a watershed-accumulation framework (upstream pairs)
#'
#' Builds a staged upstream-pairs dataset (stored as Arrow parquet files) from a
#' directed from-to topology table. Optionally uses a preset
#' (e.g., MR NHDPlus CONUS) that provides standard column names,
#' translation tables, and default edge cleaning.
#'
#' The resulting framework can be used with [sr_accumulate_ws()] to accumulate
#' zonal statistics upstream across the network.
#'
#' @param out_dir Output directory for the staged upstream-pairs dataset.
#' @param preset Preset configuration. Use `"mr_nhdplus"` for MR NHDPlus
#' conventions, or `"none"` for a generic from-to table.
#' @param nhdplus_mr Logical; if TRUE, equivalent to `preset="mr_nhdplus"`.
#' If FALSE, equivalent to `preset="none"`. If NULL, uses `preset`.
#' @param data_dir Directory containing preset parquet inputs (only for
#' presets).
#'   If NULL and a preset is used, data will be fetched on-demand.
#' @param from_to From-to topology table (data.frame/Arrow Table/Dataset)
#' or a file path. Required if `preset="none"`. Ignored for
#' `preset="mr_nhdplus"` unless explicitly provided.
#' @param edge_filter Optional function that takes a data.frame of edges and
#' returns a filtered data.frame. Applied after preset default cleaning.
#' @param include_self Logical; if TRUE, include each ID as upstream of itself
#' (common for accumulation).
#' @param overwrite Logical; if TRUE, overwrite existing `out_dir`.
#' @param ... Additional arguments passed to internal builders.
#'
#' @return An object of class `"sr_accum_framework"` containing `pairs_dir`,
#' `preset`, and schema metadata (including id columns used).
#' @export
sr_streamcat_framework <- function(
    out_dir,
    preset = c("none", "mr_nhdplus"),
    data_dir = NULL,           # can be path or a list with $paths like your streamcat_files
    from_to = NULL,
    edge_filter = NULL,
    include_self = TRUE,
    overwrite = FALSE,
    # NEW:
    from_col = NULL,
    to_col   = NULL,
    id_col   = NULL,
    up_col   = NULL,
    preset_paths = NULL,       # optional named list: from_to, ftype, translation, special_handling
    nhdplus_mr = NULL,         # deprecated
    ...
) {
  preset <- match.arg(preset, choices = c("none", "mr_nhdplus"))
  if (!is.null(nhdplus_mr)) {
    warning("`nhdplus_mr` is deprecated; use preset='mr_nhdplus' instead.")
    preset <- if (isTRUE(nhdplus_mr)) "mr_nhdplus" else "none"
  }

  if (dir.exists(out_dir)) {
    if (!overwrite) stop("`out_dir` exists. Use overwrite=TRUE to replace: ",
                         out_dir, call. = FALSE)
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (preset == "mr_nhdplus") {
    # Resolve preset paths:
    paths <- NULL

    # 1) explicit paths override
    if (is.list(preset_paths) && !is.null(preset_paths$from_to)) {
      paths <- preset_paths
    }

    # 2) user provided `data_dir` is a list (like your streamcat_files)
    if (is.null(paths) && is.list(data_dir) && !is.null(data_dir$paths)) {
      paths <- data_dir$paths
      if (is.null(paths$from_to) || is.null(paths$ftype) ||
          is.null(paths$translation) || is.null(paths$special_handling)) {
        stop("`data_dir` (list) must include $paths with from_to, ftype, translation, special_handling.")
      }
    }

    # 3) fall back to fetch or a directory
    if (is.null(paths)) {
      if (is.null(data_dir)) {
        info <- sr_fetch_accum_data("mr_nhdplus_conus", quiet = TRUE)
        data_dir <- info$data_dir
      }
      paths <- .sr_accum_expected_paths("mr_nhdplus_conus", data_dir)$paths
    }

    # Allow override of from_to file
    if (!is.null(from_to)) paths$from_to <- from_to

    # Prepare edges using preset defaults
    edges <- .sr_read_edges(paths$from_to, from_col = "FROMCOMID", to_col = "TOCOMID")
    edges <- .sr_prepare_edges_mr_nhdplus(
      edges,
      translation      = paths$translation,
      special_handling = paths$special_handling,
      ftype            = paths$ftype
    )

    if (is.function(edge_filter)) {
      edges <- edge_filter(edges)
    }

    # Names for output columns (COMID / UPCOMIDS unless overridden)
    id_col_final <- id_col %||% "COMID"
    up_col_final <- up_col %||% paste0("UP", id_col_final, "S")

    schema <- list(
      preset = "mr_nhdplus",
      from_col = "FROMCOMID",
      to_col   = "TOCOMID",
      id_col   = id_col_final,
      up_col   = up_col_final,
      gridcode_col   = "GRIDCODE",
      featureid_col  = "FEATUREID",
      translation_path = paths$translation
    )

  } else {
    # preset == "none" (generic)
    if (is.null(from_to)) stop("For preset='none', `from_to` must be provided.", call. = FALSE)

    # Column names (explicit or defaults)
    from_col <- from_col %||% "FROM"
    to_col   <- to_col   %||% "TO"

    edges <- .sr_read_edges(from_to, from_col = from_col, to_col = to_col)

    # Generic cleaning
    drop_from0 <- isTRUE(list(...)[["drop_from0"]] %||% TRUE)
    edges <- .sr_prepare_edges_generic(edges, from_col = from_col,
                                       to_col = to_col, drop_from0 = drop_from0)

    if (is.function(edge_filter)) {
      edges <- edge_filter(edges)
    }

    # Focal ID base name
    id_col_final <- id_col %||% .sr_infer_id_col(from_col, to_col)
    up_col_final <- up_col %||% paste0("UP", id_col_final, "S")

    schema <- list(
      preset = "none",
      from_col = from_col,
      to_col   = to_col,
      id_col   = id_col_final,
      up_col   = up_col_final
    )
  }

  # Build and write upstream pairs with the chosen names
  .sr_build_upstream_pairs_arrow(
    edges = edges,
    out_dir = out_dir,
    from_col = schema$from_col,
    to_col = schema$to_col,
    include_self = include_self,
    id_col_name = schema$id_col,   # NEW
    up_col_name = schema$up_col,   # NEW
    ...
  )

  structure(
    list(
      preset = schema$preset,
      pairs_dir = normalizePath(out_dir, winslash = "/", mustWork = FALSE),
      schema = schema
    ),
    class = "sr_accum_framework"
  )
}

# A tiny helper (place near your other helpers)
`%||%` <- function(x, y) if (!is.null(x)) x else y

# Infer an ID base from FROM*/TO* if id_col is not provided
.sr_infer_id_col <- function(from_col, to_col) {
  # Prefer stripping leading FROM*
  m <- sub("^FROM", "", from_col, ignore.case = TRUE)
  if (nzchar(m) && m != from_col) return(m)
  # Else try TO*
  m2 <- sub("^TO", "", to_col, ignore.case = TRUE)
  if (nzchar(m2) && m2 != to_col) return(m2)
  # Fallback
  "ID"
}
