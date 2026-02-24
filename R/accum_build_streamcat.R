#' Build a watershed-accumulation framework (upstream pairs)
#'
#' Builds a staged upstream-pairs dataset (stored as Arrow parquet files) from a
#' directed from-to topology table. Optionally uses a preset (e.g., MR NHDPlus CONUS)
#' that provides standard column names, translation tables, and default edge cleaning.
#'
#' The resulting framework can be used with [sr_accumulate_ws()] to accumulate zonal
#' statistics upstream across the network.
#'
#' @param out_dir Output directory for the staged upstream-pairs dataset.
#' @param preset Preset configuration. Use `"mr_nhdplus"` for MR NHDPlus conventions,
#'   or `"none"` for a generic from-to table.
#' @param MRNHDPlus Logical; if TRUE, equivalent to `preset="mr_nhdplus"`. If FALSE,
#'   equivalent to `preset="none"`. If NULL, uses `preset`.
#' @param data_dir Directory containing preset parquet inputs (only for presets).
#'   If NULL and a preset is used, data will be fetched on-demand.
#' @param from_to From-to topology table (data.frame/Arrow Table/Dataset) or a file path.
#'   Required if `preset="none"`. Ignored for `preset="mr_nhdplus"` unless explicitly provided.
#' @param edge_filter Optional function that takes a data.frame of edges and returns a filtered
#'   data.frame. Applied after preset default cleaning.
#' @param include_self Logical; if TRUE, include each ID as upstream of itself (common for accumulation).
#' @param overwrite Logical; if TRUE, overwrite existing `out_dir`.
#' @param ... Additional arguments passed to internal builders.
#'
#' @return An object of class `"sr_accum_framework"` containing `pairs_dir`, `preset`,
#'   and schema metadata (including id columns used).
#' @export
sr_streamcat_framework <- function(
    out_dir,
    preset = c("none", "mr_nhdplus"),
    MRNHDPlus = NULL,
    data_dir = NULL,
    from_to = NULL,
    edge_filter = NULL,
    include_self = TRUE,
    overwrite = FALSE,
    ...
) {
  preset <- match.arg(preset, choices = c(
    "none",
    "mr_nhdplus"
  ))
  if (!is.null(MRNHDPlus)) preset <- if (isTRUE(MRNHDPlus)) "mr_nhdplus" else "none"

  if (dir.exists(out_dir)) {
    if (!overwrite) stop("`out_dir` exists. Use overwrite=TRUE to replace: ", out_dir, call. = FALSE)
    unlink(out_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (preset == "mr_nhdplus") {
    # Locate preset inputs (fetch if needed)
    if (is.null(data_dir)) {
      info <- sr_fetch_accum_data("mr_nhdplus_conus", quiet = TRUE)
      data_dir <- info$data_dir
    }
    paths <- .sr_accum_expected_paths("mr_nhdplus_conus", data_dir)$paths

    # Allow override: user-supplied from_to
    if (is.null(from_to)) from_to <- paths$from_to

    # Prepare edges using preset defaults
    edges <- .sr_read_edges(from_to, from_col = "FROMCOMID", to_col = "TOCOMID")
    edges <- .sr_prepare_edges_mr_nhdplus(
      edges,
      translation      = paths$translation,
      special_handling = paths$special_handling,
      ftype            = paths$ftype
    )

    if (is.function(edge_filter)) {
      edges <- edge_filter(edges)
    }

    schema <- list(
      preset = "mr_nhdplus",
      from_col = "FROMCOMID",
      to_col   = "TOCOMID",
      id_col   = "COMID",
      gridcode_col = "GRIDCODE",
      featureid_col = "FEATUREID",
      translation_path = paths$translation
    )
  } else {
    if (is.null(from_to)) stop("For preset='none', `from_to` must be provided.", call. = FALSE)

    # Generic: user must supply (or accept default) column names via ...
    # For now, default to FROM/TO but allow overrides via ...
    dots <- list(...)
    from_col <- dots$from_col %||% "FROM"
    to_col   <- dots$to_col   %||% "TO"

    edges <- .sr_read_edges(from_to, from_col = from_col, to_col = to_col)

    # Generic cleaning: drop missing, drop self loops, drop from==0 if requested
    drop_from0 <- isTRUE(dots$drop_from0 %||% TRUE)
    edges <- .sr_prepare_edges_generic(edges, from_col = from_col, to_col = to_col, drop_from0 = drop_from0)

    if (is.function(edge_filter)) {
      edges <- edge_filter(edges)
    }

    schema <- list(
      preset = "none",
      from_col = from_col,
      to_col   = to_col
    )
  }

  # Build and write upstream pairs dataset
  .sr_build_upstream_pairs_arrow(
    edges = edges,
    out_dir = out_dir,
    from_col = schema$from_col,
    to_col = schema$to_col,
    include_self = include_self,
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
