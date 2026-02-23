#' @keywords internal
#' @noRd
.sr_accum_cache_dir <- function(dataset) {
  # Prefer R's user cache dir if available
  if (requireNamespace("tools", quietly = TRUE) && "R_user_dir" %in% getNamespaceExports("tools")) {
    base <- tools::R_user_dir("StreamCatR", which = "cache")
  } else if (requireNamespace("rappdirs", quietly = TRUE)) {
    base <- rappdirs::user_cache_dir("StreamCatR")
  } else {
    base <- file.path(tempdir(), "StreamCatR-cache")
  }
  file.path(base, "accum", dataset)
}

#' @keywords internal
#' @noRd
.sr_accum_expected_paths <- function(dataset, data_dir) {
  if (dataset != "mr_nhdplus_conus") stop("Unknown dataset: ", dataset, call. = FALSE)
  paths <- list(
    from_to          = file.path(data_dir, "full_from_to_conus.parquet"),
    ftype            = file.path(data_dir, "ftype_conus.parquet"),
    translation      = file.path(data_dir, "gridcode_comid_translation_conus.parquet"),
    special_handling = file.path(data_dir, "special_comid_handling.parquet")
  )
  list(dataset = dataset, data_dir = data_dir, paths = paths)
}

#' @keywords internal
#' @noRd
.sr_read_table <- function(x) {
  if (is.character(x) && length(x) == 1L) {
    if (!file.exists(x)) stop("Missing file: ", x, call. = FALSE)
    return(arrow::read_parquet(x, as_data_frame = TRUE))
  }
  if (inherits(x, "ArrowTabular")) {
    return(as.data.frame(x))
  }
  if (is.data.frame(x)) return(x)
  stop("Unsupported table input type.", call. = FALSE)
}

#' @keywords internal
#' @noRd
.sr_read_edges <- function(from_to, from_col, to_col) {
  df <- .sr_read_table(from_to)
  if (!all(c(from_col, to_col) %in% names(df))) {
    stop("from_to must contain columns: ", from_col, ", ", to_col, call. = FALSE)
  }
  df
}

#' @keywords internal
#' @noRd
.sr_prepare_edges_generic <- function(edges, from_col, to_col, drop_from0 = TRUE) {
  edges <- edges[, c(from_col, to_col)]
  edges <- edges[!is.na(edges[[from_col]]) & !is.na(edges[[to_col]]), , drop = FALSE]
  if (drop_from0) edges <- edges[edges[[from_col]] != 0, , drop = FALSE]
  edges <- edges[edges[[from_col]] != edges[[to_col]], , drop = FALSE]
  edges
}

#' @keywords internal
#' @noRd
.sr_prepare_edges_mr_nhdplus <- function(edges, special_handling, ftype) {
  edges <- .sr_prepare_edges_generic(edges, "FROMCOMID", "TOCOMID", drop_from0 = TRUE)

  # Drop special FROMCOMID edges
  sh <- .sr_read_table(special_handling)
  if ("removeFROMCOMID" %in% names(sh)) {
    edges <- edges[!edges$FROMCOMID %in% sh$removeFROMCOMID, , drop = FALSE]
  }

  edges
}
