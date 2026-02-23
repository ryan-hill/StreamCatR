#' @importFrom rlang .data
#' @importFrom stats setNames
NULL

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
.sr_build_upstream_pairs_arrow <- function(
    edges,
    out_dir,
    from_col,
    to_col,
    include_self = TRUE,
    vertex_batch_size = 5000L,
    max_rows_per_part = 2e7,
    compression_codec = "zstd",
    compression_level = 4L,
    use_dictionary = TRUE
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Need igraph for building upstream pairs.", call. = FALSE)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Build graph
  e <- edges[, c(from_col, to_col)]
  names(e) <- c("FROM", "TO")

  # Vertex set: all IDs seen (excluding 0)
  vids <- sort(unique(c(e$FROM, e$TO)))
  vids <- vids[!is.na(vids) & vids != 0]
  vertices <- data.frame(name = as.character(vids), stringsAsFactors = FALSE)

  g <- igraph::graph_from_data_frame(
    d = data.frame(from = as.character(e$FROM), to = as.character(e$TO)),
    directed = TRUE,
    vertices = vertices
  )

  num_names <- as.integer(igraph::V(g)$name)
  if (any(is.na(num_names))) stop("NA while converting vertex names to int.", call. = FALSE)

  mindist_val <- if (isTRUE(include_self)) 0L else 1L
  order_val   <- igraph::vcount(g)

  part_counter <- 0L

  write_part <- function(COMIDs, UPCOMIDs) {
    part_counter <<- part_counter + 1L
    df <- data.frame(COMID = COMIDs, UPCOMIDS = UPCOMIDs)
    arrow::write_parquet(
      df,
      sink = file.path(out_dir, sprintf("part-%06d.parquet", part_counter)),
      compression = compression_codec,
      compression_level = compression_level,
      use_dictionary = use_dictionary
    )
  }

  n <- igraph::vcount(g)

  for (start in seq(1L, n, by = vertex_batch_size)) {
    idxs <- start:min(start + vertex_batch_size - 1L, n)

    anc_batch <- igraph::ego(
      g, order = order_val, nodes = igraph::V(g)[idxs],
      mode = "in", mindist = mindist_val
    )
    rows_per_vertex <- lengths(anc_batch)

    i <- 1L
    while (i <= length(idxs)) {
      rows_acc <- 0L
      sub_start <- i
      while (i <= length(idxs) && rows_acc + rows_per_vertex[i] <= max_rows_per_part) {
        rows_acc <- rows_acc + rows_per_vertex[i]
        i <- i + 1L
      }
      sub_end <- i - 1L
      if (sub_end < sub_start) next

      sub_anc  <- anc_batch[sub_start:sub_end]
      sub_idxs <- idxs[sub_start:sub_end]
      sub_rows <- rows_per_vertex[sub_start:sub_end]

      COMIDs   <- rep(num_names[sub_idxs], sub_rows)
      UPCOMIDs <- unlist(lapply(sub_anc, function(vs) num_names[as.integer(vs)]), use.names = FALSE)

      write_part(COMIDs, UPCOMIDs)
      rm(COMIDs, UPCOMIDs)
      gc(FALSE)
    }
  }

  invisible(out_dir)
}

#' @keywords internal
#' @noRd
.sr_prepare_edges_mr_nhdplus <- function(edges, special_handling, ftype, translation) {
  edges <- .sr_prepare_edges_generic(edges, "FROMCOMID", "TOCOMID", drop_from0 = TRUE)

  # Coastline COMIDs (FTYPE == COASTLINE)
  ft <- .sr_read_table(ftype)
  names(ft) <- toupper(names(ft))
  coast_ids <- integer(0)
  if (all(c("COMID","FTYPE") %in% names(ft))) {
    coast_ids <- unique(as.integer(ft$COMID[toupper(as.character(ft$FTYPE)) == "COASTLINE"]))
  }

  # Sink COMIDs: negative FEATUREID
  tr <- .sr_read_table(translation)
  names(tr) <- toupper(names(tr))
  sink_ids <- integer(0)
  if ("FEATUREID" %in% names(tr)) {
    vals <- as.integer(tr$FEATUREID)
    sink_ids <- unique(vals[!is.na(vals) & vals < 0L & vals != 0L])
  }

  # Special removals
  sh <- .sr_read_table(special_handling)
  names(sh) <- toupper(names(sh))
  special_remove <- integer(0)
  if ("REMOVEFROMCOMID" %in% names(sh)) {
    special_remove <- unique(as.integer(sh$REMOVEFROMCOMID[!is.na(sh$REMOVEFROMCOMID) & sh$REMOVEFROMCOMID != 0L]))
  }

  remove_from <- unique(c(coast_ids, special_remove))

  # Apply MR rules
  edges <- edges[!edges$FROMCOMID %in% remove_from, , drop = FALSE]
  edges <- edges[!edges$TOCOMID %in% sink_ids, , drop = FALSE]

  # (Optional) return coast_ids/sink_ids as attributes for debugging
  attr(edges, "coast_ids") <- coast_ids
  attr(edges, "sink_ids") <- sink_ids
  attr(edges, "special_remove") <- special_remove
  edges
}

#' @keywords internal
#' @noRd
.sr_accumulate_with_pairs <- function(cat_tbl, pairs, id_col, sum_col, count_col) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Need dplyr.", call. = FALSE)
  if (!requireNamespace("arrow", quietly = TRUE)) stop("Need arrow.", call. = FALSE)

  # Ensure we have a data.frame for catchment metrics
  cat_tbl <- as.data.frame(cat_tbl)

  # Standardize names to stable outputs
  names(cat_tbl)[names(cat_tbl) == id_col] <- "COMID"
  names(cat_tbl)[names(cat_tbl) == sum_col] <- "cat_sum"
  names(cat_tbl)[names(cat_tbl) == count_col] <- "cat_count"

  # Build watershed sums by joining UPCOMIDS -> COMID on catchment table
  ds_pairs <- pairs

  ws <- ds_pairs |>
    dplyr::transmute(
      COMID    = as.integer(.data$COMID),
      UPCOMIDS = as.integer(.data$UPCOMIDS)
    ) |>
    dplyr::left_join(
      arrow::Table$create(cat_tbl),
      by = c("UPCOMIDS" = "COMID")
    ) |>
    dplyr::group_by(.data$COMID) |>
    dplyr::summarise(
      ws_sum   = sum(.data$cat_sum, na.rm = TRUE),
      ws_count = sum(.data$cat_count, na.rm = TRUE),
      .groups = "drop"
    )

  ws <- as.data.frame(ws)

  # Join cat + ws
  out <- dplyr::left_join(
    cat_tbl,
    ws,
    by = "COMID"
  )

  data.table::as.data.table(out)
}
