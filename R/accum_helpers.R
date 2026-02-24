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

    # If NULL, filled from .sr_choose_params(profile=...)
    vertex_batch_size = NULL,
    max_rows_per_part = NULL,
    compression_codec = NULL,

    # Keep explicit defaults (not tuned)
    compression_level = 4L,
    use_dictionary = TRUE,

    # Auto-tuning controls (02-style)
    profile = c("auto", "small", "large"),
    ram_override_gb = NULL,
    cores_override = NULL,

    # Optional mini-batch mode (02-style)
    use_mini_batches = NULL,
    mini_batch_size = NULL
) {
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Need igraph for building upstream pairs.", call. = FALSE)
  }
  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Need arrow for writing parquet.", call. = FALSE)
  }

  profile <- match.arg(profile)

  # ---- choose params (02 behavior) ----
  p <- .sr_choose_params(
    profile         = profile,
    ram_override_gb = ram_override_gb,
    cores_override  = cores_override
  )

  # Fill ONLY when user didn't explicitly set them
  if (is.null(vertex_batch_size)) vertex_batch_size <- p$vertex_batch_size
  if (is.null(max_rows_per_part)) max_rows_per_part <- p$max_rows_per_part
  if (is.null(compression_codec)) compression_codec <- p$compression_codec

  if (is.null(use_mini_batches)) use_mini_batches <- isTRUE(p$use_mini_batches)
  if (is.null(mini_batch_size))  mini_batch_size  <- p$mini_batch_size

  if (isTRUE(use_mini_batches)) {
    if (is.null(mini_batch_size) || !is.finite(mini_batch_size) || mini_batch_size <= 0) {
      stop("`use_mini_batches=TRUE` requires a positive integer `mini_batch_size`.", call. = FALSE)
    }
    mini_batch_size <- as.integer(mini_batch_size)
  }

  vertex_batch_size <- as.integer(vertex_batch_size)
  if (!is.finite(vertex_batch_size) || vertex_batch_size <= 0L) {
    stop("`vertex_batch_size` must be a positive integer.", call. = FALSE)
  }

  # `max_rows_per_part` can be large; keep as numeric but validate
  if (!is.finite(max_rows_per_part) || max_rows_per_part <= 0) {
    stop("`max_rows_per_part` must be a positive number.", call. = FALSE)
  }

  # ---- output dir ----
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ---- normalize/clean edges ----
  e <- edges[, c(from_col, to_col)]
  names(e) <- c("FROM", "TO")

  # Match 02: drop NAs and 0 endpoints (common outlet encoding)
  e <- e[!is.na(e$FROM) & !is.na(e$TO), , drop = FALSE]
  e <- e[e$FROM != 0 & e$TO != 0, , drop = FALSE]

  if (nrow(e) == 0L) {
    stop("No edges remain after filtering NA/0 endpoints. Check inputs.", call. = FALSE)
  }

  # igraph likes character vertex names
  e_from <- as.character(e$FROM)
  e_to   <- as.character(e$TO)

  vids <- sort(unique(c(e_from, e_to)))
  vertices <- data.frame(name = vids, stringsAsFactors = FALSE)

  g <- igraph::graph_from_data_frame(
    d = data.frame(from = e_from, to = e_to, stringsAsFactors = FALSE),
    directed = TRUE,
    vertices = vertices
  )

  # Optional: store numeric COMIDs for compactness; validate conversion
  num_names <- suppressWarnings(as.integer(igraph::V(g)$name))
  if (anyNA(num_names)) {
    stop(
      "Failed converting vertex names to integer COMIDs. ",
      "This indicates non-numeric IDs in the graph.",
      call. = FALSE
    )
  }

  mindist_val <- if (isTRUE(include_self)) 0L else 1L
  order_val   <- igraph::vcount(g)
  n           <- order_val

  # ---- writer ----
  part_counter <- 0L

  write_part <- function(COMIDs, UPCOMIDs) {
    part_counter <<- part_counter + 1L
    df <- data.frame(COMID = COMIDs, UPCOMIDS = UPCOMIDs)
    arrow::write_parquet(
      df,
      sink = file.path(out_dir, sprintf("part-%06d.parquet", part_counter)),
      compression       = compression_codec,
      compression_level = compression_level,
      use_dictionary    = use_dictionary
    )
    invisible(NULL)
  }

  # ---- helper: split an index vector into chunks ----
  .chunk_idxs <- function(x, chunk_size) {
    split(x, ceiling(seq_along(x) / chunk_size))
  }

  # ---- main loop: vertex batches, optional mini-batches ----
  for (start in seq.int(1L, n, by = vertex_batch_size)) {
    idxs <- start:min(start + vertex_batch_size - 1L, n)

    # Further split for low-RAM machines to reduce peak size of ego() result
    idx_groups <- if (isTRUE(use_mini_batches)) {
      .chunk_idxs(idxs, mini_batch_size)
    } else {
      list(idxs)
    }

    for (grp in idx_groups) {
      grp <- as.integer(grp)

      anc_batch <- igraph::ego(
        g,
        order   = order_val,
        nodes   = igraph::V(g)[grp],
        mode    = "in",
        mindist = mindist_val
      )

      rows_per_vertex <- lengths(anc_batch)

      # Partition this ego result into parquet parts by `max_rows_per_part`
      i <- 1L
      m <- length(grp)

      while (i <= m) {
        # Single-vertex overflow guard: prevent infinite loop
        if (rows_per_vertex[i] > max_rows_per_part) {
          # Write this vertex alone (still might be huge; but we at least make progress)
          vs <- anc_batch[[i]]
          sub_idx <- grp[i]

          COMIDs   <- rep(num_names[sub_idx], length(vs))
          UPCOMIDs <- num_names[as.integer(vs)]

          write_part(COMIDs, UPCOMIDs)
          rm(COMIDs, UPCOMIDs); gc(FALSE)

          i <- i + 1L
          next
        }

        rows_acc  <- 0
        sub_start <- i

        while (i <= m && (rows_acc + rows_per_vertex[i]) <= max_rows_per_part) {
          rows_acc <- rows_acc + rows_per_vertex[i]
          i <- i + 1L
        }

        sub_end <- i - 1L
        if (sub_end < sub_start) next

        sub_anc  <- anc_batch[sub_start:sub_end]
        sub_idxs <- grp[sub_start:sub_end]
        sub_rows <- rows_per_vertex[sub_start:sub_end]

        COMIDs <- rep(num_names[sub_idxs], sub_rows)

        # UPCOMIDs: concatenate ancestors per vertex in the same order
        UPCOMIDs <- unlist(
          lapply(sub_anc, function(vs) num_names[as.integer(vs)]),
          use.names = FALSE
        )

        write_part(COMIDs, UPCOMIDs)
        rm(COMIDs, UPCOMIDs); gc(FALSE)
      }

      rm(anc_batch, rows_per_vertex); gc(FALSE)
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

#' @keywords internal
#' @noRd
.sr_choose_params <- function(
    profile = c("auto", "small", "large"),
    target_frac_ram = NULL,
    default_ram_gb = 8,
    row_bytes_est = 24L,
    ram_override_gb = NULL,
    cores_override  = NULL
) {
  profile <- match.arg(profile)

  avail_gb <- tryCatch(ps::ps_virtual_memory()$available / 2^30, error = function(e) NA_real_)
  cores    <- tryCatch(future::availableCores(), error = function(e) NA_integer_)

  if (!is.null(ram_override_gb)) avail_gb <- ram_override_gb
  if (!is.null(cores_override))  cores    <- cores_override
  if (is.na(avail_gb)) avail_gb <- default_ram_gb
  if (is.na(cores)) cores <- 4L

  if (is.null(target_frac_ram)) {
    target_frac_ram <- switch(
      profile,
      small = 0.10,
      large = 0.30,
      auto  = if (avail_gb >= 64) 0.25 else if (avail_gb >= 32) 0.20 else 0.15
    )
  }

  budget_gb         <- max(0.5, avail_gb * target_frac_ram)
  budget_bytes      <- budget_gb * 2^30
  max_rows_per_part <- as.integer(floor(budget_bytes / row_bytes_est))

  if (profile == "small") {
    vertex_batch_size <- 5000L
    compression_codec <- if (cores <= 4) "snappy" else "zstd"
    use_mini_batches  <- TRUE
    mini_batch_size   <- 2000L
  } else if (profile == "large") {
    vertex_batch_size <- 150000L
    compression_codec <- "zstd"
    use_mini_batches  <- FALSE
    mini_batch_size   <- NA_integer_
  } else {
    vertex_batch_size <- if (avail_gb < 16) 10000L else if (avail_gb < 64) 50000L else 100000L
    compression_codec <- if (cores <= 4 && avail_gb < 16) "snappy" else "zstd"
    use_mini_batches  <- avail_gb < 16
    mini_batch_size   <- if (use_mini_batches) 2000L else NA_integer_
  }

  list(
    avail_gb          = avail_gb,
    cores             = cores,
    budget_gb         = budget_gb,
    vertex_batch_size = vertex_batch_size,
    max_rows_per_part = max_rows_per_part,
    compression_codec = compression_codec,
    compression_level = 4L,
    use_mini_batches  = use_mini_batches,
    mini_batch_size   = mini_batch_size
  )
}
