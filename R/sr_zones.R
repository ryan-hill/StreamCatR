#' Construct sr_zones object (auto-discovers files)
#' @export
sr_zones <- function(
    region_id,
    zone_dir,
    grid_index_path = NULL,
    zidx_path       = NULL,
    wins_path       = NULL,
    blocksize       = NULL
) {
  stopifnot(is.character(region_id), length(region_id) == 1L, nzchar(region_id))
  stopifnot(is.character(zone_dir),  length(zone_dir) == 1L,  nzchar(zone_dir))

  guessed <- .sr_guess_zone_paths(region_id, zone_dir)

  grid_index_path <- grid_index_path %||% guessed$grid_index_path
  zidx_path       <- zidx_path       %||% guessed$zidx_path
  wins_path       <- wins_path       %||% guessed$wins_path
  blocksize       <- blocksize       %||% guessed$blocksize

  if (is.null(grid_index_path) || !file.exists(grid_index_path)) {
    stop("grid index not found for region '", region_id, "' in '", zone_dir, "'. ",
         "Looked for: ", guessed$grid_index_path %||% "<none>", call. = FALSE)
  }
  if (is.null(zidx_path) || !file.exists(zidx_path)) {
    stop("zone index (zidx) not found for region '", region_id, "' in '", zone_dir, "'. ",
         "Looked for: ", guessed$zidx_path %||% "<none>", call. = FALSE)
  }

  # Infer blocksize from windows parquet when available
  if ((is.null(blocksize) || is.na(blocksize) || blocksize <= 0L) &&
      !is.null(wins_path) && file.exists(wins_path)) {
    blocksize <- .sr_blocksize_from_windows(wins_path)
  }

  structure(
    list(
      grid_index_path = grid_index_path,
      zidx_path       = zidx_path,
      wins_path       = wins_path,
      region_id       = region_id,
      blocksize       = blocksize,
      zone_dir        = zone_dir
    ),
    class = "sr_zones"
  )
}
