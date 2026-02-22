#' Construct sr_zones object
#'
#' @keywords internal
#' @noRd
sr_zones <- function(
    grid_index_path,
    zidx_path,
    wins_path = NULL,
    region_id = NULL,
    blocksize = NULL,
    zone_dir  = NULL
) {
  stopifnot(
    is.character(grid_index_path), length(grid_index_path) == 1L,
    is.character(zidx_path),       length(zidx_path) == 1L
  )

  if (!is.null(wins_path)) stopifnot(is.character(wins_path), length(wins_path) == 1L)
  if (!is.null(region_id)) stopifnot(is.character(region_id), length(region_id) == 1L)
  if (!is.null(blocksize)) blocksize <- as.integer(blocksize)
  if (!is.null(zone_dir))  stopifnot(is.character(zone_dir), length(zone_dir) == 1L)

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

#' @export
print.sr_zones <- function(x, ...) {
  cat("<sr_zones>\n")
  if (!is.null(x$region_id)) cat("  Region:     ", x$region_id, "\n")
  if (!is.null(x$blocksize)) cat("  Blocksize:  ", x$blocksize, "\n")
  if (!is.null(x$zone_dir))  cat("  Zone dir:   ", x$zone_dir, "\n")
  cat("  Zidx:       ", x$zidx_path, "\n")
  cat("  Grid index: ", x$grid_index_path, "\n")
  cat("  Windows:    ", if (is.null(x$wins_path)) "<not built>" else x$wins_path, "\n")
  invisible(x)
}
