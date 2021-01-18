#' Threshold interactions at a given precision cutoff
#'
#' @param interactions the ranked list of interactions output by
#'   \code{\link{score_interactions}}, including a \code{precision}
#'   column as calculated after the fact by \code{\link{calculate_precision}}
#' @param threshold the minimum precision of the unweighted interaction
#'   network to return
#'
#' @return the subset of the original ranked list of interactions including 
#'   the greatest number of interactions that yields a network at the given
#'   precision
#'
#' @export
threshold_precision = function(interactions, threshold) {
  up_to = which(interactions$precision >= threshold)
  if (!length(up_to)) {
    stop("no interactions at precision >= ", threshold)
  } else {
    interactions[seq_len(max(up_to, na.rm = T)), ]
  }
}