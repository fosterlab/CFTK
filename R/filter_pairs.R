#' Filter protein pairs from a chromatogram score matrix
#' 
#' Replace pairs of proteins in a matrix output by \link{score_pairs} with 
#' \code{NA}s, if the two proteins were not jointly quantified in some minimum
#' number of overlapping fractions.
#' 
#' @param pairs a square matrix scoring the similarity of each pair of 
#'   chromatograms, as returned by \link{score_pairs}
#' @param mat The original CF-MS chromatogram matrix, with proteins in rows and 
#'   fractions in columns
#' @param min_pairs minimum number of fractions in which a pair of proteins must 
#'   both be detected to consider scoring the candidate interaction
#' 
#' @return a filtered version of the input square matrix, with a subset of 
#'    protein pairs replaced by \code{NA}s
#' 
#' @export
filter_pairs = function(pairs, mat, min_pairs = 0) {
  # check dimensions
  if (nrow(pairs) != ncol(pairs))
    stop("pairs matrix must be square")
  if (nrow(pairs) != nrow(mat))
    stop("dimensions of pairs and chromatogram matrices do not match")
  
  # sanitize missing values, in case this hasn't been done yet
  mat = sanitize_missing_values(mat)
  # count the number of missing values
  n_pairs = crossprod(!is.na(t(mat)) & t(mat) != 0)
  # replace pairs with <min_n with NA
  pairs[n_pairs <= min_pairs] = NA
  return(pairs)
}
