#' Filter low-quality chromatograms from a CF-MS chromatogram matrix
#' 
#' Remove chromatograms from a matrix if the protein in question was not 
#' quantified in some minimum number of fractions. 
#' 
#' @param mat a CF-MS chromatogram matrix, with proteins in rows and fractions 
#'   in columns
#' @param min_fractions minimum number of fractions for a 'good-quality' 
#'   chromatogram
#' 
#' @return a filtered matrix with low-quality chromatograms removed
#' 
#' @export
filter_matrix = function(mat, min_fractions = 0) {
  # sanitize missing values, in case this hasn't been done yet
  mat = sanitize_missing_values(mat)
  # count number of fractions
  n_fractions = rowSums(!is.na(mat))
  # discard chromatograms with less than that number
  mat = mat[n_fractions >= min_fractions, ]
  return(mat)
}
