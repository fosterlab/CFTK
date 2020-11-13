#' Log-transform or quantile normalize a CF-MS matrix
#' 
#' Transform a CF-MS matrix by either log-transformation or quantile 
#' normalization. Note that quantile normalization relies on the implementation
#' in the \code{preprocessCore} Bioconductor package, which is not installed 
#' alongside the CFTK pacakge, and must be installed from Bioconductor 
#' separately by the user.
#' 
#' @param mat a CF-MS chromatogram matrix, with proteins in rows and fractions 
#'   in columns
#' @param mode the transformation to apply; one of \code{'quantile'} or 
#'   \code{'log-transform'}
#' 
#' @export
transform_matrix = function(mat, mode = c('quantile', 'log-transform')) {
  mode = match.arg(mode)
  if (mode == 'quantile') {
    dims = dimnames(mat)
    mat %<>% preprocessCore::normalize.quantiles()
    dimnames(mat) = dims
  } else if (mode == 'log-transform') {
    mat %<>% log()
  } else {
    stop("not sure what to do with mode '", mode, "'")
  }
  return(mat)
}