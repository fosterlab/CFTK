#' Replace missing values in a CF-MS chromatogram matrix
#' 
#' Sanitize different representations of missing values (e.g., zeroes, 
#' \code{NaN}s, or infinite values) in a chromatogram matrix, then optionally
#' replace them all with zeroes or with near-zero noise. 
#'
#' @param mat a CF-MS chromatogram matrix, with proteins in rows and fractions 
#'   in columns
#' @param missing one of \code{'NA'}, \code{'zero'}, or \code{'noise'}. 
#'   If \code{'zero'}, missing values will be replaced with zeroes.
#'   If \code{'noise'}, missing values will be replaced with random noise, 
#'   sampled from a uniform distribution with a maximum value of 
#'   \code{noise_floor}. If \code{'NA'}, missing values will be left as 
#'   \code{NA}s.
#' @param noise_floor upper limit of the uniform distribution for 
#'   \code{missing = 'noise'}
#' 
#' @return a copy of the input matrix with missing values optionally replaced
#' 
#' @examples 
#' mat = matrix(rnorm(100), ncol = 10, nrow = 10,
#'              dimnames = list(paste0('protein_', LETTERS[1:10]),
#'                              paste0('fraction_', seq_len(10))))
#' mat[sample(seq_along(mat), 5)] = 0
#' mat[sample(seq_along(mat), 8)] = NA
#' mat[sample(seq_along(mat), 2)] = Inf
#' sanitized = replace_missing_values(mat, missing = 'zero')
#' table(sanitized == 0)
#' 
#' @importFrom stats runif
#' 
#' @export
replace_missing_values = function(mat, missing = c('NA', 'zero', 'noise'),
                                 noise_floor = 0.001) {
  missing = match.arg(missing)
  
  # first, standardize the missing values
  mat = sanitize_missing_values(mat)
  
  # next, handle them
  if (missing == 'zero') {
    mat[is.na(mat)] = 0
  } else if (missing == 'noise') {
    is_missing = is.na(mat)
    mat[is_missing] = runif(sum(is_missing), min = 0, max = noise_floor)
  } else {
    # do nothing
  }
  
  return(mat)
}

#' Sanitize missing values in a CF-MS chromatogram matrix
#' 
#' Standardize the way missing values are represented in a chromatogram matrix
#' by replacing zeroes, infinite values, and \code{NaN}s with \code{NA}. 
#' 
#' @param mat a CF-MS chromatogram matrix, with proteins in rows and fractions 
#'   in columns
#' 
#' @return a sanitized matrix with all zeroes, infinite values, and \code{NaN}s
#'   replaced with \code{NA}
#' 
#' @export
sanitize_missing_values = function(mat) {
  mat[!is.finite(mat)] = NA
  mat[is.nan(mat)] = NA
  mat[mat == 0] = NA
  return(mat)
}
