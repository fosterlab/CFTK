#' Combine multiple CF-MS matrices into a single matrix
#' 
#' Combine fractions from multiple CF-MS matrices into a single matrix, 
#' containing either the union of all proteins found in any individual matrix,
#' or the intersection of proteins found in all matrices.
#' 
#' @param mats a list of matrices
#' @param mode one of \code{'union'} or \code{'intersect'}; the proteins to
#'   include in the merged matrix
#' 
#' @return a merged matrix
#' 
#' @importFrom purrr map
#' @importFrom magrittr %>% %<>% extract
#' 
#' @export
combine_matrices = function(mats, mode = c('union', 'intersect')) {
  mode = match.arg(mode)
  if (mode == 'union') {
    fn = union
  } else {
    fn = intersect
  }
  all_genes = map(mats, rownames) %>% Reduce(union, .)
  mats0 = map(mats, ~ {
    mat1 = .
    missing = setdiff(all_genes, rownames(mat1))
    mat2 = matrix(NA, nrow = length(missing), ncol = ncol(mat1),
                  dimnames = list(missing, colnames(mat1)))
    rbind(mat1, mat2) %>%
      extract(all_genes, )
  })
  mat = do.call(cbind, mats0)
  return(mat)
}
