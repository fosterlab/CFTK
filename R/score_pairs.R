#' Score the similarity of all pairs of chromatograms
#' 
#' Compute an index measuring the similarity of chromatogram profiles for every 
#' pair of proteins detected in a CF-MS experiment. This function implements
#' a total of 24 such indices, or measures of association. The full list 
#' can be accessed using the \link{metrics} function. 
#'
#' @param mat a CF-MS chromatogram matrix, with proteins in rows and fractions 
#'   in columns
#' @param metric the measure of association to use in scoring protein pairs
#' @param seed seed for the random number generator, used to ensure the 
#'   reproducibility some metrics, such as GENIE3
#' 
#' @return a matrix of dimensions (# of proteins) x (# of proteins), scoring
#'   every possible protein pair, in which higher values reflect more similar
#'   pairs
#' 
#' @importFrom purrr map2_dbl
#' @importFrom stats cor
#' @importFrom WGCNA bicor mutualInfoAdjacency
#' @importFrom lsa cosine
#' @importFrom propr phis perb
#' @importFrom treeClust treeClust treeClust.control
#' @importFrom GENIE3 GENIE3
#' @importFrom Pigengene dcor.matrix
#' @importFrom gtools permutations
#' @importFrom vegan vegdist
#' @importFrom wccsom wcc
#' 
#' @export
score_pairs = function(mat, metric = metrics(), seed = 0) {
  # set seed for reproducibilty of some metrics, e.g., GENIE3
  set.seed(seed)
  
  # transpose the matrix
  mat = t(mat)
  
  if (metric %in% c("pearson", "spearman", "kendall")) {
    cor = cor(mat, method = metric, use = 'pairwise.complete.obs')
  } else if (metric == 'bicor') {
    cor = bicor(mat, use = 'pairwise.complete.obs')
  } else if (metric == 'binomial') {
    cor = binomial(mat)
  } else if (metric == 'cosine') {
    cor = cosine(mat)
  } else if (metric == 'jaccard') {
    cor = jaccard(mat)
  } else if (metric == 'canberra') {
    cor = -1.0 * as.matrix(dist(t(mat), method = 'canberra'))
  } else if (metric == 'euclidean') {
    cor = -1.0 * as.matrix(dist(t(mat), method = 'euclidean'))
  } else if (metric == 'manhattan') {
    cor = -1.0 * as.matrix(dist(t(mat), method = 'manhattan'))
  } else if (metric == 'weighted_rank') {
    cor = wtd_rank(mat)
  } else if (metric == 'hamming') {
    cor = -1.0 * hamming(mat)
  } else if (metric == 'dice') {
    cor = -1.0 * dice(mat)
  } else if (metric == 'phi_s') {
    cor = -1.0 * phis(mat, select = colnames(mat))@matrix
  } else if (metric == 'rho_p') {
    cor = perb(mat, select = colnames(mat))@matrix
  } else if (metric == 'zi_kendall') {
    cor = kendall_zi(mat)
  } else if (metric == 'MI') {
    cor = mutualInfoAdjacency(mat)$AdjacencyUniversalVersion1
    rownames(cor) = colnames(cor) = colnames(mat)
  } else if (metric == 'bayes_cor') {
    cor = Bayes_Corr_Prior3(t(mat))
  } else if (metric == 'bray_curtis') {
    cor = -1.0 * as.matrix(vegdist(t(mat), method = 'bray', na.rm = T))
  } else if (metric == 'wccor') {
    permutations = permutations(n = ncol(mat), r = 2, v = seq_len(ncol(mat)))
    cor = matrix(NA, nrow = ncol(mat), ncol = ncol(mat))
    cor[permutations] = map2_dbl(permutations[, 1], permutations[, 2],
                                 ~ wcc(mat[, .x], mat[, .y], trwdth = 1))
    rownames(cor) = colnames(cor) = colnames(mat)
  } else if (metric == 'distance_cor') {
    cor = dcor.matrix(mat)
  } else if (metric == 'profile_cor') {
    cor1 = cor(mat, method = 'p', use = 'pairwise.complete.obs')
    cor = cor(cor1, use = 'pairwise.complete.obs')
  } else if (metric == 'GENIE3') {
    cor = GENIE3(t(mat))
  } else if (metric == 'treeClust') {
    input = mat %>%
      t() %>%
      as.data.frame() %>%
      set_colnames(paste0('fraction', seq_len(ncol(.))))
    dist = treeClust(dfx = input, d.num = 2, 
                     control = treeClust.control(return.dists = T))$dists
    cor = -1.0 * as.matrix(dist)
  } else {
    stop("don't know what to do for metric: ", metric)
  }
  
  return(cor)
}
