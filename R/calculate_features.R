#' Calculate and merge features for a list of CF-MS matrices
#' 
#' Score each protein pair in a series of CF-MS experiments, using a measure of
#' association of choice, then merge the calculated features into a single
#' data frame that contains scores for each pair across all of the input
#' experiments. This data frame can then be provided to a classifier as input,
#' alongside a set of 'gold-standard' interacting pairs, in order to score
#' interactions in a manner that integrates data from multiple CF-MS replicates.
#' Note that the final data frame may contain missing values; to replace them, 
#' use the \link{impute_missing_features} function.
#' 
#' @param mats a list of CF-MS matrices, with proteins in rows and fractions in 
#'   columns
#' @param metric the measure of association to use in scoring protein pairs
#' 
#' @return a data frame containing features for all protein pairs across all 
#'   replicates 
#' 
#' @importFrom purrr map
#' @importFrom tidyr drop_na
#' @importFrom reshape2 melt
#' 
#' @export 
calculate_features = function(mats, metric = metrics()) {
  metric = match.arg(metric)
  
  # first, score the pairs
  pairs = map(mats, score_pairs, metric = metric)
  
  # then convert each pair to a data frame
  feature_dfs = map(pairs, ~ melt(.x,
                                  varnames = c('protein_A', 'protein_B'),
                                  value.name = metric,
                                  as.is = TRUE) %>%
                      drop_na() %>%
                      filter(protein_A < protein_B))
  
  # then merge the features
  feature_df = merge_features(feature_dfs)
}

#' Merge features across multiple replicates
#' 
#' Merge features extracted from multiple replicates into a single data 
#' frame for input to a classifier. Note that the merged data frame may contain 
#' missing values; to replace them, use the \link{impute_missing_features}
#' function.
#' 
#' @param feature_dfs a list of feature data frames, each of which has protein
#'   pairs in the first two columns. Alternatively, a list of square matrices 
#'   can be provided as input, and will be coerced into a list of feature
#'   data frames.
#' 
#' @return a data frame containing features for all protein pairs across all 
#'   replicates 
#' 
#' @importFrom purrr map map_lgl
#' @importFrom reshape2 melt
#' @importFrom dplyr full_join filter
#' @importFrom tidyr drop_na
#' @importFrom tester is_square_matrix
#' 
#' @export
merge_features = function(feature_dfs) {
  # catch square matrices and convert them to pairwise data frames
  is_square_mat = map_lgl(feature_dfs, tester::is_square_matrix)
  if (all(is_square_mat)) {
    feature_dfs %<>% map(~ melt(.x,
                                varnames = c('protein_A', 'protein_B'),
                                value.name = 'feature',
                                as.is = TRUE) %>%
                           drop_na() %>%
                           filter(protein_A < protein_B))
  }
  
  features = Reduce(function(x, y) full_join(x, y, by = colnames(x)[c(1, 2)]),
                    feature_dfs)
  return(features)
}

#' Impute missing features with the median, plus or minus some random noise 
#' 
#' Replace missing data within each numeric column of a data frame with 
#' the column median, plus or minus some random noise, in order to train 
#' classifiers that do not readily ignore missing data (e.g. random forests or
#' support vector machines).
#' 
#' @param dat the feature data frame for which to replace missing data
#' @param noise_pct the standard deviation of the random normal 
#'   distribution from which to draw added noise, expressed as a 
#'   percentage of the standard deviation of the non-missing values in each 
#'   column
#' 
#' @return a data frame with missing values in each numeric column replaced
#' by the column median, plus or minus some random noise
#' 
#' @importFrom stats rnorm median sd
#' 
#' @export
impute_missing_features = function(dat, noise_pct = 0.05) {
  for (col_name in colnames(dat)) {
    column = dat[[col_name]]
    if (!is.numeric(column))
      next
    ## first, replace infinite values
    infinite = is.infinite(column)
    column[infinite] = NA
    ## second, replace other missing values
    missing = !is.finite(column)
    dat[[col_name]][missing] = 
      median(column, na.rm = TRUE) + rnorm(sum(missing)) * 
      sd(column, na.rm = TRUE) * noise_pct 
  }
  return(dat)
}
