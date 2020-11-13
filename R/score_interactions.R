#' Score potential interactions in cross-validation
#' 
#' Use a machine-learning approach to integrate data from across multiple CF-MS
#' replicates, or multiple features within a single replicate. This function
#' takes a data frame of features as input, alongside a set of 'gold-standard'
#' reference interactions. The gold standard is split into a user-specified
#' number of folds, and a classifier is trained on the reference interactions
#' after leaving out each fold in turn. Each classifier is then applied to
#' predict interactions in the entire feature data frame, minus the protein
#' pairs that overlap with the training interactions for that fold. The mean
#' classifier score across all folds is calculated for each protein pair,
#' and the proteins are sorted in descending order by their mean score.
#' 
#' @param features a data frame containing features for all protein pairs 
#'   across all replicates, containing columns \code{protein_A} and 
#'   \code{protein_B}, as returned by \link{calculate_features} 
#' @param gold_standard a data frame with columns \code{protein_A} and 
#'   \code{protein_B}, as returned by \link{to_pairwise_df}, containing 
#'   'gold standard' interacting protein pairs
#' @param classifier the classifier to use; one of \code{'RF'} (random forest),
#'   \code{'NB'} (naive Bayes), \code{'SVM'} (support vector machine), or
#'   \code{'LR'} (logistic regression)
#' @param split_by the mechanism by which to split the gold standard into
#'   cross-validation folds; either by protein complex subunits 
#'   (\code{'proteins'}) or by pairwise interactions between those subunits 
#'   (\code{'pairs'})
#' @param n_folds the number of folds of cross-validation to perform
#' @param verbose set to \code{FALSE} to disable messages from the function
#' 
#' @import naivebayes
#' @import LiblineaR
#' @import ranger
#' @import speedglm
#' @importFrom dplyr left_join pull select mutate arrange desc
#' @importFrom magrittr %>%
#' @importFrom stats predict setNames
#' 
#' @export
score_interactions = function(features, gold_standard,
                              classifier = c('RF', 'NB', 'SVM', 'LR'),
                              split_by = c('proteins', 'pairs'),
                              n_folds = 10,
                              verbose = TRUE) {
  classifier = match.arg(classifier)
  split_by = match.arg(split_by)
  
  # get splits
  splits = get_splits(gold_standard, split_by = split_by, n_folds = n_folds)
  
  # set up CV scores holder
  n_interactions = nrow(features)
  interaction_names = with(features, paste0(protein_A, '|', protein_B))
  clf_scores = matrix(NA, ncol = n_folds, nrow = n_interactions,
                      dimnames = list(interaction_names))
  
  # iterate through splits
  for (split_idx in seq_along(splits)) {
    if (verbose)
      message("working on split ", split_idx, " of ", length(splits), " ...")
    split = splits[[split_idx]]
    train_labels = split$train
    test_labels = split$test
    
    # make labels
    train_labels = features %>%
      left_join(train_labels, by = c('protein_A', 'protein_B')) %>% 
      pull(label)
    test_labels = features %>%
      left_join(test_labels, by = c('protein_A', 'protein_B')) %>% 
      pull(label)
    labels = list(train_label = train_labels, test_label = test_labels)
    
    # train classifier on labelled training data
    clf_idxs = which(!is.na(labels$train_label))
    clf_labels = labels$train_label[clf_idxs] %>%
      as.factor()
    clf_data = features[clf_idxs, ] %>%
      dplyr::select(-protein_A, -protein_B)
    clf_data_labeled = clf_data %>%
      mutate(label = clf_labels)
    clf = switch(classifier,
                 NB = naive_bayes(clf_data, clf_labels),
                 SVM = LiblineaR(clf_data, clf_labels, type = 2),
                 RF = ranger(data = clf_data_labeled, 
                             dependent.variable.name = "label",
                             probability = TRUE,
                             num.trees = 100,
                             num.threads = n_threads),
                 LR = speedglm(label ~ ., clf_data_labeled, 
                               family = stats::binomial()))
    
    # predict interactions on the withheld set
    test_data = features[-clf_idxs, ] %>%
      dplyr::select(-protein_A, -protein_B)
    predictions = switch(
      classifier, 
      NB = predict(clf, test_data, type = 'prob', threshold = 5e-324),
      SVM = predict(clf, test_data, decisionValues = TRUE),
      RF = predict(clf, test_data, num.threads = n_threads),
      LR = predict(clf, test_data, type = 'response'))
    
    ## extract predictions as numeric vector
    predictions = switch(
      classifier,
      NB = predictions[, "1"],
      SVM = -1.0 * predictions$decisionValues[, "0"],
      RF = predictions[[1]][, "1"],
      LR = predictions)
    # assign scores
    clf_scores[-clf_idxs, split_idx] = predictions
  }
  
  # average over folds 
  mean_scores = setNames(rowMeans(clf_scores, na.rm = TRUE),
                         rownames(clf_scores))
  
  # collate the labels from the last train/test split
  merged_labels = pmin(labels$train_label, labels$test_label, na.rm = TRUE)
  
  # create ranked data frame
  interactions = features %>%
    dplyr::select(protein_A, protein_B) %>%
    mutate(score = mean_scores, label = merged_labels) %>%
    arrange(desc(score)) 
  
  return(interactions)
}

#' Split a set of reference protein-protein interactions into cross-validation
#' folds, for use in \link{score_interactions}. 
#' 
#' @param gold_standard a data frame with columns \code{protein_A} and 
#'   \code{protein_B}, as returned by \link{to_pairwise_df}, containing 
#'   'gold standard' interacting protein pairs
#' @param n_folds the number of folds of cross-validation to perform
#' @param split_by the mechanism by which to split the gold standard into
#'   cross-validation folds; either by protein complex subunits 
#'   (\code{'proteins'}) or by pairwise interactions between those subunits 
#'   (\code{'pairs'})
#'   
#' @return a list of length \code{n_folds}, in which each item contains 
#'   \code{"train"} and \code{"test"} data frames with the columns 
#'   \code{protein_A}, \code{protein_B}, and \code{label}
#' 
#' @importFrom purrr map
#' @importFrom dplyr bind_rows distinct
#' @importFrom magrittr extract %>% %<>%
#' @importFrom reshape2 melt
#' @importFrom stats setNames
#' 
#' @export
get_splits = function(gold_standard,
                      n_folds = 10,
                      split_by = c('proteins', 'pairs')) {
  split_by = match.arg(split_by)
  percent_in_train = (n_folds - 1) / n_folds
  
  # convert to adjacency matrix
  adj = gold_standard %>%
    distinct(protein_A, protein_B) %>%
    to_adjacency_matrix()
  
  # split by gold standard into folds
  if (split_by == 'proteins') {
    # split the matrix into folds
    borders = floor(seq(1, nrow(adj), nrow(adj) / n_folds))
    fold_idxs = seq_len(nrow(adj)) %>%
      sample() %>%
      split(borders)
    folds = map(fold_idxs, ~ extract(adj, ., .)) %>%
      setNames(seq_along(.))
    
    # convert back to long
    folds %<>%
      map(~ reshape2::melt(., varnames = c('protein_A', 'protein_B'), 
                           value.name = 'label', as.is = TRUE) %>% 
            filter(protein_A < protein_B))
  } else if (split_by == 'pairs') {
    # convert back to long
    pairs = adj %>%
      melt(varnames = c('protein_A', 'protein_B'), 
           value.name = 'label', as.is = TRUE) %>% 
      filter(protein_A < protein_B)
    
    # now, split the pairs into folds
    borders = floor(seq(1, nrow(pairs), nrow(pairs) / n_folds))
    fold_idxs = seq_len(nrow(pairs)) %>%
      sample() %>%
      split(borders)
    folds = map(fold_idxs, ~ extract(pairs, ., )) %>%
      setNames(seq_along(.))
  } else {
    stop("not sure what to do with split_by: '", split_by, "'")
  }
  
  # get labels for each fold
  splits = list()
  for (test_idx in seq_along(folds)) {
    test = folds[[test_idx]]
    train = folds[-test_idx] %>% bind_rows()
    # add to split
    split = list(train = train, test = test)
    splits[[test_idx]] = split
  }
  splits %<>% setNames(seq_along(.))
  
  return(splits)
}
