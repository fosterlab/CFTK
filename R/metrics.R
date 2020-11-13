#' Metrics available in the CFTK
#' 
#' Returns a named vector of measures of association available for scoring 
#' chromatogram similarity in CFTK, where the names of the vector represent
#' human-readable names and the entries themselves are inputs to 
#' \link{score_pairs}.
#'
#' @return a vector enumerating the measures of association implemented in CFTK
#' 
#' @export
metrics = function() {
  c("Bayesian correlation" = "bayes_cor",
    "Biweight midcorrelation" = "bicor",
    "Co-dependency index" = "binomial",
    "Bray-Curtis distance" = "bray_curtis",
    "Canberra distance" = "canberra",
    "Cosine distance" = "cosine",
    "Dice coefficient" = "dice",
    "Distance correlation" = "distance_cor",
    "Euclidean distance" = "euclidean",
    "GENIE3" = "GENIE3",
    "Hamming distance" = "hamming",
    "Jaccard index" = "jaccard",
    "Kendall correlation" = "kendall",
    "Manhattan distance" = "manhattan",
    "Mutual information" = "MI",
    "Pearson correlation" = "pearson",
    "Proportionality (phi)" = "phi_s",
    "Profile correlation" = "profile_cor",
    "Proportionality (rho)" = "rho_p",
    "Spearman correlation" = "spearman",
    "treeClust" = "treeClust",
    "Weighted cross-correlation" = "wccor",
    "Weighted rank correlation" = "weighted_rank",
    "Zero-inflated Kendall correlation" = "zi_kendall"
  )
}
