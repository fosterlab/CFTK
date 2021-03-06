% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/to_adjacency_matrix.R
\name{to_adjacency_matrix}
\alias{to_adjacency_matrix}
\title{Create an adjacency matrix from a pairwise data frame}
\usage{
to_adjacency_matrix(dat, symmetric = TRUE, node_columns = c(1, 2))
}
\arguments{
\item{dat}{a data frame containing pairwise interactions}

\item{symmetric}{if true, interactions in both directions will be 
added to the adjacency matrix}

\item{node_columns}{a vector of length two, denoting either the indices 
(integer vector) or column names (character vector) of the columns within 
the data frame containing the nodes participating in pairwise interactions;
defaults to the first two columns of the data frame (\code{c(1, 2)})}
}
\value{
an adjacency matrix between all interacting proteins
}
\description{
Convert a data frame containing pairwise interactions into an
adjacency matrix. The resulting square adjacency matrix contains ones for 
interacting proteins and zeroes for non-interacting proteins.
}
\examples{
ppi = data.frame(protein_A = paste0("protein", seq_len(10)),
                 protein_B = paste0("protein", c(rep(3, 2), rep(5, 5), 
                                    rep(7, 3))), stringsAsFactors = FALSE)
adj = to_adjacency_matrix(ppi)
}
