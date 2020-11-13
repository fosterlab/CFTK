#' Convert a data frame to a list of gene sets
#' 
#' Convert a data frame with gene sets and gene names in separate columns
#' into a list of gene sets, in which the names of the list are the gene sets
#' and the entries are vectors of genes.
#'
#' @param annotations a data frame containing some annotations
#' @param set_col the name of the column that contains the gene sets
#'   (e.g. Gene Ontology terms)
#' @param gene_col the name of the column that contains the annotated targets 
#' (e.g. UniProt accessions)
#' @return a named list where the names are found in keyCol and the entries are 
#' unique items from termCol
#' 
#' @importFrom utils unstack
#' 
#' @export
#' 
#' @examples 
#' # set up a gene set data frame
#' ann_df = data.frame(pathway = c('a', 'a', 'a', 'b', 'c'),
#'                     gene = c('gene_1', 'gene_2', 'gene_3', 'gene_1', 
#'                              'gene_3'))
#' ann = as_annotation_list(ann_df, 'pathway', 'gene')
#' str(ann)
#' 
#' # go = read_gpa("goa_human.gpa.gz")
#' # ann = as_annotation_list(go, "GO ID", "DB Object ID")
as_annotation_list = function(annotations, set_col, gene_col) {
  ann_list = unstack(annotations[, c(gene_col, set_col)])
  # get unique items only
  sapply(ann_list, unique)
}
