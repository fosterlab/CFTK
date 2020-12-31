#' Read a GAF file 
#' 
#' Read a GO annotation file in GAF format. Further information about the GAF
#' format is available from the 
#' \href{http://geneontology.org/page/go-annotation-file-gaf-format-21}{Gene Ontology Consortium}. 
#' 
#' @param filepath the location of the GAF file 
#' @param filter_NOT if true, filter annotations with the qualifier NOT
#' @param filter_evidence optionally, specify evidence codes to filter. By 
#'   default, evidence codes ND, IPI, IEA and NAS are filtered. 
#' @param ontology optionally, provide the Gene Ontology itself as an object
#'   read in by the \code{ontologyIndex} package from an OBO file
#' @param propagate if true, and an ontology file is provided, all ancestors of 
#'   a given term are associated with each protein.
#' 
#' @return a data.frame containing the filtered GAF file  
#' 
#' @importFrom readr read_tsv
#' @importFrom dplyr filter
#' 
#' @export
#' 
#' @examples 
#' # read human GOA and filtering IEA annotations
#' # go = read_gaf("goa_human.gaf.gz", filter_evidence = "IEA")
#' # read mouse GOA wthout filtering any annotations
#' # go = read_gaf("goa_mouse.gaf.gz", filter_evidence = NULL)
read_gaf = function(filepath,
                    filter_NOT = TRUE, 
                    filter_evidence = NULL,
                    ontology = NULL, 
                    propagate = TRUE) {
  goa = suppressWarnings(suppressMessages(
    readr::read_tsv(filepath, comment = "!", col_names = gaf_colnames())))
  
  # optionally, filter annotations with the qualifier NOT
  if (filter_NOT) 
    goa %<>% filter(!grepl("NOT", Qualifier))
  
  # optionally, filter out evidence 
  if (!is.null(filter_evidence) && !is.na(filter_evidence) &&
      length(filter_evidence) > 0) 
    goa %<>% filter(!`Evidence Code` %in% filter_evidence)
  
  # read the ontology
  if (!is.null(ontology) & propagate) {
    if (!"ontology_index" %in% class(ontology))
      stop("Ontology must be of class ontology_index")
    goa %<>%
      mutate(ancestors = ontology$ancestors[`GO ID`]) %>%
    # filter out terms missing ancestors (these are deprecated)
      filter(lengths(ancestors) > 0) %>%
      unnest(ancestors)
    
    # replace column
    goa[["GO ID"]] = goa[["ancestors"]]
    goa = goa[, -ncol(goa)]
  }
  return(goa)
}

#' Column names for a GAF file
#' 
#' Column names for a tab-delimited annotation file in GAF format. 
#' Further information about the GAF format is available from the 
#' \href{http://geneontology.org/page/go-annotation-file-gaf-format-21}{
#' Gene Ontology Consortium}. 
#' 
#' @export
gaf_colnames = function() {
  # define column names for a GOA file in GAF format
  c('DB',
    'DB Object ID',
    'DB Object Symbol',
    'Qualifier',
    'GO ID',
    'DB:Reference',
    'Evidence Code',
    'With (or) From',
    'Aspect',
    'DB Object Name',
    'DB Object Synonym',
    'DB Object Type', 
    'Taxon',
    'Date',
    'Assigned By',
    'Annotation Extension',
    'Gene Product Form ID')
}

#' Roots of the Gene Ontology
#' 
#' Return a vector containing the Gene Ontology codes for 'biological process',
#' 'molecular function', and 'cellular compartment'. 
#' 
#' @export
go_roots = function() {
  c(BP = "GO:0008150", 
    CC = "GO:0005575",
    MF = "GO:0003674")
}
