#' Read complexes from the CORUM database
#' 
#' Utility function to read complexes from the 'coreComplexes.txt' file 
#' distributed through the CORUM database into a two-column data frame, mapping
#' each complex (column 1) to its subunits (column 2). 
#'
#' @param filepath path to the 'coreComplexes.txt' file
#' @param identifiers the set of protein identifiers to parse from the file;
#'   one of \code{'uniprot'}, \code{'entrez'}, or \code{'gene_name'}
#' 
#' @return a data frame with protein complex names in column 1 and individual
#'   subunits in column 2
#' 
#' @importFrom dplyr select mutate distinct group_by filter ungroup
#' @importFrom tidyr unnest 
#' @importFrom magrittr set_colnames %>%
#' 
#' @export
read_corum = function(filepath, 
                      identifiers = c('uniprot', 'entrez', 'gene_name')) {
  identifiers = match.arg(identifiers)
  
  # file should be named coreComplexes.txt 
  # (might be compressed)
  if (!grepl("coreComplexes\\.txt", basename(filepath)))
    warning("reading file: ", basename(filepath), 
            ", was expecting coreComplexes.txt")
  
  # read data
  dat = read.delim(filepath, stringsAsFactors = FALSE)
  # tag identifier
  column = switch(identifiers,
                  uniprot = 'subunits.UniProt.IDs.',
                  entrez = 'subunits.Entrez.IDs.',
                  gene_name = 'subunits.Gene.name.',
                  stop("not sure what to do with identifiers: ", identifiers))
  dat$identifiers = dat[[column]]
  
  # unnest subunits
  complexes = dat %>%
    dplyr::select(ComplexName, identifiers) %>%
    set_colnames(c('complex', 'subunit')) %>%
    mutate(subunit = strsplit(subunit, ';')) %>%
    unnest(subunit) %>%
    distinct() %>%
    # remove complexes with only one protein
    group_by(complex) %>%
    filter(n_distinct(subunit) > 1) %>%
    ungroup() %>%
    # remove empty strings
    filter(subunit != "")
  
  return(complexes)
}
