#' Read CF-MS data from MaxQuant
#' 
#' Read in a 'proteinGroups.txt' file output by MaxQuant, and convert it to a 
#' CF-MS chromatogram matrix, in which each protein is a row and each fraction 
#' is a column.
#' 
#' @param filepath path to the proteinGroups.txt file
#' @param identifiers one of \code{'protein groups'} or \code{'genes'}; if the
#'   latter, protein groups will be matched to gene names (using the 
#'   'Gene names' column in the proteinGroups file), keeping only the best row
#'   per gene when there are redundant gene names)
#' @param quant_mode the protein quantitations to use; one of \code{'iBAQ'}, 
#'   \code{'MS1 intensity'}, \code{'MaxLFQ'}, \code{'spectral counts'}, or 
#'   \code{'ratio'}. Note that in ratiometric experiments with more than one
#'   isotopic label, the output matrix may need to be manually split further
#'   (e.g. to separate heavy and medium over light ratios).
#'   
#' @return a numeric matrix with proteins in rows and fractions in columns,
#'   and a \code{'metadata'} attribute containing metadata from the 
#'   proteinGroups file for each row
#'  
#' @importFrom dplyr filter select starts_with first bind_rows first
#' @importFrom tidyr replace_na
#' @importFrom stringr str_count
#' @importFrom magrittr %>% %<>% extract set_colnames
#' @importFrom utils read.delim
#' 
#' @export
read_maxquant = function(filepath, identifiers = c('protein groups', 'genes'),
                         quant_mode = c('iBAQ', 'MS1 intensity', 'MaxLFQ', 
                                        'spectral counts', 'ratio')) {
  identifiers = match.arg(identifiers)
  quant_mode = match.arg(quant_mode)
  
  # read data
  dat = read.delim(filepath, stringsAsFactors = FALSE) %>%
    # remove contaminants, reverse hits
    replace_na(list(Only.identified.by.site = "",
                    Reverse = "",
                    Potential.contaminant = "")) %>%
    filter(Only.identified.by.site != '+',
           Reverse != '+',
           Potential.contaminant != '+')
  if (nrow(dat) == 0)
    stop('something went wrong: no proteins detected in file ', filepath)
  
  # extract metadata
  meta = dat %>%
    dplyr::select(Protein.IDs:Unique.peptides)
  if (!"Gene.names" %in% colnames(meta) & identifiers == 'genes') {
    stop("can't map to genes: column \"Gene.names\" does not exist")
  }
  
  # extract chromatograms
  if (quant_mode == 'iBAQ') {
    mat = dat %>%
      dplyr::select(starts_with('iBAQ.')) %>% 
      as.matrix()
  } else if (quant_mode == 'MS1 intensity') {
    mat = dat %>%
      dplyr::select(starts_with('Intensity.')) %>% 
      as.matrix()
  } else if (quant_mode == 'MaxLFQ') {
    mat = dat %>%
      dplyr::select(starts_with('LFQ.')) %>% 
      as.matrix()
  } else if (quant_mode == 'spectral counts') {
    mat = dat %>%
      dplyr::select(starts_with('MS.MS.count.')) %>% 
      as.matrix() 
  } else if (quant_mode == 'ratio') {
    mat = dat %>%
      dplyr::select(starts_with('Ratio.')) %>% 
      extract(, !grepl("iso|type|norm|var|count", colnames(.))) %>%
      as.matrix()
    keep = str_count(colnames(mat), '\\.') >= 3
    mat %<>% extract(, keep)
  }
  ## handle ratiometric datasets
  keep = grepl("\\.L\\.|\\.M\\.|\\.H\\.", colnames(mat))
  if (sum(keep) > 0) {
    mat %<>% extract(, keep)
  }
  
  # set rownames
  rownames(mat) = meta$Majority.protein.IDs
  
  # optionally, map protein groups to genes
  if (identifiers == 'genes') {
    # keep one row per gene
    gene_map = meta %>%
      dplyr::select(Majority.protein.IDs, Gene.names) %>%
      set_colnames(c("protein_group", "gene")) %>%
      mutate(gene = strsplit(gene, ';')) %>%
      unnest(gene) %>%
      drop_na()
    genes = unique(gene_map$gene)
    gene_mat = matrix(NA, nrow = length(genes), ncol = ncol(mat),
                      dimnames = list(genes, colnames(mat)))
    n_fractions = rowSums(!is.na(mat) & is.finite(mat) & mat != 0)
    out_map = data.frame()
    for (gene in genes) {
      protein_groups = gene_map$protein_group[gene_map$gene == gene]
      # pick the best protein for this replicate
      n_fractions0 = n_fractions[protein_groups]
      best = names(which(n_fractions0 == max(n_fractions0))) %>% dplyr::first()
      gene_mat[gene, ] = mat[best, ]
      # save the mapping
      out_map %<>% bind_rows(data.frame(gene = gene, protein_group = best))
    }
    
    ## now drop duplicated rows
    mat = gene_mat
    order = rownames(mat) %>% order()
    mat %<>% extract(order, )
    out_map %<>% extract(order, )
    drop = which(duplicated(mat))
    mat %<>% extract(-drop, )
    out_map %<>% extract(-drop, )
    
    # set attribute on the gene matrix
    attr(mat, 'gene_map') = out_map
  } else {
    # set metadata as attributes
    attr(mat, 'metadata') = meta
  }
  
  return(mat)
}
