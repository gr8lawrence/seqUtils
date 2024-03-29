#' Download and process the microarray data from Gene Expression Omnibus (GEO) by its GSE accession number
#'
#' This function downloads and processes microarray gene expression data stored in the GSE format from GEO into a count matrix using
#' its accession number. This function makes use of the GEOquery package built by Sean Davis.
#' 
#' @param acc_number the GEO accession number of the microarray data.
#' @param take_log whether we want to log2-transform the expression matrix. The default is `FALSE`.
#' @param ... other parameters to be passed to [GEOquery::getGEO].
#' 
#' @return a list containing two fields: `data_matrix` gives the genotype expression data in a matrix form; 
#' `p_data` gives the phenotype data in a data frame.
#' 
#' @importFrom GEOquery getGEO
#' @importFrom GEOquery Meta
#' @importFrom GEOquery GSMList
#' @importFrom GEOquery GPLList
#' @importFrom GEOquery Table
#' @importFrom BiocGenerics Filter
#' 
#' @examples 
#' ls <- get.GSE.microarray('GSE61754')
#' ls
#' 
#' @export

get.GSE.microarray <- function(acc_number, take_log = FALSE, ...) {
  
  ## download the data and the series matrix
  gse <- GEOquery::getGEO(acc_number, GSEMatrix = FALSE, ...)
  gse_series <- GEOquery::getGEO(acc_number, GSEMatrix = TRUE, AnnotGPL = TRUE, ...)[[paste(acc_number, 'series_matrix.txt.gz', sep = '_')]]
  
  ## return an error is there is more than one platforms in this GSE platform
  pl_list <- GEOquery::GPLList(gse)
  if (length(names(pl_list)) > 1) stop('More than one platforms found in the GSE data.')

  ## only take GSMs from one platform
  gse_gsmlist <- BiocGenerics::Filter(function(gsm) {GEOquery::Meta(gsm)$platform_id == names(pl_list)[1]}, GEOquery::GSMList(gse))

  ## get the active probesets (for microarrays)
  gse_probesets <- GEOquery::Table(pl_list[[1]])$ID
  gse_probe_subsets <- intersect(gse_probesets, GEOquery::Table(gse_gsmlist[[1]])$ID_REF)
  
  ## make the data matrix
  data_matrix <- do.call('cbind', lapply(gse_gsmlist, 
                                         function(x) {tab = GEOquery::Table(x)
                                         mymatch = match(gse_probe_subsets, tab$ID_REF)
                                         return(tab$VALUE[mymatch])}
  ))
  data_matrix <- apply(data_matrix, 2, function(x) {as.numeric(as.character(x))})
  data_matrix <- data_matrix[, which(!is.na(colSums(data_matrix)))]
  rownames(data_matrix) <- gse_probe_subsets
  if (take_log) data_matrix <- log2(data_matrix) # log2-transform if take_log is TRUE
  
  ## getting the phenotype data
  p_data_names <- names(gse_series @ phenoData @ data)
  p_data <- (gse_series @ phenoData @ data)[, grep(':ch1', p_data_names)]
  
  ## return a list that contains both genotype and phenotype data
  data_ls <- list(data_matrix = data_matrix, p_data = p_data)
  data_ls
}