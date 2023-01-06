#' Download and process the GSE data from Gene Expression Omnibus (GEO) by its accession number
#'
#' This function downloads and process gene expression data stored in the GSE format from GEO into a count matrix using
#' its accession number. This function makes use of the GEOquery package built by Sean Davis.
#' 
#' 

get_GSE <- function(acc_number, take_log = FALSE) {
  
  ## download the data and the series matrix
  gse <- getGEO(acc_number, GSEMatrix = FALSE)
  gse_series <- getGEO(acc_number, GSEMatrix = TRUE, AnnotGPL = TRUE)[[paste(acc_number, 'series_matrix.txt.gz', sep = '_')]]
  
  ## inspect GSM and GPL elements
  # names(GEOquery::GSMList(gse))
  # names(GEOquery::GPLList(gse))
  ## return an error is there is more than one platforms
  pl_list <- GEOquery::GPLList(gse)
  if (length(names(pl_list)) > 1) stop('More than one platforms found in the GSE data.')
  
  ## get the platforms of each GSM (since each GSM contains one platform)
  # gse_gsmplatforms <- lapply(GSMList(gse), function(x) {Meta(x)$platform_id})

  ## only take GSMs from one platform
  gse_gsmlist <- BiocGenerics::Filter(function(gsm) {Meta(gsm)$platform_id == names()[1]}, GSMList(gse))
  # GEOquery::Table(gse_gsmlist[[1]])[1:5, ] # we can see the value of microarray expression values
  
  ## get the probeset
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
  if (take_log) data_matrix <- log2(data_matrix) # log2-transform if take_log is TRUE
  
  
}