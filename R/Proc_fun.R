#' RNA Sequencing Data Manual Filtering
#' 
#' Manual filtering of unqualified rows from an RNA-seq dataset. The parameters are set
#' for bulk RNA seq data. Four filtering steps are performed to remove unqualified rows:
#' * Remove non-genes (i.e. rows for summary statistics, etc.)
#' * Remove genes with low expression across samples (sum of expression < `min_sum`)
#' * Remove genes with low variance across samples (variance < `min_var`)
#' * Remove duplicate rows (i.e. rows with duplicate gene names, etc.)
#' 
#' Although software like DESeq2 has steps like independent filtering that allows users
#' to skip manual pre-filtering, it is still recommended to do so sometimes.
#' 
#' @param raw raw sequencing data in the wide format 
#' @param min_sum minimal sum of gene expression over all samples to be allowed in the data; default value is `10`
#' @param min_var minimal gene expression variance over all samples to be allowed in the data; default value is `1e-4`
#' @param verbose logical value to indicate whether detailed information of filtering should be printed; default value is `TRUE`
#' 
#' @return a filtered raw sequencing data based on the specified values
#' 
#' @examples 
#' data(example)
#' example[1, ] = runif(ncol(example), 0, 0.01) # to create an example gene with low expression
#' example[2, ] = 15 # an example of gene with no variation
#' d1 = man.filter(example) # default setting
#' d2 = man.filter(example, min_sum = 0) # keep the new row 1
#' d3 = man.filter(example, min_var = 0) # keep the new row 2
#' d4 = man.filter(example, verbose = FALSE) # same as d1 in the output but no information printed
#' 
#' @export

man.filter <- function(raw, min_sum = 10, min_var = 1e-4, verbose = TRUE) {
  all_ids <- NULL
  ## remove non-genes in the rawrix
  bad_name_ids <- which(substr(rownames(raw), 1, 1) == "_")
  all_ids <- union(all_ids, bad_name_ids)
  
  ## keep genes with more than min_sum reads
  no_read_ids <- which(rowSums(raw) < min_sum)
  all_ids <- union(all_ids, no_read_ids)
  
  ## remove genes with no variance
  no_var_ids <- which(rowSums((raw - rowMeans(raw))^2) < min_var)
  all_ids <- union(all_ids, no_var_ids)
  
  ## remove duplicate rows
  dup_ids <- which(duplicated.array(raw))
  all_ids <- union(all_ids, dup_ids)
  
  ## remove grand totals
  gt_ids <- which(rownames(raw) == "Grand Total")
  all_ids <- union(all_ids, gt_ids)
  
  ## filter rows
  if (length(all_ids) > 0) raw <- raw[-all_ids, ]
  
  ## print messages
  if (verbose) {
    if (length(all_ids) > 0) {
      message(paste0("A total of ", length(all_ids), " rows are removed."))
      message("Detailed removal info is below. There maybe rows that fit multiple criteria below:")
      if (length(bad_name_ids) > 0) message(paste0(length(bad_name_ids), " non-gene rows are removed."))
      if (length(no_read_ids)  > 0) message(paste0(length(no_read_ids), " rows with small sum of counts (<", min_sum , ") are removed."))
      if (length(no_var_ids)  > 0) message(paste0(length(no_var_ids), " rows with small variance of counts (<", min_var , ") are removed."))
      if (length(dup_ids) > 0) message(paste0(length(dup_ids), " duplicate rows are removed."))
      if (length(gt_ids) > 0) message(paste0(length(gt_ids), " row for the grand totals are removed."))
    }
  }
  
  return(raw)
}
