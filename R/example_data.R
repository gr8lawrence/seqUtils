#' Toy Example Expression Value Data
#' 
#' This toy example is for demonstrating the format in which expression values 
#' should be presented for using functions in this package. It is generated using
#' the following code:
#' 
#' \code{example <- matrix(round(rchisq(1000, 10^3)), 100, 10)}
#' \code{rownames(example) <- paste('gene', seq_len(100), sep = '_')}
#' \code{colnames(example) <- paste('sample', seq_len(10), sep = '_')}
#' 
#' @docType data
#' 
#' @usage data(example)
#' 
#' @format an object of class \code{"matrix"} and \code{"array"}
#' 
#' @keywords example dataset
#' 
#' @examples 
#' data(example)
#' dim(example)
#' rownames(example)[1:10]
#' colnames(example)[1:10]

"example"