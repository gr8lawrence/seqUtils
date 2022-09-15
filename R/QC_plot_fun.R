#' Plot Expression Value Boxplots for Samples
#'
#' @param Y a matrix of expression values with row names corresponding to genes and column names to sample names
#' @param n the number of samples used to plot the 
#' 
#' @return boxplots for the expression values with each representing those of a sample
#' 
#' @examples 
#' data(example)
#' plot.bp(Y = example)
#' plot.bp(Y = example, n = 5)
#' 
#' @importFrom dplyr as_tibble
#' @importFrom tidyr pivot_longer
#' @export

plot.bp <- function(Y, n = 10) {
  if (n > 10) {
    message("It is recommended to only plot the first few samples (n <= 10) for best visualization.")
  }
  df = dplyr::as_tibble(Y[, seq_len(n)])
  df$Gene = rownames(Y)
  long_df = tidyr::pivot_longer(data = df, 
                                cols = seq(1, ncol(df) - 1), 
                                names_to = "samples", 
                                values_to = "nc")
  long_df$samples = factor(long_df$samples, levels = colnames(Y), ordered = TRUE)
  long_df = long_df[which(long_df$nc >= 0), ] # such intensities could be negative if using microarray
  long_df$log10nc = log(long_df$nc + 1, base = 10) 
  p = ggplot2::ggplot(data = long_df, aes(x = samples, y = log10nc)) +
    geom_boxplot() +
    labs(y = bquote("log"[10] ~ "(count + 1)"),
         x = 'Samples')
  
  return(p)
}

#' Plot Pairwise Scatterplot for Samples
#'
#' @param Y a matrix of expression values with row names corresponding to genes and column names to sample names
#' @param notes additional notes for file names, such as generations, platforms, sample information, etc. Default is \code{NULL}
#' @param dir_name name of the directory to which the output jpeg files are saved. Default is the working directory (\code{getwd()})
#' @param par_param a vector of 2 specifying the number of rows and columns of the output plots. Default is \code{c(3, 3)}
#'
#' @return output no direct return but save the pairwise scatterplots in a jpeg file to the named directory
#' 
#' @examples 
#' data(example)
#' pw.scatter(Y = example[, 1:5])
#' pw.scatter(Y = example[, 1:5], notes = 'example')
#' pw.scatter(Y = example[, 1:5], notes = 'example', par_param = c(2, 3))
#'
#' @export

pw.scatter <- function(Y, notes = NULL, dir_name = '.', par_param = c(3, 3)) {
  n = ncol(Y)
  if (is.null(notes)) {
    jpeg(paste(dir_name, "scatterplots%03d.jpeg", sep = '/'),
         width = 1080, height = 1080, units = "px", pointsize = 24)
  } else {
    jpeg(paste(dir_name, paste(notes, "scatterplots%03d.jpeg", sep = '_'), sep = '/'),
         width = 1080, height = 1080, units = "px", pointsize = 24)
  }
  par(mfrow = par_param)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      col1 = Y[, i] 
      col2 = Y[, j]
      name1 = colnames(Y)[i]
      name2 = colnames(Y)[j]
      p_cor = round(cor(col1, col2, method = "pearson"), 3)
      plot(log10(col1 + 1), log10(col2 + 1), col = "blue", xlab = name1, ylab = name2, 
           main = paste(name1, "vs.", name2), sub = paste("Pearson's cor =", p_cor, "for original counts"), col.sub = "red",
           xlim = c(0, 5), ylim = c(0, 5))
      abline(a = 0, b = 1, col = "orange")
    }
  }
  dev.off()
}