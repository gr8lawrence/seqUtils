#' Plot Expression Value Boxplots for Samples
#'
#' @param Y a matrix of expression values with row names corresponding to genes and column names to sample names
#' @param n the number of samples used to plot the boxplots. Default is 10
#' 
#' @return boxplots for the expression values with each representing those of a sample
#' 
#' @examples 
#' data(example)
#' expr.bp(Y = example)
#' expr.bp(Y = example, n = 5)
#' 
#' @importFrom dplyr as_tibble
#' @importFrom tidyr pivot_longer
#' @import ggplot2
#' @export

expr.bp <- function(Y, n = 10) {
  if (n > 10) {
    message("It is recommended to only plot the first few samples (n <= 10) for best visualization.")
  }
  df <- dplyr::as_tibble(Y[, seq_len(n)])
  df$Gene <- rownames(Y)
  long_df <- tidyr::pivot_longer(data = df, 
                                cols = seq(1, ncol(df) - 1), 
                                names_to = "samples", 
                                values_to = "nc")
  long_df$samples <- factor(long_df$samples, levels = colnames(Y), ordered = TRUE)
  long_df <- long_df[which(long_df$nc >= 0), ] # such intensities could be negative if using microarray
  long_df$log10nc <- log(long_df$nc + 1, base = 10) 
  p <- ggplot2::ggplot(data = long_df, ggplot2::aes(x = samples, y = log10nc)) +
    ggplot2::geom_boxplot() +
    ggplot2::labs(y = bquote("log"[10] ~ "(count + 1)"),
         x = 'Samples')
  
  return(p)
}

#' Plot Pairwise Scatterplot for Samples
#'
#' @param Y a matrix of expression values with row names corresponding to genes and column names to sample names
#' @param notes additional notes for file names, such as generations, platforms, sample information, etc. Default is \code{NULL}
#' @param dir_name name of the directory to which the output jpeg files are saved. Default is the working directory (\code{getwd()})
#' @param par_param a vector of 2 specifying the number of rows and columns of the output plots. Default is \code{c(3, 3)}
#' @param ... other parameters for the [base::plot] function
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

pw.scatter <- function(Y, notes = NULL, dir_name = '.', par_param = c(3, 3), ...) {
  n <- ncol(Y)
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
      col1 <- Y[, i] 
      col2 <- Y[, j]
      name1 <- colnames(Y)[i]
      name2 <- colnames(Y)[j]
      p_cor <- round(cor(col1, col2, method = "pearson"), 3)
      plot(log10(col1 + 1), log10(col2 + 1), col = "blue", xlab = name1, ylab = name2, 
           main = paste(name1, "vs.", name2), sub = paste("Pearson's cor =", p_cor), col.sub = "red",
           ...)
      abline(a = 0, b = 1, col = "orange")
    }
  }
  dev.off()
}

#' Plot a Heat map for Pairwise Correlation
#' 
#' Plot the pairwise correlations of columns of a given matrix in a heat map.
#' 
#' @param Y the matrix whose pairwise correlations between columns we are interested in.
#' @param inds the array of column indices of `Y` we wish to plot. The default is all columns of `Y`.
#' @import ggplot2
#' 
#' @examples
#' Y <- matrix(c(1:4, 1:4, 4:1), 4, 3)
#' p <- pw.cor.heatmap(Y)
#' p2 <- pw.cor.heatmap(Y, c(1, 3))
#' p
#' p2
#' 
#' @export

pw.cor.heatmap <- function(Y, inds = seq_len(ncol(Y))) {
  cor_mat <- cor(Y[, inds])
  rownames(cor_mat) <- colnames(cor_mat) <- colnames(Y)
  cor_df <- reshape2::melt(t(cor_mat))
  cor_plot <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Var1, y = Var2, fill = value)) + 
    ggplot2::geom_tile() +
    ggplot2::theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) +
    ggplot2::scale_fill_gradient(low = 'yellow', high = 'red')
  cor_plot
}

#' Plot a Heatmap for an Expression Value Matrix
#' 
#' This function will plot a heatmap for an expression value matrix. It will automatically group genes
#' that are highly expressed for each column so the diagonal blocks of the heatmap show the highest expression values.
#' This function will be particularly helpful for constructing a heatmap of a signature matrix or a differential expression map.
#' 
#' 
#' @param sig_mat The matrix for which the heatmap is plotted.
#' @param plot_title Custom title for thie plot. Default is 'Signature Matrix Heatmap'.
#' @param names The names of the columns of the final heatmap. Default is `colnames(sig_name)`.
#' 
#' @export

make.heatmap <- function(sig_mat, plot_title = 'Signature Matrix Heatmap', names = colnames(sig_name)) {
  sig_ordered <- rearrange.matrix(sig_mat) 
  p_sig <- make.theta.heat.map(sig_ordered, plot_title)
  p_sig
}

#' Function to Rearrange The Expression Matrix
#'

rearrange.matrix <- function(Theta_star) {
  inds <- apply(Theta_star, 1, function(x) which(x == max(x)))
  genes_ordered <- names(inds)[order(inds, decreasing = TRUE)]
  Theta_ordered <- Theta_star[genes_ordered, ]
  Theta_ordered
}

#' Function to Produce The Heatmap from a Matrix Using `ggplot2`
#' 
#' @importFrom reshape2 melt
#' @import ggplot2

make.theta.heat.map <- function(Theta_star, title, all_names) {
  heat_df <- reshape2::melt(Theta_star) 
  colnames(heat_df) <- c('gene', 'ct', 'value')
  p <- ggplot(heat_df, aes(x = ct, y = gene, fill = log10(value + 1))) +
    geom_tile() +
    labs(title = title,
         x = 'Cell type',
         y = 'Gene feature') +
    theme(axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = 18)) +
    scale_fill_distiller(palette = 'YlGnBu')
  p
}


