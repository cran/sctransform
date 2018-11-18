#' Denoise data by setting all latent factors to their median values and reversing the regression model
#'
#' @param x A list that provides model parameters and optionally meta data; use output of vst function
#' @param data The name of the entry in x that holds the data
#' @param cell_attr Provide cell meta data holding latent data info
#' @param do_round Round the result to integers
#' @param do_pos Set negative values in the result to zero
#' @param show_progress Whether to print progress bar
#'
#' @return De-noised data as UMI counts
#'
#' @importFrom stats median model.matrix as.formula
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
#'
#' @examples
#' vst_out <- vst(pbmc, return_cell_attr = TRUE)
#' umi_denoised <- denoise(vst_out)
#'
denoise <- function(x, data = 'y', cell_attr = x$cell_attr, do_round = TRUE, do_pos = TRUE,
                    show_progress = TRUE) {
  if (is.character(data)) {
    data <- x[[data]]
  }
  # when denoising, set all latent variables to median values
  cell_attr[, x$arguments$latent_var] <- apply(cell_attr[, x$arguments$latent_var, drop=FALSE], 2, function(x) rep(median(x), length(x)))
  regressor_data <- model.matrix(as.formula(gsub('^y', '', x$model_str)), cell_attr)

  genes <- rownames(data)
  bin_size <- x$arguments$bin_size
  bin_ind <- ceiling(x = 1:length(x = genes) / bin_size)
  max_bin <- max(bin_ind)
  if (show_progress) {
    message('Computing de-noised UMI count matrix')
    pb <- txtProgressBar(min = 0, max = max_bin, style = 3)
  }
  denoised_data <- matrix(NA, length(genes), nrow(regressor_data), dimnames = list(genes, rownames(regressor_data)))
  for (i in 1:max_bin) {
    genes_bin <- genes[bin_ind == i]
    pearson_residual <- data[genes_bin, ]
    coefs <- x$model_pars_fit[genes_bin, -1]
    theta <- x$model_pars_fit[genes_bin, 1]
    mu <- exp(tcrossprod(coefs, regressor_data))
    variance <- mu + mu^2 / theta
    denoised_data[genes_bin, ] <- mu + pearson_residual * sqrt(variance)
    if (show_progress) {
      setTxtProgressBar(pb, i)
    }
  }
  if (show_progress) {
    close(pb)
  }

  if (do_round) {
    denoised_data <- round(denoised_data, 0)
  }
  if (do_pos) {
    denoised_data[denoised_data < 0] <- 0
  }
  return(denoised_data)
}
