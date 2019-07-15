#' Simulate data from an efast matrix
#'
#' Simulate a dataset as in the simulations of Van Kesteren & Kievit (2019).
#' There are 17 regions of interest, measured in both the left and right
#' hemisphere. These ROIs have a predefined amount of correlation over and
#' above that expected by only the underlying factors.
#'
#' @param N <int> Sample size
#' @param lam_lat <numeric> factor loading for the lateralised factor
#' @param lam_bil <numeric> factor loading for the bilateral factors
#' @param psi_cov <numeric> covariances of latent variables in (0, 1)
#' @param cor_uniq <numeric> residual correlation
#'
#' @return data frame with 17 regions of interest, bilaterally measured with
#' 4 underlying factors and contralateral homology.
#'
#' @references Van Kesteren, E. J., & Kievit, R. K. (2019) Exploratory factor
#'     analysis with structured residuals applied to brain morphology.
#'
#' @export
simulate_efast <- function(N = 650L, lam_lat = 0.595, lam_bil = 0.7,
                           psi_cov = 0.5, cor_uniq = 0.4) {
  # 17 regions of interest per hemisphere, 5 latent variables
  lat_names <- c("LH", "F1", "F2", "F3")
  obs_names <- c(paste0("lh_ROI", 1:17), paste0("rh_ROI", 1:17))

  # lambda
  lambda <- matrix(0, nrow = 34, ncol = 4)
  colnames(lambda) <- lat_names
  rownames(lambda) <- obs_names

  lambda[1:8,   1] <- lam_lat

  lambda[9:11,  2] <- lam_bil
  lambda[12:14, 3] <- lam_bil
  lambda[15:17, 4] <- lam_bil

  lambda[18:20, 2] <- lam_bil
  lambda[21:23, 3] <- lam_bil
  lambda[24:25, 4] <- lam_bil

  lambda[26:28, 2] <- lam_bil
  lambda[29:31, 3] <- lam_bil
  lambda[32:34, 4] <- lam_bil

  # psi
  psi <- matrix(psi_cov, 4, 4)
  diag(psi) <- 1
  colnames(psi) <- rownames(psi) <- lat_names

  trait_matrix <- lambda %*% psi %*% t(lambda)
  diag(trait_matrix) <- 1 # uniquenesses

  methd_matrix <- matrix(0, 34, 34)
  diag(methd_matrix[1:17, 18:34]) <- cor_uniq
  diag(methd_matrix[18:34, 1:17]) <- cor_uniq

  dat <- MASS::mvrnorm(N, mu = rep(0, 34), Sigma = trait_matrix + methd_matrix)
  colnames(dat) <- obs_names
  as.data.frame(dat)
}
