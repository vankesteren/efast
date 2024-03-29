
#' Check whether an object is an efast model
#'
#' @param x any object
#'
#' @export
is_efast <- function(x) {
  if (!inherits(x, "lavaan")) return(FALSE)
  if (is.null(x@external$type)) return(FALSE)
  if (substr(x@external$type, 1, 5) == "EFAST") return(TRUE)
  FALSE
}

#' Check whether an object is an efast-hemi model
#'
#' @param x any object
#'
#' @export
is_efast_hemi <- function(x) {
  if (!inherits(x, "lavaan")) return(FALSE)
  if (is.null(x@external$type)) return(FALSE)
  if (x@external$type == "EFAST_HEMI") return(TRUE)
  FALSE
}

#' Check whether an object is an efast-efa model
#'
#' @param x any object
#'
#' @export
is_efast_efa <- function(x) {
  if (!inherits(x, "lavaan")) return(FALSE)
  if (is.null(x@external$type)) return(FALSE)
  if (x@external$type == "EFAST_EFA") return(TRUE)
  FALSE
}

#' Get loadings
#'
#' @param fit efa or efast model
#' @param symmetry whether to display the left and right hemisphere loadings side-by-side
#'
#' @export
efast_loadings <- function(fit, symmetry = FALSE) {
  if (!is_efast(fit)) stop("Enter an efast model.")
  if (symmetry && !is_efast_hemi(fit))
    stop("Enter an efast-hemi model to use the symmetry option.")

  glist         <- fit@Model@GLIST
  efa_idx       <- fit@Model@lv.efa.idx[[1]]$efa1
  M             <- length(efa_idx)
  lam           <- glist$lambda[, efa_idx, drop = FALSE]
  colnames(lam) <- paste0("F", 1:M)
  rownames(lam) <- fit@Model@dimNames[[1]][[1]]
  if (symmetry) {
    # find the number of variables per hemisphere
    nn     <- fit@Model@nvar / 2
    lh_lam <- lam[1:nn,]
    rh_lam <- lam[(nn + 1):(2*nn),]
    lam    <- matrix(0, nn, M*2)
    for (i in 1:M) {
      lam[,i*2 - 1] <- lh_lam[,i]
      lam[,i*2]     <- rh_lam[,i]
    }
    colnames(lam) <- paste0(
      rep(paste0("F", 1:M), each = 2),
      rep(c("_lh", "_rh"), M)
    )
    rownames(lam) <- fit@Model@dimNames[[3]][[2]][-efa_idx]
  }
  return(structure(lam, class = "loadings"))
}

#' Get lateralization per ROI
#'
#' The lateralization index (LI) is calculated as (sum of the uniquenesses for
#' the contralateral homologues) / (total residual variance after accounting for
#' the exploratory factors). Thus, if there is no symmetry at all, the LI is 1,
#' and if all the residual variance is explained by symmetry, the LI is 0.
#'
#' @param fit an efast model
#'
#' @return parameter estimates table of lateralization per ROI
#'
#' @export
lateralization <- function(fit) {
  # This function computes the per-region lateralisation index
  # It is computed as sum(residual variance of homologous regions) /
  # (sum(total variance of hr) - sum(trait variance of hr))
  # it is thus the variance that cannot be explained by symmetry
  pe     <- lavaan::parameterestimates(fit)
  li_idx <- stringr::str_detect(pe$label, "LI_")
  li_pe  <- pe[li_idx, -1:-4]
  li_pe$label <- stringr::str_extract(li_pe$label, "(?<=_).*")
  colnames(li_pe)[1] <- "region"
  rownames(li_pe) <- 1:nrow(li_pe)
  return(li_pe)
}

#' Get correlation matrix decomposition
#'
#' @param fit an efast model
#'
#' @return The decomposition returns four components: the observed correlation
#' matrix, the factor-implied covaraiance, the residual variance, and the
#' structural covariance matrix. You can plot these matrices using corrplot.
#'
#' @export
decomposition <- function(fit) {
  if (!is_efast(fit)) stop("Enter an EFAST model.")
  glist   <- fit@Model@GLIST
  efa_idx <- fit@Model@lv.efa.idx[[1]]$efa1

  # factor implied
  H     <- fit@Model@H[[1]][efa_idx, efa_idx]
  lam_e <- glist$lambda[, efa_idx]
  psi_e <- crossprod(H)
  factor_implied <- tcrossprod(lam_e %*% psi_e, lam_e)

  # residual
  theta <- glist$theta
  residual_variance <- diag(diag(theta))

  # structure
  structure_component <- theta - residual_variance


  return(list(observed  = fit@SampleStats@cov[[1]],
              factor    = factor_implied,
              residual  = residual_variance,
              structure = structure_component))
}
