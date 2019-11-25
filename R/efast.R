#' Estimate an EFAST-hemi model
#'
#' This function estimates efast models with covariance due to hemispheric
#' symmetry.
#'
#' @param data <data.frame> the dataset
#' @param M <numeric> How many factors, minimum 2
#' @param lh_idx <numeric> column numbers of left hemisphere variables
#' @param rh_idx <numeric> column numbers of right hemisphere variables
#' @param roi_names <character> optional names of rois
#' @param constrain <bool> whether to constrain the symmetry (see details)
#' @param auto.fix.first <bool> see lavaan
#' @param auto.var <bool> see lavaan
#' @param auto.efa <bool> see lavaan
#' @param information <character> see lavaan
#' @param std.ov <bool> see lavaan
#' @param ... other arguments passed to lavaan
#'
#' @details The constrained model constrains the residual covariance to be
#'   equal across the different ROIs.
#'
#' @examples
#' \dontrun{
#' # create a test dataset
#' test_data <- simulate_efast()
#' fit_efast <- efast_hemi(simdat, M = 4, 1:17, 18:34)
#' summary(fit_efast)
#' }
#'
#' @importFrom lavaan lavaan summary
#'
#' @export
efast_hemi <- function(data, M, lh_idx, rh_idx, roi_names, constrain = FALSE,
                       auto.fix.first = FALSE, auto.var = TRUE, auto.efa = TRUE,
                       information = "observed", std.ov = TRUE, ...) {
  if (!is.data.frame(data))
    stop("data argument should be a data.frame")
  if (!is.numeric(M) || M < 2)
    stop("We need at least 2 factors (M) for EFA.")
  if (!is.numeric(lh_idx) || !is.numeric(rh_idx))
    stop("lh.idx and rh.idx should be numeric vectors.")
  if (length(lh_idx) != length(rh_idx))
    stop("RH and LH should have an equal number of variables.")
  if (length(intersect(lh_idx, rh_idx)) != 0)
    stop("There should be no overlap between lh.idx and rh.idx.")

  if (missing(roi_names)) {
    roi_names <- paste0("ROI", seq_along(lh_idx))
  } else {
    if (length(roi_names) != length(lh_idx))
      stop("roi_names should be the same length as idx.")
  }

  df           <- data[, c(lh_idx, rh_idx)]
  colnames(df) <- c(paste0("lh_", roi_names), paste0("rh_", roi_names))

  efa_mod  <- efa_model(df, M)
  hemi_mod <- hemi_model(roi_names, var_const = constrain)
  lat_idx  <- ifelse(constrain, "", lat_index(roi_names))

  model <- paste(efa_mod, hemi_mod, lat_idx, sep = "\n")

  res <- lavaan(
    model          = model,
    data           = df,
    auto.fix.first = auto.fix.first,
    auto.var       = auto.var,
    auto.efa       = auto.efa,
    information    = information,
    std.ov         = std.ov,
    ...
  )
  res@external$type   <- "EFAST_hemi"
  res@external$syntax <- model
  res
}


#' Estimate an EFA model in lavaan
#'
#' This function estimates efa models in the same way that efast models are
#' estimated: using lavaan.
#'
#' @param data <data.frame> the dataset
#' @param M <numeric> How many factors, minimum 2
#' @param auto.fix.first <bool> see lavaan
#' @param auto.var <bool> see lavaan
#' @param auto.efa <bool> see lavaan
#' @param information <character> see lavaan
#' @param std.ov <bool> see lavaan
#' @param ... other arguments passed to lavaan
#'
#' @details The constrained model constrains the residual covariance to be
#'   equal across the different ROIs.
#'
#' @examples
#' \dontrun{
#' # create a test dataset
#' test_data <- simulate_efast()
#' fit_efa <- efa_esem(simdat, M = 4)
#' summary(fit_efa)
#' }
#'
#' @importFrom lavaan lavaan
#'
#' @export
efa_esem <- function(data, M, auto.fix.first = FALSE, auto.var = TRUE,
                     auto.efa = TRUE, information = "observed", std.ov = TRUE,
                     ...) {
  if (!is.data.frame(data))
    stop("data argument should be a data.frame")
  if (!is.numeric(M) || M < 2)
    stop("We need at least 2 factors (M) for EFA.")
  efa_mod  <- efa_model(data, M)
  res <- lavaan(
    model          = efa_mod,
    data           = data,
    auto.fix.first = auto.fix.first,
    auto.var       = auto.var,
    auto.efa       = auto.efa,
    information    = information,
    std.ov         = std.ov,
    ...
  )
  res@external$type   <- "EFA"
  res@external$syntax <- efa_mod
  res
}

#' model_syntax
#' @rdname model_syntax
#' @keywords internal
efa_model <- function(dat, M = 3) {
  # input: dataset
  # output: efa model with M factors
  nm <- colnames(dat)
  mod <- paste(
    "# Exploratory factor model\n# ------------------------\n",
    paste0("efa('efa1')*F", 1:M, collapse = " +\n"),
    "=~",
    paste(nm, collapse = " + ")
  )
  # put relaxed constraints on the residual variances
  mod <- paste0(mod, "\n\n# Uniqueness constraints\n# ----------------------\n")
  for (p in 1:ncol(dat)) {
    var_p <- stats::var(dat[,p])
    nm_p  <- nm[p]
    mod <- paste0(mod, nm_p, " ~~ v_", p, "*", nm_p, "\n",
                  "v_", p, " > ", -var_p * 0.05, "\n",
                  "v_", p, " < ",  var_p * 1.05, "\n")
  }
  return(mod)
}

#' hemi_model
#' @rdname model_syntax
#' @keywords internal
hemi_model  <- function(rois, var_const = FALSE) {
  # input: roi names
  # output: lavaan structure model
  mod <- "# Hemisphere model\n# ----------------\n"
  for (roi in rois) {
    mod <- paste0(mod, roi, " =~ 1*lh_", roi, " + 1*rh_", roi, "\n")
  }
  mod <- paste0(mod, "\n# Hemisphere covariances\n# ----------------------\n")
  for (roi in rois) {
    if (var_const) {
      mod <- paste0(mod, roi, " ~~ cov_lat*", roi, "\n")
      mod <- paste0(mod, "cov_lat > 0\n")
    } else {
      mod <- paste0(mod, roi, " ~~ cov_", roi, "*", roi, "\n")
      mod <- paste0(mod, "cov_", roi, " > 0\n")
    }
  }

  mod
}

#' lat_index stuff
#' @rdname model_syntax
#' @keywords internal
lat_index <- function(rois) {
  # input: roi names
  # output: lateralisation index per region code
  nroi <- length(rois)

  mod <- "# Lateralization index code\n# -----------------------------------\n"

  for (i in 1:nroi) {
    mod <- paste0(mod, "usum_", rois[i], " := v_", i, " + v_", i + nroi, "\n")
    mod <- paste0(mod, "LI_", rois[i], " := usum_", rois[i],
                  " / (usum_", rois[i], " + 2*cov_", rois[i], ")\n")
  }
  mod
}
