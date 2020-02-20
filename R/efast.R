#' Estimate an EFAST model
#'
#' This function estimates efa models with residual covariance.
#'
#' @param data <data.frame> the dataset
#' @param M <numeric> How many factors, minimum 2
#' @param res_struct <list> residual structure (see details)
#' @param auto.fix.first <bool> see lavaan
#' @param auto.var <bool> see lavaan
#' @param auto.efa <bool> see lavaan
#' @param information <character> see lavaan
#' @param std.ov <bool> see lavaan
#' @param ... other arguments passed to lavaan
#'
#' @details
#'
#' @examples
#' \dontrun{
#' # Use a lavaan test dataset
#' test_data <- lavaan::HolzingerSwineford1939[,7:15]
#'
#' # create an EFA model
#' test_efa <- efast(test_data, 3)
#'
#' # create a (simple) residual structure
#' res_struct <- list(
#'   c("x4", "x7"),
#'   c("x5", "x9")
#' )
#'
#' # create an efast model
#' test_efast <- efast(test_data, 3, res_struct)
#'
#' compare the models
#' lavaan::lavTestLRT(test_efa, test_efast)
#' }
#'
#' @importFrom lavaan lavaan summary
#'
#' @export
efast <- function(data, M, rstruct, auto.fix.first = FALSE, auto.var = TRUE,
                  auto.efa = TRUE, information = "observed",
                  std.ov = TRUE, ...) {
  if (!is.data.frame(data))
    stop("data argument should be a data.frame")
  if (!is.numeric(M) || M < 2)
    stop("We need at least 2 factors (M) for EFA.")

  if (missing(rstruct) || is.null(rstruct) || length(rstruct) == 0) {
    cat("No residual structure specified, performing an EFA.\n")
    res <- efast_efa(data, M, auto.fix.first = auto.fix.first,
                     auto.var = auto.var, auto.efa = auto.efa,
                     information = information, std.ov = std.ov, ...)
    return(res)
  }


  # generate the model syntax
  efa_mod <- efa_model(data, M)
  str_mod <- str_model(data, rstruct)
  model   <- paste(efa_mod, str_mod, sep = "\n")

  # fit the model
  res <- lavaan(
    model          = model,
    data           = data,
    auto.fix.first = auto.fix.first,
    auto.var       = auto.var,
    auto.efa       = auto.efa,
    information    = information,
    std.ov         = std.ov,
    ...
  )
  res@external$type   <- "EFAST"
  res@external$syntax <- model
  res
}


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
  res@external$type   <- "EFAST_HEMI"
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
#' fit_efa <- efast_efa(simdat, M = 4)
#' summary(fit_efa)
#' }
#'
#' @importFrom lavaan lavaan
#'
#' @export
efast_efa <- function(data, M, auto.fix.first = FALSE, auto.var = TRUE,
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
  res@external$type   <- "EFAST_EFA"
  res@external$syntax <- efa_mod
  res
}


#' Synthesised volume data for DKT-atlas ROIs.
#'
#' This data has been synthesised on the basis of real-world data from the
#' Cam-CAN cohort. It contains 68 grey matter volume measurements (34 LH and 34
#' RH ROIs) for 647 participants.
#'
#' @name roi_volume
#' @usage data(roi_volume)
#' @docType data
#' @references \url{cam-can.org}
#' @keywords datasets
#' @format A data frame with 647 rows and 68 columns
NULL

