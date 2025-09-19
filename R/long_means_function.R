#' Calculate margins and their differences from longitudinal mixed-effects or GEE
#' regression models
#'
#' @param object An object of type "MixMod" that is obtained from the 
#' GLMMadaptive package
#' @param tx.contrast0 A vector specifying fixed values of covariates and variables
#' to be averaged over when (NA) calculating a margin for treatment group at baseline.
#' Must be the same length as the coefficient vector from the object model.
#' @param tx.contrast1 A vector to calculate the margin for the treatment group at follow-up.
#' @param ctrl.contrast0 A vector to calculate the margin for the control group at baseline.
#' @param ctrl.contrast1 A vector to calculate the margin for the control group at follow-up.
#'
#' @return A data frame of margins, their standard errors, z-statistics, and p-values.
#' If a contrast is specified for both baseline and follow-up, their difference is
#' also calculated. If all 4 contrasts are specified, the difference in differences
#' is also calculated. 
#' @export
#'
#' @examples
long_means <- function(object, tx.contrast0=NULL, tx.contrast1=NULL,
                      ctrl.contrast0=NULL, ctrl.contrast1=NULL) {
  
  # Error handling
  if (class(object) != "MixMod"){
    stop("Currently, only objects of type `MixMod' from the
    GLMMadaptive package are supported")
  } else {
    output <- mix.means(object=object, tx.contrast0=tx.contrast0, tx.contrast1=tx.contrast1,
                          ctrl.contrast0=ctrl.contrast0, ctrl.contrast1=ctrl.contrast1)
  }
  output
}
  