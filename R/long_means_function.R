#' Calculate margins and their differences from longitudinal mixed-effects or GEE
#' regression models
#'
#' @param object An object of type "MixMod" that is obtained from the
#' GLMMadaptive package or an object of type "geeglm" that is obtained
#' from the gee package.
#' @param tx.contrast0 A vector specifying fixed values of covariates and variables
#' to be averaged over (indicated by NA) when calculating a margin for treatment group at baseline.
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
#' # Fit a mixed-effects model using the GLMMadaptive backage
#' library(GLMMadaptive)
#' fm1 <- mixed_model(fixed = quit ~ group + time + grouptime + race + tv + manual,
#'                    random = ~ time | id, data = gruder, family = binomial())
#'
#' # Specify contrasts for means
#' tx.contrast0=c(1, 1, 0, 0, NA, NA, NA) # tx time
#' tx.contrast1=c(1, 1, 4, 4, NA, NA, NA) # tx time 4
#'
#' ctrl.contrast0=c(1, 0, 0, 0, NA, NA, NA) #ctrl time 0
#' ctrl.contrast1=c(1, 0, 4, 0, NA, NA, NA) #ctrl time 4
#'
#' # Enter model object and contrasts into the long_means function
#' output <- long_means(fm1, tx.contrast0=tx.contrast0, tx.contrast1=tx.contrast1,
#'                     ctrl.contrast0=ctrl.contrast0, ctrl.contrast1=ctrl.contrast1)
#'
#' # Using a GEE model with the same contrasts
#' if (requireNamespace("geepack", quietly = TRUE)) {
#' library(geepack)
#' mod1 <- geeglm(quit ~ group + time + grouptime + race + tv + manual, id=id,
#'                       data=gruder, corstr="unstructured", family=binomial())
#'
#' output <- long_means(fm1, tx.contrast0=tx.contrast0, tx.contrast1=tx.contrast1,
#'                     ctrl.contrast0=ctrl.contrast0, ctrl.contrast1=ctrl.contrast1)
#' output
#' }

long_means <- function(object, tx.contrast0=NULL, tx.contrast1=NULL,
                      ctrl.contrast0=NULL, ctrl.contrast1=NULL) {

  # Error handling
  if ("MixMod" %in% class(object)){
    output <- mix_means(object=object, tx.contrast0=tx.contrast0, tx.contrast1=tx.contrast1,
                        ctrl.contrast0=ctrl.contrast0, ctrl.contrast1=ctrl.contrast1)

  } else if ("geeglm" %in% class(object)){
    output <- gee_means(object=object, tx.contrast0=tx.contrast0, tx.contrast1=tx.contrast1,
                        ctrl.contrast0=ctrl.contrast0, ctrl.contrast1=ctrl.contrast1)

  } else {
    stop("Currently, only objects of type `MixMod' from the
    GLMMadaptive package and `geeglm' from the geepack package are supported")
  }
  output
}
