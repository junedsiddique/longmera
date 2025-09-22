#' Calculate margins and their differences from longitudinal mixed-effects regression models
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
mix_means <- function(object, tx.contrast0=NULL, tx.contrast1=NULL,
                      ctrl.contrast0=NULL, ctrl.contrast1=NULL) {

  # Error handling
  if (is.null(tx.contrast0) & is.null(tx.contrast1)
      & is.null(ctrl.contrast0) & is.null(ctrl.contrast1)){
    stop("At least 1 contrast vector must be specified.")
  }

  # Extract link function
  model.link <- object$family$link

  if (model.link != "logit" & model.link != "log"){
    stop("Currently, only logit and log links are supported")
  }

  # Obtain marginalized coefficients
  marginals <- GLMMadaptive::marginal_coefs(object, std_errors = TRUE, cores = 1L)

  # Extract marginalized coefficients and their covariance
  betas     <- marginals$betas
  var_betas <- marginals$var_betas

  # Extract the data used in the model
  # Note that the order of the variables in this dataset
  # is the same order as the variables in the formula
  # statement. VERY CONVENIENT!!
  X <- GLMMadaptive::model.matrix(object)

  # Calculate group effects based on contrasts

  # Treatment contrasts
  # Calculate tx change if at least one tx contrast
  if (!is.null(tx.contrast0) | !is.null(tx.contrast1)){
    tx.output <- group_means(tx.contrast0, tx.contrast1,
                             X=X, betas=betas, var_betas=var_betas, link=model.link)
    tx.result <- data.frame()
  }

  # Calculate test stats and label tx results
  if(sum(!is.null(tx.contrast1), !is.null(tx.contrast0))==2){
    tx.result[1,1:4] <- calculate_z_statistics(tx.output$c1$mean.prob,
                                               tx.output$c1$se)

    tx.result[2,1:4] <- calculate_z_statistics(tx.output$c0$mean.prob,
                                               tx.output$c0$se)

    tx.result[3,1:4] <- calculate_z_statistics(tx.output$diff.prob,
                                               tx.output$diff.se)
    row.names(tx.result) <- c("Tx Contrast 1", "Tx Contrast 0", "Tx Change")

  } else if (!is.null(tx.contrast0)) {
    tx.result <- calculate_z_statistics(tx.output$c0$mean.prob,
                                        tx.output$c0$se)
    row.names(tx.result) <- c("Tx Contrast 0")
  } else if (!is.null(tx.contrast1)) {
    tx.result <- calculate_z_statistics(tx.output$c1$mean.prob,
                                        tx.output$c1$se)
    row.names(tx.result) <- c("Tx Contrast 1")
  } else {
    tx.result <- NULL
  }

  # Control contrasts
  # Calculate ctrl change if at least one ctrl contrast
  if (!is.null(ctrl.contrast0) | !is.null(ctrl.contrast1)){
    ctrl.output <- group_means(ctrl.contrast0, ctrl.contrast1,
                               X=X, betas=betas, var_betas=var_betas, link=model.link)
    ctrl.result <- data.frame()
  }

  # calculate test stats and label control results
  if(sum(!is.null(ctrl.contrast1), !is.null(ctrl.contrast0))==2){
    ctrl.result[1,1:4] <- calculate_z_statistics(ctrl.output$c1$mean.prob,
                                                 ctrl.output$c1$se)

    ctrl.result[2,1:4] <- calculate_z_statistics(ctrl.output$c0$mean.prob,
                                                 ctrl.output$c0$se)

    ctrl.result[3,1:4] <- calculate_z_statistics(ctrl.output$diff.prob,
                                                 ctrl.output$diff.se)
    row.names(ctrl.result) <- c("Ctrl Contrast 1", "Ctrl Contrast 0", "Ctrl Change")
  } else if (!is.null(ctrl.contrast0)) {
    ctrl.result <- calculate_z_statistics(ctrl.output$c0$mean.prob,
                                          ctrl.output$c0$se)
    row.names(ctrl.result) <- c("Ctrl Contrast 0")
  } else if (!is.null(ctrl.contrast1)) {
    ctrl.result <- calculate_z_statistics(ctrl.output$c1$mean.prob,
                                          ctrl.output$c1$se)
    row.names(ctrl.result) <- c("Ctrl Contrast 1")
  } else {
    ctrl.result <- NULL
  }

  # Calculate tx minus control at same time point
  if (!is.null(tx.contrast0) & !is.null(ctrl.contrast0)){
    base.output <- group_means(ctrl.contrast0, tx.contrast0,
                               X=X, betas=betas, var_betas=var_betas, link=model.link)
    base.result <- data.frame()

    base.result[1,1:4] <- calculate_z_statistics(base.output$diff.prob,
                                                 base.output$diff.se)
    row.names(base.result) <- c("Time 0: Tx - Ctrl")
  } else {
    base.result <- NULL
  }

  # Calculate tx minus control at same time point
  if (!is.null(tx.contrast1) & !is.null(ctrl.contrast1)){
    time1.output <- group_means(ctrl.contrast1, tx.contrast1,
                                X=X, betas=betas, var_betas=var_betas, link=model.link)
    time1.result <- data.frame()

    time1.result[1,1:4] <- calculate_z_statistics(time1.output$diff.prob,
                                                  time1.output$diff.se)
    row.names(time1.result) <- c("Time 1: Tx - Ctrl")
  } else {
    time1.result <- NULL
  }


  # Calculate difference in differences when
  # Two contrasts specified for each group
  if (sum(!is.null(tx.contrast0), !is.null(tx.contrast1),
          !is.null(ctrl.contrast0),!is.null(ctrl.contrast1))==4) {

    diff.prob <- tx.output$diff.prob - ctrl.output$diff.prob
    diff.j <- tx.output$diff.j - ctrl.output$diff.j # Jacobian

    # Delta method
    diff.se <- sqrt(diff.j%*%var_betas%*%t(diff.j))
    diff <- calculate_z_statistics(diff.prob, diff.se)
    row.names(diff) <- c("Diff in Change")
  } else {
    diff <- NULL
  }

  output <- rbind(tx.result, ctrl.result, base.result, time1.result, diff)

  output
}
