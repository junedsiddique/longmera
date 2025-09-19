
#' Calculate means and their difference for two contrasts
#'
#' @param contrast0 A vector specifying fixed values of covariates and variables
#' to be averaged over when (NA) calculating a margin for treatment group at baseline
#' @param contrast1 A contrast for a second time point
#' @param X The matrix of covariates used in the model
#' @param betas Marginalized regression coefficients
#' @param var_betas Variance covariance matrix of betas
#' @param link Link function. Extracted from object
#'
#' @return A list of means, SEs, Jacobians, for the two contrasts and their difference
#' @export
#'
group_means <- function(contrast0=NULL, contrast1=NULL,
                            X=X, betas=betas, var_betas=var_betas, link=link) {

  # Error handling
    if (!is.null(contrast0) & length(contrast0) != length(betas)){
      stop("Length of contrast 0 vector should be equal to ",
           length(betas), " which is the length of the coefficient vector.")
    }

  if (!is.null(contrast1) & length(contrast1) != length(betas)){
    stop("Length of contrast 1 vector should be equal to ",
         length(betas), " which is the length of the coefficient vector.")
  }

  # Only first contrast specified
  if (is.null(contrast1) & !is.null(contrast0)){
  c0 <- contrast_mean(contrast0, X=X, betas=betas, var_betas=var_betas, link=link)
  output <- list(c0=c0)

  # Only second contrast specified
  } else if (is.null(contrast0) & !is.null(contrast1)){
    c1 <- contrast_mean(contrast1, X=X, betas=betas, var_betas=var_betas, link=link)

    output <- list(c1=c1)

  # Both contrasts specified
  } else {
  c0 <- contrast_mean(contrast0, X=X, betas=betas, var_betas=var_betas, link=link)
  c1 <- contrast_mean(contrast1, X=X, betas=betas, var_betas=var_betas, link=link)

  diff.prob <- c1$mean.prob - c0$mean.prob
  diff.j <- c1$j - c0$j # Jacobian

  # Delta method
  diff.se <- sqrt(diff.j%*%var_betas%*%t(diff.j))
  diff <- calculate_z_statistics(diff.prob, diff.se)

  output <- list(c0=c0, c1=c1, diff.prob=diff.prob,
                 diff.se=diff.se, diff.j=diff.j)

  }

  output
}


#' Calculate a mean for a single contrast
#'
#' @param contrast A vector specifying fixed values of covariates and variables
#' to be averaged over when (NA) calculating a margin
#' @param X The matrix of covariates used in the model
#' @param betas Marginalized regression coefficients
#' @param var_betas Variance covariance matrix of betas
#' @param link Link function. Extracted from object
#'
#' @return A list containing the mean, its SE, and a vector of Jacobians
#' @export
#'
contrast_mean <- function(contrast, X=X, betas=betas,
                          var_betas=var_betas, link=link) {


  # Create counterfactual dataset based on contrasts
  # the NAs in the contrast statement mean that
  # the values for those variables should not change
  # Use this counterfactual dataset to estimate
  # probabilities and their standard errors

  # Create the counterfactual dataset based on the contrast
  X[, which(!is.na(contrast))] <-
  matrix(rep(contrast[!is.na(contrast)], nrow(X)),
         nrow = nrow(X), byrow = TRUE)

  pred <- X %*% betas              # Linear predictor
  result <- calculate_derivative(link=link, pred=pred) # Derivative
  prob <-  result$prob  # Probability
  deriv <- result$deriv

  j <- (t(deriv)%*% X) / (nrow(X)) # Jacobian

  # Delta method
  se <- sqrt(j%*%var_betas%*%t(j))
  mean.prob <- mean(prob)

  output <- list(mean.prob=mean.prob, se=se, j=j)

  output

}

#' Calculate derivatives
#'
#' @param link Link function. Extracted from object
#' @param pred Linear predictor X%*%betas
#'
#' @returns A vector of derivatives
#' @export
#'
calculate_derivative <- function(link, pred) {

  if (link == "logit") {
    prob <-  1/ (1 + exp(0-(pred)))  # Probability
    deriv <- prob*(1-prob)           # Derivative

  } else if (link == "log") {
    prob <- exp(pred)
    deriv <- exp(pred)
  }
  result <- list(prob=prob, deriv=deriv)
  result
}



#' Calculate z-statistics
#'
#' @param mean Mean value
#' @param se Standard Error
#'
#' @return Vector of mean, SE, z statistic, and p-value
#' @export
#'
#' @keywords internal
calculate_z_statistics <- function(mean, se) {

  z <- mean/se

  pval=2*stats::pnorm(-abs(z))
  nice.pval <- format.pval(pval, eps = 0.001, digits = 3)

  cmat <- data.frame(cbind(round(mean,4), round(se,4), round(z,4), nice.pval))
  colnames(cmat) <- c("Estimate", "Std.Err", "z-value", "p-value")

  cmat
}







