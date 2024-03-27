#' @name maxnet.coefficients
#' @title A table of 'features' and lambdas
#' @author Michelle Fink
#' 
#' @description
#' Similar to what was produced in the java Maxent \code{lambdas} file. Note that at this time, hinge features display min/max values in the name, not in the output \code{minval} and \code{maxval} columns.
#' 
#' @importFrom utils installed.packages
#' 
#' @export
#' @param envnames character vector of model input names. This is required because the maxnet model object only keeps non-zero coefficients, and so may drop original inputs.
#' @param object model object of type \code{maxnet}.
#' 
#' @returns a data.frame or data.table (if \pkg{data.table} installed) with fields \code{input, minval, maxval, penalty, coeff, abscoef}, sorted descending by absolute value of coefficient.
maxnet.coefficients <- 
  function(envnames, object) {
    if('data.table' %in% rownames(installed.packages())) {
      outable <- data.table::data.table(input=envnames)
      xtrafeatures <- setdiff(names(object$betas), outable$input)
      outable <- rbind(outable, xtrafeatures, use.names=FALSE)
      data.table::setkey(outable, input)
      tempdf1 <- data.table::SJ(input = names(object$betas), coeff = object$betas)
      outable <- tempdf1[outable, on="input"]
      outable <- outable[is.na(coeff), coeff := 0]
      tempdf2 <- data.table::SJ(input = names(object$featuremins), 
                                minval = object$featuremins, 
                                maxval = object$featuremaxs, 
                                penalty = object$penalty.factor)
      outable <- tempdf2[outable, on="input"][, abscoef := abs(coeff)]

    } else {
      outable <- data.frame(input=envnames, stringsAsFactors = FALSE)
      xfeatures <- which(names(object$betas) %in% outable$input)
      outable$input <- rbind(outable$input, names(object$betas[!xfeatures]),
                             use.names=FALSE)
      bs <- which(outable$input %in% names(object$betas))
      outable$coeff = rep(0, length(outable))
      outable$coeff[bs] = object$betas[xfeatures]
      fmn <- which(names(object$featuremins) %in% outable$input)
      outable$minval = object$featuremins[fmn]
      outable$maxval = object$featuremaxs[fmn]
      #outable$devexp = object$dev.ratio[fmn] * 100
      outable$penalty = object$penalty.factor[fmn]
      outable$abscoef = abs(coeff)
    }
    
    return(outable)
  }
