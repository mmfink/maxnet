#' @importFrom utils installed.packages
#' @export
maxnet.coefficients <- 
  function(envnames, object) {
    if('data.table' %in% rownames(installed.packages())) {
      #' @import data.table
      `%nchin%` <- Negate(data.table::`%chin%`)
      outable <- data.table::data.table(input=envnames)
      xtrafeatures <- names(object$betas) %nchin% outable[, input]
      outable <- rbind(outable, names(object$betas[xtrafeatures]),
                       use.names=FALSE)
      data.table::setkey(outable, input)
      bs <- which(outable[, input] %in% names(object$betas))
      outable[bs, coeff := object$betas][!bs, coeff := 0]
      
      fmn <- which(names(object$featuremins) %in% outable[, input])
      outable[, minval := object$featuremins[fmn]]
      outable[, maxval := object$featuremaxs[fmn]]
      outable[, penalty := object$penalty.factor[fmn]]
      outable[, abscoef := abs(coeff)]
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
      outable$penalty = object$penalty.factor[fmn]
      outable$abscoef = abs(coeff)
    }
    
    return(outable)
  }
