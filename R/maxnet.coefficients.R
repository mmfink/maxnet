#' @importFrom utils installed.packages
#' @export
maxnet.coefficients <- 
  function(envnames, object) {
    if('data.table' %in% rownames(installed.packages())) {
      #' @import data.table
      #`%nchin%` <- Negate(data.table::`%chin%`)
      outable <- data.table::data.table(input=envnames)
      xtrafeatures <- setdiff(names(object$betas), outable$input)
      outable <- rbind(outable, xtrafeatures, use.names=FALSE)
      #xtrafeatures <- names(object$betas) %nchin% outable[, input]
      #outable <- rbind(outable, names(object$betas[xtrafeatures]),
      #                 use.names=FALSE)
      data.table::setkey(outable, input)
      #xtra <- which(outable[, input] %nchin% names(object$betas))
      tempdf1 <- data.table::SJ(input = names(object$betas), coeff = object$betas)
      outable <- tempdf1[outable, on="input"]
      outable <- outable[is.na(coeff), coeff := 0]
      tempdf2 <- data.table::SJ(input = names(object$featuremins), 
                               minval = object$featuremins, 
                               maxval = object$featuremaxs, 
                               penalty = object$penalty.factor)
      outable <- tempdf2[outable, on="input"][, abscoef := abs(coeff)]
      # fmn <- which(names(object$featuremins) %in% outable[, input])
      # outable[, minval := object$featuremins[fmn]]
      # outable[, maxval := object$featuremaxs[fmn]]
      # #outable[, devexp := object$dev.ratio[fmn]] * 100
      # outable[, penalty := object$penalty.factor[fmn]]
      # outable[, abscoef := abs(coeff)]
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
