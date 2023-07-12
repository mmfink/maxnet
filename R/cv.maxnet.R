#' Cross-validation for maxnet
#' 
#' Runs a k-fold cross-validation set for maxnet and plots the resulting ROC curve
#' 
#' @export
#' @param p numeric, a vector of 1 (for presence) or 0 (for background)
#' @param data a matrix or data.frame of predictor variables
#' @param f formula, determines the features to be used
#' @param regmult numeric, a constant to adjust regularization
#' @param regfun function, computes regularization constant for each feature
#' By default it uses the [maxnet.default.regularization()] function to calculate
#' the base regularization constant for each feature. The output vector is then multiplied 
#' by `regmult` to create the vector of `penalty.factor`s to pass on to the [glmnet()] call.
#' @param addsamplestobackground logical, if TRUE then add to the background any
#'   presence sample that is not already there
#' @param nfolds integer, number of folds - default is 10. Allowable range of values
#'   is 3 to length(p), however, a large number of folds is not recommended for large datasets.
#' @return returns an object of class \code{cv.maxnet}, which is a list
#'   consisting of a cv.glmnet model with the following elements added:
#'  * `betas` nonzero coefficients of the fitted model 
#'  * `alpha` constant offset making the exponential model sum to one over the background data 
#'  * `entropy` entropy of the exponential model 
#'  * `penalty.factor` the regularization constants used for each feature 
#'  * `featuremins` minimum of each feature, to be used for clamping 
#'  * `featuremaxs` maximum of each feature, to be used for clamping 
#'  * `varmin` minimum of each predictor, to be used for clamping 
#'  * `varmax` maximum of each predictor, to be used for clamping 
#'  * `samplemeans` mean of each predictor over samples (majority for factors) 
#'  * `levels` levels of each predictor that is a factor
#'  * `cmx` confusion matrix
#'  * `AUC` Area under the Receiver Operator Characteristic (ROC) curve
#'  * `PCC` Percent Correctly Classified
#'  
#' @seealso [cv.glmnet()]

cv.maxnet <-
  function(p, data, f=maxnet.formula(p, data), regmult=1.0, 
           regfun=maxnet.default.regularization, addsamplestobackground=T,
           nfolds=10)
  {
    if (anyNA(data)) stop("NA values in data table. Please remove them and rerun.")
    if (!is.vector(p)) {
      p <- as.vector(p)
    }
    if (addsamplestobackground) {
      pdata <- data[p==1, , drop = FALSE]
      ndata <- data[p==0, , drop = FALSE]
      # add to the background any presence data that isn't already in the background
      # TODO: investigate using future.apply here (upstream issue #18)
      toadd <- dplyr::setdiff(pdata, ndata)
      p <- c(p, rep(0, nrow(toadd)))
      data <- rbind(data, toadd)
    }   
    mm <- model.matrix(f, data)
    reg <- regfun(p,mm) * regmult
    weights <- p+(1-p)*100
    ny <- as.factor(p)
    glmnet::glmnet.control(pmin=1.0e-8, fdev=0)  
    model <- glmnet::cv.glmnet(x=mm, y=ny, family="binomial", 
                               lambda=10^(seq(4,0,length.out=200))*sum(reg)/length(reg)*sum(p)/sum(weights), 
                               weights=weights,
                               type.measure="deviance",
                               nfolds=nfolds,
                               keep=TRUE,
                               standardize=F, 
                               penalty.factor=reg)
    class(model) <- c("cv.maxnet", class(model))
    if (length(model$lambda) < 200) {
      msg <- "Error: glmnet failed to complete regularization path.  Model may be infeasible."
      if (!addsamplestobackground) 
        msg <- paste(msg, " Try re-running with addsamplestobackground=T.")
      stop(msg)
    }
    id1se = match(model$lambda.1se, model$lambda)
    bb <- model$glmnet.fit$beta[,id1se]
    model$betas <- bb[bb!=0]
    model$alpha <- 0
    rr <- predict.maxnet(model, data[p==0, , drop = FALSE], type="exponent", clamp=F)
    raw <- rr / sum(rr)
    model$entropy <- -sum(raw * log(raw))
    model$alpha <- -log(sum(rr))
    model$penalty.factor <- reg
    model$featuremins <- apply(mm, 2, min)
    model$featuremaxs <- apply(mm, 2, max)
    vv <- (sapply(data, class)!="factor")
    model$varmin <- apply(data[,vv, drop = FALSE], 2, min)
    model$varmax <- apply(data[,vv, drop = FALSE], 2, max)
    means <- apply(data[p==1,vv, drop = FALSE], 2, mean)
    majorities <- sapply(names(data)[!vv], 
                         function(n) which.max(table(data[p==1,n, drop = FALSE])))
    names(majorities) <- names(data)[!vv]
    model$samplemeans <- unlist(c(means, majorities))
    model$levels <- lapply(data, levels)
    
    model$cmx <- glmnet::confusion.glmnet(model$fit.preval, 
                                          newy=ny, 
                                          family="binomial")[[id1se]]
    metrics <- glmnet::assess.glmnet(model$fit.preval, newy=ny, family="binomial")
    model$AUC <- metrics$auc[[id1se]]
    model$PCC <- 1 - metrics$class[[id1se]]
    roc <- glmnet::roc.glmnet(model$fit.preval, newy=ny)
    
    plot(roc[[id1se]], xlab="False Positive Rate", ylab="True Positive Rate",
         main="ROC Curve", xlim=c(0,1),
         type="l", col="red", asp=1)
    lines(x=c(0,1), y=c(0,1), lty=3, col="lightgray")
    invisible(NULL)
    
    model
  }