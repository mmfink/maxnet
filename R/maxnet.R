#' @name maxnet
#' @aliases maxnet.default.regularization maxnet.formula
#' @title Maxent over glmnet
#'
#' @description 
#' Maxent species distribution modeling using glmnet for model fitting
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
#' @param m a matrix of feature values
#' @param classes charcater, feature classes, with `l` = "linear", `q` = "quadratic",
#' `p` = "product", `h` = "hinge", and `t` = "threshold". Features are simple
#' mathematical transformations of the input continuous predictors. Categorical predictors
#' are used as-is. The default classes used are based on sample size of presence points (`np`) 
#' as follows:
#'   
#' * `(np < 10) classes <- "l"`
#' * `(np < 15) classes <- "lq"`
#' * `(np < 80) classes <- "lqh"`
#' * `else classes <- "lqph"`
#' 
#'   [maxnet.formula()] is used to calculate the above default classes. To override this,
#'either designate the classes to use with `classes` and set `f = NULL`, or set
#'`classes = NULL` and set `f` to your own formula.
#' @param ... not used
#'
#' @details
#' Using `lp` for the linear predictor and `entropy` for the entropy
#' of the exponential model over the background data, the values plotted on
#' the y-axis of the response curves are:
#'   
#'   * `lp` if `type` is "link"
#'   * `exp(lp)` if `type` is "exponential"     
#'   * `1-exp(-exp(entropy+lp))` if `type` is "cloglog"
#'   * `1/(1+exp(-entropy-lp))` if `type` is "logistic"
#'
#' @return returns an object of class \code{maxnet}, which is a list
#'   consisting of a glmnet model with the following elements added:
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
#' @author Steve Phillips
#' @examples
#' ```
#'   library(maxnet)
#'   data(bradypus)
#'   p <- bradypus$presence
#'   data <- bradypus[,-1]
#'   mod <- maxnet(p, data)
#'   plot(mod, type="cloglog")
#'   mod <- maxnet(p, data, maxnet.formula(p, data, classes="lq"))
#'   plot(mod, "tmp6190_ann")
#' ```
maxnet <-
function(p, data, f=maxnet.formula(p, data), regmult=1.0, 
         regfun=maxnet.default.regularization, addsamplestobackground=T, ...)
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
   glmnet::glmnet.control(pmin=1.0e-8, fdev=0)  
   model <- glmnet::glmnet(x=mm, y=as.factor(p), family="binomial", 
                           standardize=F, penalty.factor=reg, 
                           lambda=10^(seq(4,0,length.out=200))*sum(reg)/length(reg)*sum(p)/sum(weights), 
                           weights=weights)
   
   class(model) <- c("maxnet", class(model))
   if (length(model$lambda) < 200) {
        msg <- "Error: glmnet failed to complete regularization path.  Model may be infeasible."
        if (!addsamplestobackground) 
           msg <- paste(msg, " Try re-running with addsamplestobackground=T.")
        stop(msg)
   }
   bb <- model$beta[,200]
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
   model
}
