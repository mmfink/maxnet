% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maxnet.R, R/maxnet.default.regularization.R,
%   R/maxnet.formula.R
\name{maxnet}
\alias{maxnet}
\alias{maxnet.default.regularization}
\alias{maxnet.formula}
\title{Maxent over glmnet}
\usage{
maxnet(
  p,
  data,
  f = maxnet.formula(p, data),
  regmult = 1,
  regfun = maxnet.default.regularization,
  addsamplestobackground = T,
  ...
)

maxnet.default.regularization(p, m)

maxnet.formula(p, data, classes = "default")
}
\arguments{
\item{p}{numeric, a vector of 1 (for presence) or 0 (for background)}

\item{data}{a matrix or data.frame of predictor variables}

\item{f}{formula, determines the features to be used}

\item{regmult}{numeric, a constant to adjust regularization}

\item{regfun}{function, computes regularization constant for each feature
By default it uses the \code{\link[=maxnet.default.regularization]{maxnet.default.regularization()}} function to calculate
the base regularization constant for each feature. The output vector is then multiplied
by \code{regmult} to create the vector of \code{penalty.factor}s to pass on to the \code{\link[=glmnet]{glmnet()}} call.}

\item{addsamplestobackground}{logical, if TRUE then add to the background any
presence sample that is not already there}

\item{...}{not used}

\item{m}{a matrix of feature values}

\item{classes}{charcater, feature classes, with \code{l} = "linear", \code{q} = "quadratic",
\code{p} = "product", \code{h} = "hinge", and \code{t} = "threshold". Features are simple
mathematical transformations of the input continuous predictors. Categorical predictors
are used as-is. The default classes used are based on sample size of presence points (\code{np})
as follows:
\itemize{
\item \verb{(np < 10) classes <- "l"}
\item \verb{(np < 15) classes <- "lq"}
\item \verb{(np < 80) classes <- "lqh"}
\item \verb{else classes <- "lqph"}

\code{\link[=maxnet.formula]{maxnet.formula()}} is used to calculate the above default classes. To override this,
either designate the classes to use with \code{classes} and set \code{f = NULL}, or set
\code{classes = NULL} and set \code{f} to your own formula.
}}
}
\value{
returns an object of class \code{maxnet}, which is a list
consisting of a glmnet model with the following elements added:
\itemize{
\item \code{betas} nonzero coefficients of the fitted model
\item \code{alpha} constant offset making the exponential model sum to one over the background data
\item \code{entropy} entropy of the exponential model
\item \code{penalty.factor} the regularization constants used for each feature
\item \code{featuremins} minimum of each feature, to be used for clamping
\item \code{featuremaxs} maximum of each feature, to be used for clamping
\item \code{varmin} minimum of each predictor, to be used for clamping
\item \code{varmax} maximum of each predictor, to be used for clamping
\item \code{samplemeans} mean of each predictor over samples (majority for factors)
\item \code{levels} levels of each predictor that is a factor
}
}
\description{
Maxent species distribution modeling using glmnet for model fitting
}
\details{
Using \code{lp} for the linear predictor and \code{entropy} for the entropy
of the exponential model over the background data, the values plotted on
the y-axis of the response curves are:
\itemize{
\item \code{lp} if \code{type} is "link"
\item \code{exp(lp)} if \code{type} is "exponential"
\item \code{1-exp(-exp(entropy+lp))} if \code{type} is "cloglog"
\item \code{1/(1+exp(-entropy-lp))} if \code{type} is "logistic"
}
}
\examples{
```
  library(maxnet)
  data(bradypus)
  p <- bradypus$presence
  data <- bradypus[,-1]
  mod <- maxnet(p, data)
  plot(mod, type="cloglog")
  mod <- maxnet(p, data, maxnet.formula(p, data, classes="lq"))
  plot(mod, "tmp6190_ann")
```
}
\author{
Steve Phillips
}
