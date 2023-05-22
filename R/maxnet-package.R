#' @name maxnet-package
#' @aliases Maxent
#' @docType package
#' @title maxnet: Maxent over glmnet
#' @description Procedures to fit species distributions models from occurrence records and environmental variables, using 'glmnet' for model fitting. Model structure is the same as for the 'Maxent' Java package, version 3.4.0, with the same feature types and regularization options.  See the 'Maxent' website <http://biodiversityinformatics.amnh.org/open_source/maxent> for more details.
#' @references
#' Phillips, S. and Dudik, M. 2008. Modeling of species distributions with MaxEnt: new extensions and a comprehensive evaluation. Ecography 31:161.
#' Fithian, W. and Hastie, T. 2012. Finite-sample equivalence of several statistical models for presence-only data. <http://arxiv.org/abs/1207.6950v1>.
#' @author Steve Phillips
#' @seealso [url(https://github.com/mrmaxent/maxnet)] for source code, and [glmnet()] for base functionality.
#' @importFrom stats approx formula model.matrix predict sd setNames
#' @keywords internals

