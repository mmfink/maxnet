#' @name plot
#' @title Create response plots for user selected predictors in a maxnet model
#' 
#' @seealso [response.plot()]
#' @export
#' @param x	an object of class maxnet, i.e., a fitted model.
#' @param vars character, vector of predictors for which response plots are desired.
#' @param common.scale logical,	if true, all plots use the same scale on the y-axis.
#' @param type character, type of response to plot on y-axis.
#' @param ylab character, label for y-axis
#' @param plot logical, if TRUE render a plot, if FALSE return a list of 
#'   data frames with variable and response columns
#' @param mar numeric, 4 element value for margins (lines, in order of bottom, left, top, right)
#'   See `par` for details.
#' @param N numeric, the number of intervals over which to sample the response
#' @param ... other arguments passed to [plot()] or [barplot()]
#' @return if `plot` is `FALSE` then return a list of data frames 
#'   that contain variable and response columns otherwise `NULL` invisibly
#' @examples  
#' ```
#' # To generate the same type of plots that the Maxent java program does
#'   library(maxnet)
#'   data(bradypus)
#'   p <- bradypus$presence
#'   data <- bradypus[,-1]
#'   mod <- maxnet(p, data)
#'   plot(mod, type="cloglog", asp=0.5)
#' ```
plot.maxnet <-
function(x, vars=names(x$samplemeans), 
         common.scale=TRUE, 
         type=c("link","exponential","cloglog","logistic"), 
         ylab=NULL, 
         plot = TRUE,
         mar = c(4,3,3,1) + 0.1,
         N = 100, 
         ...)
{
   type <- match.arg(type)
   ylim=NULL
   if (common.scale && (type=="link" || type=="exponential")) {
        vals <- do.call(c, lapply(vars, function(v) 
          response.plot(x, v, type, plot=F)$pred))
        ylim=c(min(vals), max(vals))
   }
   if (type=="cloglog" || type=="logistic") ylim=c(0,1)
   
   if (plot){
     nc <- ceiling(sqrt(length(vars)))
     nr <- ceiling(length(vars)/nc)
     graphics::par(mfrow=c(nr,nc), mar=mar,
                   cex.axis=1.8, cex.lab=2.2, lwd=2)
     for (v in vars) response.plot(x, v, type, ylim=ylim, ylab=ylab, 
                                   N = N, col="red")
   } else {
     return(sapply(vars, 
                   function(v){
                     response.plot(x, v, type, ylim=ylim, ylab=ylab, 
                                   N = N, plot = FALSE)
                   }, simplify = FALSE))
   }
   invisible(NULL)
}
