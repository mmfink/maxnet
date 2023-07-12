#' Variable importance for a maxnet model
#'
#' Computes variable importance scores for all variables present in a maxnet model object. Scores are indicative and should be interpreted with caution as variable interaction is not considered in the computation of the scores.
#' @export
#' @param modobj Object of class maxnet
#' @param resp_table Data.frame. Response table. First column is assumed to be 'response'
#' @param numReplicates Integer. Number of permutations performed to compute importance.
#'
#' @details {
#' The method used to compute variable importance follows that used in R packages \pkg{biomod2} and \pkg{ecospat}. First, model predictions are made for each row of the environmental data table (or 'response table') \emph{resp_table}. 
#'
#' For each variable in the maxnet model object, values for the variable are permuted between rows and a model prediction made for each row using the permuted or shuffled data table. The permutation is performed \emph{numReplicates} times for each variable.
#'
#' At each permutation, a Pearson correlation is computed between reference predictions and the predicted values from the shuffled table. The importance score is 1 - correlation coefficient.
#'
#' A vector of mean scores for each variable expressed as a percentage of the sum of all mean scores is returned.
#' }
#'
#' @return A plot and a named vector of percent importance scores for the variables present in the maxnet model object sorted from highest to lowest.
#' @references based on peterbat1/fitMaxnet varImportance.R [https://rdrr.io/github/peterbat1/fitMaxnet/src/R/varImportance.R]

maxnet.varImp <- function(modobj, resp_table, numReplicates = 5) {
  varList <- names(modobj$samplemeans)
  importance <- vector("numeric", length(varList))
  names(importance) <- varList
  envData <- resp_table[, -1]
  fullModelVals <- predict(modobj, envData)
  
  for(thisVar in varList)
  {
    correls <- vector("numeric", numReplicates)
    origVals <- envData[[thisVar]]
    
    for (thisRep in 1:numReplicates)
    {
      permInd <- sample(1:nrow(envData), nrow(envData))
      envData[[thisVar]] <- origVals[permInd]
      newModelVals <- predict(modobj, envData)
      correls[thisRep] <- cor(fullModelVals, newModelVals)
    }
    
    # Re-instate original values for the variable
    envData[, thisVar] <- origVals
    
    # Compute the mean correlation for this variable
    importance[thisVar] <- mean(correls)
  }
  
  importance <- 1 - importance
  out <- round(100*importance/sum(importance), 2)
  
  par(mar=c(4,7,4,1) + 0.1)
  plotvar <- out[order(out)]
  b <- barplot(height=plotvar, width=length(plotvar), 
               col=hcl.colors(length(plotvar), palette="PuBu"),
               yaxt="n", horiz=TRUE, main="Variable Importance")
  mtext(names(plotvar), side=2, line=0.1, las=1, at=b, cex=0.8)
  invisible(NULL)
  
  return(out)
}