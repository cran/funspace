#' Functional space GAM
#'
#' Mapping response variables in a functional space
#'
#'@details
#'
#'Different response variables can be mapped onto a functional space. In \code{funspace}, we follow the approach by Carmona et al. (2021), in which a generalized additive model is estimated across the bidimensional functional space. The resulting models show the predicted values of the response variable at each position of the portion of the functional space that is defined in the TPD of the global set of observations or of individual groups.
#'
#'@references
#'
#'CP Carmona, et al. (2021). Erosion of global functional diversity across the tree of life. Science Advances eabf2675
#'
#'@param y vector including the variable to be mapped inside the functional space. There must be a correspondence between the elements of y and the observations used to make the PCA (contained in 'pca.object'), both in the number of elements and in their order.
#'@param funspace An object of class \code{funspace} providing the functional space to be considered. See function \code{funspace}
#'@param family A family object specifying the distribution and link to use in the gam model. Defaults to "gaussian". See package \code{mgcv} for more details.
#'@param minObs minimum number of observations needed in a  group to make a model (defaults to 30).

#'@return The function returns an object of class \code{funspace} containing the functional space, trait probability distributions, and the fitted gam models. The \code{funspace} class has specific methods exists for the generic functions \code{plot} and \code{summary}.
#'
#'
#' @examples
#'
#'# 1. GAM on a space based on a PCA
#'x <- princomp(GSPFF)
#'funtest <- funspace(x = x, PCs = c(1, 2), threshold = 0.95)
#'y <- abs(x$scores[, 1] * x$scores[, 2]) + rnorm(nrow(GSPFF), mean = 0, sd = 1)
#'funtestGAM <- funspaceGAM(y = y, funspace = funtest)
#'plot(funtestGAM, quant.plot = TRUE, quant.col = "grey90")
#'summary(funtestGAM)
#'
#'
#' @import mgcv
#' @export

funspaceGAM <- function(y, funspace, family = "gaussian", minObs = 30) {
  #1. Checkings.
  # 1.1 funspace.object must be a pca object of class "funspace"
  if (!inherits(funspace, "funspace")){
    stop("'funspace' must be a pca object of class 'funspace'")
  }
  # length of y same as number of points in trait space
  if (length(y) != nrow(funspace$parameters$pcs)){
    stop(" number of observations differ between 'y' and 'funspace")
  }
  results <- funspace
  # FIT A SINGLE GLOBAL MODEL (for all points, regardless of groups): ### IMPORT FROM PREVIOUS
  #7. GAM model
  PC1 <- results$parameters$pcs[, 1]
  PC2 <- results$parameters$pcs[, 2]
  results$global$gam <- mgcv::gam(y ~ s(PC1, PC2), family = family, method = "REML")
  #8. Predict across the whole trait space:
  predAux <- stats::predict(object = results$global$gam,
                     newdata =  results$parameters$evaluation_grid,
                     type = "response",
                     se.fit = TRUE)
  results$global$predicted <- cbind(results$parameters$evaluation_grid,
                                    predicted = predAux$fit,
                                    se = predAux$se.fit)
  imagePred <- imageSEPred <- results$global$images$TPD.quantiles
  imagePred[ , ] <- imageSEPred[ , ] <- NA
  for(k in 1:length(unique(results$parameters$evaluation_grid[, 2]))){
    colAux <- subset(results$global$predicted,
                     results$global$predicted[,2] == unique(results$parameters$evaluation_grid[, 2])[k])
    imagePred[, k] <- colAux$predicted
    imageSEPred[, k] <- colAux$se
  }
  ### apply mask
  maskQuant <- results$global$images$TPD.quantiles
  maskQuant[!is.na(maskQuant)] <- 1
  imagePredMask <- imagePred * maskQuant
  imageSEPredMask <- imageSEPred * maskQuant

  results$global$images$predicted <- imagePredMask
  results$global$images$SE.predicted <- imageSEPredMask

  ########## IF there are groups, a model for each group:
  if(!is.null(results$parameters$group.vec)){
    group.vec <- results$parameters$group.vec
    df <- data.frame(pcs = results$parameters$pcs,
                     y = y,
                     group.vec = results$parameters$group.vec)
    colnames(df) <- c("PC1", "PC2", "y", "group.vec")
    for (i in 1:nlevels(df$group.vec)) {
      targetGroup <- levels(df$group.vec)[i]
      dfGroup <- subset(df, df$group.vec == targetGroup)
      if(nrow(dfGroup) < minObs){
        warning(c("Group '", targetGroup, "' was filtered out" ))
        results$groups[[targetGroup]]$gam <- list()
        results$groups[[targetGroup]]$gam$exist <- FALSE
      } else{
        # GAM model
          PC1 <- dfGroup$PC1
          PC2 <- dfGroup$PC2
          y_aux <- dfGroup$y
          results$groups[[targetGroup]]$gam <- mgcv::gam(y_aux ~ s(PC1, PC2), family = family, method = "REML")
          results$groups[[targetGroup]]$gam$exist <- TRUE
          #8. Predict across the whole trait space:
          predAux <- stats::predict(object = results$groups[[targetGroup]]$gam,
                                    newdata = results$parameters$evaluation_grid,
                                    type = "response", se.fit = TRUE)
          results$groups[[targetGroup]]$predicted <- cbind(results$parameters$evaluation_grid,
                                                           predicted = predAux$fit,
                                                           se = predAux$se.fit)
          imagePred <- imageSEPred <- results$groups[[targetGroup]]$images$TPD.quantiles
          imagePred[ , ] <- imageSEPred[ , ] <- NA
          for(k in 1:length(unique(results$parameters$evaluation_grid[, 2]))){
            colAux <- subset(results$groups[[targetGroup]]$predicted,
                             results$groups[[targetGroup]]$predicted[,2] ==
                               unique(results$parameters$evaluation_grid[, 2])[k])
            imagePred[, k] <- colAux$predicted
            imageSEPred[, k] <- colAux$se
          }
          ### apply mask
          maskQuant <- results$groups[[targetGroup]]$images$TPD.quantiles
          maskQuant[!is.na(maskQuant)] <- 1
          imagePredMask <- imagePred * maskQuant
          imageSEPredMask <- imageSEPred * maskQuant
          results$groups[[targetGroup]]$images$predicted <- imagePredMask
          results$groups[[targetGroup]]$SE.predicted <- imageSEPredMask
      }
    }
  }
  results$parameters$type <- "gam"
  class(results) <- "funspace"
  return(results)
}

