#' Null models in functional space
#'
#' Comparing the amount of occupied functional space against null models
#'
#'@details
#'
#'\code{funspaceNull} The function tests for the statistical difference between the size (functional richness) of the considered TPD, obtained using the \code{funspace} function, against a vector of functional richness values generated using null models (see below) across a user-defined number of iterations. Two null models are currently available for testing. One generates data with a multivariate normal distribution, creating a dataset with normally distributed variables having the same mean and covariance than the observations used to build the functional space (see Carmona et al. 2021). This null model returns a theoretical TPD where some trait combinations (those around the mean of the trait space axes, thus towards the center of the null trait space) are more likely than others (i.e., this null model resembles an ellipse). The other null model generates a dataset with variables following a uniform distribution (see null model 1 in Diaz et al. 2016), creating a distribution where all trait combinations within the range of the original observations are equally possible (i.e., the approximate shape of this null model is a rectangle).
#' Note that the function does not work for funspace objects that are based on a TPDs object created using the package \code{TPD}
#'
#'@references
#'
#'CP Carmona, et al. (2021). Fine-root traits in the global spectrum of plant form and function. Nature 597, 683–687
#'S Diaz, et al. (2016). The global spectrum of plant form and function. Nature 529, 167–171
#'
#'@param funspace An object of class \code{funspace}
#'@param nrep \code{numeric}The number of generated null surfaces
#'@param alter \code{character}. The hypothesis to be tested when comparing the observed trait space against the null model. Options are 'greater', 'less', and 'two-sided'. See specification of the \code{as.randtest} function in the \code{ade4} R package
#'@param null.distribution \code{character}. Data distribution for null model building. Available options are 'multnorm' and 'uniform' to generate data with a multivariate normal or uniform distribution, respectively.
#'@param verbose \code{logical}. Do you want to  information about the progress of the null model to be written to the console?

#'
#'@return \code{funspaceNull} The function returns the list containing all the simulated datasets, the area of the observed trait space, the mean value of the area for the null model (calculated across iterations), the p-value of the difference between observed and simulated trait space, as well as a standardized effect size of the difference between observed trait space and mean null model areas. This output is reported together with the output of  \code{funspace}.
#'
#'@examples
#'
#'# 1. PCA space, multivariate model (see Carmona et al. 2021, Nature)
#'x <- princomp(GSPFF)
#'funtest <- funspace(x = x, PCs = c(1, 2), threshold = 0.95)
#'\donttest{funtestNull <- funspaceNull(funtest, null.distribution = 'multnorm', nrep = 1000)}
#'\donttest{summary(funtestNull)}
#'
#'#'# 2. Two raw traits and uniform distribution (see Diaz et al. 2016, Nature)
#'x <- GSPFF[, c("ph", "sla")]
#'funtest <- funspace(x = x, threshold = 0.95)
#'\donttest{funtestNull <- funspaceNull(funtest, null.distribution = 'uniform', nrep = 1000)}
#'\donttest{summary(funtestNull)}
#'
#'@import ade4
#' @export
funspaceNull <- function(funspace, nrep = 100,
                         alter = 'greater', null.distribution = 'multnorm',
                         verbose = TRUE
                         ){
  #1. Checkings.
  # 1.1 funspace.object must be a pca object of class "funspace"
  if (!inherits(funspace, "funspace")){
    stop("'funspace' must be an object of class 'funspace'")
  }
  if (!(null.distribution %in% c("multnorm", "uniform"))){
    stop("'null.distribution' must be either 'multnorm' or 'uniform'")
  }
  if (funspace$parameters$objectClass == "TPDs"){
    stop("'funspaceNull' does not support funspace objects based on a TPDs function")
  }
  n_divisions <- sqrt(nrow(funspace$parameters$evaluation_grid))

  Null_results <- list()

  # Getting FRich from funspace object
  target.FR <- funspace$FD$global$FRich

  # Setting same mean and cov as original data
  mean.aux <- base::colMeans(funspace$parameters$pcs[, c(1,2)])
  cov.aux  <- stats::cov(funspace$parameters$pcs[, c(1,2)])

  null.FRic = c()
  data.aux.list = list()

  lengthCellX <- unique(funspace$parameters$evaluation_grid[, 1])[2] - unique(funspace$parameters$evaluation_grid[, 1])[1]
  lengthCellY <- unique(funspace$parameters$evaluation_grid[, 2])[2] - unique(funspace$parameters$evaluation_grid[, 2])[1]
  volCell <- lengthCellX * lengthCellY

  for(i in 1:nrep){

  if(verbose){
    cat(paste0("\r REPETITION: ",i,"/",nrep,"\r"))
  }

  # Target distribution (multivariate normal and uniform distribution)

  if(null.distribution != 'multnorm'){

    data.aux <- data.frame(PCx = stats::runif(nrow(funspace$parameters$pcs),
                                       min = min(funspace$parameters$pcs[,1]),
                                       max = max(funspace$parameters$pcs[,1])),
                           PCy = stats::runif(nrow(funspace$parameters$pcs),
                                       min = min(funspace$parameters$pcs[,2]),
                                       max = max(funspace$parameters$pcs[,2])))
  } else {

  data.aux <- data.frame(MASS::mvrnorm(nrow(funspace$parameters$pcs),
                                 mu = mean.aux,
                                 Sigma = cov.aux))
  colnames(data.aux) <- c("PCx", "PCy")
  }

  # kde estimation in the null grid
  est.aux <- ks::kde(x = data.aux,
               H = funspace$parameters$bandwidth,
               eval.points = funspace$parameters$evaluation_grid)

  # Re-escaling of kde to integrate to 1:
  est.aux$estimate <- est.aux$estimate /sum(est.aux$estimate)

  # Retrieving TPD quantiles using the same threshold as the input funspace object
  TPD.quants.aux <- imageTPD(est.aux, thresholdPlot = funspace$parameters$threshold)

  null.nCells <- sum(TPD.quants.aux > 0, na.rm = TRUE)
  null.FRic[i] <- null.nCells * volCell

  data.aux.list[[i]] <- data.aux

}

# Output

  test.diff <- ade4::as.randtest(null.FRic, target.FR, alter = alter)

  Null_results$dataNull <- data.aux.list
  Null_results$Distr.Test <- null.distribution
  Null_results$ObsFric <- test.diff$obs
  Null_results$NullFric <- test.diff$expvar[2]
  Null_results$p_value <- test.diff$pvalue
  Null_results$StdEffSize <- (test.diff$obs - mean(null.FRic)) / stats::sd(null.FRic)
  Null_results$parameters$type <- "null"

  class(Null_results) <- "funspace"
  return(Null_results)
}
