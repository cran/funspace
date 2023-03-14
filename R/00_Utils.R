# \code{new_var_loading} Calculates variable loadings on selected principal components.
#
# @param pca.object a pca object created with the \code{princomp} function.
# @param PCs a vector of length 2, indicating which two principal components will be selected.
# @param trait_ranges the range of values to consider for each of the dimensions.


new_var_loading <- function(pca.object, PCs, trait_ranges){

  new.loadings <- t(pca.object$sdev * t(pca.object$loadings[,]))
  pcs <- pca.object$scores[, c(PCs[1], PCs[2])]
  fit <- new.loadings[, c(PCs[1], PCs[2])]

  angle <- length <- numeric()
  for(i in 1:nrow(fit)){
    angle[i] <- atan2(fit[i, 2], fit[i, 1])
    length[i] <- sqrt(fit[i, 2]**2 + fit[i, 1]**2)
  }
  for(i in 1:nrow(fit)){
    meanRange <- mean(unlist(lapply(trait_ranges, abs)))
    fit[i, 1] <- meanRange * cos(angle[i]) * length[i]
    fit[i, 2] <- meanRange * sin(angle[i]) * length[i]
  }
  load.fit <- list(arrows = fit, loadings = new.loadings)
  return(load.fit)
}

# \code{find_nrow_ncol} Returns the number of rows and columns necessary to accommodate all the levels of the grouping variable to be plotted
#
# @param group.vec a vector containing the levels of the variable

find_nrow_ncol <- function(group.vec){
  group.vec <- as.factor(group.vec) #Added this to make it work when groups vector is not a factor
  npanels <- nlevels(group.vec)
  nrows <- floor(sqrt(npanels))
  ncols <- ceiling(npanels / nrows)
  return(c(nrows, ncols))
}


# \code{limits_plot} Function to set the limits of the functional space evaluated.
#
# @param pcs the scores of observations in the selected functional space.
# @param H bandwidth use in the kernel density estimation.
# @param bufferSD size of the buffer around the most extreme estimations, expressed as number of sd of the bandwidth.
# @param trait_ranges the range of values to consider for each of the dimensions.


limits_plot <- function(pcs, H, bufferSD = 4, trait_ranges = NULL){
  if (is.null(trait_ranges)) {
    trait_ranges <- rep (15, 2)
  }
  if (!inherits(trait_ranges, "list")) {
    trait_ranges <- list()
    for (dimens in 1:2) {
      min_aux <- min(pcs[, dimens]) - bufferSD * sqrt(H[dimens, dimens])
      max_aux <- max(pcs[, dimens]) + bufferSD * sqrt(H[dimens, dimens])
      trait_ranges[[dimens]] <- c(min_aux, max_aux)
    }
  }
  return(trait_ranges)
}


# \code{grid_creation} Function to create grid where kde is evaluated.
#
# @param trait_ranges the range of values to consider for each of the dimensions.
# @param n_divisions number of divisions per dimension.


grid_creation <- function(trait_ranges, n_divisions){
  grid_evaluate<-list()
  for (dimens in 1:2){
    grid_evaluate[[dimens]] <- seq(from = trait_ranges[[dimens]][1],
                                   to = trait_ranges[[dimens]][2],
                                   length=n_divisions)
  }
  evaluation_grid <- expand.grid(grid_evaluate)
  return(evaluation_grid)
}


# \code{imageTPD} Function to map a TPD function in the functional space, and express it as quantiles.
#
# @param x a TPD function (estimated with ks::kde).
# @param thresholdPlot probability threshold to apply (quantiles larger than the threshold will be removed).


imageTPD <- function(x, thresholdPlot = 0.99){
  percentile <- rep(NA, length(x$estimate))
  TPDList <- cbind(index = 1:length(x$estimate), prob = x$estimate, percentile)
  orderTPD <- order(TPDList[,"prob"], decreasing = TRUE)
  TPDList <- TPDList[orderTPD,]
  TPDList[,"percentile"] <- cumsum(TPDList[,"prob"])
  TPDList <- TPDList[order(TPDList[,"index"]),]
  imageTPD <- TPDList

  spacePercentiles <- matrix(data=0, nrow = nrow(x$eval.points), ncol = 1)
  trait1Edges <- unique(x$eval.points[,1])
  trait2Edges <- unique(x$eval.points[,2])
  imageMat <- matrix(NA, nrow = length(trait1Edges), ncol = length(trait2Edges),
                     dimnames = list(trait1Edges, trait2Edges))
  percentileSpace <- x$eval.points
  percentileSpace$percentile <- rep(NA, nrow(percentileSpace))
  percentileSpace[, "percentile"] <- imageTPD[,"percentile"]
  for(i in 1:length(trait2Edges)){
    colAux <- subset(percentileSpace, percentileSpace[,2] == trait2Edges[i])
    imageMat[, i] <- colAux$percentile
  }
  imageMat[imageMat > thresholdPlot] <- NA
  return(imageMat)
}



# \code{FRichnessGlobal} Function to estimate functional richness.
#
# @param x a funspace object

FRichnessGlobal <- function(x){
  evaluation_grid <- x$parameters$evaluation_grid
  lengthCellX <- unique(evaluation_grid[, 1])[2] - unique(evaluation_grid[, 1])[1]
  lengthCellY <- unique(evaluation_grid[, 2])[2] - unique(evaluation_grid[, 2])[1]
  volCell <- lengthCellX * lengthCellY
  nCells <- sum(x$global$images$TPD.quantiles > 0, na.rm = TRUE)
  FRic <- nCells * volCell
  return(FRic)
}

# \code{FRichnessGlobal} Function to estimate functional richness for groups of observations.
#
# @param x a funspace object

FRichnessGroups <- function(x){
  FRic <- numeric()
  evaluation_grid <- x$parameters$evaluation_grid
  lengthCellX <- unique(evaluation_grid[, 1])[2] - unique(evaluation_grid[, 1])[1]
  lengthCellY <- unique(evaluation_grid[, 2])[2] - unique(evaluation_grid[, 2])[1]
  volCell <- lengthCellX * lengthCellY
  for(i in 1:length(x$groups)){
    nCells <- sum(x$groups[[i]]$images$TPD.quantiles > 0, na.rm = TRUE)
    FRic[i] <- nCells * volCell
    names(FRic)[i] <- names(x$groups)[i]
  }
  return(FRic)
}

# \code{FRichnessGlobal} Function to estimate functional divergence.
#
# @param x a funspace object

FDivergenceGlobal <- function(x){
  #### 1. Find center of gravity of occupied space
  COG <- numeric()
  colsOcc <- colSums(x$global$images$TPD.quantiles, na.rm = TRUE)
  colMin <- which(colsOcc > 0)[1]
  colMax <- which(colsOcc > 0)[length(which(colsOcc > 0))]
  rowsOcc <- rowSums(x$global$images$TPD.quantiles, na.rm = TRUE)
  rowMin <- which(rowsOcc > 0)[1]
  rowMax <- which(rowsOcc > 0)[length(which(rowsOcc > 0))]
  trimmedSpace <- x$global$images$TPD.quantiles[rowMin:rowMax, colMin:colMax]

  rowsN <- as.numeric(rownames(trimmedSpace))
  rowsN <- (rowsN - min(rowsN)) / (max(rowsN) - min(rowsN))
  colsN <- as.numeric(colnames(trimmedSpace))
  colsN <- (colsN - min(colsN)) / (max(colsN) - min(colsN))

  occupiedCells <- replace(trimmedSpace, trimmedSpace > 0, 1)
  rownames(occupiedCells) <- rowsN
  colnames(occupiedCells) <- colsN
  weightX <- rowMeans(occupiedCells, na.rm = TRUE)
  coordsX <- rowsN
  COG[1] <- stats::weighted.mean(coordsX[!is.na(weightX)], weightX[!is.na(weightX)])
  weightY <- colMeans(occupiedCells, na.rm = TRUE)
  coordsY <- colsN
  COG[2] <- stats::weighted.mean(coordsY[!is.na(weightY)], weightY[!is.na(weightY)])

  # 2. Calculate the distance of each point in the occupied functional space to the  COG:
  pointsOcc <- matrix(NA, nrow = sum(occupiedCells, na.rm = TRUE), ncol = 3,
                      dimnames = list(1:sum(occupiedCells, na.rm = TRUE),
                                      c("PC1", "PC2", "TPD")))
  index <- 1
  for(row in 1:nrow(occupiedCells)){
    if(rowSums(occupiedCells, na.rm = TRUE)[row] > 0){
      occupiedCols <- which(occupiedCells[row, ] > 0)
      for(col in occupiedCols){
        pointsOcc[index, 1:2] <- c(as.numeric(rownames(occupiedCells)[row]),
                                   as.numeric(colnames(occupiedCells)[col]))
        pointsOcc[index, "TPD"] <- 1 - trimmedSpace[row, col]
        index <- index + 1
      }
    }
  }
  pointsOcc[, "TPD"] <- pointsOcc[, "TPD"] / sum(pointsOcc[, "TPD"])

  dist_COG <- function(x, COG) {
    result_aux<-stats::dist(rbind(x, COG))
    return(result_aux)
  }
  COGDist <- apply(pointsOcc[, 1:2], 1, dist_COG, COG)
  pointsOcc <- cbind(pointsOcc, COGDist)


  # 3. Calculate the mean of the COGDist's
  meanCOGDist <- mean(pointsOcc[, "COGDist"])
  # 4. Calculate the sum of the abundance-weighted deviances for distaces
  #   from the COG (AWdistDeviances) and the absolute abundance-weighted
  #   deviances:
  distDeviances <- pointsOcc[, "COGDist"] - meanCOGDist
  AWdistDeviances <- sum(pointsOcc[, "TPD"] * distDeviances)
  absdistDeviances <- abs( pointsOcc[, "COGDist"] - meanCOGDist)
  AWabsdistDeviances <- sum(pointsOcc[, "TPD"] * absdistDeviances)
  #Finally, calculate FDiv:
  FDiv <- (AWdistDeviances + meanCOGDist) / (AWabsdistDeviances +  meanCOGDist)
  return(FDiv)
}

# \code{FRichnessGlobal} Function to estimate functional divergence for groups of observations.
#
# @param x a funspace object

FDivergenceGroups <- function(x){
  FDiv <- numeric()
  for(i in 1:length(x$groups)){

    #### 1. Find center of gravity of occupied space
    COG <- numeric()
    colsOcc <- colSums(x$groups[[i]]$images$TPD.quantiles, na.rm = TRUE)
    colMin <- which(colsOcc > 0)[1]
    colMax <- which(colsOcc > 0)[length(which(colsOcc > 0))]
    rowsOcc <- rowSums(x$groups[[i]]$images$TPD.quantiles, na.rm = TRUE)
    rowMin <- which(rowsOcc > 0)[1]
    rowMax <- which(rowsOcc > 0)[length(which(rowsOcc > 0))]
    trimmedSpace <- x$groups[[i]]$images$TPD.quantiles[rowMin:rowMax, colMin:colMax]

    rowsN <- as.numeric(rownames(trimmedSpace))
    rowsN <- (rowsN - min(rowsN)) / (max(rowsN) - min(rowsN))
    colsN <- as.numeric(colnames(trimmedSpace))
    colsN <- (colsN - min(colsN)) / (max(colsN) - min(colsN))

    occupiedCells <- replace(trimmedSpace, trimmedSpace > 0, 1)
    rownames(occupiedCells) <- rowsN
    colnames(occupiedCells) <- colsN
    weightX <- rowMeans(occupiedCells, na.rm = TRUE)
    coordsX <- rowsN
    COG[1] <- stats::weighted.mean(coordsX[!is.na(weightX)], weightX[!is.na(weightX)])
    weightY <- colMeans(occupiedCells, na.rm = TRUE)
    coordsY <- colsN
    COG[2] <- stats::weighted.mean(coordsY[!is.na(weightY)], weightY[!is.na(weightY)])

    # 2. Calculate the distance of each point in the occupied functional space to the  COG:
    pointsOcc <- matrix(NA, nrow = sum(occupiedCells, na.rm = TRUE), ncol = 3,
                        dimnames = list(1:sum(occupiedCells, na.rm = TRUE),
                                        c("PC1", "PC2", "TPD")))
    index <- 1
    for(row in 1:nrow(occupiedCells)){
      if(rowSums(occupiedCells, na.rm = TRUE)[row] > 0){
        occupiedCols <- which(occupiedCells[row, ] > 0)
        for(col in occupiedCols){
          pointsOcc[index, 1:2] <- c(as.numeric(rownames(occupiedCells)[row]),
                                     as.numeric(colnames(occupiedCells)[col]))
          pointsOcc[index, "TPD"] <- 1 - trimmedSpace[row, col]
          index <- index + 1
        }
      }
    }
    pointsOcc[, "TPD"] <- pointsOcc[, "TPD"] / sum(pointsOcc[, "TPD"])

    dist_COG <- function(x, COG) {
      result_aux<-stats::dist(rbind(x, COG))
      return(result_aux)
    }
    COGDist <- apply(pointsOcc[, 1:2], 1, dist_COG, COG)
    pointsOcc <- cbind(pointsOcc, COGDist)


    # 3. Calculate the mean of the COGDist's
    meanCOGDist <- mean(pointsOcc[, "COGDist"])
    # 4. Calculate the sum of the abundance-weighted deviances for distaces
    #   from the COG (AWdistDeviances) and the absolute abundance-weighted
    #   deviances:
    distDeviances <- pointsOcc[, "COGDist"] - meanCOGDist
    AWdistDeviances <- sum(pointsOcc[, "TPD"] * distDeviances)
    absdistDeviances <- abs( pointsOcc[, "COGDist"] - meanCOGDist)
    AWabsdistDeviances <- sum(pointsOcc[, "TPD"] * absdistDeviances)
    #Finally, calculate FDiv:
    FDiv[i] <- (AWdistDeviances + meanCOGDist) / (AWabsdistDeviances +  meanCOGDist)
    names(FDiv)[i] <- names(x$groups)[i]
  }
  return(FDiv)
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
