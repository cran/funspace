#' Functional space
#'
#' Defines the functional structure of a set of species
#'
#'@details
#'
#'The functional structure of a set of organisms refers to how these organisms are distributed within a functional space (a space defined by traits). Functional structure can be expressed in probabilistic terms using trait probability density functions (TPD). TPD functions reflect how densely the organisms occupy the different parts of the functional space, and are implemented in the package \code{TPD} (Carmona et al. 2019).
#'
#'\code{funspace} allows the user to define functional structure in a two-dimensional functional space created using a PCA, other ordination methods, or raw traits. The function automatically estimates the probability of occurrence of trait combinations within the space using kernel density estimation with unconstrained bandwidth using the functions from the \code{ks} R package  (Duong, 2007). Contour lines can be drawn at any quantile of the probability distribution. Colored areas, corresponding to the target quantiles, visually summarize the probability of occurrence of certain trait combinations within the trait space.
#'
#'@references
#'
#'CP Carmona, F de Bello, NWH Mason, J Leps (2019). Trait Probability Density (TPD): measuring functional diversity across scales based on trait probability density with R. Ecology e02876.
#'T Duong, T., (2007). ks: Kernel Density Estimation and Kernel Discriminant Analysis for Multivariate Data in R. J. Stat. Softw. 21(7), 1-16.
#'
#'@param x Data to create the functional space. It can be either a PCA object obtained using the \code{princomp} function, a PCoA obtained using the \code{capscale} function from the \code{vegan} package, an NMDS generated with the \code{metaMDS} or \code{monoMDS} functions from \code{vegan}, a TPDs object generated with the \code{TPD} package or a matrix or data frame with at least two columns (representing two dimensions which can be either traits or ordination scores obtained with other methods).
#'@param PCs A vector specifying the Principal Components to be considered (e.g. choosing \code{PCs = c(1,2)} would lead to to consider the first and the second principal components). Only applies if \code{x} contains a PCA. Defaults to \code{c(1, 2)}, which selects the first two principal components.
#'@param group.vec An object of class factor specifying the levels of the grouping variable.
#'@param fixed.bw Logical indicating whether the same bandwidth that is used in the kde estimation for the whole dataset should also be used for the kde estimation of individual groups of observations (\code{fixed.bw = T}), or if a different bandwidth has to be estimated for each group ((\code{fixed.bw = F})). Defaults to TRUE, which makes the most extreme quantiles of the individual groups to coincide with those of the global distribution, and allows for more meaningful comparisons of the amount of functional space occupied by groups (functional richness).
#'@param n_divisions The number of equal-length parts in which each principal component  should be divided to calculate the grid in which calculations are based. Higher values of n_divisions will result in larger computation times, but also more smooth graphics. Defaults to 100.
#'@param trait_ranges A list indicating the range of values that will be considered in the calculations for each of the considered PCA components. The list should contain the range (minimum and maximum) of  values that will be considered. Each element of the list corresponds with one PCA component. The order of the components must be the same as the order provided in \code{PCs}. Defaults to	NULL, in which case ranges are automatically calculated to ensure the functional space considered is sufficiently large to encompass the whole TPD function.
#'@param threshold The probability threshold to consider to estimate the TPD function. TPD functions are positive across the whole trait space; \code{threshold} defines boundaries beyond which the TPD function is set to 0 (see Carmona et al. 2016; 2019 for more information). Defaults to 0.999.

#'@return \code{funspace} The function returns an object of class \code{funspace} containing characteristics of the functional space and the trait probability distributions. The object includes estimations of functional richness and functional divergence for all observations taken together (global) and for each individual group (if groups are provided). The \code{funspace} class has specific methods exists for the generic functions \code{plot} and \code{summary}.
#'
#' @examples
#'
#'# 1. Plotting a space based on a PCA
#'x <- princomp(GSPFF)
#'funtest <- funspace(x = x, PCs = c(1, 2), threshold = 0.95)
#'summary(funtest)
#'plot(funtest, type = "global")
#'
#'#2. To include groups, let's consider two major families.
#'# We will use two raw traits, ph and sla:
#'selFam <- c("Pinaceae", "Fabaceae")
#'selRows <- which(GSPFF_tax$family %in% selFam)
#'GSPFF_subset <- GSPFF[selRows, c("ph", "sla")]
#'tax_subset <- droplevels(GSPFF_tax[selRows, ])
#'funtest <- funspace(x = GSPFF_subset, threshold = 0.95, group.vec = tax_subset$family)
#'summary(funtest)
#'plot(funtest, type = "global")
#'plot(funtest, type = "groups", axis.title.x = "Plant height",
#'    axis.title.y = "Specific leaf area",
#'    quant.plot = TRUE, pnt = TRUE, pnt.cex = 0.5,
#'    pnt.col = rgb(0, 1, 1, alpha = 0.2))
#'
#' @import ks
#' @import vegan
#' @export

funspace <- function(x,
                     PCs = c(1, 2),
                     group.vec = NULL,
                     fixed.bw = TRUE,
                     n_divisions = 100,
                     trait_ranges = NULL,
                     threshold = 0.999
                     ){

  args <- list(x = x, PCs = PCs, group.vec = group.vec, n_divisions = n_divisions,
               trait_ranges = trait_ranges, threshold = threshold)
  #1. Creating list to store results
  results <- list()
  #2. Definition of limits and grid for kde estimation
  results$parameters <- list()
  results$originalX <- x
  sType  <- do.call(spaceType, args)
  results$parameters$objectClass <- sType$objectClass
  results$parameters$PCs <- sType$PCs
  results$parameters$pcs <- sType$pcs
  results$parameters$bandwidth <-sType$bandwidth
  results$parameters$threshold <- sType$threshold
  results$parameters$group.vec <- sType$group.vec
  results$parameters$trait_ranges <- sType$trait_ranges
  results$parameters$evaluation_grid <- sType$evaluation_grid

  #### EXCEPT WHEN A TPD OBJECT IS PROVIDED, MAKE A SINGLE GLOBAL PLOT:
  #4. kde estimation
  if(results$parameters$objectClass != "TPDs"){
    est <- ks::kde(x = results$parameters$pcs,
                   H = results$parameters$bandwidth,
                   compute.cont = TRUE,
                   eval.points = results$parameters$evaluation_grid)

    #5. Re-escaling of kde to integrate to 1:
    est$estimate <- est$estimate /sum(est$estimate)
    #6. Creating and plotting image matrix with quantiles and threshold
    results$global <- list() # To store the global TPD
    results$global$images <- list()
    results$global$images$noThreshold.quantiles <- imageTPD(est, thresholdPlot = 1)
    results$global$images$TPD.quantiles <- imageTPD(est, thresholdPlot = results$parameters$threshold)
  }
  if(results$parameters$objectClass == "TPDs"){
    results$global <- list()
    results$global$images <- list()
    results$global$images$noThreshold.quantiles <- NA
    results$global$images$TPD.quantiles <- NA
  }
  #### IF THERE ARE ALSO GROUPS:
  if(!is.null(results$parameters$group.vec)){
    results$parameters$group.vec <- as.factor(results$parameters$group.vec)
    results$groups <- list()
    if(results$parameters$objectClass != "TPDs"){
      pcs.gr <- data.frame(results$parameters$pcs, results$parameters$group.vec)
      colnames(pcs.gr) <- c(colnames(results$parameters$pcs), "group.vec")
      pcs.gr.list <- split(pcs.gr, pcs.gr$group.vec)
      for(i in 1:length(pcs.gr.list)){
        results$groups[[i]] <- list()
        names(results$groups)[i] <- names(pcs.gr.list)[i]
        results$groups[[i]]$parameters <- list()
        results$groups[[i]]$parameters$points <- pcs.gr.list[[i]][,c(1,2)]
        if(fixed.bw == FALSE){
          results$groups[[i]]$parameters$bandwidth <- ks::Hpi(x = results$groups[[i]]$parameters$points)      # optimal bandwidth estimation at group level
        } else{
          results$groups[[i]]$parameters$bandwidth <- results$parameters$bandwidth      # Use global bandwidth
        }
        kde.gr <- ks::kde(x = results$groups[[i]]$parameters$points,
                          H = results$groups[[i]]$parameters$bandwidth,
                          compute.cont=TRUE,
                          eval.points = results$parameters$evaluation_grid)     # kernel density estimation
        #8. Re-escaling of kde to integrate to 1:
        kde.gr$estimate <- kde.gr$estimate / sum(kde.gr$estimate)
        #9. Creating and plotting image matrix with quantiles and threshold
        results$groups[[i]]$images <- list()
        results$groups[[i]]$images$noThreshold.quantiles <- imageTPD(kde.gr, thresholdPlot = 1)
        results$groups[[i]]$images$TPD.quantiles  <- imageTPD(kde.gr, thresholdPlot = results$parameters$threshold)
      }
    }
    if(results$parameters$objectClass == "TPDs"){
      pcs.gr <- data.frame(results$parameters$pcs, results$parameters$group.vec)
      colnames(pcs.gr) <- c(colnames(results$parameters$pcs), "group.vec")
      pcs.gr.list <- split(pcs.gr, pcs.gr$group.vec)
      for(i in 1:length(pcs.gr.list)){
        results$groups[[i]] <- list()
        names(results$groups)[i] <- names(pcs.gr.list)[i]
        results$groups[[i]]$parameters <- list()
        results$groups[[i]]$parameters$points <- pcs.gr.list[[i]][,c(1,2)]
        results$groups[[i]]$images <- list()
        results$groups[[i]]$images$noThreshold.quantiles <- imageTPDs(x = x, targetTPD = names(results$groups)[i], thresholdPlot = 1)
        results$groups[[i]]$images$TPD.quantiles  <- imageTPDs(x = x, targetTPD = names(results$groups)[i],
                                                               thresholdPlot = results$parameters$threshold)
      }
    }
  }
  #### In the PCA case, some information about the variance explained by the selected components for each trait
  if(results$parameters$objectClass == "princomp"){
    results$PCAInfo <- list()
    results$PCAInfo$pca.object <- x
    results$PCAInfo$fit <- new_var_loading(pca.object = x,
                                           PCs = results$parameters$PCs,
                                           trait_ranges = results$parameters$trait_ranges)

    var.explained.out <- (results$PCAInfo$fit[[2]]**2)[1:nrow(results$PCAInfo$fit[[1]]), c(PCs[1], PCs[2])]
    Overall_explained <- rowSums(var.explained.out)
    results$PCAInfo$varExpl <- data.frame(var.explained.out, Overall_explained)
  }
  #### Estimation of functional richness and divergence
  results$FD <- list()
  if(results$parameters$objectClass != "TPDs"){
      results$FD$global <-list()
      results$FD$global$FRich <- FRichnessGlobal(results)
      results$FD$global$FDiv <- FDivergenceGlobal(results)
    if(!is.null(group.vec)){
      results$FD$groups <-list()
      results$FD$groups$FRich <- FRichnessGroups(results)
      results$FD$groups$FDiv <- FDivergenceGroups(results)
    }
  }
  if(results$parameters$objectClass == "TPDs"){
    results$FD$global <-list()
    results$FD$global$FRich <- NA
    results$FD$global$FDiv <- NA
    results$FD$groups <-list()
    results$FD$groups$FRich <- FRichnessGroups(results)
    results$FD$groups$FDiv <- FDivergenceGroups(results)
  }

  results$parameters$type <- "standard"
  class(results) <- "funspace"

  return(results)
}

