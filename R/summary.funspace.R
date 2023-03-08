#'Summarizing Functional Spaces
#'
#'\code{summary} method for class \code{funspace}"
#'
#'@param object A \code{funspace} object produced by \code{funspace()},
#'\code{funspaceGAM()}, or \code{funspaceNull()}.
#'@param ... Other arguments
#'
#'@details Produces default summary. If the input object was generated with \code{funspace()}, the summary includes information about the characteristics of the functional space (particularly if it derives from a PCA), along with functional diversity indicators (functional richness and functional divergence) for the whole set of observations and for each group (in case groups are specified). If the input object was generated with \code{funspaceGAM()}, the function returns the summary for the GAM models for the whole set of observations and individual groups. If the input was generated  with \code{funspaceNull()}, the function returns tests exploring the difference between the observed functional richness and the null model functional richness.
#'
#'@return No return value. This function is called for its side effect: summarizing objects of class \code{"funspace"}.
#'
#' @examples
#'x <- princomp(GSPFF)
#'funtest <- funspace(x = x, PCs = c(1, 2), threshold = 0.95)
#'summary(funtest)
#'
#'
#'@export

summary.funspace <- function(object, ...){
  ############################################################################
  ### CASE 1. object contains TPD ----
  ### Summary returns info on the functional space characteristics, plus functional richness and divergence
  if(object$parameters$type == "standard"){
    if(object$parameters$princomp == TRUE){
      cat(paste0("\nFunctional space based on a PCA with ", ncol(object$PCAInfo$fit$loadings), " dimensions"))
      cat(paste0("\nDimensions ", object$parameters$PCs[1], " and ", object$parameters$PCs[2],
                 " are considered in analyses\n"))
      cat(paste0("\nLoadings:\n"))
      print(round(object$PCAInfo$fit$loadings, 3))

      cat(paste0("\nPercentage of variance explained for each trait:\n"))
      print(round(100 * object$PCAInfo$varExpl, 2))
    }
    if(object$parameters$princomp == FALSE){
      cat(paste0("\nFunctional space not based on PCA\n"))
    }

    cat(paste0("\n------------------------------------------------------------\n"))

    cat(paste0("\nFunctional diversity indicators:\n"))

    cat(paste0("\n ---> For the global set of species:\n"))
    cat(paste0("\nFunctional richness (",
               100 * object$parameters$threshold, "% probability threshold) = ",
               round(object$FD$global$FRich, 2), "\n"))
    cat(paste0("\nFunctional divergence  = ",
               round(object$FD$global$FDiv, 2), "\n"))

    if(!is.null(object$groups)){
      cat(paste0("\n ---> For each group:\n"))
      cat(paste0("\nFunctional richness (",
                 100 * object$parameters$threshold, "% probability threshold; ",
                 names(object$FD$groups$FRich), ") = ", round(object$FD$groups$FRich, 2)))
      cat(paste0("\n"))
      cat(paste0("\nFunctional divergence (", names(object$FD$groups$FDiv), ") = ",
                 round(object$FD$groups$FDiv, 2)),"\n")

    }
  }
  ############################################################################
  ### CASE 2. object contains GAM ----
  ### Summary returns info on the GAM models
  if(object$parameters$type == "gam"){
    if(object$parameters$princomp == TRUE){
      cat(paste0("\nFunctional space based on a PCA with ", ncol(object$PCAInfo$fit$loadings), " dimensions"))
      cat(paste0("\nDimensions ", object$parameters$PCs[1], " and ", object$parameters$PCs[2],
               " are considered in analyses\n"))
    }
    if(object$parameters$princomp == FALSE){
      cat(paste0("\nFunctional space not based on PCA\n"))
    }
    cat(paste0("\n------------------------------------------------------------\n"))
    cat(paste0("---> GAM model for the global set of species:\n"))
    print(summary(object$global$gam))

    if(!is.null(object$groups)){
      cat(paste0("\n------------------------------------------------------------\n"))
      cat(paste0("---> GAM models for individual groups of species:\n"))
      for(i in 1:length(object$groups)){
        cat(paste0("\n\n -------------------> ", names(object$groups)[i] ,":\n"))
        if(object$groups[[i]]$gam$exist == FALSE){
          cat(paste0("NOT ENOUGH OBSERVATIONS TO FIT MODEL\n"))
        } else{
        print(summary(object$groups[[i]]$gam))
        }
      }
    }
  }
  ############################################################################
  ### CASE 3. object contains Null model testing ----
  ### Summary returns info on the null vs. observed space
  if(object$parameters$type == "null"){


    cat(paste0("\n------------------------------------------------\n"))


    cat(paste0("\nDifference between observed and null trait space\n"))

    cat(paste0("\n--->  Null distribution:"," ",object$Distr.Test, "\n"))

    cat(paste0("\nObserved trait space area: ", round(object$ObsFric, 3), "\n"))
    cat(paste0("\nNull trait space area: ", round(object$NullFric, 3), "\n"))
    cat(paste0("\np-value: ", round(object$p_value, 4), "\n"))
    cat(paste0("\nStand. Effect Size: ", round(object$StdEffSize, 2), "\n"))
  }
}
