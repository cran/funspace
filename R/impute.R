#' Imputing Trait Information
#'
#' Imputing incomplete trait information, with the possibility of using phylogenetic information
#'
#'@details
#'
#' \code{impute} imputes trait values in trait matrices with incomplete trait information. It uses the Random Forest approach implemented in the \code{missForest} package. Phylogenetic information can be incorporated in the imputation in the form of a phylogenetic tree, from which a number of phylogenetic eigenvectors are added to the trait matrix.
#'
#'@param traits A matrix or data.frame containing trait information with missing values. The rows correspond to observations (generally species) and the columns to the variables (generally traits). Traits can be continuous and/or categorical. Row names of the \code{traits} object must contain the names of the species. We recommend writing species name in the format "Genus_species" or "Genus species".
#'@param phylo (optional) A phylogenetic tree (an object of class "phylo") containing the evolutionary relationships between species. \code{phylo} is used to estimate phylogenetic eigenvectors that are added to the \code{traits} matrix. Not all species in \code{traits} need to be necessarily included in \code{phylo}, despite this is highly recommended. Note that in order to assign phylogenetic information to species reliably, the names in \code{phylo$tip.label} must be exactly the same as \code{row.names(traits)}, although not necessarily in the same order. Note that computing cophenetic distances for very large trees (ca. 30,000 species) can result in memory allocation problems.
#'@param addingSpecies Logical, defaults to FALSE. Should species present in the trait matrix but not in the phylogeny be added to the phylogeny? If TRUE, the \code{phytools::add.species.to.genus} function is used to add species to the root of the genus (in case there are any other congeneric species in the tree). Note that \code{phytools::add.species.to.genus} has other arguments that provide more flexibility, but those are not considered here for simplicity; users who want to make use of those options can instead modify their phylogenetic tree beforehand.
#'@param nEigen The number of phylogenetic eigenvectors to be considered. Defaults to 10.
#'@param messages Logical, defaults to TRUE. Should the function return messages?

#'@return The function returns a list containing both the original trait data (incomplete) and the imputed trait data.
#'
#'
#' @examples
#'
#' # GSPFF_missing dataset includes >10,000 species.
#' # Preparing and imputing this data takes very long time.
#' # Let's select a small random subset:
#'selectSPS <- 200
#'set.seed(2)
#'subset_traits <- GSPFF_missing[sample(1:nrow(GSPFF_missing), selectSPS), ]
#'deleteTips <- setdiff(phylo$tip.label, rownames(subset_traits))
#'subset_phylo <- ape::drop.tip(phylo, tip = deleteTips)
#'GSPFF_subset <- impute(traits = subset_traits, phylo = subset_phylo, addingSpecies = TRUE)
#'pca <- princomp(GSPFF_subset$imputed)
#'funtest <- funspace(pca)
#'plot(funtest, pnt = TRUE, pnt.cex = 0.2, arrows = TRUE)
#'summary(funtest)
#'
#'
#' @import  ape
#'          missForest
#'          phytools
#' @export

impute <- function(traits, phylo = NULL, addingSpecies = FALSE, nEigen = 10,
                   messages = TRUE){
  # 1. INITIAL CHECKS:----
	# 1.1 traits must be a matrix or data.frame.
	if (!is.matrix(traits) & !is.data.frame(traits)){
		stop("'traits' must be of class 'matrix' or 'data.frame'")
	}
  # 1.2. traits must contain missing information
  if (sum(is.na(traits)) == 0){
    stop("'traits' must contain missing information")
  }
	# Creation of lists to store results:
	results <- list()
	# 2. PREPARATION OF PHYLOGENETIC EIGENVECTORS:----
	# 2.1 phylo must be an object of 'phylo' class.
	if (!is.null(phylo)) {
	  if (!inherits(phylo, "phylo")){
	    stop("'phylo' must be of class 'phylo'")
	  }

  	# 2.2 Keep in phylo only species present in the trait matrix
  	tipsToDrop <- which(! phylo$tip.label %in% rownames(traits))
  	droptreetips <- phylo$tip.label[tipsToDrop]
  	phylogenyTraits <- ape::drop.tip(phylo, droptreetips)

  	# 2.3 Add species to phylogeny
  	if(addingSpecies){
  	  tipsToAdd <- which(!rownames(traits) %in% phylogenyTraits$tip.label)
  	  namesToAdd <- rownames(traits)[tipsToAdd]
  	  phylogenyTraitsAdd <- phylogenyTraits
  	  if(length(namesToAdd) > 0){
  	    for(i in 1:length(namesToAdd)){
  	      if(messages){
  	        cat(paste("\r adding species", i ,"/", length(namesToAdd), "to phylogeny\r"))
  	      }
  	      phylogenyTraitsAdd <- phytools::add.species.to.genus(tree = phylogenyTraitsAdd,
  	                                                           species = namesToAdd[i],
  	                                                           genus=NULL, where="root")
  	    }
  	    if(messages){message(cat(paste0("\n", length(phylogenyTraitsAdd$tip.label) -
  	                                      length(phylogenyTraits$tip.label),
  	                                    " species added to phylogeny")))
  	    }
  	  }
  	  phylDiss <- sqrt(stats::cophenetic(phylogenyTraitsAdd))
  	}
  	# 2.4 Estimate phylogenetic PCoA
  	phylDiss <- sqrt(stats::cophenetic(phylogenyTraits))
  	pcoaPhyl <- stats::cmdscale(phylDiss, k = nEigen)
  	# 2.5 Add phylogenetic eigenvectors to traits
  	colnames(pcoaPhyl) <- paste0("Eigen.", 1:ncol(pcoaPhyl))
  	missingInPhyl <- setdiff(rownames(traits), rownames(pcoaPhyl))
  	spNoPhyl <- matrix(NA, nrow = length(missingInPhyl), ncol = nEigen,
  	                   dimnames = list(missingInPhyl, colnames(pcoaPhyl)))
  	pcoaPhyl <- rbind(pcoaPhyl, spNoPhyl)
  	pcoaPhyl <- pcoaPhyl[rownames(traits), ]
  	traitsToImpute <- cbind(traits, pcoaPhyl)
	} else{
	  traitsToImpute <- traits
	}
	# 3. TRAIT IMPUTATION:----
	if(messages) imputed <- missForest::missForest(xmis= traitsToImpute)$ximp[, colnames(traits)]
	if(! messages) imputed <- quiet(missForest::missForest(xmis= traitsToImpute)$ximp[, colnames(traits)])

	# 4. RETURNING RESULTS
	results$imputed <- imputed
	results$originalTraits <- traitsToImpute
  return(results)
}

