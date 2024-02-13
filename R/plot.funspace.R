#'Functional space plotting
#'
#'Takes a \code{funspace} object produced by \code{funspace()} or \code{funspaceGAM()} and plots the trait probability distribution (TPD) or the map of the response variable (depending of which kind of funspace object is provided) in a functional space.
#'
#'
#'@param x A \code{funspace} object produced by \code{funspace()} or \code{funspaceGAM()}.
#'@param type character indicating whether the plots should represent the global distribution of observations (\code{type = "global"}), or be separated by the groups (\code{type = "groups"}) provided when the \code{funspace} object was created. Defaults to \code{"global"}. In the case of funspace objects based on a TPD function created with the \code{TPD} package, only groups are plotted (there is no "global" distribution).
#'@param which.group when plotting groups, either a character or a number indicating the name (character) or position (number) of a single group to be plotted individually.
#'@param quant.plot Logical, Default is \code{TRUE}. Should contour lines representing quantiles (specified in \code{quant}) be plotted?
#'@param quant A vector specifying the quantiles to be plotted (in case \code{quant.plot} is set to \code{TRUE}. In case a TPD function is plotted, the quantiles represent the quantiles of the trait probability density function (lower quantiles indicate areas with higher probability density). In case a GAM object is plotted, the quantiles of the fitted response variable are plotted. In case a TPD function is plotted, default quantiles are 0.99 (or the selected threshold if it is lower), 0.5 and 0.25. In the GAM alternative, default quantiles are 0.99, 0.5 and 0.25.
#'@param quant.lty type of line to be used to represent quantiles. See \code{lty} argument in \code{graphics::par()}.
#'@param quant.col Color to be used in the quantile lines. Defaults to \code{"grey30"}.
#'@param quant.lwd Line width to be used in the quantile lines. Defaults to 1.
#'@param quant.labels Logical, Default is \code{TRUE}. Should labels be added to quantile lines?
#'@param colors A vector defining the colors of plotted quantiles in the TPD case. Only two colors need to be specified. The first color is automatically assigned to the highest quantile in \code{quantiles} (e.g. 0.99), the second color is assigned to the lowest quantile. These colors are then used to automatically generate a gradient from the greatest to the lowest quantile. Any color is admitted. Default is \code{NULL}, in which case \code{c("yellow", "red")} is used in case a trait probability density function is plotted and to \code{viridis::viridis(5)} in the GAM case.
#'@param ncolors number of colors to include in the color gradients set by \code{colors}. Defaults to 100.
#'@param pnt Logical, defaults to \code{FALSE}. Should data points be added to the functional space?
#'@param pnt.pch Numerical. Graphical parameter to select the type of point to be drawn. Default is set to 19. See \code{pch} argument in \code{graphics::par()}..
#'@param pnt.cex Numerical. Graphical parameter to set the size of the points. Default is 0.5. See \code{cex} argument in \code{graphics::par()}.
#'@param pnt.col Graphical parameter to set the points color. Default is \code{"grey80"}.
#'@param arrows Logical, defaults to \code{FALSE}. In case the functional space is based on a PCA, should the loadings of the original traits be represented by arrows in the functional space?
#'@param arrows.length Numerical. Graphical parameter to set the length of the arrow  (see \code{arrows}). Lower values lead to shorter arrows, which can help to make arrows fit within the represented functional space. Defaults to 1.
#'@param arrows.head Numerical. Graphical parameter to set the length of the arrow head (see \code{arrows}). Defaults to 0.08.
#'@param arrows.col Graphical parameter to set the arrows color (see \code{arrows}). Default is \code{"black"}.
#'@param arrows.label.col Graphical parameter to set the color of the arrows labels  color. Default is \code{"black"}.
#'@param arrows.label.pos Numerical. Graphical parameter to set the position of the arrow labels with respect to the arrow heads. Default is 1.1, which draws arrow labels slightly beyond the arrow heads. A value of 1 means drawing labels on top of arrow heads.
#'@param arrows.label.cex Numerical. Graphical parameter to set the size of arrow labels. Defaults to 1.
#'@param axis.title Logical. Default is \code{TRUE}. Should axes titles be plotted?
#'@param axis.title.x Character. The title to be plotted in the x axis if \code{axis.title} is set to \code{TRUE}. If not specified, a default axis title is plotted.
#'@param axis.title.y Character. The title to be plotted in the y axis if \code{axis.title} is set to \code{TRUE}. If not specified, a default axis title is plotted.
#'@param axis.title.cex Numerical. Graphical parameter to set the size of the axes titles. Default is 1.
#'@param axis.title.line Numerical. Graphical parameter to set the on which margin line to plot axes titles. Default is 2.
#'@param axis.cex Numerical. Graphical parameter to set the size of the axes annotation. Default is 1.
#'@param globalContour Logical, Default is \code{TRUE}. Should a contour line representing the global distribution be plotted when \code{type} is set to \code{"groups"}. Adding a global contour lines provides a common reference for all groups and makes comparisons easier.
#'@param globalContour.quant A vector specifying the quantiles to be plotted (in case \code{globalContour} is set to \code{TRUE}. Defaults to the threshold selected when the provided \code{funspace} object was originally created.
#'@param globalContour.lwd Line width to be used in the global contour lines. Defaults to 3.
#'@param globalContour.lty type of line to be used to represent the global contour lines. See \code{lty} argument in \code{graphics::par()}. Defaults to 1 (a continuous line).
#'@param globalContour.col Graphical parameter to set the color of the global contour lines. Default is \code{"grey50"}.
#'@param xlim the x limits (x1, x2) of the plot.
#'@param ylim the y limits (y1, y2) of the plot.
#'@param ... Other arguments
#'
#'@details Produces default plots. If the input object was generated with \code{funspace()}, the plot shows a bivariate functional trait space displaying trait probability densities (for single or multiple groups). If the input object was  generated with \code{funspaceGAM}, the plot shows a heatmap depicting how a target variable is distributed within the functional trait space (for single or multiple groups).
#'
#'@return No return value. This function is called for its side effect: generating plots.
#'
#' @examples
#'x <- princomp(GSPFF)
#'funtest <- funspace(x = x, PCs = c(1, 2), threshold = 0.95)
#'plot(funtest, type = "global", quant.plot = TRUE, quant.lwd = 2, pnt = TRUE, pnt.cex = 0.1,
#'    pnt.col = rgb(0.1, 0.8, 0.2, alpha = 0.2), arrows = TRUE, arrows.length = 0.7)
#'
#'@import viridis
#'        graphics
#'@export

plot.funspace <- function(x = NULL,
                          type = "global",
                          which.group = NULL,
                          quant.plot = FALSE,
                          quant = NULL,
                          quant.lty = 1,
                          quant.col = "grey30",
                          quant.lwd = 1,
                          quant.labels = TRUE,
                          colors = NULL,
                          ncolors = 100,
                          pnt = FALSE,
                          pnt.pch = 19,
                          pnt.cex = 0.5,
                          pnt.col = "grey80",
                          arrows = FALSE,
                          arrows.length = 1,
                          arrows.head = 0.08,
                          arrows.col = "black",
                          arrows.label.col = "black",
                          arrows.label.pos = 1.1,
                          arrows.label.cex = 1,
                          axis.title = TRUE,
                          axis.title.x = NULL,
                          axis.title.y = NULL,
                          axis.title.cex = 1,
                          axis.title.line = 2,
                          axis.cex = 1,
                          globalContour = TRUE,
                          globalContour.quant = NULL,
                          globalContour.lwd = 3,
                          globalContour.lty = 1,
                          globalContour.col = "grey50",
                          xlim = NULL,
                          ylim = NULL,
                          ...) {
  #1. Set plotting limits
  if(is.null(xlim)){
    xlim <- range(x$parameters$evaluation_grid[, 1])
  }
  if(is.null(ylim)){
    ylim <- range(x$parameters$evaluation_grid[, 2])
  }
  aspUse <- 1
  if(x$parameters$objectClass == "princomp"){
    aspUse <- NA
  }

  ############################################################################
  ############################################################################
  ### CASE 1. Plot TPD ----
  if(x$parameters$type == "standard"){
    if(quant.plot == TRUE){
      if(is.null(quant)){
        quant <- c(min(x$parameters$threshold, 0.99), 0.5, 0.25)
      } else{
        if(max(quant) > x$parameters$threshold){
          warning(paste0("The maximum quantile for contour plotting is higher than the selected threshold."))
        }
      }
    }
    if(is.null(colors)){
      colors <- c("yellow", "red")
    }

    ###     ###     ###     ###     ###
    ### CASE 1.1 Single plot ----
    if(type == "global" & x$parameters$objectClass != "TPDs"){
      im <- x$global$images$TPD.quantiles
      imNoThreshold <- x$global$images$noThreshold.quantiles
      #Color palette
      threshold <- x$parameters$threshold
      gradientColors <- grDevices::colorRampPalette(colors, space = "Lab")
      ColorRamp <- rev(gradientColors(ncolors))
      ColorRamp_ex <- ColorRamp[round(1 + (min(im, na.rm=T)) * ncolors / (threshold)) :
                                  round((max(im, na.rm=T)) * ncolors / (threshold))]
      # make base plot
      graphics::image(x = unique(x$parameters$evaluation_grid[,1]),
            y = unique(x$parameters$evaluation_grid[,2]),
            z = im,
            col = ColorRamp_ex, xaxs="r", yaxs="r", xlab = "", ylab = "",
            axes=F, xlim = xlim, ylim = ylim, main = "", asp = aspUse)
      graphics::box(which="plot")
      graphics::axis(1, tcl=0.3, lwd=0.8, cex = axis.cex)
      graphics::axis(2, las=1, tcl=0.3, lwd=0.8, cex = axis.cex)
      # Adding points
      if(pnt == TRUE){
        graphics::points(x$parameters$pcs[, ], pch = pnt.pch, cex = pnt.cex, col = pnt.col)
      }
      # Adding contour lines
      if(quant.plot == TRUE){
        graphics::contour(x = unique(x$parameters$evaluation_grid[,1]),
                y = unique(x$parameters$evaluation_grid[,2]),
                z = imNoThreshold,
                levels = quant, lwd = quant.lwd, lty = quant.lty, col = quant.col,
                drawlabels = quant.labels, add = TRUE)
      }
      # Plotting arrows
      if(arrows == TRUE & x$parameters$objectClass == "princomp"){
        graphics::arrows(x0 = 0, y0 = 0,
                         x1 = x$PCAInfo$fit[[1]][, 1] * arrows.length,
                         y1 = x$PCAInfo$fit[[1]][, 2] * arrows.length,
                         length = arrows.head, col = arrows.col)
        graphics::text(x = x$PCAInfo$fit[[1]][, 1] * arrows.length * arrows.label.pos,
                       y = x$PCAInfo$fit[[1]][, 2] * arrows.length * arrows.label.pos,
                       labels = rownames(x$PCAInfo$fit[[1]]),
                       col = arrows.label.col, cex = arrows.label.cex)
      }
      # Plotting axes titles
      if(axis.title == TRUE){
        axistTitle <- axisTitles(x, axis.title.x = axis.title.x, axis.title.y = axis.title.y )
        graphics::mtext(axistTitle[1], cex = axis.title.cex, side = 1, line = axis.title.line, adj = 0.5)
        graphics::mtext(axistTitle[2], cex = axis.title.cex, side = 2, line = axis.title.line, adj = 0.5)
      }
    }

    ###     ###     ###     ###     ###
    ### CASE 1.2 Multiple plots ----
    if(type == "groups" | x$parameters$objectClass == "TPDs"){
      #Color palette
      rangeVals <- numeric()
      for(i in 1:length(x$groups)){
        rangeVals <- range(rangeVals, x$groups[[i]]$images$TPD.quantiles, na.rm = TRUE)
      }
      threshold <- x$parameters$threshold
      gradientColors <- grDevices::colorRampPalette(colors, space = "Lab")
      ColorRamp <- rev(gradientColors(ncolors))
      ColorRamp_ex <- ColorRamp[round(1 + (min(rangeVals)) * ncolors / (threshold)) :
                                  round((max(rangeVals)) * ncolors / (threshold))]
      # Plot each group
      groupsPlot <- x$groups
      if(!is.null(which.group)){
        if(length(which.group) != 1){
          stop("'which.group' must have length 1.")
        }
        if(inherits(which.group, "character")){
          if(which.group %in% names(x$groups)){
            x$groups <- x$groups[which.group]
          } else{
            stop("When 'which.group' is a a character, it must be one of the groups contained in x$groups")
          }
        }
        if(inherits(which.group, "numeric")){
          if(which.group %in% 1:length(names(x$groups))){
            x$groups <- x$groups[which.group]
          } else{
            stop("When 'which.group' is a a number, it must be a positive integer not larger than the number of groups contained in x$groups")
          }
        }
      }
      for(i in 1:length(x$groups)){
        groupName <- names(x$groups)[i]
        im <- x$groups[[i]]$images$TPD.quantiles
        imNoThreshold <- x$groups[[i]]$images$noThreshold.quantiles
        # make base plot
        graphics::image(x = unique(x$parameters$evaluation_grid[,1]),
              y = unique(x$parameters$evaluation_grid[,2]),
              z = im,
              col = ColorRamp_ex, xaxs="r", yaxs="r", xlab = "", ylab = "",
            axes=F, xlim = xlim, ylim = ylim, main = "", asp = aspUse)
        graphics::box(which="plot")
        graphics::axis(1,tcl=0.3,lwd=0.8)
        graphics::axis(2, las=1, tcl=0.3,lwd=0.8)
        # Adding global contour lines
        if(globalContour == TRUE & x$parameters$objectClass != "TPDs"){
          if(is.null(globalContour.quant)){
            globalContour.quant <- x$parameters$threshold
          }
          graphics::contour(x = unique(x$parameters$evaluation_grid[,1]),
                  y = unique(x$parameters$evaluation_grid[,2]),
                  z = x$global$images$noThreshold.quantiles,
                  levels = globalContour.quant, lwd = globalContour.lwd, lty = globalContour.lty,
                  col = globalContour.col,
                  drawlabels = F, add=T)
        }
        # Adding points
        if(pnt == TRUE){
          pointsInGroup <- which(x$parameters$group.vec == groupName)
          graphics::points(x$parameters$pcs[pointsInGroup, ], pch = pnt.pch, cex = pnt.cex, col = pnt.col)
        }
        # Adding contour lines
        if(quant.plot == TRUE){
          graphics::contour(x = unique(x$parameters$evaluation_grid[,1]),
                  y = unique(x$parameters$evaluation_grid[,2]),
                  z = imNoThreshold,
                  levels = quant, lwd = quant.lwd, lty = quant.lty, col = quant.col,
                  drawlabels = quant.labels, add=T)
        }
        # Plotting arrows
        if(arrows == TRUE & x$parameters$objectClass == "princomp"){
          graphics::arrows(x0 = 0, y0 = 0,
                           x1 = x$PCAInfo$fit[[1]][, 1] * arrows.length, y1 = x$PCAInfo$fit[[1]][, 2] * arrows.length,
                           length = arrows.head, col = arrows.col)
          graphics::text(x = x$PCAInfo$fit[[1]][, 1] * arrows.length * arrows.label.pos,
                         y = x$PCAInfo$fit[[1]][, 2] * arrows.length * arrows.label.pos,
                         labels = rownames(x$PCAInfo$fit[[1]]),
                         col = arrows.label.col, cex = arrows.label.cex)
        }
        # Adding group names
        graphics::legend("topleft", legend = groupName, bty = "n")
        # Plotting axes titles
        if(axis.title == TRUE){
          axistTitle <- axisTitles(x, axis.title.x = axis.title.x, axis.title.y = axis.title.y )
          graphics::mtext(axistTitle[1], cex = axis.title.cex, side = 1, line = axis.title.line, adj = 0.5)
          graphics::mtext(axistTitle[2], cex = axis.title.cex, side = 2, line = axis.title.line, adj = 0.5)
        }
      }
    }
  }

  ############################################################################
  ############################################################################
  ### CASE 2. Plot GAM ----
  if(x$parameters$type == "gam"){
    if(quant.plot == TRUE){
      if(is.null(quant)){
        quant <- c(0.99, 0.5, 0.25)
      }
    }

    ###     ###     ###     ###     ###
    ### CASE 2.1 Single plot ----
    if(type == "global"){
      if(x$parameters$objectClass != "TPDs"){
        im <- x$global$images$predicted
        #Color palette
        if(is.null(colors)){
          colors <- viridis::viridis(5)
        }
        gradientColors <- grDevices::colorRampPalette(colors, space = "Lab")
        ColorRamp <- gradientColors(ncolors)

        # make base plot
        graphics::image(x = unique(x$parameters$evaluation_grid[,1]),
              y = unique(x$parameters$evaluation_grid[,2]),
              z = im,
              col = ColorRamp, xaxs="r", yaxs="r", xlab = "", ylab = "",
              axes = FALSE, xlim = xlim, ylim = ylim, main = "", asp = aspUse)
        graphics::box(which="plot")
        graphics::axis(1, tcl=0.3, lwd=0.8, cex = axis.cex)
        graphics::axis(2, las=1, tcl=0.3, lwd=0.8, cex = axis.cex)
        # Adding points
        if(pnt == TRUE){
          graphics::points(x$parameters$pcs[, ], pch = pnt.pch, cex = pnt.cex, col = pnt.col)
        }
        # Adding contour lines
        if(quant.plot == TRUE){
          graphics::contour(x = unique(x$parameters$evaluation_grid[,1]),
                  y = unique(x$parameters$evaluation_grid[,2]),
                  z = im,
                  levels = stats::quantile(im, na.rm=T, probs = quant),
                  labels = round(stats::quantile(im, na.rm=T, probs = quant), 2),
                  lwd = quant.lwd, lty = quant.lty, col = quant.col,
                  drawlabels = quant.labels, add = TRUE)
        }
        # Plotting arrows
        if(arrows == TRUE & x$parameters$objectClass == "princomp"){
          graphics::arrows(x0 = 0, y0 = 0,
                           x1 = x$PCAInfo$fit[[1]][, 1] * arrows.length, y1 = x$PCAInfo$fit[[1]][, 2] * arrows.length,
                           length = arrows.head, col = arrows.col)
          graphics::text(x = x$PCAInfo$fit[[1]][, 1] * arrows.length * arrows.label.pos,
                         y = x$PCAInfo$fit[[1]][, 2] * arrows.length * arrows.label.pos,
                         labels = rownames(x$PCAInfo$fit[[1]]),
                         col = arrows.label.col, cex = arrows.label.cex)
        }
        # Plotting axes titles
        if(axis.title == TRUE){
          axistTitle <- axisTitles(x, axis.title.x = axis.title.x, axis.title.y = axis.title.y )
          graphics::mtext(axistTitle[1], cex = axis.title.cex, side = 1, line = axis.title.line, adj = 0.5)
          graphics::mtext(axistTitle[2], cex = axis.title.cex, side = 2, line = axis.title.line, adj = 0.5)
        }
      }
    }
    ###     ###     ###     ###     ###
    ### CASE 2.2 Multiple plots ----
    if(type == "groups"| x$parameters$objectClass == "TPDs"){
      #Color palette
      rangeVals <- numeric()
      for(i in 1:length(x$groups)){
        if(x$groups[[i]]$gam$exist){
          rangeVals <- range(rangeVals, x$groups[[i]]$images$predicted, na.rm = TRUE)
        }
      }
      if(is.null(colors)){
        colors <- viridis::viridis(5)
      }
      gradientColors <- grDevices::colorRampPalette(colors, space = "Lab")
      ColorRamp <- gradientColors(ncolors)
      # Plot each group
      if(!is.null(which.group)){
        if(length(which.group) != 1){
          stop("'which.group' must have length 1.")
        }
        if(inherits(which.group, "character")){
          if(which.group %in% names(x$groups)){
            x$groups <- x$groups[which.group]
          } else{
            stop("When 'which.group' is a a character, it must be one of the groups contained in x$groups")
          }
        }
        if(inherits(which.group, "numeric")){
          if(which.group %in% 1:length(names(x$groups))){
            x$groups <- x$groups[which.group]
          } else{
            stop("When 'which.group' is a a number, it must be a positive integer not larger than the number of groups contained in x$groups")
          }
        }
      }
      for(i in 1:length(x$groups)){
        if(x$groups[[i]]$gam$exist){
          groupName <- names(x$groups)[i]
          im <- x$groups[[i]]$images$predicted
          range_im <- range(x$groups[[i]]$images$predicted, na.rm = TRUE)
          colSeqTotal <- seq(min(rangeVals), max(rangeVals), length = ncolors)
          ColorRamp_ex <- ColorRamp[which.min(abs(colSeqTotal - range_im[1])) :
                                      which.min(abs(colSeqTotal - range_im[2]))]
          # make base plot
          graphics::image(x = unique(x$parameters$evaluation_grid[,1]),
                y = unique(x$parameters$evaluation_grid[,2]),
                z = im,
                col = ColorRamp_ex, xaxs="r", yaxs="r", xlab = "", ylab = "",
                axes=F, xlim = xlim, ylim = ylim, main = "", asp = 1)
          graphics::box(which="plot")
          graphics::axis(1,tcl=0.3,lwd=0.8)
          graphics::axis(2, las=1, tcl=0.3,lwd=0.8)
          # Adding global contour lines
          if(globalContour == TRUE & x$parameters$objectClass != "TPDs"){
            if(is.null(globalContour.quant)){
              globalContour.quant <- x$parameters$threshold
            }
            graphics::contour(x = unique(x$parameters$evaluation_grid[,1]),
                    y = unique(x$parameters$evaluation_grid[,2]),
                    z = x$global$images$noThreshold.quantiles,
                    levels = globalContour.quant, lwd = globalContour.lwd, lty = globalContour.lty,
                    col = globalContour.col,
                    drawlabels = FALSE, add = TRUE)
          }
          # Adding points
          if(pnt == TRUE){
            pointsInGroup <- which(x$parameters$group.vec == groupName)
            graphics::points(x$parameters$pcs[pointsInGroup, ], pch = pnt.pch, cex = pnt.cex, col = pnt.col)
          }
          # Adding contour lines
          if(quant.plot == TRUE){
            graphics::contour(x = unique(x$parameters$evaluation_grid[,1]),
                    y = unique(x$parameters$evaluation_grid[,2]),
                    z = im,
                    levels = stats::quantile(im, na.rm = TRUE, probs = quant),
                    labels = round(stats::quantile(im, na.rm = TRUE, probs = quant), 2),
                    lwd = quant.lwd, lty = quant.lty, col = quant.col,
                    drawlabels = quant.labels, add = TRUE)
          }
          # Plotting arrows
          if(arrows == TRUE & x$parameters$objectClass == "princomp"){
            graphics::arrows(x0 = 0, y0 = 0,
                             x1 = x$PCAInfo$fit[[1]][, 1] * arrows.length, y1 = x$PCAInfo$fit[[1]][, 2] * arrows.length,
                             length = arrows.head, col = arrows.col)
            graphics::text(x = x$PCAInfo$fit[[1]][, 1] * arrows.length * arrows.label.pos,
                           y = x$PCAInfo$fit[[1]][, 2] * arrows.length * arrows.label.pos,
                           labels = rownames(x$PCAInfo$fit[[1]]),
                           col = arrows.label.col, cex = arrows.label.cex)
          }
          # Adding group names
          graphics::legend("topleft", legend = groupName, bty = "n")
          # Plotting axes titles
          if(axis.title == TRUE){
            axistTitle <- axisTitles(x, axis.title.x = axis.title.x, axis.title.y = axis.title.y )
            graphics::mtext(axistTitle[1], cex = axis.title.cex, side = 1, line = axis.title.line, adj = 0.5)
            graphics::mtext(axistTitle[2], cex = axis.title.cex, side = 2, line = axis.title.line, adj = 0.5)
          }
        }else{
          groupName <- names(x$groups)[i]
          graphics::plot(0, type = "n", xlab = "", ylab = "", axes = FALSE)
          graphics::legend("center", legend = "No model", bty = "n")
          graphics::legend("topleft", legend = groupName, bty = "n")
        }
      }
    }
  }
}
