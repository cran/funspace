% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.funspace.R
\name{plot.funspace}
\alias{plot.funspace}
\title{Functional space plotting}
\usage{
\method{plot}{funspace}(
  x = NULL,
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
  ...
)
}
\arguments{
\item{x}{A \code{funspace} object produced by \code{funspace()} or \code{funspaceGAM()}.}

\item{type}{character indicating whether the plots should represent the global distribution of observations (\code{type = "global"}), or be separated by the groups (\code{type = "groups"}) provided when the \code{funspace} object was created. Defaults to \code{"global"}. In the case of funspace objects based on a TPD function created with the \code{TPD} package, only groups are plotted (there is no "global" distribution).}

\item{which.group}{when plotting groups, either a character or a number indicating the name (character) or position (number) of a single group to be plotted individually.}

\item{quant.plot}{Logical, Default is \code{TRUE}. Should contour lines representing quantiles (specified in \code{quant}) be plotted?}

\item{quant}{A vector specifying the quantiles to be plotted (in case \code{quant.plot} is set to \code{TRUE}. In case a TPD function is plotted, the quantiles represent the quantiles of the trait probability density function (lower quantiles indicate areas with higher probability density). In case a GAM object is plotted, the quantiles of the fitted response variable are plotted. In case a TPD function is plotted, default quantiles are 0.99 (or the selected threshold if it is lower), 0.5 and 0.25. In the GAM alternative, default quantiles are 0.99, 0.5 and 0.25.}

\item{quant.lty}{type of line to be used to represent quantiles. See \code{lty} argument in \code{graphics::par()}.}

\item{quant.col}{Color to be used in the quantile lines. Defaults to \code{"grey30"}.}

\item{quant.lwd}{Line width to be used in the quantile lines. Defaults to 1.}

\item{quant.labels}{Logical, Default is \code{TRUE}. Should labels be added to quantile lines?}

\item{colors}{A vector defining the colors of plotted quantiles in the TPD case. Only two colors need to be specified. The first color is automatically assigned to the highest quantile in \code{quantiles} (e.g. 0.99), the second color is assigned to the lowest quantile. These colors are then used to automatically generate a gradient from the greatest to the lowest quantile. Any color is admitted. Default is \code{NULL}, in which case \code{c("yellow", "red")} is used in case a trait probability density function is plotted and to \code{viridis::viridis(5)} in the GAM case.}

\item{ncolors}{number of colors to include in the color gradients set by \code{colors}. Defaults to 100.}

\item{pnt}{Logical, defaults to \code{FALSE}. Should data points be added to the functional space?}

\item{pnt.pch}{Numerical. Graphical parameter to select the type of point to be drawn. Default is set to 19. See \code{pch} argument in \code{graphics::par()}..}

\item{pnt.cex}{Numerical. Graphical parameter to set the size of the points. Default is 0.5. See \code{cex} argument in \code{graphics::par()}.}

\item{pnt.col}{Graphical parameter to set the points color. Default is \code{"grey80"}.}

\item{arrows}{Logical, defaults to \code{FALSE}. In case the functional space is based on a PCA, should the loadings of the original traits be represented by arrows in the functional space?}

\item{arrows.length}{Numerical. Graphical parameter to set the length of the arrow  (see \code{arrows}). Lower values lead to shorter arrows, which can help to make arrows fit within the represented functional space. Defaults to 1.}

\item{arrows.head}{Numerical. Graphical parameter to set the length of the arrow head (see \code{arrows}). Defaults to 0.08.}

\item{arrows.col}{Graphical parameter to set the arrows color (see \code{arrows}). Default is \code{"black"}.}

\item{arrows.label.col}{Graphical parameter to set the color of the arrows labels  color. Default is \code{"black"}.}

\item{arrows.label.pos}{Numerical. Graphical parameter to set the position of the arrow labels with respect to the arrow heads. Default is 1.1, which draws arrow labels slightly beyond the arrow heads. A value of 1 means drawing labels on top of arrow heads.}

\item{arrows.label.cex}{Numerical. Graphical parameter to set the size of arrow labels. Defaults to 1.}

\item{axis.title}{Logical. Default is \code{TRUE}. Should axes titles be plotted?}

\item{axis.title.x}{Character. The title to be plotted in the x axis if \code{axis.title} is set to \code{TRUE}. If not specified, a default axis title is plotted.}

\item{axis.title.y}{Character. The title to be plotted in the y axis if \code{axis.title} is set to \code{TRUE}. If not specified, a default axis title is plotted.}

\item{axis.title.cex}{Numerical. Graphical parameter to set the size of the axes titles. Default is 1.}

\item{axis.title.line}{Numerical. Graphical parameter to set the on which margin line to plot axes titles. Default is 2.}

\item{axis.cex}{Numerical. Graphical parameter to set the size of the axes annotation. Default is 1.}

\item{globalContour}{Logical, Default is \code{TRUE}. Should a contour line representing the global distribution be plotted when \code{type} is set to \code{"groups"}. Adding a global contour lines provides a common reference for all groups and makes comparisons easier.}

\item{globalContour.quant}{A vector specifying the quantiles to be plotted (in case \code{globalContour} is set to \code{TRUE}. Defaults to the threshold selected when the provided \code{funspace} object was originally created.}

\item{globalContour.lwd}{Line width to be used in the global contour lines. Defaults to 3.}

\item{globalContour.lty}{type of line to be used to represent the global contour lines. See \code{lty} argument in \code{graphics::par()}. Defaults to 1 (a continuous line).}

\item{globalContour.col}{Graphical parameter to set the color of the global contour lines. Default is \code{"grey50"}.}

\item{xlim}{the x limits (x1, x2) of the plot.}

\item{ylim}{the y limits (y1, y2) of the plot.}

\item{...}{Other arguments}
}
\value{
No return value. This function is called for its side effect: generating plots.
}
\description{
Takes a \code{funspace} object produced by \code{funspace()} or \code{funspaceGAM()} and plots the trait probability distribution (TPD) or the map of the response variable (depending of which kind of funspace object is provided) in a functional space.
}
\details{
Produces default plots. If the input object was generated with \code{funspace()}, the plot shows a bivariate functional trait space displaying trait probability densities (for single or multiple groups). If the input object was  generated with \code{funspaceGAM}, the plot shows a heatmap depicting how a target variable is distributed within the functional trait space (for single or multiple groups).
}
\examples{
x <- princomp(GSPFF)
funtest <- funspace(x = x, PCs = c(1, 2), threshold = 0.95)
plot(funtest, type = "global", quant.plot = TRUE, quant.lwd = 2, pnt = TRUE, pnt.cex = 0.1,
   pnt.col = rgb(0.1, 0.8, 0.2, alpha = 0.2), arrows = TRUE, arrows.length = 0.7)

}
