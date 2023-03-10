% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/funspaceDim.R
\name{funspaceDim}
\alias{funspaceDim}
\title{Dimensionality of a trait space}
\usage{
funspaceDim(data)
}
\arguments{
\item{data}{A \code{data.frame} or \code{matrix} containing trait data}
}
\value{
\code{funspaceDim} returns the number of dimensions to be retained. The output is stored and printed out in the R console as well.
}
\description{
Calculating the dimensionality of a functional space based on PCA
}
\details{
\code{funspaceDim} allows the user to identify the number of dimensions that are needed to build a trait space. The identified dimensions are those that minimize redundancy while maximizing the information contained in the trait data. The number of significant PCA axes to be retained is determined by using the \code{paran()} function of the R package \code{paran} (Dinno, 2018). \code{paran()} is based on the method proposed by Horn (1965), which involves contrasting the eigenvalues produced through PCAs run on (30 * (number of variables)) random datasets with the same number of variables and observations as the input dataset. Eigenvalues > 1 are retained in the adjustment.
}
\examples{

# Dimensionality of the GSPFF
funspaceDim(GSPFF)

}
\references{
Horn, J.L. (1965). A rationale and test for the number of factors in factor analysis. Psychometrika 30: 179-185.

Dinno, A. (2018). paran: Horn's test of principal components/factors. R package version 1.5.2.
}
