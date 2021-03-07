\name{gbmt-package}
\alias{gbmt-package}
\docType{package}
\title{Group-Based Multi-Trajectory Modeling}
\description{Estimation and analysis of group-based multi-trajectory models.}
\details{
  \tabular{ll}{
  Package: \tab gbmt\cr
  Type: \tab Package\cr
  Version: \tab 0.0\cr
  Date: \tab 2021-03-04\cr
  License: \tab GPL-2
  }
Group-based trajectory modeling is a statistical method to create groups of units based on the similarity in the multivariate time series of a set of variables.
It is a special case of latent class growth curves where the units in the same group have the same trajectory (Nagin, 2005), but it assumes a multivariate Normal linear regression within each group, instead of a univariate one, to account for multiple indicators (Nagin et al., 2018; Magrini, 2021).
The main functions of the package are:
\itemize{
\item{\link{gbmt}: }{to estimate a group-based multi-trajectory model through the Expectation-Maximization algorithm;}
\item{\link{plot.gbmt}:} {to obtain the graphics of group-based trajectories;}
\item{\link{predict.gbmt}:} {to make prediction on group-based trajectories.}
}
}
\author{Alessandro Magrini <alessandro.magrini@unifi.it>}
\references{
A. Magrini (2021). Assessment of agricultural sustainability in European Union countries: a group-based multi-trajectory approach. \emph{Under review}.

D. S. Nagin, B. L. Jones, V. L. Passos and R. E. Tremblay (2018). Group-based multi-trajectory modeling. \emph{Statistical Methods in Medical Research}, 27(7): 2015-2023.

D. S. Nagin (2005). Group-based modeling of development. Harvard University Press, Cambridge, US-MA.
}
