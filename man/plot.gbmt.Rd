\name{plot.gbmt}
\alias{plot.gbmt}
\title{Graphics for a group-based multi-trajectory model}
\description{Obtain the graphics of group-based trajectories from a group-based multi-trajectory model.}
\usage{
\method{plot}{gbmt}(x, group, equal.scale=FALSE, ylim=NULL, trim=0,
  titles=NULL, mfrow=NULL, mar=c(5.1,4.1,4.1,2.1), xlab="", ylab="", ...)
}
\arguments{
  \item{x}{An object of class \code{gbmt}.}
  \item{group}{Numeric value indicating the group for which the trajectories should be displayed.}
  \item{equal.scale}{Logical indicating whether each indicator should have the same scale across the groups. Default is \code{TRUE}.}
  \item{ylim}{vector of length 2 indicating the limits of the y-axis, which will be applied to all the indicators. If \code{NULL} (the default), it will be determined independently for each indicator based on data.}
  \item{trim}{Numeric value between 0 and 1 indicating the order of extreme quantiles to be trimmed in the graphical representation. Default is 0.}
  \item{titles}{vector of titles for the indicators. If \code{NULL}, the name of the indicators is used as title.}
  \item{mfrow}{vector of length 2 indicating the partition into panels of the graphical area. If \code{NULL} (the default), it will be automatically determined based on the number of indicators.}
  \item{mar}{vector of length 4 indicating the margin size in the order: bottom, left, top, right, which will be applied to all the indicators. Default is \code{c(5.1,4.1,4.1,2.1)}.}
  \item{xlab}{label for the x-axis, which will be applied to all the indicators. Default is empty string.}
  \item{ylab}{label for the y-axis, which will be applied to all the indicators. Default is empty strings.}
  \item{...}{Further graphical parameters.}
}
%\value{}
%\details{}
%\note{}
%\references{}
\seealso{\link{gbmt}.}
\examples{
data(agrisus)

# names of indicators
varNames <- c("TFP_2005_CAP","NetCapital_GVA","Manager_ratio",
  "FactorIncome_paid_2010","EntrIncome_unpaid_2010",
  "Income_rur","Unempl_rur","Poverty_rur",
  "RenewProd","Organic_p","GHG_UAA","GNB_UAA")

# model with 2 polynomial degrees and 4 groups
m4_2 <- gbmt(x.names=varNames, unit="Country", time="Year", d=2, ng=4 ,data=agrisus)

# graphics of group-based trajectories
plot(m4_2, group=1)
plot(m4_2, group=2)
plot(m4_2, group=3)
plot(m4_2, group=4)
}
