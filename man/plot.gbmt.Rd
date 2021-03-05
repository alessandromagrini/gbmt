\name{plot.gbmt}
\alias{plot.gbmt}
\title{Graphics for a group-based multi-trajectory model}
\description{Obtain the graphics of group-based trajectories from a group-based multi-trajectory model.}
\usage{
\method{plot}{gbmt}(x, group, equal.scale=FALSE, ylim=NULL, trim=0,
  titles=NULL, mfrow=NULL, mar=c(5.1,4.1,4.1,2.1), xlab="", ylab="", ...)
}
\arguments{
  \item{x}{An object of class \code{cibn}.}
  \item{group}{Numeric value indicating the group for which the trajectories should be displayed.}
  \item{equal.scale}{...}
  \item{ylim}{...}
  \item{trim}{...}
  \item{titles}{...}
  \item{mfrow}{...}
  \item{mar}{...}
  \item{xlab}{...}
  \item{ylab}{...}
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

# model with 2 degrees and 4 groups
m4_2 <- gbmt(x.names=varNames, unit="Country", time="Year", d=2, ng=4 ,data=agrisus)

# graphics of group-based trajectories
plot(m4_2, group=1)
plot(m4_2, group=2)
plot(m4_2, group=3)
plot(m4_2, group=4)
}
