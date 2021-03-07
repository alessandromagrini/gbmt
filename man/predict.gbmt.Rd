\name{predict.gbmt}
\alias{predict.gbmt}
\title{Prediction of group-based trajectories}
\description{Make prediction on group-based trajectories for a group-based multi-trajectory model.}
\usage{
\method{predict}{gbmt}(object, n.ahead=NULL, ...)
}
\arguments{
  \item{object}{An object of class \code{gbmt}.}
  \item{n.ahead}{Numeric value indicating the number of steps ahead to perform prediction. If \code{NULL} (the default), group-based trajectories in the observed period are returned.}
  \item{...}{Further arguments for the generic \code{predict} method.}
}
\value{A list with one component for each group. Each component includes an object of class \code{data.frame} with one entry for each step ahead at which prediction is made.}
%\details{}
\note{Predictions are made on the original scales of indicators.}
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

# group-based trajectories in the observed period
predict(m4_2)

# 3-step ahead prediction og group-based trajectories
predict(m4_2, n.ahead=3)
}
