\name{gbmt}
\alias{gbmt}
\title{Estimation of a group-based multi-trajectory model}
\description{Estimation of a group-based multi-trajectory model through the Expectation-Maximization algorithm.}
\usage{gbmt(x.names, unit, time, ng, d=2, data, tol=1e-6, maxit=500,
  nstart=NULL, quiet=FALSE)}
\arguments{
  \item{x.names}{Character vector including the names of the indicators, that must be strictly positive variables.}
  \item{unit}{Character indicating the name of the variable identifying the units.}
  \item{time}{Character indicating the name of the variable identifying the time.}
  \item{ng}{The number of groups to create.}
  \item{d}{The polynomial degree of the trajectories. Default is 2.}
  \item{data}{Object of class \code{data.frame} containing variables in \code{x.names}, \code{unit} and \code{time}. They cannot contain missing values.}
  \item{tol}{The tolerance of the Expectation-Maximization algorithm. Default is 1e-6.}
  \item{maxit}{The maximum number of iterations of the Expectation-Maximization algorithm. Default is 500.}
  \item{nstart}{The number of random restarts of the Expectation-Maximization algorithm. If \code{NULL} (the default), the Expectation-Maximization algorithm is run with several different starting values, each obtained by applying agglomerative hiearchical clustering with Ward linkage for each time point.}
  \item{quiet}{Logical value indicating whether prompt messages should be displayed during the execution. Effective only if \code{nstart=NULL}. Default is \code{TRUE}.}
  }
\value{An object of class \code{gbmt}, including the following components:
\itemize{
  \item{\code{ng}: }{the number of groups as provided to argument \code{ng}.}
  \item{\code{d}: }{the polynomial degree of the trajectories as provided to argument \code{d}.}
  \item{\code{nt}: }{the number of time points.}
  \item{\code{pi}: }{prior probabilities.}
  \item{\code{beta}: }{estimated regression coefficients for each group.}
  \item{\code{Sigma}: }{estimated covariance matrix of random errors for each group.}
  \item{\code{predicted}: }{estimated group-based trajectories.}
  \item{\code{reg}: }{regression models.}
  \item{\code{posterior}: }{posterior probabilities}
  \item{\code{assign}: }{assignments to the groups.}
  \item{\code{assign.list}: }{assignments to the groups.}
  \item{\code{Rsq}: }{R-squared values.}
  \item{\code{logLik}: }{log-likelihood of the model.}
  \item{\code{ic}: }{information criteria for the model.}
  \item{\code{em}: }{infomration on each run of the Expectation-Maximization algorithm.}
  \item{\code{data.orig}: }{data provided to argument \code{data}.}
  \item{\code{data.scaled}: }{data after rescaling to unit mean.}
  }
}
\details{S3 methods \code{print}, \code{summary} and \code{plot} are available for class \code{gbmt}.
}
%\note{}
\references{
A. Magrini (2021). Assessment of agricultural sustainability in European Union countries: a group-based multi-trajectory approach. \emph{Under revision}.

D. S. Nagin, B. L. Jones, V. L. Passos and R. E. Tremblay (2018). Group-based multi-trajectory modeling. \emph{Statistical Methods in Medical Research}, 27(7): 2015-2023.

D. S. Nagin (2005). Group-based modeling of development. Harvard University Press, Cambridge, US-MA.
}
\seealso{\link{plot.gbmt}.}
\examples{
data(agrisus)

# names of indicators
varNames <- c("TFP_2005_CAP","NetCapital_GVA","Manager_ratio",
  "FactorIncome_paid_2010","EntrIncome_unpaid_2010",
  "Income_rur","Unempl_rur","Poverty_rur",
  "RenewProd","Organic_p","GHG_UAA","GNB_UAA")

# model with 2 degrees and 4 groups
m4_2 <- gbmt(x.names=varNames, unit="Country", time="Year", d=2, ng=4 ,data=agrisus)

# resulting groups
m4_2$assign.list

# regressions by group
m4_2$reg

## NOT RUN: summary of regressions by group
#summary(m4_2)

## NOT RUN: test models combining 2 and 3 degrees with a number of groups from 3 to 8
#d <- 2:3
#ng <- 3:8
#modList <- list()
#mnam <- c()
#for(i in 1:length(d)) {
#  for(j in 1:length(ng)) {
#    ijnam <- paste("m",ng[j],"_",d[i],sep="")
#    mnam <- c(mnam,ijnam)
#    ijm <- gbmt(x.names=varNames, unit="Country", time="Year",
#                d=d[i], ng=ng[j], data=agrisus)
#    modList <- c(modList,list(ijm))
#    }
#  }
#names(modList) <- mnam
#
## comparison among models based on BIC
#sort(sapply(modList,function(x){x$ic["bic"]}))
}
