\name{gbmt}
\alias{gbmt}
\title{Estimation of a group-based multi-trajectory model}
\description{Estimation of a group-based multi-trajectory model through the Expectation-Maximization algorithm.}
\usage{gbmt(x.names, unit, time, ng, d=2, data, tol=1e-6, maxit=500,
  nstart=NULL, quiet=FALSE)}
\arguments{
  \item{x.names}{Character vector including the names of the indicators.}
  \item{unit}{Character indicating the name of the variable identifying the units.}
  \item{time}{Character indicating the name of the variable identifying the time.}
  \item{ng}{The number of groups to create.}
  \item{d}{The polynomial degree of the trajectories. Default is 2.}
  \item{data}{Object of class \code{data.frame} containing the variables indicated in arguments \code{x.names}, \code{unit} and \code{time}. Missing values for these variables are not allowed. Variables indicated in argument \code{x.names} must be numerical with strictly positive values.}
  \item{tol}{The tolerance of the Expectation-Maximization algorithm. Default is 1e-6.}
  \item{maxit}{The maximum number of iterations of the Expectation-Maximization algorithm. Default is 500.}
  \item{nstart}{The number of random restarts of the Expectation-Maximization algorithm. If \code{NULL} (the default), the Expectation-Maximization algorithm is run with several different starting values, each obtained by applying agglomerative hiearchical clustering with Ward linkage for each time point.}
  \item{quiet}{Logical value indicating whether prompt messages should be displayed during the execution (only if \code{nstart=NULL}). Default is \code{TRUE}.}
  }
\value{An object of class \code{gbmt}, including the following components:
\itemize{
  \item{\code{ng}: }{the number of groups as provided to argument \code{ng}.}
  \item{\code{d}: }{the polynomial degree of the trajectories as provided to argument \code{d}.}
  \item{\code{nt}: }{the number of time points.}
  \item{\code{pi}: }{vector including prior probabilities.}
  \item{\code{beta}: }{list of matrices, one for each group, including the estimated regression coefficients.}
  \item{\code{Sigma}: }{list of matrices, one for each group, including the estimated covariance matrix of random errors.}
  \item{\code{fitted}: }{list of matrices, one for each group, including the estimated group-based trajectories.}
  \item{\code{reg}: }{list of objects of class \code{mlm}, one for each group, including the fitted multivariate linear regressions.}
  \item{\code{posterior}: }{matrix including posterior probabilities.}
  \item{\code{assign}: }{vector indicating the group to which each unit is assigned.}
  \item{\code{assign.list}: }{list indicating the units assigned to each group.}
  \item{\code{Rsq}: }{list of matrices, one for each group, including R-squared values.}
  \item{\code{logLik}: }{log-likelihood of the model.}
  \item{\code{ic}: }{information criteria for the model.}
  \item{\code{em}: }{matrix including information on number of iteration and convergence (1=yes, 0=no) for each run of the Expectation-Maximization algorithm.}
  \item{\code{data.orig}: }{data provided to argument \code{data}.}
  \item{\code{data.scaled}: }{data after rescaling to unit mean.}
  }
}
\details{The considered indicators must have strictly positive values, and they are rescaled through division by the respective sample mean.

S3 methods available for class \code{gbmt} include:
\itemize{
  \item{\code{print}: to see regression coefficients for each group}{;}
  \item{\code{summary}: to obtain the summary of the multivariate linear regression for each group}{;}
  \item{\code{plot}: }{to display the group-based trajectories. See \link{plot.gbmt} for details;}
  \item{\code{fitted}: }{to obtain the fitted values for each group;}
  \item{\code{residuals}: }{to obtain the residuals for each group;}
  \item{\code{predict}: }{to make prediction on the group-based trajectories. See \link{predict.gbmt} for details.}
  }
}
%\note{}
\references{
A. Magrini (2021). Assessment of agricultural sustainability in European Union countries: a group-based multi-trajectory approach. \emph{Under revision}.

D. S. Nagin, B. L. Jones, V. L. Passos and R. E. Tremblay (2018). Group-based multi-trajectory modeling. \emph{Statistical Methods in Medical Research}, 27(7): 2015-2023.

D. S. Nagin (2005). Group-based modeling of development. Harvard University Press, Cambridge, US-MA.
}
\seealso{\link{plot.gbmt}; \link{predict.gbmt}.}
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

# group-based trajectories
m4_2$fitted

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
