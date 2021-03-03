\name{gbmt}
\alias{gbmt}
\title{Estimation of a group-based multi-trajectory model}
\description{Estimation of a group-based multi-trajectory model through the Expectation-Maximization algorithm.}
\usage{gbmt(x.names, unit, time, ng, d=2, data, tol=1e-6, maxit=500, nstart=NULL, quiet=FALSE)}
\arguments{
  \item{x.names}{...}  
  \item{unit}{...}  
  \item{time}{...}  
  \item{ng}{...}
  \item{d}{...}
  \item{data}{...}
  \item{tol}{...}
  \item{maxit}{...}
  \item{nstart}{...}
  \item{quiet}{...}  
  }
\value{An object of class \code{gbmt}.}
\details{

S3 methods \code{print}, \code{summary} and \code{plot} are available for class \code{cibn}.
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

# models
d <- 2:3
ng <- 3:8
modList <- list()
mnam <- c()
for(i in 1:length(d)) {
  for(j in 1:length(ng)) {
    ijnam <- paste("m",ng[j],"_",d[i],sep="")
    mnam <- c(mnam,ijnam)
    ijm <- gbmt(x.names=varNames,unit="Country",time="Year",d=d[i],ng=ng[j],data=agrisus)
    modList <- c(modList,list(ijm))
    }
  }
names(modList) <- mnam

# comparison among models based on BIC
sort(sapply(modList,function(x){x$ic["bic"]}))

# groups identified by the final model
modList$m4_2$assign.list

# graphics for the final model
plot(modList$m4_2, group=1)
plot(modList$m4_2, group=2)
plot(modList$m4_2, group=3)
plot(modList$m4_2, group=4)
}
