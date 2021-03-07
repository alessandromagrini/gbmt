# gbmt
__Estimation of group-based multi-trajectory models__

`gbmt` is an R package implementing several functionalities for group-based multi-trajectory models, including maximum likelihood estimation through the Expectation-Maximization algorithm, graphical visualization and prediction.
The reference papers are:

A. Magrini (2021). Assessment of agricultural sustainability in European Union countries: a group-based multi-trajectory approach. _Under review_.

D. S. Nagin, B. L. Jones, V. L. Passos and R. E. Tremblay (2018). Group-based multi-trajectory modeling. _Statistical Methods in Medical Research_, 27(7): 2015-2023.


R (The R Project for Statistical Computing) needs to be installed on your system in order
to use the `gbmt` package. R can be downloaded from https://www.r-project.org/.

To install the `gbmt` package, open the console of R and type:
```
install.packages("devtools")  ## do not run if package 'devtools' is already installed
library(devtools)
install_github("alessandromagrini/gbmt")
```

For any request or feedback, please write to <alessandro.magrini@unifi.it> (Alessandro Magrini)

Below, you find some examples of use of the package.
_________________________________________________________________

Load data on indicators of agricultural sustainability in EU countries, 2004-2018:
```
data(agrisus)
?agrisus          ## data description
summary(agrisus)  ## data summaries
```
Define the indicators:
```
varNames <- c("TFP_2005_CAP","NetCapital_GVA","Manager_ratio",
  "FactorIncome_paid_2010","EntrIncome_unpaid_2010",
  "Income_rur","Unempl_rur","Poverty_rur",
  "RenewProd","Organic_p","GHG_UAA","GNB_UAA")
```
Fit a group-based multi-trajectory model with 2 polynomial degrees and 4 groups:
```
m4_2 <- gbmt(x.names=varNames, unit="Country", time="Year", d=2, ng=4 ,data=agrisus)
```
See the resulting groups:
```
m4_2$assign.list
```
See the group-based trajectories:
```
m4_2$fitted
```
Graphics of the group-based trajectories:
```
plot(m4_2, group=1)  ## group 1
plot(m4_2, group=2)  ## group 2
plot(m4_2, group=3)  ## group 3
plot(m4_2, group=4)  ## group 4
```
Prediction of group-based trajectories:
```
predict(m4_2, n.ahead=3)  ## 3-step ahead
```
The following code allows to fit all models combining 2 and 3 polynomial degrees with a number of groups from 3 to 8:
```
d <- 2:3
ng <- 3:8
modList <- list()
mnam <- c()
for(i in 1:length(d)) {
  for(j in 1:length(ng)) {
    ijnam <- paste("m",ng[j],"_",d[i],sep="")
    mnam <- c(mnam,ijnam)
    ijm <- gbmt(x.names=varNames, unit="Country", time="Year",
                d=d[i], ng=ng[j], data=agrisus)
    modList <- c(modList,list(ijm))
    }
  }
names(modList) <- mnam

# comparison among the models based on BIC
sort(sapply(modList,function(x){x$ic["bic"]}))  ## the lowest BIC is for d=2 and ng=4
```
