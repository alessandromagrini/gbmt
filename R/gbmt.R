# linear interpolation (auxiliary)
linImp <- function(x) {
  res <- x
  auxNA <- which(is.na(x))
  if(length(auxNA)>0) {
    naL <- split(auxNA,cumsum(c(1,diff(auxNA)!=1)))
    for(i in 1:length(naL)) {
      ina <- naL[[i]]
      x1 <- min(ina)-1
      x2 <- max(ina)+1
      if(x1>=1 & x2<=length(x)) {
        y1 <- x[x1]
        y2 <- x[x2]
        b <- (y2-y1)/(x2-x1)
        a <- y1-b*x1
        res[ina] <- a+b*ina
        }
      }
    }
  res
  }

# univariate imputation (auxiliary)
imputFun <- function(x) {
  log.scale <- sum(x<=0,na.rm=T)==0
  ind <- which(is.na(x))
  if(length(ind)>0) {
    if(log.scale) {
      exp(spline(log(x), xout=1:length(x), method="natural")$y)
      } else {
      spline(x, xout=1:length(x), method="natural")$y
      }
    } else {
    x  
    }
  }

# univariate imputation with polynomials (auxiliary)
imputFun2 <- function(x,max.deg=1,sqrt.scale=F) {
  #if(log.scale) ls <- sum(x<=0,na.rm=T)==0 else ls <- F
  if(sqrt.scale) ls <- sum(x<0,na.rm=T)==0 else ls <- F
  ind <- which(is.na(x))
  if(length(ind)>0) {
    max.deg <- max(min(max.deg,length(x)-length(ind)-2),1)
    xt <- 1:length(x)
    currentAIC <- Inf
    for(i in 1:max.deg) {
      if(ls) {
        #iform <- paste("log(x)~poly(xt,",i,",raw=T)",sep="")
        iform <- paste("sqrt(x)~poly(xt,",i,",raw=T)",sep="")
        } else {
        iform <- paste("x~poly(xt,",i,",raw=T)",sep="")
        }
      im <- lm(formula(iform))
      iaic <- AIC(im)
      if(iaic<currentAIC) {
        mOK <- im
        currentAIC <- iaic
        }
      }
    options(warn=-1)
    if(ls) {
      #x[ind] <- exp(predict(mOK,data.frame(xt=xt[ind]))+0.5*summary(mOK)$sigma^2)
      x[ind] <- predict(mOK,data.frame(xt=xt[ind]))^2
      } else {
      x[ind] <- predict(mOK,data.frame(xt=xt[ind]))
      }
    options(warn=0)
    }
  x
  }

# random assignment (auxiliary)
randomCl <- function(n,ng) {
  n0 <- round(n/ng+0.5)
  x <- rep(1:ng,each=n0)
  u <- runif(length(x))
  res <- x[order(u)]
  res[1:n]
  }

# multinormal density (auxiliary)
logdmN <- function(x,mu,S) {
  k <- ncol(S)
  cost <- -k/2*log(2*pi)-0.5*log(det(S))
  ll <- 0
  for(i in 1:nrow(x)){  
    if(sum(is.na(x[i,]))==0) {
      idx <- x[i,]-mu[i,]
      ill <- c(cost-0.5*t(idx)%*%solve(S)%*%(idx))
      ll <- ll+ill
      }
    }
  ll
  }

# fit gbmt
gbmt <- function(x.names,unit,time,ng,d=2,data,tol=1e-6,maxit=500,nstart=NULL,quiet=FALSE) {
  #
  ## check arguments
  #
  ## gestire sbilanciamento degli istanti temporali
  ## check dati positivi
  ## check missing
  #
  dataOK <- do.call(rbind, lapply(split(data, data[,unit]), function(x){
    cbind(x[,setdiff(colnames(x),x.names)],apply(x[,x.names],2,function(y){y/mean(y,na.rm=T)}))
    }))
  if(!is.null(nstart)) {
    if(quiet==F) {
      cat("\r","Restart 1/",nstart,sep="")
      flush.console()
      }
    mOK <- gbmtFit(x.names=x.names,unit=unit,time=time,ng=ng,d=d,data=dataOK,tol=tol,maxit=maxit,quiet=T,init=NULL)
    res <- c(logLik=mOK$logLik,mOK$em)
    for(i in 2:nstart) {
      if(quiet==F) {
        cat("\r","Restart ",i,"/",nstart,sep="")
        flush.console()
        }
      m0 <- gbmtFit(x.names=x.names,unit=unit,time=time,ng=ng,d=d,data=dataOK,tol=tol,maxit=maxit,quiet=T,init=NULL)
      res <- rbind(res,(c(logLik=m0$logLik,m0$em)))
      if(m0$logLik>mOK$logLik) mOK <- m0
      }
    if(quiet==F) cat("\n")
    } else {
    dataI <- dataOK
    cou <- levels(factor(dataOK[,unit]))
    for(i in 1:length(cou)) {
      ind <- which(dataOK[,unit]==cou[i])
      dataI[ind,x.names] <- apply(dataOK[ind,x.names,drop=F],2,imputFun)
      }
    #
    iniList <- list()
    dataList <- split(dataI, dataI[,time])
    for(i in 1:length(dataList)) {
      idat <- dataList[[i]][,x.names,drop=F]
      rownames(idat) <- dataList[[i]][,unit]
      if(sum(apply(idat,2,var)==0)==0) {
        #iini <- cutree(hclust(dist(scale(idat)),method="ward.D"),ng)
        iini <- kmeans(scale(idat),ng,iter.max=maxit,nstart=10)$cluster
        #iini <- Mclust(scale(idat),G=ng)$classification
        } else {
        iini <- NULL
        }
      iniList <- c(iniList, list(iini))
      }
    #
    #dataList <- split(dataI, dataI[,time])
    #cou <- dataList[[1]][,unit]
    #dataList_s <- lapply(dataList, function(x){x[which(x[,unit]%in%cou),x.names]})
    #fulldat <- scale(do.call(cbind,dataList_s))
    #rownames(fulldat) <- cou
    ##iniList <- list(cutree(hclust(dist(fulldat),method="ward.D"),ng))
    #iniList <- list(kmeans(fulldat,ng,iter.max=maxit,nstart=100)$cluster)
    #
    mOK <- list(logLik=-Inf)
    res <- c()
    #if(quiet==F) cat("|")
    for(i in 1:length(iniList)) {
      m0 <- gbmtFit(x.names=x.names,unit=unit,time=time,ng=ng,d=d,data=dataOK,tol=tol,maxit=maxit,quiet=T,init=iniList[[i]])
      #if(quiet==F) cat("*")
      res <- rbind(res,(c(logLik=m0$logLik,m0$em)))
      if(m0$logLik>mOK$logLik) mOK <- m0
      }
    #if(quiet==F) cat("|","\n")
    }
  rownames(res) <- NULL
  mOK$em <- res
  mOK$mean <- apply(data[,x.names,drop=F],2,mean,na.rm=T)
  mOK$data.orig <- data[,c(unit,time,x.names)] 
  mOK$data.scaled <- dataOK[,c(unit,time,x.names)]
  mOK
  }

# fit function (auxiliary)
gbmtFit <- function(x.names,unit,time,ng,d=2,data,tol=1e-6,maxit=150,quiet=F,init=NULL) {
  data[,unit] <- factor(data[,unit])
  tnam <- levels(factor(data[,time]))
  tnum <- as.numeric(factor(data[,time]))
  dataList <- split(data,data[,unit])
  cou <- names(dataList)
  n <- length(dataList)
  np <- length(x.names)
  nt <- length(tnam)
  # inizializzo Z
  if(is.null(init)) {
    Z <- randomCl(n,ng)
    } else {
    Z <- init
    }
  # inizializzo parametri
  Xmat <- c()
  for(i in 0:d) Xmat <- cbind(Xmat,(1:nt)^i)
  rownames(Xmat) <- tnam
  beta <- S <- vector("list",length=ng)
  postd <- matrix(nrow=n,ncol=ng)
  reg <- vector("list",length=ng)
  ll_old <- -Inf
  fine <- nch0 <- count <- 0
  mstr <- ""
  while(fine==0) {
    if(count>0 && quiet==F) {
      mstr <- paste("EM iteration ",count,". Log lik: ",round(ll_old,4),sep="")
      nch <- nchar(mstr)
      if(nch<nch0) {
        estr <- rep(" ",nch0-nch)
        } else {
        estr <- ""
        }
      cat('\r',mstr,estr)
      flush.console()
      nch0 <- nch
      }
    count <- count+1
    # stima pi dato Z
    Zaux <- factor(Z,levels=1:ng)
    pi <- prop.table(table(Zaux))
    # stima beta e S dato Z
    for(i in 1:ng) {
      iOK <- which(Z==i)
      if(length(iOK)>0) {
        iind <- c()
        for(j in 1:length(iOK)) {
          iind <- c(iind,which(data[,unit]==cou[iOK[j]]))
          }
        idat <- data.frame(data[iind,x.names,drop=F],t=tnum[iind])
        iform <- formula(paste("cbind(",paste(x.names,collapse=","),") ~ ",paste("I(t^",1:d,")",collapse="+"),sep=""))
        imod <- lm(iform,data=idat)
        imod$call$formula <- iform
        imod$call$data <- NULL
        #deparse(substitute(data))
        beta[[i]] <- imod$coefficients
        if(np>1) {
          iS <- cov(imod$residuals)
          diag(iS) <- sapply(summary(imod),function(g){g$sigma^2})
          #itry <- try(as.matrix(nearPD(iS)$mat))
          #if(inherits(itry, "try-error")) {
          #  iSok <- diag(rep(1,nrow(iS)))
          #  } else {
          #  iSok <- itry
          #  }
          #S[[i]] <- iSok
          S[[i]] <- as.matrix(nearPD(iS)$mat)
          } else {
          iS <- summary(imod)$sigma^2
          S[[i]] <- iS
          }
        reg[[i]] <- imod
        }
      }
    # stima postd dato pi, beta e sigma
    for(j in 1:ng) {
      ijmu <- Xmat%*%beta[[j]]
      for(i in 1:n) {
        idat <- as.matrix(dataList[[i]][,x.names,drop=F])
        if(np>1) {
          ijD <- logdmN(idat,ijmu,S[[j]])
          } else {
          ijD <- sum(dnorm(idat,ijmu,sqrt(S[[j]]),log=T))
          }
        if(is.na(ijD)) postd[i,j] <- -Inf else postd[i,j] <- ijD+log(pi[j])
        }
      }
    # determino Z
    Z <- apply(postd,1,which.max)
    pOK <- as.numeric(names(table(Z)))
    # check convergence
    ll_new <- sum(log(apply(exp(postd[,pOK]),1,sum)))
    if(ll_new-ll_old<tol | count>=maxit) {
      fine <- 1
      } else {
      ll_old <- ll_new
      }
    }
  if(quiet==F) {
    cat("\n")
    if(count>=maxit) {
      cat("Maximum number of iterations exceeded","\n")
      } else {
      cat("Converged","\n")
      }
    }
  #if(count>=maxit) warning("Maximum number of iterations exceeded")
  conv <- ifelse(count>=maxit,0,1)
  postp <- postd
  for(i in 1:n) postp[i,] <- round(exp(postd[i,])/sum(exp(postd[i,])),4)
  names(Z) <- rownames(postp) <- names(dataList)
  assL <- vector("list",length=ng)
  pred <- list()
  for(i in 1:ng) {
    assL[[i]] <- names(which(Z==i))
    pred[[i]] <- Xmat%*%beta[[i]]
    }
  names(pi) <- names(reg) <- names(assL) <- colnames(postp) <- names(beta) <- names(S) <- names(pred) <- 1:ng
  npar <- ng*(1+(d+1)*np+np*(np+1)/2)-1
  datOK <- data[,c(unit,time,x.names)]
  rownames(datOK) <- NULL
  if(np>1) {
    Rsq <- rsqCalc(reg)
    } else {
    Rsq <- rbind(sapply(reg,function(x){summary(x)$r.squared}))
    rownames(Rsq) <- x.names
    beta <- lapply(beta,function(x){m<-rbind(x);rownames(m)<-x.names;m})
    S <- lapply(S,function(x){m<-matrix(x,nrow=1,ncol=1);rownames(m)<-colnames(m)<-x.names;m})
    pred <- lapply(pred,function(x){colnames(x)<-x.names;x})
    }
  ic <- c(aic=-2*ll_new+2*npar,
          aicc=-2*ll_new+2*npar*(1+(npar+1)/(n*nt-npar-1)),
          bic=-2*ll_new+npar*log(n*nt))
  res <- list(ng=ng, d=d, nt=nt, pi=pi, beta=lapply(beta,t), Sigma=S,
              fitted=pred, reg=reg, posterior=postp, assign=Z, assign.list=assL,
              Rsq=Rsq, logLik=ll_new, ic=ic,
              em=c(n.iter=count,converged=conv))
  class(res) <- "gbmt"
  res
  }

# compute R-squared (auxiliary)
rsqCalc <- function(reg) {
  rsq <- c()
  for(i in 1:length(reg)) {
    isumm <- summary(reg[[i]])
    irsq <- c()
    for(j in 1:length(isumm)) {
      irsq[j] <- isumm[[j]]$r.squared
      }
    rsq <- cbind(rsq,irsq)
    }
  colnames(rsq) <- 1:length(reg)
  rownames(rsq) <- colnames(reg[[1]]$coefficients)
  rsq
  }

# compute confidence bands (auxiliary)
confBands <- function(x,conf=0.95) {
  reg <- x$reg
  nomiVar <- colnames(reg[[1]]$coefficients)
  X <- c()
  for(i in 0:x$d) X <- cbind(X,(1:x$nt)^i)
  mu <- lower <- upper <- vector("list",length=x$ng)
  names(mu) <- names(lower) <- names(upper) <- 1:x$ng
  for(i in 1:x$ng) {
    mu[[i]] <- lower[[i]] <- upper[[i]] <- matrix(nrow=x$nt,ncol=length(nomiVar))
    colnames(mu[[i]]) <- colnames(lower[[i]]) <- colnames(upper[[i]]) <- nomiVar
    }
  for(i in 1:length(reg)) {
    im <- reg[[i]]
    isumm <- summary(im)
    itval <- -qt((1-conf)/2,nobs(im)-x$d-2)
    for(j in 1:length(isumm)) {
      iSb <- isumm[[j]]$cov.unscaled
      ise <- isumm[[j]]$sigma*sqrt(1+diag(X%*%iSb%*%t(X)))
      imu <- X%*%isumm[[j]]$coefficients[,1]
      mu[[i]][,j] <- imu
      lower[[i]][,j] <- imu-itval*ise
      upper[[i]][,j] <- imu+itval*ise
      }
    }
  res <- list(fit=mu,lower=lower,upper=upper)
  attr(res,"conf") <- conf
  res
  }

# print method
print.gbmt <- function(x, ...) {
  gdis <- table(factor(x$assign,levels=1:x$ng))
  for(i in 1:x$ng) {
    cat("Group ",i,": ",gdis[i]," units",sep="","\n")
    print(cbind(x$beta[[i]],'(Rsq)'=x$Rsq[,i]))
    if(i<length(gdis)) cat("\n")
    }
  }

# summary method
summary.gbmt <- function(object, ...) {
  lapply(object$reg,summary,...)
  }

# plot method
plot.gbmt <- function(x, group, equal.scale=FALSE, ylim=NULL, trim=0, titles=NULL, mfrow=NULL, mar=c(5.1,4.1,4.1,2.1), xlab="", ylab="", ...) {
  #
  ## check arguments
  #
  mu <- x$fitted
  np <- ncol(mu[[1]])
  xnam <- colnames(mu[[1]])
  tnam <- rownames(mu[[1]])
  gnam <- x$assign.list[[group]]
  ##conf <- 0.95 
  #confMat <- confBands(c,conf=conf)
  ind <- which(x$data.scaled[,1]%in%gnam)
  if(is.null(mfrow)) mfrow <- n2mfrow(length(xnam))
  par(mar=mar,mfrow=mfrow)
  for(i in 1:length(xnam)) {
    itit <- ifelse(is.null(titles),xnam[i],titles[i])
    if(equal.scale==T) {
      ily <- quantile(x$data.scaled[,xnam[i]],prob=c(trim/2,1-trim/2),na.rm=T)
      } else {
      if(is.null(ylim)) {
        ily <- range(x$data.scaled[ind,xnam[i]],na.rm=T)
        } else {
        if(!is.matrix(ylim)) ily <- ylim else ily <- ylim[i,]  
        }
      }
    plot(0,type="n",ylim=ily,xlim=c(1,x$nt),xaxt="n",xlab=xlab,ylab=ylab,main=itit,...)
    grid(NA,NULL,lty=3,col="grey80")
    abline(v=1:x$nt,lty=3,col="grey80")
    axis(1, at=1:x$nt, labels=tnam, las=2, ...)
    for(j in 1:length(gnam)) {
      ijind <- which(x$data.scaled[,1]==gnam[j])
      lines(1:x$nt,x$data.scaled[ijind,xnam[i]],col="grey70")
      }
    #
    #lines(1:x$nt,confMat$lower[[group]][,xnam[i]])
    #lines(1:x$nt,confMat$upper[[group]][,xnam[i]])
    #
    lines(1:x$nt,mu[[group]][,xnam[i]],lwd=2)
    box()
    }
  par(mar=c(5.1,4.1,4.1,2.1),mfrow=c(1,1))
  }

# residuals method
residuals.gbmt <- function(object, ...) {
  lapply(object$reg, residuals)
  }

# fitted method
fitted.gbmt <- function(object, ...) {
  lapply(object$reg, fitted)
  }

# predict method
predict.gbmt <- function(object, n.ahead=NULL, ...) {
  mu <- object$mean
  if(is.null(n.ahead)) {
    res <- object$fitted
    } else {
    tdat <- data.frame(t=(object$nt+1):(object$nt+n.ahead))
    res <- lapply(object$reg, predict, newdata=tdat, ...)
    }
  resOK <- res
  for(i in 1:length(res)) {
    for(j in 1:ncol(res[[i]])) {
      resOK[[i]][,j] <- res[[i]][,j]*mu[colnames(res[[i]])[j]]
      }
    }
  resOK
  }
