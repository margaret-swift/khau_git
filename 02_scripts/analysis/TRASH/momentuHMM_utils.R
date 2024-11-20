momentuHMMdists<-sort(c('gamma','weibull','exp','lnorm','beta','pois','wrpcauchy','vm','norm','bern','vmConsensus','mvnorm2','mvnorm3','rw_mvnorm2','rw_mvnorm3','rw_norm','cat','negbinom','logis','t'))
moveHMMdists<-sort(c('gamma','weibull','exp','lnorm','wrpcauchy','vm'))
angledists<-sort(c('wrpcauchy','vm','vmConsensus'))
stepdists<-sort(c('gamma','weibull','exp','lnorm'))
singleParmdists<-sort(c('exp','pois','bern'))
nonnegativedists<-sort(c('gamma','weibull','exp','lnorm','pois','negbinom'))
zeroInflationdists<-sort(c('gamma','weibull','exp','lnorm','beta'))
oneInflationdists<-sort(c('beta'))
integerdists<-sort(c('bern','pois','cat','negbinom'))
mvndists <- c('mvnorm2','mvnorm3','rw_mvnorm2','rw_mvnorm3')
rwdists <- c('rw_norm','rw_mvnorm2','rw_mvnorm3')
splineList<-c("bs","ns","bSpline","mSpline","cSpline","iSpline")
meansList<-c("matrix","numeric","integer","logical","Date","POSIXlt","POSIXct","difftime")
meansListNoTime<-c("numeric","integer","logical")
plotArgs <- c("cex","cex.main","cex.lab","cex.axis","cex.legend","lwd","asp","legend.pos")
fitMethods<-c("nlm","Nelder-Mead","SANN")
badNames <- c("beta", "delta", "pi", "g0", "theta")

#' @importFrom stats dbinom
dbern <- function (x, prob, log = FALSE) 
{
  return(stats::dbinom(x, 1, prob, log))
}

#' @importFrom stats pbinom
pbern <- function (q, prob, lower.tail = TRUE, log.p = FALSE) 
{
  return(stats::pbinom(q, 1, prob, lower.tail, log.p))
}

#' @importFrom stats rbinom
rbern <- function (n, prob) 
{
  return(stats::rbinom(n, 1, prob))
}

#' @importFrom stats dnbinom
dnegbinom <- function (x, mu, size, log = FALSE) 
{
  return(stats::dnbinom(x, size = size, mu = mu, log = log))
}

#' @importFrom stats pnbinom
pnegbinom <- function (q, mu, size, lower.tail = TRUE, log.p = FALSE) 
{
  return(stats::pnbinom(q, size = size, mu = mu, lower.tail = lower.tail, log.p = log.p))
}

#' @importFrom stats rnbinom
rnegbinom<- function (n, mu, size) 
{
  return(stats::rnbinom(n, size = size, mu = mu))
}

dmvnorm2 <- dmvnorm3 <- drw_mvnorm2 <- drw_mvnorm3 <- function(x,mean,sigma){
  dmvnorm_rcpp(x,mean,sigma)
}

rmvnorm2 <- rmvnorm3 <- rrw_mvnorm2 <- rrw_mvnorm3 <- mvtnorm::rmvnorm

drw_norm <- stats::dnorm
rrw_norm <- stats::rnorm

RWdata <- function(dist,data,knownStates){
  distnames <- names(dist)
  if(any(unlist(dist) %in% rwdists)){
    newdata <- NULL
    colInd <- NULL
    if(length(knownStates)){
      if("knownStates" %in% colnames(data)) stop("data cannot include a column named 'knownStates'")
      data$knownStates <- knownStates
    }
    ID <- unique(data$ID)
    for(j in ID){
      jInd <- which(data$ID==j)
      for(i in distnames){
        if(dist[[i]] %in% rwdists){
          tmpdata <- ldata <- data[jInd,,drop=FALSE]
          lInd <- 1:nrow(tmpdata)
          if(dist[[i]] %in% mvndists){
            if(inherits(data,"hierarchical")){
              iLevel <- attr(data,"coordLevel")
              lInd <- which(tmpdata$level==iLevel)
              ldata <- tmpdata[lInd,]
              colInd <- NULL
            }
            tmpdata[[paste0(i,".x_tm1")]] <- tmpdata[[paste0(i,".x")]]
            tmpdata[[paste0(i,".y_tm1")]] <- tmpdata[[paste0(i,".y")]]
            ldata[[paste0(i,".x_tm1")]] <- ldata[[paste0(i,".x")]]
            ldata[[paste0(i,".y_tm1")]] <- ldata[[paste0(i,".y")]]
            if(dist[[i]]=="rw_mvnorm2"){
              colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% c(paste0(i,".x"),paste0(i,".y"))))]))
            } else if(dist[[i]]=="rw_mvnorm3"){
              tmpdata[[paste0(i,".z_tm1")]] <- tmpdata[[paste0(i,".z")]]
              ldata[[paste0(i,".z_tm1")]] <- ldata[[paste0(i,".z")]]
              colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% c(paste0(i,".x"),paste0(i,".y"),paste0(i,".z"))))]))
            }
          } else {
            if(inherits(data,"hierarchical")){
              iLevel <- attr(data,"coordLevel")
              lInd <- which(tmpdata$level==iLevel)
              ldata <- tmpdata[lInd,]
              colInd <- NULL
            }
            tmpdata[[paste0(i,"_tm1")]] <- tmpdata[[i]]
            ldata[[paste0(i,"_tm1")]][lInd] <- ldata[[i]]
            colInd <- unique(c(colInd,colnames(tmpdata)[which(!(colnames(tmpdata) %in% distnames))]))
          }
          if(inherits(data,"hierarchical")){
            ldata[,colInd] <- rbind(rep(NA,length(colInd)),ldata[-nrow(ldata),colInd])
            ldata <- ldata[-1,,drop=FALSE]
            tmpdata[lInd,colInd] <- rbind(rep(NA,length(colInd)),tmpdata[lInd[-length(lInd)],colInd])
            tmpdata <- tmpdata[-lInd[1],,drop=FALSE]
            tmpdata[lInd[-1]-1,colnames(tmpdata)] <- ldata[,colnames(tmpdata)]
            tmpdata[which(tmpdata$level!=iLevel),paste0(colnames(tmpdata)[!colnames(tmpdata) %in% colInd],"_tm1")] <- 0 # can't have NAs in covariates
          }
        }
      }
      if(!inherits(data,"hierarchical")){
        lastRow <- tmpdata[nrow(tmpdata),]
        tmpdata[,colInd] <- rbind(rep(NA,length(colInd)),tmpdata[-nrow(tmpdata),colInd])
        tmpdata <- tmpdata[-1,,drop=FALSE]
        tmpdata <- rbind(tmpdata,lastRow)
        tmpdata[nrow(tmpdata),colnames(tmpdata)[which(!(colnames(tmpdata) %in% colInd))]] <- NA
      }
      newdata <- rbind(newdata,tmpdata)
    }
    class(newdata) <- class(data)
  } else newdata <- data
  newdata
}

# @importFrom dplyr lag
crw <- function(x_tm1,lag=1){
  for(pkg in c("dplyr")){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop("Package \"",pkg,"\" needed for crw function to work. Please install it.",
           call. = FALSE)
    }
  }
  dplyr::lag(x_tm1,n=lag-1,default=x_tm1[1])-dplyr::lag(x_tm1,n=lag,default=x_tm1[1])
}

radian <- function(degree) 
{
  radian <- degree * (pi/180)
  radian
}

# startup message
#' @importFrom utils packageDescription available.packages
print.momentuHMM.version <- function()
{ pkgDescr <- utils::packageDescription("momentuHMM")
  hello <- paste("momentuHMM ",pkgDescr$Version," (",pkgDescr$Date,")",sep="")
  curVersion <- tryCatch(suppressWarnings(utils::available.packages(repos = "http://cran.us.r-project.org")["momentuHMM","Version"]),error=function(e) e)
  packageStartupMessage(hello)
  if(!inherits(curVersion,"error")){
    if(pkgDescr$Version<curVersion) warning("  A newer version (",curVersion,") is available from CRAN")
  }
}

.onAttach <- function(...) { 
  print.momentuHMM.version()
}

# suppress RNG warning when using %dorng%
muffleRNGwarning <- function(w) {
  if(any(grepl("Foreach loop \\(doParallelMC\\) had changed the current RNG type: RNG was restored to same type, next state",w))
     | any(grepl("already exporting variable\\(s\\):",w)))
    invokeRestart("muffleWarning")
}

# .combine function for multiple rbinds in foreach
comb <- function(x, ...) {  
    mapply(rbind,x,...,SIMPLIFY=FALSE)
}

# #' @importFrom doFuture registerDoFuture
#' @importFrom doRNG %dorng%
#' @importFrom foreach %dopar% foreach
# #' @importFrom future multisession plan
# #' @importFrom iterators icount
progBar <- function(kk, N, per = 1) {
  if (kk %in% seq(1, N, per)) {
    x <- round(kk * 100 / N)
    message("[ ", 
            paste(rep("=", x), collapse = ""),
            paste(rep("-", 100 - x), collapse = ""), 
            " ] ", x, "%", "\r",
            appendLF = FALSE)
    if (kk == N) message("\n")
  }
}

installDataTree <- function(){
  if (!requireNamespace("data.tree", quietly = TRUE)) {
    stop("Package \"data.tree\" needed for this function to work. Please install it.",
         call. = FALSE)
  }
}

#' @importFrom MASS ginv
# this function maintains backwards compatibility with momentuHMM versions <1.4.0 (workBounds), <1.4.3 (betaCons), and <1.5.0 (mixtures)
delta_bc <- function(m){
  
  if(is.momentuHMM(m) | is.miSum(m)){
    if(!is.null(m$conditions$fit)){
      if(!m$conditions$fit) warning("The given model hasn't been fitted.")
    } else m$conditions$fit <- TRUE
    if(is.null(m$conditions$workBounds)){
      distnames <- names(m$conditions$dist)
      
      parCount<- lapply(m$conditions$fullDM,ncol)
      for(i in distnames[!unlist(lapply(m$conditions$circularAngleMean,isFALSE))]){
        parCount[[i]] <- length(unique(gsub("cos","",gsub("sin","",colnames(m$conditions$fullDM[[i]])))))
      }
      parindex <- c(0,cumsum(unlist(parCount))[-length(m$conditions$fullDM)])
      names(parindex) <- distnames
      
      workBounds <- vector('list',length(distnames))
      names(workBounds) <- distnames
      if(is.miSum(m)){
        beta <- m$Par$beta$beta$est
        delta <- m$Par$beta$delta$est
      } else {
        beta <- m$CIbeta$beta$est
        delta <- m$CIbeta$delta$est
      }
      beta <- list(beta=beta,g0=m$mle$g0,theta=m$mle$theta)
      m$conditions$workBounds <- getWorkBounds(workBounds,distnames,m$mod$estimate,parindex,parCount,m$conditions$DM,beta,delta)
    }
    if(length(m$stateNames)>1 && is.null(m$conditions$betaCons)){
      if(is.miSum(m) & !is.null(m$Par$beta$beta)) m$conditions$betaCons <- matrix(1:length(m$Par$beta$beta$est),nrow(m$Par$beta$beta$est),ncol(m$Par$beta$beta$est))
      else if(is.momentuHMM(m) & !is.null(m$mle$beta)) m$conditions$betaCons <- matrix(1:length(m$mle$beta),nrow(m$mle$beta),ncol(m$mle$beta))
    }
    if(is.null(m$conditions$betaRef)) m$conditions$betaRef <- as.integer(1:length(m$stateNames))
    if(is.momentuHMM(m)){
      if(is.null(m$mod$wpar)) m$mod$wpar <- m$mod$estimate
      if(is.null(m$mod$Sigma) & !is.null(m$mod$hessian)) m$mod$Sigma <- MASS::ginv(m$mod$hessian)
    } else {
      ####### compatability hack for change to MIcombine in momentuHMM >= 1.4.3 ######
      if(is.null(m$conditions$optInd)){
        for(i in names(m$conditions$dist)){
          m$conditions$workBounds[[i]]<-matrix(c(-Inf,Inf),nrow(m$conditions$workBounds[[i]]),2,byrow=TRUE)
        }
      }
      ################################################################################
    }
    if(is.null(m$conditions$mixtures)) m$conditions$mixtures <- 1
    if(is.null(m$covsPi)) m$covsPi <- matrix(1,length(unique(m$data$ID)),1)
    if(is.null(attr(m$data,"coords")) & !is.null(m$data$x) & !is.null(m$data$y)) attr(m$data,"coords") <- c("x","y")
  } else if(!is.miHMM(m) & any(unlist(lapply(m,is.momentuHMM)))){
    m <- HMMfits(m)
  }
  m
}

checkCrawlArgs <- function(ids,nbAnimals,ind_data,obsData,Time.name,retryFits,attempts,timeStep,mov.model,err.model,activity,drift,theta,fixPar,retrySD,constr,prior,proj,predTime){
  
  if(is.null(mov.model)){
    mov.model<-list()
    for(i in ids){
      mov.model[[i]] <- ~1
    }
  } else if(is.formula(mov.model)){
    tmpmov.model<-mov.model
    mov.model<-list()
    for(i in ids){
      mov.model[[i]] <- tmpmov.model
    }    
  }
  if(!is.null(names(mov.model))) {
    if(!all(names(mov.model) %in% ids)) stop("mov.model names must match obsData$ID")
    mov.model <- mov.model[ids]
  } else names(mov.model) <- ids
  
  if(is.null(err.model)){
    err.model<-vector('list',nbAnimals)
    names(err.model)<-ids
  } else if(is.list(err.model)){
    if(all(unlist(lapply(err.model,is.formula)))){
      tmperr.model<-err.model
      err.model<-list()
      for(i in ids){
        err.model[[i]] <- tmperr.model
      } 
    }
  }
  if(!is.null(names(err.model)))  {
    if(!all(names(err.model) %in% ids)) stop("err.model names must match obsData$ID")
    err.model <- err.model[ids]
  } else names(err.model) <- ids
  
  if(is.null(activity)){
    activity<-vector('list',nbAnimals)
    names(activity)<-ids
  } else if(is.formula(activity)){
    tmpactivity<-activity
    activity<-list()
    for(i in ids){
      activity[[i]] <- tmpactivity
    }    
  }
  if(!is.null(names(activity))) {
    if(!all(names(activity) %in% ids)) stop("activity names must match obsData$ID")
    activity <- activity[ids]
  } else names(activity) <- ids
  
  if(is.null(drift)){
    drift<-list()
    for(i in ids){
      drift[[i]] <- FALSE
    }
  } else if(is.logical(drift)){
    tmpdrift<-drift
    drift<-list()
    for(i in ids){
      drift[[i]] <- tmpdrift
    }    
  }
  if(!is.null(names(drift))) {
    if(!all(names(drift) %in% ids)) stop("drift names must match obsData$ID")
    drift <- drift[ids]
  } else names(drift) <- ids
  
  if(!missing(fixPar)){
    if(!is.list(fixPar)){
      tmpfixPar<-fixPar
      fixPar<-list()
      for(i in ids){
        fixPar[[i]] <- tmpfixPar
      } 
    }
    if(!is.null(names(fixPar))) {
      if(!all(names(fixPar) %in% ids)) stop("fixPar names must match obsData$ID")
      fixPar <- fixPar[ids]
      for(i in ids){
        if(is.null(fixPar[[i]])){
          fixPar[[i]] <- crawl::displayPar(mov.model =  mov.model[[i]], err.model = err.model[[i]], activity = activity[[i]], drift = drift[[i]], data = ind_data[[i]], Time.name = Time.name)$fixPar
        }
      }
    } else names(fixPar) <- ids
  } else {
    fixPar <- list()
    for(i in ids){
      fixPar[[i]] <- crawl::displayPar(mov.model =  mov.model[[i]], err.model = err.model[[i]], activity = activity[[i]], drift = drift[[i]], data = ind_data[[i]], Time.name = Time.name)$fixPar
    } 
  }
  
  if(!missing(theta)){
    if(!is.list(theta)){
      tmptheta<-theta
      theta<-list()
      for(i in ids){
        theta[[i]] <- tmptheta
      } 
    }
    if(!is.null(names(theta))) {
      if(!all(names(theta) %in% ids)) stop("theta names must match obsData$ID")
      theta <- theta[ids]
      for(i in ids){
        if(is.null(theta[[i]])){
          theta[[i]] <- rep(0,sum(is.na(fixPar[[i]])))
        }
      }
    } else names(theta) <- ids
  } else {
    theta <- list()
    for(i in ids){
      theta[[i]] <- rep(0,sum(is.na(fixPar[[i]])))
    }     
  }
  
  if(retryFits>0 | attempts>1){
    if(!is.list(retrySD)){
      tmpretrySD <- retrySD
      retrySD<-list()
      for(i in ids){
        if(length(tmpretrySD)>1){
          if(length(theta[[i]])!=length(tmpretrySD)) stop("retrySD is not the correct length for individual ",i)
          retrySD[[i]] <- tmpretrySD
        } else {
          retrySD[[i]] <- rep(tmpretrySD,length(theta[[i]]))
        }
      }
    } else {
      if(!is.null(names(retrySD))) {
        if(!all(names(retrySD) %in% ids)) stop("retrySD names must match obsData$ID")
        retrySD <- retrySD[ids]
      } else {
        if(length(retrySD)!=nbAnimals) stop('when no list object names are provided, retrySD must be a list of length ',nbAnimals)
        names(retrySD) <- ids
      }
      for(i in ids){
        if(is.null(retrySD[[i]])){
          retrySD[[i]] <- rep(1,length(theta[[i]]))
        } else {
          tmpretrySD <- retrySD[[i]]
          if(length(tmpretrySD)>1){
            if(length(theta[[i]])!=length(tmpretrySD)) stop("retrySD is not the correct length for individual ",i)
          } else {
            retrySD[[i]] <- rep(tmpretrySD,length(theta[[i]]))
          }
        }
      }
    }
    retrySD <- retrySD[ids]
  } else {
    retrySD <- vector('list',nbAnimals)
    names(retrySD) <- ids
    for(i in ids){
      retrySD[[i]] <- rep(1,length(theta[[i]]))
    }
  }
  
  if(is.null(constr)){
    constr<-list()
    for(i in ids){
      constr[[i]]<-list(lower = -Inf,upper = Inf)
    }
  } else if(is.list(constr)){
    if(!is.null(names(constr))){
      if(all(names(constr) %in% c("upper","lower"))){
        tmpconstr<-constr
        constr<-list()
        for(i in ids){
          constr[[i]] <- tmpconstr
        } 
      }
    }
  }
  if(!is.null(names(constr))) {
    if(!all(names(constr) %in% ids)) stop("constr names must match obsData$ID")
    constr <- constr[ids]
  } else names(constr) <- ids
  
  if(is.null(prior)){
    prior<-vector('list',nbAnimals)
    names(prior)<-ids
  } else if(is.function(prior)){
    tmpprior<-prior
    prior<-list()
    for(i in ids){
      prior[[i]] <- tmpprior
    }    
  }
  if(!is.null(names(prior))) {
    if(!all(names(prior) %in% ids)) stop("prior names must match obsData$ID")
    prior <- prior[ids]
  } else names(prior) <- ids
  
  if(is.null(proj)){
    proj<-vector('list',nbAnimals)
    names(proj)<-ids
  } else if(!is.list(proj)){
    tmpproj<-proj
    proj<-list()
    for(i in ids){
      proj[[i]] <- tmpproj
    }    
  }
  if(!is.null(names(proj))) {
    if(!all(names(proj) %in% ids)) stop("proj names must match obsData$ID")
    proj <- proj[ids]
  } else names(proj) <- ids
  
  if(is.null(predTime)){
    predTime<-list()
  }
  if(!is.list(predTime)){
    tmpPredTime<-predTime
    predTime<-list()
    for(i in ids){
      predTime[[i]]<-tmpPredTime
    }
  }
  for(i in ids){
    if(is.null(predTime[[i]])){
      iTime <- range(ind_data[[i]][[Time.name]])
      if(inherits(obsData[[Time.name]],"POSIXct")){
        tzone <- attributes(obsData[[Time.name]])$tzone
        predTime[[i]] <- as.POSIXct(seq(iTime[1],iTime[2],timeStep),tz=tzone)
      } else {
        tzone <- NULL
        predTime[[i]] <- seq(iTime[1],iTime[2],timeStep)
      }
      attributes(predTime[[i]])$tzone <- tzone
    }
  }
  if(!is.null(names(predTime))) {
    if(!all(names(predTime) %in% ids)) stop("predTime names must be character strings that match elements of obsData$ID. See example in ?crawlWrap.")
    predTime <- predTime[ids]
  } else names(predTime) <- ids
  
  return(list(mov.model=mov.model,err.model=err.model,activity=activity,drift=drift,
              theta=theta,fixPar=fixPar,retrySD=retrySD,constr=constr,prior=prior,proj=proj,predTime=predTime))
}
check_fit <- function(mle) {
  if (length(mle) == 1) {
    message('MLE does not exist...')
    return( TRUE )
  } else if (mle$convergence > 0 ) {
    message("MLE exists but DID NOT converge...")
    return( TRUE )
  } 
  message("MLE exists and DID converge. Attempting to solve...")
  C.tmp <- try(2 * solve(mle$hessian), silent = TRUE)
  if (inherits(C.tmp, "try-error")) {
    message("System is exactly singular :(")
    return( TRUE )
  } else {
    message('If diagonals are negative, try again...')
    return( any(diag(C.tmp) <= 0) )
  }
}
crwData <- function(m)
{
  if(is.null(m$crwFits) | is.null(m$crwPredict))
    stop("Can't construct crwData object: fields are missing")
  
  if(any(!unlist(lapply(m$crwFits,function(x) inherits(x,"crwFit")))))
    stop("Can't construct crwData object: crwFits must be a list of crwFit objects")
  
  if(!inherits(m$crwPredict,"crwPredict")) stop("Can't construct crwData object: crwPredict must be a crwPredict object")
  
  obj <- m
  
  class(obj) <- append("crwData",class(obj))
  return(obj)
}

#' Is crwData
#'
#' Check that an object is of class \code{\link{crwData}}. Used in \code{\link{MIfitHMM}}.
#'
#' @param x An R object
#'
#' @return \code{TRUE} if \code{x} is of class \code{\link{crwData}}, \code{FALSE} otherwise.

is.crwData <- function(x)
  inherits(x,"crwData")

is.formula <- function(x)
  tryCatch(inherits(x,"formula"),error= function(e) {FALSE})