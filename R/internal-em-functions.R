## E-STEP FUNCTIONS

## prob of ND given Z
cpNDz <- function(z0, mu, pyfit){
  predict(pyfit, newdata=data.frame(gavg=z0+mu), type="response")
}

## joint prob of ND,Z
pNDZ <- function(z0, mu, s2, pyfit){
  cpNDz(z0, mu, pyfit)*dnorm(z0, 0, sqrt(s2))
}

## marginal prob of ND
pND <- function(mu, s2, pyfit){
  f <- function(z0) pNDZ(z0, mu, s2, pyfit)
  integrate(f, lower=-Inf, upper=Inf)$value
}

## E[Z]
EZ <- function(mu, s2, pyfit){
  f <- function(z0) (z0)*pNDZ(z0, mu, s2, pyfit)/pND(mu, s2, pyfit)
  integrate(f, lower=-Inf, upper=Inf)$value + mu
}

## E[Z^2]
EZ2 <- function(mu, s2, pyfit, EZ){
  f <- function(z0) ((z0)^2)*pNDZ(z0, mu, s2, pyfit)/pND(mu, s2, pyfit)
  integrate(f, lower=-Inf, upper=Inf)$value + 2*mu*EZ - mu^2
}

#######################

## M-STEP FUNCTIONS

updateS2 <- function(Ct, thetaVec, dj, ez, ez2, i.nd, ngene, nperts){
  s2Vec <- vector(length=length(Ct))
  s2Mat <- matrix(nrow=max(ngene), ncol=max(nperts))
  for(i in 1:max(ngene)){
    for(j in 1:max(nperts)){
      ind <- which(ngene == i)
      indD <- intersect(ind, which(!i.nd))
      indND <- intersect(ind, which(i.nd))
      p1 <- sum((Ct[indD]-(thetaVec[indD]+dj[indD]))^2)
      p2 <- sum(ez2[indND]-2*ez[indND]*(thetaVec[indND]+dj[indND])
                +((thetaVec[indND]+dj[indND])^2))
      s2Vec[ind] <- s2Mat[i,j] <- (p1+p2)/length(ind)
    }
  }
  return(list(s2Vec, s2Mat))
}

updateTheta <- function(Y, ez, dj, i.nd, ngene, nperts, Design){
  thetaVec <- vector(length=length(Y))
  thetaMat <- matrix(nrow=max(ngene), ncol=max(nperts))
  Y[which(i.nd==TRUE)]<-ez[which(i.nd==TRUE)]-dj[which(i.nd==TRUE)]
  #DesLM=model.matrix(formula,as.data.frame(nrep))
  fit.lmFit=lmFit(Y, design=Design)

  for(i in 1:max(ngene)){
    for(j in 1:max(nperts)){
      ind2 <- intersect(which(ngene == i), which(nperts == j))
      #thetaVec[ind2] <- thetaMat[i,j] <- mean(Ct[ind2]-dj[ind2])
      #thetaVec[ind2] <- thetaMat[i,j] <- mean(Ct[ind2])
      thetaVec[ind2]<-thetaMat[i,j]<-fit.lmFit$coefficients[i,j]
    }
  }

  return(list(thetaVec, thetaMat, fit.lmFit$sigma, fit.lmFit$cov.coefficients))
}

# updateTheta <- function(Ct, ez, dj, i.nd, ngene, nperts){
#   thetaVec <- vector(length=length(Ct))
#   thetaMat <- matrix(nrow=max(ngene), ncol=max(nperts))
#   for(i in 1:max(ngene)){
#     for(j in 1:max(nperts)){
#       ind <- intersect(which(ngene==i), which(nperts==j))
#       indD <- intersect(ind, which(!i.nd))
#       indND <- intersect(ind, which(i.nd))
#       num <- sum(Ct[indD]-dj[indD])+sum(ez[indND]-dj[indND])
#       denom <- length(ind)
#       thetaVec[ind] <- thetaMat[i,j] <- num/denom
#     }
#   }
#   return(list(thetaVec, thetaMat))
# }


logLik <- function(Ct, ez, ez2, s2Vec, thetaVec, dj, i.nd){
  p1 <- -sum(log(2*pi*s2Vec)/2)
  p2 <- -sum(((Ct[which(!i.nd)]-(thetaVec[which(!i.nd)]+dj[which(!i.nd)]))^2)
             / (2*s2Vec[which(!i.nd)]))
  p3 <- -sum((ez2[which(i.nd)]-2*ez[which(i.nd)]
              *(thetaVec[which(i.nd)]+dj[which(i.nd)])+
              (thetaVec[which(i.nd)]+dj[which(i.nd)])^2)
             / (2*s2Vec[which(i.nd)]))
  p1+p2+p3
}


#########################################
## Multiple imputation
multy <- function(pyfit, numsam, params.new, Ct, Y, dj, ez, ez2,
                  i.nd, ngene, ntype, DesLM, iterMax, tol,
                  vary_fit, vary_model, add_noise)
{
  #set.seed(10291986)
  cat("vary model=",vary_model,"vary_fit=",vary_fit,"add_noise=",add_noise)
  ezInit <- ez
  multyfit<-pyfit
  Sigma<-summary(pyfit)$cov.unscaled
  Betas<-pyfit$coefficients
  multylist <- list()
  for(k in 1:numsam)
  {
   cat("creating data set ",k,"\n")
   ez <- ezInit
   params<-params.new
   if (vary_fit)
   {
    ## varying beta0, beta1 for the fit curve
    multyfit$coefficients<-rmvnorm(1, mean=Betas, sigma=Sigma)

    ll <- vector(length=iterMax)
    iter <- 1
    cond <- TRUE

    ## Update estimates of Theta, sigma2, and missing values
    while(cond){
      print(paste(iter, "/", iterMax))
#      params.old <- params

      ## E-step
      ## Expactation of missing data points
      ez <- rep(NA, length=length(Ct))
      for(i in 1:length(Ct))
      {
        if(i.nd[i])
        {
          xi <- params$thetaVec[i]+dj[i]
          ez[i] <- EZ(xi, params$s2Vec[i], multyfit)
        }
      }

      ## M-Step -- theta
      thetas <- updateTheta(Y, ez, dj, i.nd, ngene, ntype, DesLM)
      params$thetaVec   <- thetas[[1]]
      params$thetaMat   <- thetas[[2]]
      params$sigma      <- thetas[[3]]
      params$cov.matrix <- thetas[[4]]

      ## E-Step
      ## missing values
      ez <- ez2 <- rep(NA, length=length(Ct))
      for(i in 1:length(Ct)){
        if(i.nd[i]){
          xi <- params$thetaVec[i]+dj[i]
          ez[i] <- EZ(xi, params$s2Vec[i], multyfit)
          ez2[i] <- EZ2(xi, params$s2Vec[i], multyfit, ez[i])
        }
      }

      ## M-Step -- sigma2
      sigmas <- updateS2(Ct, params$thetaVec, dj, ez, ez2,
                         i.nd, ngene, ntype)
      params$s2Vec <- sigmas[[1]]
      params$s2Mat <- sigmas[[2]]

      ## log-likelihoog
      ll[iter] <- logLik(Ct, ez, ez2, params$s2Vec,
                         params$thetaVec, dj, i.nd)
      message(ll[iter])

      if(iter>1) cond <- (abs(ll[iter]-ll[iter-1]) > tol) & (iter < iterMax)
      iter <- iter+1
    }
   }

   if (vary_model)
   {
    ## Valying the Thetas
    error.sd <- params.new$sigma # use sigma from EM
    #    error.sd <- params$sigma # use sigma from updated fit
    j <- which(!is.na(ez))
    gtmp <- unique(ngene[j])
    ## Generate theta vector for genes with missing values
    for(i in 1:length(gtmp))
      {
        VarMat <- (error.sd[gtmp[i]])^2*params$cov.matrix
        params$thetaMat[gtmp[i],] <- rmvnorm(1, mean=params$thetaMat[gtmp[i],], sigma=VarMat)
      }
    for(i in 1:max(ngene))
      {
       for(j in 1:max(ntype))
        {
         ind2 <- intersect(which(ngene == i), which(ntype == j))
         params$thetaVec[ind2]<-params$thetaMat[i,j]
        }
      }
    ez <- rep(NA, length=length(Ct))
    ## Update E(z) based on new thetas
    for(i in 1:length(Ct))
     {
      if(i.nd[i])
       {
        xi <- params$thetaVec[i]+dj[i]
        ez[i] <- EZ(xi, params$s2Vec[i], multyfit)
       }
     }
   }


if (add_noise)
{
  ## Adding random noise
  epsilonVec <- as.vector(c(rep(0 ,length(Ct))))
  epsilonMat <- matrix(0, nrow=max(ngene), ncol=max(ntype))
  error.sd <- params.new$sigma # use sigma from EM
  #    error.sd <- params$sigma # use sigma from updated fit

  j <- which(!is.na(ez))
  gtmp <- unique(ngene[j])
  ## Generate epsilon vector for genes with missing values
  for(i in 1:length(gtmp))
  {
    epsilonMat[gtmp[i],] <- rnorm(max(ntype), 0, error.sd[gtmp[i]])
  }
  for(i in 1:max(ngene))
  {
    for(j in 1:max(ntype))
    {
      ind2 <- intersect(which(ngene == i), which(ntype == j))
      epsilonVec[ind2]<-epsilonMat[i,j]
    }
  }

  ## Update ez based on new errors
  for(i in 1:length(Ct))
  {
    if(i.nd[i])
    {
     # xi <- params$thetaVec[i]+dj[i]
      ez[i] <- ez[i] + epsilonVec[i]
    }
  }

}
    ## Crerate objects with Imputed values
    Ct[which(as.logical(i.nd))] <- ez[which(as.logical(i.nd))]
    ind <- grep("target", featureType(object), ignore.case=TRUE)
    exprs(object)[ind,] <- Ct

    fc <- as.matrix(featureCategory(object))
    fc[which(fc == "Undetermined", arr.ind=TRUE)] <- "Imputed"
    featureCategory(object) <- as.data.frame(fc)

    if (nrow(getCtHistory(object)) == 0){
      setCtHistory(object) <- data.frame(
        history = "Manually created qPCRset object.",
        stringsAsFactors = FALSE)
    }
    setCtHistory(object) <- rbind(getCtHistory(object),
                                  capture.output(match.call(multy)))
    ## Save objects to the list
    multylist[[k]]<-object
  }
multylist[[k+1]]<-params.new
multylist[[k+2]]<-params
return(multylist)
}
