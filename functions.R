#'
#' Code for computing covariate-adjusted treatment effect estimators
#'
#' Provided by Elizabeth Colantuoni
#'
#' The key functions in this file are "unadjusted", "PLEASE", "Rotnitzky", "run_analysis"
#' There are many convenience/helper functions, unused test functions, and other related functions included as well.


# Write the unadjusted estimator 
unadjusted <- function(y,trt) {
  list(mean(y[trt==1],na.rm=T),mean(y[trt==0],na.rm=T))
}

# Write the code for the IPTW estimator
IPTW <- function(y,trt,piX,ctrl=list(epsilon=1e-7,maxit=50)){
  if( !all (trt %in% c(0,1) ) ) stop( "'TRT' should be binary." )
  piFIT <- try(glm(trt~as.matrix(piX),family=quasibinomial,control=ctrl,na.action=na.omit),TRUE)
  failed.pi <- inherits(piFIT,"try-error")
  if(failed.pi) piFIT <- suppressWarnings(glm(trt~as.matrix(piX),family=quasibinomial,control=ctrl,na.action=na.omit))
  mu.hat.1 <- sum(y*trt/piFIT$fitted)/sum(trt/piFIT$fitted)
  mu.hat.0 <- sum(y*(1-trt)/(1-piFIT$fitted))/sum((1-trt)/(1-piFIT$fitted))
  list(mu.hat.1,mu.hat.0)
}




# Write the code for the Joffe
Joffe.calcs <- function(y,trt,piX,phiX,ss=1,famY,ctrl){
  if(ss==1) trtnew=trt else trtnew = (1-trt)
  piFIT <- try(glm(trtnew~as.matrix(piX),family=binomial,control=ctrl,na.action=na.omit),TRUE)
  failed.pi <- inherits(piFIT,"try-error")
  if(failed.pi) piFIT <- supressWarnings(glm(trtnew~as.matrix(piX),family=binomial,control=ctrl,na.action=na.omit))
  wt.alpha<-piFIT$fit
  check.wt<-max(1/piFIT$fit[trt==ss])
  ave.wt<-mean(trtnew/piFIT$fit)
  phiFIT <- try(glm(y~as.matrix(phiX),family=famY,weights=1/piFIT$fit,subset=trt==ss,control=ctrl,na.action=na.omit),TRUE)
  failed.phi <- inherits(phiFIT,"try-error")
  if(failed.phi) phiFIT <- supressWarnings(glm(y~as.matrix(phiX),family=famY,weights=1/piFIT$fit,subset=trt==ss,control=ctrl,na.action=na.omit))
  Q0 <- predict(phiFIT,as.data.frame(phiX),type="response")
  mu.hat <- mean(Q0)
  list(wt.alpha=wt.alpha,Q0=Q0,mu.hat=mu.hat,check.wt=check.wt,ave.wt=ave.wt,failed.pi=failed.pi,failed.phi=failed.phi)
}


Joffe <- function(y,trt,piX,phiX,ctrl=list(epsilon=1e-8,maxit=50)) {
  if( !all (trt %in% c(0,1) ) ) stop( "'TRT' should be binary." )
  YisBin <- all(y[!is.na(y)] %in% c(0,1))
  famY <- ifelse(YisBin,"quasibinomial","gaussian")
  trt1 <- Joffe.calcs(y,trt,piX,phiX,ss=1,famY,ctrl)
  trt0 <- Joffe.calcs(y,trt,piX,phiX,ss=0,famY,ctrl)
  mu.hat.1 <- trt1$mu.hat
  mu.hat.0 <- trt0$mu.hat
  list(mu.hat.1,mu.hat.0,trt1$check.wt,trt0$check.wt,trt1$ave.wt,trt0$ave.wt,trt1$failed.pi,trt1$failed.phi,trt0$failed.pi,trt1$failed.phi)
}

# Write the code for the Rotnitzky estimator

# Write the code for the Rotnitzky estimator

Rotnitzky.calcs <- function(y,trt,piX,phiX,initEtacontrol=NULL,initEtatrt=NULL,ctrl=list(epsilon=1e-8,maxit=50),famY,YisBin){
  #delta <- ifelse(is.na(y),0,1)
  #y <- y*delta
  # if(ss==1) trtnew = trt*delta
  # if(ss==0) trtnew = (1-trt)*delta
  options(warn=2)
  piFIT <- try(glm(trt ~ as.matrix(piX),family=binomial,control=ctrl,na.action=na.omit),TRUE)
  failed.pi <- inherits(piFIT,"try-error")
  if(failed.pi) piFIT <- suppressWarnings(glm(trt ~ as.matrix(piX),family=binomial,control=ctrl,na.action=na.omit))
  piMLE <- fitted(piFIT)
  trtdPi <- trt/piMLE
  trtdPi[ trt==0 ] <- 0
  ytrtdPi <- y * trtdPi
  ytrtdPi[ trt==0] <- 0
  controldPi <- (1-trt)/(1-piMLE) # check this
  controldPi[ (1-trt)==0 ] <- 0
  ycontroldPi <- y * controldPi
  ycontroldPi[ (1-trt)==0] <- 0
  
  phiX1 <- cbind(1,phiX)
  if( is.vector(phiX) ) nEta <- 2  else nEta <- dim(phiX)[[2]]+1
  if( is.vector(piX) ) nAlph <- 2  else nAlph <- dim(piX)[[2]]+1
  Salph <- (trt - piMLE)* cbind(1,piX)
  piXtpiX <- t( apply( cbind(1,piX),1,function(x) outer(x,x)))
  VSalph <- matrix(colMeans( piMLE*(1-piMLE)*piXtpiX), nAlph)
  ViSalph <- qr.solve(VSalph, t(Salph))

    fitY <- function(eta) {
      etaX <- c( cbind(1, as.matrix(phiX)) %*% eta)
      if(!YisBin) return(etaX)
      return( c(1/(1+exp(-etaX))))
    }

    if(is.null(initEtatrt)) initEtatrt <- coef(suppressWarnings(glm(y~as.matrix(phiX),subset=(trt==1),family=famY,na.action=na.omit)))  else {
      if(length(initEtatrt) != nEta) stop("'initEtatrt' should be a vector of length ", nEta)
      if(!is.numeric(initEtatrt)) stop("'initEtatrt' should be a numeric vector")}
    if(is.null(initEtacontrol)) initEtacontrol <- coef(suppressWarnings(glm(y~as.matrix(phiX),subset=(trt==0),family=famY,na.action=na.omit)))  else {
      if(length(initEtacontrol) != nEta) stop("'initEtacontrol' should be a vector of length ", nEta)
      if(!is.numeric(initEtacontrol)) stop("'initEtacontrol' should be a numeric vector")
	}    
	initEta <- c(initEtacontrol,initEtatrt)
	
  # Step 1: Function "bounded" returns nu_OR, y_hat and beta_OR
  bounded <- function(Pi){
    phiFITtrt <- try(glm(y~as.matrix(phiX),family=famY,weights=(trt/Pi),control=ctrl,na.action=na.omit),TRUE)
    failed.phitrt <- inherits(phiFITtrt,"try-error")
    if(failed.phitrt) phiFITtrt <- suppressWarnings(glm(y~as.matrix(phiX),family=famY,weights=(trt/Pi),control=ctrl,na.action=na.omit))
    etaBtrt <- phiFITtrt$coef
    phiBtrt <- fitY(etaBtrt)
    betaBtrt <- mean(phiBtrt)

    phiFITcontrol <- try(glm(y~as.matrix(phiX),family=famY,weights=(1-trt)/(1-Pi),control=ctrl,na.action=na.omit),TRUE)
    failed.phicontrol <- inherits(phiFITcontrol,"try-error")
    if(failed.phicontrol) phiFITcontrol <- suppressWarnings(glm(y~as.matrix(phiX),family=famY,weights=(1-trt)/(1-Pi),control=ctrl,na.action=na.omit))
    etaBcontrol <- phiFITcontrol$coef
    phiBcontrol <- fitY(etaBcontrol)
    betaBcontrol <- mean(phiBcontrol)
	list(etaBtrt=etaBtrt,etaBcontrol=etaBcontrol,phiBtrt=phiBtrt,phiBcontrol=phiBcontrol,betaBcontrol=betaBcontrol,betaBtrt=betaBtrt,failed.phitrt=failed.phitrt,failed.phicontrol=failed.phicontrol)
  }

  OR <- bounded(piMLE)
  phiOR <- list(control=OR$phiBcontrol,trt=OR$phiBtrt)
  betaOR <- list(control=OR$betaBcontrol,trt=OR$betaBtrt)
  failed.phi <- list(control=OR$failed.phicontrol,trt=OR$failed.phitrt)
  
  
  # Step 2: Calculate nu_t minimizing the estimated variance of betaOR
  eff <- function(initEta,BETA) {
    effObj <- function(eta, BETA) {
	  etacontrol <- eta[1:nEta]
	  etatrt <- eta[(nEta+1):(2*nEta)]
      phihatcontrol <- fitY(etacontrol)
      psicontrol <- ycontroldPi - BETA$control - (controldPi - 1)*phihatcontrol
      covarcontrol <- colMeans(psicontrol * Salph)
      
      phihattrt <- fitY(etatrt)
      psitrt <- ytrtdPi - BETA$trt - (trtdPi - 1)*phihattrt
      covartrt <- colMeans(psitrt * Salph)

      sum((
      (-1/mean(controldPi))*(psicontrol-c(covarcontrol %*% ViSalph))
      +(1/mean(trtdPi))*(psitrt-c(covartrt %*% ViSalph))
      )^2)
    }
    effgr <- function(eta,BETA){
	  etacontrol <- eta[1:nEta]
	  etatrt <- eta[(nEta+1):(2*nEta)]

      phihatcontrol <- fitY(etacontrol)
      psicontrol <- ycontroldPi - BETA$control - (controldPi - 1)*phihatcontrol
      covarcontrol <- colMeans(psicontrol * Salph)
      CFcontrol = 2*(psicontrol - c(covarcontrol %*% ViSalph))
      EXcontrol = c(exp(-as.matrix(cbind(1,phiX)) %*% etacontrol))
      dEXcontrol = 1 / (1/EXcontrol + 2 + EXcontrol)
      gr1control <- - (controldPi - 1) * dEXcontrol * phiX1
      prjcontrol <- function(i) colMeans(-(controldPi-1)*dEXcontrol*phiX1[,i]*Salph)%*%ViSalph
      gr2control <- sapply(1:nEta,prjcontrol)
           
      phihattrt <- fitY(etatrt)
      psitrt <- ytrtdPi - BETA$trt - (trtdPi - 1)*phihattrt
      covartrt <- colMeans(psitrt * Salph)
      CFtrt = 2*(psitrt - c(covartrt %*% ViSalph))
      EXtrt = c(exp(-as.matrix(cbind(1,phiX)) %*% etatrt))
      dEXtrt = 1 / (1/EXtrt + 2 + EXtrt)
      gr1trt <- - (trtdPi - 1) * dEXtrt * phiX1
      prjtrt <- function(i) colMeans(-(trtdPi-1)*dEXtrt*phiX1[,i]*Salph)%*%ViSalph
      gr2trt <- sapply(1:nEta,prjtrt)

      c(colSums( ((-1/mean(controldPi))*CFcontrol+(1/mean(trtdPi))*CFtrt) * (-1/mean(controldPi))*(gr1control +gr2control)),
      colSums( ((-1/mean(controldPi))*CFcontrol+(1/mean(trtdPi))*CFtrt) * (1/mean(trtdPi))*(gr1trt+gr2trt)))
    }
    ### NOTE: I need to add in the flag for convergence for this optimization
    optim.out <- optim(initEta, effObj, gr=effgr, BETA=BETA, method="BFGS")

    return(c(optim.out$par,optim.out$convergence))
  }

  #effcont <- function(BETA){
  #  B <-  (trtdPI -1) * phiX1
   # newY <- ytrtdPi - BETA
    #newY <- newY - c(colMeans(newY*Salph) %*% ViSalp)
    #newX2 <- sapply(1:nEta,function(i) colMeans(B[,i]*Salph)%*%ViSalph)
    #newX <- B - newX2
    #coef( glm(newY~0 + newX))
  #}

    if(YisBin) {
      out <- eff(initEta, betaOR)
      etaE <- out[1:(length(out)-1)]
      optim.conv <- out[length(out)]
    }
    else {return(0); print("Binary outcome required")
      #etaE <- effcont(betaOR)
      #optim.conv <- 0
    }
   
  # Step 3: fit the extended model for treatment assignment
  # and repeat step 1 with fitted treatment assignment
  # probability from extended treatment assignment model

	etaEcontrol <- etaE[1:nEta]
 	etaEtrt <- etaE[(nEta+1):(2*nEta)]

  phiE <- list(control=fitY(etaEcontrol),trt= fitY(etaEtrt))
  betaE <- list(control=mean(ycontroldPi - (controldPi -1)*phiE$control),trt=mean(ytrtdPi - (trtdPi -1)*phiE$trt))
  V1 <- 1/(1-piMLE) * (phiOR$control - betaOR$control)
  V2 <- 1/(1-piMLE) * (phiE$control - betaE$control)
  V3 <- 1/piMLE * (phiOR$trt - betaOR$trt)
  V4 <- 1/piMLE * (phiE$trt - betaE$trt)
  extX <- cbind(piX,V1,V2,V3,V4)
  write.table(cbind(y,trt,extX),file="test.csv",sep=",",row.names=FALSE)
  piFIText <- try(glm(trt~as.matrix(extX),family=binomial,control=ctrl,na.action=na.omit),TRUE)
  failed.piDR <- inherits(piFIText,"try-error")
  if(failed.piDR) piFIText <- piFIT
  #piFIText <- suppressWarnings(glm(trt~as.matrix(extX),family=binomial,control=ctrl,na.action=na.omit))
  piext <- piFIText$fitted
  DR <- bounded(piext)
  betaDR <- list(control=DR$betaBcontrol,trt=DR$betaBtrt)
  etaDR <- c(DR$etaBcontrol,DR$etaBtrt)
  failed.phiDR <- list(control=DR$failed.phicontrol,trt=DR$failed.phitrt)

  max.wt.ext<-max(max(1/piext),max(1/(1-piext)))
  ave.wt.ext<-mean(c(mean(trt/piext),mean((1-trt)/(1-piext))))
  list(betaDR=betaDR,etaDR=etaDR,max.wt.ext=max.wt.ext,
       ave.wt.ext=ave.wt.ext,failed.pi=failed.pi,failed.phi=failed.phi,failed.piDR=failed.piDR,failed.phiDR=failed.phiDR,optim.conv=optim.conv)
}

Rotnitzky <- function(y,trt,piX,phiX,initEtacontrol=NULL,initEtatrt=NULL,ctrl=list(epsilon=1e-7,maxit=50)) {
  if( !all (trt %in% c(0,1) ) ) stop( "'TRT' should be binary." )
  YisBin <- all(y[!is.na(y)] %in% c(0,1))
  famY <- ifelse(YisBin,"quasibinomial","gaussian")
  rout <- Rotnitzky.calcs(y,trt,piX,phiX,initEtacontrol,initEtatrt,ctrl,famY,YisBin)
  list(mu.hat.1=rout$betaDR$trt,mu.hat.0=rout$betaDR$control,
  	   rout$max.wt.ext,rout$ave.wt.ext,rout$failed.pi,rout$failed.phi$trt,rout$failed.piDR,rout$failed.phiDR$trt,
       rout$failed.phi$control,rout$failed.piDR,rout$failed.phiDR$control,
       rout$optim.conv)
}

Rotnitzky.calcs.old <- function(y,trt,piX,phiX,initEta=NULL,ctrl=list(epsilon=1e-8,maxit=50),ss=1,famY,YisBin){
  delta <- ifelse(is.na(y),0,1)
  y <- y*delta
  if(ss==1) trtnew = trt*delta
  if(ss==0) trtnew = (1-trt)*delta
  options(warn=2)
  piFIT <- try(glm(trtnew ~ as.matrix(piX),family=quasibinomial,control=ctrl,na.action=na.omit),TRUE)
  failed.pi <- inherits(piFIT,"try-error")
  if(failed.pi) piFIT <- suppressWarnings(glm(trtnew ~ as.matrix(piX),family=quasibinomial,control=ctrl,na.action=na.omit))
  piMLE <- fitted(piFIT)
  trtdPi <- trtnew/piMLE
  trtdPi[ trtnew==0 ] <- 0
  ytrtdPi <- y * trtdPi
  ytrtdPi[ trtnew==0] <- 0
  phiX1 <- cbind(1,phiX)
  if( is.vector(phiX) ) nEta <- 2  else nEta <- dim(phiX)[[2]]+1
  if( is.vector(piX) ) nAlph <- 2  else nAlph <- dim(piX)[[2]]+1
  Salph <- (trtnew - piMLE)* cbind(1,piX)
  piXtpiX <- t( apply( cbind(1,piX),1,function(x) outer(x,x)))
  VSalph <- matrix(colMeans( piMLE*(1-piMLE)*piXtpiX), nAlph)
  ViSalph <- qr.solve(VSalph, t(Salph))

    fitY <- function(eta) {
      etaX <- c( cbind(1, as.matrix(phiX)) %*% eta)
      if(!YisBin) return(etaX)
      return( c(1/(1+exp(-etaX))))
    }

    if(is.null(initEta)) initEta <- coef(suppressWarnings(glm(y~as.matrix(phiX),family=famY,subset=trt==ss,na.action=na.omit)))  else {
      if(length(initEta) != nEta) stop("'initEta' should be a vector of length ", nEta)
      if(!is.numeric(initEta)) stop("'initEta' should be a numeric vector")
    }    

  # Step 1: Function "bounded" returns nu_OR, y_hat and beta_OR
  bounded <- function(Pi){
    phiFIT <- try(glm(y~as.matrix(phiX),family=famY,weights=1/Pi,subset=trt==ss,control=ctrl,na.action=na.omit),TRUE)
    failed.phi <- inherits(phiFIT,"try-error")
    if(failed.phi) phiFIT <- suppressWarnings(glm(y~as.matrix(phiX),family=famY,weights=1/Pi,subset=trt==ss,control=ctrl,na.action=na.omit))
    etaB <- phiFIT$coef
    phiB <- fitY(etaB)
    betaB <- mean(phiB)
    list(etaB=etaB,phiB=phiB,betaB=betaB,failed.phi=failed.phi)
  }

  OR <- bounded(piMLE)
  phiOR <- OR$phiB
  betaOR <- OR$betaB
  failed.phi <- OR$failed.phi
  
  # Step 2: Calculate nu_t minimizing the estimated variance of betaOR
  eff <- function(initEta,BETA) {
    effObj <- function(eta, BETA) {
      phihat <- fitY(eta)
      psi <- ytrtdPi - BETA - (trtdPi - 1)*phihat
      covar <- colMeans(psi * Salph)
      sum((psi-c(covar %*% ViSalph))^2)
    }
    effgr <- function(eta,BETA){
      phihat <- fitY(eta)
      psi <- ytrtdPi - BETA - (trtdPi - 1)*phihat
      covar <- colMeans(psi * Salph)
      CF = 2*(psi - c(covar %*% ViSalph))
      EX = c(exp(-as.matrix(cbind(1,phiX)) %*% eta))
      dEX = 1 / (1/EX + 2 + EX)
      gr1 <- - (trtdPi - 1) * dEX * phiX1
      prj <- function(i) colMeans(-(trtdPi-1)*dEX*phiX1[,i]*Salph)%*%ViSalph
      gr2 <- sapply(1:nEta,prj)
      colSums( CF * (gr1 + gr2))
    }
    ### NOTE: I need to add in the flag for convergence for this optimization
    optim.out <- optim(initEta, effObj, gr=effgr, BETA=BETA, method="BFGS")
    return(c(optim.out$par,optim.out$convergence))
  }

  effcont <- function(BETA){
    B <-  (trtdPI -1) * phiX1
    newY <- ytrtdPi - BETA
    newY <- newY - c(colMeans(newY*Salph) %*% ViSalp)
    newX2 <- sapply(1:nEta,function(i) colMeans(B[,i]*Salph)%*%ViSalph)
    newX <- B - newX2
    coef( glm(newY~0 + newX))
  }

    if(YisBin) {
      out <- eff(initEta, betaOR)
      etaE <- out[1:(length(out)-1)]
      optim.conv <- out[length(out)]
    }
    else {
      etaE <- effcont(betaOR)
      optim.conv <- 0
    }
   
    
  # Step 3: fit the extended model for treatment assignment
  # and repeat step 1 with fitted treatment assignment
  # probability from extended treatment assignment model

  phiE <- fitY(etaE)
  betaE <- mean(ytrtdPi - (trtdPi -1)*phiE)
  V1 <- 1/piMLE * (phiOR - betaOR)
  V2 <- 1/piMLE * (phiE - betaE)
  extX <- cbind(piX,V1,V2)
  piFIText <- try(glm(trtnew~as.matrix(extX),family=quasibinomial),TRUE)
  failed.piDR <- inherits(piFIText,"try-error")
  if(failed.piDR) piFIText <- suppressWarnings(glm(trtnew~as.matrix(extX),family=quasibinomial))
  piext <- piFIText$fitted
  DR <- bounded(piext)
  betaDR <- DR$betaB
  etaDR <- DR$etaB
  failed.phiDR <- DR$failed.phi
  max.wt.ext<-max(1/piext[trt==ss])
  ave.wt.ext<-mean(trtnew/piext)
  list(betaDR=betaDR,etaDR=etaDR,max.wt.ext=max.wt.ext,
       ave.wt.ext=ave.wt.ext,failed.pi=failed.pi,failed.phi=failed.phi,failed.piDR=failed.piDR,failed.phiDR=failed.phiDR,optim.conv=optim.conv)
}

Rotnitzky.old <- function(y,trt,piX,phiX,initEta=NULL,ctrl=list(epsilon=1e-7,maxit=50)) {
  if( !all (trt %in% c(0,1) ) ) stop( "'TRT' should be binary." )
  YisBin <- all(y[!is.na(y)] %in% c(0,1))
  famY <- ifelse(YisBin,"quasibinomial","gaussian")
  trt1 <- Rotnitzky.calcs(y,trt,piX,phiX,initEta,ctrl,ss=1,famY,YisBin)
  trt0 <- Rotnitzky.calcs(y,trt,piX,phiX,initEta,ctrl,ss=0,famY,YisBin)
  mu.hat.1 <- trt1$betaDR
  mu.hat.0 <- trt0$betaDR
  list(mu.hat.1,mu.hat.0,trt1$max.wt,trt1$max.wt.ext,
       trt0$max.wt,trt0$max.wt.ext,trt1$ave.wt.ext,trt0$ave.wt.ext,
       trt1$failed.pi,trt1$failed.phi,trt1$failed.piDR,trt1$failed.phiDR,
       trt0$failed.pi,trt0$failed.phi,trt0$failed.piDR,trt0$failed.phiDR,
       trt1$optim.conv,trt0$optim.conv)
}

## Gruber estimator

bound <- function(x,bounds) {
  x[x<min(bounds)] <- min(bounds)
  x[x>max(bounds)] <- max(bounds)
  return(x)
}

tmle <- function(y,trt,g1W,Q,f,gbds=c(10^-8,1-10^-8)) {
  g1W <- bound(g1W,gbds)
  eps1 <- eps2 <- Inf
  epsilon <- 0.00001
  maxIter <- 30
  iterations <- 0
  options(warn=2)
  while((any(abs(c(eps1,eps2))>epsilon)) & iterations <= maxIter){
    iterations <- iterations + 1
    h <- cbind(trt/g1W - (1-trt)/(1-g1W),1/g1W,-1/(1-g1W))
    m <- try(glm(y~-1+offset(Q[,"QAW"]) + h[,1],family=quasibinomial,na.action=na.omit),TRUE)
    failed.m <- inherits(m,"try-error")
    if(failed.m) m <- suppressWarnings(glm(y~-1+offset(Q[,"QAW"]) + h[,1],family=quasibinomial,na.action=na.omit))
    eps1 <- coef(m)
    Q <- Q + eps1*h
    h2 <- plogis(Q[,"Q1W"])/g1W + plogis(Q[,"Q0W"])/(1-g1W)
    h3 <- f/(g1W * (1-g1W))
    g <- try(glm(trt~-1 + offset(qlogis(g1W)) + h2 + h3,family=quasibinomial,na.action=na.omit),TRUE)
    failed.g <- inherits(g,"try-error")
    if(failed.g) g <- suppressWarnings(glm(trt~-1 + offset(qlogis(g1W)) + h2 + h3,family=quasibinomial,na.action=na.omit))
    g1W <- bound(predict(g,type="response"),gbds)
    eps2 <- coef(g)
  }
    Q <- as.data.frame(plogis(as.matrix(Q)))
    names(Q) <- c("QAW","Q0W","Q1W")
    psi.en <- mean(Q[,"Q1W"]-Q[,"Q0W"])
    psi.IPTW <- mean((trt/g1W - (1-trt)/(1-g1W))*y)
    psi.AIPTWQstargstar <- mean((trt/g1W - (1-trt)/(1-g1W))*y - (Q[,"Q1W"]/g1W - Q[,"Q0W"]/(1-g1W))*(trt-g1W))
    psi.AIPTWQegstar <- mean((trt/g1W - (1-trt)/(1-g1W))*y - f/(g1W*(1-g1W))*(trt-g1W))
    mu.hat.1 <- mean(trt/g1W*y - f/(g1W*(1-g1W))*(trt-g1W))
    mu.hat.0 <- mean((1-trt)/(1-g1W)*y - f/(g1W*(1-g1W))*(trt-g1W))
  return(c(mu.hat.1,mu.hat.0,psi.en,psi.IPTW,psi.AIPTWQstargstar,psi.AIPTWQegstar,failed.m,failed.g,iterations))
}

ff <- function(x,y,wt,W){
  beta.W <- NULL
  for(i in 1:ncol(as.matrix(W))) beta.W <- cbind(beta.W,x[i+2]*as.matrix(W)[,i])
  sum(wt^2*(y-(x[1]+exp(x[2]+beta.W)/(1+exp(x[2]+beta.W))))^2)
}

Gruber <- function(y,trt,piX,phiX,ctrl=list(epsilon=1e-8,maxit=50)) {
  if( !all (trt %in% c(0,1) ) ) stop( "'TRT' should be binary." )
  YisBin <- all(y[!is.na(y)] %in% c(0,1))
  famY <- ifelse(YisBin,"quasibinomial","gaussian")
  options(warn=2)
  piFIT <- try(glm(trt~as.matrix(piX),family=quasibinomial,control=ctrl,na.action=na.omit),TRUE)
  failed.pi <- inherits(piFIT,"try-error")
  if(failed.pi) piFIT <- suppressWarnings(glm(trt~as.matrix(piX),family=quasibinomial,control=ctrl,na.action=na.omit))
  wt = (2*trt-1)/bound(piFIT$fitted,c(10^-8,1-10^-8))
  phiFIT <- try(glm(y~as.matrix(cbind(trt,phiX)),family=famY,control=ctrl,na.action=na.omit),TRUE)
  failed.phi <- inherits(phiFIT,"try-error")
  if(failed.phi) phiFIT <- suppressWarnings(glm(y~as.matrix(cbind(trt,phiX)),family=famY,control=ctrl,na.action=na.omit))
  Q <- as.data.frame(cbind(predict(phiFIT,as.data.frame(cbind(trt,phiX)),type="link"),
               predict(phiFIT,as.data.frame(cbind(trt=rep(0,length(trt)),phiX)),type="link"),
               predict(phiFIT,as.data.frame(cbind(trt=rep(1,length(trt)),phiX)),type="link")))
  names(Q) = c("QAW","Q0W","Q1W")
  ## Now get the solution to "f"
  get.f <- try(nlm(ff,rep(0,(ncol(as.matrix(phiX))+2)),y,wt,phiX),TRUE)
  failed.f <- inherits(get.f,"try-error")
  if(!failed.f) {
        nlm.code <- get.f$code
        expit <- function(x) exp(x) / (1+exp(x))
        beta.W <- apply(sweep(as.matrix(phiX),2,get.f$estimate[3:length(get.f$estimate)],"*"),1,sum)
        ## Replace the NaN with 1 since the issue is the large value within the expit
        test <- expit(get.f$estimate[2]+beta.W)
        #print(summary(test))
        test[is.nan(test)] <- 1
        f <- get.f$estimate[1] + test
        grub.out <- tmle(y,trt,piFIT$fitted,Q,f,gbds=c(10^-9,1-10^-9))
      }
      if(failed.f) {grub.out <- rep(NA,9);nlm.code<-NA}
    return(c(grub.out,failed.pi,failed.phi,failed.f,nlm.code))
}

TAN <- function(y,trt,piX,phiX,ctrl=list(epsilon=1e-8,maxit=50)) {
  if( !all (trt %in% c(0,1) ) ) stop( "'TRT' should be binary." )
  YisBin <- all(y[!is.na(y)] %in% c(0,1))
  famY <- ifelse(YisBin,"quasibinomial","gaussian")
  piFIT <- try(glm(trt~as.matrix(piX),family=binomial,control=ctrl,na.action=na.omit),TRUE)
  failed.pi <- inherits(piFIT,"try-error")
  if(failed.pi) piFIT <- supressWarnings(glm(trt~as.matrix(piX),family=binomial,control=ctrl))
  X <- model.matrix(piFIT)
  phiFIT <- try(glm(y~as.matrix(phiX),family=famY,subset=trt==1,control=ctrl,na.action=na.omit),TRUE)
  failed.phi1 <- inherits(phiFIT,"try-error")
  if(failed.phi1) phiFIT <- supressWarnings(glm(y~as.matrix(phiX),family=famY,subset=trt==1,control=ctrl,na.action=na.omit))
  Q0.1 <- predict(phiFIT,as.data.frame(phiX),type="response")
  phiFIT <- try(glm(y~as.matrix(phiX),family=famY,subset=trt==0,control=ctrl,na.action=na.omit),TRUE)
  failed.phi0 <- inherits(phiFIT,"try-error")
  if(failed.phi0) phiFIT <- supressWarnings(glm(y~as.matrix(phiX),family=famY,subset=trt==0,control=ctrl,na.action=na.omit))
  Q0.0 <- predict(phiFIT,as.data.frame(phiX),type="response")
  out <- try(ate.clik(y,trt,piFIT$fitted,g0=cbind(1,Q0.0),g1=cbind(1,Q0.1),X,inv="ginv"),TRUE)
  failed.tan <- inherits(out,"try-error")
  if(failed.tan) mu.hat.1 <- mu.hat.0 <- v.diff <- conv1 <- conv0 <- NA
  if(!failed.tan) {
    mu.hat.1 <- out$mu[1]
    mu.hat.0 <- out$mu[2]
    v.diff <- out$v.diff
    conv1 <-out$conv[1]
    conv0 <- out$conv[2]
  }
  list(mu.hat.1=mu.hat.1,mu.hat.0=mu.hat.0,v.diff=v.diff,conv1=conv1,conv0=conv0,failed.pi=failed.pi,failed.phi1=failed.phi1,failed.phi0=failed.phi0)
}

# PLEASE.calcs computes the estimated mean for a given treatment group
# This function takes the outcome (y), treatment indicator (trt),
# matrix of covariates defining the missing data model (dX),
# matrix of baseline variables defining the propensity score model (piX),
# matrix of baseline variables defining the outcome regression model (phiX),
# treatment group value (ss), and distribution of Y
PLEASE.calcs <- function(y,trt,dX,piX,phiX,ss=1,famY,ctrl){
  if(ss==1) trtnew=trt else trtnew=(1-trt)
  delta <- ifelse(is.na(y),0,1)
  if(sum(delta)==length(y)) h.wt <- rep(1,length(y))
  if(sum(delta)<length(y)) {
      hFIT <- try(glm(delta~as.matrix(cbind(trt,dX)),family=binomial,control=ctrl,na.action=na.omit),TRUE)
      failed.h <- inherits(hFIT,"try-error")
      if(failed.h) hFIT <- suppressWarnings(glm(delta~as.matrix(cbind(trt,dX)),family=binomial,control=ctrl,na.action=na.omit))
      h.wt <- hFIT$fitted           
  }
  piFIT <- try(glm(trtnew~as.matrix(piX),family=binomial,control=ctrl,na.action=na.omit),TRUE)
  failed.pi <- inherits(piFIT,"try-error")
  if(failed.pi) piFIT <- suppressWarnings(glm(trtnew~as.matrix(piX)[,1:ncol(phiX)],family=binomial,control=ctrl,na.action=na.omit))
  wt.alpha<-piFIT$fit
  check.wt<-max(1/piFIT$fit[trt==ss & delta==1])
  ave.wt<-mean(trtnew/piFIT$fit)
  phiFIT <- try(glm(y~as.matrix(phiX),family=famY,weights=1/(h.wt*piFIT$fit),subset=(trt==ss & delta==1),control=ctrl,na.action=na.omit),TRUE)
  failed.phi <- inherits(phiFIT,"try-error")
  if(failed.phi) phiFIT <- suppressWarnings(glm(y~as.matrix(phiX),family=famY,weights=1/(h.wt*piFIT$fit),subset=(trt==ss & delta==1),control=ctrl,na.action=na.omit))
  Q0 <- predict(phiFIT,as.data.frame(phiX),type="response")
  mu.hat <- mean(Q0)
  list(wt.alpha=wt.alpha,Q0=Q0,mu.hat=mu.hat,check.wt=check.wt,ave.wt=ave.wt,failed.pi=failed.pi,failed.phi=failed.phi,h.wt=h.wt)
}


# PLEASE takes the outcome (y), treatment indicator (trt), 
# baseline covariates that form the missing data model (dX),
# the propensity score model (piX), and the outcome regression
# model (phiX)
PLEASE <- function(y,trt,dX,piX,phiX,ctrl=list(epsilon=1e-7,maxit=50)) {
  if( !all (trt %in% c(0,1) ) ) stop( "'TRT' should be binary." )
  YisBin <- all(y[!is.na(y)] %in% c(0,1))
  famY <- ifelse(YisBin,"binomial","gaussian")
  # Get the Joffe estimators and save Pr(A=1/0|W) and Pr(Y=1|A,W)
  drwls1 <- PLEASE.calcs(y,trt,dX,piX,phiX,ss=1,famY,ctrl)
  drwls0 <- PLEASE.calcs(y,trt,dX,piX,phiX,ss=0,famY,ctrl)
  # Create the new augmenting variables for the missingness model
  newd1 <- trt/(drwls1$h.wt)*(drwls1$Q0 - drwls1$mu.hat)
  newd0 <- (1-trt)/(drwls0$h.wt)*(drwls0$Q0 - drwls0$mu.hat)
  dXnew <- cbind(dX,newd1,newd0)
  # Create the new predictor for the propensity score model
  newX1 <- drwls1$Q0 - drwls1$mu.hat
  newX0 <- drwls0$Q0 - drwls0$mu.hat
  piXnew <- cbind(piX,newX1,newX0)
  # Get the updated estimate using piX and newX1
  # in the propensity score model
  trt1 <- PLEASE.calcs(y,trt,dXnew,piXnew,phiX,ss=1,famY,ctrl)
  trt0 <- PLEASE.calcs(y,trt,dXnew,piXnew,phiX,ss=0,famY,ctrl)
  mu.bar.1 <- trt1$mu.hat
  mu.bar.0 <- trt0$mu.hat
  max.wt.1 <- trt1$check.wt
  max.wt.0 <- trt0$check.wt
  ave.wt.1 <- trt1$ave.wt
  ave.wt.0 <- trt0$ave.wt
  failed.piDR.1 <- trt1$failed.pi
  failed.phiDR.1 <- trt1$failed.phi
  failed.piDR.0 <- trt0$failed.pi
  failed.phiDR.0 <- trt0$failed.phi
  list(mu.bar.1,mu.bar.0,drwls1$failed.pi,drwls1$failed.phi,drwls0$failed.pi,drwls0$failed.phi,max.wt.1,max.wt.0,ave.wt.1,ave.wt.0,failed.piDR.1,failed.phiDR.1,failed.piDR.0,failed.phiDR.0)
}

upper.fun <- function(x,data,xs) {
  apply(x,2,my.sim,data,xs)
}



upper.fun.ci <- function(x,data,xs,k) {
  apply(x,2,my.sim,data,xs,k)
}

my.boot <- function(d,index,xs) {
  d2 <- d[index,]
  y <- d2$y
  trt <- d2$trt
  piX <- as.matrix(d2[,xs])
  phiX <- as.matrix(d2[,xs])
  una <- unlist(unadjusted(y,trt))[1:2]
  ros <- unlist(Rosenblum(y,trt,piX,phiX,dX=piX)[1:14])[1:2]
  rot <- unlist(Rotnitzky(y,trt,piX,phiX)[1:14])[1:2]
  out <- c(una[1]-una[2],ros[1]-ros[2],rot[1]-rot[2],una[1],una[2],ros[1],ros[2],rot[1],rot[2])
  out
}

my.boot_full <- function(f, index2, xs, nboot){
	f2 <- f[index2,]
	boot.out <- boot(f2, my.boot, R=nboot, xs=xs)
	c((sd(boot.out$t[,1]) - sd(boot.out$t[,2]))/sd(boot.out$t[,1]), (sd(boot.out$t[,1]) - sd(boot.out$t[,3]))/sd(boot.out$t[,1]))
}

#####################################
#####################################


run_analysis_boot <- function(data, xs, nboot=100){

	## Load the boot library
	library(boot)

	## Define y, trt, piX, phiX
	y <- data$y
	trt <- data$trt
	#xs <- c("age","node","grade","size")
	# NOTE:  when you add your subtype variables to xs you will need to create
	# dummy variables since that variable is a factor
	piX <- as.matrix(data[,xs])
	phiX <- as.matrix(data[,xs])

	## Call the functions and bootstrap
	# Call unadjusted function
	una <- unadjusted(y,trt)
	# Call the Rosenblum estimator
	adj <- Rosenblum(y,trt,piX,phiX,dX=piX)
	# Call the Rotnitzky estimoator
	adj2 <- Rotnitzky(y,trt,piX,phiX)
	# Set the seed
	set.seed(321)
	# Call the bootstrap
	#boot.out <- boot(data,my.boot,R=nboot,xs=xs)
	# Note: this is going to take FOREVER
	boot.out_full <- boot(data, my.boot_full, R=nboot, xs=xs, nboot=nboot)
	# Call on the boot.ci to calculate confidence intervals
	#b1 <- boot.ci(boot.out,index=1,type="bca")
	#b2 <- boot.ci(boot.out,index=2,type="bca")
	#b3 <- boot.ci(boot.out,index=3,type="bca")
	b1 <- boot.ci(boot.out_full, index=1, type="bca")
	b2 <- boot.ci(boot.out_full, index=2, type="bca")

	#list(una, adj, b1, b2)

	#list("una_est"=boot.out$t0[1], "ros_est"=boot.out$t0[2], "rot_est"=boot.out$t0[3],
	 #    "una_est_sd"=sd(boot.out$t[,1]), "ros_est_sd"=sd(boot.out$t[,2]), "rot_est_sd"=sd(boot.out$t[,3]),
	  #   "una_ci"=b1$bca[4:5], "ros_ci"=b2$bca[4:5], "rot_ci"=b3$bca[4:5]) #"fit"=adj[[15]], "fit2"=adj2[[15]])

	#boot.out
	list(una, adj, adj2, boot.out_full$t0[1], boot.out_full$t0[2], b1$bca[4:5], b2$bca[4:5])
}

# Jackknife instead of bootstrapping
run_analysis_jk <- function(data, xs){
	y <- data$y
        trt <- data$trt
        #xs <- c("age","node","grade","size")
        # NOTE:  when you add your subtype variables to xs you will need to create
        # dummy variables since that variable is a factor
        piX <- as.matrix(data[,xs])
        phiX <- as.matrix(data[,xs])

	# We will do leave-one-out
	nsample <- length(y)

	unavec <- adjvec <- adj2vec <- vector("numeric", nsample)

	for(i in 1:nsample){
		# Call unadjusted function
        	una <- unadjusted(y[-i],trt[-i])
	        # Call the Rosenblum estimator
	        adj <- Rosenblum(y[-i],trt[-i],as.matrix(piX[-i,]),as.matrix(phiX[-i,]),dX=as.matrix(piX[-i,]))
	        # Call the Rotnitzky estimoator
	        adj2 <- Rotnitzky(y[-i],trt[-i],as.matrix(piX[-i,]),as.matrix(phiX[-i,]))
		unavec[i] <- una[[1]] - una[[2]]
		adjvec[i] <- adj[[1]] - adj[[2]]
		adj2vec[i] <- adj2[[1]] - adj2[[2]]

	}

	juna <- mean(unavec)
	jadj <- mean(adjvec)
	jadj2 <- mean(adj2vec)
	vuna <- ((nsample - 1)/nsample)*sum((unavec - juna)^2)
	vadj <- ((nsample - 1)/nsample)*sum((adjvec - jadj)^2)
	vadj2 <- ((nsample - 1)/nsample)*sum((adj2vec - jadj2)^2)


	list(juna, vuna, jadj, vadj, jadj2, vadj2)

}

# No bootstrap or jackknife
run_analysis <- function(data, xs){
	options(warn = -1)
        y <- data$y
        trt <- data$trt
        #xs <- c("age","node","grade","size")
        # NOTE:  when you add your subtype variables to xs you will need to create
        # dummy variables since that variable is a factor
        piX <- as.matrix(data[,xs])
        phiX <- as.matrix(data[,xs])

	una <- unadjusted(y,trt)
        # Call the Colantuoni estimator
        adj2 <- tryCatch(PLEASE(y,trt,as.matrix(piX),as.matrix(phiX),dX=as.matrix(piX)), error=function(e) e)
        # Call the Rotnitzky estimoator
        adj <- tryCatch(Rotnitzky(y,trt,as.matrix(piX),as.matrix(phiX)), error=function(e) e)

	if(inherits(adj2, "error")){
		adj2 <- una
	}

	if(inherits(adj, "error")){
		adj <- una
	}

        # Output order is unadjusted, rotnitzky, colantuoni
	list(una[[1]]-una[[2]], adj[[1]]-adj[[2]], adj2[[1]]-adj2[[2]])
}

