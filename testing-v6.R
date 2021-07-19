library("swCRTdesign")

################################################functions
swGlmPwr_new <- function(design,distn,n,fixed.intercept,fixed.treatment.effect,fixed.time.effect,tau=0,eta=0,rho=0,gamma=0,alpha=0.05,retDATA=F)
{
  #help functions:
  logit <- function(x){log(x/(1 - x))}
  expit <- function(x){exp(x)/(1 + exp(x))}
  if (alpha <= 0 | alpha >= 1) {
    stop("Alpha must be strictly between 0 and 1.")
  }
  if (!all(n%%1 == 0)) {
    warning("n (either scalar, vector, or matrix) should consist only of integers.")
  }
  if (rho < -1 | rho > 1) {
    stop("The correlation between the random cluster effect and the random treatment effect must be a numeral between -1 and 1.")
  }
  if (tau < 0 | eta < 0 | eta < 0) {
    stop("Standard deviations of the random effects must be non-negative numerals.")
  }
  
  if(distn =="binomial")
  { 
    phi <- 1
    a <- 1
    gprime <- function(x){1/(x*(1-x))}
    v <- function(x){x*(1-x)}
    h <- expit
    g <- logit
    print("Outcome type: binomial")
  } 
  if(distn == "poisson")
  {
    phi <- 1
    #take 1 for now
    a <- 1 
    gprime <- function(x){1/x}
    v <- function(x){x}
    h <- function(x){exp(x)}
    g <- function(x){log(x)}
    print("Outcome type: poisson")
  } 
  if(distn != "binomial" && distn != "poisson")
  {
    stop("Valid outcome distributions are 'binomial' and 'poisson'")
  }
  #check inputs of fixed and random effect components  
  if(is.null(fixed.intercept) | is.null(fixed.time.effect) | is.null(fixed.treatment.effect))
  {
    stop("Parameters for fixed effects (intercept, treatment effect, and time effect) 
           must be specified.")
  }
  if(is.null(tau) | is.null(gamma) | is.null(eta))
  {
    warning("Standard deviations of all random effects are needed. If a random effect does not exist, its standard deviation should be set to 0. The default standard deviations are set to 0. ")
    if(is.null(tau)) {tau <- 0}
    if(is.null(gamma)) {gamma <- 0}
    if(is.null(eta)) {eta <- 0}
  }
  if(is.null(rho))
  {
    rho <- 0 
    if(eta*gamma != 0)
    {
      warning("The correlation between the random cluster effect and the random treatment effect (rho) is set to 0.")
    }
  }
  if(eta*gamma == 0 & (!is.null(rho)) & (rho!=0))
  {
    warning("The correlation between the random cluster effect and the random treatment effect (rho) is not used.")
  }
  
  power <- function(var.alt,var.null,betap,alpha)
  {
    p1 <- pnorm((betap-qnorm(1-alpha/2)*sqrt(var.null))/sqrt(var.alt))
    p2 <- pnorm((betap+qnorm(1-alpha/2)*sqrt(var.null))/sqrt(var.alt))
    return(p1 + 1 - p2)
  }
  
  #design related
  Nseq <- design$n.waves
  cldistninseq <- design$clusters
  design.mtrx <- design$swDsn.unique.clusters
  Ntime <- design$total.time
  time.mt <- matrix(rep(1:Ntime,each=Nseq),ncol=Ntime)
  
  #check if n is correctly specified  
  if (length(n) == 1) 
  {
    size.matrix <- matrix(n,nrow = sum(cldistninseq),ncol = Ntime)
  }
  else 
  {
    if ((is.vector(n) & length(n) > 1)) 
    {
      if (length(n) != design$n.clusters) 
      {
        stop("The number of clusters in 'design' (design$n.clusters) and 'n' (length(n)) do not match.")
      }
      size.matrix <- matrix(rep(n,Ntime),ncol = Ntime)
    }
    else if (is.matrix(n)) 
    {
      if ((nrow(n) != design$n.clusters) | (ncol(n) != design$total.time)) 
      {
        stop("The number of clusters and/or time steps in 'design' (design$n.clusters and design$total.time) and 'n' (number of rows and columns, respectively) do not match.")
      }
      size.matrix <- n
    }
  }
  
  ##fixed.effects
  #a vector of fixed effect under the alternative. c(fixed.intercept, fixed.time effect, fixed.treatment effect)
  if(length(fixed.time.effect)==1)
  {
    fixed.time.effect <- rep(fixed.time.effect,Ntime - 1)
  }
  else
  {
    if((is.vector(fixed.time.effect)) & length(fixed.time.effect) > 1)
    {
      if(length(fixed.time.effect) != Ntime - 1)
      {
        stop(paste0("Time effect must be specified for time 2 through time", Ntime, ", as it is coded as dummy variables with time 1 being the reference group."))
      }
    }
  }
  fixed.effects <- c(fixed.intercept, fixed.time.effect,fixed.treatment.effect)
  beta0.vector <- beta1.vector <- fixed.effects
  beta_x1 <- beta1.vector[length(beta0.vector)]
  beta0.vector[length(beta0.vector)] <- 0
  
  #approximation
  Info1 <- Info0 <- 0
  for(seq in 1:Nseq)
  {
    if(seq ==1)
    {
      cluster.index <- 1: cldistninseq[seq]
    }
    if(seq > 1)
    {
      cluster.index <- (sum(cldistninseq[1:(seq-1)])+1):sum(cldistninseq[1:(seq)])
    }
    
    clustersize <- matrix(size.matrix[cluster.index,],ncol=Ntime)
    treatment.temp <- as.matrix(design.mtrx[seq,])
    colnames(treatment.temp) <- "treatment.i"
    time.temp <- time.mt[seq,]
    time.temp.mtrx <- model.matrix(~factor(time.temp))
    time.temp.mtrx.z <- diag(1,nrow = Ntime,ncol = Ntime)
    
    X.seq.i  <- cbind(time.temp.mtrx, treatment.temp)
    mu.null.i <- h(X.seq.i  %*% beta0.vector)
    mu.alt.i <- h(X.seq.i  %*% beta1.vector)
    
    w.null.i <- (phi*a*v(mu.null.i)*(gprime(mu.null.i))^2)
    w.alt.i <- (phi*a*v(mu.alt.i)*(gprime(mu.alt.i))^2)
    
    varfun <- function(type,clustersize,treatment.temp,Ntime,time.temp.mtrx,tau,eta,gamma,rho)
    {
      if(type == "null")
      {
        w <- w.null.i 
      }
      if(type == "alt")
      {
        w <- w.alt.i
      }
      info.seq <- 0
      for(i in 1:nrow(clustersize))
      {
        tr <- treatment.temp
        LxTLx <- LxTtr <- trTtr <-  0
        for(j in 1:Ntime)
        {
          kx.j <-  as.matrix(time.temp.mtrx[j,])
          LxTLx <-  LxTLx + clustersize[i,j]*(kx.j%*%t(kx.j))/w[j]
          LxTtr <- LxTtr + clustersize[i,j]*(kx.j*tr[j])/w[j]
          trTtr <- trTtr + clustersize[i,j]*(tr[j])^2/w[j]
        }
        
        xTx <- rbind(cbind(LxTLx,LxTtr),cbind(t(LxTtr),trTtr))
        
        if(tau==0)
        {
          stop("The cluster random effect cannot be 0.")
        }
        #cluster + treatment 
        if(eta!=0 & gamma==0)
        {
          dimD <- 2
          D.inv <- matrix(0,dimD,dimD)
          D.inv <- (tau*tau*(1 - rho^2))^(-1)*matrix(c(eta^2,-tau*eta*rho,-tau*eta*rho,tau^2),2,2)
          
          LxTLz <- LzTtr <- LzTLz <- 0
          for(j in 1:Ntime)
          {
            kx.j <-  as.matrix(time.temp.mtrx[j,])
            kz.j <-  1
            LxTLz <-  LxTLz + clustersize[i,j]*(kx.j)/w[j]
            LzTtr <- LzTtr + clustersize[i,j]*(tr[j])/w[j]
            LzTLz <- LzTLz + clustersize[i,j]/w[j]
          }
          
          zTz <- rbind(cbind(LzTLz,LzTtr),cbind(t(LzTtr),trTtr))
          xTz <- rbind(cbind(LxTLz,LxTtr),cbind(t(LzTtr),trTtr))
          increase <- (xTx - xTz%*%solve(D.inv + zTz)%*%t(xTz))
        }
        #cluster + time + treatment
        if(eta!=0 & gamma!=0)
        {
          dimD <- 2+Ntime
          D.inv <- matrix(0,dimD,dimD)
          D.inv[1:2,1:2] <- (tau*tau*(1 - rho^2))^(-1)*matrix(c(eta^2,-tau*eta*rho,-tau*eta*rho,tau^2),2,2)
          Dtime.inv <- diag(rep(gamma^(-2),Ntime))
          D.inv[3:dimD,3:dimD] <- Dtime.inv
          
          LxTLz <- LzTtr <- LzTLz <- 0
          for(j in 1:Ntime)
          {
            kx.j <-  as.matrix(time.temp.mtrx[j,])
            kz.j <-  as.matrix(c(1,time.temp.mtrx.z[j,]))
            LxTLz <-  LxTLz + clustersize[i,j]*(kx.j %*% t(kz.j))/w[j]
            LzTtr <- LzTtr + clustersize[i,j]*(kz.j*tr[j])/w[j]
            LzTLz <- LzTLz + clustersize[i,j]*(kz.j %*% t(kz.j))/w[j]
          }
          
          zTz <- rbind(cbind(LzTLz,LzTtr),cbind(t(LzTtr),trTtr))
          xTz <- rbind(cbind(LxTLz,LxTtr),cbind(t(LzTtr),trTtr))
          increase <- (xTx - xTz%*%solve(D.inv + zTz)%*%t(xTz))
        }
        #cluster + time
        if(eta==0 & gamma!=0)
        {
          dimD <- 1 + Ntime
          D.inv <- matrix(0,dimD,dimD)
          D.inv[1,1] <- tau^(-2)
          Dtime.inv <- diag(rep(gamma^(-2),Ntime))
          D.inv[2:dimD,2:dimD] <- Dtime.inv
          
          zTz <- LxTLz <- LzTtr <- 0
          for(j in 1:Ntime)
          {
            kx.j <-  as.matrix(time.temp.mtrx[j,])
            kz.j <-  as.matrix(c(1,time.temp.mtrx.z[j,]))
            zTz <-  zTz + clustersize[i,j]*(kz.j %*% t(kz.j))/w[j]
            LxTLz <-  LxTLz + clustersize[i,j]*(kx.j %*% t(kz.j))/w[j]
            LzTtr <-  LzTtr + clustersize[i,j]*(kz.j * tr[j])/w[j]
          }
          xTz <- rbind(LxTLz,t(LzTtr))
          increase <- (xTx - xTz%*%solve(D.inv + zTz)%*%t(xTz))
        }
        #only random cluster effect
        if(eta==0 & gamma==0)
        {
          D.inv <- tau^(-2)
          
          zTz <- LxTLz <-LzTtr <-  0
          for(j in 1:Ntime)
          {
            kx.j <-  as.matrix(time.temp.mtrx[j,])
            kz.j <-  1
            zTz <-  zTz + clustersize[i,j]/w[j]
            LxTLz <-  LxTLz + clustersize[i,j]*(kx.j)/w[j]
            LzTtr <-  LzTtr + clustersize[i,j]*(tr[j])/w[j]
          }
          xTz <- rbind(LxTLz,t(LzTtr))
          increase <- (xTx - xTz%*%t(xTz)*(D.inv + zTz)^(-1))
        }
        info.seq <- info.seq + increase
      }
      info.seq
    }
    
    Info0 <- Info0 +   varfun("null",clustersize,treatment.temp,Ntime,time.temp.mtrx,tau,eta,gamma,rho)
    Info1 <- Info1 +   varfun("alt",clustersize,treatment.temp,Ntime,time.temp.mtrx,tau,eta,gamma,rho)
  }
  var.theta.alt <- solve(Info1)[Ntime+1,Ntime+1]
  var.theta.null <- solve(Info0)[Ntime+1,Ntime+1]
  
  pwrGLM <- power(var.theta.alt,var.theta.null,beta_x1,alpha)
  
  
  
  if(retDATA==T)
  {
    return(list(design, distn,n,fixed.intercept,fixed.treatment.effect,fixed.time.effect,tau,eta,rho,gamma,alpha,
                var.theta.null,var.theta.alt,pwrGLM))
  }
  
  if(retDATA==F)
  {
    return(pwrGLM)
  }
  
  
}

swGlmPwr <- function(design,distn,n,fixed.intercept,fixed.treatment.effect,fixed.time.effect,tau=0,eta=0,rho=0,gamma=0,alpha=0.05,retDATA=F)
{
  #help functions:
  logit <- function(x){log(x/(1 - x))}
  expit <- function(x){exp(x)/(1 + exp(x))}
  if (alpha <= 0 | alpha >= 1) {
    stop("Alpha must be strictly between 0 and 1.")
  }
  if (!all(n%%1 == 0)) {
    warning("n (either scalar, vector, or matrix) should consist only of integers.")
  }
  if (rho < -1 | rho > 1) {
    stop("The correlation between the random cluster effect and the random treatment effect must be a numeral between -1 and 1.")
  }
  if (tau < 0 | eta < 0 | eta < 0) {
    stop("Standard deviations of the random effects must be non-negative numerals.")
  }
  
  if(distn =="binomial")
  { 
    phi <- 1
    a <- 1
    gprime <- function(x){1/(x*(1-x))}
    v <- function(x){x*(1-x)}
    h <- expit
    g <- logit
    print("Outcome type: binomial")
  } 
  if(distn == "poisson")
  {
    phi <- 1
    #take 1 for now
    a <- 1 
    gprime <- function(x){1/x}
    v <- function(x){x}
    h <- function(x){exp(x)}
    g <- function(x){log(x)}
    print("Outcome type: poisson")
  } 
  if(distn != "binomial" && distn != "poisson")
  {
    stop("Valid outcome distributions are 'binomial' and 'poisson'")
  }
  #check inputs of fixed and random effect components  
  if(is.null(fixed.intercept) | is.null(fixed.time.effect) | is.null(fixed.treatment.effect))
  {
    stop("Parameters for fixed effects (intercept, treatment effect, and time effect) 
           must be specified.")
  }
  if(is.null(tau) | is.null(gamma) | is.null(eta))
  {
    warning("Standard deviations of all random effects are needed. If a random effect does not exist, its standard deviation should be set to 0. The default standard deviations are set to 0. ")
    if(is.null(tau)) {tau <- 0}
    if(is.null(gamma)) {gamma <- 0}
    if(is.null(eta)) {eta <- 0}
  }
  if(is.null(rho))
  {
    rho <- 0 
    if(eta*gamma != 0)
    {
      warning("The correlation between the random cluster effect and the random treatment effect (rho) is set to 0.")
    }
  }
  if(eta*gamma == 0 & (!is.null(rho)) & (rho!=0))
  {
    warning("The correlation between the random cluster effect and the random treatment effect (rho) is not used.")
  }
  
  power <- function(var.alt,var.null,betap,alpha)
  {
    p1 <- pnorm((betap-qnorm(1-alpha/2)*sqrt(var.null))/sqrt(var.alt))
    p2 <- pnorm((betap+qnorm(1-alpha/2)*sqrt(var.null))/sqrt(var.alt))
    return(p1 + 1 - p2)
  }
  
  #design related
  Nseq <- design$n.waves
  cldistninseq <- design$clusters
  design.mtrx <- design$swDsn.unique.clusters
  Ntime <- design$total.time
  time.mt <- matrix(rep(1:Ntime,each=Nseq),ncol=Ntime)
  
  #check if n is correctly specified  
  if (length(n) == 1) 
  {
    size.matrix <- matrix(n,nrow = sum(cldistninseq),ncol = Ntime)
  }
  else 
  {
    if ((is.vector(n) & length(n) > 1)) 
    {
      if (length(n) != design$n.clusters) 
      {
        stop("The number of clusters in 'design' (design$n.clusters) and 'n' (length(n)) do not match.")
      }
      size.matrix <- matrix(rep(n,Ntime),ncol = Ntime)
    }
    else if (is.matrix(n)) 
    {
      if ((nrow(n) != design$n.clusters) | (ncol(n) != design$total.time)) 
      {
        stop("The number of clusters and/or time steps in 'design' (design$n.clusters and design$total.time) and 'n' (number of rows and columns, respectively) do not match.")
      }
      size.matrix <- n
    }
  }
  
  ##fixed.effects
  #a vector of fixed effect under the alternative. c(fixed.intercept, fixed.time effect, fixed.treatment effect)
  if(length(fixed.time.effect)==1)
  {
    fixed.time.effect <- rep(fixed.time.effect,Ntime - 1)
  }
  else
  {
    if((is.vector(fixed.time.effect)) & length(fixed.time.effect) > 1)
    {
      if(length(fixed.time.effect) != Ntime - 1)
      {
        stop(paste0("Time effect must be specified for time 2 through time", Ntime, ", as it is coded as dummy variables with time 1 being the reference group."))
      }
    }
  }
  fixed.effects <- c(fixed.intercept, fixed.time.effect,fixed.treatment.effect)
  beta0.vector <- beta1.vector <- fixed.effects
  beta_x1 <- beta1.vector[length(beta0.vector)]
  beta0.vector[length(beta0.vector)] <- 0
  
  #approximation
  Info1 <- Info0 <- 0
  for(seq in 1:Nseq)
  {
    if(seq ==1)
    {
      cluster.index <- 1: cldistninseq[seq]
    }
    if(seq > 1)
    {
      cluster.index <- (sum(cldistninseq[1:(seq-1)])+1):sum(cldistninseq[1:(seq)])
    }
    
    clustersize <- matrix(size.matrix[cluster.index,],ncol=Ntime)
    treatment.temp <- design.mtrx[seq,]
    time.temp <- time.mt[seq,]
    for(i in 1:nrow(clustersize))
    {
      treatment.i <- rep(treatment.temp,clustersize[i,])
      time.i <- rep(time.temp,clustersize[i,])
      time.i.f <- factor(time.i)
      dummy.time <- model.matrix(~time.i.f)
      
      X.cluster.i <- as.matrix(cbind(dummy.time,treatment.i))
      
      
      if(tau==0)
      {
        stop("The cluster random effect cannot be 0.")
      }
      if(eta==0)
      {
        if(gamma==0)
        {
          D.inv <- tau^(-2)
          #D <- tau^2
          Buse <- 0
          Z.cluster.i <- matrix(1,ncol=1,nrow=sum(clustersize[i,]))
        }
        if(gamma != 0)
        {
          dimD <- 1 + Ntime
          D.inv <- matrix(0,dimD,dimD)
          D.inv[1,1] <- tau^(-2)
          Dtime.inv <- diag(rep(gamma^(-2),Ntime))
          D.inv[2:dimD,2:dimD] <- Dtime.inv
          
          # D <- matrix(0,dimD,dimD)
          # D[1,1] <- tau^2
          # D[2:dimD,2:dimD] <- diag(rep(gamma^(2),Ntime))
          
          Buse <- matrix(0,ncol = sum(clustersize[i,]),nrow = Ntime)
          for(j in 1:(Ntime))
          {
            if(j==1)
            {
              Buse[j,1:clustersize[i,j]] <- 1
            }
            if(j>1)
            {
              Buse[j,(sum(clustersize[i,1:(j-1)])+1):sum(clustersize[i,1:j])] <- 1
            }
          }
          Buse <- t(Buse)
          Z.cluster.i <- as.matrix(cbind(1,Buse))
        }
      }
      if(eta*gamma != 0)
      {
        dimD <- 2+Ntime
        
        D.inv <- matrix(0,dimD,dimD)
        D.inv[1:2,1:2] <- (tau*tau*(1 - rho^2))^(-1)*matrix(c(eta^2,-tau*eta*rho,-tau*eta*rho,tau^2),2,2)
        Dtime.inv <- diag(rep(gamma^(-2),Ntime))
        D.inv[3:dimD,3:dimD] <- Dtime.inv
        
        # D <- matrix(0,dimD,dimD)
        # D[1:2,1:2] <- matrix(c(tau^2,tau*eta*rho,tau*eta*rho,eta^2),2,2)
        # D[3:dimD,3:dimD] <- diag(rep(gamma^(2),Ntime))
        
        
        Buse <- matrix(0,ncol = sum(clustersize[i,]),nrow = Ntime)
        for(j in 1:(Ntime))
        {
          if(j==1)
          {
            Buse[j,1:clustersize[i,j]] <- 1
          }
          if(j>1)
          {
            Buse[j,(sum(clustersize[i,1:(j-1)])+1):sum(clustersize[i,1:j])] <- 1
          }
        }
        Buse <- t(Buse)
        Z.cluster.i <- as.matrix(cbind(1,Buse,treatment.i))
      }
      if(eta!=0 & gamma==0)
      {
        dimD <- 2
        
        D.inv <- matrix(0,dimD,dimD)
        D.inv <- (tau*tau*(1 - rho^2))^(-1)*matrix(c(eta^2,-tau*eta*rho,-tau*eta*rho,tau^2),2,2)
        
        Z.cluster.i <- as.matrix(cbind(1,treatment.i))
      }
      
      null.eta.approx.i <- X.cluster.i  %*% beta0.vector
      null.mu.approx.i <- h(null.eta.approx.i)
      alt.eta.approx.i <- X.cluster.i  %*% beta1.vector
      alt.mu.approx.i <- h(alt.eta.approx.i)
      
      #Woodbury  
      w.inv.null <- diag(as.vector((phi*a*v(null.mu.approx.i)*(gprime(null.mu.approx.i))^2)^(-1)))
      v.inv.null <-  w.inv.null - w.inv.null %*% Z.cluster.i %*% solve(D.inv +  t(Z.cluster.i)%*%w.inv.null%*%Z.cluster.i)%*%t(Z.cluster.i)%*%w.inv.null
      Info0 <- Info0 + t(X.cluster.i) %*% v.inv.null %*% X.cluster.i
      
      w.inv.alt <- diag(as.vector((phi*a*v(alt.mu.approx.i)*(gprime(alt.mu.approx.i))^2)^(-1)))
      v.inv.alt <-  w.inv.alt - w.inv.alt %*% Z.cluster.i %*% solve(D.inv +  t(Z.cluster.i)%*%w.inv.alt%*%Z.cluster.i)%*%t(Z.cluster.i)%*%w.inv.alt
      Info1 <- Info1 + t(X.cluster.i) %*% v.inv.alt %*% X.cluster.i
      
    }
  }
  var.theta.alt <- solve(Info1)["treatment.i","treatment.i"]
  var.theta.null <- solve(Info0)["treatment.i","treatment.i"]
  
  pwrGLM <- power(var.theta.alt,var.theta.null,beta_x1,alpha)
  
  
  if(retDATA==T)
  {
    return(list(design, distn,n,fixed.intercept,fixed.treatment.effect,fixed.time.effect,tau,eta,rho,gamma,alpha,
                var.theta.null,var.theta.alt,pwrGLM))
  }
  
  if(retDATA==F)
  {
    return(pwrGLM)
  }
  
  
}
#source("eff_swGlmPwr.R")
###helper
expit <- function(x){exp(x)/(1+exp(x))}
logit <- function(x){log(x/(1-x))}


################################################tests
#change scale (swPwr is in linear scale)
mu0 <- .05/6
mu1 <- -.025/6 + .05/6
swPwr(design = swDsn(rep(1,23)),distn="binomial",n=50,mu0=mu0,
      mu1= mu1,tau=.01/6, eta=0, rho=0, gamma=0,  alpha=.05)

beta0 <- logit(mu0)
beta1 <- logit(mu1) - beta0
mubar=(mu0+mu1)/2

tauprime = (.01/6)/(mubar*(1-mubar)) # based on delta method
swGlmPwr(design = swDsn(rep(1,23)),distn="binomial",n=50,fixed.intercept=beta0,fixed.time.effect=0,
             fixed.treatment.effect= beta1,tau=tauprime, eta=0, rho=0, gamma=0,  alpha=.05)
swGlmPwr_new(design = swDsn(rep(1,23)),distn="binomial",n=50,fixed.intercept=beta0,fixed.time.effect=0,
             fixed.treatment.effect= beta1,tau=tauprime, eta=0, rho=0, gamma=0,  alpha=.05)


#test-case large clusters
size = c(35219,53535,63785,456132,128670,96673,51454,156667,127440,68615,56502,17719,75931,58655,52874,75936)

mu0=log(28.62/(2*100000))

sdid = .31

sdtime= .15

delta = log(.6)

sdelta = abs(0.4*delta)
#size = c(3,5,6,4,12,9,5,15,12,6,5,17,7,5,5,7)


swGlmPwr(design=swDsn(c(4,3,5,4)),distn="binomial",n=size,fixed.intercept=mu0,fixed.time.effect = 1,fixed.treatment.effect = delta,tau=sdid,
         eta=sdelta,rho=0,gamma=sdtime,alpha=.05)

swGlmPwr_new(design=swDsn(c(4,3,5,4)),distn="binomial",n=size,fixed.intercept=mu0,fixed.time.effect = 1,fixed.treatment.effect = delta,tau=sdid,
             eta=sdelta,rho=0,gamma=sdtime,alpha=.05)

#for nongaussian, make sure it works, with n being a scalar, a vector, and a matrix
test.fun.nongaussian <- function(distn = "binomial",design.test,fixed.intercept=1,fixed.time.effect=1,fixed.treatment.effect=1,tau=1,gamma=1,eta=1,rho=0.1)
{
  res <- NULL
  #same size
  n=100
  res1 <- swGlmPwr(design=design.test,distn = distn, n = n,fixed.intercept = fixed.intercept,fixed.treatment.effect=fixed.treatment.effect,fixed.time.effect = fixed.time.effect,tau= tau,eta=eta,
                rho= rho,gamma= gamma,alpha = 0.05)
  
  res1.new <- swGlmPwr_new(design=design.test,distn = distn, n = n,fixed.intercept = fixed.intercept,fixed.treatment.effect=fixed.treatment.effect,fixed.time.effect = fixed.time.effect,tau= tau,eta=eta,
                       rho= rho,gamma= gamma,alpha = 0.05)
  
  ##different size for different clusters
  n <- rep(10*(1:design.test$n.clusters))
  res2 <- swGlmPwr(design=design.test,distn = distn, n = n,fixed.intercept = fixed.intercept,fixed.treatment.effect=fixed.treatment.effect,fixed.time.effect = fixed.time.effect,tau= tau,eta=eta,
                rho= rho,gamma= gamma,alpha = 0.05)
  
  res2.new <- swGlmPwr_new(design=design.test,distn = distn, n = n,fixed.intercept = fixed.intercept,fixed.treatment.effect=fixed.treatment.effect,fixed.time.effect = fixed.time.effect,tau= tau,eta=eta,
                       rho= rho,gamma= gamma,alpha = 0.05)
  
  ##different size for different cluster/time.
  n <- matrix(rep(10*(1:design.test$total.time),design.test$n.clusters),ncol = design.test$total.time) + diag(1,ncol = design.test$total.time,nrow = design.test$n.clusters)
  res3 <- swGlmPwr(design=design.test,distn = distn, n = n,fixed.intercept = fixed.intercept,fixed.treatment.effect=fixed.treatment.effect,fixed.time.effect = fixed.time.effect,tau= tau,eta=eta,
                rho= rho,gamma= gamma,alpha = 0.05)
  
  res3.new <- swGlmPwr(design=design.test,distn = distn, n = n,fixed.intercept = fixed.intercept,fixed.treatment.effect=fixed.treatment.effect,fixed.time.effect = fixed.time.effect,tau= tau,eta=eta,
                       rho= rho,gamma= gamma,alpha = 0.05)
  
  res <- cbind(res1,res2,res3)
  rownames(res) <- c("swGlmPwr")
  res.new <- cbind(res1.new,res2.new,res3.new)
  rownames(res.new) <- c("swGlmPwr.new")
  
  return(list(res,res.new))
}


#Different random effect combinations with different designs
#only random cluster effect
gamma.temp = 0 
eta.temp = 0
rho.temp = 0
#Unequal clusters per sequence
design.test <- swDsn(c(2,3,4))
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Extra time points between crossing over
design.test <- swDsn(c(1,0,1,0,1,0))
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Extra time in the beginning and end
design.test <- swDsn(c(0,0,0,1,1,1),extra.time = 3)
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Full treatment effect is never achieved?
#swDsn(c(1,1,1),tx.effect.frac=c(0.1,0.2,0.3))

#Full treatment effect is never achieved
design.test <- swDsn(c(1,1,1),all.ctl.time0 = FALSE)
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Some clusters don't get any treatment time
design.test <- swDsn(c(1,1,1),tx.effect.frac = c(0,1),all.ctl.time0 = FALSE)
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#######
#random cluster + random treatment
gamma.temp = 0 
eta.temp = 0.1
rho.temp = 0 
#Unequal clusters per sequence
design.test <- swDsn(c(2,3,4))
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Extra time points between crossing over
design.test <- swDsn(c(1,0,1,0,1,0))
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Extra time in the beginning and end
design.test <- swDsn(c(0,0,0,1,1,1),extra.time = 3)
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Full treatment effect is never achieved?
#swDsn(c(1,1,1),tx.effect.frac=c(0.1,0.2,0.3))

#Full treatment effect is never achieved
design.test <- swDsn(c(1,1,1),all.ctl.time0 = FALSE)
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Some clusters don't get any treatment time
design.test <- swDsn(c(1,1,1),tx.effect.frac = c(0,1),all.ctl.time0 = FALSE)
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)


#######
#random cluster effect and random time effect
gamma.temp = 0.1 
eta.temp = 0
rho.temp = 0 
#Unequal clusters per sequence
design.test <- swDsn(c(2,3,4))
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Extra time points between crossing over
design.test <- swDsn(c(1,0,1,0,1,0))
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Extra time in the beginning and end
design.test <- swDsn(c(0,0,0,1,1,1),extra.time = 3)
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Full treatment effect is never achieved?
#swDsn(c(1,1,1),tx.effect.frac=c(0.1,0.2,0.3))

#Full treatment effect is never achieved
design.test <- swDsn(c(1,1,1),all.ctl.time0 = FALSE)
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Some clusters don't get any treatment time
design.test <- swDsn(c(1,1,1),tx.effect.frac = c(0,1),all.ctl.time0 = FALSE)
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#####
#random cluster effect + random treatment effect + random time effect
gamma.temp = 0.1 
eta.temp = 0.1
rho.temp = 0.1
#Unequal clusters per sequence
design.test <- swDsn(c(2,3,4))
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Extra time points between crossing over
design.test <- swDsn(c(1,0,1,0,1,0))
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Extra time in the beginning and end
design.test <- swDsn(c(0,0,0,1,1,1),extra.time = 3)
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Full treatment effect is never achieved?
#swDsn(c(1,1,1),tx.effect.frac=c(0.1,0.2,0.3))

#Full treatment effect is never achieved
design.test <- swDsn(c(1,1,1),all.ctl.time0 = FALSE)
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)

#Some clusters don't get any treatment time
design.test <- swDsn(c(1,1,1),tx.effect.frac = c(0,1),all.ctl.time0 = FALSE)
#test.fun.guassian(distn = "gaussian",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "binomial",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)
test.fun.nongaussian(distn = "poisson",design = design.test,gamma = gamma.temp, eta = eta.temp, rho = rho.temp)






