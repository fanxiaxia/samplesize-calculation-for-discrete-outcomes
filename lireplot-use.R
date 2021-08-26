library(MASS)
library("lme4")
library(dummies)
library(beepr)


#####parameter combos
Ntime <- 6
type <- "Binomial"

#target Ncluster
#Ncluster <- 20
Totalindv <- 4*6
Nindv <- Totalindv




Npcp <- 3500/100
#185
#random effect
clinic.int <- 0.011
tx <- 0.0015
pcp.int <- 0.0054
#fixed effect
beta0 <- -1.43
beta_x1 <- -0.055
beta_tm <- -0.124


#help functions:
logit <- function(x){log(x/(1 - x))}
expit <- function(x){exp(x)/(1 + exp(x))}

#data.X <- design.fixed.fun(Nseq,Ntime,Ncluster,Nindv)[[1]]
#Time.vname <- design.fixed.fun(Nseq,Ntime,Ncluster,Nindv)[[2]]

#for different types of link functions, the corresponding W has different forms. 
# <<- is global variable 

W.parameter <- function(type,sigma=NA)
{
  if(type=="Binomial")
  { 
    phi <<- 1
    a <<- 1
    gprime <<- function(x){1/(x*(1-x))}
    v <<- function(x){x*(1-x)}
    h <<- expit
    g <<- logit
  } 
  if(type == "Normal")
  { 
    phi <<- sigma^2
    a <<- 1
    gprime <<- function(x){rep(1, times=length(x))}
    v <<- function(x){rep(1, times=length(x))}
    h <<- function(x){x}
    g <<- h
  } 
  if(type == "Poisson")
  {
    phi <<- 1
    #take 1 for now
    a <<- 1 
    gprime <<- function(x){1/x}
    v <<- function(x){x}
    h <<- function(x){exp(x)}
    g <<- function(x){log(x)}
  } 
}   

#When doing all clusters together, the covariance matrix is big and the computation is slow.
#We look at each cluster (\Sum X_iV_iX_i)^(-1)
Nseq <- Ntime-1
design.mtrx <- matrix(0, nrow = Nseq, ncol = Ntime)
design.mtrx[upper.tri(design.mtrx)] <- 1
time.mt <- matrix(rep(1:Ntime,each=Nseq),ncol=Ntime)

beta0.vector <- as.matrix(c(beta0,beta_tm,0))
beta1.vector<- as.matrix(c(beta0,beta_tm,beta_x1))



Nullhyp <- function(isnull,beta0=NA,beta1=NA)
{
  if(isnull)
  {
    beta.vector <<- beta0
    return("Under Null hypothesis")
  }
  if(!isnull)
  {
    beta.vector <<- beta1
    return("Under Alternative hypothesis")
  }
}

#Nindv <- 6
bycluster.proposed.fun <- function(isnull,Ncluster)
{
  dimD <- 2+Npcp
  #For each cluster, we have a random intercept, and a random treatment effect, and a random pcp effect. 
  #The covariance matrix of the random effects is
  Dc <- clinic.int^2
  Dtx <- tx^2
  
  Dpcp <- diag(rep(pcp.int^2,Npcp))
  D <- matrix(0,dimD,dimD)
  D[1,1] <- Dc
  D[2,2] <- Dtx
  D[3:dimD,3:dimD] <- Dpcp
  
  D.inv <- matrix(0,dimD,dimD)
  D.inv[1,1] <- 1/Dc
  D.inv[2,2] <- 1/Dtx
  Dpcp.inv <- diag(rep(pcp.int^(-2),Npcp))
  D.inv[3:dimD,3:dimD] <- Dpcp.inv
  
  Buse <- matrix(0,ncol = Nindv*Npcp,nrow=Npcp)
  for(i in 1:(Npcp))
  {
    Buse[i,((i-1)*Nindv+1):(i*Nindv)] <- 1
  }
  Buse <- t(Buse)
  
  
  #The fixed effect design matrix is (1 (for beta_0), Time1, Time2, Time3, Treatment)
  #The random effect design matrix is (1 (for usi), Treatment)
  
  #set the parameters
  Null.info <- Nullhyp(isnull,beta0.vector,beta1.vector)
  W.parameter(type,sigma)
  
  #Calculate X_iV_iX_i for cluster i 
  Info <- 0
  for(sequence in 1:Nseq)
  {
    treatment.temp <- design.mtrx[sequence,]
    time.temp <- time.mt[sequence,]
    for (i in 1:Ncluster)
    {
      treatment.i <- rep(treatment.temp,Npcp*Nindv/Ntime)
      time.i <- rep(time.temp,Npcp*Nindv/Ntime)
      X.cluster.i <- as.matrix(cbind(1,time.i,treatment.i))
      Z.cluster.i <- as.matrix(cbind(1,treatment.i,Buse))
      
      
      eta.approx.i <- X.cluster.i  %*% beta.vector
      mu.approx.i <- h(eta.approx.i)
      #for normal case it is actually not a function of mu.approx.i
      W.cluster.i <- diag(as.vector(phi*a*v(mu.approx.i)*(gprime(mu.approx.i))^2))
      V.cluster.i <- W.cluster.i + Z.cluster.i %*% D %*% t(Z.cluster.i)
      print(paste0("solve",i))
      
      #woodbury  
      w.inv <- diag(as.vector((phi*a*v(mu.approx.i)*(gprime(mu.approx.i))^2)^(-1)))
      v.inv <-  w.inv - w.inv %*% Z.cluster.i %*% solve(D.inv +  t(Z.cluster.i)%*%w.inv%*%Z.cluster.i)%*%t(Z.cluster.i)%*%w.inv
      
      Info <- Info + t(X.cluster.i) %*% v.inv %*% X.cluster.i
    }
  }
  
  #This result has been checked by three different ways of coding!!!
  Var.treatment <- solve(Info)["treatment.i","treatment.i"]
  
  return(data.frame(Nindv,Var.treatment))
}
power <- function(var.alt,var.null,betap,alpha=0.05)
{
  p1 <- pnorm((betap-qnorm(1-alpha/2)*sqrt(var.null))/sqrt(var.alt))
  p2 <- pnorm((-betap-qnorm(1-alpha/2)*sqrt(var.null))/sqrt(var.alt))
  if (p1 > alpha & p2 > alpha)
  {
    return("both large")
  }
  if (p1 > alpha & p2 < alpha)
  {
    return(p1)
  }
  if (p1 < alpha & p2 > alpha)
  {
    return(p2)
  }
  if (p1 < alpha & p2 < alpha)
  {
    return("both small")
  }
}

#Nindv <- 72

#real:238886/3761
#use 3500 (closer to the real number 3761) then Nindv is 54.


#var.proposed <- bycluster.proposed.fun(isnull,beta0,beta1,type,sigma,eta,tau,rho,data.X,Time.vname)
#Now we generate the data and estimate the variance empirically. 


power.cal <- function(Ncluster)
{
  var.proposed <- bycluster.proposed.fun(F,Ncluster)
  var.proposed.null <- bycluster.proposed.fun(T,Ncluster)
  
  powervalue <- power(var.proposed$Var.treatment,var.proposed.null$Var.treatment,beta_x1)
  beep()
  print(powervalue)
}
power.value <- rep(NA,20)
taketime <- rep(NA,20)
for(i in 1:20)
{
  start_time <- Sys.time()
  power.value[i] <- power.cal(i*2)
  end_time <- Sys.time()
  
  taketime[i] <- difftime(end_time,start_time,units = "secs")
  print(i)
}
write.csv(taketime,"taketime6.csv")
beep(sound = 5)


#i <- 20
library(ggplot2)
#plot(1:35, power.value,main="LIRE: Power VS Number of Clinics",xlab = "Number of Clinics",ylab = "Power")
clinic <- (1:20)*5*2
power.clinic <- power.value
clinicplot <- data.frame(clinic,power.clinic)
clinicplot4 <- clinicplot
clinicplot.final <- rbind(clinicplot4,clinicplot5,clinicplot6)
clinicplot.final$nindv <- c(rep(4*35,20),rep(5*35,20),rep(6*35,20))
clinicplot.final$nindv <- as.factor(clinicplot.final$nindv)
write.csv(clinicplot.final,"clinicplot.csv")
p <- ggplot(clinicplot.final, aes(clinic, power.clinic,colour=nindv))
p + geom_smooth(method = "loess")+
  ggtitle("LIRE: Power VS Number of Clinics") +
  xlab("Total Number of Clinics") + ylab("Power")+
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    axis.title.x = element_text(color="black", size=14, face="bold"),
    axis.title.y = element_text(color="black", size=14, face="bold")
  )+
  scale_x_continuous(breaks=clinic)+geom_hline(yintercept = 0.8,linetype="dashed",color="black",size=0.5)+ylim(0,1)+
  geom_segment(x = 135,y=0,xend=135,yend=0.8,linetype="dashed",color="black",size=0.5)+
  geom_segment(x = 160,y=0,xend=160,yend=0.8,linetype="dashed",color="black",size=0.5)+
  geom_segment(x = 200,y=0,xend=200,yend=0.8,linetype="dashed",color="black",size=0.5)+theme(legend.position="bottom")+
   scale_colour_discrete(name  ="Number of individuals per clinic-period")
  

