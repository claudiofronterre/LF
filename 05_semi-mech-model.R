lf <- readRDS(file = "data/lf_clean.rds")


subsamp <- sample(names(table(lf$Country)), 2)

lf <- lf[lf$Country %in% subsamp, ]


indICT <- which(lf$Method_1=="Serological")
indMF <- which(lf$Method_1=="Parasitological")
lf <- lf[c(indICT,indMF), ]

n1 <- length(indICT)
n2 <- length(indMF)

# data$Blood_volume_cat <- as.character(data$Blood_volume)
# data$Blood_volume_cat[data$Blood_volume == 120 | data$Blood_volume == 100] <- ">= 100"
# data$Blood_volume_cat <- factor(data$Blood_volume_cat,levels=c("20","60",">= 100"))
coords <- lf[, c("Longitude", "Latitude")]

library(PrevMap)
coords <- jitterDupCoords(coords, max = 5)

U <- as.matrix(dist(coords))
U2 <- as.matrix(dist(coords[(n1 + 1):(n1 + n2), ]))
hist(U)

units.m <- lf$Examined
y <- lf$Positive
n <- length(y)

D1 <- cbind(rep(1, n))
D2 <- cbind(rep(1, n2))

p <- ncol(D1)
q <- ncol(D2)

load("data/estim_LF_mech.RData")
beta1 <- estim$par[1:p]
beta2 <- estim$par[(p + 1):(p + q)]

sigma2 <- exp(estim$par[p + q + 1])
phi <- exp(estim$par[p + q + 2])
nu2 <- exp(estim$par[p + q + 3])

Sigma <- sigma2 * exp(-U / phi)
diag(Sigma) <- diag(Sigma) + sigma2 * nu2
Sigma.inv <- solve(Sigma)

n <- n1 + n2
indICT <- 1:n1
indMF <- (n1 + 1):(n1 + n2)
indICT.S <- 1:n
indMF.S <- (n + 1):(n2 + n)

mu1 <- as.numeric(D1 %*% beta1)
mu2 <- as.numeric(D2 %*% beta2)
p.sens <- 0.97

r0 <- exp(mu2)
S <- mu1
integrand <- function(S) {
  m <- exp(S)
  pICT <- p.sens * (1 - exp(-m[indICT])) 
  pMF <- 1 - exp(-m[indMF] * (1 - exp(-r0)))
  
  diff.S <- S - mu1
  q.f_S <- t(diff.S) %*% Sigma.inv %*% (diff.S)
  
  llik <- sum(y[indICT] * log(pICT / (1 - pICT)) + units.m[indICT] * log(1 - pICT)) + 
    sum(y[indMF] * log(pMF / (1 - pMF)) + units.m[indMF] * log(1 - pMF))
  
  as.numeric(-0.5 * q.f_S + llik)
}

grad.integrand <- function(S) {
  m <- exp(S)
  pICT <- p.sens * (1 - exp(-m[indICT]))
  pMF <- 1 - exp(-m[indMF] * (1 - exp(-r0)))
  der.pICT <- p.sens * m[indICT] * exp(-m[indICT])  
  der.pMF.S <- (1 - exp(- r0)) * exp(-m[indMF] * (1 - exp(-r0))) * m[indMF]
  
  auxMF <- (y[indMF]/(pMF*(1-pMF))-units.m[indMF]/(1-pMF))
  
  diff.S <- S-mu1
  
  der.q.f_S <- -Sigma.inv%*%diff.S
  der.llik <- c((y[indICT]/(pICT*(1-pICT))-units.m[indICT]/(1-pICT))*der.pICT,
                auxMF*der.pMF.S)
  out <- der.q.f_S+der.llik
  as.numeric(out)
}

hessian.integrand <- function(S) {
  m <- exp(S)
  pICT <- p.sens*(1-exp(-m[indICT]))
  pMF <- 1-exp(-m[indMF]*(1-exp(-r0)))
  der.pICT <- p.sens*m[indICT]*exp(-m[indICT])  
  der2.pICT <- p.sens*exp(-m[indICT])*m[indICT]*(1-m[indICT])
  der.pMF.S <- (1-exp(-r0))*exp(-m[indMF]*(1-exp(-r0)))*m[indMF]
  der2.pMF.S <- der.pMF.S-((1-exp(-r0))^2)*exp(-m[indMF]*(1-exp(-r0)))*(m[indMF]^2)
  
  
  auxMF1 <- (y[indMF]/(pMF*(1-pMF))-units.m[indMF]/(1-pMF))
  auxMF2 <- (y[indMF]*((2*pMF-1)/((pMF-1)^2*pMF^2))+
               -units.m[indMF]/((1-pMF)^2))
  
  
  hess.llik1 <- c((y[indICT]*((2*pICT-1)/((pICT-1)^2*pICT^2))+
                     -units.m[indICT]/((1-pICT)^2))*(der.pICT)^2+
                    (y[indICT]/(pICT*(1-pICT))-units.m[indICT]/(1-pICT))*der2.pICT,
                  auxMF2*(der.pMF.S^2)+auxMF1*der2.pMF.S)
  
  res <- -Sigma.inv
  diag(res) <- diag(res)+hess.llik1
  
  return(res)
}

S.estim <- maxBFGS(
  integrand,
  grad.integrand,
  hessian.integrand,
  start=mu1,
  print.level = 1, iterlim = 1000)

S.estim$Sigma.tilde <- solve(-S.estim2$hessian)
Sigma.sroot <- t(chol(S.estim$Sigma.tilde))
A <- solve(Sigma.sroot)
library(Matrix)
Sigma.W.inv <- solve(A%*%Sigma%*%t(A))
mu.W <- as.numeric(A%*%(mu1-S.estim$estimate))

n.sim <- 10000
burnin <- 2000
thin <- 8
n.samples <- (n.sim-burnin)/thin
h <- 1.65/(n^(1/6))
c1.h <- 0.01
c2.h <- 1e-04

cond.dens.W <- function(W,S) {
  diff.W <- W-mu.W
  
  m <- exp(S)
  pICT <- p.sens*(1-exp(-m[indICT]))
  pMF <- 1-exp(-m[indMF]*(1-exp(-r0)))
  
  
  llik <- sum(y[indICT]*log(pICT/(1-pICT))+units.m[indICT]*log(1-pICT))+
    sum(y[indMF]*log(pMF/(1-pMF))+units.m[indMF]*log(1-pMF))
  
  as.numeric(-0.5*as.numeric(t(diff.W)%*%Sigma.W.inv%*%diff.W))+llik
}

lang.grad <- function(W,S) {
  diff.W <- W-mu.W
  
  m <- exp(S)
  pICT <- p.sens*(1-exp(-m[indICT]))
  pMF <- 1-exp(-m[indMF]*(1-exp(-r0)))
  der.pICT <- p.sens*m[indICT]*exp(-m[indICT])  
  der.pMF.S <- (1-exp(-r0))*exp(-m[indMF]*(1-exp(-r0)))*m[indMF]
  
  auxMF <- (y[indMF]/(pMF*(1-pMF))-units.m[indMF]/(1-pMF))
  
  diff.S <- S-mu1
  der.q.f_S <- -Sigma.inv%*%diff.S
  der.llik <- c((y[indICT]/(pICT*(1-pICT))-units.m[indICT]/(1-pICT))*der.pICT,
                auxMF*der.pMF.S)
  out <- as.numeric(-Sigma.W.inv%*%diff.W+t(Sigma.sroot)%*%der.llik)
}

S.estim$mode <- S.estim$estimate
W.curr <- rep(0,n)
S.curr <- as.numeric(Sigma.sroot%*%W.curr+S.estim$mode)
mean.curr <- as.numeric(W.curr + (h^2/2)*lang.grad(W.curr,S.curr))
lp.curr <- cond.dens.W(W.curr,S.curr)

acc <- 0
sim <- matrix(NA,nrow=n.samples,ncol=n)
h.vec <- rep(NA,n.sim)


for(i in 1:n.sim) {
  W.prop <- mean.curr+h*rnorm(n)
  S.prop <-  as.numeric(Sigma.sroot%*%W.prop+S.estim$mode)
  mean.prop <- as.numeric(W.prop + (h^2/2)*lang.grad(W.prop,S.prop))
  lp.prop <- cond.dens.W(W.prop,S.prop)
  
  dprop.curr <- -sum((W.prop-mean.curr)^2)/(2*(h^2))
  dprop.prop <- -sum((W.curr-mean.prop)^2)/(2*(h^2))
  
  log.prob <- lp.prop+dprop.prop-lp.curr-dprop.curr
  
  if(log(runif(1)) < log.prob) {
    acc <- acc+1
    W.curr <- W.prop
    S.curr <- S.prop
    lp.curr <- lp.prop
    mean.curr <- mean.prop
  }
  
  if( i > burnin & (i-burnin)%%thin==0) {
    sim[(i-burnin)/thin,] <- S.curr
  }
  
  h.vec[i] <- h <- max(0,h + c1.h*i^(-c2.h)*(acc/i-0.57))
  cat("Iteration",i,"out of",n.sim,"\r")
  flush.console()
}

if(TRUE) {
  acf.plot <- acf(sim[,1],plot=FALSE)
  plot(acf.plot$lag,acf.plot$acf,type="l",xlab="lag",ylab="autocorrelation",
       ylim=c(-0.1,1),main="Autocorrelogram of the simulated samples")
  for(i in 2:ncol(sim)) {
    acf.plot <- acf(sim[,i],plot=FALSE)
    lines(acf.plot$lag,acf.plot$acf)
  }
  abline(h=0,lty="dashed",col=2)
}

log.integrand <- function(S,val) {
  m <- exp(S)
  r0 <- exp(val$mu2)
  pICT <- p.sens*(1-exp(-m[indICT]))
  pMF <- 1-exp(-m[indMF]*(1-exp(-r0)))
  
  diff.S <- S-val$mu1
  q.f_S <- -0.5*(n*log(val$sigma2)+val$l.det.R+t(diff.S)%*%
                   val$R.inv%*%(diff.S)/val$sigma2)
  
  llik <- sum(y[indICT]*log(pICT/(1-pICT))+units.m[indICT]*log(1-pICT))+
    sum(y[indMF]*log(pMF/(1-pMF))+units.m[indMF]*log(1-pMF))
  
  as.numeric(q.f_S+llik)
}

compute.log.f <- function(par,l.det.R=NA,R.inv=NA) {
  beta1 <- par[1:p]
  beta2 <- par[(p+1):(p+q)]
  phi <- exp(par[p+q+2])
  nu2 <- exp(par[p+q+3])
  
  
  val <- list()
  val$mu1 <- as.numeric(D1%*%beta1)
  val$mu2 <- as.numeric(D2%*%beta2)
  val$sigma2 <- exp(par[p+q+1])
  if(is.na(l.det.R) & is.na(as.numeric(R.inv)[1])) {
    R <- exp(-U/phi)
    diag(R) <- diag(R) + nu2
    val$l.det.R <- determinant(R)$modulus
    val$R.inv <- solve(R)
  } else {
    val$l.det.R <- l.det.R
    val$R.inv <- R.inv
  }
  sapply(1:n.samples,function(i) log.integrand(sim[i,],val))
}

par0 <- c(beta1,beta2,log(c(sigma2,phi,nu2)))

rm(beta1,beta2,sigma2,phi,nu2,r0,mu1,mu2)


log.f.tilde <- compute.log.f(par0)

MC.log.lik <- function(par) {
  log(mean(exp(compute.log.f(par)-log.f.tilde)))
}


grad.MC.log.lik <- function(par) {
  beta1 <- par[1:p]; mu1 <- D1%*%beta1
  beta2 <- par[(p+1):(p+q)]; mu2 <- D2%*%beta2
  sigma2 <- exp(par[p+q+1])
  phi <- exp(par[p+q+2])
  nu2 <- exp(par[p+q+3])
  
  R <- exp(-U/phi)
  diag(R) <- diag(R)+nu2
  R.inv <- solve(R)
  l.det.R <- determinant(R)$modulus
  
  exp.fact <- exp(compute.log.f(par,l.det.R,R.inv)-log.f.tilde)
  L.m <- sum(exp.fact)
  exp.fact <- exp.fact/L.m
  
  
  R1.phi <- (U*exp(-U/phi))/phi
  m1.phi <- R.inv%*%R1.phi
  t1.phi <- -0.5*sum(diag(m1.phi))
  m2.phi <- m1.phi%*%R.inv
  
  t1.nu2 <- -0.5*sum(diag(R.inv))*nu2
  m2.nu2 <- R.inv%*%R.inv*nu2
  
  r0 <- exp(mu2)
  
  gradient.S <- function(S) {
    diff.S <- S-mu1
    m.MF <-exp(S[indMF])
    pMF <- 1-exp(-m.MF*(1-exp(-r0)))  
    auxMF <- (y[indMF]/(pMF*(1-pMF))-units.m[indMF]/(1-pMF))
    
    der.pMF.mu2 <- exp(-m.MF*(1-exp(-r0)))*(m.MF*exp(-r0))*r0
    
    q.f <- t(diff.S)%*%R.inv%*%diff.S
    grad.beta1 <-  t(D1)%*%R.inv%*%(diff.S)/sigma2
    grad.beta2 <-  as.numeric(t(D2)%*%(auxMF*der.pMF.mu2))
    
    grad.log.sigma2 <- (-n/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2
    
    grad.log.phi <- (t1.phi+0.5*as.numeric(t(diff.S)%*%m2.phi%*%(diff.S))/sigma2)
    
    grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(diff.S)%*%m2.nu2%*%(diff.S))/sigma2)
    
    out <- c(grad.beta1,grad.beta2,grad.log.sigma2,grad.log.phi,grad.log.nu2)
    
    return(out)
  }
  out <- rep(0,length(par))
  for(i in 1:n.samples) {
    out <- out + exp.fact[i]*gradient.S(sim[i,])
  }
  out
}

hess.MC.log.lik <- function(par) {
  beta1 <- par[1:p]; mu1 <- D1%*%beta1
  beta2 <- par[(p+1):(p+q)]; mu2 <- D2%*%beta2
  sigma2 <- exp(par[p+q+1])
  phi <- exp(par[p+q+2])
  nu2 <- exp(par[p+q+3])
  
  R <- exp(-U/phi)
  diag(R) <- diag(R)+nu2
  R.inv <- solve(R)
  l.det.R <- determinant(R)$modulus
  
  exp.fact <- exp(compute.log.f(par,l.det.R,R.inv)-log.f.tilde)
  L.m <- sum(exp.fact)
  exp.fact <- exp.fact/L.m
  
  
  R1.phi <- (U*exp(-U/phi))/phi
  m1.phi <- R.inv%*%R1.phi
  t1.phi <- -0.5*sum(diag(m1.phi))
  m2.phi <- m1.phi%*%R.inv
  
  t1.nu2 <- -0.5*sum(diag(R.inv))*nu2
  m1.nu2 <- R.inv*nu2
  m2.nu2 <- R.inv%*%R.inv*nu2
  
  R2.phi <- R1.phi+(U*(U-2*phi)*exp(-U/phi))/phi^2
  t2.phi <- -0.5*(sum(R.inv*R2.phi)-sum(m1.phi*t(m1.phi)))
  n2.phi <- R.inv%*%(2*R1.phi%*%m1.phi-R2.phi)%*%R.inv
  
  t2.nu2 <- -0.5*(sum(diag(R.inv)*nu2)-sum(diag(m2.nu2)*nu2))
  nu2.aux <- 2*m1.nu2*nu2
  diag(nu2.aux) <- diag(nu2.aux)-nu2
  n2.nu2 <- R.inv%*%nu2.aux%*%R.inv
  
  t2.phi.nu2 <- 0.5*sum(m1.phi*t(R.inv))*nu2
  n2.phi.nu2 <- R.inv%*%(nu2*m1.phi+nu2*R1.phi%*%R.inv)%*%R.inv
  
  ind.beta1 <- 1:p
  ind.beta2 <- (p+1):(p+q)
  ind.sigma2 <- p+q+1
  ind.phi <- p+q+2
  ind.nu2 <- p+q+3
  
  H <- matrix(0,nrow=length(par),ncol=length(par))
  H[ind.beta1,ind.beta1] <- -t(D1)%*%R.inv%*%D1/sigma2
  
  r0 <- exp(mu2)
  hessian.S <- function(S,ef) {
    diff.S <- S-mu1
    m.MF <-exp(S[indMF])
    pMF <- 1-exp(-m.MF*(1-exp(-r0)))  
    auxMF <- (y[indMF]/(pMF*(1-pMF))-units.m[indMF]/(1-pMF))
    
    der.pMF.mu2 <- exp(-m.MF*(1-exp(-r0)))*(m.MF*exp(-r0))*r0
    
    q.f <- t(diff.S)%*%R.inv%*%diff.S
    grad.beta1 <-  t(D1)%*%R.inv%*%(diff.S)/sigma2
    grad.beta2 <-  as.numeric(t(D2)%*%(auxMF*der.pMF.mu2))
    
    grad.log.sigma2 <- (-n/(2*sigma2)+0.5*q.f/(sigma2^2))*sigma2
    
    grad.log.phi <- (t1.phi+0.5*as.numeric(t(diff.S)%*%m2.phi%*%(diff.S))/sigma2)
    
    grad.log.nu2 <-  (t1.nu2+0.5*as.numeric(t(diff.S)%*%m2.nu2%*%(diff.S))/sigma2)
    
    g <- c(grad.beta1,grad.beta2,grad.log.sigma2,grad.log.phi,grad.log.nu2)
    
    # exp(-m.MF*(1-exp(-r0)))*(m.MF*exp(-r0))*r0
    der2.pMF.mu2 <- der.pMF.mu2+
      (-(exp(-m.MF*(1-exp(-r0))))*((m.MF*exp(-r0))^2)+
         -exp(-m.MF*(1-exp(-r0)))*(m.MF*exp(-r0)))*(r0^2)
    
    auxMF1 <- (y[indMF]/(pMF*(1-pMF))-units.m[indMF]/(1-pMF))
    auxMF2 <- (y[indMF]*((2*pMF-1)/((pMF-1)^2*pMF^2))+
                 -units.m[indMF]/((1-pMF)^2))
    
    der2.pMF.mu2.tot <- auxMF2*(der.pMF.mu2^2)+auxMF1*der2.pMF.mu2
    
    H[ind.beta1,ind.sigma2] <-
      H[ind.sigma2,ind.beta1] <- -t(D1)%*%R.inv%*%(diff.S)/sigma2
    H[ind.beta2,ind.beta2] <- t(D2)%*%(D2*der2.pMF.mu2.tot)
    
    H[ind.beta1,ind.phi] <-
      H[ind.phi,ind.beta1] <- -as.numeric(t(D1)%*%m2.phi%*%(diff.S))/sigma2
    
    H[ind.beta1,ind.nu2] <-
      H[ind.nu2,ind.beta1] <- -as.numeric(t(D1)%*%m2.nu2%*%(diff.S))/sigma2
    
    H[ind.sigma2,ind.sigma2] <- (n/(2*sigma2^2)-q.f/(sigma2^3))*sigma2^2+
      grad.log.sigma2
    
    H[ind.sigma2,ind.phi] <-
      H[ind.phi,ind.sigma2] <- -(grad.log.phi-t1.phi)
    
    H[ind.sigma2,ind.nu2] <-
      H[ind.nu2,ind.sigma2] <- -(grad.log.nu2-t1.nu2)
    
    H[ind.phi,ind.phi] <- t2.phi-0.5*t(diff.S)%*%n2.phi%*%(diff.S)/sigma2
    
    H[ind.nu2,ind.nu2] <- t2.nu2-0.5*t(diff.S)%*%n2.nu2%*%(diff.S)/sigma2
    H[ind.phi,ind.nu2] <-
      H[ind.nu2,ind.phi] <- (t2.phi.nu2-0.5*t(diff.S)%*%n2.phi.nu2%*%(diff.S)/sigma2)
    
    
    
    out <- list()
    out$mat1<- ef*(g%*%t(g)+H)
    out$g <- g*ef
    out
  }
  
  a <- rep(0,length(par))
  A <- matrix(0,length(par),length(par))
  for(i in 1:n.samples) {
    out.i <- hessian.S(sim[i,],exp.fact[i])
    a <- a+out.i$g
    A <- A+out.i$mat1
  }
  (A-a%*%t(a))
}


estim <- nlminb(par0,
                function(x) -MC.log.lik(x),
                function(x) -grad.MC.log.lik(x),
                function(x) -hess.MC.log.lik(x),
                control=list(trace=1))
estim$gradient <- grad.MC.log.lik(estim$par)
estim$hessian <- hess.MC.log.lik(estim$par)
#save(estim,file="estim_LF_mech2.RData")
