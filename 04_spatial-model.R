if (!require("pacman")) install.packages("pacman")
pkgs = c("spBayes", "geoR", "PrevMap") # package names
pacman::p_load(pkgs, character.only = T)

rmvn <- function(n, mu = 0, V = matrix(1)){
  p <- length(mu)
  if(any(is.na(match(dim(V), p)))){stop("Dimension problem!")}
  D <- chol(V)
  t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}

set.seed(1)

##generate some data
n <- 100 ##number of locations
q <- 3 ##number of outcomes at each location
nltr <- q*(q+1)/2 ##number of triangular elements in the cross-covariance matrix

coords <- cbind(runif(n, 0, 1), runif(n, 0, 1))

##parameters for generating a multivariate spatial GP covariance matrix
theta <- rep(3 / 0.5, q) ##spatial decay

A <- matrix(0, q, q)
A[lower.tri(A, TRUE)] <- c(1, 1, -1, 1, 0.5, 0.25)
K <- A %*% t(A)
K ##spatial cross-covariance
cov2cor(K) ##spatial cross-correlation
Psi <- diag(0.01,q)
C <- mkSpCov(coords, K, Psi, theta, cov.model = "exponential")

w <- rmvn(1, rep(0, nrow(C)), C) ##spatial random effects

w.a <- w[seq(1, length(w), q)]
w.b <- w[seq(2, length(w), q)]
w.c <- w[seq(3, length(w), q)]

##covariate portion of the mean
x.a <- cbind(1, rnorm(n))
x.b <- cbind(1, rnorm(n))
x.c <- cbind(1, rnorm(n))
x <- mkMvX(list(x.a, x.b, x.c))

B.1 <- c(1, -1)
B.2 <- c(-1, 1)
B.3 <- c(-1, -1)
B <- c(B.1, B.2, B.3)

Psi <- c(0.1, 0.1, 0.1) ##non-spatial residual variance, i.e., nugget

beta.0 <- c(1, 5, 10)

y <- rmvn(1, rep(beta.0,n), C)
y <- rnorm(n*q, x %*% B + w, rep(sqrt(Psi), n))

y.a <- y[seq(1, length(y), q)]
y.b <- y[seq(2, length(y), q)]
y.c <- y[seq(3, length(y), q)]

##subsample to make spatially misaligned data
sub.1 <- 1:50
y.1 <- y.a[sub.1]
w.1 <- w.a[sub.1]
coords.1 <- coords[sub.1, ]
x.1 <- x.a[sub.1, ]

sub.2 <- 25:75
y.2 <- y.b[sub.2]
w.2 <- w.b[sub.2]
coords.2 <- coords[sub.2, ]
x.2 <- x.b[sub.2, ]

sub.3 <- 50:100
y.3 <- y.c[sub.3]
w.3 <- w.c[sub.3]
coords.3 <- coords[sub.3, ]
x.3 <- x.c[sub.3, ]

##call spMisalignLM
q <- 3
A.starting <- diag(1, q)[lower.tri(diag(1, q), TRUE)]
n.samples <- 5000

starting <- list("phi" = rep(3/0.5, q), "A" = A.starting, "Psi" = rep(1, q))
tuning <- list("phi" = rep(0.5, q), "A" = rep(0.01, length(A.starting)), "Psi" = rep(0.1, q))
priors <- list("phi.Unif" = list(rep(3 / 0.75, q), rep(3 / 0.25, q)),
               "K.IW" = list(q + 1, diag(0.1, q)), "Psi.ig" = list(rep(2, q), rep(0.1, q)))

m.1 <- spMisalignLM(list(y.1~ 1, y.2~ 1, y.3~ 1), 
                    coords=list(coords.1, coords.2, coords.3),
                    starting = starting, tuning = tuning, priors = priors, 
                    n.samples = n.samples, cov.model = "exponential", n.report = 100)

burn.in <- floor(0.75*n.samples)

x11()
plot(m.1$p.theta.samples, density=FALSE)

##recover regression coefficients and random effects
m.1 <- spRecover(m.1, start = burn.in)

round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.1$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)

##predict for all locations, i.e., observed and not observed
pred.covars <- list(as.matrix(rep(1, n)), as.matrix(rep(1, n)),
                    as.matrix(rep(1, n)))

out <- spPredict(m.1, start=burn.in, thin=10, pred.covars=pred.covars, 
                 pred.coords=list(coords, coords, coords))

##summary and check
quants <- function(x){quantile(x, prob=c(0.5,0.025,0.975))}

y.hat <- apply(out$p.y.predictive.samples, 1, quants)

##unstack and plot
y.a.hat <- y.hat[,1:n]
y.b.hat <- y.hat[,(n+1):(2*n)]
y.c.hat <- y.hat[,(2*n+1):(3*n)]

par(mfrow=c(1,3))
plot(y.a, y.a.hat[1,], xlab="Observed y.a", ylab="Fitted & predicted y.a",
     xlim=range(y), ylim=range(y.hat), main="")
arrows(y.a[-sub.1], y.a.hat[1,-sub.1], y.a[-sub.1], y.a.hat[2,-sub.1], length=0.02, angle=90)
arrows(y.a[-sub.1], y.a.hat[1,-sub.1], y.a[-sub.1], y.a.hat[3,-sub.1], length=0.02, angle=90)
lines(range(y.a), range(y.a))

plot(y.b, y.b.hat[1,], xlab="Observed y.b", ylab="Fitted & predicted y.b",
     xlim=range(y), ylim=range(y.hat), main="")
arrows(y.b[-sub.2], y.b.hat[1,-sub.2], y.b[-sub.2], y.b.hat[2,-sub.2], length=0.02, angle=90)
arrows(y.b[-sub.2], y.b.hat[1,-sub.2], y.b[-sub.2], y.b.hat[3,-sub.2], length=0.02, angle=90)
lines(range(y.b), range(y.b))

plot(y.c, y.c.hat[1,], xlab="Observed y.c", ylab="Fitted & predicted y.c",
     xlim=range(y), ylim=range(y.hat), main="")
arrows(y.c[-sub.3], y.c.hat[1,-sub.3], y.c[-sub.3], y.c.hat[2,-sub.3], length=0.02, angle=90)
arrows(y.c[-sub.3], y.c.hat[1,-sub.3], y.c[-sub.3], y.c.hat[3,-sub.3], length=0.02, angle=90)
lines(range(y.c), range(y.c))

# Apply to the dataset ------------------------------------------------

lf <- readRDS(file = "data/lf_clean.rds")
ict <- which(lf$Method_1 == "Serological")
mf <- which(lf$Method_1 == "Parasitological")

logit <- log((lf$Positive + 0.5)/(lf$Examined - lf$Positive + 0.5))
lf$logit <- logit
hist(lf$logit[ict])
hist(lf$logit[mf])

coords <- lf[, c("Longitude", "Latitude")]

vari_ict <- variog(coords = coords[ict, ], data = logit[ict], max.dist = max(dist(coords[ict, ])) * 2/3)
plot(vari_ict)
vari.fit_ict <- variofit(vari_ict, ini.cov.pars = c(2, 500), cov.model = "exponential", fix.nugget = FALSE, nugget = 0)
vari.fit_ict

vari_mf <- variog(coords = coords[mf, ], data = logit[mf], max.dist = max(dist(coords[mf, ])) * 2/3)
plot(vari_mf)
vari.fit_mf <- variofit(vari_mf, ini.cov.pars = c(2, 500), cov.model = "exponential", fix.nugget = FALSE, nugget = 0)
vari.fit_mf

par(mfrow = c(2, 1))
plot(vari_ict, main = "Fitted variogram ICT")
lines(vari.fit_ict)

plot(vari_mf, main = "Fitted variogram MF")
lines(vari.fit_mf)

coords_jt <- jitterDupCoords(coords, max = 0.1)
lf$Longitudej <- coords_jt[, 1]
lf$Latitudej <- coords_jt[, 2]

f <- as.formula(paste("logit", 1, sep = "~"))
fit.MLE_ict <- linear.model.MLE(formula = f, coords = ~ Longitudej + Latitudej, 
                                data = lf[ict, ], start.cov.pars = c(100, 0.5), kappa = 0.5, 
                                messages = TRUE, method = "nlminb")

fit.MLE_mf <- linear.model.MLE(formula = f, coords = ~ Longitudej + Latitudej, 
                               data = lf[mf, ], start.cov.pars = c(100, 0.5), kappa = 0.5, 
                               messages = TRUE, method = "nlminb")

f2 <- as.formula(paste("logit", "Method_1", sep = "~"))
fit.MLE <- linear.model.MLE(formula = f2, coords = ~ Longitudej + Latitudej, 
                            data = lf, start.cov.pars = c(100, 0.5), kappa = 0.5, 
                            messages = TRUE, method = "nlminb")


summary(fit.MLE_ict, log.cov.pars = FALSE)
summary(fit.MLE_mf, log.cov.pars = FALSE)
summary(fit.MLE, log.cov.pars = FALSE)

q <- 2
A.starting <- diag(2, q)[lower.tri(diag(1, q), TRUE)]
#A.starting <- c(3, 0.6, 1.8)
n.samples <- 5000

starting <- list("phi" = rep(1/0.01, q), "A" = A.starting, "Psi" = rep(1, q))
tuning <- list("phi"=rep(0.001,q), "A"=rep(0.001,length(A.starting)), "Psi"=rep(0.001,q))
priors <- list("phi.Unif" = list(rep(1/0.02, q), rep(1/0.0067, q)),
               "K.IW"=list(q+1, diag(c(15,1.5))), "Psi.ig"=list(rep(2,q), c(5,1)))

system.time(m.1 <- spMisalignLM(list(logit[ict] ~ 1, logit[mf] ~ 1), 
                    coords = list(coords_jt[ict, ], coords_jt[mf, ]),
                    starting = starting, tuning = tuning, priors = priors, 
                    n.samples = n.samples, cov.model = "exponential", n.report = 100))
save(m.1, file = "output/sp_misalign_model.rds")

burn.in <- floor(0.75*n.samples)

x11()
plot(m.1$p.theta.samples, density=FALSE)

##recover regression coefficients and random effects
m.1 <- spRecover(m.1, start=burn.in)

round(summary(m.1$p.theta.recover.samples)$quantiles[, c(3, 1, 5)], 2)
round(summary(m.1$p.beta.recover.samples)$quantiles[, c(3, 1, 5)], 2)

K2Cor <- function(K, m){
  tmp <- matrix(0, m, m)
  tmp[lower.tri(tmp, TRUE)] <- K
  tmp[upper.tri(tmp, FALSE)] <- t(tmp)[lower.tri(tmp, FALSE)]
  (cov2cor(tmp))[lower.tri(tmp, TRUE)]
}

cor.K.hat <- mcmc(t(apply(m.1$p.theta.recover.samples[,grep("K",colnames(m.1$p.theta.recover.samples))], 1, K2Cor, q)))

phi.hat <- m.1$p.theta.recover.samples[,grep("phi",colnames(m.1$p.theta.recover.samples))]

Psi.hat <- m.1$p.theta.recover.samples[,grep("Psi",colnames(m.1$p.theta.recover.samples))]

head(cor.K.hat)
summary(cor.K.hat)
             