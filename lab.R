if (!require("pacman")) install.packages("pacman")
pkgs = c("spBayes", "geoR", "PrevMap", "MBA", "sp", "mapview", "tmap") # package names
pacman::p_load(pkgs, character.only = T)


lf <- readRDS(file = "data/lf_clean.rds")

table(lf$Method_1, lf$Country)

country <- "Benin"

lf <- lf[lf$Country == country, ]

ict <- which(lf$Method_1 == "Serological")
mf <- which(lf$Method_1 == "Parasitological")

test <- sample(mf, floor(length(mf) * 0.3))
mf <- mf[!mf %in% test]

logit <- log((lf$Positive + 0.5)/(lf$Examined - lf$Positive + 0.5))
lf$logit <- logit
hist(lf$logit[ict])
hist(lf$logit[mf])

coords <- lf[, c("Longitude", "Latitude")]

lf_sp <- lf
coordinates(lf_sp) <- ~ Longitude + Latitude
proj4string(lf_sp) <- CRS("+init=epsg:3857 +units=km")
#mapView(lf_sp, zcol = "Method_1", legend = T, burst = T, hide = T)

africa <- rgdal::readOGR("data/Africa.shp")
proj4string(africa) <- CRS("+init=epsg:4326") 
africa <- spTransform(africa, proj4string(lf_sp))

lf_sp$Validation <- 1 
lf_sp$Validation[test] <- "Test"
lf_sp$Validation[-test] <- "Training"

map <- tm_shape(africa, bbox = country, is.master = T) +
  tm_borders("black") +
  tm_fill("grey85") + 
  tm_shape(africa[africa$COUNTRY == country, ]) +
  tm_fill("white") +
  tm_borders("black") +
  tm_shape(lf_sp) +
  tm_symbols(col = "Method_1", title.col = "Survey data by\ndiagnostic method", shape = "Validation", shapes = c(4, 21), 
             shapes.legend = 21, shapes.legend.fill = "red",
             size = .2, border.col = "black", palette = c("red", "blue")) +
  tm_compass(position = "left") +
  tm_layout(bg.color = "lightblue", design.mode = F, legend.bg.color = "white", legend.frame = "black") +
  tm_scale_bar() 
map      

# tmap_mode(mode = "view")
# map
# ttm()

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
n.samples <- 20000

starting <- list("phi" = c(1 / 34^-1, 1 / 9^-1), "A" = A.starting, "Psi" = rep(1, q))
tuning <- list("phi" = rep(0.5, q), "A" = c(0.01, 0.5, 0.01), "Psi" = rep(0.1, q))
priors <- list("phi.Unif" = list(c(1 / 5^-1, 1 / 2^-1), c(1 / 100^-1, 1 / 50^-1)),
               "K.IW" = list(q + 1, diag(c(1, 1))), "Psi.ig" = list(rep(2, q), c(.1, .1)))

system.time(m.1 <- spMisalignLM(list(logit[ict] ~ 1, logit[mf] ~ 1), 
                                coords = list(coords_jt[ict, ], coords_jt[mf, ]),
                                starting = starting, tuning = tuning, priors = priors, 
                                n.samples = n.samples, cov.model = "exponential", n.report = 1000))
#save(m.1, file = "output/sp_misalign_model.rds")

burn.in <- floor(0.75 * n.samples)

x11()
plot(m.1$p.theta.samples, density=FALSE)
plot(mcmc(m.1$p.theta.samples[burn.in:n.samples, ]), density=FALSE)

##recover regression coefficients and random effects
m.1 <- spRecover(m.1, start = burn.in)
#save(m.1, file = "output/sp_misalign_model_recovered.rds")

round(summary(m.1$p.theta.recover.samples)$quantiles[,c(3,1,5)],2)
round(summary(m.1$p.beta.recover.samples)$quantiles[,c(3,1,5)],2)

# Predict mf only at test points
n <- nrow(coords[mf, ])
pred.covars <- list(as.matrix(1), as.matrix(rep(1, n)))
pred.coords <- list(as.matrix(coords[1, ]), as.matrix(coords[mf, ]))

# Predict mf everywhere but not ict
n <- nrow(coords)
pred.covars <- list(as.matrix(1), as.matrix(rep(1, n)))
pred.coords <- list(as.matrix(coords[1, ]), as.matrix(coords))

# Predict both everywhere
n <- nrow(coords)
pred.covars <- list(as.matrix(rep(1, n)), as.matrix(rep(1, n)))
pred.coords <- list(as.matrix(coords), as.matrix(coords))

# Predict at new points
coords_new <- coordinates(spsample(africa[africa$COUNTRY == country, ], n = 1000, type = "regular"))
n <- nrow(coords_new)
pred.covars <- list(as.matrix(rep(1, n)), as.matrix(rep(1, n)))
pred.coords <- list(as.matrix(coords_new), as.matrix(coords_new))

pred <- spPredict(m.1, pred.covars = pred.covars, pred.coords = pred.coords, start = burn.in, thin = 10)


##summary and check
quants <- function(x){quantile(x, prob = c(0.5, 0.025, 0.975))}

y.hat <- apply(pred$p.y.predictive.samples, 1, quants)
y.hat <- apply(pred$p.y.predictive.samples, 1, median)
p.hat <- apply(psych::logistic(pred$p.y.predictive.samples), 1, median)


##unstack and plot
y.mf.hat <- y.hat[, 2:(n + 1)]

y.ict.hat <- y.hat[, 1:n]
y.mf.hat <- y.hat[, (n + 1):(2 * n)]

cbind(logit[ict], t(y.ict.hat[,ict]))
cbind(logit[mf], t(y.mf.hat[,mf]))

plot(logit[ict], y.ict.hat[1, ict], xlab = "Observed MF", ylab = "Predicted MF",
     main = "")
plot(logit[mf], y.mf.hat[1, mf], xlab = "Observed MF", ylab = "Predicted MF",
     main = "")
arrows(logit[mf], y.mf.hat[1, ], logit[mf], y.mf.hat[2, ], length = 0.02, angle = 90)
arrows(y.a[-sub.1], y.a.hat[1,-sub.1], y.a[-sub.1], y.a.hat[3,-sub.1], length = 0.02, angle = 90)
abline(0, 1, col = "red")


p.ict.hat <- p.hat[1:n]
p.mf.hat <- p.hat[(n + 1):(2 * n)]
surf_ict <- mba.surf(cbind(coords_new, p.ict.hat), no.X = 200, no.Y = 200)$xyz.est
fields::image.plot(surf_ict)
points(coords[ict, ])

surf_mf <- mba.surf(cbind(coords_new, p.mf.hat), no.X = 200, no.Y = 200)$xyz.est
fields::image.plot(surf_mf, col = )

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

# INLA ------------------------------------------------------------------------------------

### R code for Chapter 8
### Last update: 18/08/2014

###################################################
### Set working directory and load packages
###################################################
remove(list=ls())
my.dir <- paste(getwd(),"/",sep="")

require(INLA)
inla.setOption(scale.model.default=TRUE)

require(splancs)
require(lattice)
require(fields)
require(maptools)
require(spdep)


###################################################
### Code for Section 8.1
###################################################
n_y <- 15
n_x <- 20 

set.seed(2)
loc_y <- cbind(runif(n_y), runif(n_y))
set.seed(1)
loc_x <- cbind(runif(n_x), runif(n_x))

# *** Code for Figure 8.1
plot(loc_x, xlim=0:1, ylim=0:1, pch=19, xlab="", ylab="", asp=1, cex=1.5)
points(loc_y, pch=21, bg="white", col=1, lwd=2, cex=1.5)
# ***

###################################################
### Code for Section 8.1.1 
###################################################
#--- Data simulation ---#
kappa_xi <- 5
sigma2_xi <- 0.5
kappa_u <- 7
sigma2_u <- 0.3

simulate_GF <- function(coords, kappa, variance, lambda=1) {
  # Compute the number of locations
  n <- nrow(coords)
  # Compute the distance matrix  
  dist.m <- as.matrix(dist(coords))
  # Compute the Matern correlation matrix
  cor.m <- 2^(1-lambda)/gamma(lambda)*(kappa*dist.m)^lambda*
    besselK(x=dist.m*kappa, nu=lambda)
  diag(cor.m) <- 1
  # Compute the covariance matrix
  Sigma <- variance * cor.m
  # Simulate date using standard algorithm based on Cholesky fact.
  c(chol(Sigma) %*% rnorm(n=n,mean=0,sd=1))
} 

set.seed(223)
u <- simulate_GF(coords=loc_y, kappa=kappa_u, variance=sigma2_u)
length(u)

set.seed(233)
xi <- simulate_GF(coords=rbind(loc_x, loc_y), kappa=kappa_xi, variance=sigma2_xi)
length(xi)

b0 <- 10
beta1 <- -0.5
sigma2_e <- 0.16
sigma2_x <- 0.25

set.seed(4455)
x <- xi[1:n_x] + rnorm(n=n_x, mean=0, sd=sqrt(sigma2_x))
set.seed(5544)
y <- b0 + beta1*xi[n_x + 1:n_y] + u + rnorm(n=n_y, mean=0, sd=sqrt(sigma2_e))

#--- Model fitting ---#
mesh <- inla.mesh.2d(loc=rbind(loc_x, loc_y), max.edge=0.15, 
                     cutoff=0.03, offset=0.1)

# *** Code for Figure 8.2
plot(mesh, main="", asp=1)
points(loc_x, pch=19, cex=1.5)
points(loc_y, pch=21, bg="white", col=1, lwd=2, cex=1.5)
# ***

spde <- inla.spde2.matern(mesh=mesh, alpha=2)

A_x <- inla.spde.make.A(mesh=mesh, loc=loc_x)
dim(A_x)
A_y <- inla.spde.make.A(mesh=mesh, loc=loc_y)
dim(A_y)

stk_x <- inla.stack(data=list(y=cbind(x, NA)),  
                    effects=list(xi.field=1:spde$n.spde), 
                    A=list(A_x), 
                    tag="est.x")

stk_y <- inla.stack(data=list(y=cbind(NA, y)),  
                    effects=list(list(u.field=1:spde$n.spde, x.field=1:spde$n.spde), list(intercept=rep(1,n_y))),
                    A=list(A_y, 1),
                    tag="est.y")

stk <- inla.stack(stk_x, stk_y)


formula <- y ~ -1 + f(xi.field, model=spde) + 
  intercept + f(u.field, model=spde) +
  f(x.field, copy="xi.field", fixed=FALSE, hyper=list(theta=list(param=c(-1, 10))))

precprior <- list(theta=list(param=c(1, 0.1)))
output <- inla(formula, family=c("gaussian","gaussian"),
               data=inla.stack.data(stk),  
               control.predictor=list(compute=TRUE, A=inla.stack.A(stk)),
               control.family=list(list(hyper=precprior), list(hyper=precprior)))

#---- Getting the results ---#
round(cbind(True=b0,output$summary.fixed[1,1:5,drop=FALSE]), 4)
round(cbind(True=1/c(Prec.x=sigma2_x, Prec.e=sigma2_e), output$summary.hyperpar[1:2,1:5]), 4)
round(cbind(True=beta1, output$summary.hyperpar[7,1:5,drop=FALSE]), 4)

# *** Code for Figure 8.4
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(2,0.7,0))
plot(inla.smarginal(output$marginals.fixed[[1]]), type='l', xlab=expression(b[0]), ylab='Density')
abline(v=b0, col=gray(.4))

plot(inla.smarginal(output$marginals.hyperpar[[7]]), type='l', xlab=expression(beta[1]), ylab='Density')
abline(v=beta1, col=gray(.4))

plot.default(inla.tmarginal(function(x) 1/x, output$marginals.hyperpar[[1]]), type='l', xlab=expression(sigma[x]^2), ylab='Density')
abline(v=sigma2_x, col=gray(.4))

plot.default(inla.tmarginal(function(x) 1/x, output$marginals.hyperpar[[2]]), type='l', xlab=expression(sigma[e]^2), ylab='Density')
abline(v=sigma2_e, col=gray(.4))
# ***

xi_field <- inla.spde2.result(output, name="xi.field", spde)
u_field <- inla.spde2.result(output, name="u.field", spde)

round(cbind(True=c(sigma2.xi=sigma2_xi, sigma2.u=sigma2_u), 
            mean=c(inla.emarginal(function(x) x, xi_field$marginals.var[[1]]), 
                   inla.emarginal(function(x) x, u_field$marginals.var[[1]])), 
            rbind(inla.hpdmarginal(.95, xi_field$marginals.var[[1]]), 
                  inla.hpdmarginal(.95, u_field$marginals.var[[1]]))), 4)

round(cbind(True=c(range.xi=sqrt(8)/kappa_xi, range.u=sqrt(8)/kappa_u), 
            mean=c(inla.emarginal(function(x) x, xi_field$marginals.range[[1]]), 
                   inla.emarginal(function(x) x, u_field$marginals.range[[1]])), 
            rbind(inla.hpdmarginal(.95, xi_field$marginals.range[[1]]), 
                  inla.hpdmarginal(.95, u_field$marginals.range[[1]]))), 4)

# *** Code for Figure 8.5
par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(2,0.7,0))
plot.default(xi_field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[xi]^2), ylab='Density')
abline(v=sigma2_xi, col=gray(.4))

plot.default(u_field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[u]^2), ylab='Density')
abline(v=sigma2_u, col=gray(.4))

plot.default(xi_field$marginals.range.nominal[[1]], type='l', xlab=expression(r[xi]), ylab='Density')
abline(v=sqrt(8)/kappa_xi, col=gray(.4))

plot.default(u_field$marginals.range.nominal[[1]], type='l', xlab=expression(r[u]), ylab='Density')
abline(v=sqrt(8)/kappa_u, col=gray(.4))
# ***

###################################################
### Code for Section 8.1.2
###################################################
#--- Data simulation ---#
set.seed(134)
E <- rgamma(n=n_x, shape=10, rate=10) 
rho <- exp(xi[1:n_x])
set.seed(14)
x <- rpois(n=n_x, lambda=E*rho)

b0 <- 1
set.seed(19)
ntrials <- 5 + rpois(n=n_y, lambda=10)

eta_y <- b0 + beta1 * xi[n_x + 1:n_y] + u 
set.seed(553)
y <- rbinom(n=n_y, size=ntrials, prob=exp(eta_y)/(1 + exp(eta_y)) )

#--- Fitting the model ---#
stk_x <- inla.stack(data=list(y=cbind(x, NA), E=E, link="log"), 
                    effects=list(xi.field=1:spde$n.spde), 
                    A=list(A_x),
                    tag="est.x")

stk_y <- inla.stack(data=list(y=cbind(NA, y), Ntrials=ntrials, link="logit"), 
                    effects=list(list(x.field=1:spde$n.spde, u.field=1:spde$n.spde), list(intercept=rep(1,length(y)))), 
                    A=list(A_y, 1),
                    tag="est.y")

stk <- inla.stack(stk_x, stk_y)

output2 <- inla(formula,
                data=inla.stack.data(stk), 
                family=c("poisson", "binomial"), 
                E=inla.stack.data(stk)$E, 
                Ntrials=inla.stack.data(stk)$Ntrials, 
                control.predictor=list(compute=TRUE, A=inla.stack.A(stk)))

# *** Code for Figure 8.6
par(mfrow=c(1,2), mar=c(3,3,0.5,0.5), mgp=c(2,0.7,0))
plot(inla.smarginal(output2$marginals.fixed[[1]]), type='l', xlab=expression(b[0]), ylab='Density')
abline(v=b0, col=gray(.4))

plot(inla.smarginal(output2$marginals.hyperpar[[5]]), type='l', xlab=expression(beta[1]), ylab='Density')
abline(v=beta1, col=gray(.4))
# ***

# *** Code for Figure 8.7
xi_field <- inla.spde2.result(output2, "xi.field", spde)
u_field <- inla.spde2.result(output2, "u.field", spde)

par(mfrow=c(2,2), mar=c(3,3,0.5,0.5), mgp=c(2,0.7,0))
plot.default(xi_field$marginals.variance.n[[1]], type='l', 
             xlab=expression(sigma[xi]^2), ylab='Density')
abline(v=sigma2_xi, col=gray(.4))

plot.default(u_field$marginals.variance.n[[1]], type='l', xlab=expression(sigma[u]^2), ylab='Density')
abline(v=sigma2_u, col=gray(.4))

plot.default(xi_field$marginals.range.n[[1]], type='l', xlab=expression(r[xi]), ylab='Density')
abline(v=sqrt(8)/kappa_xi, col=gray(.4))

plot.default(u_field$marginals.range.n[[1]], type='l', xlab=expression(r[u]), ylab='Density')
abline(v=sqrt(8)/kappa_u, col=gray(.4))
# ***

###################################################
### Code for Section 8.2
###################################################
library(splancs) 
library(lattice)

#--- Parana' state rainfall data ---#
data(PRprec) 
dim(PRprec)

PRprec[1:3, 1:10] 
apply(X=PRprec[,4:10], MARGIN=2, FUN=summary) 
colSums(PRprec[,4:10]>0, na.rm=TRUE) 

jday <- 6 #columns in which data are stored
z <- as.numeric(PRprec[,jday] > 0) 
sum(z==0,na.rm=T) #stations with no rain
y <- ifelse(PRprec[,jday] > 0, PRprec[,jday], NA) 

# *** Code for Figure 8.8
data(PRborder)
par(mfrow=c(1,2), mar=c(0,0,0,0))
color.map <- gray(c(0.9, 0.7, 0.5, 0.3, 0))
iz <- which(z==1)
iz.o <- iz[order(y[iz])]
rain.factor <- cut(y[iz.o], quantile(y[iz.o], 0:5/5, na.rm=TRUE))
rcol <- color.map[unclass(rain.factor)]
# Points with some rain
plot(PRprec[,1:2], axes=FALSE, xlab="",ylab="",asp=1, cex=0.5)
points(PRprec[iz.o,1:2], pch=21, bg=rcol, cex=0.5+unclass(rain.factor)/5)
# Points with no rain
lines(PRborder)
legend("topright", pch=c(1,rep(21,5)), pt.cex=0.5+0:5/5, bty="n", 
       pt.bg=c(1, color.map), legend=c("0", levels(rain.factor)))

par(mar=c(2,4,1,0.5), mgp=c(2,0.7,0))
hist(y, xlab=y, col=gray(.9), main='')
# ***

#--- Fitting the model ---#
prec.coords <- cbind(PRprec$Longitude, PRprec$Latitude)
boundary <- inla.nonconvex.hull(points=prec.coords, convex=0.2, concave=0.2)
mesh <- inla.mesh.2d(loc=prec.coords, boundary=boundary, max.edge=c(0.3, 0.8), cutoff=0.1)

# *** Code for Figure 8.9
par(mar=c(0,0,0,0))
plot(mesh, asp=1,main="")
points(PRprec$Longitude, PRprec$Latitude, pch=19, cex=0.5)
# ***

spde <- inla.spde2.matern(mesh=mesh, alpha=2) 
A <- inla.spde.make.A(mesh=mesh, loc=prec.coords) 
dim(A)

stk.y <- inla.stack(data=list(amount=y, #for the single model
                              alldata=cbind(y, NA)), #for the joint model
                    A=list(A, 1), 
                    effects=list(list(y.field=1:spde$n.spde), list(y.intercept=rep(1,length(y)))),
                    tag="est.y") 

stk.z <- inla.stack(data=list(occurence=z, #for the single model
                              alldata=cbind(NA, z)), #for the joint model
                    A=list(A, 1), 
                    effects=list(list(z.field=1:spde$n.spde, zc.field=1:spde$n.spde), list(z.intercept=rep(1,length(z)))),
                    tag="est.z") 

formula.y <- amount ~ -1 +  y.intercept + f(y.field, model=spde)
out.y <- inla(formula.y, family="gamma", data=inla.stack.data(stk.y),
              control.predictor=list(A=inla.stack.A(stk.y)), 
              control.compute=list(dic=TRUE), 
              control.inla=list(strategy="laplace"))

formula.z <- occurence ~ -1 + z.intercept + f(z.field, model=spde)
out.z <- inla(formula.z, family="binomial", data=inla.stack.data(stk.z),
              control.predictor=list(A=inla.stack.A(stk.z)), 
              control.compute=list(dic=TRUE), 
              control.inla=list(strategy="laplace"))

stk.yz <- inla.stack(stk.y, stk.z) 

formula.yz <- alldata ~ -1 + y.intercept + z.intercept + 
  f(y.field, model=spde) + f(z.field, model=spde) + 
  f(zc.field, copy="y.field", fixed=FALSE)
out.yz <- inla(formula.yz, family=c("gamma", "binomial"), 
               data=inla.stack.data(stk.yz), 
               control.predictor=list(A=inla.stack.A(stk.yz)),
               control.compute=list(dic=TRUE, config=TRUE), 
               control.inla=list(strategy="laplace"))

rbind(separate=c(y=out.y$dic$dic, z=out.z$dic$dic), 
      joint=tapply(out.yz$dic$local.dic, out.yz$dic$family, sum))

idy <- which(PRprec[,jday]>0)
exp.y <- sapply(out.yz$marginals.linear.pred[idy], function(m) inla.emarginal(exp, m))
idz <- which(!is.na(z))
exp.z <- sapply(out.yz$marginals.linear.pred[nrow(PRprec) + idz], function(m) inla.emarginal(inla.link.invlogit, m))

c(yPositive.obs=mean(PRprec[which(PRprec[,jday]>0),jday], na.rm=TRUE), yPositive.pred=mean(exp.y))
c(z.obs=mean(z, na.rm=TRUE), z.pred=mean(exp.z))

# *** Code for Figure 8.10
par(mfrow=c(2,2), mar=c(3,3,1,1), mgp=2:0)
plot(inla.smarginal(out.yz$marginals.fixed[[1]]), type='l', ylab='Density', 
     ylim=c(0,max(out.yz$marginals.fixed[[1]][,2], out.y$marginals.fixed[[1]][,2])), xlab=expression(b[0]^y))
lines(inla.smarginal(out.y$marginals.fixed[[1]]), lty=2)
legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

plot(inla.smarginal(out.yz$marginals.fixed[[2]]), type='l', ylab='Density', 
     ylim=c(0,max(out.yz$marginals.fixed[[2]][,2],out.z$marginals.fixed[[1]][,2])), xlab=expression(b[0]^z))
lines(inla.smarginal(out.z$marginals.fixed[[1]]), lty=2)
legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

plot(inla.smarginal(out.yz$marginals.hyperpar[[1]]), type='l', ylab='Density', 
     ylim=c(0,max(out.yz$marginals.hy[[1]][,2], out.y$marginals.hy[[1]][,2])), xlab=expression(tau))
lines(inla.smarginal(out.y$marginals.hyperpar[[1]]),lty=2)
legend("topright", lty=c(1,2), legend=c("joint", "separate"), box.lty=0)

plot(inla.smarginal(out.yz$marginals.hyperpar[[6]]), type='l', ylab='Density',  
     ylim=c(0,max(out.yz$marginals.hy[[6]][,2])),  xlab=expression(beta[1]))
# ***

# *** Code for Figure 8.11
data(PRborder)
in.pr <- which(inout(mesh$loc, PRborder))
par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(2,1,0))
ordxy <- intersect(order(out.yz$summary.random$y.field$mean), in.pr)
plot(out.yz$summary.random$y.field$mean[ordxy], type="l", ylab=expression(xi[i]), ylim=range(out.yz$summary.random$y.field[, 4:6]))
for (i in c(4,6)) lines(out.yz$summary.random$y.field[ordxy,i], lty=2) 
abline(h=0, lty=3) 
ordxz <- intersect(order(out.yz$summary.random$z.field$mean), in.pr)
plot(out.yz$summary.random$z.field$mean[ordxz], type="l", ylab=expression(u[i]), ylim=range(out.yz$summary.random$z.field[, 4:6]))
for (i in c(4,6)) lines(out.yz$summary.random$z.field[ordxz,i], lty=2) 
abline(h=0, lty=3) 
# ***

#--- Significance of the spatial effects ---#
c(n.nodes.in.pr=length(in.pr), expected.out=0.05*length(in.pr), 
  observed.out=sum(out.yz$summary.random$y.field[in.pr,4]>0 | out.yz$summary.random$y.field[in.pr,6]<0))

form.nospatial <- alldata ~ -1 + z.intercept + y.intercept
out.nospatial <- inla(form.nospatial, family=c("gamma", "binomial"), 
                      data=inla.stack.data(stk.yz),
                      control.predictor=list(A=inla.stack.A(stk.yz)),
                      control.compute=list(dic=TRUE), 
                      control.inla=list(strategy="laplace"))

form.oneshared <- alldata ~ -1 + y.intercept + z.intercept + 
  f(y.field, model=spde) + f(z.field, model=spde) 
out.oneshared <- inla(form.oneshared, family=c("gamma", "binomial"), 
                      data=inla.stack.data(stk.yz), 
                      control.predictor=list(A=inla.stack.A(stk.yz)),
                      control.compute=list(dic=TRUE), 
                      control.inla=list(strategy="laplace"))

rbind(nospatial=tapply(out.nospatial$dic$local.dic, out.nospatial$dic$family, sum), 
      one.shared=tapply(out.oneshared$dic$local.dic, out.oneshared$dic$family, sum),
      separated=c(out.y$dic$dic, out.z$dic$dic),
      joint=tapply(out.yz$dic$local.dic, out.yz$dic$family, sum))

#--- The spatial random effects ---#
y.field <- inla.spde2.result(out.yz, "y.field", spde) 
z.field <- inla.spde2.result(out.yz, "z.field", spde) 

# *** Code for Figure 8.12
par(mfrow=c(2,2))
plot.default(y.field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[xi]^2), ylab='Density')
plot.default(y.field$marginals.range.nominal[[1]], type='l', xlab=expression(r[xi]), ylab='Density')
plot.default(z.field$marginals.variance.nominal[[1]], type='l', xlab=expression(sigma[u]^2), ylab='Density')
plot.default(z.field$marginals.range.nominal[[1]], type='l', xlab=expression(r[u]), ylab='Density')
# ***

data(PRborder) # state boundary
nxy <- round(c(diff(range(PRborder[,1])), diff(range(PRborder[,2]))) /.02) 
nxy # grid size

projgrid <- inla.mesh.projector(mesh, xlim=range(PRborder[,1]), ylim=range(PRborder[,2]), dims=nxy) 
rf.grid <- list(xi.mean = inla.mesh.project(projgrid, out.yz$summary.random$y.field$mean)) 
rf.grid$xi.sd <- inla.mesh.project(projgrid, out.yz$summary.random$y.field$sd) 
rf.grid$u.mean <- inla.mesh.project(projgrid, out.yz$summary.random$z.field$mean) 
rf.grid$u.sd <- inla.mesh.project(projgrid, out.yz$summary.random$z.field$sd) 

# *** Code for Figure 8.13
xy.in <- inout(projgrid$lattice$loc, PRborder)
for (i in 1:4) rf.grid[[i]][!xy.in] <- NA 

library(fields)
par(mfcol=c(2,2), mar=c(0,0,0,0))
for (i in 1:4) {
  image.plot(list(x=projgrid$x, y=projgrid$y, z=rf.grid[[i]]), 
             asp=1, axes=FALSE, nlevel=101, col=gray(100:0/100), 
             legend.mar=7, legend.args=list(text=names(rf.grid)[i], line=0))
  lines(PRborder)
}
# ***

#--- Prediction of the responses ---#
s1 <- inla.posterior.sample(n=1,result=out.yz)
names(s1[[1]])
grep("y.intercept", rownames(s1[[1]]$latent), fixed=TRUE)

ids <- lapply(c('y.intercept', 'z.intercept', 'y.field', 'z.field', 'zc.field'), function(x) grep(x, rownames(s1[[1]]$latent), fixed=TRUE))
pred.y.f <- function(s) exp(s$latent[ids[[1]], 1] + s$latent[ids[[3]], 1])
pred.z.f <- function(s) 1/(1 + exp(-(s$latent[ids[[2]], 1] + s$latent[ids[[4]], 1] + s$latent[ids[[5]], 1])))

s1000 <- inla.posterior.sample(1000, out.yz)

prd.y <- sapply(s1000, pred.y.f)
prd.z <- sapply(s1000, pred.z.f)

prd <- list(y.mean=inla.mesh.project(projgrid, field=rowMeans(prd.y)))
prd$y.sd <- inla.mesh.project(projgrid, field=apply(prd.y, 1, sd))
prd$z.mean <- inla.mesh.project(projgrid, field=rowMeans(prd.z))
prd$z.sd <- inla.mesh.project(projgrid, field=apply(prd.z, 1, sd))
for (j in 1:4) prd[[j]][!xy.in] <- NA

# *** Code for Figure 8.14
par(mfcol=c(2,2), mar=c(0,0,0,0))
for (i in 1:4) {
  image.plot(list(x=projgrid$x, y=projgrid$y, z=prd[[i]]), 
             asp=1, axes=FALSE, nlevel=101, col=gray(100:0/100), 
             legend.mar=7, legend.args=list(text=names(prd)[i], line=0))
  lines(PRborder)
}
# ***

###################################################
### Code for Section 8.3.1
###################################################
#--- Data simulation ---#
library(spdep)
file <- system.file("etc/shapes/eire.shp", package="spdep")[1]
require(maptools)
pol <- readShapePoly(file)

nbl <- poly2nb(pol)
n <- length(nbl)
R <- Matrix(nb2mat(nbl, style="B"))

mu.xi <- c(-10, 0.5);        tau <- c(0.5, 0.3) 
ddiag <- c(0.3, 0.5);       gamma <- c(.8, .6) 

n <- 26
nnb <- colSums(R) # number of neighbors
n <- length(nnb)

Q1bp <- tau[1]*(diag(ddiag[1],n)+diag(nnb)-R)
Q2bp <- tau[2]*(diag(ddiag[2],n)+diag(nnb)-R)

inla.setOption(inla.call=loc.call)
T <- 30
xi1 <- omega1 <- inla.qsample(n=T, Q=Q1bp, seed=1)
xi2 <- omega2 <- inla.qsample(n=T, Q=Q2bp, seed=2)
xi1[,1] <- sqrt(1-gamma[1]^2)*omega1[,1]
xi2[,1] <- sqrt(1-gamma[2]^2)*omega2[,1]
for (j in 2:T) {
  xi1[, j] <- gamma[1]*xi1[, j-1] + sqrt(1-gamma[1]^2)*omega1[, j]
  xi2[, j] <- gamma[2]*xi2[, j-1] + sqrt(1-gamma[2]^2)*omega2[, j]
}

# Covariate matrix
set.seed(3)
x2 <- matrix(runif(T*n), nrow=n)

# Expected values
set.seed(4)
ee <- outer(rgamma(n, 5000000, scale=0.1), rgamma(T, 100, scale=0.01)) * rgamma(n*T, 100, scale=0.01)

# Observations
set.seed(5)
yy <- rpois(n=n*T, lambda=ee*exp((xi1+mu.xi[1]) + (xi2+mu.xi[2])*x2))

#--- Estimation throught the group feature ---#
dat.bp <- data.frame(y=yy, x=as.vector(x2),  e=as.vector(ee), xi1=rep(1:n, T), index.t=rep(1:T, each=n)) 
dat.bp$xi2 <- dat.bp$xi1
dat.bp$index.t2 <- dat.bp$index.t 

form.bp <- y ~ x + 
  f(xi1, model="besagproper", graph=R, group=index.t, control.group=list(model="ar1")) +
  f(xi2, x, model="besagproper", graph=R, group=index.t2, control.group=list(model="ar1"))

res.bp <- inla(form.bp, family="poisson", data=dat.bp, E=dat.bp$e, control.predictor=list(compute=TRUE))

# *** Code for Figure 8.15
par(mfrow=c(4,2), mar=c(3,3,1,0.5), mgp=c(2,0.7,0), las=1)
plot(res.bp$marginals.fixed$'(Intercept)', type='l', ylab='Density', xlab=expression(mu[xi[1]]))
abline(v=mu.xi[1], col=gray(.4))
plot(res.bp$marginals.fixed$'x', type='l', ylab='Density', xlab=expression(mu[xi[2]]))
abline(v=mu.xi[2], col=gray(.4))
names <- c(expression(tau[1]),expression(d[1]),expression(gamma[1]),
           expression(tau[2]),expression(d[2]),expression(gamma[2]))
for (j in c(1,4,2,5,3,6)) {
  plot(res.bp$marginals.hyperpar[[j]], type='l', ylab='Density', xlab=names[j])
  abline(v=c(tau, ddiag, gamma)[c(1,3,5,2,4,6)][j], col=gray(.4))
}
# ***

# *** Code for Figure 8.16
set.seed(1); id.areas.sel <- sample(1:n,3) 
par(mfcol=c(3,3), mar=c(1.5,2,1.5,.5), mgp=c(1,.5,0), las=1)
for (i in id.areas.sel) {
  i.aux <- which(dat.bp$xi1==i)
  yran <- range(res.bp$summary.random$xi1[i.aux,4:6])
  plot.ts(xi1[i, ], lwd=2, col=gray(.4), ylim=yran, ylab='', xlab='', main=paste('xi1 in area', i))
  for (j in c(2, 4, 6))
    lines(matrix(res.bp$summary.random$xi1[,j], n)[i,], lty=c(1,2,2)[j/2])
  yran <- range(res.bp$summary.random$xi2[i.aux,4:6])
  plot.ts(xi2[i, ], lwd=2, col=gray(.4), ylim=yran, ylab='', xlab='', main=paste('xi2 in area', i))
  for (j in c(2, 4, 6))
    lines(matrix(res.bp$summary.random$xi2[,j], n)[i,], lty=c(1,2,2)[j/2])
  yran <- range(res.bp$summary.fitted.val[i.aux,4:6]*dat.bp$e, yy[i.aux])
  plot.ts(matrix(yy,n)[i, ], lwd=2, col=gray(.4), ylim=yran, ylab='', xlab='', main=paste('y in area', i))
  for (j in c(1, 3, 5))
    lines(matrix(res.bp$summary.fitted.val[,j]*dat.bp$e, n)[i,], lty=c(1,2,2)[(j+1)/2])
}
# ***

# *** Code for Figure 8.17
T.2 <- round(T/2);   j0 <- ncol(pol@data)
colnames(xi1) <- colnames(xi2) <- 1:T 
pol@data <- data.frame(
  pol@data, xi1=xi1,  xi1pm=
    matrix(res.bp$summary.random$xi1$mean, n), 
  xi2=xi2, xi2pm=matrix(res.bp$summary.random$xi2$mean, n))
j.aux <- j0 + rep(c(1, T.2, T), each=4) + rep(0:3,3)*T

spplot(pol, names(pol)[j.aux], col.regions=grey(seq(0,1,0.01)), names.attr=list("xi1,time1,true", "xi1,time1,est","xi2,time1,true","xi2,time1,est","xi1,time15,true","xi1,time15,est","xi2,time15,true","xi2,time15,est", "xi1,time30,true","xi1,time30,est","xi2,time30,true","xi2,time30,est")) 
# ***

c(xi1=cor(as.vector(xi1), res.bp$summary.random$xi1$mean), 
  xi2=cor(as.vector(xi2), res.bp$summary.random$xi2$mean))

###################################################
### Code for Section 8.3.2
###################################################
library(spdep)
require(maptools)

set.seed(1)
id.areas.sel <- sample(1:n,3) 

file <- system.file("etc/shapes/eire.shp", package="spdep")[1]
pol <- readShapePoly(file)
require(spdep)
nbl <- poly2nb(pol)
n <- length(nbl)
R <- Matrix(nb2mat(nbl, style="B"))
mu.xi <- c(-10, 0.5);   gamma <- c(.8, .6) 
nnb <- colSums(R) # number of neighbors
n <- length(nnb); T <- 30

# Covariate matrix
set.seed(3)
x2 <- matrix(runif(T*n), nrow=n)

# Expected values
set.seed(4)
ee <- outer(rgamma(n, 5000000, 10), rgamma(T, 100, 100)) *rgamma(n*T, 100, 100)

phi <- c(0.95, 0.85)
tau <- c(3, 5)

lamb.max <- max(eigen(R, only.values=TRUE)$values)
Q1g1 <- tau[1]*(diag(n)-phi[1]/lamb.max*R)
Q2g1 <- tau[2]*(diag(n)-phi[2]/lamb.max*R)

xi1 <- omega1 <- inla.qsample(n=T, Q=Q1g1, seed=1)
xi2 <- omega2 <- inla.qsample(n=T, Q=Q2g1, seed=2)
for (j in 2:T) {
  xi1[, j] <- gamma[1]*xi1[, j-1] + omega1[, j]
  xi2[, j] <- gamma[2]*xi2[, j-1] + omega2[, j]
}
set.seed(5)
yy <- rpois(n=n*T,lambda=ee*exp((mu.xi[1]+xi1) + (mu.xi[2]+xi2)*x2))

#--- The data augmentation strategy ---#
stc.mat <- kronecker(diag(T), R)
dat.aug <- list(e=c(as.vector(ee), rep(NA, n*T*2)), 
                x1=rep(1:0, c(n*T, n*T*2)), 
                x2=c(as.vector(x2), rep(0, n*T*2)))
dat.aug$Y <- cbind(c(yy, rep(NA, n*T*2)), rep(c(NA, 0), c(n*T, n*T*2)))
dat.aug$o1 <- c(rep(NA, n*T), 1:(n*T), rep(NA, n*T)) 
dat.aug$o2 <- c(rep(NA, n*T*2), 1:(n*T)) 

dat.aug$xi1 <- c(1:(n*T), 1:(n*T), rep(NA, n*T))
dat.aug$xi1d <- c(rep(NA, n*T+n), 1:(n*(T-1)), rep(NA, n*T))
dat.aug$xi2 <- c(1:(n*T), rep(NA, n*T), 1:(n*T))
dat.aug$xi2d <- c(rep(NA, n*T*2+n), 1:(n*(T-1)))

dat.aug$wo1 <- dat.aug$wxi1d <- rep(c(0,-1,0), c(n*T, n*T, n*T)) 
dat.aug$wo2 <- dat.aug$wxi2d <- rep(c(0,0,-1), c(n*T, n*T, n*T)) 

dat.aug$wxi1 <- rep(1:0, c(n*T*2, n*T))
dat.aug$wxi2 <- c(as.vector(x2), rep(0:1, c(n*T, n*T)))

form.aug <- Y ~ -1 + x1 + x2 + 
  f(o1, wo1, model="generic1", Cmatrix=stc.mat) + 
  f(o2, wo2, model="generic1", Cmatrix=stc.mat) + 
  f(xi1, wxi1, model="iid", hyper=list(theta=list(initial=-10, fixed=TRUE))) + 
  f(xi1d, wxi1d, copy="xi1", range=c(0,1), hyper=list(beta=list(fixed=FALSE, param=c(0, 1.5)))) + 
  f(xi2, wxi2, model="iid", hyper=list(theta=list(initial=-10, fixed=TRUE))) + 
  f(xi2d, wxi2d, copy="xi2", range=c(0,1), hyper=list(beta=list(fixed=FALSE, param=c(0, 1.5)))) 

# *** Code for Figure 8.19
theta.prior.sample <- rnorm(10000, 0, 1.5)
beta.prior.sample <- 1/(1+exp(-theta.prior.sample))
par(mar=c(3,3,0.5,0.5), mgp=c(2,1,0))
hist(beta.prior.sample, main='', ylab='', freq=F, xlab=expression(beta)) 
box(); axis(1)
# ***

res.aug <- inla(form.aug, c("poisson", "gaussian"), 
                data=dat.aug, E=dat.aug$e, 
                control.family=list(list(), list(hyper=list(theta=list(initial=10, fixed=TRUE)))), 
                control.predictor=list(compute=TRUE, link=1))

# *** Code for Figure 8.20
par(mfrow=c(4,2), mar=c(3,3,1,0.5), mgp=c(2,0.7,0), las=1)
plot(res.aug$marginals.fixed$'x1', type='l', xlab=expression(mu[xi[1]]))
abline(v=mu.xi[1], col=gray(.4))
plot(res.aug$marginals.fixed$'x2', type='l', xlab=expression(mu[xi[2]]))
abline(v=mu.xi[2], col=gray(.4))
h.names <- c(expression(tau[1]),expression(tau[2]),expression(phi[1]),
             expression(phi[2]),expression(gamma[1]),expression(gamma[2]))
for (j in 1:6) {
  plot(res.aug$marginals.hyperpar[[c(1,3,2,4,5,6)[j]]], 
       type='l', xlab=h.names[j])
  abline(v=c(tau, phi, gamma)[j], col=gray(.4))
}
# ***

# *** Code for Figure 8.21
par(mfcol=c(3,3), mar=c(1.5,2,1.5,.5), mgp=c(1,.5,0), las=1)
for (i in id.areas.sel) {
  i.aux <- which(rep(1:n, T)==i)
  yran <- range(res.aug$summary.random$xi1[i.aux,4:6], xi1[i,])
  plot.ts(xi1[i, ], lwd=2, col=gray(.4), ylim=yran, ylab='', xlab='', main=paste('xi1 in area', i))
  for (j in c(2, 4, 6))
    lines(matrix(res.aug$summary.random$xi1[,j], n)[i,], lty=c(1,2,2)[j/2])
  yran <- range(res.aug$summary.random$xi2[i.aux,4:6], xi2[i,])
  plot.ts(xi2[i, ], lwd=2, col=gray(.4), ylim=yran, ylab='', xlab='', main=paste('xi2 in area', i))
  for (j in c(2, 4, 6))
    lines(matrix(res.aug$summary.random$xi2[,j], n)[i,], lty=c(1,2,2)[j/2])
  yran <- range(res.aug$summary.fitted.val[i.aux,4:6]*ee[i, ], matrix(yy, n)[i,])
  plot.ts(matrix(yy,n)[i, ], lwd=2, col=gray(.4), ylim=yran, ylab='', xlab='', main=paste('y in area', i))
  for (j in c(1, 3, 5))
    lines(matrix(res.aug$summary.fitted.val[,j]*ee[i, ], n)[i,], lty=c(1,2,2)[(j+1)/2])
}
# ***

# *** Code for Figure 8.22
colnames(xi1) <- colnames(xi2) <- 1:T 
pol@data <- data.frame(
  pol@data[, 1:j0], xi1=xi1,  x1pm=
    matrix(res.aug$summary.random$xi1$mean, n), 
  xi2=xi2, x2pm=
    matrix(res.aug$summary.random$xi2$mean, n))
j.aux <- j0 + rep(c(1, T.2, T), each=4) + rep(0:3,3)*T
require(lattice)
spplot(pol, names(pol)[j.aux], col.regions=grey(seq(0,1,0.01)), names.attr=list("xi1,time1,true", "xi1,time1,est","xi2,time1,true","xi2,time1,est","xi1,time15,true","xi1,time15,est","xi2,time15,true","xi2,time15,est", "xi1,time30,true","xi1,time30,est","xi2,time30,true","xi2,time30,est")) 
# ***

c(xi1=cor(as.vector(xi1), res.aug$summary.random$xi1$mean), 
  xi2=cor(as.vector(xi2), res.aug$summary.random$xi2$mean))

###################################################
### Code for Section 8.4
###################################################
obst <- c(1:2,4:7,10)
n <- length(obst)
knots <- c(2, 5, 8)
k <- length(knots)

mesh1d <- inla.mesh.1d(loc=knots)    
wmat <- inla.spde.make.A(mesh=mesh1d, loc=obst)
wmat

# *** Code for Figure 8.23
par(mar=c(3,3,0.5,0.5), mgp=c(2,.7,0), las=1)
plot(knots, rep(max(obst),k), xlim=c(min(obst), max(obst+0.3)), ylim=c(min(obst), max(obst+1)), 
     type='h', axes=FALSE, col=1,  ylab='knots', xlab='observed time points')
axis(1, obst); axis(2, knots)
for (j in 1:k) 
  text(obst, knots[j], format(wmat[,j],,2), col=1*(wmat[,j]>0))
# ***

n <- 1000
obst <- 1:n
k <- 200
knots <- seq(n/k, n, n/k) - 0.5*n/k

mesh1d <- inla.mesh.1d(loc=knots)
wmat <- inla.spde.make.A(mesh=mesh1d, loc=obst)
dim(wmat)

taue <- 5;     taux <- 3;     gamma <- 0.7
set.seed(1)
x <- as.vector(arima.sim(n=k, model=list(ar=gamma), sd=sqrt(1/taux)))
e <- rnorm(n=n, mean=0, sd=sqrt(1/taue))
y <- drop(wmat%*%x) + e

formula <- y ~ -1 + f(i, model="ar1")
ar1 <- inla(formula, data=list(y=y, i=1:n), 
            control.predictor=list(compute=TRUE),
            control.compute=list(dic=TRUE))

id.mid.knots <- findInterval(obst, c(0.5, knots+n/k))
ar1k <- inla(formula, data=list(y=y, i=id.mid.knots), 
             control.predictor=list(compute=TRUE),
             control.compute=list(dic=TRUE))

ar1w <- inla(formula, data=list(y=y, i=1:k), 
             control.predictor=list(A=wmat, compute=TRUE),
             control.compute=list(dic=TRUE))

c(ar1$dic$dic, ar1k$dic$dic, ar1w$dic$dic)

# *** Code for Figure 8.24
plot(knots, x, type='l', lwd=5, xlim=c(1,50), ylim=range(y[1:50]), col=gray(.6),ylab="")
lines(obst, y, lwd=3)
for (j in 3:5) {
  lines(ar1$summary.fitted.values[,j], lwd=1, lty=2)
  lines(ar1k$summary.fitted.values[,j], lwd=2, lty=2)
  lines(ar1w$summary.fitted.values[,j], lwd=3, lty=2)
}
legend('bottomleft', c("x_true", "y_obs", "ar1","ar1k","ar1w"), 
       lwd=c(5,3,1:3), lty=c(1,1,2,2,2), bty='n', col=c(gray(.6), 1,1,1,1))
# ***

###################################################
### Code for Section 8.4.1
###################################################
file <- system.file("etc/shapes/eire.shp", package="spdep")[1]
pol <- readShapePoly(file)
nbl <- poly2nb(pol) 
nbmat <- Matrix(nb2mat(nbl, style="B"))
n <- nrow(nbmat)
T <- 20
nnb <- colSums(nbmat)
tau <- c(1, 5)
gamma <- 0.7
Qbp <- tau[1]*(diag(nnb+1e-5) - nbmat) 
xi <- omega <- inla.qsample(n=T, Q=Qbp, seed=1, constr=list(A=matrix(1,1,n), e=0))
xi[,1] <- sqrt(1-gamma^2)*omega[,1]
for (j in 2:T)  
  xi[, j] <- gamma*xi[, j-1] + sqrt(1-gamma^2)*omega[, j]

knots <- seq(from=2, to=T*4-2, by=4)
length(knots)
obst <- 1:(T*4)
k <- length(obst)
mesh1d <- inla.mesh.1d(loc=knots)
wmat <- inla.spde.make.A(mesh=mesh1d, loc=obst)
dim(wmat)

yy <- matrix(NA, nrow=nrow(xi), ncol=length(obst))
set.seed(2)
for (i in 1:n) 
  yy[i, ] <- drop(wmat %*% xi[i,]) + rnorm(n=k, mean=0, sd=sqrt(1/tau[2]))

ldat <- list(y.resp=as.vector(t(yy)), xi.idx=rep(1:n, each=T), time=rep(1:T, n))
length(ldat$y.resp)
length(ldat$xi.idx)
length(ldat$time)

formula <- y.resp ~ -1 + 
  f(xi.idx, model="besag", graph=nbmat, scale.model=FALSE,group=time, control.group=list(model="ar1")) 

A <- kronecker(Diagonal(n), wmat)
out.st.low <- inla(formula, data=ldat, control.predictor=list(A=A, compute=TRUE))

# *** Code for Figure 8.25
par(mfrow=c(1,3), mar=c(2.5, 2.5, 0.5, 0.5), mgp=c(1.7,.5,0))
j=1
plot(out.st.low$marginals.hy[[j]], type='l', xlab=expression(tau[e]), ylab='Density')
abline(v=c(tau[2:1], gamma)[j], col=gray(.4))
j=2
plot(out.st.low$marginals.hy[[j]], type='l', xlab=expression(tau[x]), ylab='Density')
abline(v=c(tau[2:1], gamma)[j], col=gray(.4))
j=3
plot(out.st.low$marginals.hy[[j]], type='l', xlab=expression(rho), ylab='Density')
abline(v=c(tau[2:1], gamma)[j], col=gray(.4))
#***

# *** Code for Figure 8.26
par(mfcol=c(2,3), mar=c(2,2,2,.5), mgp=c(1,.5,0), las=1)
for (i in sample(1:n,3)) {
  yran <- range(matrix(out.st.low$summary.random$xi[, 4], n)[i, ], matrix(out.st.low$summary.random$xi[, 6], n)[i, ])
  plot.ts(xi[i, ], lwd=2, col=gray(.4), ylim=yran, ylab='', xlab='', main=paste('xi at area', i))
  for (j in c(2, 4, 6))
    lines(matrix(out.st.low$summary.random$xi[,j], n)[i,], lty=c(2,2,2)[j/2])
  i.aux <- which(rep(1:n, each=k)==i)
  yran <- range(out.st.low$summary.fitted.val[i.aux, 3:5], matrix(yy,n)[i,])
  plot.ts(matrix(yy,n)[i, ], lwd=2, col=gray(.4), ylim=yran, ylab='', xlab='', main=paste('y at area', i))
  for (j in c(1, 3, 5))
    lines(out.st.low$summary.fitted.val[(1:nrow(A))[i.aux], j], lty=c(2,2,2)[(j+1)/2])
}
# ***

cor(as.vector(xi), out.st.low$summary.random$xi$mean)

# *** Code for Figure 8.27
T.2 <- round(T/2);   j0 <- ncol(pol@data);   colnames(xi) <- 1:T
pol@data <- data.frame(
  pol@data, x=xi,  xpm=matrix(out.st.low$summary.random$xi$mean, n))
j.aux <- j0 + c(1, T.2, T) + rep(0:1, each=3)*T

spplot(pol, names(pol)[j.aux],col.regions=grey(seq(0,1,0.01))) 
# ***