library(zipfR)   # package for incomplete gamma function
library(nloptr)  # package for optimization
library(KMsurv)  # package for breast-fed data
library(survidm) # package for German Breast Cancer Study Data.
library(Deriv)   # package for derivative
library(prodlim)
library(pec) # package For ipcw
library(stats)
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
## basic properties for the reference distributions ##
##+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##
# Standard Normal
f0_AN <- function(s){1/sqrt(2*pi)*exp(-s^2/2)} # density function of N(0,1)
F0_AN <- function(s){pnorm(s, mean = 0,sd = 1)} # distribution function of N(0,1)
QF_AN <- function(tau){qnorm(tau, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)}
S0_AN <- function(tau){qnorm(tau, mean = 0, sd = 1, lower.tail = FALSE, log.p = FALSE)}
# Standard Laplace
f0_ALa <- function(s){0.5*exp(-abs(s))} # density function of Laplace(0,1)
F0_ALa <- function(s){0.5+0.5*sign(s)*(1-exp(-abs(s)))} # distribution function of Laplace(0,1)
QF_ALa <- function(tau){-sign(tau-0.5)*log(1-2*abs(tau-0.5))} # quantile
S0_ALa <- function(s){ifelse(s<0, 1-0.5*exp(s), 0.5*exp(-s))}# survival
#Standard Logistic
f0_AL <- function(s){exp(-s)/(1+exp(-s))^2}
F0_AL <- function(s){1/(1+exp(-s))}
QF_AL <- function(tau){log(tau/(1-tau))}
S0_AL <- function(s){exp(-s)/(1+exp(-s))}
# Standard student's -t distribution
f0_st <- function(s, nu){dt(s,df = nu)}
F0_st <- function(s, nu){pt(s, df = nu, lower.tail = TRUE)}
QF_st <- function(tau, nu){qt(tau, df = nu)}
# Link functions and its inverse
glog <- function(s){log(s)}
glog.inv <- function(s){exp(s)}
glogit <- function(s, lambda){log(exp(lambda*s)-1)}# lambda=1
glogit.inv <- function(s, lambda){(1/lambda)*log(exp(s)+1)}

# density
dgad <-function(y,eta=1,phi=1,alpha=0.5,lambda=NULL, f0,g)
{
   g.prime <- Deriv::Deriv(g, "s")
   if(is.null(lambda)){
      coeff <- 2*alpha*(1-alpha)*g.prime(y)/phi
      den <-   ifelse(y<eta,(coeff*f0((1-alpha)*(g(eta)-g(y))/phi)),
                      (coeff*f0(alpha*(g(y)-g(eta))/phi)))
   } else{
      coeff <- 2*alpha*(1-alpha)*g.prime(y, lambda)/phi
      den <-   ifelse(y < eta, (coeff*f0((1-alpha)*(g(eta, lambda)-g(y, lambda))/phi)),
                      (coeff*f0(alpha*(g(y, lambda)-g(eta, lambda))/phi)))
   }
   return(den=den)
}
# CDF and Survival

pgad <- function(y, eta=1, phi=1, alpha=0.5, lambda = NULL,
                 F0, g, lower.tail = TRUE) {
   if(is.null(lambda)) {
      p <- ifelse(y < eta, (2*alpha*F0((1 - alpha)*(g(y) - g(eta))/phi)),
                  (2*alpha - 1 + 2*(1 - alpha)*F0(alpha*(g(y) - g(eta))/phi)))
   } else {
      p <- ifelse(y < eta, (2*alpha*F0((1 - alpha)*(g(y, lambda) - g(eta, lambda))/phi)),
                  (2*alpha - 1 + 2*(1 - alpha)*F0(alpha*(g(y, lambda) - g(eta, lambda))/phi)))
   }
   ifelse(test = lower.tail == TRUE, yes = return(p), no = return(1-p))
}
# quantile function
qgad <- function(tau, eta=1, phi=1, alpha=0.5, lambda = NULL,
                 F0, g, g.inv, QF = NULL){
   if (is.null(QF)){
      QF<-GoFKernel::inverse(F0, lower =0, upper = Inf)
   }
   if(is.null(lambda)){

      q<-ifelse(tau< alpha, g.inv(g(eta) + (phi/(1-alpha))*QF(tau/(2*alpha))),
                g.inv(g(eta) + (phi/alpha)*QF((1+tau-2*alpha)/(2*(1-alpha)))))
   } else {

      q<-ifelse(tau < alpha, g.inv((g(eta, lambda)+(phi/(1-alpha))*QF(tau/(2*alpha))), lambda),
                g.inv((g(eta, lambda)+(phi/alpha)*QF((1+tau-2*alpha)/(2*(1-alpha)))), lambda))
   }
   return(q=q)
}


# Random number generation
rgad<- function(n,eta=1,phi=1,alpha=0.5,lambda = NULL, F0,QF, g, g.inv){
   u <- runif(n, min = 0, max = 1)
   if (is.null(lambda)){
      r <- ifelse(u< alpha, g.inv(g(eta)+(phi/(1-alpha))*QF(u/(2*alpha))),
                  g.inv(g(eta)+(phi/alpha)*QF((1+u-2*alpha)/(2*(1-alpha)))))
   } else{

      r <- ifelse(u< alpha, g.inv((g(eta, lambda) + (phi/(1-alpha))*QF(u/(2*alpha))), lambda),
                  g.inv((g(eta, lambda)+(phi/alpha)*QF((1+u-2*alpha)/(2*(1-alpha)))), lambda))
   }

   return(r=r)
}
# GQBA for student's t family is defined separately.

# density for GQBA student's-t disribution
dgad_st <- function(y,eta = 1, phi = 1, alpha = 0.5,lambda = NULL, nu = 2, g)
{
   g.prime <- Deriv::Deriv(g, "s")
   if(is.null(lambda)){
      coeff <- 2*alpha*(1-alpha)*g.prime(y)/phi
      den <-   ifelse(y < eta, (coeff*f0_st(((1-alpha)*(g(eta)-g(y))/phi),nu)),
                      (coeff*f0_st((alpha*(g(y)-g(eta))/phi),nu)))
   } else{
      coeff <- 2*alpha*(1-alpha)*g.prime(y, lambda)/phi
      den <-   ifelse(y < eta, (coeff*f0_st(((1-alpha)*(g(eta, lambda)-g(y, lambda))/phi),nu)),
                      (coeff*f0_st((alpha*(g(y, lambda)-g(eta, lambda))/phi),nu)))
   }
   return(den=den)
}
# CDF and Survival
pgad_st <- function(y, eta = 1, phi = 1, alpha = 0.5, lambda = NULL, nu = 2,
                    g, lower.tail = TRUE) {
   if(is.null(lambda)) {
      p = ifelse(y < eta, (2*alpha*F0_st(((1 - alpha)*(g(y) - g(eta))/phi),nu)),
                 (2*alpha - 1 + 2*(1 - alpha)*F0_st((alpha*(g(y) - g(eta))/phi),nu)))
   } else {
      p <- ifelse(y < eta, (2*alpha*F0_st(((1 - alpha)*(g(y, lambda) - g(eta, lambda))/phi),nu)),
                  (2*alpha - 1 + 2*(1 - alpha)*F0_st((alpha*(g(y, lambda) - g(eta, lambda))/phi),nu)))
   }
   ifelse(test = lower.tail == TRUE, yes = return(p), no = return(1-p))
}
# quantile function
qgad_st <- function(tau, eta=1, phi=1, alpha=0.5, lambda = NULL, nu = 2,
                    g, g.inv ){
   if(is.null(lambda)){

      q<-ifelse(tau< alpha, g.inv(g(eta) + (phi/(1-alpha))*QF_st((tau/(2*alpha)),nu)),
                g.inv(g(eta) + (phi/alpha)*QF_st(((1+tau-2*alpha)/(2*(1-alpha))),nu)))
   } else {

      q<-ifelse(tau < alpha, g.inv((g(eta, lambda) + (phi/(1-alpha))*QF_st((tau/(2*alpha)),nu)), lambda),
                g.inv((g(eta, lambda) + (phi/alpha)*QF_st(((1+tau-2*alpha)/(2*(1-alpha))),nu)), lambda))
   }
   return(q=q)
}


# Random number generation from GQBA student's-t disribution
rgad_st <- function(n,eta=1,phi=1,alpha=0.5,lambda=NULL, nu = 2, g, g.inv){
   u <- runif(n, min = 0, max = 1)
   if (is.null(lambda)){
      r <- ifelse(u< alpha, g.inv(g(eta)+(phi/(1-alpha))*QF_st(u/(2*alpha),nu)),
                  g.inv(g(eta)+(phi/alpha)*QF_st((1+u-2*alpha)/(2*(1-alpha)),nu)))
   } else{

      r <- ifelse(u< alpha, g.inv((g(eta, lambda) + (phi/(1-alpha))*QF_st(u/(2*alpha),nu)), lambda),
                  g.inv((g(eta, lambda)+(phi/alpha)*QF_st((1+u-2*alpha)/(2*(1-alpha)),nu)), lambda))
   }

   return(r=r)
}

# log-likelihood function for  GQBA families: AL, AN, ALa
# Implementing our estimation algorithm used in the real data analysis
nll_GQBA<- function(theta, thetahat = NULL, y, d, link = c("log", "logit"), f0, F0){
   y <- as.matrix(y) ; d <- as.matrix(d)
   if(link == "log" & !is.null(thetahat)) {
      alpha <- theta[1] # Step 2
      den <- dgad(y, eta = thetahat[1], phi = thetahat[2], alpha = alpha, f0 = f0, g=glog)
      sur <- pgad(y, eta = thetahat[1], phi = thetahat[2], alpha = alpha, F0 = F0, g=glog, lower.tail = FALSE)
   } else if(link == "log" & is.null(thetahat)) {
      eta <- theta[1]; phi  <-  theta[2] ; alpha <- theta[3] # step 3
      den <- dgad(y, eta = eta, phi = phi, alpha = alpha, f0 = f0, g = glog)
      sur <- pgad(y, eta = eta, phi = phi, alpha = alpha, F0 = F0, g = glog, lower.tail = FALSE)
   } else if(link == "logit" & !is.null(thetahat)) {
      lambda <- theta[1] # step 4
      den <- dgad(y, eta = thetahat[1], phi = thetahat[2], alpha = thetahat[3],
                  lambda = lambda, f0 = f0, g = glogit)
      sur <- pgad(y, eta = thetahat[1], phi = thetahat[2], alpha = thetahat[3],
                  lambda = lambda, F0 = F0, g = glogit, lower.tail = FALSE)
   } else {
      # final step 5
      eta <-  theta[1]; phi <- theta[2]; alpha <-  theta[3]; lambda <-  theta[4]
      den <- dgad(y, eta = eta, phi = phi, alpha = alpha,
                  lambda = lambda, f0 = f0, g = glogit)
      sur <- pgad(y, eta = eta, phi = phi, alpha = alpha,
                  lambda = lambda, F0 = F0, g = glogit, lower.tail = FALSE)
   }
   LL <- suppressWarnings((d*log(den) + (1 - d)*log(sur)))
   return(-sum(LL[!is.infinite(LL)]))
}

# likelihood function for GQBA student's-t family: ASt
nll_student <- function(theta, thetahat = NULL, nu = NULL, y, d,  link = c("log", "logit")){
   y <- as.matrix(y)
   d <- as.matrix(d)

   if(link == "log" & !is.null(thetahat) & !is.null(nu)) {
      alpha <- theta[1]     # estimate alpha based on step 2
      den <- dgad_st(y, eta = thetahat[1], phi = thetahat[2], alpha = alpha, nu = nu, g = glog)
      sur <- pgad_st(y, eta = thetahat[1], phi = thetahat[2], alpha = alpha, nu = nu, g = glog, lower.tail = FALSE)

   }  else if(link == "log" & !is.null(thetahat) & is.null(nu)) {
      nu <- theta[1] # estimate v based on step 2
      den <- dgad_st(y, eta = thetahat[1], phi = thetahat[2], alpha = thetahat[3], nu = nu, g = glog)
      sur <- pgad_st(y, eta = thetahat[1], phi = thetahat[2], alpha = thetahat[3], nu = nu, g = glog, lower.tail = FALSE)

   }  else if(link == "log" & is.null(thetahat) & !is.null(nu)) {
      # final estimate for log link with known v
      eta <- theta[1]; phi  <-  theta[2] ; alpha <- theta[3]
      den <- dgad_st(y, eta = eta, phi = phi, alpha = alpha, nu = nu, g = glog)
      sur <- pgad_st(y, eta = eta, phi = phi, alpha = alpha, nu = nu, g = glog, lower.tail = FALSE)

   } else if(link == "log" & is.null(thetahat) & is.null(nu)) {
      # final estimate for log link with unknown v
      eta <- theta[1]; phi  <-  theta[2]; alpha <- theta[3]; nu = theta[4]
      den <- dgad_st(y, eta = eta, phi = phi, alpha = alpha, nu = nu, g = glog)
      sur <- pgad_st(y, eta = eta, phi = phi, alpha = alpha, nu = nu, g = glog, lower.tail = FALSE)

   }   else if(link == "logit" & !is.null(thetahat) & !is.null(nu)) {
      lambda <- theta[1] # estimate lambda based on step 4
      den <- dgad_st(y, eta = thetahat[1], phi = thetahat[2], alpha = thetahat[3],
                     lambda = lambda, nu = nu, g = glogit)
      sur <- pgad_st(y, eta = thetahat[1], phi = thetahat[2], alpha = thetahat[3],
                     lambda = lambda, nu = nu, g = glogit, lower.tail = FALSE)

   }  else if(link == "logit" & is.null(thetahat) & !is.null(nu)) {
      # final estimate for logit type link with known v
      eta <-  theta[1]; phi <- theta[2]; alpha <-  theta[3]; lambda <-  theta[4]
      den <- dgad_st(y, eta = eta, phi = phi, alpha = alpha,
                     lambda = lambda, nu = nu, g = glogit)
      sur <- pgad_st(y, eta = eta, phi = phi, alpha = alpha,
                     lambda = lambda, nu = nu, g = glogit, lower.tail = FALSE)

   } else {       # final estimate for logit type link with unknown v
      eta <-  theta[1]; phi <- theta[2]; alpha <-  theta[3]; lambda <-  theta[4]; nu <- theta[5]
      den <- dgad_st(y, eta = eta, phi = phi, alpha = alpha,
                     lambda = lambda, nu = nu, g = glogit)
      sur <- pgad_st(y, eta = eta, phi = phi, alpha = alpha,
                     lambda = lambda, nu = nu, g = glogit, lower.tail = FALSE)
   }
   LL <- suppressWarnings((d*log(den) + (1 - d)*log(sur)))
   return(-sum(LL[!is.infinite(LL)]))
}
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Example 1: Times to weaning of breast-fed data
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
#gbc <- data.frame("time"=gbcsIDM$rectime/365.25, "status"=gbcsIDM$censrec)
#TAM <- subset(gbcsIDM, gbcsIDM$hormone==2) # hormone therapy treatment group
#NOTAM <- subset(gbcsIDM, gbcsIDM$hormone==1)# untreated group
#data_gbc <-  gbc
#data_tam <-  data.frame("time"=TAM$rectime/365.25, "status"=TAM$censrec)
#data_notam <-  data.frame("time"=NOTAM$rectime/365.25, "status"=NOTAM$censrec)

data_feed <-  data.frame("time"=bfeed$duration, "status"=bfeed$delta)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Log-logistic classical distribution
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
library(FAdist) # for initail value
llogistic_classic <- function(y, d, ...){
   y <- as.matrix(y)
   d <- as.matrix(d)
   if (nrow(y)[1] != nrow(d)) {
      stop(paste("y and d should have equal  dimension."))
   }
   nlogL <- function(theta) {
      eta <- theta[1]
      mu <- log(eta)
      phi <- theta[2]
      den <- dllog(y,  scale = mu, shape = phi)
      sur <- pllog(y,  scale = mu, shape = phi, lower.tail = FALSE)
      LL <- (d*log(den) + (1 - d)*log(sur))
      return(-sum(LL[!is.infinite(LL)]))
   }
   # starting values
   data<- as.data.frame(cbind(y, d))
   y.sort <- sort(data[,1][data[,2]==1])
   n.c <- length(y.sort)
   W <- numeric(0)
   for (j in 1:n.c) {
      W[j]<-y.sort[j]*((j-1)/(n.c-1))
   }
   W0 <- mean(y.sort)
   W1 <- sum(W)/n.c
   phi.pwmom <- (2*W1-W0)/W0
   eta.pwmom <- (1/phi.pwmom)*(W0*sin(pi*phi.pwmom)/pi)

   theta.start <- c(eta.pwmom, phi.pwmom)

   opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD",
                "xtol_rel" = 1.0e-10,  "print_level"=0,
                "maxeval" = 1000)
   fit_nm <- nloptr(x0 = theta.start, eval_f =  nlogL,
                    lb = c(1.0e-10, 1.0e-10),
                    ub = c(Inf,  Inf), opts = opts)
   estimates = c(fit_nm$solution[1], fit_nm$solution[2])
   nll = fit_nm$objective
   AIC = 2*2 + 2*nll
   BIC = 2*log(nrow(y)) + 2*nll
   return(list(MLE = estimates, nll = nll, AIC = AIC, BIC = BIC))
}

fit_llogistic <- llogistic_classic(y = data_feed$time, d = data_feed$status)
par_llogistic <-  round(fit_llogistic$MLE, 3)
nll_llogistic <-  fit_llogistic$nll
AIC_llogistic <-  fit_llogistic$AIC
BIC_llogistic <-  fit_llogistic$BIC

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Log-normal classical distribution
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
library(EnvStats) # for initial value

lnormal_classic <- function(y, d, ...){
   y <- as.matrix(y)
   d <- as.matrix(d)
   nlogL <- function(theta) {
      eta <- theta[1]
      phi <- theta[2]
      mu <- log(eta)
      den <- dlnorm(y, meanlog = mu,  sdlog = phi)
      sur <- plnorm(y, meanlog = mu,  sdlog = phi, lower.tail = FALSE)
      LL <- (d*log(den) + (1 - d)*log(sur))
      return(-sum(LL[!is.infinite(LL)]))
   }
   # starting values
   data <- as.data.frame(cbind(y, d))
   x = as.vector(data[,1])
   theta.start <- as.numeric(elnorm(x= x, method =  "mle/mme")$parameters)

   opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD",
                "xtol_rel" = 1.0e-10,
                "maxeval" = 1000)
   fit_nm <- nloptr(x0 = theta.start , eval_f =  nlogL,
                    lb = c(1.0e-10, 1.0e-10),
                    ub = c(Inf,  Inf), opts = opts)
   estimates = c(fit_nm$solution[1], fit_nm$solution[2])
   nll = fit_nm$objective
   AIC = 2*2 + 2*nll
   BIC = 2*log(nrow(y)) + 2*nll
   return(list(MLE = estimates, nll = nll, AIC = AIC, BIC = BIC))
}

fit_lnormal <- lnormal_classic(y = data_feed$time, d = data_feed$status)
par_lnormal  <- fit_lnormal$MLE
nll_lnormal  <- fit_lnormal$nll
AIC_lnormal  <- fit_lnormal$AIC
BIC_lnormal  <- fit_lnormal$BIC

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Log-laplace classical distribution
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
library(LaplacesDemon)# for log-laplace disrtibtuon
library(ExtDist) # for laplace distribution

llaplace_classic <- function(y, d, ...){
   y <- as.matrix(y)
   d <- as.matrix(d)

   nlogL <- function(theta) {
      eta <- theta[1]
      phi <- theta[2]
      den <- dllaplace(y, location = eta,  scale  = phi)
      sur <- 1 - pllaplace(y, location = eta,  scale = phi)
      LL <- (d*log(den) + (1 - d)*log(sur))
      return(-sum(LL[!is.infinite(LL)]))
   }

   # starting values
   data <- as.data.frame(cbind(y, d))
   x = as.vector(data[,1])
   fit.lap <- eLaplace(x, method = c("analytic.MLE"))

   theta.start <- c(fit.lap$mu, fit.lap$b)

   opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD",
                "xtol_rel" = 1.0e-8,
                "maxeval" = 10000)
   fit_nm <- nloptr(x0 = theta.start , eval_f =  nlogL,
                    lb = c(1.0e-8, 1.0e-8),
                    ub = c(Inf,  Inf), opts = opts)
   estimates = c(exp(fit_nm$solution[1])/2, fit_nm$solution[2])
   nll = fit_nm$objective
   AIC = 2*2 + 2*nll
   BIC = 2*log(nrow(y)) + 2*nll
   return(list(MLE = estimates, nll = nll, AIC = AIC, BIC = BIC))
}
# check VGAM package
fit_llaplace <- llaplace_classic(y = data_feed$time, d = data_feed$status)
par_llaplace <- fit_llaplace$MLE
nll_llaplace <- fit_llaplace$nll
AIC_llaplace <- fit_llaplace$AIC
BIC_llaplace <- fit_llaplace$BIC
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Weibull distribution
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
library(ForestFit)

weibull <- function(y, d,...){
   y <- as.matrix(y)
   d <- as.matrix(d)
   nlogL <- function(theta) {
      eta <- theta[1]
      phi <- theta[2]
      den <- dweibull(y, shape = phi, scale = eta)
      sur <- pweibull(y, shape = phi, scale = eta, lower.tail = FALSE)
      LL <- (d*log(den) + (1 - d)*log(sur))
      return(-sum(LL[!is.infinite(LL)]))
   }
   data <- as.data.frame(cbind(y, d))
   theta.start <- fitWeibull(data = data[,1][data[,2] == 1],
                             location =  FALSE, method =  "moment",
                             starts=c(1,1,0))$estimate
   theta.start <- c(theta.start[2], theta.start[1])
   opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD",
                "xtol_rel" = 1.0e-10,
                "maxeval" = 1000)
   fit_nm <- nloptr(x0 = as.numeric(theta.start) , eval_f =  nlogL,
                    lb = c(1.0e-10, 1.0e-10),
                    ub = c(Inf,  Inf), opts = opts)
   estimates = fit_nm$solution
   nll = fit_nm$objective
   AIC = 2*2 + 2*nll
   BIC = 2*log(nrow(y)) + 2*nll
   return(list(MLE = estimates, nll = nll, AIC = AIC, BIC = BIC))
}

fit_weibull <-  weibull(y=data_feed$time, d=data_feed$status)
par_weibull <- fit_weibull$MLE
nll_weibull <- fit_weibull$nll
AIC_weibull <- fit_weibull$AIC
BIC_weibull <- fit_weibull$BIC

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA normal with log link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# first estimate alpha
library(QBAsyDist)# initial value
int_alpha <- mleAND(log(data_feed$time[data_feed$status==1]))$alpha.MLE
alphahat <- optim(par = int_alpha, fn  = nll_GQBA, method = "BFGS",
                  thetahat = par_lnormal, link = "log", f0=f0_AN, F0 =F0_AN,
                  y = data_feed$time, d = data_feed$status)$par
#
theta.start <- c(par_lnormal, alphahat)
opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"=1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level"=1)
fit_AN_log<-  nloptr(x0 = theta.start, eval_f   = nll_GQBA, opts = opts,
                     thetahat = NULL, y = data_feed$time, d = data_feed$status, link ="log",
                     f0 = f0_AN, F0 = F0_AN,
                     lb = rep(1.0e-10, length(theta.start)),
                     ub = c(Inf, Inf, 0.9999))
par_AN_log <- fit_AN_log$solution
nll_AN_log <- fit_AN_log$objective
AIC_AN_log <- 2*3 + 2*nll_AN_log
BIC_AN_log <- 2*log(nrow(data_feed)) + 2*nll_AN_log

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA normal with logit-type link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# first estimate lambda
opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel"=1.0e-4,
             "xtol_abs" = 1.0e-4,"maxeval" = 10000, "print_level"=1)

lambdahat <- nloptr(x0 = 0.0102, eval_f  = nll_GQBA, opts = opts,
                    thetahat = par_AN_log, link = "logit", f0 = f0_AN, F0 = F0_AN,
                    y = data_feed$time, d = data_feed$status,
                    lb = 0, ub = Inf)$solution
#
theta.start <- c(par_lnormal, alphahat, lambdahat)

opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"= 1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level"=1) #  NLOPT_LN_SBPLX takes morethan 10,000 iteration
fit_AN_logit<-  nloptr(x0 = theta.start, eval_f   = nll_GQBA, opts = opts,
                       thetahat = NULL, y = data_feed$time, d = data_feed$status, link ="logit",
                       f0 = f0_AN, F0 = F0_AN,
                       lb = rep(1.0e-10, length(theta.start)),
                       ub = c(Inf, Inf, 0.9999, Inf))
par_AN_logit <- fit_AN_logit$solution
nll_AN_logit <- fit_AN_logit$objective
AIC_AN_logit <- 2*4 + 2*nll_AN_logit
BIC_AN_logit <- 2*log(nrow(data_feed)) + 2*nll_AN_logit

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA logistic with log link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# first estimate alpha
library(QBAsyDist)# initial value
int_alpha <- mleALoD(log(data_feed$time[data_feed$status==1]))$alpha.MLE
alphahat <- optim(par = int_alpha, fn  = nll_GQBA, method = "BFGS",
                  thetahat = par_llogistic, link = "log", f0 = f0_AL, F0 =F0_AL,
                  y = data_feed$time, d = data_feed$status)$par
#
theta.start <- c(par_llogistic, alphahat)
opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"=1.0e-4,
             "xtol_abs" = 1.0e-4,
             "maxeval" = 10000, "print_level"=1)
fit_AL_log<-  nloptr(x0 = theta.start, eval_f   = nll_GQBA, opts = opts,
                     thetahat = NULL, y = data_feed$time, d = data_feed$status, link ="log",
                     f0 = f0_AL, F0 = F0_AL,
                     lb = rep(1.0e-10, length(theta.start)),
                     ub = c(Inf, Inf, 0.9999))
par_AL_log <- fit_AL_log$solution
nll_AL_log <- fit_AL_log$objective
AIC_AL_log <- 2*3 + 2*nll_AL_log
BIC_AL_log <- 2*log(nrow(data_feed)) + 2*nll_AL_log

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA logistic with logit-type link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# first estimate lambda
opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel"=1.0e-4,
             "xtol_abs" = 1.0e-4,"maxeval" = 10000, "print_level"=1)

lambdahat <- nloptr(x0 = 0.015, eval_f  = nll_GQBA, opts = opts,
                    thetahat = par_AL_log, link = "logit", f0 = f0_AL, F0 = F0_AL,
                    y = data_feed$time, d = data_feed$status,
                    lb = 0, ub = Inf)$solution
#theta.start <- c(par_AN_log, lambdahat) # takes longer time
theta.start <- par_AN_logit
opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel"= 1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level"=1) #  NLOPT_LN_SBPLX takes morethan 10,000 iteration
fit_AL_logit <-  nloptr(x0 = theta.start, eval_f   = nll_GQBA, opts = opts,
                        thetahat = NULL, y = data_feed$time, d = data_feed$status, link ="logit",
                        f0 = f0_AL, F0 = F0_AL,
                        lb = rep(0.00001, length(theta.start)),
                        ub = c(Inf, Inf, 0.9999, Inf))
par_AL_logit <- fit_AL_logit$solution
nll_AL_logit <- fit_AL_logit$objective
AIC_AL_logit <- 2*4 + 2*nll_AL_logit
BIC_AL_logit <- 2*log(nrow(data_feed)) + 2*nll_AL_logit

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA Laplace with log link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# first estimate alpha
int_alpha <- mleALaD(log(data_feed$time[data_feed$status==1]))$alpha.ALaD
alphahat <- optim(par = int_alpha, fn  = nll_GQBA, method = "BFGS",
                  thetahat = par_llogistic, link = "log", f0 = f0_ALa, F0 = F0_ALa,
                  y = data_feed$time, d = data_feed$status)$par
#
theta.start <- c(par_llaplace, alphahat)
#theta.start <- par_AN_log
opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"=1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level"=1)
fit_ALa_log <-  nloptr(x0 = theta.start, eval_f   = nll_GQBA, opts = opts,
                       thetahat = NULL, y = data_feed$time, d = data_feed$status, link ="log",
                       f0 = f0_ALa, F0 = F0_ALa,
                       lb = rep(1.0e-10, length(theta.start)),
                       ub = c(Inf, Inf, 0.9999))
par_ALa_log <- fit_ALa_log$solution
nll_ALa_log <- fit_ALa_log$objective
AIC_ALa_log <- 2*3 + 2*nll_ALa_log
BIC_ALa_log <- 2*log(nrow(data_feed)) + 2*nll_ALa_log

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA Laplace with logit-type link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# first estimate lambda
opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel"=1.0e-8,
             "xtol_abs" = 1.0e-8,"maxeval" = 10000, "print_level"=1)

lambdahat <- nloptr(x0 = 0.015, eval_f  = nll_GQBA, opts = opts,
                    thetahat = par_ALa_log, link = "logit", f0 = f0_ALa, F0 = F0_ALa,
                    y = data_feed$time, d = data_feed$status,
                    lb = 0, ub = Inf)$solution
#theta.start <- c(par_ALa_log, 0.012)
#theta.start <- par_AN_logit
theta.start <- c(par_ALa_log, lambdahat)
opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel"= 1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level"=1) #  NLOPT_LN_SBPLX takes morethan 10,000 iteration
fit_ALa_logit <-  nloptr(x0 = theta.start, eval_f   = nll_GQBA, opts = opts,
                         thetahat = NULL, y = data_feed$time, d = data_feed$status, link ="logit",
                         f0 = f0_ALa, F0 = F0_ALa,
                         lb = rep(1.0e-8, length(theta.start)),
                         ub = c(Inf, Inf, 0.9999, Inf))
par_ALa_logit <- fit_ALa_logit$solution
nll_ALa_logit <- fit_ALa_logit$objective
AIC_ALa_logit <- 2*4 + 2*nll_ALa_logit
BIC_ALa_logit <- 2*log(nrow(data_feed)) + 2*nll_ALa_logit

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA student's t with log link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
mom <- momATD(log(data_feed$time[data_feed$status==1]))
vhat <- mom$nu.MoM
alphamom <- mom$alpha.MoM
opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"=1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level" = 1)
alphahat_st <- optim(par = alphamom, fn  = nll_student, method = "BFGS",
                     thetahat = par_lnormal, nu = vhat, link = "log",
                     y = data_feed$time, d = data_feed$status)$par

thetahat_st <- c(par_lnormal, alphahat_st)
vhat_st <- optim(par = vhat, fn = nll_student, method = "BFGS",
                 thetahat = thetahat_st, nu = NULL, link = "log",
                 y = data_feed$time, d = data_feed$status)$par
# known nu
par_st_log <-  nloptr(x0 = thetahat_st, eval_f   = nll_student, opts = opts,
                      thetahat = NULL, y = data_feed$time, d = data_feed$status,
                      nu = vhat, link ="log",
                      lb = rep(1.0e-10, length(thetahat_st)),
                      ub = c(Inf, Inf, 0.9999))$solution
# unknown nu
#thetahat_st_nu <- c(par_st_log, vhat)
thetahat_st_nu <- c(par_st_log, vhat)
opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"=1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level" = 1)
fit_Ast_log <-  nloptr(x0 = thetahat_st_nu, eval_f  = nll_student, opts = opts,
                       thetahat = NULL, y = data_feed$time, d = data_feed$status,
                       nu = NULL,  link ="log",
                       lb = c(1.0e-10, 1.0e-10, 1.0e-10, 1),
                       ub = c(Inf, Inf, 0.9999, Inf))
par_Ast_log <- fit_Ast_log$solution
nll_Ast_log <- fit_Ast_log$objective
AIC_Ast_log <- 2*4 + 2*nll_Ast_log
BIC_Ast_log <- 2*log(nrow(data_feed)) + 2*nll_Ast_log

#thetahat <- par_st_log_nu$solution[1:3]
lambdahat_st <- optim(par = 0.015, fn  = nll_student,
                      thetahat = par_Ast_log, y = data_feed$time, d = data_feed$status,
                      nu = vhat,   link = "logit")$par
theta.start_st <- c(par_AN_logit, par_Ast_log[4])
#t <- c(19.865, 0.486, 0.697, 0.03573719, 80)
fit_Ast_logit <-  nloptr(x0 = theta.start_st, eval_f   = nll_student, opts = opts,
                         thetahat = NULL, y = data_feed$time, d = data_feed$status,
                         nu = NULL,  link ="logit",
                         lb = c(rep(1.0e-10, 4), 1),
                         ub = c(Inf, Inf, 0.99999, Inf, Inf))
par_Ast_logit <- fit_Ast_logit$solution
nll_Ast_logit <- fit_Ast_logit$objective
AIC_Ast_logit <- 2*5 + 2*nll_Ast_logit
BIC_Ast_logit <- 2*log(nrow(data_feed)) + 2*nll_Ast_logit

#
nll <- rbind(nll_llogistic, nll_lnormal, nll_llaplace, nll_weibull,
             nll_AL_log, nll_AN_log, nll_ALa_log, nll_Ast_log,
             nll_AL_logit, nll_AN_logit, nll_ALa_logit, nll_Ast_logit)
AIC <- rbind(AIC_llogistic, AIC_lnormal, AIC_llaplace, AIC_weibull,
             AIC_AL_log, AIC_AN_log, AIC_ALa_log, AIC_Ast_log,
             AIC_AL_logit, AIC_AN_logit, AIC_ALa_logit, AIC_Ast_logit)
BIC <- rbind(BIC_llogistic, BIC_lnormal, BIC_llaplace, BIC_weibull,
             BIC_AL_log, BIC_AN_log, BIC_ALa_log, BIC_Ast_log,
             BIC_AL_logit, BIC_AN_logit, BIC_ALa_logit, BIC_Ast_logit)
Estimate <- rbind(par_llogistic, par_lnormal, par_llaplace, par_weibull,
                  par_AL_log, par_AN_log, par_ALa_log, par_Ast_log,
                  par_AL_logit, par_AN_logit, par_ALa_logit, par_Ast_logit)
breast_par_nll_AIC_BIC <- cbind(Estimate, nll, AIC, BIC)
model <- c("llogistic", "lnormal", "llaplace", "weibull",
           "AL_log", "AN_log", "ALa_log", "Ast_log",
           "AL_logit", "AN_logit", "ALa_logit", "Ast_logit")
colN <- c("etahat","phihat", "alphahat", "lambdahat", "nuhat","nll", "AIC","BIC")
colnames(breast_par_nll_AIC_BIC) <- paste("",colN)
rownames(breast_par_nll_AIC_BIC) <- paste("",model)
# exporting the summarized results
#write.csv(breast_par_nll_AIC_BIC, "breast_par_nll_AIC_BIC_updated.csv")

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Example 2: German breast cancer data , all patients
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
data_gbc <- data.frame("time"=gbcsIDM$rectime/365.25, "status"=gbcsIDM$censrec)
#TAM <- subset(gbcsIDM, gbcsIDM$hormone==2) # hormone therapy treatment group
#NOTAM <- subset(gbcsIDM, gbcsIDM$hormone==1)# untreated group
#data_gbc <-  gbc
#data_tam <-  data.frame("time"=TAM$rectime/365.25, "status"=TAM$censrec)
#data_notam <-  data.frame("time"=NOTAM$rectime/365.25, "status"=NOTAM$censrec)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# Log-logistic classical distribution
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
fit_llogistic <- llogistic_classic(y = data_gbc$time, d = data_gbc$status)
par_llogistic <-  round(fit_llogistic$MLE, 3)
nll_llogistic <-  fit_llogistic$nll
AIC_llogistic <-  fit_llogistic$AIC
BIC_llogistic <-  fit_llogistic$BIC
# Log-normal classical distribution
fit_lnormal <- lnormal_classic(y = data_gbc$time, d = data_gbc$status)
par_lnormal  <- fit_lnormal$MLE
nll_lnormal  <- fit_lnormal$nll
AIC_lnormal  <- fit_lnormal$AIC
BIC_lnormal  <- fit_lnormal$BIC
# Log-laplace classical distribution
fit_llaplace <- llaplace_classic(y = data_gbc$time, d = data_gbc$status)
par_llaplace <- fit_llaplace$MLE
nll_llaplace <- fit_llaplace$nll
AIC_llaplace <- fit_llaplace$AIC
BIC_llaplace <- fit_llaplace$BIC

# Weibull distribution
weibull <- function(y, d,...){
   y <- as.matrix(y)
   d <- as.matrix(d)
   nlogL <- function(theta) {
      eta <- theta[1]
      phi <- theta[2]
      den <- dweibull(y, shape = phi, scale = eta)
      sur <- pweibull(y, shape = phi, scale = eta, lower.tail = FALSE)
      LL <- (d*log(den) + (1 - d)*log(sur))
      return(-sum(LL[!is.infinite(LL)]))
   }
   data <- as.data.frame(cbind(y, d))
   theta.start <- fitWeibull(data = data[,1][data[,2] == 1],
                             location =  FALSE, method =  "moment",
                             starts=c(1,1,0))$estimate
   theta.start <- c(theta.start[2], theta.start[1])
   opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD",
                "xtol_rel" = 1.0e-10,
                "maxeval" = 1000)
   fit_nm <- nloptr(x0 = as.numeric(theta.start) , eval_f =  nlogL,
                    lb = c(1.0e-10, 1.0e-10),
                    ub = c(Inf,  Inf), opts = opts)
   estimates = fit_nm$solution
   nll = fit_nm$objective
   AIC = 2*2 + 2*nll
   BIC = 2*log(nrow(y)) + 2*nll
   return(list(MLE = estimates, nll = nll, AIC = AIC, BIC = BIC))
}

fit_weibull <-  weibull(y=data_gbc$time, d=data_gbc$status)
par_weibull <- fit_weibull$MLE
nll_weibull <- fit_weibull$nll
AIC_weibull <- fit_weibull$AIC
BIC_weibull <- fit_weibull$BIC

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA normal with log link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# first estimate alpha
int_alpha <- mleAND(log(data_gbc$time[data_gbc$status==1]))$alpha.MLE
alphahat <- optim(par = int_alpha, fn  = nll_GQBA, method = "BFGS",
                  thetahat = par_lnormal, link = "log", f0=f0_AN, F0 =F0_AN,
                  y = data_gbc$time, d = data_gbc$status)$par
#
theta.start <- c(par_lnormal, alphahat)
opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"=1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level"=1)
fit_AN_log<-  nloptr(x0 = theta.start, eval_f   = nll_GQBA, opts = opts,
                     thetahat = NULL, y = data_gbc$time, d = data_gbc$status, link ="log",
                     f0 = f0_AN, F0 = F0_AN,
                     lb = rep(1.0e-10, length(theta.start)),
                     ub = c(Inf, Inf, 0.9999))
par_AN_log <- fit_AN_log$solution
nll_AN_log <- fit_AN_log$objective
AIC_AN_log <- 2*3 + 2*nll_AN_log
BIC_AN_log <- 2*log(nrow(data_gbc)) + 2*nll_AN_log

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA normal with logit-type link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# first estimate lambda
opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel"=1.0e-4,
             "xtol_abs" = 1.0e-4,"maxeval" = 10000, "print_level"=1)

lambdahat <- nloptr(x0 = 0.0102, eval_f  = nll_GQBA, opts = opts,
                    thetahat = par_AN_log, link = "logit", f0 = f0_AN, F0 = F0_AN,
                    y = data_gbc$time, d = data_gbc$status,
                    lb = 0, ub = Inf)$solution
#
theta.start <- c(par_lnormal, alphahat, lambdahat)

opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"= 1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level"=1) #  NLOPT_LN_SBPLX takes morethan 10,000 iteration
fit_AN_logit<-  nloptr(x0 = theta.start, eval_f   = nll_GQBA, opts = opts,
                       thetahat = NULL, y = data_gbc$time, d = data_gbc$status, link ="logit",
                       f0 = f0_AN, F0 = F0_AN,
                       lb = rep(1.0e-10, length(theta.start)),
                       ub = c(Inf, Inf, 0.9999, Inf))
par_AN_logit <- fit_AN_logit$solution
nll_AN_logit <- fit_AN_logit$objective
AIC_AN_logit <- 2*4 + 2*nll_AN_logit
BIC_AN_logit <- 2*log(nrow(data_gbc)) + 2*nll_AN_logit

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA logistic with log link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

# first estimate alpha
int_alpha <- mleALoD(log(data_gbc$time[data_gbc$status==1]))$alpha.MLE
alphahat <- optim(par = int_alpha, fn  = nll_GQBA, method = "BFGS",
                  thetahat = par_llogistic, link = "log", f0 = f0_AL, F0 =F0_AL,
                  y = data_gbc$time, d = data_gbc$status)$par
#
theta.start <- c(par_llogistic, alphahat)
opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"=1.0e-4,
             "xtol_abs" = 1.0e-4,
             "maxeval" = 10000, "print_level"=1)
fit_AL_log<-  nloptr(x0 = theta.start, eval_f   = nll_GQBA, opts = opts,
                     thetahat = NULL, y = data_gbc$time, d = data_gbc$status, link ="log",
                     f0 = f0_AL, F0 = F0_AL,
                     lb = rep(1.0e-10, length(theta.start)),
                     ub = c(Inf, Inf, 0.9999))
par_AL_log <- fit_AL_log$solution
nll_AL_log <- fit_AL_log$objective
AIC_AL_log <- 2*3 + 2*nll_AL_log
BIC_AL_log <- 2*log(nrow(data_gbc)) + 2*nll_AL_log

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA logistic with logit-type link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# first estimate lambda
opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel"=1.0e-4,
             "xtol_abs" = 1.0e-4,"maxeval" = 10000, "print_level"=1)

lambdahat <- nloptr(x0 = 0.015, eval_f  = nll_GQBA, opts = opts,
                    thetahat = par_AL_log, link = "logit", f0 = f0_AL, F0 = F0_AL,
                    y = data_gbc$time, d = data_gbc$status,
                    lb = 0, ub = Inf)$solution
#theta.start <- c(par_AN_log, lambdahat) # takes longer time
theta.start <- par_AN_logit
opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel"= 1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level"=1) #  NLOPT_LN_SBPLX takes morethan 10,000 iteration
fit_AL_logit <-  nloptr(x0 = theta.start, eval_f   = nll_GQBA, opts = opts,
                        thetahat = NULL, y = data_gbc$time, d = data_gbc$status, link ="logit",
                        f0 = f0_AL, F0 = F0_AL,
                        lb = rep(0.00001, length(theta.start)),
                        ub = c(Inf, Inf, 0.9999, Inf))
par_AL_logit <- fit_AL_logit$solution
nll_AL_logit <- fit_AL_logit$objective
AIC_AL_logit <- 2*4 + 2*nll_AL_logit
BIC_AL_logit <- 2*log(nrow(data_gbc)) + 2*nll_AL_logit

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA Laplace with log link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# first estimate alpha
int_alpha <- mleALaD(log(data_gbc$time[data_gbc$status==1]))$alpha.ALaD
alphahat <- optim(par = int_alpha, fn  = nll_GQBA, method = "BFGS",
                  thetahat = par_llogistic, link = "log", f0 = f0_ALa, F0 = F0_ALa,
                  y = data_gbc$time, d = data_gbc$status)$par
#
theta.start <- c(par_llaplace, alphahat)
#theta.start <- par_AN_log
opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"=1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level"=1)
fit_ALa_log <-  nloptr(x0 = theta.start, eval_f   = nll_GQBA, opts = opts,
                       thetahat = NULL, y = data_gbc$time, d = data_gbc$status, link ="log",
                       f0 = f0_ALa, F0 = F0_ALa,
                       lb = rep(1.0e-10, length(theta.start)),
                       ub = c(Inf, Inf, 0.9999))
par_ALa_log <- fit_ALa_log$solution
nll_ALa_log <- fit_ALa_log$objective
AIC_ALa_log <- 2*3 + 2*nll_ALa_log
BIC_ALa_log <- 2*log(nrow(data_gbc)) + 2*nll_ALa_log

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA Laplace with logit-type link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# first estimate lambda
opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel"=1.0e-8,
             "xtol_abs" = 1.0e-8,"maxeval" = 10000, "print_level"=1)

lambdahat <- nloptr(x0 = 0.015, eval_f  = nll_GQBA, opts = opts,
                    thetahat = par_ALa_log, link = "logit", f0 = f0_ALa, F0 = F0_ALa,
                    y = data_gbc$time, d = data_gbc$status,
                    lb = 0, ub = Inf)$solution
#theta.start <- c(par_ALa_log, 0.012)
#theta.start <- par_AN_logit
theta.start <- c(par_ALa_log, lambdahat)
opts <- list("algorithm" = "NLOPT_LN_SBPLX", "xtol_rel"= 1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level"=1) #  NLOPT_LN_SBPLX takes morethan 10,000 iteration
fit_ALa_logit <-  nloptr(x0 = theta.start, eval_f   = nll_GQBA, opts = opts,
                         thetahat = NULL, y = data_gbc$time, d = data_gbc$status, link ="logit",
                         f0 = f0_ALa, F0 = F0_ALa,
                         lb = rep(1.0e-8, length(theta.start)),
                         ub = c(Inf, Inf, 0.9999, Inf))
par_ALa_logit <- fit_ALa_logit$solution
nll_ALa_logit <- fit_ALa_logit$objective
AIC_ALa_logit <- 2*4 + 2*nll_ALa_logit
BIC_ALa_logit <- 2*log(nrow(data_gbc)) + 2*nll_ALa_logit

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
# GQBA student's t with log link function
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#
mom <- momATD(log(data_gbc$time[data_gbc$status==1]))
vhat <- mom$nu.MoM
alphamom <- mom$alpha.MoM
opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"=1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level" = 1)
alphahat_st <- optim(par = alphamom, fn  = nll_student, method = "BFGS",
                     thetahat = par_lnormal, nu = vhat, link = "log",
                     y = data_gbc$time, d = data_gbc$status)$par

thetahat_st <- c(par_lnormal, alphahat_st)
vhat_st <- optim(par = vhat, fn = nll_student, method = "BFGS",
                 thetahat = thetahat_st, nu = NULL, link = "log",
                 y = data_gbc$time, d = data_gbc$status)$par
# known nu
par_st_log <-  nloptr(x0 = thetahat_st, eval_f   = nll_student, opts = opts,
                      thetahat = NULL, y = data_gbc$time, d = data_gbc$status,
                      nu = vhat, link ="log",
                      lb = rep(1.0e-10, length(thetahat_st)),
                      ub = c(Inf, Inf, 0.9999))$solution
# unknown nu
#thetahat_st_nu <- c(par_st_log, vhat)
thetahat_st_nu <- c(par_st_log, vhat)
opts <- list("algorithm" = "NLOPT_LN_NELDERMEAD", "xtol_rel"=1.0e-8,
             "xtol_abs" = 1.0e-8,
             "maxeval" = 10000, "print_level" = 1)
fit_Ast_log <-  nloptr(x0 = thetahat_st_nu, eval_f  = nll_student, opts = opts,
                       thetahat = NULL, y = data_gbc$time, d = data_gbc$status,
                       nu = NULL,  link ="log",
                       lb = c(1.0e-10, 1.0e-10, 1.0e-10, 1),
                       ub = c(Inf, Inf, 0.9999, Inf))
par_Ast_log <- fit_Ast_log$solution
nll_Ast_log <- fit_Ast_log$objective
AIC_Ast_log <- 2*4 + 2*nll_Ast_log
BIC_Ast_log <- 2*log(nrow(data_gbc)) + 2*nll_Ast_log

#thetahat <- par_st_log_nu$solution[1:3]
lambdahat_st <- optim(par = 0.015, fn  = nll_student,method =  "Brent",
                      thetahat = par_Ast_log, y = data_gbc$time, d = data_gbc$status,
                      nu = vhat,   link = "logit", lower = 0, upper = 1)$par
theta.start_st <- c(par_AN_logit, par_Ast_log[4])
#t <- c(19.865, 0.486, 0.697, 0.03573719, 80)
fit_Ast_logit <-  nloptr(x0 = theta.start_st, eval_f   = nll_student, opts = opts,
                         thetahat = NULL, y = data_gbc$time, d = data_gbc$status,
                         nu = NULL,  link ="logit",
                         lb = c(rep(1.0e-10, 4), 1),
                         ub = c(Inf, Inf, 0.99999, Inf, Inf))
par_Ast_logit <- fit_Ast_logit$solution
nll_Ast_logit <- fit_Ast_logit$objective
AIC_Ast_logit <- 2*5 + 2*nll_Ast_logit
BIC_Ast_logit <- 2*log(nrow(data_gbc)) + 2*nll_Ast_logit

#
nll <- rbind(nll_llogistic, nll_lnormal, nll_llaplace, nll_weibull,
             nll_AL_log, nll_AN_log, nll_ALa_log, nll_Ast_log,
             nll_AL_logit, nll_AN_logit, nll_ALa_logit, nll_Ast_logit)
AIC <- rbind(AIC_llogistic, AIC_lnormal, AIC_llaplace, AIC_weibull,
             AIC_AL_log, AIC_AN_log, AIC_ALa_log, AIC_Ast_log,
             AIC_AL_logit, AIC_AN_logit, AIC_ALa_logit, AIC_Ast_logit)
BIC <- rbind(BIC_llogistic, BIC_lnormal, BIC_llaplace, BIC_weibull,
             BIC_AL_log, BIC_AN_log, BIC_ALa_log, BIC_Ast_log,
             BIC_AL_logit, BIC_AN_logit, BIC_ALa_logit, BIC_Ast_logit)
Estimate <- rbind(par_llogistic, par_lnormal, par_llaplace, par_weibull,
                  par_AL_log, par_AN_log, par_ALa_log, par_Ast_log,
                  par_AL_logit, par_AN_logit, par_ALa_logit, par_Ast_logit)
gbc_par_nll_AIC_BIC <- cbind(Estimate, nll, AIC, BIC)
model <- c("llogistic", "lnormal", "llaplace", "weibull",
           "AL_log", "AN_log", "ALa_log", "Ast_log",
           "AL_logit", "AN_logit", "ALa_logit", "Ast_logit")
colN <- c("etahat","phihat", "alphahat", "lambdahat", "nuhat","nll", "AIC","BIC")
colnames(gbc_par_nll_AIC_BIC) <- paste("",colN)
rownames(gbc_par_nll_AIC_BIC) <- paste("",model)
gbc_par_nll_AIC_BIC[1:4,3:5] <- NA
gbc_par_nll_AIC_BIC[5:7,4:5] <- NA
gbc_par_nll_AIC_BIC[8,4] <- NA
gbc_par_nll_AIC_BIC[8,5] <- par_Ast_log[4]
gbc_par_nll_AIC_BIC[9:11,5] <- NA

#write.csv(gbc_par_nll_AIC_BIC, "breast_par_nll_AIC_BIC_updated.csv")
