rm(list = ls())
########################################################################################################
###                               Bivariate logistic spectral density                                ###
########################################################################################################

########################################################################################################
### Log-likelihood and its 1st and 2nd derivatives (wrt its parameter \alpha) : Needed for the NR step
########################################################################################################

#Below, x is the pseudo-observation (in the simplex) and y is the dependence parameter

#log-likelihood
bilogistic  <- function(x,y) {
  log((x - x^2)^(-1 - 1/y)*((1 - x)^(-y^(-1)) + x^(-y^(-1)))^(-2 + y)*(-1 + 1/y))
}

#1st derivative of the log-likelihood
dbilogistic <- function(x,y) {
  (x^(1/y)*(2 - 3*y + y^2)*log(1 - x) + (1 - x)^(1/y)*(2 - 3*y + y^2)*log(x) + 
     ((1 - x)^(1/y) + x^(1/y))*((-1 + y)*log((-(-1 + x))*x) + 
                                  y*(1 + (-1 + y)*y*log((1 - x)^(-y^(-1)) + x^(-y^(-1))))))/
    (((1 - x)^(1/y) + x^(1/y))*(-1 + y)*y^2)
}

#2nd derivative of the log-likelihood
d2bilogistic <- function(x,y) {
  (((-(-1 + x))*x)^(1/y)*(-2 + y)*(-1 + y)^2*log(1 - x)^2 + 
     4*((1 - x)^(2/y) + ((-(-1 + x))*x)^(1/y))*(-1 + y)^2*y*log(x) + 
     ((-(-1 + x))*x)^(1/y)*(-2 + y)*(-1 + y)^2*log(x)^2 + 
     2*(-1 + y)^2*log(1 - x)*(2*(x^(2/y) + ((-(-1 + x))*x)^(1/y))*y - 
                                ((-(-1 + x))*x)^(1/y)*(-2 + y)*log(x)) - 
     ((1 - x)^(2/y) + x^(2/y) + 2*((-(-1 + x))*x)^(1/y))*y*
     (y*(-1 + 2*y) + 2*(-1 + y)^2*log((-(-1 + x))*x)))/
    (((1 - x)^(1/y) + x^(1/y))^2*(-1 + y)^2*y^4)
}

#Maximize the log-likelihood (needed to set the initial parameters)

maxLikelihood.logistic.spec <- function(data, init = NULL, maxit = 500, hess = F) {
  lhood <- function(alpha) {
    out <- NULL
    for(i in 1:nrow(data)){
      out[i] <- bilogistic(data[i,1],alpha)
    }
    return(-sum(out))
  }
  if(is.null(init))
    init <- 0.1
  
  min <- 1e-6
  max <- 1 - 1e-6
  opt <- optim(init, lhood, method = "L-BFGS-B",lower = min, upper = max, control = list(maxit = maxit), hessian = hess)
  
  rtn <- list()
  rtn$par <- opt$par
  rtn$value <- opt$value
  rtn$convergence <- opt$convergence
  rtn$w <- data
  if(hess == F) {
    rtn$hessian <- 1
  }
  else {
    rtn$hessian <- opt$hessian
  }
  return(rtn)
}

########################################################################################################
### VGAM family
########################################################################################################

log2d.vgam <- function (lshape1 = "logit", ishape1 = NULL, parallel = FALSE, zero = NULL,
                        tol12 = 1e-04) {
  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")
  if (length(ishape1) && (!is.Numeric(ishape1, length.arg = 1, 
                                      positive = TRUE))) 
    stop("bad input for argument 'ishape1'")
  if (!is.Numeric(tol12, length.arg = 1, positive = TRUE)) 
    stop("bad input for argument 'tol12'")
  
  new("vglmff", blurb = c("logistic_2d distribution\n\n", "Links:    ", 
                          namesof("shape1", tag = FALSE), "\n"),
      constraints = eval(substitute(expression({
        constraints <- cm.VGAM(matrix(1, M, 1), x = x, bool = .parallel, 
                               constraints = constraints)
        constraints <- cm.zero.VGAM(constraints, x = x, .zero, 
                                    M = M, predictors.names = predictors.names, M1 = 1)
      }), list(.parallel = parallel, .zero = zero))),
      infos = eval(substitute(function(...) {
        list(M1 = 1, Q1 = 1, expected = FALSE, multipleResponses = TRUE, 
             parameters.names = c("shape1"), lshape1 = .lshape1, parallel = .parallel, zero= .zero)
      }, list(.lshape1 = lshape1, .parallel = parallel, .zero = zero))), 
      initialize = eval(substitute(expression({
        checklist <- w.y.check(w = w, y = y, Is.positive.y = TRUE, 
                               ncol.w.max = Inf, ncol.y.max = Inf, out.wy = TRUE, 
                               colsyperw = 1, maximize = TRUE)
        w <- checklist$w
        y <- checklist$y
        if (any((y <= 0) | (y >= 1))) stop("the response must be in (0, 1)")
        extra$multiple.responses <- TRUE
        extra$ncoly <- ncoly <- ncol(y)
        extra$M1 <- M1 <- 1
        M <- M1 * ncoly
        mynames1 <- param.names("shape1", ncoly)
        predictors.names <- c(namesof(mynames1,tag = FALSE))[interleave.VGAM(M,M1 = M1)]
        if (!length(etastart)) {
          shape1.init <- .ishape1 
          shape1.init <- matrix(shape1.init, n, ncoly, 
                                byrow = TRUE)
          etastart <- cbind(theta2eta(shape1.init, .lshape1, 
                                      earg = .eshape1))[, interleave.VGAM(M, M1 = M1)]
        }
      }), list(.lshape1 = lshape1, .ishape1 = ishape1, .eshape1 = eshape1))), 
      linkinv = eval(substitute(function(eta,extra = NULL) {
        shape1 <- eta2theta(eta, .lshape1, 
                            earg = .eshape1)
        shape1
      }, list(.lshape1 = lshape1, .eshape1 = eshape1))), last = eval(substitute(expression({
        misc$link <- c(rep_len(.lshape1, ncoly))[interleave.VGAM(M, M1 = M1)]
        temp.names <- c(mynames1)[interleave.VGAM(M, M1 = M1)]
        names(misc$link) <- temp.names
        misc$expected    <- FALSE
        misc$multiple.responses <- TRUE
        misc$earg <- vector("list", M)
        names(misc$earg) <- temp.names
        for (ii in 1:ncoly) {
          misc$earg[[M1 * ii]] <- .eshape1
        }
      }), list(.lshape1 = lshape1, .eshape1 = eshape1))),
      loglikelihood = eval(substitute(function(mu, y, w, residuals = FALSE, eta, extra = NULL, summation = TRUE) {
        shape1 <- eta2theta(eta, .lshape1, 
                            earg = .eshape1)
        ll.elts <- c(w) * bilogistic(y, shape1)
        if (summation) sum(ll.elts) else ll.elts
      }, list(.lshape1 = lshape1, .eshape1 = eshape1))),
      vfamily = c("logistic2d"), validparams = eval(substitute(function(eta, y, extra = NULL) {
        shape1 <- eta2theta(eta, .lshape1, 
                            earg = .eshape1)
        okay1 <- all(is.finite(shape1)) && all(0 < (shape1)) && all((shape1) <= 1)
        okay1
      }, list(.lshape1 = lshape1, .eshape1 = eshape1))), 
      simslot = eval(substitute(function(object, nsim) {
        eta <- predict(object)
        shape1 <- eta2theta(eta, .lshape1, 
                            earg = .eshape1)
        rmevspec(n=nsim*length(shape1), d=2, param=shape1, model="log")[,1]
      }, list(.lshape1 = lshape1, .eshape1 = eshape1))), deriv = eval(substitute(expression({
        shape1 <- eta2theta(eta, .lshape1, 
                            earg = .eshape1)
        dshape1.deta <- dtheta.deta(shape1, link = .lshape1, 
                                    earg = .eshape1)
        d2shape1.deta2 <- d2theta.deta2(shape1, link = .lshape1, 
                                        earg = .eshape1)
        dl.dshape1 <- dbilogistic(y,shape1)
        dl.deta <- c(w) * cbind(dl.dshape1 * dshape1.deta)
        dl.deta
      }), list(.lshape1 = lshape1, .eshape1 = eshape1))), 
      weight = eval(substitute(expression({
        ned2l.dshape11 <- -d2bilogistic(y,shape1)
        wz <- array(c(c(w) * cbind((ned2l.dshape11* dshape1.deta^2)-
                                     (dbilogistic(y,shape1)*d2shape1.deta2))), 
                    dim = c(n, M/M1, 1))
        wz <- arwz2wz(wz, M = M, M1 = M1)
        wz
      }), list(.lshape1 = lshape1, .eshape1 = eshape1, 
               .tol12 = tol12))))
}

########################################################################################################
### Examples of use
########################################################################################################
if(!require(VGAM)) install.packages("VGAM") 
if(!require(mev)) install.packages("mev") 

library("VGAM")
library("mev")

n <- 600
t <- seq(0.1,2,length=n)
#dependence parameter for the first subset (below y1)
alpha.1 <- VGAM::logit(t,inverse = TRUE)
#dependence parameter for the second subset (below y2)
alpha.2 <- VGAM::logit(t-0.5,inverse=TRUE)

#simulate from the logistic spectral density
z <- matrix(0,ncol=2,nrow=2*n)
set.seed(22)
for (i in 1:n)
{
  z[i,] <- rmevspec(n=1, d=2, param=alpha.1[i], model="log")
}
set.seed(22)
for (i in (n+1):(2*n))
{
  z[i,] <- rmevspec(n=1, d=2, param=alpha.2[i-n], model="log")
}

db       <- data.frame("z1"=z[,1],"z2"=z[,2],"t"=t)
z.data   <- data.frame("y1"=db$z1[1:n],"y2"=db$z1[(n+1):(2*n)],"t"=t)

#initial value = MLE 
par.init <- maxLikelihood.logistic.spec(cbind(z.data$y1,1-z.data$y1))$par

#fit a VGAM to the first subset "y1" with an "intercept-only" parameter (should be equal to the MLE)
log.vgam.0 <- vgam(y1~1,family=log2d.vgam(ishape1=par.init), data=z.data, trace=TRUE)
unique(fitted(log.vgam.0))

# fit a VGAM to the first subset "y1"
log.vgam.1 <- vgam(y1~sm.ps(t),family=log2d.vgam(ishape1=par.init),
                 data=z.data, trace=TRUE)

### compare the fitted values with the true ones

quantile.fun <- function(x) quantile(x, c((1 - conf)/2, 1 - (1 - conf)/2))
myCI         <- function(x) t(apply(x, 2, quantile.fun))
conf         <- .95

Xp <- model.matrix(log.vgam.1,type="lm")
b <- coef(log.vgam.1)
Vp <- vcov(log.vgam.1)
br <- mvrnorm(10000, b, Vp)
calib <- br %*% t(Xp)
calib.CI.b <- myCI(VGAM::logit(calib,inverse=TRUE))

plot(t,alpha.1,type="l",ylim=c(0,1),ylab=expression(alpha[1]),col=2);
lines(t,VGAM::logit(log.vgam.1@predictors,inverse = TRUE))
lines(t,calib.CI.b[,1],lty=2,lwd=2)
lines(t,calib.CI.b[,2],lty=2,lwd=2)
legend("bottomright", c("TRUE", "VGAM","CI"), col = c(2, 1, 1),
       lty = c(1,1,2),lwd=c(1,1,2),bg = "gray90")

# fit a VGAM to both subsets "y1" and "y2"
log.vgam.2 <- vgam(cbind(y1,y2)~sm.ps(t),family=log2d.vgam(ishape1=par.init),
                 data=z.data, trace=TRUE)

plot(alpha.1,fitted(log.vgam.2)[,1],xlab=expression(alpha[1]),ylab="Fitted values")
abline(a=0,b=1,col=2)

plot(alpha.2,fitted(log.vgam.2)[,2],xlab=expression(alpha[2]),ylab="Fitted values")
abline(a=0,b=1,col=2)

# fit a VGAM to both subsets "y1" and "y2" AND imposing the same effect of t on both subsets
log.vgam.3 <- vgam(cbind(y1,y2)~sm.ps(t),family=log2d.vgam(ishape1=par.init,parallel=TRUE),
                   data=z.data, trace=TRUE)

plot(alpha.1,fitted(log.vgam.3)[,1],type="l",xlab=expression(alpha[1]),ylab="Fitted values")
abline(a=0,b=1,col=2)

plot(alpha.2,fitted(log.vgam.3)[,2],type="l",xlab=expression(alpha[2]),ylab="Fitted values")
abline(a=0,b=1,col=2)

# fit a VGAM to both subsets "y1" and "y2" AND imposing an "intercept-only" dependence parameter in
# the second subset
log.vgam.4 <- vgam(cbind(y1,y2)~sm.ps(t),family=log2d.vgam(ishape1=par.init,zero=2),
                   data=z.data, trace=TRUE)

plot(alpha.1,fitted(log.vgam.4)[,1],type="l")
abline(a=0,b=1,col=2)
#the estimation in the first subset is not affected

plot(t,alpha.2,ylab=expression(alpha[2]),type="l")
lines(t,fitted(log.vgam.4)[,2],col=2)




