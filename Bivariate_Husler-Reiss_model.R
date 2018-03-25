rm(list = ls())
########################################################################################################
###                               Bivariate Husler-Reiss spectral density                            ###
########################################################################################################

########################################################################################################
### Log-likelihood and its 1st and 2nd derivatives (wrt its parameter \lambda) : 
### Needed for the NR step
########################################################################################################

#Below, w is the pseudo-observation (in the simplex) and alpha and beta the dependence parameters

#log-likelihood
log.h.HR <- function(x, lambda){
  (-(1/2))*log(2*pi) + log(lambda/(2*(-1 + x)^2*x)) - (2 + lambda^2*log(x/(1 - x)))^2/
    (8*lambda^2)
}

#1st derivative of the log-likelihood
#wrt \lambda 
HR.d1 <- function(x,lambda){
  1/lambda^3 + 1/lambda - (1/4)*lambda*log(x/(1 - x))^2
}

#2nd derivative of the log-likelihood
#wrt \lambda^2
HR.d2 <- function(x,lambda){
  -((3 + lambda^2)/lambda^4) - (1/4)*log(x/(1 - x))^2
}

#Maximize the log-likelihood (needed to set the initial parameter)
maxLikelihood.HR.spec <- function(data, init = NULL, maxit = 500, hess = F)
{
  
  lhood <- function(alpha)
  {
    out <- NULL
    for(i in 1:nrow(data)){
      out[i] <- log.h.HR(data[i,1],alpha)
    }
    return(-sum(out))
  }
  if(is.null(init))
  {
    init <- 0.5
  }
  
  min <- 1e-6
  max <- 1e6
  opt <- optim(init, lhood, method = "L-BFGS-B",lower = min, upper = max, control = list(maxit = maxit), hessian = hess)
  
  
  rtn <- list()
  rtn$par <- opt$par
  rtn$value <- opt$value
  rtn$convergence <- opt$convergence
  rtn$w <- data
  if(hess == F)
  {
    rtn$hessian <- 1
  }
  else
  {
    rtn$hessian <- opt$hessian
  }
  
  return(rtn)
}

########################################################################################################
### VGAM family
########################################################################################################

hr2d.vgam <- function (lshape1 = "loge", ishape1 = NULL, 
                       tol12 = 1e-04, parallel = FALSE, 
                       zero = NULL) 
{
  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")
  if (length(ishape1) && (!is.Numeric(ishape1))) 
    stop("bad input for argument 'ishape1'")
  if (!is.Numeric(tol12, length.arg = 1, positive = TRUE)) 
    stop("bad input for argument 'tol12'")
  
  new("vglmff", blurb = c("HR_2d distribution\n\n", "Links:    ", 
                          namesof("shape1", tag = FALSE), "\n"), 
      constraints = eval(substitute(expression({
        constraints <- cm.VGAM(matrix(1, M, 1), x = x, bool = .parallel, 
                               constraints = constraints)
        constraints <- cm.zero.VGAM(constraints, x = x, .zero, 
                                    M = M, predictors.names = predictors.names, M1 = 1)
      }), list(.parallel = parallel,.zero = zero))), infos = eval(substitute(function(...) {
        list(M1 = 1, Q1 = 1, expected = FALSE, multipleResponses = TRUE, 
             parameters.names = c("shape1"), lshape1 = .lshape1, 
             parallel = .parallel, zero= .zero)
      }, list(.parallel = parallel, .zero = zero, .lshape1 = lshape1))), 
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
        predictors.names <- c(namesof(mynames1, .lshape1, earg = .eshape1, tag = FALSE))[interleave.VGAM(M,M1 = M1)]
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
        misc$earg        <- vector("list", M)
        names(misc$earg) <- temp.names
        misc$expected    <- FALSE
        misc$multiple.responses <- TRUE
        for (ii in 1:ncoly) {
          misc$earg[[M1 * ii]] <- .eshape1
        }
      }), list(.lshape1 = lshape1, .eshape1 = eshape1))),
      loglikelihood = eval(substitute(function(mu, y, w, residuals = FALSE, eta, extra = NULL, summation = TRUE) {
        shape1  <- eta2theta(eta, .lshape1, earg = .eshape1)
        ll.elts <- c(w) * log.h.HR(y, shape1)
        if (summation) sum(ll.elts) else ll.elts
      }, list(.lshape1 = lshape1, .eshape1 = eshape1))),
      vfamily = c("hr2d"), validparams = eval(substitute(function(eta, y, extra = NULL) {
        shape1 <- eta2theta(eta, .lshape1, earg = .eshape1)
        okay1 <- all(is.finite(shape1)) && all(0 < (shape1))
        okay1
      }, list(.lshape1 = lshape1, .eshape1 = eshape1))), 
      simslot = eval(substitute(function(object, nsim) {
        eta    <- predict(object)
        shape1 <- eta2theta(eta, .lshape1, earg = .eshape1)
        rmevspec(n=nsim*length(shape1), d=2, sigma=matrix(c(0,shape1^2,shape1^2,0),byrow=TRUE,ncol=2), 
                 model="hr")[,1]
      }, list(.lshape1 = lshape1, .eshape1 = eshape1))), deriv = eval(substitute(expression({
        shape1 <- eta2theta(eta, .lshape1, earg = .eshape1)
        dshape1.deta <- dtheta.deta(shape1, link = .lshape1, earg = .eshape1)
        d2shape1.deta2 <- d2theta.deta2(shape1, link = .lshape1, earg = .eshape1)
        dl.dshape1 <- HR.d1(y,shape1)
        dl.deta <- c(w) * cbind(dl.dshape1 * dshape1.deta)
        dl.deta
      }), list(.lshape1 = lshape1, .eshape1 = eshape1))), 
      weight = eval(substitute(expression({
        ned2l.dshape11 <- -HR.d2(y,shape1)
        wz <- array(c(c(w) * cbind((ned2l.dshape11* dshape1.deta^2)-
                                     (HR.d1(y,shape1)*d2shape1.deta2))), 
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

n      <- 400
t      <- seq(0.1,2,length=n)
lambda <- exp(t)

#simulate from the Husler-Reiss spectral density
z      <- matrix(0,ncol=2,nrow=n)
set.seed(22)
for (i in 1:n)
{
  z[i,] <- rmevspec(n=1, d=2, sigma=matrix(c(0,1/lambda[i]^2,1/lambda[i]^2,0),byrow=TRUE,ncol=2),model="hr")
}

db       <- data.frame("z1"=z[,1],"z2"=z[,2],"t"=t)
z.data   <- data.frame("y"=db$z1/(db$z1+db$z2), "t"=t)

#initial value = MLE 
par.init <- maxLikelihood.HR.spec(cbind(db$z1,db$z2))$par

#fit a VGAM with an "intercept-only" parameter (should be equal to the MLE)
hr.vgam.0 <- vgam(y~1,family=hr2d.vgam(ishape1=par.init), data=z.data, trace=TRUE)
unique(fitted(hr.vgam.0))

# fit a VGAM to the response y (angular observations)
hr.vgam.1 <- vgam(y~sm.ps(t),family=hr2d.vgam(ishape1=par.init),
                   data=z.data, trace=TRUE)

#-------------------------------------------------------------------------------------------------------
#Generate n.simplex values from the simplex (d=2)
unit.simplex <- function(p,n)
{
  if(p<=0){stop("Impossible value for dimension p")}
  if(p==1){stop("Copulas not defined for univariate data")}
  if(p==2){
    output <- cbind(seq(0,1,length.out=n),1 - seq(0,1,length.out=n))
  }
  if (p==3){
    output <- matrix(NA,nrow=n,ncol=p)
    X <- matrix(rexp(n*p),nrow=n,ncol=p)
    
    for(j in 1:n){
      output[j,] <- X[j,]/sum(X[j,])
    }
  }
  return(output)
}
n.simplex <- 300
Omega     <- unit.simplex(2,n.simplex+2)
Omega     <- Omega[-c(1,nrow(Omega)),]

#compare the true and the estimated spectral density conditional on a fixed value of t
t.fixed <- 280
library("evd")
CE.hr.h <- NULL
CE.hr.h <- hbvevd(Omega[,1], dep=loge(hr.vgam.1@predictors[t.fixed], inverse=TRUE),
                    model="hr",plot=FALSE,half=TRUE)

CE.hr.true.h <- NULL
CE.hr.true.h <- hbvevd(Omega[,1], dep=lambda[t.fixed],
                         model="hr",plot=FALSE, half=TRUE)

#compute asymptotic confidence intervals for the spectral density conditional on the fixed value of t
Xp <- model.matrix(hr.vgam.1,type="lm")
b  <- coef(hr.vgam.1)
Vp <- vcov(hr.vgam.1)
#the vector \beta that contains the B-splines coefficients is asymptotically normal
br <- mvrnorm(10000, b, Vp)
#recover the the predictor \eta
calib <- br %*% t(Xp)

calib.CI.b.L <- matrix(NA,nrow=1,ncol=n.simplex)
calib.CI.b.U <- matrix(NA,nrow=1,ncol=n.simplex)

#this will take a while as the confidence intervals are computed pointwise
for(i in 1:n.simplex){
  h.asy      <- apply(calib,1,function(x){hbvevd(Omega[i,1], 
                                                  dep=VGAM::loge(x[t.fixed],inverse=TRUE), 
                                                  model="hr", half=TRUE)})
  CI.asy           <- quantile.fun(h.asy)
  calib.CI.b.L[,i] <- CI.asy[1]
  calib.CI.b.U[,i] <- CI.asy[2]
}

plot(Omega[,1], CE.hr.h,lwd=3,ylab=expression(paste(h[t],"(",w,")")),
     xlab="w",ylim=c(0,5),main=bquote("t" == .(round(db$t[t.fixed],2))), col=1,lty=2,type="n")
lines(Omega[,1],CE.hr.true.h,lwd=2)

polygon(c(Omega[,1],rev(Omega[,1])),
        c(calib.CI.b.L[1:n.simplex],rev(calib.CI.b.U[1:n.simplex])),col="lightgrey",border=NA)
lines(Omega[,1], CE.hr.true.h,lwd=2)
lines(Omega[,1], CE.hr.h,lwd=3,col=1,lty=2)


