rm(list = ls())
########################################################################################################
###                               Bivariate Dirichlet spectral density                               ###
########################################################################################################

########################################################################################################
### Log-likelihood and its 1st and 2nd derivatives (wrt its parameters \alpha and \beta) : 
### Needed for the NR step
########################################################################################################

#Below, w is the pseudo-observation (in the simplex) and alpha and beta the dependence parameters

#log-likelihood
log.h.dirichlet <- function(w, alpha, beta) {
  log((alpha*beta*(alpha*w)^(-1 + alpha)*(beta - beta*w)^(-1 + beta)*
         (beta + alpha*w - beta*w)^(-1 - alpha - beta)*gamma(1 + alpha + beta))/
        (gamma(alpha)*gamma(beta)))
}

#1st derivative of the log-likelihood
#wrt \alpha 
diri.d1.a <- function(w,alpha,beta){
  (1/(beta + alpha*w - beta*w))*(beta - w - 2*beta*w + beta*log(alpha*w) + 
                                   alpha*w*log(alpha*w) - beta*w*log(alpha*w) - beta*log(beta + alpha*w - beta*w) - 
                                   alpha*w*log(beta + alpha*w - beta*w) + beta*w*log(beta + alpha*w - beta*w) + 
                                   (beta*(-1 + w) - alpha*w)*digamma(alpha) + (beta + alpha*w - beta*w)*
                                   digamma(1 + alpha + beta))
}

#wrt \beta
diri.d1.b <- function(w,alpha,beta){
  (1/(beta*(-1 + w) - alpha*w))*(1 + alpha - w - 2*alpha*w - beta*log(beta - beta*w) - 
                                   alpha*w*log(beta - beta*w) + beta*w*log(beta - beta*w) + 
                                   beta*log(beta + alpha*w - beta*w) + alpha*w*log(beta + alpha*w - beta*w) - 
                                   beta*w*log(beta + alpha*w - beta*w) + (beta + alpha*w - beta*w)*digamma(beta) + 
                                   (beta*(-1 + w) - alpha*w)*digamma(1 + alpha + beta))
}

#2nd derivative of the log-likelihood
#wrt \alpha^2
diri.d2.a <- function(w,alpha,beta){
  (beta^2 - 2*beta^2*w + alpha*w^2 + alpha*beta*w^2 + beta^2*w^2 - 
     alpha*(beta + alpha*w - beta*w)^2*trigamma(alpha) + 
     alpha*(beta + alpha*w - beta*w)^2*trigamma(1 + alpha + beta))/
    (alpha*(beta + alpha*w - beta*w)^2)
}

#wrt \beta^2
diri.d2.b <- function(w,alpha,beta){
  (beta + alpha*beta - 2*beta*w - 2*alpha*beta*w + alpha^2*w^2 + beta*w^2 + 
     alpha*beta*w^2 - beta*(beta*(-1 + w) - alpha*w)^2*trigamma(beta) + 
     beta*(beta + alpha*w - beta*w)^2*trigamma(1 + alpha + beta))/
    (beta*(beta*(-1 + w) - alpha*w)^2)
}

#cross derivative
diri.d2.a.b <- function(w,alpha,beta){
  ((-beta)*(-1 + w)^2 - w*(-1 + w + alpha*w) + (beta*(-1 + w) - alpha*w)^2*
     trigamma(1 + alpha + beta))/(beta*(-1 + w) - alpha*w)^2
}

#Maximize the log-likelihood (needed to set the initial parameters)
maxLikelihood.dirichlet <- function(data, init = NULL, maxit = 500, hess = F)
{
  lhood <- function(alpha)
  {
    out <- NULL
    for(i in 1:nrow(data)){
      out[i] <- log.h.dirichlet(data[i,1]/(data[i,1]+data[i,2]),alpha[1],alpha[2])
    }
    return(-sum(out))
  }
  if(is.null(init)) {
    init <- c(1,1)
  }
  
  opt <- optim(init, lhood, method = "BFGS", control = list(maxit = maxit), hessian = hess)
  
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

dirichlet2d.vgam <- function (lshape1 = "loge", lshape2 = "loge", ishape1 = NULL, 
                              ishape2=NULL, tol12 = 1e-04, parallel = FALSE,
                              zero = NULL) 
{
  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")
  lshape2 <- as.list(substitute(lshape2))
  eshape2 <- link2list(lshape2)
  lshape2 <- attr(eshape2, "function.name")
  
  if (length(ishape1) && (!is.Numeric(ishape1))) 
    stop("bad input for argument 'ishape1'")
  if (length(ishape2) && !is.Numeric(ishape2)) 
    stop("bad input for argument 'ishape2'")
  new("vglmff", blurb = c("Dirichlet_2d distribution\n\n", "Links:    ", 
                          namesof("shape1", lshape1, eshape1, tag = FALSE), ", ", 
                          namesof("shape2", lshape2, eshape2, tag = FALSE), "\n"), 
      constraints = eval(substitute(expression({
        constraints <- cm.VGAM(matrix(1, M, 1), x = x, bool = .parallel, 
                               constraints = constraints)
        constraints <- cm.zero.VGAM(constraints, x = x, .zero, 
                                    M = M, predictors.names = predictors.names, M1 = 2)
      }), list(.parallel = parallel, .zero = zero))), infos = eval(substitute(function(...) {
        list(M1 = 2, Q1 = 1, expected = FALSE, multipleResponses = TRUE, 
             parameters.names = c("shape1", "shape2"), lshape1 = .lshape1, 
             lshape2 = .lshape2, parallel = .parallel, zero = .zero)
      }, list(.parallel = parallel, .zero = zero, .lshape1 = lshape1, .lshape2 = lshape2))), 
      initialize = eval(substitute(expression({
        checklist <- w.y.check(w = w, y = y, Is.positive.y = TRUE, 
                               ncol.w.max = Inf, ncol.y.max = Inf, out.wy = TRUE, 
                               colsyperw = 1, maximize = TRUE)
        w <- checklist$w
        y <- checklist$y
        if (any((y <= 0) | (y >= 1))) stop("the response must be in (0, 1)")
        extra$multiple.responses <- TRUE
        extra$ncoly <- ncoly <- ncol(y)
        extra$M1 <- M1 <- 2
        M <- M1 * ncoly
        mynames1 <- param.names("shape1", ncoly)
        mynames2 <- param.names("shape2", ncoly)
        predictors.names <- c(namesof(mynames1, .lshape1, earg = .eshape1, tag = FALSE), 
                              namesof(mynames2, .lshape2, earg = .eshape2, tag = FALSE)
        )[interleave.VGAM(M, M1 = M1)]
        if (!length(etastart)) {
          shape1.init <- .ishape1
          shape1.init <- matrix(shape1.init, n, ncoly, 
                                byrow = TRUE)
          shape2.init <- .ishape2
          shape2.init <- matrix(shape2.init, n, ncoly, 
                                byrow = TRUE)
          etastart <- cbind(theta2eta(shape1.init, .lshape1, earg = .eshape1), 
                            theta2eta(shape2.init, .lshape2, earg = .eshape2))[, interleave.VGAM(M, M1 = M1)]
        }
      }), list(.lshape1 = lshape1, .lshape2 = lshape2, .ishape1 = ishape1, 
               .ishape2 = ishape2, .eshape1 = eshape1, .eshape2 = eshape2 
      ))), linkinv = eval(substitute(function(eta, extra = NULL) {
        shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1, 
                            earg = .eshape1)
        shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2, 
                            earg = .eshape2)
        c(shape1,shape2)
      }, list(.lshape1 = lshape1, .lshape2 = lshape2, .eshape1 = eshape1, 
              .eshape2 = eshape2))), last = eval(substitute(expression({
                misc$link <- c(rep_len(.lshape1, ncoly), rep_len(.lshape2, 
                                                                 ncoly))[interleave.VGAM(M, M1 = M1)]
                temp.names <- c(mynames1, mynames2)[interleave.VGAM(M, M1 = M1)]
                names(misc$link) <- temp.names
                misc$earg        <- vector("list", M)
                names(misc$earg) <- temp.names
                misc$expected    <- FALSE
                misc$multiple.responses <- TRUE
                for (ii in 1:ncoly) {
                  misc$earg[[M1 * ii - 1]] <- .eshape1
                  misc$earg[[M1 * ii]]     <- .eshape2
                }
              }), list(.lshape1 = lshape1, .lshape2 = lshape2, .eshape1 = eshape1, 
                       .eshape2 = eshape2))), 
      loglikelihood = eval(substitute(function(mu, y, w, residuals = FALSE, eta, extra = NULL, summation = TRUE) {
        shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1, 
                            earg = .eshape1)
        shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2, 
                            earg = .eshape2)
        ll.elts <- c(w) * log.h.dirichlet(y, shape1, shape2)
        if (summation) sum(ll.elts) else ll.elts
      }, list(.lshape1 = lshape1, .lshape2 = lshape2, .eshape1 = eshape1, 
              .eshape2 = eshape2))), vfamily = c("dirichlet2d"), validparams = eval(substitute(function(eta, y, extra = NULL) {
                shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1, 
                                    earg = .eshape1)
                shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2, 
                                    earg = .eshape2)
                okay1 <- all(is.finite(shape1)) && all(0 < shape1) && 
                  all(is.finite(shape2)) && all(0 < shape2)
                okay1
              }, list(.lshape1 = lshape1, .lshape2 = lshape2, .eshape1 = eshape1, 
                      .eshape2 = eshape2))), simslot = eval(substitute(function(object,  nsim) {
                        eta <- predict(object)
                        shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1, 
                                            earg = .eshape1)
                        shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2, 
                                            earg = .eshape2)
                        rmevspec(n=nsim*length(shape1), d=2, param=c(shape1,shape2), model="ct")
                      }, list(.lshape1 = lshape1, .lshape2 = lshape2, .eshape1 = eshape1, 
                              .eshape2 = eshape2))), deriv = eval(substitute(expression({
                                shape1 <- eta2theta(eta[, c(TRUE, FALSE)], .lshape1, 
                                                    earg = .eshape1)
                                shape2 <- eta2theta(eta[, c(FALSE, TRUE)], .lshape2, 
                                                    earg = .eshape2)
                                dshape1.deta <- dtheta.deta(shape1, link = .lshape1, 
                                                            earg = .eshape1)
                                d2shape1.deta2 <- d2theta.deta2(shape1, link = .lshape1, 
                                                                earg = .eshape1)
                                dshape2.deta <- dtheta.deta(shape2, link = .lshape2, 
                                                            earg = .eshape2)
                                d2shape2.deta2 <- d2theta.deta2(shape2, link = .lshape2, 
                                                                earg = .eshape2)
                                dl.dshape1 <- diri.d1.a(y,shape1,shape2)
                                dl.dshape2 <- diri.d1.b(y,shape1,shape2)
                                
                                dl.deta <- c(w) * cbind(dl.dshape1 * dshape1.deta, 
                                                        dl.dshape2 * dshape2.deta)
                                dl.deta[, interleave.VGAM(M, M1 = M1)]
                              }), list(.lshape1 = lshape1, .lshape2 = lshape2, .eshape1 = eshape1, 
                                       .eshape2 = eshape2))), weight = eval(substitute(expression({
                                         ned2l.dshape11 <- (-diri.d2.a(y,shape1,shape2))
                                         ned2l.dshape22 <- (-diri.d2.b(y,shape1,shape2))
                                         ned2l.dshape12 <- (-diri.d2.a.b(y,shape1,shape2))
                                         
                                         wz <- array(c(c(w) * ((ned2l.dshape11 * dshape1.deta^2)-
                                                                 (diri.d1.a(y,shape1,shape2)*d2shape1.deta2)),
                                                       c(w) * ((ned2l.dshape22 * dshape2.deta^2)-
                                                                 (diri.d1.b(y,shape1,shape2)*d2shape2.deta2)),
                                                       c(w) * ned2l.dshape12 * dshape1.deta * dshape2.deta),
                                                     dim = c(n, M/M1, 3))
                                         wz <- arwz2wz(wz, M = M, M1 = M1)
                                         wz
                                       }), list(.lshape1 = lshape1, .lshape2 = lshape2, .eshape1 = eshape1, 
                                                .eshape2 = eshape2, .tol12 = tol12))))
}

########################################################################################################
### Examples of use
########################################################################################################
if(!require(VGAM)) install.packages("VGAM") 
if(!require(mev)) install.packages("mev") 

library("VGAM")
library("mev")

n       <- 600
t       <- seq(0.8,3.3,length=n)
alpha.1 <- exp(t)
alpha.2 <- exp(0.3*t^2)
beta    <- rep(3,length(t))

#simulate from the logistic spectral density
z      <- matrix(0,ncol=2,nrow=2*n)
set.seed(22)
for (i in 1:n)
{
  z[i,] <- rmevspec(n=1, d=2, param=c(alpha.1[i],beta[i]), model="ct")
}
set.seed(22)
for (i in (n+1):(2*n))
{
  z[i,] <- rmevspec(n=1, d=2, param=c(alpha.2[i-n],beta[i-n]), model="ct")
}

db       <- data.frame("z1"=z[,1],"z2"=z[,2],"t"=t)
z.data   <- data.frame("y1"=db$z1[1:n],"y2"=db$z1[(n+1):(2*n)],"t"=t)

#initial value = MLE 
par.init <- maxLikelihood.dirichlet(cbind(z.data$y1,1-z.data$y1))$par

#fit a VGAM to the first subset "y1" with "intercept-only" parameters (should be equal to the MLE)
diri.vgam.0 <- vgam(y1~1,family=dirichlet2d.vgam(ishape1 = par.init[1], 
                                                 ishape2 = par.init[2]), 
                    data=z.data, trace=TRUE)
unique(loge(unique(diri.vgam.0@predictors),inverse=TRUE))

# fit a VGAM to the first subset "y1" with an "intercept-only" beta parameter
diri.vgam.1 <- vgam(y1~sm.ps(t),family=dirichlet2d.vgam(ishape1 = par.init[1], 
                                                        ishape2 = par.init[2],
                                                        zero = 2),
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
CE.diri.h <- NULL
CE.diri.h <- hbvevd(Omega[,1], alpha=loge(diri.vgam.1@predictors[t.fixed,1], inverse=TRUE),
                    beta=loge(diri.vgam.1@predictors[t.fixed,2],inverse=TRUE),
                    model="ct",plot=FALSE,half=TRUE)

CE.diri.true.h <- NULL
CE.diri.true.h <- hbvevd(Omega[,1], alpha=alpha.1[t.fixed], beta=unique(beta),
                         model="ct",plot=FALSE, half=TRUE)

plot(Omega[,1], CE.diri.h,lwd=3,ylab=expression(paste(h[t],"(",w,")")),
     xlab="w",ylim=c(0,4),main=bquote("t" == .(round(db$t[t.fixed],2))), col=1,lty=2,type="l")
lines(Omega[,1],CE.diri.true.h,lwd=2)
#-------------------------------------------------------------------------------------------------------

# fit a VGAM to both subsets "y1" and "y2" with "intercept-only" \beta parameters
diri.vgam.2 <- vgam(cbind(y1,y2)~sm.ps(t),family=dirichlet2d.vgam(ishape1 = par.init[1], 
                                                                  ishape2 = par.init[2],
                                                                  zero=c(2,4)),
                   data=z.data, trace=TRUE)
par(mfrow=c(1,2))
plot(alpha.1,loge(diri.vgam.2@predictors[,1],inverse=TRUE),type="l",xlab=expression(alpha[1]),
     ylab=expression(hat(alpha)[1]))
abline(a=0,b=1,col=2)
plot(alpha.2,loge(diri.vgam.2@predictors[,3],inverse=TRUE),type="l",xlab=expression(alpha[2]),
     ylab=expression(hat(alpha)[2]))
abline(a=0,b=1,col=2)


# fit a VGAM to both subsets "y1" and "y2" with "intercept-only" \beta parameters 
# AND parallelism on the \alpha parameters
diri.vgam.3 <- vgam(cbind(y1,y2)~sm.ps(t),family=dirichlet2d.vgam(ishape1 = par.init[1], 
                                                                  ishape2 = par.init[2],
                                                                  parallel=TRUE,zero=c(2,4)),
                    data=z.data, trace=TRUE)
par(mfrow=c(1,2))
plot(alpha.1,loge(diri.vgam.3@predictors[,1],inverse=TRUE),type="l",xlab=expression(alpha[1]),
     ylab=expression(hat(alpha)[1]))
abline(a=0,b=1,col=2)
plot(alpha.2,loge(diri.vgam.3@predictors[,3],inverse=TRUE),type="l",xlab=expression(alpha[2]),
     ylab=expression(hat(alpha)[2]))
abline(a=0,b=1,col=2)

#the fit diri.vgam.3 is poorer than diri.vgam.2 as we are imposing the parallelism on the \alpha 
#parameters which is clearly a wrong assumption!

## An example of the use of the constraints argument + the parallelism on only one term
# fit a VGAM to both subsets "y1" and "y2" with parallelism only on the \alpha parameters
z.data.4   <- data.frame("y1"=db$z1[1:n],"y2"=db$z1[(n+1):(2*n)],"t1"=t, "t2"=t, "t3"=t)

diri.vgam.4 <- vgam(cbind(y1,y2)~sm.ps(t1)+sm.ps(t2),family=dirichlet2d.vgam(ishape1 = par.init[1], 
                                                                             ishape2 = par.init[2],
                                                                             parallel = TRUE ~ sm.ps(t1)-1),
                    data=z.data.4, trace=TRUE, constraints = list("(Intercept)" = diag(4),
                                                                  "sm.ps(t1)" = rbind(1,0,1,0), 
                                                                  "sm.ps(t2)" = cbind(c(0,1,0,0),c(0,0,0,1))))
