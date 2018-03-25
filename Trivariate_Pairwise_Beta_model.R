rm(list = ls())
########################################################################################################
###                            Trivariate Pairwise Beta spectral density                             ###
########################################################################################################

########################################################################################################
### Log-likelihood and its 1st and 2nd derivatives : Needed for the NR step
########################################################################################################

#Below, (w1,w2,1-w1-w2) is the pseudo-observation (in the simplex)
# alpha is the global dependence parameter
# beta1, beta2, and beta3 are the pairwise dependence parameters

#log-likelihood
pBeta.3d <- function(w1,w2,alpha,beta1,beta2,beta3){
  log(-((gamma(1 + 3*alpha)*((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*
                               (w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
                               gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                                                 (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                                                 (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                                                 gamma(beta2)^2*gamma(2*beta3))))/(w1*w2*(-1 + w1 + w2)*gamma(alpha)*
                                                                                     gamma(1 + 2*alpha)*gamma(beta1)^2*gamma(beta2)^2*gamma(beta3)^2)))
}

#1st derivative of the log-likelihood
#wrt alpha
pBeta.3d.d1.alpha <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
     gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2*(log(1 - w1 - w2) + 2*log(w1 + w2) - 
                                                     digamma(alpha) - 2*digamma(1 + 2*alpha) + 3*digamma(1 + 3*alpha)) + 
     gamma(beta1)^2*((1 - w1)^(2*alpha)*(-1 + w1)*w1^alpha*
                       (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)*
                       (2*log(1 - w1) + log(w1) - digamma(alpha) - 2*digamma(1 + 2*alpha) + 
                          3*digamma(1 + 3*alpha)) + (1 - w2)^(2*alpha)*(-1 + w2)*w2^alpha*
                       (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2*
                       (2*log(1 - w2) + log(w2) - digamma(alpha) - 2*digamma(1 + 2*alpha) + 
                          3*digamma(1 + 3*alpha))))/
    ((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
       gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 + 
       gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                         (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                         (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                         gamma(beta2)^2*gamma(2*beta3)))
}

#wrt beta1
pBeta.3d.d1.beta1 <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*
     gamma(beta2)^2*gamma(beta3)^2*(log((w1*w2)/(w1 + w2)^2) - 2*digamma(beta1) + 
                                      2*digamma(2*beta1)))/((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*
                                                              (w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
                                                              gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                                                                                (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                                                                                (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                                                                                gamma(beta2)^2*gamma(2*beta3)))
}

#wrt beta2
pBeta.3d.d1.beta2 <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((1 - w2)^(2*alpha)*(-1 + w2)*w2^alpha*(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*
     gamma(beta1)^2*gamma(2*beta2)*gamma(beta3)^2*
     (log(-((w1*(-1 + w1 + w2))/(-1 + w2)^2)) - 2*digamma(beta2) + 
        2*digamma(2*beta2)))/((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*
                                (w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 + 
                                gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                                                  (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                                                  (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                                                  gamma(beta2)^2*gamma(2*beta3)))
}

#wrt beta3
pBeta.3d.d1.beta3 <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((1 - w1)^(2*alpha)*(-1 + w1)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
     gamma(beta1)^2*gamma(beta2)^2*gamma(2*beta3)*
     (log(-((w2*(-1 + w1 + w2))/(-1 + w1)^2)) - 2*digamma(beta3) + 
        2*digamma(2*beta3)))/((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*
                                (w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 + 
                                gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                                                  (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                                                  (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                                                  gamma(beta2)^2*gamma(2*beta3)))
}

#2nd derivative of the log-likelihood
pBeta.3d.d2.alpha <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((-((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
        gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2*log(1 - w1 - w2) + 
        gamma(beta1)^2*((1 - w1)^(2*alpha)*(-1 + w1)*w1^alpha*
                          (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)*
                          (2*log(1 - w1) + log(w1)) + (1 - w2)^(2*alpha)*(-1 + w2)*w2^alpha*
                          (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2*
                          (2*log(1 - w2) + log(w2))) - 2*(1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*
        (w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2*
        log(w1 + w2)))*((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*
                          (w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2*
                          (log(1 - w1 - w2) + 2*log(w1 + w2) - digamma(alpha) - 
                             2*digamma(1 + 2*alpha) + 3*digamma(1 + 3*alpha)) + 
                          gamma(beta1)^2*((1 - w1)^(2*alpha)*(-1 + w1)*w1^alpha*
                                            (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)*
                                            (2*log(1 - w1) + log(w1) - digamma(alpha) - 2*digamma(1 + 2*alpha) + 
                                               3*digamma(1 + 3*alpha)) + (1 - w2)^(2*alpha)*(-1 + w2)*w2^alpha*
                                            (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2*
                                            (2*log(1 - w2) + log(w2) - digamma(alpha) - 2*digamma(1 + 2*alpha) + 
                                               3*digamma(1 + 3*alpha)))) + 
     ((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
        gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 + 
        gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                          (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                          (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                          gamma(beta2)^2*gamma(2*beta3)))*
     ((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
        gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2*log(1 - w1 - w2)*
        (log(1 - w1 - w2) + 2*log(w1 + w2) - digamma(alpha) - 
           2*digamma(1 + 2*alpha) + 3*digamma(1 + 3*alpha)) - 
        2*(1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
        gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2*log(w1 + w2)*
        (log(1 - w1 - w2) + 2*log(w1 + w2) - digamma(alpha) - 
           2*digamma(1 + 2*alpha) + 3*digamma(1 + 3*alpha)) + 
        gamma(beta1)^2*(-2*(1 - w1)^(1 + 2*alpha)*w1^alpha*
                          (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)*
                          log(1 - w1)*(2*log(1 - w1) + log(w1) - digamma(alpha) - 
                                         2*digamma(1 + 2*alpha) + 3*digamma(1 + 3*alpha)) + 
                          (1 - w1)^(2*alpha)*(-1 + w1)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                          gamma(beta2)^2*gamma(2*beta3)*log(w1)*(2*log(1 - w1) + log(w1) - 
                                                                   digamma(alpha) - 2*digamma(1 + 2*alpha) + 
                                                                   3*digamma(1 + 3*alpha)) - 2*(1 - w2)^(1 + 2*alpha)*w2^alpha*
                          (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2*
                          log(1 - w2)*(2*log(1 - w2) + log(w2) - digamma(alpha) - 
                                         2*digamma(1 + 2*alpha) + 3*digamma(1 + 3*alpha)) + 
                          (1 - w2)^(2*alpha)*(-1 + w2)*w2^alpha*(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*
                          gamma(2*beta2)*gamma(beta3)^2*log(w2)*(2*log(1 - w2) + log(w2) - 
                                                                   digamma(alpha) - 2*digamma(1 + 2*alpha) + 
                                                                   3*digamma(1 + 3*alpha)) + (1 - w2)^(1 + 2*alpha)*w2^alpha*
                          (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2*
                          (trigamma(alpha) + 4*trigamma(1 + 2*alpha) - 
                             9*trigamma(1 + 3*alpha)) + (1 - w1)^(1 + 2*alpha)*w1^alpha*
                          (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)*
                          (trigamma(alpha) + 4*trigamma(1 + 2*alpha) - 
                             9*trigamma(1 + 3*alpha))) + (1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*
        (w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2*
        (trigamma(alpha) + 4*trigamma(1 + 2*alpha) - 
           9*trigamma(1 + 3*alpha))))/
    ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
       gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
       gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                         (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                         (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                         gamma(beta2)^2*gamma(2*beta3)))^2
}

pBeta.3d.d2.beta1 <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*
     gamma(beta2)^2*gamma(beta3)^2*
     (((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
         gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - gamma(beta1)^2*
         ((-(1 - w2)^(1 + 2*alpha))*w2^alpha*(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*
            gamma(2*beta2)*gamma(beta3)^2 - (1 - w1)^(1 + 2*alpha)*w1^alpha*
            (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)))*
        log((w1*w2)/(w1 + w2)^2)*(log((w1*w2)/(w1 + w2)^2) - 2*digamma(beta1) + 
                                    2*digamma(2*beta1)) + 2*((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*
                                                               (w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
                                                               gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                                                                                 (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                                                                                 (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                                                                                 gamma(beta2)^2*gamma(2*beta3)))*digamma(2*beta1)*
        (log((w1*w2)/(w1 + w2)^2) - 2*digamma(beta1) + 2*digamma(2*beta1)) - 
        (log((w1*w2)/(w1 + w2)^2) - 2*digamma(beta1) + 2*digamma(2*beta1))*
        (-2*gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                              (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                              (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                              gamma(beta2)^2*gamma(2*beta3))*digamma(beta1) + 
           (1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
           gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2*(log((w1*w2)/(w1 + w2)^2) + 
                                                           2*digamma(2*beta1))) + ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*
                                                                                     (w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
                                                                                     gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                                                                                                       (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                                                                                                       (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                                                                                                       gamma(beta2)^2*gamma(2*beta3)))*(-2*trigamma(beta1) + 
                                                                                                                                          4*trigamma(2*beta1))))/
    ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
       gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
       gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                         (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                         (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                         gamma(beta2)^2*gamma(2*beta3)))^2
}

pBeta.3d.d2.beta2 <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((1 - w2)^(2*alpha)*(-1 + w2)*w2^alpha*(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*
     gamma(beta1)^2*gamma(2*beta2)*gamma(beta3)^2*
     (((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
         gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 + gamma(beta1)^2*
         ((-(1 - w2)^(1 + 2*alpha))*w2^alpha*(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*
            gamma(2*beta2)*gamma(beta3)^2 - (1 - w1)^(1 + 2*alpha)*w1^alpha*
            (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)))*
        log(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))*(log(-((w1*(-1 + w1 + w2))/(-1 + w2)^2)) - 
                                                   2*digamma(beta2) + 2*digamma(2*beta2)) + 
        2*((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
             gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 + gamma(beta1)^2*
             ((-(1 - w2)^(1 + 2*alpha))*w2^alpha*(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*
                gamma(2*beta2)*gamma(beta3)^2 - (1 - w1)^(1 + 2*alpha)*w1^alpha*
                (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)))*
        digamma(2*beta2)*(log(-((w1*(-1 + w1 + w2))/(-1 + w2)^2)) - 
                            2*digamma(beta2) + 2*digamma(2*beta2)) - 
        (log(-((w1*(-1 + w1 + w2))/(-1 + w2)^2)) - 2*digamma(beta2) + 
           2*digamma(2*beta2))*(-2*(1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*
                                  (w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2*
                                  digamma(beta2) + gamma(beta1)^2*(-2*(1 - w1)^(1 + 2*alpha)*w1^alpha*
                                                                     (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)*
                                                                     digamma(beta2) + (1 - w2)^(2*alpha)*(-1 + w2)*w2^alpha*
                                                                     (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2*
                                                                     (log(-((w1*(-1 + w1 + w2))/(-1 + w2)^2)) + 2*digamma(2*beta2)))) + 
        ((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
           gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 + gamma(beta1)^2*
           ((-(1 - w2)^(1 + 2*alpha))*w2^alpha*(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*
              gamma(2*beta2)*gamma(beta3)^2 - (1 - w1)^(1 + 2*alpha)*w1^alpha*
              (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)))*
        (-2*trigamma(beta2) + 4*trigamma(2*beta2))))/
    ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
       gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
       gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                         (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                         (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                         gamma(beta2)^2*gamma(2*beta3)))^2
}

pBeta.3d.d2.beta3 <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((1 - w1)^(2*alpha)*(-1 + w1)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
     gamma(beta1)^2*gamma(beta2)^2*gamma(2*beta3)*
     (((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
         gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 + gamma(beta1)^2*
         ((-(1 - w2)^(1 + 2*alpha))*w2^alpha*(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*
            gamma(2*beta2)*gamma(beta3)^2 - (1 - w1)^(1 + 2*alpha)*w1^alpha*
            (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)))*
        log(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))*(log(-((w2*(-1 + w1 + w2))/(-1 + w1)^2)) - 
                                                   2*digamma(beta3) + 2*digamma(2*beta3)) + 
        2*((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
             gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 + gamma(beta1)^2*
             ((-(1 - w2)^(1 + 2*alpha))*w2^alpha*(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*
                gamma(2*beta2)*gamma(beta3)^2 - (1 - w1)^(1 + 2*alpha)*w1^alpha*
                (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)))*
        digamma(2*beta3)*(log(-((w2*(-1 + w1 + w2))/(-1 + w1)^2)) - 
                            2*digamma(beta3) + 2*digamma(2*beta3)) - 
        (log(-((w2*(-1 + w1 + w2))/(-1 + w1)^2)) - 2*digamma(beta3) + 
           2*digamma(2*beta3))*(-2*(1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*
                                  (w1 + w2)^(1 + 2*alpha)*gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2*
                                  digamma(beta3) + gamma(beta1)^2*(-2*(1 - w2)^(1 + 2*alpha)*w2^alpha*
                                                                     (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2*
                                                                     digamma(beta3) + (1 - w1)^(2*alpha)*(-1 + w1)*w1^alpha*
                                                                     (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)*
                                                                     (log(-((w2*(-1 + w1 + w2))/(-1 + w1)^2)) + 2*digamma(2*beta3)))) + 
        ((-(1 - w1 - w2)^alpha)*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
           gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 + gamma(beta1)^2*
           ((-(1 - w2)^(1 + 2*alpha))*w2^alpha*(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*
              gamma(2*beta2)*gamma(beta3)^2 - (1 - w1)^(1 + 2*alpha)*w1^alpha*
              (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*gamma(beta2)^2*gamma(2*beta3)))*
        (-2*trigamma(beta3) + 4*trigamma(2*beta3))))/
    ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
       gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
       gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                         (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                         (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                         gamma(beta2)^2*gamma(2*beta3)))^2
}

pBeta.3d.d2.alpha.beta1 <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*gamma(beta1)^2*
     gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2*
     ((1 - w1)^(2*alpha)*(-1 + w1)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
        gamma(beta2)^2*gamma(2*beta3)*(2*log(1 - w1) + log(w1) - log(1 - w1 - w2) - 
                                         2*log(w1 + w2)) + (1 - w2)^(2*alpha)*(-1 + w2)*w2^alpha*
        (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2*
        (2*log(1 - w2) - log(1 - w1 - w2) + log(w2) - 2*log(w1 + w2)))*
     (log((w1*w2)/(w1 + w2)^2) - 2*digamma(beta1) + 2*digamma(2*beta1)))/
    ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
       gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
       gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                         (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                         (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                         gamma(beta2)^2*gamma(2*beta3)))^2
}

pBeta.3d.d2.beta1.beta2 <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((1 - w2)^(2*alpha)*(1 - w1 - w2)^alpha*(-1 + w2)*w2^alpha*
     (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*((w1*w2)/(w1 + w2)^2)^beta1*
     (w1 + w2)^(1 + 2*alpha)*gamma(beta1)^2*gamma(2*beta1)*gamma(beta2)^2*gamma(2*beta2)*
     gamma(beta3)^4*(log((w1*w2)/(w1 + w2)^2) - 2*digamma(beta1) + 
                       2*digamma(2*beta1))*(log(-((w1*(-1 + w1 + w2))/(-1 + w2)^2)) - 
                                              2*digamma(beta2) + 2*digamma(2*beta2)))/
    ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
       gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
       gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                         (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                         (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                         gamma(beta2)^2*gamma(2*beta3)))^2
}

pBeta.3d.d2.beta2.beta3 <- function(w1,w2,alpha,beta1,beta2,beta3){
  -(((1 - w1)^(1 + 2*alpha)*w1^alpha*(1 - w2)^(1 + 2*alpha)*w2^alpha*
       (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^
       beta3*gamma(beta1)^4*gamma(beta2)^2*gamma(2*beta2)*gamma(beta3)^2*gamma(2*beta3)*
       (log(-((w1*(-1 + w1 + w2))/(-1 + w2)^2)) - 2*digamma(beta2) + 
          2*digamma(2*beta2))*(log(-((w2*(-1 + w1 + w2))/(-1 + w1)^2)) - 
                                 2*digamma(beta3) + 2*digamma(2*beta3)))/
      ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
         gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
         gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                           (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                           (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                           gamma(beta2)^2*gamma(2*beta3)))^2)
}

pBeta.3d.d2.alpha.beta2 <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((1 - w2)^(1 + 2*alpha)*w2^alpha*(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*
     gamma(beta1)^2*gamma(beta2)^2*gamma(2*beta2)*gamma(beta3)^2*
     ((1 - w1)^(2*alpha)*(-1 + w1)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
        gamma(beta1)^2*gamma(2*beta3)*(2*log(1 - w1) + log(w1) - 2*log(1 - w2) - log(w2)) + 
        (1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
        gamma(2*beta1)*gamma(beta3)^2*(2*log(1 - w2) - log(1 - w1 - w2) + log(w2) - 
                                         2*log(w1 + w2)))*(log(-((w1*(-1 + w1 + w2))/(-1 + w2)^2)) - 
                                                             2*digamma(beta2) + 2*digamma(2*beta2)))/
    ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
       gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
       gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                         (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                         (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                         gamma(beta2)^2*gamma(2*beta3)))^2
}

pBeta.3d.d2.beta1.beta3 <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((1 - w1)^(2*alpha)*(-1 + w1)*w1^alpha*(1 - w1 - w2)^alpha*
     (-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*((w1*w2)/(w1 + w2)^2)^beta1*
     (w1 + w2)^(1 + 2*alpha)*gamma(beta1)^2*gamma(2*beta1)*gamma(beta2)^4*gamma(beta3)^2*
     gamma(2*beta3)*(log((w1*w2)/(w1 + w2)^2) - 2*digamma(beta1) + 
                       2*digamma(2*beta1))*(log(-((w2*(-1 + w1 + w2))/(-1 + w1)^2)) - 
                                              2*digamma(beta3) + 2*digamma(2*beta3)))/
    ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
       gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
       gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                         (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                         (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                         gamma(beta2)^2*gamma(2*beta3)))^2
}

pBeta.3d.d2.alpha.beta3 <- function(w1,w2,alpha,beta1,beta2,beta3){
  ((1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
     gamma(beta1)^2*gamma(beta2)^2*gamma(beta3)^2*gamma(2*beta3)*
     ((1 - w2)^(1 + 2*alpha)*w2^alpha*(-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*
        gamma(beta1)^2*gamma(2*beta2)*(2*log(1 - w1) + log(w1) - 2*log(1 - w2) - log(w2)) + 
        (1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
        gamma(2*beta1)*gamma(beta2)^2*(2*log(1 - w1) + log(w1) - log(1 - w1 - w2) - 
                                         2*log(w1 + w2)))*(log(-((w2*(-1 + w1 + w2))/(-1 + w1)^2)) - 
                                                             2*digamma(beta3) + 2*digamma(2*beta3)))/
    ((1 - w1 - w2)^alpha*((w1*w2)/(w1 + w2)^2)^beta1*(w1 + w2)^(1 + 2*alpha)*
       gamma(2*beta1)*gamma(beta2)^2*gamma(beta3)^2 - 
       gamma(beta1)^2*((-(1 - w2)^(1 + 2*alpha))*w2^alpha*
                         (-((w1*(-1 + w1 + w2))/(-1 + w2)^2))^beta2*gamma(2*beta2)*gamma(beta3)^2 - 
                         (1 - w1)^(1 + 2*alpha)*w1^alpha*(-((w2*(-1 + w1 + w2))/(-1 + w1)^2))^beta3*
                         gamma(beta2)^2*gamma(2*beta3)))^2
}

#Maximize the log-likelihood (needed to set the initial parameters)
maxLikelihood.pb <- function(data, model = c("pb", "diri"), init = NULL, maxit = 500, hess = F)
{
  p = dim(data)[2]
  
  lhood <- function(vec)
  {
    out <- dpairbeta(data,par = c(vec[1],vec[2:(choose(p,2)+1)]),log = TRUE)
    return(-sum(out))
  }
  if(is.null(init))
  {
    init <- c(1, rep(1, choose(p, 2)))
  }
  
  min <- rep(1e-5, choose(p,2)+1)
  max <- rep(1e6, length(init))
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
  
  #likelihood contributions
  par <- opt$par
  likeContrib <- dpairbeta(data,par = c(par[1],par[2:(choose(p,2)+1)]),log = TRUE)
  rtn$likeContrib <- likeContrib
  
  return(rtn)
}

########################################################################################################
### VGAM family
########################################################################################################

pBeta3d.vgam <- function(lshape1 = "loge", lshape2 = "loge", lshape3 = "loge", lshape4 = "loge", 
                         ishape1 = NULL, ishape2=NULL, ishape3=NULL, ishape4=NULL, tol12 = 1e-04,
                         zero = NULL, parallel = FALSE) 
{
  lshape1 <- as.list(substitute(lshape1))
  eshape1 <- link2list(lshape1)
  lshape1 <- attr(eshape1, "function.name")
  lshape2 <- as.list(substitute(lshape2))
  eshape2 <- link2list(lshape2)
  lshape2 <- attr(eshape2, "function.name")
  lshape3 <- as.list(substitute(lshape3))
  eshape3 <- link2list(lshape3)
  lshape3 <- attr(eshape3, "function.name")
  lshape4 <- as.list(substitute(lshape4))
  eshape4 <- link2list(lshape4)
  lshape4 <- attr(eshape4, "function.name")
  if (length(ishape1) && (!is.Numeric(ishape1))) 
    stop("bad input for argument 'ishape1'")
  if (length(ishape2) && !is.Numeric(ishape2)) 
    stop("bad input for argument 'ishape2'")
  if (length(ishape3) && !is.Numeric(ishape3)) 
    stop("bad input for argument 'ishape3'")
  if (length(ishape3) && !is.Numeric(ishape4)) 
    stop("bad input for argument 'ishape4'")
  new("vglmff", blurb = c("pBeta_3d distribution\n\n", "Links:    ", 
                          namesof("shape1", lshape1, eshape1, tag = FALSE), ", ", 
                          namesof("shape2", lshape2, eshape2, tag = FALSE),", ", 
                          namesof("shape3", lshape3, eshape3, tag = FALSE),", ",
                          namesof("shape4", lshape4, eshape4, tag = FALSE), "\n"), 
      constraints = eval(substitute(expression({
        constraints <- cm.VGAM(matrix(1, M, 1), x = x, bool = .parallel, 
                               constraints = constraints)
        constraints <- cm.zero.VGAM(constraints, x = x, .zero, 
                                    M = M, predictors.names = predictors.names, M1 = 4)
      }), list(.parallel = parallel, .zero = zero))), infos = eval(substitute(function(...) {
        list(M1 = 4, Q1 = 2, expected = FALSE, multipleResponses = TRUE, 
             parameters.names = c("shape1", "shape2", "shape3", "shape4"), lshape1 = .lshape1, 
             lshape2 = .lshape2, lshape3 = .lshape3, lshape4 = .lshape4, zero = .zero, parallel = .parallel)
      }, list(.parallel = parallel, .zero = zero, .lshape1 = lshape1, .lshape2 = lshape2, .lshape3 = lshape3, .lshape4 = lshape4))), 
      initialize = eval(substitute(expression({
        Q1        <- 2
        checklist <- w.y.check(w = w, y = y, Is.positive.y = TRUE, 
                               ncol.w.max = Inf, ncol.y.min = Q1, ncol.y.max = Inf, out.wy = TRUE, 
                               colsyperw = Q1, maximize = TRUE)
        w <- checklist$w
        y <- checklist$y
        if (any((y <= 0) | (y >= 1))) stop("the response must be in (0, 1)")
        extra$multiple.responses <- TRUE
        extra$ncoly <- ncoly <- ncol(y)
        extra$M1 <- M1 <- 1+choose(Q1+1,2)
        M <- M1 * (ncoly/Q1)
        mynames1 <- param.names("shape1", ncoly/Q1)
        mynames2 <- param.names("shape2", ncoly/Q1)
        mynames3 <- param.names("shape3", ncoly/Q1)
        mynames4 <- param.names("shape4", ncoly/Q1)
        predictors.names <- c(namesof(mynames1, .lshape1, earg = .eshape1, tag = FALSE), 
                              namesof(mynames2, .lshape2, earg = .eshape2, tag = FALSE),
                              namesof(mynames3, .lshape3, earg = .eshape3, tag = FALSE),
                              namesof(mynames4, .lshape4, earg = .eshape4, tag = FALSE))[interleave.VGAM(M, M1 = M1)]
        if (!length(etastart)) {
          shape1.init <- .ishape1
          shape1.init <- matrix(shape1.init, n, ncoly/Q1, 
                                byrow = TRUE)
          shape2.init <- .ishape2
          shape2.init <- matrix(shape2.init, n, ncoly/Q1, 
                                byrow = TRUE)
          shape3.init <- .ishape3
          shape3.init <- matrix(shape3.init, n, ncoly/Q1, 
                                byrow = TRUE)
          shape4.init <- .ishape4
          shape4.init <- matrix(shape4.init, n, ncoly/Q1, 
                                byrow = TRUE)
          etastart <- cbind(theta2eta(shape1.init, .lshape1, earg = .eshape1), 
                            theta2eta(shape2.init, .lshape2, earg = .eshape2),
                            theta2eta(shape3.init, .lshape3, earg = .eshape3),
                            theta2eta(shape4.init, .lshape4, earg = .eshape4))[, interleave.VGAM(M, M1 = M1)]
        }
      }), list(.lshape1 = lshape1, .lshape2 = lshape2, .lshape3 = lshape3, .lshape4 = lshape4, 
               .ishape1 = ishape1, .ishape2 = ishape2, .ishape3 = ishape3, .ishape4 = ishape4, 
               .eshape1 = eshape1, .eshape2 = eshape2, .eshape3 = eshape3, .eshape4 = eshape4))), 
      linkinv = eval(substitute(function(eta, extra = NULL) {
        shape1 <- eta2theta(eta[, c(TRUE, FALSE, FALSE, FALSE)], .lshape1, 
                            earg = .eshape1)
        shape2 <- eta2theta(eta[, c(FALSE, TRUE, FALSE, FALSE)], .lshape2, 
                            earg = .eshape2)
        shape3 <- eta2theta(eta[, c(FALSE, FALSE, TRUE, FALSE)], .lshape3, 
                            earg = .eshape3)
        shape4 <- eta2theta(eta[, c(FALSE, FALSE, FALSE, TRUE)], .lshape4, 
                            earg = .eshape4)
        c(shape1,shape2,shape3,shape4)
      }, list(.lshape1 = lshape1, .lshape2 = lshape2, .lshape3 = lshape3, .lshape4 = lshape4, 
              .eshape1 = eshape1, .eshape2 = eshape2, .eshape3 = eshape3, .eshape4 = eshape4))), 
      last = eval(substitute(expression({
        misc$link <- c(rep_len(.lshape1, ncoly/Q1), rep_len(.lshape2,ncoly/Q1),
                       rep_len(.lshape3, ncoly/Q1), rep_len(.lshape4, ncoly/Q1))[interleave.VGAM(M, M1 = M1)]
        temp.names <- c(mynames1, mynames2, mynames3, mynames4)[interleave.VGAM(M, 
                                                                                M1 = M1)]
        names(misc$link) <- temp.names
        misc$earg        <- vector("list", M)
        names(misc$earg) <- temp.names
        misc$expected    <- FALSE
        misc$multiple.responses <- TRUE
        for (ii in 1:(ncoly/Q1)) {
          misc$earg[[M1 * ii - 3]] <- .eshape1
          misc$earg[[M1 * ii - 2]] <- .eshape2
          misc$earg[[M1 * ii - 1]] <- .eshape3
          misc$earg[[M1 * ii]]     <- .eshape4
        }
      }), list(.lshape1 = lshape1, .lshape2 = lshape2, .lshape3 = lshape3, .lshape4 = lshape4, 
               .eshape1 = eshape1, .eshape2 = eshape2, .eshape3 = eshape3, .eshape4 = eshape4))), 
      loglikelihood = eval(substitute(function(mu, y, w, residuals = FALSE, eta, extra = NULL, summation = TRUE) {
        shape1 <- eta2theta(eta[, c(TRUE, FALSE, FALSE, FALSE)], .lshape1, 
                            earg = .eshape1)
        shape2 <- eta2theta(eta[, c(FALSE, TRUE, FALSE, FALSE)], .lshape2, 
                            earg = .eshape2)
        shape3 <- eta2theta(eta[, c(FALSE, FALSE, TRUE, FALSE)], .lshape3, 
                            earg = .eshape3)
        shape4 <- eta2theta(eta[, c(FALSE, FALSE, FALSE, TRUE)], .lshape4, 
                            earg = .eshape4)
        
        ll.elts <- c(w) * pBeta.3d(y[,1],y[,2],shape1,shape2,shape3,shape4)
        if (summation) sum(ll.elts) else ll.elts
      }, list(.lshape1 = lshape1, .lshape2 = lshape2, .lshape3 = lshape3, .lshape4 = lshape4,
              .eshape1 = eshape1, .eshape2 = eshape2, .eshape3 = eshape3, .eshape4 = eshape4))), 
      vfamily = c("pBeta3d"), validparams = eval(substitute(function(eta, 
                                                                     y, extra = NULL) {
        shape1 <- eta2theta(eta[, c(TRUE, FALSE, FALSE, FALSE)], .lshape1, 
                            earg = .eshape1)
        shape2 <- eta2theta(eta[, c(FALSE, TRUE, FALSE, FALSE)], .lshape2, 
                            earg = .eshape2)
        shape3 <- eta2theta(eta[, c(FALSE, FALSE, TRUE, FALSE)], .lshape3, 
                            earg = .eshape3)
        shape4 <- eta2theta(eta[, c(FALSE, FALSE, FALSE, TRUE)], .lshape4, 
                            earg = .eshape4)
        okay1 <- all(is.finite(shape1)) && all(0 < shape1) && 
          all(is.finite(shape2)) && all(0 < shape2) &&
          all(is.finite(shape3)) && all(0 < shape3) &&
          all(is.finite(shape4)) && all(0 < shape4)
        okay1
      }, list(.lshape1 = lshape1, .lshape2 = lshape2, .lshape3 = lshape3, .lshape4 = lshape4, 
              .eshape1 = eshape1, .eshape2 = eshape2, .eshape3 = eshape3, .eshape4 = eshape4))), 
      simslot = eval(substitute(function(object, nsim) {
        eta <- predict(object)
        shape1 <- eta2theta(eta[, c(TRUE, FALSE, FALSE, FALSE)], .lshape1, 
                            earg = .eshape1)
        shape2 <- eta2theta(eta[, c(FALSE, TRUE, FALSE, FALSE)], .lshape2, 
                            earg = .eshape2)
        shape3 <- eta2theta(eta[, c(FALSE, FALSE, TRUE, FALSE)], .lshape3, 
                            earg = .eshape3)
        shape4 <- eta2theta(eta[, c(FALSE, FALSE, FALSE, TRUE)], .lshape4, 
                            earg = .eshape4)
        rpairbeta(n = nsim*length(shape1), dimData = 3, par = c(shape1,shape2,shape3,shape4))[,c(1,2)]
      }, list(.lshape1 = lshape1, .lshape2 = lshape2, .lshape3 = lshape3, .lshape4 = lshape4, 
              .eshape1 = eshape1, .eshape2 = eshape2, .eshape3 = eshape3, .eshape4 = eshape4))), 
      deriv = eval(substitute(expression({
        shape1 <- eta2theta(eta[, c(TRUE, FALSE, FALSE, FALSE)], .lshape1, 
                            earg = .eshape1)
        shape2 <- eta2theta(eta[, c(FALSE, TRUE, FALSE, FALSE)], .lshape2, 
                            earg = .eshape2)
        shape3 <- eta2theta(eta[, c(FALSE, FALSE, TRUE, FALSE)], .lshape3, 
                            earg = .eshape3)
        shape4 <- eta2theta(eta[, c(FALSE, FALSE, FALSE, TRUE)], .lshape4, 
                            earg = .eshape4)
        dshape1.deta <- dtheta.deta(shape1, link = .lshape1, 
                                    earg = .eshape1)
        dshape2.deta <- dtheta.deta(shape2, link = .lshape2, 
                                    earg = .eshape2)
        dshape3.deta <- dtheta.deta(shape3, link = .lshape3, 
                                    earg = .eshape3)
        dshape4.deta <- dtheta.deta(shape4, link = .lshape4, 
                                    earg = .eshape4)
        d2shape1.deta2 <- d2theta.deta2(shape1, link = .lshape1, 
                                        earg = .eshape1)
        d2shape2.deta2 <- d2theta.deta2(shape2, link = .lshape2, 
                                        earg = .eshape2)
        d2shape3.deta2 <- d2theta.deta2(shape3, link = .lshape3, 
                                        earg = .eshape3)
        d2shape4.deta2 <- d2theta.deta2(shape4, link = .lshape4, 
                                        earg = .eshape4)
        dl.dshape1 <- pBeta.3d.d1.alpha(y[,1],y[,2],shape1,shape2,shape3,shape4)
        dl.dshape2 <- pBeta.3d.d1.beta1(y[,1],y[,2],shape1,shape2,shape3,shape4)
        dl.dshape3 <- pBeta.3d.d1.beta2(y[,1],y[,2],shape1,shape2,shape3,shape4)
        dl.dshape4 <- pBeta.3d.d1.beta3(y[,1],y[,2],shape1,shape2,shape3,shape4)
        
        dl.deta <- c(w) * cbind(dl.dshape1 * dshape1.deta, 
                                dl.dshape2 * dshape2.deta,
                                dl.dshape3 * dshape3.deta,
                                dl.dshape4 * dshape4.deta)
        dl.deta[, interleave.VGAM(M, M1 = M1)]
      }), list(.lshape1 = lshape1, .lshape2 = lshape2, .lshape3 = lshape3, .lshape4 = lshape4,  
               .eshape1 = eshape1, .eshape2 = eshape2,.eshape3 = eshape3,.eshape4 = eshape4))),
      weight = eval(substitute(expression({
        ned2l.dshape11 <- (-pBeta.3d.d2.alpha(y[,1],y[,2],shape1,shape2,shape3,shape4))
        ned2l.dshape22 <- (-pBeta.3d.d2.beta1(y[,1],y[,2],shape1,shape2,shape3,shape4))
        ned2l.dshape33 <- (-pBeta.3d.d2.beta2(y[,1],y[,2],shape1,shape2,shape3,shape4))
        ned2l.dshape44 <- (-pBeta.3d.d2.beta3(y[,1],y[,2],shape1,shape2,shape3,shape4))
        ned2l.dshape12 <- (-pBeta.3d.d2.alpha.beta1(y[,1],y[,2],shape1,shape2,shape3,shape4))
        ned2l.dshape23 <- (-pBeta.3d.d2.beta1.beta2(y[,1],y[,2],shape1,shape2,shape3,shape4))
        ned2l.dshape34 <- (-pBeta.3d.d2.beta2.beta3(y[,1],y[,2],shape1,shape2,shape3,shape4))
        ned2l.dshape13 <- (-pBeta.3d.d2.alpha.beta2(y[,1],y[,2],shape1,shape2,shape3,shape4))
        ned2l.dshape24 <- (-pBeta.3d.d2.beta1.beta3(y[,1],y[,2],shape1,shape2,shape3,shape4))
        ned2l.dshape14 <- (-pBeta.3d.d2.alpha.beta3(y[,1],y[,2],shape1,shape2,shape3,shape4))
        
        wz <- array(c(c(w) * ((ned2l.dshape11 * dshape1.deta^2)-
                                (pBeta.3d.d1.alpha(y[,1],y[,2],shape1,shape2,shape3,shape4)*d2shape1.deta2)),
                      c(w) * ((ned2l.dshape22 * dshape2.deta^2)-
                                (pBeta.3d.d1.beta1(y[,1],y[,2],shape1,shape2,shape3,shape4)*d2shape2.deta2)),
                      c(w) * ((ned2l.dshape33 * dshape3.deta^2)-
                                (pBeta.3d.d1.beta2(y[,1],y[,2],shape1,shape2,shape3,shape4)*d2shape3.deta2)),
                      c(w) * ((ned2l.dshape44 * dshape4.deta^2)-
                                (pBeta.3d.d1.beta3(y[,1],y[,2],shape1,shape2,shape3,shape4)*d2shape4.deta2)),
                      c(w) * ned2l.dshape12 * dshape1.deta * dshape2.deta,
                      c(w) * ned2l.dshape23 * dshape2.deta * dshape3.deta,
                      c(w) * ned2l.dshape34 * dshape3.deta * dshape4.deta,
                      c(w) * ned2l.dshape13 * dshape1.deta * dshape3.deta,
                      c(w) * ned2l.dshape24 * dshape2.deta * dshape4.deta,
                      c(w) * ned2l.dshape14 * dshape1.deta * dshape4.deta),
                    dim = c(n, M/M1, 10))
        wz <- arwz2wz(wz, M = M, M1 = M1)
        wz
      }), list(.lshape1 = lshape1, .lshape2 = lshape2, .lshape3 = lshape3, .lshape4 = lshape4,  
               .eshape1 = eshape1, .eshape2 = eshape2, .eshape3 = eshape3, .eshape4 = eshape4,
               .tol12 = tol12))))
}

########################################################################################################
### Examples of use
########################################################################################################
if(!require(VGAM)) install.packages("VGAM") 
if(!require(mev)) install.packages("BMAmevt") 

library("VGAM")
library("BMAmevt")
#needed for dpairbeta

B <- 300
t <- seq(0.8,3.3,length=B)
alpha <- exp(t)
beta1 <- t
beta2 <- log(t+1)
beta3 <- log(t+2)

z <- matrix(0,ncol=3,nrow=B)
set.seed(22)
for (i in 1:B)
{
  z[i,] <- rpairbeta(n = 1, dimData = 3, par = c(alpha[i],beta1[i],beta2[i],beta3[i]))
}

db <- data.frame("z1"=z[,1],"z2"=z[,2],"z3"=z[,3],"t"=t)

t             <- db$t
z.data        <- data.frame("y1"=db$z1,"y2"=db$z2, "t"=t)

bla         <- maxLikelihood.pb(cbind(db$z1,db$z2,db$z3),model="pb")
alpha.init  <- bla$par[1]
beta1.init  <- bla$par[2]
beta2.init  <- bla$par[3]
beta3.init  <- bla$par[4]

#fit a VGAM to the bivariate random vector (y1,y2) with "intercept-only" parameters 
#(should be equal to the MLE)
pbeta.vgam.0 <- vgam(cbind(y1,y2)~1,family=pBeta3d.vgam(ishape1 = alpha.init, 
                                                        ishape2 = beta1.init,
                                                        ishape3 = beta2.init,
                                                        ishape4 = beta3.init), 
                    data=z.data, trace=TRUE)
head(loge(pbeta.vgam.0@predictors,inverse = TRUE))

# fit a VGAM with no constraints
pbeta.vgam.1 <- vgam(cbind(y1,y2)~sm.ps(t),family=pBeta3d.vgam(ishape1 = alpha.init, 
                                                             ishape2 = beta1.init,
                                                             ishape3 = beta2.init,
                                                             ishape4 = beta3.init),
                               data=z.data,trace=TRUE)

########### Asymptotic CI

quantile.fun <- function(x) quantile(x, c((1 - conf)/2, 1 - (1 - conf)/2))
myCI <- function(x) t(apply(x, 2, quantile.fun))
conf <- .95

par(mfrow=c(2,2))

k    <- 1
Xp <- model.matrix(pbeta.vgam.1,linpred.index=k,type="lm")
b  <- coef(pbeta.vgam.1)[seq(k,40,by=4)]
Vp <- vcov(pbeta.vgam.1)[seq(k,40,by=4),seq(k,40,by=4)]
br <- mvrnorm(10000, b, Vp)
calib <- br %*% t(Xp)

calib.CI.b <- myCI(VGAM::loge(calib,inverse=TRUE))

par(mar=c(3,3,1.5,0.5),mgp=c(1.6,0.5,0),font.main=1,cex=1.5,cex.main=1)
plot(t,alpha,type="l",ylab=expression(alpha),col=2,ylim=c(range(calib.CI.b[,1])[1],range(calib.CI.b[,2])[2]));
lines(t,VGAM::loge(pbeta.vgam.1@predictors[,k],inverse=TRUE))
lines(t,calib.CI.b[,1],lty=2,lwd=2)
lines(t,calib.CI.b[,2],lty=2,lwd=2)

k    <- 2
Xp <- model.matrix(pbeta.vgam.1,linpred.index=k,type="lm")
b  <- coef(pbeta.vgam.1)[seq(k,40,by=4)]
Vp <- vcov(pbeta.vgam.1)[seq(k,40,by=4),seq(k,40,by=4)]
br <- mvrnorm(10000, b, Vp)
calib <- br %*% t(Xp)

calib.CI.b <- myCI(VGAM::loge(calib,inverse=TRUE))

par(mar=c(3,3,1.5,0.5),mgp=c(1.6,0.5,0),font.main=1,cex=1.5,cex.main=1)
plot(t,beta1,type="l",ylab=expression(beta[1]),col=2,ylim=c(range(calib.CI.b[,1])[1],range(calib.CI.b[,2])[2]));
lines(t,VGAM::loge(pbeta.vgam.1@predictors[,k],inverse=TRUE))
lines(t,calib.CI.b[,1],lty=2,lwd=2)
lines(t,calib.CI.b[,2],lty=2,lwd=2)

k    <- 3
Xp <- model.matrix(pbeta.vgam.1,linpred.index=k,type="lm")
b  <- coef(pbeta.vgam.1)[seq(k,40,by=4)]
Vp <- vcov(pbeta.vgam.1)[seq(k,40,by=4),seq(k,40,by=4)]
br <- mvrnorm(10000, b, Vp)
calib <- br %*% t(Xp)

calib.CI.b <- myCI(VGAM::loge(calib,inverse=TRUE))

par(mar=c(3,3,1.5,0.5),mgp=c(1.6,0.5,0),font.main=1,cex=1.5,cex.main=1)
plot(t,beta2,type="l",ylab=expression(beta[2]),col=2,ylim=c(range(calib.CI.b[,1])[1],range(calib.CI.b[,2])[2]));
lines(t,VGAM::loge(pbeta.vgam.1@predictors[,k],inverse=TRUE))
lines(t,calib.CI.b[,1],lty=2,lwd=2)
lines(t,calib.CI.b[,2],lty=2,lwd=2)

k    <- 4
Xp <- model.matrix(pbeta.vgam.1,linpred.index=k,type="lm")
b  <- coef(pbeta.vgam.1)[seq(k,40,by=4)]
Vp <- vcov(pbeta.vgam.1)[seq(k,40,by=4),seq(k,40,by=4)]
br <- mvrnorm(10000, b, Vp)
calib <- br %*% t(Xp)

calib.CI.b <- myCI(VGAM::loge(calib,inverse=TRUE))

par(mar=c(3,3,1.5,0.5),mgp=c(1.6,0.5,0),font.main=1,cex=1.5,cex.main=1)
plot(t,beta3,type="l",ylab=expression(beta[3]),col=2,ylim=c(range(calib.CI.b[,1])[1],range(calib.CI.b[,2])[2]));
lines(t,VGAM::loge(pbeta.vgam.1@predictors[,k],inverse=TRUE))
lines(t,calib.CI.b[,1],lty=2,lwd=2)
lines(t,calib.CI.b[,2],lty=2,lwd=2)

# fit a VGAM to multiple responses (same response (y1,y2) here) with a paralellism assumption
pbeta.vgam.2 <- vgam(cbind(y1,y2,y1,y2)~sm.ps(t),family=pBeta3d.vgam(ishape1 = alpha.init, 
                                                               ishape2 = beta1.init,
                                                               ishape3 = beta2.init,
                                                               ishape4 = beta3.init,
                                                               parallel = TRUE),
                     data=z.data,trace=TRUE)

#with no surprises...
all.equal(loge(pbeta.vgam.2@predictors,inverse = TRUE)[,1],loge(pbeta.vgam.2@predictors,inverse = TRUE)[,5])
all.equal(loge(pbeta.vgam.2@predictors,inverse = TRUE)[,2],loge(pbeta.vgam.2@predictors,inverse = TRUE)[,6])
all.equal(loge(pbeta.vgam.2@predictors,inverse = TRUE)[,3],loge(pbeta.vgam.2@predictors,inverse = TRUE)[,7])
all.equal(loge(pbeta.vgam.2@predictors,inverse = TRUE)[,4],loge(pbeta.vgam.2@predictors,inverse = TRUE)[,8])

