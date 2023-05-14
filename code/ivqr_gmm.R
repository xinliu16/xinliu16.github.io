
IDENTITY.FLAG <- FALSE
ZZ.FLAG <- FALSE


source("ivqr_onestep.R")
source("gmmq.R")

library("GenSA")

ivqr.gmm <- function(Y, X.excl, Z.excl, D, tau, h=0, dB, max.time, upper, lower, weight.mtx, 
                     structure=c('iid','ts'), LRV.kernel=c('QS','Bartlett','uniform'), 
                     Lambda=function(y, x, b){y-x%*%b }, Lambda.derivative=function(y,x,b){-x}, 
                     Lambda.gmmq=function(y,x,b) {y[,1]-cbind(y[,-1],x)%*%b}, 
                     Lambda.derivative.gmmq=function(y,x,b) {-cbind(y[,-1],x)}, 
                     b.init=0 ) {
    X <- cbind(D, X.excl)
    Z <- cbind(Z.excl, X.excl)
    Y.gmmq <- cbind(Y,D)
    
    n <- length(Y)
    dZ <- ncol(Z) #number of instruments
    dX <- ncol(X) #number of regressors

    Itilde <- Itilde.KS17
    Itilde.deriv <- Itilde.deriv.KS17
    
    # Add bandwidth as second argument to G() and G'() functions.
    Gfn <- function(v,h) {      Itilde.KS17(v/h)    }
    Gpfn <- function(v,h) {      Itilde.deriv.KS17(v/h)    }
    
    # 
    ret <- ivqr.onestep(Y=Y, X.excl=X.excl, Z.excl=Z.excl, D=D, tau=tau, h=h, dB=dB,
                      structure=structure, LRV.kernel=LRV.kernel, 
                      Lambda=Lambda, Lambda.derivative=Lambda.derivative,
                      Lambda.gmmq=Lambda.gmmq, Lambda.derivative.gmmq=Lambda.derivative.gmmq, 
                      b.init=b.init ) 
    
    b.initial <- ret$b.onestep
    h.initial <- ret$h 
    
    if(length(b.initial)!=length(lower)){stop("lower bound doesn't have right dimension.")}
    if(length(b.initial)!=length(upper)){stop("upper bound doesn't have right dimension.")}
    
    
    LRV.hat <- LRV.est.fn(tau=tau, Y=Y.gmmq, X=X.excl, Z=Z, 
                          Lambda=Lambda.gmmq, beta.hat=b.initial, Itilde=Itilde.KS17, h=h.initial, 
                          structure=structure, LRV.kernel=LRV.kernel) 
    if (IDENTITY.FLAG){W.hat <-diag(dZ) } else{
    if (ZZ.FLAG) {W.hat <- solve(t(Z)%*%Z/n)
    } else{   
    if (missing(weight.mtx)) {   W.hat <- solve(LRV.hat) }
                         else{   W.hat <- weight.mtx}
                        }
                                              }
      ivqr.obj <- function(b, h, weight.mtx) {
          L <- matrix(Lambda(y=Y, x=X, b=b), ncol=1)

          gni <- Z*repmat((Gfn(-L,h)-tau),1,dZ) 
          g.bar <- as.matrix(colMeans(Z*repmat((Gfn(-L,h)-tau),1,dZ))) 

          obj.fn <- t(g.bar)%*% weight.mtx %*%g.bar
          return(obj.fn)
      }

      obj.fn <- function(b) n*ivqr.obj(b=b,  h=h.initial, weight.mtx=W.hat)
 
     
       M.hat <- function(b) {
       L <- matrix(Lambda(y=Y, x=X, b=b), ncol=1)
       gni <- Z*repmat((Gfn(-L,h.initial)-tau),1,dZ) 
       g.bar <- as.matrix(colMeans(Z*repmat((Gfn(-L,h.initial)-tau),1,dZ))) 
       return(g.bar)        }
     
     
     ivqr.gmm1 <- GenSA(par=b.initial, fn=obj.fn,
                        lower=lower, upper=upper,
                  control=list(max.time=max.time))
     b <- ivqr.gmm1$par
     
     J.stat <- obj.fn(b)
 
     return(list(b=b, obj.fn=obj.fn, h=h.initial, W.hat=W.hat, M.hat=M.hat, J.stat=J.stat))
  }


#EOF