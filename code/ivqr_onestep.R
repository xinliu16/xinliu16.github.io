

source("gmmq.R")


ivqr.onestep <- function(Y, X.excl, Z.excl, D, tau, h=0, dB,
                         structure=c('iid','ts'), LRV.kernel=c('QS','Bartlett','uniform'),
                         Lambda=function(y, x, b){y-x%*%b }, Lambda.derivative=function(y,x,b) -x,
                         Lambda.gmmq=function(y,x,b) {y[,1]-cbind(y[,-1],x)%*%b}, 
                         Lambda.derivative.gmmq=function(y,x,b) {-cbind(y[,-1],x)}, 
                         b.init=0 ) {
  X <- cbind(D, X.excl)
  Z <- cbind(Z.excl, X.excl)
  Y.gmmq <- cbind(Y,D)

  if (b.init!=0 && length(b.init)!=dB) stop(sprintf("Length of b.init must equal dB, which is %d", dB))
  
  if ((h==0) && (b.init==0)  ) {ivqr.est <- gmmq(tau=tau, Y=Y.gmmq, X=X.excl, Z.excl=Z.excl, dB=dB, 
                                                 Lambda=Lambda.gmmq, 
                                                 Lambda.derivative=Lambda.derivative.gmmq, 
                                                 h=0, VERBOSE=FALSE, RETURN.Z=FALSE )     }
  
  if ((h!=0) && (b.init==0)  ) {ivqr.est <- gmmq(tau=tau, Y=Y.gmmq, X=X.excl, Z.excl=Z.excl, dB=dB, 
                                                 Lambda=Lambda.gmmq, 
                                                 Lambda.derivative=Lambda.derivative.gmmq, 
                                                 h=h, VERBOSE=FALSE, RETURN.Z=FALSE )  }
  
  if ((h==0) && (b.init!=0)  ) {ivqr.est <- gmmq(tau=tau, Y=Y.gmmq, X=X.excl, Z.excl=Z.excl, dB=dB, 
                                                 Lambda=Lambda.gmmq, 
                                                 Lambda.derivative=Lambda.derivative.gmmq, 
                                                 h=0, VERBOSE=FALSE, RETURN.Z=FALSE, b.init=b.init)  }
  
  if ((h!=0) && (b.init!=0)  ) {ivqr.est <- gmmq(tau=tau, Y=Y.gmmq, X=X.excl, Z.excl=Z.excl, dB=dB, 
                                                 Lambda=Lambda.gmmq, 
                                                 Lambda.derivative=Lambda.derivative.gmmq, 
                                                 h=h, VERBOSE=FALSE, RETURN.Z=FALSE, b.init=b.init)   }
  b.gmmq <- ivqr.est$b
  h.gmmq <- ivqr.est$h
        
  Itilde <- Itilde.KS17
  Itilde.deriv <- Itilde.deriv.KS17
  
        # Add bandwidth as second argument to G() and G'() functions.
        Gfn <- function(v,h) {      Itilde(v/h)       }
        Gpfn <- function(v,h) {     Itilde.deriv(v/h) }
        
        
        n <- length(Y) #number of observations, a.k.a. sample size
        dZ <- ncol(Z) #number of instruments
        dX <- ncol(X) #number of regressors
        
        L <- matrix(Lambda.gmmq(y=Y.gmmq,x=X.excl,b=ivqr.est$b), ncol=1)
        Lp <- Lambda.derivative.gmmq(y=Y.gmmq, x=X.excl, b=ivqr.est$b)
        
 
        gni <- Z*repmat((Gfn(-L,ivqr.est$h)-tau),1,dZ)
        g.bar <- colMeans(gni) 
       
        h.use <- max(sort(abs(L))[floor(n^(4/5))], ivqr.est$h)
        jac <- t(Z)%*%(Lp*repmat(Gpfn(-L,h.use),1,dB))/(-h.use) 
        

        if (structure=='iid') { Omega.bar <- t(gni)%*%gni/n}
        if (structure=='ts') { Omega.bar <- LRV.est.fn(tau=tau,Y=Y.gmmq,X=X.excl,Z=Z,Lambda=Lambda.gmmq,beta.hat=ivqr.est$b,Itilde=Itilde,h=ivqr.est$h,structure=structure,LRV.kernel=LRV.kernel)}
        
        b_onestep <- ivqr.est$b - solve(t(jac)%*%solve(Omega.bar)%*%jac) %*% t(jac)%*%solve(Omega.bar)%*%g.bar
        
        return(list(b.onestep=b_onestep,h=ivqr.est$h, b.gmmq=b.gmmq, h.gmmq=h.gmmq))
        
}




#EOF
