# Implementation of smoothed IV quantile regression from
# Kaplan and Sun (2016), "Smoothed estimating equations for instrumental variables quantile regression" (Econometric Theory)
# Original Matlab code by Dave Kaplan
# Ported to R by Lonnie Hofmann
# Questions? Comments? email Dave: kaplandm@missouri.edu

ivqr_see <- function(p,Y,X,Z,h_in) {

  # IVQR_SEE  Quantile regression with smoothed estimating equation, with or without instrumental variables.

  # ivqr_see(p,Y,X) returns a list where b is the estimated parameter vector for the quantile regression (QR) model Y=X*B+U where P(U<0|X)=p, using the smoothed estimating equation (SEE) estimator of Kaplan and Sun (2016). Y is a vector, and X is a matrix with n rows and d columns, where n is the sample size and d is the number of parameters.  For example, p=0.5 is median regression, a.k.a. least absolute deviation (LAD) regression.  The returned list also shows the computed plug-in bandwidth value (hhat) and the bandwidth actually used (h); if hhat is too small, a larger h is used to compute b.
  
  # ivqr_see(p,Y,X,Z) uses the instrumental variables quantile regression (IVQR) model Y=X*B+U where P(U<0|Z)=p, using the smoothed estimating equation (SEE) estimator of Kaplan and Sun (2016).  Y and X are the same as above, and Z is a matrix with n rows and at least d columns.  The columns of Z should include both exogenous regressors (columns from X) and additional instruments that replace the endogenous regressor columns of X.  Note that ivqr_see(p,Y,X,X) = ivqr_see(p,Y,X).  Return values are the same as before.
  
  # ivqr_see(p,Y,X,Z,h_in) additionally sets the bandwidth to h_in.  If h_in is too small, then in the returned list, hhat is the specified bandwidth (h_in) while h is the bandwidth actually used to estimate b; so h>=hhat.  (If h_in is not given, the data-dependent plug-in bandwidth selection procedure described in Kaplan and Sun (2016) is used.)

  # Check arguments: errors, warnings, defaults
r <- 4
if (!require(pracma)) {
		stop('Please install the R package pracma to run this code.')
}
if (missing(p)) {
	stop('Need input for p.')
}
if (missing(Y)) {
	stop('Need input for Y.')
}
if (missing(X)) {
	stop('Need input for X.')
}

if (!is.numeric(p) || !is.numeric(Y) || !is.numeric(X) || !missing(Z) && !is.numeric(Z) || !missing(h_in) && !is.numeric(h_in)) {
	stop('All inputs must be numeric types.')
}

if (!is.matrix(Y)) Y <- matrix(Y,ncol=1)
if (ncol(Y) > nrow(Y)) {
	stop('Y must be a matrix (column vector) with dim(Y) <- c(n,1), where n is the number of observations.')
} else if (ncol(X) > nrow(X)) {
	stop('X must have dim(X) <- c(n,d) where n is the number of observations, d is the number of regressors, and n>d.')
} else if (nrow(X)!=nrow(Y)) {
	stop('Must have nrow(X)=nrow(Y)=n, the number of observations.')
} else if (!missing(Z) && ncol(Z)>nrow(Z)) {
	stop('Z must have dim(Z) <- c(n,d2), where n is the number of observations, d2 is the number of instruments, and n>d2.')
} else if (!missing(Z) && nrow(Z)!=nrow(Y)) {
	stop('Z must have nrow(Z)=nrow(Y)=n, the number of observations.')
} else if (!missing(Z) && ncol(Z)<ncol(X)) {
	stop('Must have ncol(Z)>=ncol(X), where Z contains exogenous regressors from X, plus one or more instruments for each endogenous regressor in X.')
}

if (p<=0 || p>=1) {
	stop('Quantile p (first argument) must be strictly between 0 and 1.')
} else if (!missing(h_in) && h_in<0) {
	stop('Bandwidth must be nonnegative.')
}

  # Defaults, initialization
if (missing(Z)) {
	Z <- X
}
n <- length(Y) #number of observations, a.k.a. sample size
d <- ncol(X) #number of regressors
if (ncol(Z)>ncol(X)) {
	Z <- Z%*%solve(t(Z)%*%Z)%*%t(Z)%*%X #Z%*%(solve(Z,X)) 
} #transform, per Kaplan and Sun (2012, Eqn. 1)
  #Z*inv(Z'*Z)*Z'*X in Matlab
  
h_flag<-FALSE
if (!missing(h_in)) {
	h_flag <- TRUE
}

dgev.ls.fn <- function(x,xi,mu=0,sig=1) {
  if (sig<=0) {
    return(rep(0,length(x)))
  } else if (xi!=0) {
		tmp <- (1 + (xi * (x - mu))/sig)
		return(ifelse(tmp>0,(tmp^(-1/xi - 1)*exp(-tmp^(-1/xi)))/sig,0))
	} else {
		tmp <- exp(-(x-mu)/sig)
		return((tmp^(xi+1)*exp(-tmp))/sig)
	}
} #GEV pdf

  # Setup for smoothing function
  # For different possible parametric distributions, generate the objects needed for the bandwidth as functions of distribution parameters.
if (!h_flag) {
	fN0_fn <- function(nmu,sig) {
		dnorm(0,nmu,sig)
	}
	fN0r1_fn <- function(nmu,sig) {
			dnorm(0,nmu,sig)%*%(nmu/sig^4)%*%(-3+(nmu%*%nmu)/sig^2)
	} 
	ft0_fn <- function(v,tmu,tsig) {
		if (v>100) {
			dnorm(0,tmu,tsig)
		} else {
		pi^(-1/2)%*%((1/sqrt(v)*gamma(v*(1/2)+1/2)*((tmu^2*1/tsig^2)/v+1)^(v*(-1/2)-1/2))/(tsig*gamma(v*(1/2))))
	    }
	}
	ft0r1_fn <- function(v,tmu,tsig) {
		if (v>100) {
			fN0r1_fn(tmu,tsig)
		} else {
			pi^(-1/2)%*%((tmu*1/tsig^5*1/v^(5/2)*gamma(v*(1/2)+1/2)*(v*(1/2)+1/2)*(v*(1/2)+3/2)*((tmu^2*1/tsig^2)/v+1)^(v*(-1/2)-5/2)*-12)/gamma(v*(1/2))+(tmu^3*1/tsig^7*1/v^(7/2)*gamma(v*(1/2)+1/2)*(v*(1/2)+1/2)*(v*(1/2)+3/2)*(v*(1/2)+5/2)*((tmu^2*1/tsig^2)/v+1)^(v*(-1/2)-7/2)*8)/gamma(v*(1/2)))
	    }
	}
	fgam0_fn <- function(k,theta,Ueval) {
		dgamma(Ueval,k,scale=theta)
	}
	fgam0r1_fn <- function(k,theta,Ueval) {
			-(Ueval^(k-1)*1/theta^3*theta^(-k)*exp(-Ueval/theta))/gamma(k)+(Ueval^(k-2)*1/theta^2*theta^(-k)*exp(-Ueval/theta)*(k-1)*3)/gamma(k)-(Ueval^(k-3)*theta^(-k)*exp(-Ueval/theta)*(k-1)*(k-2)*3)/(theta*gamma(k))+(Ueval^(k-4)*theta^(-k)*exp(-Ueval/theta)*(k-1)*(k-2)*(k-3))/gamma(k)
	}
	fGEV0_fn <- function(xi,sig,gmu) {
		dgev.ls.fn(0,xi,mu,sig)
		} 
	fGEV0r1_fn <- function(xi,sig,gmu) {
		1/sig^4*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi+1)*(-(gmu*xi)/sig+1)^(-3/xi-3)-1/sig^4*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^xi*(-(gmu*xi)/sig+1)^(-3/xi-3)*(xi+1)-1/sig^4*xi*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi+1)*(1/xi+1)*(-(gmu*xi)/sig+1)^(-2/xi-3)*2-1/sig^4*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^xi*(-(gmu*xi)/sig+1)^(-1/xi-1)*(-(gmu*xi)/sig+1)^(-2/xi-2)*(xi+1)*2+1/sig^4*xi*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi-1)*(-(gmu*xi)/sig+1)^(-3/xi-3)*(xi+1)+1/sig^4*xi*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^xi*(-(gmu*xi)/sig+1)^(-2/xi-3)*(2/xi+2)*(xi+1)*2+1/sig^4*xi*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi-1)*(-(gmu*xi)/sig+1)^(-1/xi-1)*(-(gmu*xi)/sig+1)^(-2/xi-2)*(xi+1)*2-1/sig^4*xi^2*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi-1)*(1/xi+1)*(-(gmu*xi)/sig+1)^(-2/xi-3)*(xi+1)*2-1/sig^4*xi*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi+1)*(1/xi+1)*(-(gmu*xi)/sig+1)^(-1/xi-1)*(-(gmu*xi)/sig+1)^(-1/xi-2)-1/sig^4*xi*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi-2)*(-(gmu*xi)/sig+1)^(-3/xi-3)*(xi-1)*(xi+1)+1/sig^4*xi^2*exp(-(-(gmu*xi)/sig+1)^(-1/xi))*((-(gmu*xi)/sig+1)^(-1/xi))^(xi+1)*(1/xi+1)*(1/xi+2)*(-(gmu*xi)/sig+1.0)^(-1.0/xi-3.0)+1/sig^4*xi*exp(-(-(gmu*xi)/sig+1.0)^(-1.0/xi))*((-(gmu*xi)/sig+1.0)^(-1.0/xi))^xi*(1.0/xi+1.0)*(-(gmu*xi)/sig+1.0)^(-1.0/xi-1.0)*(-(gmu*xi)/sig+1.0)^(-1.0/xi-2.0)*(xi+1.0)*2.0-1/sig^4*xi^2*exp(-(-(gmu*xi)/sig+1.0)^(-1.0/xi))*((-(gmu*xi)/sig+1.0)^(-1.0/xi))^xi*(1.0/xi+1.0)*(1.0/xi+2.0)*(-(gmu*xi)/sig+1.0)^(-1.0/xi-3.0)*(xi+1.0)-1/sig^4*xi^2*exp(-(-(gmu*xi)/sig+1.0)^(-1.0/xi))*((-(gmu*xi)/sig+1.0)^(-1.0/xi))^(xi-1.0)*(1/xi+1)*(-(gmu*xi)/sig+1)^(-1/xi-1)*(-(gmu*xi)/sig+1)^(-1/xi-2)*(xi+1)
     }
	hopt_estfn1_t <- function(Gsqiv,CKv,mu,sig,v,d,n) {
		((factorial(r)^2*(1-Gsqiv)*ft0_fn(v,mu,sig)*d)/(2*r*CKv^2*ft0r1_fn(v,mu,sig)^2*n))^(1/(2*r-1))
	}
	hopt_estfn1_N <- function(Gsqiv,CKv,mu,sig,d,n) {
		((factorial(r)^2*(1-Gsqiv)*fN0_fn(mu,sig)*d)/(2*r*CKv^2*fN0r1_fn(mu,sig)^2*n))^(1/(2*r-1))
	}
	hopt_estfn1_gam <- function(Gsqiv,CKv,k,theta,Ueval,d,n) {
		((factorial(r)^2*(1-Gsqiv)*fgam0_fn(k,theta,Ueval)*d)/(2*r*CKv^2*fgam0r1_fn(k,theta,Ueval)^2*n))^(1/(2*r-1))
	}
	hopt_estfn1_GEV <- function(Gsqiv,CKv,xi,sig,mu,d,n) {
		((factorial(r)^2*(1-Gsqiv)*fGEV0_fn(xi,sig,mu)*d)/(2*r*CKv^2*fGEV0r1_fn(xi,sig,mu)^2*n))^(1/(2*r-1))
	}
}
  # Set smoothing fn G() and its derviative G'()
G_symfn <- function(u) {
  ifelse(u >= 1, 1, ifelse(u > -1, 1/2 + (105/64)*(u-(5/3)*u^3+(7/5)*u^5 -(3/7)*u^7), 0))
} 
Gp_symfn <- function(u) {
  ifelse(u > -1 & u < 1, (105/64)*(1-5*u^2+7*u^4-3*u^6), 0)
} 
Gsqint_val <- 394/429
CK_val <- -1/33

  # Add bandwidth as second argument to G() and G'() functions.
Gfn <- function(v,h) {
	G_symfn(v/h)
}
Gpfn <- function(v,h) {
	Gp_symfn(v/h)
}

  # Plug-in bandwidth as function of distributional parameters, for different possible parametric distributions.
if (!h_flag) {
	hopt_estfn_t <- function(mu,sig,v,d,n) {
		hopt_estfn1_t(Gsqint_val,CK_val,mu,sig,v,d,n)
	}
	hopt_estfn_N <- function(mu,sig,d,n) {
		hopt_estfn1_N(Gsqint_val,CK_val,mu,sig,d,n)
	}
	hopt_estfn_gam <- function(k,theta,p,d,n) {
		hopt_estfn1_gam(Gsqint_val,CK_val,k,theta,qgamma(p,k,scale=theta),d,n)
	}
	hopt_estfn_GEV <- function(xi,sig,mu,d,n) {
		hopt_estfn1_GEV(Gsqint_val,CK_val,xi,sig,mu,d,n)
	}
}

  # For solver (to find solution to estimating equations)
objfn <- function(b,h) {
	colSums(Z*repmat((Gfn(X%*%b-Y,h)-p),1,d))
} # Objective Function 
jacfn <- function(b,h) {
	t(Z)%*%(X*repmat(Gpfn(X%*%b-Y,h),1,d))/h
} # Jacobian Function
objfn1 <- function (b) {objfn(b,htmp)}
jacfn1 <- function (b) {jacfn(b,htmp)}

  # Initial Estimation
b_win <- matrix(0,d,1)
found_flag <- FALSE
h_init <- 0.3*(d*factorial(r)^2/(2*n*r))^(1/(2*r-1))
htmp <- h_init
h_pilot <- h_init
hfac <- 4
tmp_b0 <- matrix(0,d,1) #where to start fsolve
MAXITER <- 400
while (htmp>=h_init*hfac^(-0)) {
	btmp <- tryCatch(newtonsys(Ffun=objfn1,x0=tmp_b0,Jfun=jacfn1,maxiter=MAXITER),
				warning=function(w) NA,
				error=function(w) NA
				)
	exitOK1 <- !(is.na(btmp)[1] || (btmp$niter==MAXITER)[1])
	if (exitOK1) {
		h_pilot <- htmp
	}
	if (found_flag && !exitOK1) {
		break
	} else if (exitOK1[1] && (htmp==h_init)[1]) {
		found_flag <- TRUE
		b_win <- btmp$zero
		f_win <- objfn1(btmp$zero)
		htmp <- htmp/hfac
	} else if (exitOK1[1] && !found_flag[1] && (htmp>h_init)[1]) {
		b_win <- btmp$zero
		f_win <- objfn1(btmp$zero)
		break
	} else if (exitOK1[1] && found_flag[1]) {
		b_win <- btmp$zero
		f_win <- objfn1(btmp$zero)
		htmp <- htmp/hfac
	} else if (!exitOK1[1] && !found_flag[1]) {
		htmp <- htmp*hfac
	} else {
		stop('Missed case in b_pilot initialization.')
	}
}
b_pilot <- b_win

# If h_in bandwidth not provided, use plug-in procedure
dt.ls.fn <- function(x,m,sigma,df) dt(x=(x-m)/sigma,df=df)/sigma #3-parameter student's t

if (!h_flag) {
	Uhat <- Y-X%*%b_pilot
	Uhat <- Uhat-quantile(Uhat,p) # Shift Uhat dist'n such that p-quantile=0; o/w bias
	if (!require(MASS)) {
		stop('Please install the R package MASS to run this code.')
	}
	t.MLE.est <- tryCatch(fitdistr(x=Uhat,densfun=dt.ls.fn,start=list(m = mean(Uhat), sigma = sd(Uhat), df = 10),lower=c(-Inf,0,1), control=list(maxit=200)),
						warning=function(w) NA,
						error=function(w) NA
						) # If warning/error, don't use (set to Inf)
	if (is.na(t.MLE.est[1])) {
		tnlogL <- Inf
		hoptests_t <- Inf
	} else {
		tnlogL <- -t.MLE.est$loglik
		v <- max(1,round(t.MLE.est$estimate["df"]))
		mu <- t.MLE.est$estimate["m"]
		sig <- t.MLE.est$estimate["sigma"]
		hoptests_t <- hopt_estfn_t((mu+(0==mu)*.01),sig,v,d,n)
	}
	tmp.sig <- sd(Uhat)*sqrt(6)/pi
	tmp.mu <- mean(Uhat) - tmp.sig*(-digamma(1)) #approx 0.5772... (Euler's Constant)
	gev.MLE.est <- tryCatch(suppressWarnings(fitdistr(x=Uhat,densfun=dgev.ls.fn,list(xi=0, mu=tmp.mu, sig=tmp.sig), method='BFGS', control=list(maxit=200))), 
	                        #warning=function(w) NA,
	                        error=function(w) {
	                          tryCatch(suppressWarnings(fitdistr(x=Uhat,densfun=dgev.ls.fn, list(xi=0, mu=tmp.mu, sig=tmp.sig), lower=c(-Inf,-Inf,0), control=list(maxit=200))), 
	                                   error=function(w2)NA)
	                        }
	)
	if(is.na(gev.MLE.est[1]) || any(!is.finite(gev.MLE.est$estimate))) {
		gevnlogL <- Inf
		hoptests_gev <- Inf
	} else {
		gevnlogL <- -gev.MLE.est$loglik
		xi <- gev.MLE.est$estimate["xi"]
		sig <- gev.MLE.est$estimate["sig"]
		mu <- gev.MLE.est$estimate["mu"]
		hoptests_gev <- hopt_estfn_GEV(xi=xi,sig=sig,mu=mu,d=d,n=n)
	}
	gam.MLE.est <- tryCatch(suppressWarnings(fitdistr(x=Uhat-min(Uhat)+.Machine$double.eps,densfun=dgamma,start=list(shape=1,scale=sd(Uhat)),lower=.Machine$double.eps, control=list(maxit=200))),
						  # warning=function(w) NA,
						  error=function(w) NA
						  )
	if(is.na(gam.MLE.est[1]) || any(!is.finite(gam.MLE.est$estimate))) {
	  # if(is.na(gam.MLE.est[1])) {
		gamnlogL <- Inf
		hoptests_gam <- Inf
	} else {
		gamnlogL <- -gam.MLE.est$loglik
		k <- gam.MLE.est$estimate["shape"]
		theta <- gam.MLE.est$estimate["scale"]
		hoptests_gam <- hopt_estfn_gam(k,theta,p,d,n)
	}
# 	norm.MLE.est <- tryCatch(fitdistr(x=Uhat,densfun="normal"),
# 						   warning=function(w) NA,
# 						   error=function(w) NA
# 						   )
# 	if(is.na(norm.MLE.est[1])) {
# 		normnlogL <- Inf
# 		hoptests_N <- Inf
# 	} else {
# 		normnlogL <- -norm.MLE.est$loglik
# 		mu <- norm.MLE.est$estimate["mean"]
# 		sig <- norm.MLE.est$estimate["sd"]
# 		hoptests_N <- hopt_estfn_N((mu+(0==mu)*.01),sig,d,n)
# 	}
	mu <- mean(Uhat); sig <- sd(Uhat)
	hoptests_N <- hopt_estfn_N((mu+(0==mu)*.01),sig,d,n)
	h <- min(hoptests_t,hoptests_N,hoptests_gev,hoptests_gam)
} else {
	h <- h_in
}
hhat <- h # h may be reassigned below; hhat is also returned

  # Given bandwith from above, estimate parameter vector
  # Solve for b, using plug-in (feasible) h from above
objfn2 <- function (b) {objfn(b,h)}
jacfn2 <- function (b) {jacfn(b,h)}

b_est1 <- tryCatch(newtonsys(Ffun=objfn2,x0=b_pilot,Jfun=jacfn2,maxiter=MAXITER),
				warning=function(w) NA,
				error=function(w) NA
				)
exitOK2 <- !(is.na(b_est1)[1] || (b_est1$niter==MAXITER)[1])
if (!exitOK2) {
		b_est2 <- tryCatch(newtonsys(Ffun=objfn2,x0=matrix(0,d,1),Jfun=jacfn2,maxiter=MAXITER),
				warning=function(w) NA,
				error=function(w) NA
				)
		exitOK3 <- !(is.na(b_est2)[1] || (b_est2$niter==MAXITER)[1])
		if (!(is.na(b_est2))[1] && !(is.na(b_est1))[1] && colSums(as.matrix(abs(objfn2(b_est2$zero))))<colSums(as.matrix(abs(objfn2(b_est1$zero))))[1]) {
		b_est1$zero <- b_est2$zero
		if (b_est1$niter==MAXITER) {
			exitOK2 <- exitOK3
		} 
	}
	}

exitOK <- exitOK2
if (!exitOK2) {
	b_win <- b_pilot
	h_init <- h_pilot
	objfn3 <- function (b) {objfn(b,h_init)}
    jacfn3 <- function (b) {jacfn(b,h_init)}
	f_win <- colSums(t(abs(objfn3(btmp$zero))))
	h_win <- h_pilot
	htmp <- h_init
	hstep <- abs((h_init)-h)/20
	btmp <- b_pilot
	while (htmp>=h) {
		btmp1 <- tryCatch(newtonsys(Ffun=objfn1,x0=btmp,Jfun=jacfn1,maxiter=MAXITER),
				warning=function(w) NA,
				error=function(w) NA
				)
	exitOK4 <- !(is.na(btmp1)[1] || (btmp1$niter==MAXITER)[1])
	if (!exitOK4) {
		break
	} else {
	b_win <- btmp1$zero
	f_win <- objfn1(btmp1$zero)
	h_win <- htmp
	htmp <- htmp-hstep
	}
}
suppressWarnings(b_est1$zero <- b_win) #b_est1 may be NA => innocuous warning
fval1 <- objfn2(b_win)
  # ONE last try . . . 
b_est3 <- tryCatch(newtonsys(Ffun=objfn2,x0=b_win,Jfun=jacfn2,maxiter=MAXITER),
				warning=function(w) NA,
				error=function(w) NA
				)
	exitOK5 <- !(is.na(b_est3)[1] || (b_est3$niter==MAXITER)[1])
	if (exitOK5[1] && colSums(as.matrix(abs(objfn2(b_est3$zero))))<colSums(as.matrix(abs(objfn2(b_est1$zero))))[1]) {
		b_est1$zero <- b_est3$zero
		exitOK <- exitOK5
		}	else {
	h <- h_win # smallest bandwidth that fsolve converges with
	}
}
b <- b_est1$zero
return(list(b=b,h=h,hhat=hhat))
}