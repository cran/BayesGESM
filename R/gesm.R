gesm <- function(formula, data, family, eta, burn.in, post.sam.s, thin){

	if(family!="Normal" & family!="Student-t" & family!="Slash" & family!="Hyperbolic" & family!="ContNormal" & family!="Laplace")
	stop("Family of Distributions specified by the user is not supported, Check documentation!!",call.=FALSE)

	if(family=="Student-t" | family=="Slash" | family=="Hyperbolic" | family=="ContNormal"){
	if(missingArg(eta)) stop("for the Family of Distributions specified by the user an extra parameter is required!!", call.=FALSE)
	}

	if(missingArg(thin)){thin <- 1}
	else{
		if(thin<1 | thin != floor(thin))
		stop("the thin value must be a positive integer!!", call.=FALSE)
		}

	if (post.sam.s != floor(post.sam.s) | post.sam.s < 1  ){
		stop("Invalid posterior sample size value", call.=FALSE)
		}
	if (burn.in != floor(burn.in) | burn.in < 1){
		stop("Invalid burn-in value", call.=FALSE)
		}

	reem <- function(aa,b){
		 ag <- aa
		 ag <- sub("(", "", ag,fixed=TRUE)
		 ag <- sub(")", "", ag,fixed=TRUE)
		 ag <- sub(b, "", ag,fixed=TRUE)
		 ag <- strsplit(ag, ",")
		 ag <- ag[[1]][1]
		 ag
	 }

	mf <- match.call(expand.dots = FALSE)
	m <- match(c("formula"), names(mf), 0)
	mf <- mf[c(1, m)]

    if (missingArg(data)) 
	    data <- environment(formula)

	formula <- Formula(formula)
	if (length(formula)[2L] < 2L)
    	formula <- as.Formula(formula(formula), ~1)

	nl <- formula(formula, lhs=0, rhs=1)
	ll <- formula(formula, lhs=1, rhs=2)

	response <- as.matrix(eval(ll[[2]], data))
	n <- length(response)
	if(ncol(as.matrix(response)) > 1) stop("The response variable must be univariate!!", call.=FALSE)

	x <- model.matrix(nl, data)
	xa <- "bsp("
	xb <- colnames(x)
	idx <- grepl(xa,xb,fixed=TRUE)
	 
	if(sum(idx) > 1) stop("More than one nonparametric component is not supported in the location parameter!!",call.=FALSE)

	if(sum(idx) == 1){
		temp <- eval(parse(text=xb[idx]), data)
		temp.mu <- temp
		k1 <- attr(temp.mu,"kn")
		idx[1] <- TRUE
		if(length(idx) > 2){
			   X <- as.matrix(x[,!idx])
			   colnames(X) <- colnames(x)[!idx]
			   p <- sum(!idx)
			}else{p <- 0}
	 }else{
		k1 <- 0
	 	X <- as.matrix(x[,])
		colnames(X) <- colnames(x)
    	p <- ncol(X)
	 	}

	if(length(x)==0){
	  	X <- matrix(1,n,1)
	  	p <- 1
	}

	w <- model.matrix(ll, data)
	wa <- "bsp("
	wb <- colnames(w)
	idw <- grepl(wa,wb,fixed=TRUE)
	 
	if(sum(idw) > 1) stop("More than one nonparametric component is not supported in the dispersion parameter!!",call.=FALSE)

	if(sum(idw) == 1){
		temp <- eval(parse(text=wb[idw]), data)
		temp.dis <- temp
		k2 <- attr(temp.dis,"kn")
		idw[1] <- TRUE
		if(length(idw) > 2){
		   Z <- as.matrix(w[,!idw])
		   colnames(Z) <- colnames(w)[!idw]
		   q <- sum(!idw)
		}else{q <- 0}
	 }else{
	       k2 <- 0
	 	   Z <- as.matrix(w[,])
		   colnames(Z) <- colnames(w)
    	   q <- ncol(Z)
	 }


	par_ini <- list(p=0,k1=0,q=0,k2=0,family=family,n=n,y=response) 

	if((p==1 && var(X)==0) || p==0){
		  if(k1==0){
			  colnames(X) <- "Intercept"
			  par_ini$beta.i <- mean(response)
			  par_ini$p <- p
			  par_ini$k1 <- k1
			  par_ini$X <- X
			  rres <- log((response - mean(response))^2)
			  }
		  else{
			  B <- attr(temp.mu,"B")
			  k1 <- ncol(B)
			  colnames(B) <- paste("alpha",1:k1)
			  par_ini$p <- p
			  par_ini$k1 <- k1
			  par_ini$temp.mu <- temp.mu
			  par_ini$B <- B
			  par_ini$alpha.i <- solve(t(B)%*%B)%*%t(B)%*%response
			  rres <- log((response - B%*%par_ini$alpha.i)^2)
  			}
		}
	 else{
		  if(k1==0){
			 par_ini$p <- p
			 par_ini$X <- X
			 par_ini$k1 <- k1
			 par_ini$beta.i <- solve(t(X)%*%X)%*%t(X)%*%response
			 rres <- log((response -  X%*%par_ini$beta.i)^2)
		     }
     	  else{
			 nombres <- colnames(X)
			 B <- attr(temp.mu,"B")
			 k1 <- ncol(B)
			 par_ini$temp.mu <- temp.mu
			 colnames(B) <- paste("alpha",1:k1)
			 par_ini$p <- p
			 par_ini$X <- X
			 par_ini$k1 <- k1
			 par_ini$B <- B
			 X_au <- cbind(X,B)
			 b_au <- solve(t(X_au)%*%X_au)%*%t(X_au)%*%response
			 par_ini$beta.i <- b_au[1:p]
			 par_ini$alpha.i <- b_au[(p+1):(p+k1)]
			 rres <- log((response -  X_au%*%b_au)^2)
	     }  
 	}
	homosc <- 0
    if((q==1 && var(Z)==0) || q==0){
		if(k2==0){
			   q <- 1
			   colnames(Z) <- "Intercept"
			   par_ini$q <- q
			   par_ini$k2 <- k2
			   par_ini$Z <- Z
			   par_ini$gamma.i <- mean(rres)
			   homosc <- 1
	  		}
	  	else{
			   D <- attr(temp.dis,"B")
			   par_ini$temp.dis <- temp.dis
			   k2 <- ncol(D)
			   colnames(D) <- paste("lambda",1:k2)
			   par_ini$q <- q
			   par_ini$k2 <- k2
			   par_ini$D <- D
			   par_ini$lambda.i <- solve(t(D)%*%D)%*%t(D)%*%rres
			  }
		 }
	else{
		par_ini$Z <- Z
		par_ini$q <- q
			  if(k2==0){
					par_ini$k2 <- k2
					par_ini$gamma.i <- solve(t(Z)%*%Z)%*%t(Z)%*%rres
			  }
			  else{
					D <- attr(temp.dis,"B")
					par_ini$temp.dis <- temp.dis
					k2 <- ncol(D)
					colnames(D) <- paste("lambda",1:k2)
					par_ini$k2 <- k2
					par_ini$D <- D
					par_ini$q <- q
				    par_ini$Z <- Z 
					Z_au <- cbind(Z,D)
	                b_au <- solve(t(Z_au)%*%Z_au)%*%t(Z_au)%*%rres
					par_ini$gamma.i <- b_au[1:q]
					par_ini$lambda.i <- b_au[(q+1):(k2+q)]
     		  }
 	}



 
	 if(family=="Normal"){
	 eta <- 0
	    u <- function(s){
	     matrix(1,n,1)
	    }
	 pdf <- function(z){dnorm(z)}
	 cdf <- function(z){pnorm(z)}
	 }
	 if(family=="Student-t"){
	 if(eta[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
	    u <- function(s){
	     bb <- eta/2 + s/2
	 aa <- (eta + 1)/2
	 uii <- matrix(0,n,1)
	for(zz in 1:n){
	uii[zz] <- rgamma(1,shape=aa,scale=1/bb[zz])
	}
	uii
	    }
	 pdf <- function(z){dt(z,eta)}
	 cdf <- function(z){pt(z,eta)}
	 }
	 
	 if(family=="Slash"){
	 if(eta[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
	     u <- function(s){
	     bb <- s/2
	 aa <- eta/2  + 1
	  uii <- matrix(0,n,1)
	     for(zz in 1:n){
	  uii[zz] <- rtrunc(1, spec="gamma", a=0, b=1, shape=aa, scale=1/bb[zz]) 
	 }
	 uii
	    }
	    pdf <- function(z){eta*(z^2/2)^(-(eta+1/2))*gamma(eta+1/2)*pgamma(1,shape=(eta+1/2), scale=(2/z^2))/(2*pi)^(1/2)}
	    cdf <- function(z){temp <- matrix(0,length(z),1)          
	 for(i in 1:length(z)){    
	temp[i] <- integrate(pdf,-Inf,z[i])$value    
	}    
	temp    
	}   
	 }
	
	 if(family=="Laplace"){
	 eta <- 0
	    u <- function(s){
	     bb <- s
	 uii <- matrix(0,n,1)
	for(zz in 1:n){
	uii[zz] <- rgig(1,lambda=1/2,chi=bb[zz],psi=0.25)
	}
	uii
	    }
	 pdf <- function(z){dnormp(z,mu=0,sigmap=2,p=1)}
	 cdf <- function(z){pnormp(z,mu=0,sigmap=2,p=1)}
	 }
	 if(family=="Hyperbolic"){
	 if(eta[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
	    u <- function(s){
	     bb <- s + 1
	 uii <- matrix(0,n,1)
	for(zz in 1:n){
	    uii[zz] <- rgig(1,lambda=1/2,chi=bb[zz],psi=eta^2)
	}
	uii
	    }
	 dh <- function(z){exp(-eta*sqrt(1+z^2))}
	 pdf <- function(z){dh(z)/(2*integrate(dh,0,Inf)$value)}
	 cdf <- function(z){
	temp <- matrix(0,length(z),1)
	for(i in 1:length(z)){
	            temp[i] <- integrate(dh,-Inf,z[i])$value
	}
	temp/(2*integrate(dh,0,Inf)$value)
	}
	 }
	 if(family=="ContNormal"){
	 if(eta[1]<0 | eta[1]>1) stop("the extra parameter eta[1] must be within the interval (0,1)!!",call.=FALSE)
	 if(eta[2]<0 | eta[2]>1) stop("the extra parameter eta[2] must be within the interval (0,1)!!",call.=FALSE)
	     u <- function(s){
	     bb <- exp(-s*eta[2]/2)*eta[1]*sqrt(eta[2])
	 aa <- (1 - eta[1])*exp(-s/2)
	 bb <- bb/(aa + bb)
	 uu <- runif(n)
	 uii <- ifelse(uu<=bb,eta[2],1)
	uii
	    }
	 pdf <- function(z){eta[1]*dnorm(z,sd=1/sqrt(eta[2]))+(1 - eta[1])*dnorm(z)}
	 cdf <- function(z){eta[1]*pnorm(z,sd=1/sqrt(eta[2]))+(1 - eta[1])*pnorm(z)}
	 }
	 if(family=="Student-t" | family=="Slash" | family=="ContNormal"){
	 kappa <- function(u){
	 1/u  }
	 }
	 else{kappa <- function(u){
	 u  } }
	   par_ini$u <- u
	   par_ini$burn.in <- burn.in
	   par_ini$post.sam.s <- post.sam.s
	   par_ini$homosc <- homosc
	   par_ini$thin <- thin
	   par_ini$kappa <- kappa
	   par_ini$eta <- eta
	   par_ini$pdf <- pdf
	   par_ini$cdf <- cdf
	   result <- mcmc.gesm(par_ini)
	   par_ini$chains=result$chains
	   par_ini$DIC=result$DIC
	   par_ini$EAIC=result$EAIC
	   par_ini$EBIC=result$EBIC
	   par_ini$LMPL=result$LMPL
	   par_ini$res=result$res
	   par_ini$KL=result$KL
	   par_ini$X_2=result$X_2
	   class(par_ini) <- "gesm"
	   par_ini$call <- match.call()
	   par_ini
} 

