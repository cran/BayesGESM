BayesGESM <-
function(response, formula.mu, np.mu, formula.dis, np.dis, family, eta, burn.in, post.sam.s, thin){

if(family!="Normal" & family!="Student-t" & family!="Slash" & family!="Hyperbolic" & family!="ContNormal" & family!="Laplace")
stop("Family of Distributions specified by the user is not supported, Check documentation!!",call.=FALSE)

if(family=="Student-t" | family=="Slash" | family=="Hyperbolic" | family=="ContNormal"){
 if(missingArg(eta)) stop("for the Family of Distributions specified by the user an extra parameter is required!!", call.=FALSE)
 }

if(missingArg(thin)){thin <- 1}
else{if(thin<1 | thin != floor(thin))
stop("the thin value must be a positive integer!!", call.=FALSE)}
    if(ncol(as.matrix(response)) > 1) stop("The response variable must be univariate!!", call.=FALSE)
if (post.sam.s != floor(post.sam.s) | post.sam.s < 1  ){ stop("Invalid posterior sample size value", call.=FALSE)}
if (burn.in != floor(burn.in) | burn.in < 1){ stop("Invalid burn-in value", call.=FALSE)}

   n <- length(response)
par_ini <- list(p=0,k1=0,q=0,k2=0,family=family,n=n,y=response) 
     if(missingArg(formula.mu)){
  if(missingArg(np.mu)){
   p <- 1
 X <- matrix(1,n,1)
 colnames(X) <- "Intercept"
 k1 <- 0
 par_ini$beta.i <- mean(response)
 par_ini$p <- p
 par_ini$k1 <- k1
 par_ini$X <- X
 rres <- log((response - mean(response))^2)
  }
  else{
   p <- 0
 k1 <- floor(n^(1/5))
 B <- bs(np.mu,knots=quantile(np.mu,prob=seq(1,k1,by=1)/(k1+1)),intercept=TRUE)
 k1 <- ncol(B)
 colnames(B) <- paste("alpha",1:k1)
 par_ini$p <- p
 par_ini$k1 <- k1
 par_ini$B <- B
 par_ini$alpha.i <- solve(t(B)%*%B)%*%t(B)%*%response
 rres <- log((response - B%*%par_ini$alpha.i)^2)
  }
 }
 else{
  if(missingArg(np.mu)){
 mf <- model.frame(formula=formula.mu) 
 X <- as.matrix(model.matrix(attr(mf, "terms"),data=mf))
 colnames(X) <- colnames(model.matrix(attr(mf, "terms"), data=mf))
 p <- ncol(X)
 par_ini$p <- p
 par_ini$X <- X
  k1 <- 0
 par_ini$k1 <- k1
 par_ini$beta.i <- solve(t(X)%*%X)%*%t(X)%*%response
 rres <- log((response -  X%*%par_ini$beta.i)^2)
     }
     else{
  mf <- model.frame(formula=formula.mu) 
  X <- as.matrix(model.matrix(attr(mf, "terms"),data=mf))
  colnames(X) <- colnames(as.matrix(model.matrix(attr(mf, "terms"),data=mf)))
  nombres <- colnames(X)
  X <- as.matrix(X[,2:ncol(X)])
  colnames(X) <- nombres[2:length(nombres)]
  p <- ncol(X)
    k1 <- floor(n^(1/5))
  B <- bs(np.mu,knots=quantile(np.mu,prob=seq(1,k1,by=1)/(k1+1)),intercept=TRUE)
  k1 <- ncol(B)
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
     if(missingArg(formula.dis)){
  if(missingArg(np.dis)){
   q <- 1
 Z <- matrix(1,n,1)
 colnames(Z) <- "Intercept"
 k2 <- 0
 par_ini$q <- q
 par_ini$k2 <- k2
 par_ini$Z <- Z
 par_ini$gamma.i <- mean(rres)
 homosc <- 1
  }
  else{
   q <- 0
 k2 <- floor(n^(1/5))
 D <- bs(np.dis,knots=quantile(np.dis,prob=seq(1,k2,by=1)/(k2+1)),intercept=TRUE)
 k2 <- ncol(D)
 colnames(D) <- paste("lambda",1:k2)
 par_ini$q <- q
 par_ini$k2 <- k2
 par_ini$D <- D
 par_ini$lambda.i <- solve(t(D)%*%D)%*%t(D)%*%rres
  }
 }
 else{
 mf <- model.frame(formula=formula.dis)
 Z <- as.matrix(model.matrix(attr(mf, "terms"),data=mf))
 colnames(Z) <- colnames(model.matrix(attr(mf, "terms"), data=mf))
 q <- ncol(Z)
 par_ini$Z <- Z
 par_ini$q <- q
  if(missingArg(np.dis)){
 k2 <- 0
par_ini$k2 <- k2
par_ini$gamma.i <- solve(t(Z)%*%Z)%*%t(Z)%*%rres
     }
  else{
 mf <- model.frame(formula=formula.mu) 
 Z <- as.matrix(model.matrix(attr(mf, "terms"),data=mf))
      colnames(Z) <- colnames(as.matrix(model.matrix(attr(mf, "terms"),data=mf)))
        nombres <- colnames(Z)
 Z <- as.matrix(Z[,2:ncol(Z)])
 colnames(Z) <- nombres[2:length(nombres)]
 q <- ncol(Z)
 k2 <- floor(n^(1/5))
 D <- bs(np.dis,knots=quantile(np.dis,prob=seq(1,k2,by=1)/(k2+1)),intercept=TRUE)
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
 H <- function(a,x){gamma(a)*gamma_inc_P(a,x)/(x^a)}          
 z <- function(x){H(eta+1/2,x^2/2)}
 pdf <- function(z){z(z)/(2*integrate(z,0,Inf)$value)}
 cdf <- function(z){temp <- matrix(0,length(z),1)
 for(i in 1:length(z)){
temp[i] <- integrate(z,-Inf,z[i])$value
}
temp/(2*integrate(z,0,Inf)$value)
}
 }
 if(family=="Laplace"){
 eta <- 0
    u <- function(s){
     bb <- s
 uii <- matrix(0,n,1)
for(zz in 1:n){
uii[zz] <- rgig(1,lambda=1/2,chi=bb[zz],psi=1)
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
    uii[zz] <- rgig(1,lambda=1/2,chi=bb[zz],psi=eta)
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
   result <- mcmc(par_ini)
       par_ini$chains=result$chains
   par_ini$DIC=result$DIC
   par_ini$EAIC=result$EAIC
   par_ini$EBIC=result$EBIC
   par_ini$LMPL=result$LMPL
   par_ini$res=result$res
      par_ini$KL=result$KL
      par_ini$X_2=result$X_2
   class(par_ini) <- "BayesGESM"
       par_ini$call <- match.call()
   par_ini
   }
