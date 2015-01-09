fmem <-
function(formula, data, omeg, family, eta, burn.in, post.sam.s, thin){

if(family!="Normal" & family!="Student-t" & family!="Slash" & family!="Hyperbolic" & family!="ContNormal" & family!="Laplace")
stop("Family of Distributions specified by the user is not supported, Check documentation!!",call.=FALSE)

if(family=="Student-t" | family=="Slash" | family=="Hyperbolic" | family=="ContNormal"){
 if(missingArg(eta)) stop("for the Family of Distributions specified by the user an extra parameter is required!!", call.=FALSE)
 }

if(family=="Laplace" | family=="Normal"){
 if(!missingArg(eta)) stop("for the Laplace and Normal distribution must be not specify the extra parameter!!", call.=FALSE)
 }

if(missingArg(thin)){thin <- 1}
if(missingArg(omeg)){omeg <- 1}
if(omeg<=0){stop("The value of the Ratio of the error variances must be positive", call.=FALSE)}

if(thin<1 | thin != floor(thin))
stop("the thin value must be a positive integer!!", call.=FALSE)
  if (post.sam.s != floor(post.sam.s) | post.sam.s < 1  ){ stop("Invalid posterior sample size value", call.=FALSE)}
if (burn.in != floor(burn.in) | burn.in < 1){ stop("Invalid burn-in value", call.=FALSE)}

     
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
 if(length(x)==0) stop("At least one covariate with measurement error must be specified!!", call.=FALSE) 

 M <- as.matrix(x[,-1])
 colnames(M) <- colnames(x)[-1]
 xa <- "bsp("
 xb <- colnames(M)
 q <- ncol(M)
 idx <- grepl(xa,xb,fixed=TRUE)
 if(sum(idx) >= 1) stop("Nonlinear effects of  covariates with measurement error  are not supported!!",call.=FALSE)


 w <- model.matrix(ll, data)
 wa <- "bsp("
 wb <- colnames(w)
 idw <- grepl(wa,wb,fixed=TRUE)
 
 if(sum(idw) > 1) stop("More than one nonparametric component is not supported!!",call.=FALSE)

 if(sum(idw) == 1){
	temp <- eval(parse(text=wb[idw]), data)
	k1 <- attr(temp,"kn")
	idw[1] <- TRUE
	if(length(idw) > 2){
   	X <- as.matrix(w[,!idw])
   colnames(X) <- colnames(w)[!idw]
   p <- sum(!idw)
}else{p <- 0}
 }else{
       k1 <- 0
    X <- as.matrix(w[,])
   colnames(X) <- colnames(w)
       p <- ncol(X)
 }



par_ini <- list(p=0,k1=0,k2=0,q=0,family=family,n=n,y=response)
par_ini$q <- q
par_ini$M <- M

if(p==0){
  if(k1==0){
 X_au <- cbind(X,M)
 b_au <- solve(t(X_au)%*%X_au)%*%t(X_au)%*%response
 par_ini$beta.i <- b_au[1]
 par_ini$rho.i <-  b_au[2:(q+p)]
 par_ini$p <- p
 par_ini$k1 <- k1
 par_ini$X <- X
 rres <- mean((response - X_au%*%b_au)^2)
  }
  else{
   B <- attr(temp,"B")
   k1 <- ncol(B)
 colnames(B) <- paste("alpha",1:k1)

  B_au <- cbind(M,B)
 b_au <- solve(t(B_au)%*%B_au)%*%t(B_au)%*%response

 par_ini$p <- p
 par_ini$k1 <- k1
 par_ini$B <- B
 par_ini$temp.mu <- temp
 par_ini$alpha.i <- b_au[(q+1):(q+k1)]
 par_ini$rho.i <- b_au[1:q]
 rres <- mean((response - B_au%*%b_au)^2)
  }
 }
 else{
     nombres <- colnames(X)
  if(k1==0){
 X_au <- cbind(X,M)
 b_au <- solve(t(X_au)%*%X_au)%*%t(X_au)%*%response
 
 par_ini$p <- p
 par_ini$X <- X
  k1 <- 0
 par_ini$k1 <- k1
 par_ini$beta.i <- b_au[1:p]
 par_ini$rho.i <- b_au[(p+1):(q+p)]
 
 rres <- mean((response -  X_au%*%b_au)^2)
     }
     else{
	 B <- attr(temp,"B")
	 k1 <- ncol(B)	 
  colnames(B) <- paste("alpha",1:k1)
  par_ini$p <- p
  par_ini$X <- X
  par_ini$k1 <- k1
  par_ini$B <- B
  X_au <- cbind(X,B,M)
  b_au <- solve(t(X_au)%*%X_au)%*%t(X_au)%*%response
  par_ini$temp.mu <- temp
  par_ini$beta.i <- b_au[1:p]
  par_ini$alpha.i <- b_au[(p+1):(p+k1)]
  par_ini$rho.i <- b_au[(p+k1+1):(p+k1+q)]
  rres <- mean((response -  X_au%*%b_au)^2)
     }  
 }

 
 
 if(family=="Normal"){
 eta <- 0
    u <- function(s){
     1
    }
 pdf <- function(z){dnorm(z)/(2*pi)^(q/2)}
 cdf <- function(z){pnorm(z)}
 }
 if(family=="Student-t"){
 if(eta[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
    u <- function(s){
     bb <- eta/2 + s/2
 aa <- (eta + q)/2 + 1
  rgamma(1,shape=aa,scale=1/bb)
    }
 pdf <- function(z){gamma((q+1+eta)/2)*(1+z^2/eta)^(-(q+1+eta)/2)/(gamma(eta/2)*(pi*eta)^((q+1)/2))}
 cdf <- function(z){pt(z,eta)}
 }
 if(family=="Slash"){
 if(eta[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
   u <- function(s){
     bb <- s/2
 aa <- q/2 + eta + 1
     rtrunc(1, spec="gamma", a=0, b=1, shape=aa, scale=1/bb) 
    }
 pdf <- function(z){eta*(z^2/2)^(-(eta+(q+1)/2))*gamma(eta+(q+1)/2)*pgamma(1,shape=(eta+(q+1)/2), scale=(2/z^2))/(2*pi)^((q+1)/2)}
 pdf2 <- function(z){eta*(z^2/2)^(-(eta+(1)/2))*gamma(eta+(1)/2)*pgamma(1,shape=(eta+(1)/2), scale=(2/z^2))/(2*pi)^((1)/2)}

 cdf <- function(z){temp <- matrix(0,length(z),1)          
 for(i in 1:length(z)){    
temp[i] <- integrate(pdf2,-Inf,z[i])$value    
}    
temp    
}   
 }   
 if(family=="Laplace"){
 eta <- 0
    u <- function(s){
     bb <- s
 aa <- -q/2
     rgig(1,lambda=aa,chi=bb,psi=0.25)

    }
 pdf <- function(z){besselK(sqrt(z^2/4),(-(q+1)/2+1))*(z^2)^(-((q+1)-2)/4)/((pi)^((q+1)/2))*4^(-(q+2)/2)}

 cdf <- function(z){pnormp(z,mu=0,sigmap=2,p=1)}
 }
 if(family=="Hyperbolic"){
 if(eta[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
    u <- function(s){
     bb <- s + 1
 aa <- -q/2
 rgig(1,lambda=aa,chi=bb,psi=eta^2)
    }

 pdf <- function(z){besselK((sqrt(z^2+1)*eta),(-(q+1)/2+1))*eta^((q+1)/2)*(z^2+1)^((-(q+1)+2)/4)/((2*pi)^((q+1)/2)*besselK(eta,1))}
 pdf2 <- function(z){besselK((sqrt(z^2+1)*eta),(-(1)/2+1))*eta^((1)/2)*(z^2+1)^((-(1)+2)/4)/((2*pi)^((1)/2)*besselK(eta,1))}

 cdf <- function(z){
temp <- matrix(0,length(z),1)
for(i in 1:length(z)){
            temp[i] <- integrate(pdf2,-Inf,z[i])$value
}
temp/(2*integrate(pdf2,0,Inf)$value)
}
 }
 
 if(family=="ContNormal"){
 if(eta[1]<0 | eta[1]>1) stop("the extra parameter eta[1] must be within the interval (0,1)!!",call.=FALSE)
 if(eta[2]<0 | eta[2]>1) stop("the extra parameter eta[2] must be within the interval (0,1)!!",call.=FALSE)
     u <- function(s){
     bb <- exp(-s*eta[2]/2)*eta[1]*eta[2]^(q/2 + 1)
 aa <- (1 - eta[1])*exp(-s/2)
 bb <- bb/(aa + bb)
 uu <- runif(1)
 uii <- ifelse(uu<=bb,eta[2],1)
uii
    }
 pdf <- function(z){eta[2]^((q+1)/2)*eta[1]*dnorm(z*sqrt(eta[2]))/(2*pi)^(q/2) + (1-eta[1])*dnorm(z)/(2*pi)^(q/2)}  

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
   par_ini$thin <- thin
      par_ini$kappa <- kappa
   par_ini$eta <- eta
   par_ini$pdf <- pdf
      par_ini$cdf <- cdf
   par_ini$sigma2_y <- rres
   par_ini$omeg <- omeg
   result <- mcmc.fmem(par_ini)
       par_ini$chains=result$chains
   par_ini$DIC=result$DIC
   par_ini$LMPL=result$LMPL
   par_ini$res=result$residuos
   par_ini$KL=result$KL
     par_ini$X_2=result$X_2
   class(par_ini) <- "fmem"
       par_ini$call <- match.call()
   par_ini
   }
