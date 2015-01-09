mcmc.gesm <-
function(params){

rmvnorm.l <- function(mean, sigma){
d <- length(mean)
sigma2 <- chol(sigma)
t(mean + t(sigma2)%*%rnorm(d))
}


Ident <- function(k){temp <- matrix(0,k,k)
  diag(temp) <- rep(1,k)
 temp
             }

MH.gamma    <- function(alpha,beta,tau.e,lambda,gamma,uii,parms){
 Omega <- matrix(0,length(gamma),length(gamma))
for(i in 1:parms$n){
Omega <- Omega +  as.numeric((1/parms$kappa(uii[i]))*(parms$y[i]-parms$X[i,]%*%beta-parms$B[i,]%*%alpha)^2/exp(parms$Z[i,]%*%gamma+ parms$D[i,]%*%lambda))*(parms$Z[i,]%*%t(parms$Z[i,]))/2
    }
    Omega <- Omega + solve(parms$S.gamma)
sigma2.gamma <- 0.5

cont <- 0
gen <- 0
while(gen==0){
gamma.es <- rmvnorm.l(mean=gamma,sigma=sigma2.gamma*solve(Omega))
dis.es <- exp(parms$Z%*%t(gamma.es) + parms$D%*%lambda)*parms$kappa(uii)
dis <- exp(parms$Z%*%gamma + parms$D%*%lambda)*parms$kappa(uii)
a <- exp(sum(log(sqrt(dis/dis.es))))
b.es <- sum((parms$y-parms$X%*%beta-parms$B%*%alpha)^2/dis.es) + t(t(gamma.es)-parms$gamma0)%*%solve(parms$S.gamma)%*%(t(gamma.es)-parms$gamma0)
b <- sum((parms$y-parms$X%*%beta-parms$B%*%alpha)^2/dis) + t(gamma-parms$gamma0)%*%solve(parms$S.gamma)%*%(gamma-parms$gamma0)
p <- min(1,a*exp(-(b.es-b)/2))
uni <- runif(1)
if(uni<p){gen <- 1}
cont <- cont + 1
}
list(gamma=gamma.es)
}


MH.lambda <- function(alpha,beta,tau.e,gamma,lambda,uii,parms){
 Omega <- matrix(0,length(lambda),length(lambda))
for(i in 1:parms$n){
Omega <- Omega +  as.numeric((1/parms$kappa(uii[i]))*(parms$y[i]-parms$X[i,]%*%beta-parms$B[i,]%*%alpha)^2/exp(parms$Z[i,]%*%gamma+ parms$D[i,]%*%lambda))*(parms$D[i,]%*%t(parms$D[i,]))/2
    }
    Omega <- Omega + solve(tau.e*Ident(parms$k2))
sigma2.lambda <- 0.5

cont <- 0
gen <- 0
while(gen==0){
lambda.es <- rmvnorm.l(mean=lambda,sigma=sigma2.lambda*solve(Omega))
dis.es <- exp(parms$Z%*%gamma + parms$D%*%t(lambda.es))*parms$kappa(uii)
dis <- exp(parms$Z%*%gamma + parms$D%*%lambda)*parms$kappa(uii)
a <- exp(sum(log(sqrt(dis/dis.es))))
b.es <- sum((parms$y-parms$X%*%beta-parms$B%*%alpha)^2/dis.es) + t(t(lambda.es)-parms$lambda0)%*%solve(tau.e*Ident(parms$k2))%*%(t(lambda.es)-parms$lambda0)
b <-    sum((parms$y-parms$X%*%beta-parms$B%*%alpha)^2/dis) + t(lambda-parms$lambda0)%*%solve(tau.e*Ident(parms$k2))%*%(lambda-parms$lambda0)
p <- min(1,a*exp(-(b.es-b)/2))
uni <- runif(1)
if(uni<p){gen <- 1}
cont <- cont + 1
}
list(lambda=lambda.es)
}
 
 
burn.in <- params$burn.in
post.sam.s <- params$post.sam.s
thin <- params$thin
homosc <- params$homosc
y <- params$y
p <- params$p
q <- params$q
k1 <- params$k1
k2 <- params$k2
u <- params$u
pdf <- params$pdf
cdf <- params$cdf
n <- params$n
kappa <- params$kappa

total <- burn.in + post.sam.s*thin
ancho <- floor(seq(2, total, length=10))

bar <- txtProgressBar(min=0, max=ancho[10], initial=0, width=50, char="+", style=3)

if(p > 0){beta.a <- matrix(0,total,p)
  beta0 <- rep(0,p)
  S.beta <- Ident(p)*100000
  beta.a[1,] <- params$beta.i
  X <- params$X
  params$beta0 <- beta0
  params$S.beta <- S.beta}
else{beta.a <- matrix(0,total,1)
      X <- matrix(1,n,1)
   params$X <- X}

if(q > 0){gamma.a <- matrix(0,total,q)
   gamma0 <- rep(0,q)
  S.gamma <- Ident(q)*100000
  gamma.a[1,] <- params$gamma.i
  Z <- params$Z
  params$gamma0 <- gamma0
  params$S.gamma <- S.gamma}
else{gamma.a <- matrix(0,total,1)
       Z <- matrix(1,n,1)
   params$Z <- Z}

tau.a <- matrix(0,total,1)
a.tau <- 0.00001
b.tau <- 0.00001
if(k1 > 0){alpha.a <- matrix(0,total,k1)
   alpha0 <- rep(0,k1)
   alpha.a[1,] <- params$alpha.i
   B <- params$B
   params$alpha0 <- alpha0}
else{alpha.a <- matrix(0,total,1)
       B <- matrix(1,n,1)
    params$B <- B}

tau.a.l <- matrix(0,total,1)
a.tau.l <- 0.00001
b.tau.l <- 0.00001
if(k2 > 0){lambda.a <- matrix(0,total,k2)
   lambda0 <- rep(0, k2)
   lambda.a[1,] <- params$lambda.i
   D <- params$D
   params$lambda0 <- lambda0}
else{lambda.a <- matrix(0,total,1)
   D <- matrix(1,n,1)
       params$D <- D}

    cont <- 1
for(l in 2:total){

s <- (y-X%*%beta.a[l-1,]-B%*%alpha.a[l-1,])^2/(exp(Z%*%gamma.a[l-1,] + D%*%lambda.a[l-1,]))
uii <- u(s)

Sigma <- exp(Z%*%gamma.a[l-1,] + D%*%lambda.a[l-1,])*kappa(uii)

if(k1 > 0){
scale.tau <- t(alpha.a[l-1,]-alpha0)%*%(alpha.a[l-1,]-alpha0)+2*b.tau
tau.a[l] <- 1/rgamma(1,shape=(k1/2+a.tau), scale=2/scale.tau)}
else{tau.a[l] <- 0}

if(k2 > 0){
scale.tau.l <- t(lambda.a[l-1,]-lambda0)%*%(lambda.a[l-1,]-lambda0)+2*b.tau.l
tau.a.l[l] <- 1/rgamma(1,shape=(k2/2+a.tau.l), scale=2/scale.tau.l)}
else{tau.a.l[l] <- 0}

if(k1 > 0){
sig.alpha.a<- solve(Ident(k1)/tau.a[l] + t(B)%*%(matrix(1/Sigma,n,ncol(B))*B))
mu.alpha.a <- sig.alpha.a%*%(Ident(k1)%*%alpha0/tau.a[l] + t(matrix(1/Sigma,n,ncol(B))*B)%*%(y-X%*%beta.a[l-1,]))
alpha.a[l,] <- rmvnorm.l(mean=mu.alpha.a,sigma=sig.alpha.a)}
else{alpha.a[l,] <- 0}

if(p > 0){
sig.beta.a <- solve(solve(S.beta)+t(X)%*%(matrix(1/Sigma,n,ncol(X))*X))
mu.beta.a <- sig.beta.a%*%(solve(S.beta)%*%beta0+t(X*matrix(1/Sigma,n,ncol(X)))%*%(y-B%*%alpha.a[l,]))
beta.a[l,] <- rmvnorm.l(mean=mu.beta.a,sigma=sig.beta.a)}
else{beta.a[l,] <- 0}

if(k2 > 0){
salida <- MH.lambda(alpha.a[l,],beta.a[l,],tau.a.l[l],gamma.a[l-1,],lambda.a[l-1,],uii,params)
lambda.a[l,] <- salida$lambda}
else{lambda.a[l,] <- 0}

if(q > 0){
    if(homosc == 1){ gamma.a[l,] <- log(1/rgamma(1,shape=(n/2 -1), scale=(2/sum((y-X%*%beta.a[l,]-B%*%alpha.a[l,])^2/kappa(uii)))))}
else{salida2 <- MH.gamma(alpha.a[l,],beta.a[l,],tau.a.l[l],lambda.a[l,],gamma.a[l-1,],uii,params)
gamma.a[l,] <- salida2$gamma}
}
     else{gamma.a[l,] <- 0}

if(l==ancho[cont]){
 Sys.sleep(0.5);
             setTxtProgressBar(bar,ancho[cont])
 cont <- cont + 1
}

l <- l + 1
}

close(bar)

size <- seq(burn.in+thin,total,length=post.sam.s)
aa <- matrix(0,post.sam.s,p+q+k1+k2+2)
cad <- as.vector(" ")

if(p > 0){
aa[,1:p] <- beta.a[size,]
cad <- cbind(cad,t(paste("beta",1:p)))}
if(k1 > 0){
aa[,(p+1):(p+k1)] <- alpha.a[size,]
aa[,(p+k1+q+k2+1)] <- tau.a[size,]
cad <- cbind(cad,t(paste("alpha",1:k1)))}

if(q > 0){
aa[,(p+k1+1):(p+k1+q)] <- gamma.a[size,]
cad <- cbind(cad,t(paste("gamma",1:q)))}

if(k2 >0){
aa[,(p+k1+1+q):(p+k1+q+k2)] <- lambda.a[size,]
aa[,(p+k1+k2+q+2)] <- tau.a.l[size,]
    cad <- cbind(cad,t(paste("lambda",1:k2)))}

mi <- matrix(0,n,post.sam.s)
res_q <- matrix(0,n,post.sam.s)

D_bar <- 0
for(i in 1:post.sam.s){
mui <- (X%*%beta.a[(burn.in+i),]+B%*%alpha.a[(burn.in+i),])
sigma2i <- exp(Z%*%gamma.a[(burn.in+i),] + D%*%lambda.a[(burn.in+i),])
zii <- (y-mui)/sqrt(sigma2i)
D_bar <- D_bar - 2*sum((log(pdf(zii)) -log(sigma2i)/2))
      mi[,i] <- 1/(pdf(zii)/sqrt(sigma2i))
res_q[,i] <- qnorm(cdf(zii))
}
D_bar <- D_bar/post.sam.s

D_theta <- 0
mu_bar <- X%*%apply(as.matrix(beta.a[size,]),2,mean) + B%*%apply(as.matrix(alpha.a[size,]),2,mean)
sigma2_bar <- exp(Z%*%apply(as.matrix(gamma.a[size,]),2,mean) + D%*%apply(as.matrix(lambda.a[size,]),2,mean))
zii <- (y-mu_bar)/sqrt(sigma2_bar)
D_theta <- -2*sum((log(pdf(zii)) -log(sigma2_bar)/2))

DIC <- 2*D_bar - D_theta
LMPL <- -sum(log(apply(mi,1,mean)))
res <- apply(res_q,1,mean)            #residuos
KL <- log(apply(mi,1,mean))+ apply(-log(mi),1,mean) #distancia K-L
X_2 <- apply(mi^2,1,mean)/(apply(mi,1,mean))^2- 1#distancia X^2

if(k1 > 0){
H_a <- sum(diag(B%*%solve(t(B)%*%(matrix(1/sigma2_bar,n,ncol(B))*B) + apply(as.matrix(tau.a[size,]),2,mean)^(-1)*Ident(k1))%*%t(B*matrix(1/sigma2_bar,n,ncol(B)))))
if(k2 > 0){
H_l <- sum(diag(D%*%solve(t(D)%*%D + apply(as.matrix(tau.a.l[size,]),2,mean)^(-1)*Ident(k2))%*%t(D)))
EAIC <- D_bar + 2*(p + q + H_a + H_l)
EBIC <- D_bar + log(n)+(p + q + H_a + H_l)
}
else{EAIC <- D_bar + 2*(p + q + H_a)
 EBIC <- D_bar + log(n)+(p + q + H_a)}
}
else{if(k2 > 0){
H_l <- sum(diag(D%*%solve(t(D)%*%D + apply(as.matrix(tau.a.l[size,]),2,mean)^(-1)*Ident(k2))%*%t(D)))
EAIC <- D_bar + 2*(p + q  + H_l)
EBIC <- D_bar + log(n)+(p + q + H_l)
}
else{EAIC <- D_bar + 2*(p + q)
 EBIC <- D_bar + log(n)+(p + q)}
}


haver <- apply(aa,2,var)
aa2 <- aa[,haver!=0]
if(k1>0) cad <- cbind(cad,"tau_alpha  1")
if(k2>0) cad <- cbind(cad,"tau_lambda 2")
cad <- cad[2:length(cad)]
colnames(aa2) <- cad

list(DIC=DIC, EAIC=EAIC, EBIC=EBIC, LMPL=LMPL, chains=aa2, residuos=res, KL=KL, X_2=X_2)

}
