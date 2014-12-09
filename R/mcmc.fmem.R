mcmc.fmem <-
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


Wishart <- function(df,Omega){
  p <- ncol(Omega)
  x <- matrix(rnorm(df*p),df,p)
  chole <- chol(Omega)
  t(chole)%*%(t(x)%*%x)%*%chole
}  

burn.in <- params$burn.in
post.sam.s <- params$post.sam.s
thin <- params$thin
y <- params$y
p <- params$p
q <- params$q
k1 <- params$k1
u <- params$u
pdf <- params$pdf
cdf <- params$cdf
n <- params$n
kappa <- params$kappa
omeg <- 1/params$omeg

total <- burn.in + post.sam.s*thin
ancho <- floor(seq(2, total, length=10))

bar <- txtProgressBar(min=0, max=ancho[10], initial=0, width=50, char="+", style=3)

mu_m0 <- rep(0,q)
sigma2_mu0 <- Ident(q)*1000
sigma2_mu0.I <- solve(sigma2_mu0)

nu_m <- q
omega_m <- Ident(q)*1000

a.sigma2_y <- 0.0001/2
b.sigma2_y <- 0.0001/2

rho.a <- matrix(0,total,q)
rho0 <- rep(0,q)
S.rho <- Ident(q)*1000
    S.rho.I <- solve(S.rho)
rho.a[1,] <- params$rho.i
M <- params$M
params$rho0 <- rho0
params$S.rho <- S.rho

sigma2_m.a <- array(0,c(q,q,total))
sigma2_m.a2 <- array(0,c(q,q,total))
mu_m.a <- matrix(0,total,q)
m_i.a <- matrix(0,n,q)
sigma2_y.a <- matrix(0,total,1)

sigma2_m.a[,,1] <- as.matrix(var(M))
mu_m.a[1,] <- apply(M,2,mean)
m_i.a <- M
sigma2_y.a[1] <- params$sigma2_y


if(p > 0){
      beta.a <- matrix(0,total,p)
  beta0 <- rep(0,p)
  S.beta <- Ident(p)*1000
  S.beta.II <- matrix(0,p+q,p+q)
  S.beta.II[1:p,1:p] <- solve(S.beta)
  S.beta.II[(p+1):(p+q),(p+1):(p+q)] <- S.rho.I
  beta_rho0 <- matrix(0,p+q,1)
  beta_rho0[1:p] <- beta0
  beta_rho0[(p+1):(q+1)] <- rho0
  beta.a[1,] <- params$beta.i
  X <- params$X
  params$beta0 <- beta0
  params$S.beta <- S.beta}
else{beta.a <- matrix(0,total,1)
      X <- matrix(1,n,1)
   params$X <- X}

tau.a <- matrix(0,total,1)
a.tau <- 0.001
b.tau <- 0.001
if(k1 > 0){
   alpha.a <- matrix(0,total,k1)
   alpha0 <- rep(0,k1)
   alpha.a[1,] <- params$alpha.i
   B <- params$B
   params$alpha0 <- alpha0}
else{alpha.a <- matrix(0,total,1)
       B <- matrix(1,n,1)
    params$B <- B}

    cont <- 1
uii <- matrix(0,n,1)
inv.omega_m <- solve(omega_m)


for(l in 2:total){

M_m <- 0
omega_m.a <- matrix(0,ncol(M),ncol(M))

for(i in 1:n){
ss <- (y[i] - X[i,]%*%beta.a[l-1,] - m_i.a[i,]%*%rho.a[l-1,]- B[i,]%*%alpha.a[l-1,])^2/sigma2_y.a[l-1] + t(M[i,]-m_i.a[i,])%*%(M[i,]-m_i.a[i,])/(sigma2_y.a[l-1]*omeg) +
        t(m_i.a[i,]-mu_m.a[l-1,])%*%sigma2_m.a[,,l-1]%*%(m_i.a[i,]-mu_m.a[l-1,])
uii[i] <- u(ss)

sig.m_i.a <- solve((kappa(uii[i])*omeg*sigma2_y.a[l-1])^(-1)*Ident(ncol(M)) + sigma2_m.a[,,l-1]/kappa(uii[i])
 + rho.a[l-1,]%*%t(rho.a[l-1,])/(sigma2_y.a[l-1]*kappa(uii[i])))
chole <- chol(sig.m_i.a)

mu.m_i.a <- sig.m_i.a%*%(sigma2_m.a[,,l-1]%*%mu_m.a[l-1,]/kappa(uii[i]) + (kappa(uii[i])*omeg*sigma2_y.a[l-1])^(-1)*Ident(ncol(M))%*%M[i,]
+ rho.a[l-1,]*(y[i]-X[i,]%*%beta.a[l-1,] - B[i,]%*%alpha.a[l-1,])/(sigma2_y.a[l-1]*kappa(uii[i])))
m_i.a[i,] <- mu.m_i.a + t(chole)%*%rnorm(q)

omega_m.a <- omega_m.a + (m_i.a[i,]-mu_m.a[l-1,])%*%t(m_i.a[i,]-mu_m.a[l-1,])/kappa(uii[i])
M_m <- M_m + t(M[i,]- m_i.a[i,])%*%(M[i,] - m_i.a[i,])/kappa(uii[i])
}

if(p > 0){
MX <- cbind(X,m_i.a)
sig.beta.a <- solve(S.beta.II + t(MX)%*%(matrix(1/(sigma2_y.a[l-1] * kappa(uii)),n,(p+q))*MX))
mu.beta.a <-  sig.beta.a%*%(S.beta.II%*%beta_rho0 + t(MX)%*%(matrix(1/(sigma2_y.a[l-1]*kappa(uii)),n,1)*(y-B%*%alpha.a[l-1,])))
temp <- rmvnorm.l(mean=mu.beta.a, sigma=sig.beta.a)
beta.a[l,] <- temp[1:p]
rho.a[l,] <- temp[(p+1):(q+p)]
}
else{beta.a[l,] <- 0
    MX <- m_i.a
sig.beta.a <- solve(S.rho.I + t(MX)%*%(matrix(1/(sigma2_y.a[l-1] * kappa(uii)),n,(p+q))*MX))
mu.beta.a <-  sig.beta.a%*%(S.rho.I%*%rho0 + t(MX)%*%(matrix(1/(sigma2_y.a[l-1]*kappa(uii)),n,1)*(y-B%*%alpha.a[l-1,])))
rho.a[l,] <- rmvnorm.l(mean=mu.beta.a, sigma=sig.beta.a)
}


if(k1 > 0){
scale.tau <- 2*b.tau + t(alpha.a[l-1,]-alpha0)%*%(alpha.a[l-1,]-alpha0)
tau.a[l] <- 1/rgamma(1, shape= (k1/2 + a.tau), scale= 2/scale.tau)

sig.alpha.a<- solve(Ident(k1)/tau.a[l] + t(B)%*%(matrix(1/(sigma2_y.a[l-1]*kappa(uii)),n,k1)*B))
mu.alpha.a <- sig.alpha.a%*%(Ident(k1)%*%alpha0/tau.a[l] + t(B)%*%(matrix(1/(sigma2_y.a[l-1]*kappa(uii)),n,1)*(y-X%*%beta.a[l,]-m_i.a%*%rho.a[l,])))
alpha.a[l,] <- rmvnorm.l(mean=mu.alpha.a,sigma=sig.alpha.a)}

else{tau.a[l] <- 0
 alpha.a[l,] <- 0}

sig.mu_m.a <- solve(sigma2_mu0.I + sigma2_m.a[,,l-1]*n*mean(1/kappa(uii)))
mu.mu_m.a <-  sig.mu_m.a%*%(sigma2_mu0.I%*%mu_m0 + n*(sigma2_m.a[,,l-1])%*%apply(as.matrix(m_i.a)/matrix(kappa(uii),n,q),2,mean))
mu_m.a[l,] <- rmvnorm.l(mean=mu.mu_m.a, sigma=sig.mu_m.a)

sigma2_m.a[,,l] <- Wishart(n+q, solve(inv.omega_m + omega_m.a))
sigma2_m.a2[,,l] <- solve(sigma2_m.a[,,l])

scale.sigma2_y <- sum((y-X%*%beta.a[l,]-m_i.a%*%rho.a[l,]-B%*%alpha.a[l,])^2/kappa(uii)) + (1/omeg)*M_m + b.sigma2_y
sigma2_y.a[l] <- 1/rgamma(1,shape= ((n*(1+q) + a.sigma2_y)/2), scale=2/scale.sigma2_y)
 
if(l==ancho[cont]){
 Sys.sleep(0.5);
             setTxtProgressBar(bar,ancho[cont])
 cont <- cont + 1
}

l <- l + 1
}

close(bar)

size <- seq(burn.in+thin,total,length=post.sam.s)
aa <- matrix(0,post.sam.s,p+3*q+k1+2)
cad <- as.vector(" ")

if(p > 0){
aa[,1:p] <- beta.a[size,]
cad <- cbind(cad,t(paste("beta",1:p)))
}

aa[,(p+1):(p+q)] <- rho.a[size,]
cad <- cbind(cad,t(paste("rho",1:q)))

aa[,(p+q+1):(p+2*q)] <- mu_m.a[size,]
cad <- cbind(cad,t(paste("mu_m",1:q)))

for(i in 1:q){
aa[,(p+2*q+i)] <- as.vector(sigma2_m.a2[i,i,size])
}
cad <- cbind(cad,t(paste("sigma2_m",1:q)))

aa[,(p+3*q+1)] <- sigma2_y.a[size]
cad <- cbind(cad,t(paste("sigma2_y")))


if(k1 > 0){
aa[,(p+3*q+2):(p+3*q+1+k1)] <- alpha.a[size,]
aa[,(p+3*q+k1+2)] <- tau.a[size]
cad <- cbind(cad,t(paste("alpha",1:k1)))}

    mi <- matrix(0,n,post.sam.s)

D_bar <- 0
sigma_mp <- matrix(0,q+1,q+1)

for(i in 1:post.sam.s){
sigma_mp[1,1] <- sigma2_y.a[(burn.in+i)] + t(rho.a[(burn.in+i),])%*%sigma2_m.a2[,,(burn.in+i)]%*%rho.a[(burn.in+i),]
sigma_mp[1,2:(q+1)] <- rho.a[(burn.in+i),]%*%sigma2_m.a2[,,(burn.in+i)]
sigma_mp[2:(q+1),1] <- t(rho.a[(burn.in+i),]%*%sigma2_m.a2[,,(burn.in+i)])
sigma_mp[2:(q+1),2:(q+1)] <- omeg*sigma2_y.a[(burn.in+i)] + sigma2_m.a2[,,(burn.in+i)]
inv.sigma_mp <- solve(sigma_mp)
det_sigma_mp <- det(sigma_mp)
comp_s <- X%*%beta.a[(burn.in+i),] + B%*%alpha.a[(burn.in+i),] + sum(rho.a[(burn.in+i),]*mu_m.a[(burn.in+i),])
zz <- matrix(0,q+1,1)
for(j in 1:n){
zz[1] <- y[j]- comp_s[j]
zz[2:(q+1)] <- M[j,] - mu_m.a[(burn.in+i),]
ss <- sqrt(t(zz)%*%inv.sigma_mp%*%zz)
D_bar <- D_bar - 2*log(pdf(ss))
          mi[j,i] <- 1/(pdf(ss)*det_sigma_mp^(-(q+1)/2))
}
D_bar <- D_bar + n*log(det_sigma_mp)
}
D_bar <- D_bar/post.sam.s

D_theta <- 0
if(q >1){
comp_s <- X%*%apply(as.matrix(beta.a[size,]),2,mean) + B%*%apply(as.matrix(alpha.a[size,]),2,mean) + sum(apply(as.matrix(rho.a[size,]),2,mean)*apply(as.matrix(mu_m.a[(size),]),2,mean))
sigma_mp[1,1] <- mean(sigma2_y.a[size]) + t(apply(as.matrix(rho.a[size,]),2,mean))%*%apply(sigma2_m.a2[,,size],1:2,mean)%*%apply(as.matrix(rho.a[size,]),2,mean)
sigma_mp[1,2:(q+1)] <- t(apply(as.matrix(rho.a[size,]),2,mean))%*%apply(sigma2_m.a2[,,size],1:2,mean)
sigma_mp[2:(q+1),1] <- t(sigma_mp[1,2:(q+1)])
sigma_mp[2:(q+1),2:(q+1)] <- omeg*mean(sigma2_y.a[size]) + apply(sigma2_m.a2[,,size],1:2,mean)
}else{
comp_s <- X%*%apply(as.matrix(beta.a[size,]),2,mean) + B%*%apply(as.matrix(alpha.a[size,]),2,mean) + sum(apply(as.matrix(rho.a[size,]),2,mean)*apply(as.matrix(mu_m.a[(size),]),2,mean))
sigma_mp[1,1] <- mean(sigma2_y.a[size]) + t(apply(as.matrix(rho.a[size,]),2,mean))%*%mean(as.vector(sigma2_m.a2[,,size]))%*%apply(as.matrix(rho.a[size,]),2,mean)
sigma_mp[1,2:(q+1)] <- t(apply(as.matrix(rho.a[size,]),2,mean))%*%mean(as.vector(sigma2_m.a2[,,size]))
sigma_mp[2:(q+1),1] <- t(sigma_mp[1,2:(q+1)])
sigma_mp[2:(q+1),2:(q+1)] <- omeg*mean(sigma2_y.a[size]) + mean(as.vector(sigma2_m.a2[,,size]))
}
inv.sigma_mp <- solve(sigma_mp)
det_sigma_mp <- det(sigma_mp)
zz <- matrix(0,q+1,1)
for(j in 1:n){
zz[1] <- y[j]- comp_s[j]
zz[2:(q+1)] <- M[j,] - apply(as.matrix(mu_m.a[size,]),2,mean)
ss <- sqrt(t(zz)%*%inv.sigma_mp%*%zz)
D_theta <- D_theta - 2*log(pdf(ss))
}

D_theta <- D_theta + n*log(det_sigma_mp)
res_q <- qnorm(cdf((y-comp_s)/sqrt(sigma_mp[1,1])))

DIC <- 2*D_bar - D_theta
LMPL <- -sum(log(apply(mi,1,mean)))

KL <- log(apply(mi,1,mean))+ apply(-log(mi),1,mean) #distancia K-L
X_2 <- apply(mi^2,1,mean)/(apply(mi,1,mean))^2- 1#distancia X^2

haver <- apply(aa,2,var)
aa2 <- aa[,haver!=0]
if(k1>0) cad <- cbind(cad,"tau_alpha")
cad <- cad[2:length(cad)]
colnames(aa2) <- cad

list(chains=aa2, DIC=DIC, LMPL=LMPL, residuos=res_q, KL=KL, X_2=X_2)
}
