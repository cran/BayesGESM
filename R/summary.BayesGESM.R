summary.BayesGESM <-
function(object, ...){
chains <- object$chains     
p <- object$p
q <- object$q
k1 <- object$k1
k2 <- object$k2

quant005 <- function(x){quantile(x,prob=0.025)}
quant095 <- function(x){quantile(x,prob=0.975)}

cat("\n       Error distribution:", object$family,"(",object$eta,")")
cat("\n              Sample size:", length(object$res))
cat("\n Size of posterior sample:", object$post.sam.s, "\n")

cat("\n =================   Location Submodel   =================")

if(p >0){
cat("\n =========   Parametric part     \n")
a <- round(apply(as.matrix(chains[,1:p]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,1:p]), 2, median), digits=4)
b <- round(apply(as.matrix(chains[,1:p]), 2, sd), digits=4)
c <- round(apply(as.matrix(chains[,1:p]), 2, quant005), digits=4)
d <- round(apply(as.matrix(chains[,1:p]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- colnames(object$X)
printCoefmat(e)
}


if(k1 >0){
cat("\n =========   Nonparametric part   \n")
a <- round(apply(as.matrix(chains[,(p+1):(p+k1)]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,(p+1):(p+k1)]), 2, median), digits=4)
b <- round(apply(as.matrix(chains[,(p+1):(p+k1)]), 2, sd), digits=4)
c <- round(apply(as.matrix(chains[,(p+1):(p+k1)]), 2, quant005), digits=4)
d <- round(apply(as.matrix(chains[,(p+1):(p+k1)]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- colnames(object$B)
printCoefmat(e)
}

cat("\n =================  Dispersion Submodel  ================= ")

if(q >0){
cat("\n =========   Parametric part    \n")
a <- round(apply(as.matrix(chains[,(p+k1+1):(p+k1+q)]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,(p+k1+1):(p+k1+q)]), 2, median), digits=4)
b <- round(apply(as.matrix(chains[,(p+k1+1):(p+k1+q)]), 2, sd), digits=4)
c <- round(apply(as.matrix(chains[,(p+k1+1):(p+k1+q)]), 2, quant005), digits=4)
d <- round(apply(as.matrix(chains[,(p+k1+1):(p+k1+q)]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- colnames(object$Z)
printCoefmat(e)
}

if(k2 >0){
cat("\n =========   Nonparametric part   \n")
a <- round(apply(as.matrix(chains[,(p+k1+1+q):(p+k1+q+k2)]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,(p+k1+1+q):(p+k1+q+k2)]), 2, median), digits=4)
b <- round(apply(as.matrix(chains[,(p+k1+1+q):(p+k1+q+k2)]), 2, sd), digits=4)
c <- round(apply(as.matrix(chains[,(p+k1+1+q):(p+k1+q+k2)]), 2, quant005), digits=4)
d <- round(apply(as.matrix(chains[,(p+k1+1+q):(p+k1+q+k2)]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- colnames(object$D)
printCoefmat(e)
}

cat("\n ==================     Model Selection Criteria     ==================")
cat("\n DIC=" , object$DIC, "   EAIC=" , object$EAIC, "   EBIC=" , object$EBIC, "   LMPL=" , object$LMPL, "\n")

}
