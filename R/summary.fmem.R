summary.fmem <-
function(object, ...){
chains <- object$chains     
p <- object$p
q <- object$q
k1 <- object$k1

quant005 <- function(x){quantile(x,prob=0.025)}
quant095 <- function(x){quantile(x,prob=0.975)}

cat("\n       Error distribution:", object$family,"(",object$eta,")")
cat("\n              Sample size:", length(object$y))
cat("\n Size of posterior sample:", object$post.sam.s, "\n")

if(p >0){
cat("\n =========   Covariates measured without error     \n")
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

cat("\n =========   Covariates measured with error     \n")
a <- round(apply(as.matrix(chains[,(p+1):(p+q)]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,(p+1):(p+q)]), 2, median), digits=4)
b <- round(apply(as.matrix(chains[,(p+1):(p+q)]), 2, sd), digits=4)
c <- round(apply(as.matrix(chains[,(p+1):(p+q)]), 2, quant005), digits=4)
d <- round(apply(as.matrix(chains[,(p+1):(p+q)]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- colnames(object$M)
printCoefmat(e)

if(k1 >0){
cat("\n =========   Nonparametric part   \n")
a <- round(apply(as.matrix(chains[,(p+3*q+2):(p+3*q+1+k1)]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,(p+3*q+2):(p+3*q+1+k1)]), 2, median), digits=4)
b <- round(apply(as.matrix(chains[,(p+3*q+2):(p+3*q+1+k1)]), 2, sd), digits=4)
c <- round(apply(as.matrix(chains[,(p+3*q+2):(p+3*q+1+k1)]), 2, quant005), digits=4)
d <- round(apply(as.matrix(chains[,(p+3*q+2):(p+3*q+1+k1)]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- colnames(object$B)
printCoefmat(e)
}

cat("\n =========   Dispersion parameter    \n")
a <- round(apply(as.matrix(chains[,(p+3*q+1)]), 2, mean), digits=4)
a2 <- round(apply(as.matrix(chains[,(p+3*q+1)]), 2, median), digits=4)
b <- round(apply(as.matrix(chains[,(p+3*q+1)]), 2, sd), digits=4)
c <- round(apply(as.matrix(chains[,(p+3*q+1)]), 2, quant005), digits=4)
d <- round(apply(as.matrix(chains[,(p+3*q+1)]), 2, quant095), digits=4)
e <- cbind(a,a2,b,c,d)
colnames(e) <- c("Mean ","Median","SD  ","      C.I.","95%    ")
rownames(e) <- "Sigma2_y  "
printCoefmat(e)
cat("\n  Ratio of the error variances: " , (object$omeg))


cat("\n\n ==========   Model Selection Criteria   ==========")
cat("\n DIC=" , object$DIC,  "   LMPL=" , object$LMPL, "\n")
}
