bsp.graph <-
function(object,which,xlab,ylab,main){
if(which!=1 && which!=2) stop("The argument must be 1 or 2!!",call.=FALSE)
if(which ==1 && object$k1==0) stop(" A nonparametric function was not specified in the location parameter!!",call.=FALSE)
if(which ==2 && object$k2==0) stop(" A nonparametric function was not specified in the dispersion parameter!!",call.=FALSE)

   if(missingArg(xlab) || !is.character(xlab))xlab <- " "
       if(missingArg(ylab) || !is.character(ylab))ylab <- " "
       if(missingArg(main) || !is.character(main))main <- " "

quant005 <- function(x){quantile(x,prob=0.025)}
quant095 <- function(x){quantile(x,prob=0.975)}

n <- length(object$y)
chains <- object$chains
R <- nrow(chains)

if(which==1){
cov <- object$temp.mu
B <- object$B
p <- object$p
k1 <- object$k1
fs <- matrix(0,n,R)
for(i in 1:R){
fs[,i] <- B%*%chains[i,(p+1):(p+k1)]
}
li <- apply(fs, 1, quant005)
m <- apply(fs, 1, mean)
ls <- apply(fs, 1, quant095)
id <- sort(cov,index=TRUE)$ix
plot(cov[id], li[id], xlim=range(cov), ylim=range(li,ls), type="l", lty=3, xlab=xlab, ylab=ylab, main=main)
par(new=TRUE)
plot(cov[id], m[id], xlim=range(cov), ylim=range(li,ls), type="l", lty=1, xlab=xlab, ylab=ylab, main=main)
par(new=TRUE)
plot(cov[id], ls[id], xlim=range(cov), ylim=range(li,ls), type="l", lty=3, xlab=xlab, ylab=ylab, main=main)
}

if(which==2){
cov <- object$temp.dis
D <- object$D
p <- object$p
k1 <- object$k1
q <- object$q
k2 <- object$k2
fs <- matrix(0,n,R)
for(i in 1:R){
fs[,i] <- D%*%chains[i,(p+k1+1+q):(p+k1+q+k2)]
}
li <- apply(fs, 1, quant005)
m <- apply(fs, 1, mean)
ls <- apply(fs, 1, quant095)
id <- sort(cov,index=TRUE)$ix
plot(cov[id], li[id], xlim=range(cov), ylim=range(li,ls), type="l", lty=3, xlab=xlab, ylab=ylab, main=main)
par(new=TRUE)
plot(cov[id], m[id], xlim=range(cov), ylim=range(li,ls), type="l", lty=1, xlab=xlab, ylab=ylab, main=main)
par(new=TRUE)
plot(cov[id], ls[id], xlim=range(cov), ylim=range(li,ls), type="l", lty=3, xlab=xlab, ylab=ylab, main=main)
}


}
