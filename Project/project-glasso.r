setwd("/Users/syedrahman/Documents/Fall2016/SparseApprox/Project")
library(MASS)
library(lassoshooting)
library(microbenchmark)
library("EBImage")
library("wavelets")
library("spatstat")
library("wavethresh")
library("squash")
library("raster")
library("png")
library(clusterGeneration)
library(Matrix)
library(MASS)
library(ROCR)
library("BDgraph")
source("projecthelper.r")

p <- 100
spars <- 0.5
set.seed(12345)
Omega <- genPositiveDefMat(p)$Sigma
threshold <- quantile(as.vector(Omega), c(spars)) 
Omega[abs(Omega)<threshold] <- 0
###Omega2 <- rgwish(n=1,adj.g=(abs(Omega*upper.tri(abs(Omega)>0,diag=TRUE))>0))[,,1]
Omega2 <- nearPD(Omega)$mat
Y <- mvrnorm(n = 2*p, mu = rep(0,p), Sigma = solve(Omega2), tol = 1e-6, empirical = FALSE, EISPACK = FALSE)

typeof(Omega2)
 
xold <- diag(p)
rho <- 1
eps <- 10^(-3)
t <- 0.1
   
L = max(eigen(Omega2)$values)

nlambda <- seq(from = 0.0001, to = 0.3, length = 30)

fistaerr <- rep(0,length(nlambda))
k <- 1
for(lambda in nlambda){
    fistaerr[k] <- glassofista(Y,xold,1/L,lambda,eps=10^(-20),Omega)$diff
    k <- k+1
}

fistaerr3 <- rep(0,length(nlambda))
k <- 1
for(lambda in nlambda){
    fistaerr3[k] <- glassofista3(Y,xold,1/L,lambda,eps=10^(-20),Omega)$diff
    k <- k+1
}

istaerr3 <- rep(0,length(nlambda))
k <- 1
for(lambda in nlambda){
    istaerr3[k] <- glassoista3(Y,xold,1/L,lambda,eps=10^(-20),Omega)$diff
    k <- k+1
}

istaerr <- rep(0,length(nlambda))
k <- 1
for(lambda in nlambda){
    istaerr[k] <- glassoista(Y,xold,1/L,lambda,eps=10^(-20),Omega)$diff
    k <- k+1
}

subgerr3 <- rep(0,length(nlambda))
k <- 1
for(lambda in nlambda){
    subgerr3[k] <- glassosubg3(Y,xold,1/L,lambda,eps=10^(-20),Omega)$diff
    k <- k+1
}

subgerr <- rep(0,length(nlambda))
k <- 1
for(lambda in nlambda){
    subgerr[k] <- glassosubg(Y,xold,1/L,lambda,eps=10^(-20),Omega)$diff
    k <- k+1
}

admmerr <- rep(0,length(nlambda))
k <- 1
for(lambda in nlambda){
    admmerr[k] <- glassoadmm(Y,xold,rho,lambda,eps=10^(-20),Omega)$diff
    k <- k+1
}

Omegahatadmm <- glassoadmm(Y,xold,1,0.01,eps=10^(-20),Omega2)$omegahat
Omegahataubg <- glassosubg3(Y,xold,1/L,1,eps=10^(-20),Omega2)$omegahat
Omegahatista <- glassoista3(Y,xold,1/L,0.02,eps=10^(-20),Omega2)$omegahat
Omegahatfista <- glassofista3(Y,xold,1/L,0.02,eps=10^(-20),Omega2)$omegahat
est <- glasso(s=(1/(2*p))*t(Y)%*%Y,rho = 0.02)

normomega <- sqrt(sum(Omega)^2)
subgerr3[1]/normomega
istaerr3[1]/normomega
fistaerr3[1]/normomega
admmerr[30]/normomega
1798.976/normomega

library(glasso)

names(est)
sqrt(sum(est$wi - Omega)^2)/normomega

glassoista3(Y,xold,1/L,0.0001,eps=10^(-20),Omega)$diff/normomega

sum(abs(Omega2)>10^(-5))
sum(abs(Omegahatadmm)>0.0001)
sum(abs(Omegahataubg)>0.0001)
sum(abs(Omegahatista)>0.0001)
sum(abs(Omegahatfista)>0.0001)


custConcord1 = function(S,n,rmax=100,eps=10^(-5),lambda){
(p=dim(S)[1])  
(omegahatcurrent = diag(p))
(r = 1)
(converged  = FALSE)
(maxdiff = eps/10)
while((converged == FALSE)&&(r<rmax)){
    cat('Iter:',r,maxdiff,'\n')
    (maxdiff = eps/10)
    (omegahatold = omegahatcurrent)
    for (i in 1:(p-1)){
    for (j in (i+1):p){
    if(i!=j){
        #cat('r=',r,'i=',i,'j=',j,'i!=j\n')
        (x = -t(omegahatcurrent[i,!(1:p == j)])%*%(S[j,!(1:p == j)])
         - t(omegahatcurrent[j,!(1:p == i)])%*%(S[i,!(1:p == i)]))
        (omegahatcurrent[i,j] = sign(x)*max(abs(x)-(lambda/n),0))
        (omegahatcurrent[i,j] = omegahatcurrent[i,j]/(S[i,i]+S[j,j]))
        (omegahatcurrent[j,i] = omegahatcurrent[i,j])
        (maxdiff = max(maxdiff, abs(omegahatcurrent[i,j]-omegahatold[i,j])))
    }}}
for(i in 1:p)
{
#cat('r=',r,'i=',i,'i==j\n')
    (omegahatcurrent[i,i] = -t(omegahatcurrent[i,!(1:p == i)])%*%(S[i,!(1:p == i)]) + sqrt((t(omegahatcurrent[i,!(1:p == i)])%*%(S[i,!(1:p == i)]))^2+4*S[i,i]))
    (omegahatcurrent[i,i] = omegahatcurrent[i,i]/(2*S[i,i]))
    (maxdiff = max(maxdiff, abs(omegahatcurrent[i,i]-omegahatold[i,i])))
}
if(maxdiff<eps){
converged = TRUE
}
else{
    r = r+1
}
}
return(list(Omega = omegahatcurrent,Iter = r))
}

S = (1/(2*p))*t(Y)%*%Y

custConcord1(S=S,n=2*p,lambda = 0.1)


tprfpr(list(abs(Omegahatadmm)>0.00001),(abs(Omega2)>0.0001))
tprfpr(list(abs(Omegahataubg)>0.00001),(abs(Omega2)>0.0001))
tprfpr(list(abs(Omegahatista)>0.00001),(abs(Omega2)>0.0001))
tprfpr(list(abs(Omegahatfista)>0.00001),(abs(Omega2)>0.0001))
tprfpr(list(abs(est$wi)>0.00001),(abs(Omega2)>0.0001))
