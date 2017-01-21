getwd()
library("EBImage")
library("dtt")
library("MASS")
library("lassoshooting")

vni <- function(n,i,N){
    if(n==0){out <- ifelse(i<=N,sqrt(1/N),NULL)
        } else {
            out <- ifelse(i<=N,sqrt(2/N)*cos(pi*n*(i-1/2)/N),NULL)
        }
    out
}


C <- cbind(vni(0,1:5,5),vni(1,1:5,5),vni(2,1:5,5),vni(3,1:5,5),vni(4,1:5,5))
I <- diag(5)
D <- cbind(I,C)
x = runif(5)

sqrt(sum(C[,1]^2))

temp1 <- apply(I*I,1,sum)
temp2 <- apply(I*C,1,sum)
temp3 <- apply(C*C,1,sum)

mu <- max(temp1,temp2,temp3)

#### code for matrching pursuit

MP <- function(x,D,maxiter, eps){
#### to normalize D
    D <- apply(D,2,function(x) {x/sqrt(sum(x^2))})
    nk <- integer(0)
    ak <- integer(0)
#### to initialize xk and Rk
    xk <- rep(0,length(x))
    Rk <- x
    for (iter in 1:maxiter){
        nk <- rbind(nk,which.max(abs(t(D)%*%Rk)))
        ak <- rbind(ak,sum(Rk*D[,nk[iter]]))
        xk1 <- xk + sum(Rk*D[,nk[iter]])*D[,nk[iter]]
        Rk1 <- Rk - sum(Rk*D[,nk[iter]])*D[,nk[iter]]
        if(sqrt(sum((x-xk1)^2))<eps){
            break
        } else {
            xk <- xk1
            Rk <- Rk1
        }        
    }
    return(list(xk=xk1,Rk=Rk1,nk=nk,ak = ak,itr=iter,err=sqrt(sum((x-xk1)^2))))
}

MP(x,D,100,10^(-5))



mumat <- t(D)%*%D
diag(mumat) <- 0
which(mumat==max(mumat),arr.ind=TRUE)

## so d_3 and d_10 are furthest apart
## coherence in D
mu <- sum(D[,3]*D[,10])

set.seed(12345)
for(N in c(30,60,90,120)){
    d <- (N/2)
    reps <- 30
    normerrMP <- matrix(0,nrow = d, ncol = reps)
    for (sprs in 1:d){
        for(MM in 1:reps){
            nz <- sample(x=1:N,size=sprs, replace =FALSE)
            a <- rep(0,N)
            a[nz] <- rnorm(n=sprs)
            C <- vni(0,1:(N/2),(N/2))
            for(i in 1:((N/2)-1)){
                C <- cbind(C,vni(i,1:(N/2),(N/2)))
            }
            I <- diag((N/2))
            D <- cbind(I,C)
            x <- D%*%a
            ahatMP <- rep(0,N)
            mp <- MP(x = x, D = D, maxiter = 50, eps = 10^(-5))
            ahatMP[mp$nk] <- ginv(D[,mp$nk])%*%mp$xk     
            normerrMP[sprs,MM] <- sqrt(sum((ahatMP-a)^2))/sqrt(sum(a^2))
        }
    }
}




