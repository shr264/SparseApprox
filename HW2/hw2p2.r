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

OMP <- function(x, D, maxiter,eps){
#### to normalize D
    D <- apply(D,2,function(x) {x/sqrt(sum(x^2))})
    nk <- integer(0)
    ak <- integer(0)
#### to initialize xk and Rk
    xk <- rep(0,length(x))
    Rk <- x
    uk <- matrix(0,nrow = dim(D)[1], ncol = dim(D)[2])
    nk <- rbind(nk,which.max(abs(t(D)%*%Rk)))
    uk[,1] <- D[,nk[1]]
    ak <- rbind(ak,(sum(Rk*uk[,1])/sum(uk[,1]*uk[,1])))
    xk1 <- xk + (sum(Rk*uk[,1])/sum(uk[,1]*uk[,1]))*uk[,1]
    Rk1 <- Rk - (sum(Rk*uk[,1])/sum(uk[,1]*uk[,1]))*uk[,1]
    if(sqrt(sum((xk1-x)^2))<eps | maxiter ==1){
        return(list(xk=xk1,Rk=Rk1,nk = nk, ak = ak, itr=1,err=sqrt(sum((x-xk1)^2))))
    }  else    {
            xk <- xk1
            Rk <- Rk1
            #iter = 2
            for (iter in 2:min(dim(D)[2],maxiter)){
                nk <- rbind(nk,which.max(abs(t(D)%*%Rk)))
                sumterm = 0
                for(j in 1:(iter-1)){
                    sumterm = sumterm + (sum(D[,nk[iter]]*uk[,j])/sum(uk[,j]*uk[,j]))*uk[,j]
                }
                uk[,iter] <- D[,nk[iter]] - sumterm
                ak <- rbind(ak,(sum(Rk*uk[,iter])/sum(uk[,iter]*uk[,iter])))

                xk1 <- xk + (sum(Rk*uk[,iter])/sum(uk[,iter]*uk[,iter]))*uk[,iter]
                Rk1 <- Rk - (sum(Rk*uk[,iter])/sum(uk[,iter]*uk[,iter]))*uk[,iter]
                if(sqrt(sum((xk1-x)^2))<eps){
                    break
                } else {
                    xk <- xk1
                    Rk <- Rk1
                }        
            }
            return(list(xk=xk1,Rk=Rk1,nk = nk, ak = ak, itr=iter,err=sqrt(sum((x-xk1)^2))))          
        }
    }

OMP(x,D,10,10^(-5))

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
normerrMP <- normerrOMP <- normerrBP <- matrix(0,nrow = d, ncol = reps)
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
        ahatOMP <- rep(0,N)
        mp <- MP(x = x, D = D, maxiter = 50, eps = 10^(-5))
        omp <- OMP(x = x, D = D, maxiter = 50, eps = 10^(-5))
        ahatMP[mp$nk] <- ginv(D[,mp$nk])%*%mp$xk
        ahatOMP[omp$nk] <- ginv(D[,omp$nk])%*%omp$xk
        ahatBP <- BP(D,x)
        
        normerrMP[sprs,MM] <- sqrt(sum((ahatMP-a)^2))/sqrt(sum(a^2))
        normerrOMP[sprs,MM] <- sqrt(sum((ahatOMP-a)^2))/sqrt(sum(a^2))
        normerrBP[sprs,MM] <- sqrt(sum((ahatBP-a)^2))/sqrt(sum(a^2))
    }
}


filename = paste('normerr',toString(d),'.pdf', sep = '')
pdf(filename)
plot(1:d,apply(normerrOMP,1,mean), xlab = "# of non-sparse elements", ylab = "Relative Normed Error", type = 'l', col = 2)
lines(1:d,apply(normerrMP,1,mean), xlab = "# of non-sparse elements", ylab = "Relative Normed Error", type = 'l', col = 3)
lines(1:d,apply(normerrBP,1,mean), xlab = "# of non-sparse elements", ylab = "Relative Normed Error", type = 'l', col = 4)
legend("bottomright",c("OMP","MP","BP"),lty = rep(1,3),col = 2:4)
dev.off()
}

### When d = 15, both MP and OMP tend to do well when the number of non-sparse elements are very, i.e. 1. MP starts to fail pretty quickly. Even at 2 non-sparse elements, it performs poorly. However, OMP does much better. Only starts to fail after the number of non-sparse elements reaches 3 or so. 3/15 = 0.2. Hence as the number of non-sparse elements reaches 20%, OMP starts to fail. For d = 30, OMP starts to fail at 11. 11/30 = 0.37. d = 45, OMP starts to fail at about 18. 18/45 = 0.4. For d = 60, OMP starts to fail at about 20. 20/60 = 0.3.

###For the basis pursuit algorithm, we use a lasso regression to estimate a from x. So basically min_a (1/2)(x-\Phi a)_2^2 + \lambda (a)_1

cv5BP <- function(D,x,nlambda){
foldsize <- floor(length(x)/5)    
err <- matrix(0,nrow = nlambda, ncol = 5)
rndmsmpl <- sample(1:length(x))
folds <- split(rndmsmpl, ceiling(seq_along(rndmsmpl)/foldsize))
for(lambda in 1:nlambda){
    for(i in 1:5){
        a1 <- folds[[i]]
        fit <- lassoshooting(X=D[-folds[[i]],],y=x[-folds[[i]]],lambda = 0.01*lambda)
        err[lambda,i] <- sum((a1 - D[folds[[i]],]%*%fit$coeff)^2)
    }
}
0.01*which.min(apply(err,1,mean))
}
cvBP(D,x,100)

BP <- function(D,x){
    lambdastar <- cv5BP(D,x,100)
    lassoshooting(X=D,y=x,lambda = lambdastar)$coeff
}
BP(D,x)

set.seed(12345)
for(N in c(30,60,90,120)){
d <- (N/2)
reps <- 30
normerrMP <- normerrOMP <- normerrBP <- matrix(0,nrow = d, ncol = reps)
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
        ahatOMP <- rep(0,N)
        mp <- MP(x = x, D = D, maxiter = 50, eps = 10^(-5))
        omp <- OMP(x = x, D = D, maxiter = 50, eps = 10^(-5))
        for(i in 1:length(mp$nk)){
            ahatMP[mp$nk[i]] <- ahatMP[mp$nk[i]]+mp$ak[i]
        }
        ahatOMP[omp$nk] <- ginv(D[,omp$nk])%*%omp$xk
        ahatBP <- BP(D,x)
        
        normerrMP[sprs,MM] <- sqrt(sum((ahatMP-a)^2))/sqrt(sum(a^2))
        normerrOMP[sprs,MM] <- sqrt(sum((ahatOMP-a)^2))/sqrt(sum(a^2))
        normerrBP[sprs,MM] <- sqrt(sum((ahatBP-a)^2))/sqrt(sum(a^2))
    }
}


filename = paste('Method2-normerr',toString(d),'.pdf', sep = '')
pdf(filename)
plot(1:d,apply(normerrOMP,1,mean), xlab = "# of non-sparse elements", ylab = "Relative Normed Error", type = 'l', col = 2)
lines(1:d,apply(normerrMP,1,mean), xlab = "# of non-sparse elements", ylab = "Relative Normed Error", type = 'l', col = 3)
lines(1:d,apply(normerrBP,1,mean), xlab = "# of non-sparse elements", ylab = "Relative Normed Error", type = 'l', col = 4)
legend("bottomright",c("OMP","MP","BP"),lty = rep(1,3),col = 2:4)
dev.off()
}


