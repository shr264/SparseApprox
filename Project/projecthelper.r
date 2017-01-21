g <- function(betanew,y,X){(1/2)*t(y-X%*%betanew)%*%(y-X%*%betanew)}
h <- function(betanew,lambda){lambda*sum(abs(betanew))}
subgradh <- function(beta,lambda){ifelse(abs(beta)>lambda,sign(lambda),0)}
gradg <-function(y,X,beta){t(-X)%*%(y-X%*%beta)}
F <- function(betanew,y,X,lambda){g(betanew,y,X) + h(betanew,lambda)}
Q <- function(betanew,betaold,y,X,lambda,L){g(betaold,y,X) + t(betanew-betaold)%*%gradg(y,X,betaold) + (L/2)*(t(betaold-betanew)%*%(betaold-betanew)) + h(betanew,lambda)}
subgF <- function(beta,y,X,lambda){t(-X)%*%(y-X%*%beta) + lambda*subgradh(betaold,1)}
softhresh <- function(x,lambda){sign(x)*pmax(abs(x)-lambda,0)}

subgradientconst <- function(y,X,betaold,lambda,eps=10^(-3),t, betastar){
    betanew <- betaold - t*(-t(X)%*%(y-X%*%betaold) + lambda*subgradh(betaold,1))
    itr <- 1
    err <- max(abs(betanew-betaold))
    err2 <- sqrt(sum((betastar-betanew)^2))
    while(max(abs(betanew-betaold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        betanew <- betaold - t*(-t(X)%*%(y-X%*%betaold) + lambda*subgradh(betaold,1))
        err <- append(err,max(abs(betanew-betaold)))
        err2 <- append(err2,sqrt(sum((betastar-betanew)^2)))
    }
    return(list(betahat = betanew,error = err, error2 = err2, diff = sqrt(sum((betastar-betanew)^2))))
}

subgradientconst2 <- function(y,X,betaold,lambda,eps=10^(-3),t, betastar){
    betanew <- betaold - t*(-t(X)%*%(y-X%*%betaold) + lambda*subgradh(betaold,1))
    itr <- 1
    err <- max(abs(betanew-betaold))
    err2 <- sqrt(sum((betastar-betanew)^2))
    while((itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        betanew <- betaold - t*(-t(X)%*%(y-X%*%betaold) + lambda*subgradh(betaold,1))
        err <- append(err,max(abs(betanew-betaold)))
        err2 <- append(err2,sqrt(sum((betastar-betanew)^2)))
    }
    return(list(betahat = betanew,error = err, error2 = err2, diff = sqrt(sum((betastar-betanew)^2))))
}

subgradientdim <- function(y,X,eta,betaold,lambda,eps=10^(-3),initialt, betastar){
    t = initialt
    while(F(betaold-t*subgF(betaold,y,X,lambda),y,X,lambda)>F(betaold,y,X,lambda)+(t/2)*t(subgF(betaold,y,X,lambda))%*%subgF(betaold,y,X,lambda)){t = eta*t}
    betanew <- betaold - t*(-t(X)%*%(y-X%*%betaold) + lambda*subgradh(betaold))
    itr <- 1
    err <- max(abs(betanew-betaold))
    while(max(abs(betanew-betaold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        t <- initialt
    while(F(betaold-t*subgF(betaold,y,X,lambda),y,X,lambda)>F(betaold,y,X,lambda)+(t/2)*t(subgF(betaold,y,X,lambda))%*%subgF(betaold,y,X,lambda)){t = eta*t}                                                                                                                                     
        betanew <- betaold - t*(-t(X)%*%(y-X%*%betaold) + lambda*subgradh(betaold,1))
        err <- append(err,max(abs(betanew-betaold)))
    }
    return(list(betahat = betanew,error = err, error2 = err2, diff = sqrt(sum(betastar-betanew)^2)))
}

subgradientdim2 <- function(y,X,eta,betaold,lambda,eps=10^(-3),initialt, betastar){
    t = initialt
    while(F(betaold-t*subgF(betaold,y,X,lambda),y,X,lambda)>F(betaold,y,X,lambda)+(t/2)*t(subgF(betaold,y,X,lambda))%*%subgF(betaold,y,X,lambda)){t = eta*t}
    betanew <- betaold - t*(-t(X)%*%(y-X%*%betaold) + lambda*subgradh(betaold,1))
    itr <- 1
    err <- max(abs(betanew-betaold))
    while((itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        t <- initialt
        while(F(betaold-t*subgF(betaold,y,X,lambda),y,X,lambda)>F(betaold,y,X,lambda)+(t/2)*t(subgF(betaold,y,X,lambda))%*%subgF(betaold,y,X,lambda)){t = eta*t}                                                                                                                                     
        betanew <- betaold - t*(-t(X)%*%(y-X%*%betaold) + lambda*subgradh(betaold,1))
        err <- append(err,max(abs(betanew-betaold)))
    }
    return(list(betahat = betanew,error = err, diff = sqrt(sum(betastar-betanew)^2)))
}

ista <- function(y,X,lambda,betaold, eps = 10^(-3),betastar){
    eigenvalues <- eigen(t(X)%*%X)$values
    t <- 1/max(eigenvalues)
    betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t)
    err <- max(abs(betanew-betaold))
    err2 <- sqrt(sum((betastar-betanew)^2))
    itr <- 1
    while(max(abs(betanew-betaold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t)
        err <- append(err,max(abs(betanew-betaold)))
        err2 <- append(err2,sqrt(sum((betastar-betanew)^2)))
    }
    return(list(betahat = betanew,error = err, error2 = err2, diff = sqrt(sum((betastar-betanew)^2))))
}

ista2 <- function(y,X,lambda,betaold, eps = 10^(-3),betastar){
    eigenvalues <- eigen(t(X)%*%X)$values
    t <- 1/max(eigenvalues)
    betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t)
    err <- max(abs(betanew-betaold))
    err2 <- sqrt(sum((betastar-betanew)^2))
    itr <- 1
    while(itr<100){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t)
        err <- append(err,max(abs(betanew-betaold)))
        err2 <- append(err2,sqrt(sum((betastar-betanew)^2)))
    }
    return(list(betahat = betanew,error = err, error2 = err2, diff = sqrt(sum((betastar-betanew)^2))))
}

istadim <- function(y,X,lambda,betaold, eps = 10^(-3),initialL,eta,betastar){
    (L <- initialL)
    (t <- 1/L)
    (betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t))
    while(F(betanew,y,X,lambda)>Q(betanew,betaold,y,X,lambda,(1/t))){
        cat('... decreasing step-size\n')
        (t <- eta*t)
        (betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t))
    }
    (err <- max(abs(betanew-betaold)))
    (itr <- 1)
    while(max(abs(betanew-betaold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        while(F(betanew,y,X,lambda)>Q(betanew,betaold,y,X,lambda,(1/t))){
            cat('... decreasing step-size\n')
            t <- eta*t
            betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t)
        }
        err <- append(err,max(abs(betanew-betaold)))
    }
    return(list(betahat = betanew,error = err, diff = sqrt(sum(betastar-betanew)^2)))
}

istadim2 <- function(y,X,lambda,betaold, eps = 10^(-3),initialL,eta,betastar){
    (L <- initialL)
    (t <- 1/L)
    (betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t))
    while(F(betanew,y,X,lambda)>Q(betanew,betaold,y,X,lambda,(1/t))){
        cat('... decreasing step-size\n')
        (t <- eta*t)
        (betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t))
    }
    (err <- max(abs(betanew-betaold)))
    (itr <- 1)
    while((itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        while(F(betanew,y,X,lambda)>Q(betanew,betaold,y,X,lambda,(1/t))){
            cat('... decreasing step-size\n')
            t <- eta*t
            betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t)
        }
        err <- append(err,max(abs(betanew-betaold)))
    }
    return(list(betahat = betanew,error = err, diff = sqrt(sum(betastar-betanew)^2)))
}

fista <- function(y,X,lambda,betaold, eps = 10^(-3),betastar){
    eigenvalues <- eigen(t(X)%*%X)$values
    t <- 1/max(eigenvalues)
    itr <- 1
    zeta0 <- betaold
    betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
    zeta0 <- betanew + ((itr-1)/(itr+2))*(betanew-betaold)
    err <- max(abs(betanew-betaold))
    err2 <- sqrt(sum((betastar-betanew)^2))
    while(max(abs(betanew-betaold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
        zeta0 <- betanew + ((itr-1)/(itr+2))*(betanew-betaold)
        err <- append(err,max(abs(betanew-betaold)))
        err2 <- append(err2,sqrt(sum((betastar-betanew)^2)))
    }
    return(list(betahat = betanew,error = err, error2 = err2, diff = sqrt(sum((betastar-betanew)^2))))
}

fista2 <- function(y,X,lambda,betaold, eps = 10^(-3),betastar){
    eigenvalues <- eigen(t(X)%*%X)$values
    t <- 1/max(eigenvalues)
    itr <- 1
    zeta0 <- betaold
    betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
    zeta0 <- betanew + ((itr-1)/(itr+2))*(betanew-betaold)
    err <- max(abs(betanew-betaold))
    err2 <- sqrt(sum((betastar-betanew)^2))
    while((itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
        zeta0 <- betanew + ((itr-1)/(itr+2))*(betanew-betaold)
        err <- append(err,max(abs(betanew-betaold)))
        err2 <- append(err2,sqrt(sum((betastar-betanew)^2)))
    }
    return(list(betahat = betanew,error = err, error2 = err2, diff = sqrt(sum((betastar-betanew)^2))))
}

fistadim <- function(y,X,lambda,betaold, eps = 10^(-3),initialL,eta,betastar){
    L <- initialL
    t <- 1/L
    itr <- 1
    zeta0 <- betaold
    betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
    while(F(betanew,y,X,lambda)>Q(betanew,betaold,y,X,lambda,(1/t))){
        t <- eta*t
        betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
    }
    zeta0 <- betanew + ((itr-1)/(itr+2))*(betanew-betaold)
    err <- max(abs(betanew-betaold))
    while(max(abs(betanew-betaold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        while(F(betanew,y,X,lambda)>Q(betanew,betaold,y,X,lambda,(1/t))){
            t <- eta*t
            betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
        }
        zeta0 <- betanew + ((itr-1)/(itr+2))*(betanew-betaold)
        err <- append(err,max(abs(betanew-betaold)))
    }
    return(list(betahat = betanew,error = err, diff = sqrt(sum(betastar-betanew)^2)))
}

fistadim2 <- function(y,X,lambda,betaold, eps = 10^(-3),initialL,eta,betastar){
    L <- initialL
    t <- 1/L
    itr <- 1
    zeta0 <- betaold
    betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
    while(F(betanew,y,X,lambda)>Q(betanew,betaold,y,X,lambda,(1/t))){
        t <- eta*t
        betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
    }
    zeta0 <- betanew + ((itr-1)/(itr+2))*(betanew-betaold)
    err <- max(abs(betanew-betaold))
    while((itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        while(F(betanew,y,X,lambda)>Q(betanew,betaold,y,X,lambda,(1/t))){
            t <- eta*t
            betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
        }
        zeta0 <- betanew + ((itr-1)/(itr+2))*(betanew-betaold)
        err <- append(err,max(abs(betanew-betaold)))
    }
    return(list(betahat = betanew,error = err, diff = sqrt(sum(betastar-betanew)^2)))
}

admm <- function(y, X, lambda, betaold, eps = 10^(-3), rho, betastar){
    p <- length(betaold) 
    gammaold <- etaold <- betaold
    betanew <- solve(t(X)%*%X + rho*diag(p))%*%(rho*gammaold - etaold + t(X)%*%y)
    gammanew <- softhresh(betanew + etaold/rho,lambda/rho)
    etanew <- etaold + rho*(betanew-gammanew)
    itr <-  1
    err <-max(abs(betanew-betaold))
    err2 <- sqrt(sum((betastar-betanew)^2))
    while(max(abs(betanew-betaold)>eps)&&(itr<100)){
        gammaold <- gammanew
        etaold <- etanew
        betaold <- betanew
        betanew <- solve(t(X)%*%X + rho*diag(p))%*%(rho*gammaold - etaold + t(X)%*%y)
        gammanew <- softhresh(betanew + etaold/rho,lambda/rho)
        etanew <- etaold + rho*(betanew-gammanew)
        err <- append(err,max(abs(betanew-betaold)))
        err2 <- append(err2,sqrt(sum((betastar-betanew)^2)))
        itr <- itr + 1
    }
    return(list(betahat = betanew,error = err, error2 = err2, diff = sqrt(sum((betastar-betanew)^2))))
}

admm2 <- function(y, X, lambda, betaold, eps = 10^(-3), rho, betastar){
    p <- length(betaold) 
    gammaold <- etaold <- betaold
    betanew <- solve(t(X)%*%X + rho*diag(p))%*%(rho*gammaold - etaold + t(X)%*%y)
    gammanew <- softhresh(betanew + etaold/rho,lambda/rho)
    etanew <- etaold + rho*(betanew-gammanew)
    itr <-  1
    err <-max(abs(betanew-betaold))
    err2 <- sqrt(sum((betastar-betanew)^2))
    while((itr<100)){
        gammaold <- gammanew
        etaold <- etanew
        betaold <- betanew
        betanew <- solve(t(X)%*%X + rho*diag(p))%*%(rho*gammaold - etaold + t(X)%*%y)
        gammanew <- softhresh(betanew + etaold/rho,lambda/rho)
        etanew <- etaold + rho*(betanew-gammanew)
        err <- append(err,max(abs(betanew-betaold)))
        err2 <- append(err2,sqrt(sum((betastar-betanew)^2)))
        itr <- itr + 1
    }
    return(list(betahat = betanew,error = err, error2 = err2, diff = sqrt(sum((betastar-betanew)^2))))
}

psnr <- function(original,noisy){
    MSE <- sum((original-noisy)^2)/(dim(original)[1]*dim(original)[2])
    MAX <- max(original)
    out <- 10*log10((MAX^2)/MSE)
    out
}


cv5sg <- function(y,X,betaold,eps,betastar,nlambda){
    (foldsize <- floor(length(y)/5))   
    (err <- matrix(0,nrow = length(nlambda), ncol = 5))
    (rndmsmpl <- sample(1:length(y)))
    (folds <- split(rndmsmpl, ceiling(seq_along(rndmsmpl)/foldsize)))
    for(lambda in nlambda){
        for(i in 1:5){
            #i = 1
            (a1 <- folds[[i]])
            (fit <- subgradientconst2(y[-folds[[i]]],X[-folds[[i]],],betaold,lambda,eps,0.0001,betastar))
            (err[lambda,i] <- sum((a1 - X[folds[[i]],]%*%fit$betahat)^2))
        }
    }
    nlambda[which.min(apply(err,1,mean))]
}

cv5ista <- function(y,X,betaold,eps,betastar,nlambda){
    (foldsize <- floor(length(y)/5))   
    (err <- matrix(0,nrow = length(nlambda), ncol = 5))
    (rndmsmpl <- sample(1:length(y)))
    (folds <- split(rndmsmpl, ceiling(seq_along(rndmsmpl)/foldsize)))
    for(lambda in nlambda){
        for(i in 1:5){
            #i = 1
            (a1 <- folds[[i]])
            (fit <- ista2(y[-folds[[i]]],X[-folds[[i]],],lambda,betaold,eps,betastar))
            (err[lambda,i] <- sum((a1 - X[folds[[i]],]%*%fit$betahat)^2))
        }
    }
    nlambda[which.min(apply(err,1,mean))]
}

cv5fista <- function(y,X,betaold,eps,betastar,nlambda){
    (foldsize <- floor(length(y)/5))   
    (err <- matrix(0,nrow = length(nlambda), ncol = 5))
    (rndmsmpl <- sample(1:length(y)))
    (folds <- split(rndmsmpl, ceiling(seq_along(rndmsmpl)/foldsize)))
    for(lambda in nlambda){
        for(i in 1:5){
            #i = 1
            (a1 <- folds[[i]])
            (fit <- fista2(y[-folds[[i]]],X[-folds[[i]],],lambda,betaold,eps,betastar))
            (err[lambda,i] <- sum((a1 - X[folds[[i]],]%*%fit$betahat)^2))
        }
    }
    nlambda[which.min(apply(err,1,mean))]
}


cv5admm <- function(y,X,betaold,eps,betastar,rho, nlambda){
    (foldsize <- floor(length(y)/5))   
    (err <- matrix(0,nrow = length(nlambda), ncol = 5))
    (rndmsmpl <- sample(1:length(y)))
    (folds <- split(rndmsmpl, ceiling(seq_along(rndmsmpl)/foldsize)))
    for(lambda in nlambda){
        for(i in 1:5){
            #i = 1
            (a1 <- folds[[i]])
            rho <- 1
            (fit <- admm(y[-folds[[i]]],X[-folds[[i]],],lambda, betaold, eps = 10^(-3), rho, betastar))
            (err[lambda,i] <- sum((a1 - X[folds[[i]],]%*%fit$betahat)^2))
        }
    }
    nlambda[which.min(apply(err,1,mean))]
}

glassoadmm <- function(Y,xold,rho,lambda,eps=10^(-3),omegastar){
    n <- dim(Y)[1]
    p <- dim(Y)[2]
    S <- (1/n)*t(Y)%*%Y
    uold <- zold <- xold
    eig <- eigen(rho*(zold-uold)-S)
    xdiag <- (1/(2*rho))*(lambda + sqrt(lambda^2+ 4*rho))
    xnew <- Re(eig$vectors%*%(xdiag*diag(p))%*%t(eig$vectors))
    znew <- xnew + uold
    znew[row(znew)!=col(znew)] <- softhresh(znew[row(znew)!=col(znew)],lambda/rho)
    unew <- uold + (xnew-znew)
    itr <-  1
    err <-max(abs(xnew-xold))
    while(max(abs(xnew-xold)>eps)&&(itr<100)){
        cat('... itr:',itr,'\n')
        xold <- xnew
        zold <- znew
        uold <- unew
        eig <- eigen(rho*(zold-uold)-S)
        xdiag <- (1/(2*rho))*(lambda + sqrt(lambda^2+ 4*rho))
        xnew <- Re(eig$vectors%*%(xdiag*diag(p))%*%t(eig$vectors))
        znew <- xnew + uold
        znew[row(znew)!=col(znew)] <- softhresh(znew[row(znew)!=col(znew)],lambda/rho)
        unew <- uold + (xnew-znew)
        err <- append(err,max(abs(xnew-xold)))
        itr <- itr + 1
    }
    return(list(omegahat = xnew,error = err, diff = sqrt(sum(omegastar-xnew)^2)))
}

glassosubg <- function(Y,xold,t,lambda,eps=10^(-3),omegastar){
n <- dim(Y)[1]
p <- dim(Y)[2]
S <- (1/n)*t(Y)%*%Y
sgupdate <- matrix(0,p,p)
sgupdate[row(sgupdate)!=col(sgupdate)] <- subgradh(xold[row(xold)!=col(xold)],1)
xnew <- xold - t*(solve(xold)-S + lambda*sgupdate)
itr <- 1
err <- max(abs(xnew-xold))
while(max(abs(xnew-xold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        xold <- xnew
        sgupdate <- matrix(0,p,p)
        sgupdate[row(sgupdate)!=col(sgupdate)] <- subgradh(xold[row(xold)!=col(xold)],1)
        xnew <- xold - t*(solve(xold)-S + lambda*sgupdate)
        err <- append(err,max(abs(xnew-xold)))
    }
    return(list(omegahat = xnew,error = err,diff = sqrt(sum(omegastar-xnew)^2)))
}

glassoista <- function(Y,xold,t,lambda,eps=10^(-3),omegastar){
    n <- dim(Y)[1]
    p <- dim(Y)[2]
    S <- (1/n)*t(Y)%*%Y
    xnew <- xold - t*(solve(xold)-S)
    xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
    err <- max(abs(xnew-xold))
    itr <- 1
    while(max(abs(xnew-xold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        xold <- xnew
        xnew <- xold - t*(solve(xold)-S)
        xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
        err <- append(err,max(abs(xnew-xold)))
    }
return(list(omegahat = xnew,error = err,diff = sqrt(sum(omegastar-xnew)^2)))
}

glassofista <- function(Y,xold,t,lambda,eps=10^(-3),omegastar){
    n <- dim(Y)[1]
    p <- dim(Y)[2]
    S <- (1/n)*t(Y)%*%Y
    itr <- 1
    zeta0 <- xold
    xnew <- zeta0 - t*(solve(zeta0)-S)
    xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
    zeta0 <- xnew + ((itr-1)/(itr+2))*(xnew-xold)
    err <- max(abs(xnew-xold))
    while(max(abs(xnew-xold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        xold <- xnew
        zeta0 <- xold
        xnew <- zeta0 - t*(solve(zeta0)-S)
        xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
        zeta0 <- xnew + ((itr-1)/(itr+2))*(xnew-xold)
        err <- append(err,max(abs(xnew-xold)))
    }
    return(list(omegahat = xnew,error = err,diff = sqrt(sum(omegastar-xnew)^2)))
}


glassosubg2 <- function(Y,xold,t,lambda,eps=10^(-3),omegastar){
n <- dim(Y)[1]
p <- dim(Y)[2]
S <- (1/n)*t(Y)%*%Y
sgupdate <- matrix(0,p,p)
sgupdate[row(sgupdate)!=col(sgupdate)] <- subgradh(xold[row(xold)!=col(xold)],1)
xnew <- xold - t*(-solve(xold) + S - lambda*sgupdate)
itr <- 1
err <- max(abs(xnew-xold))
while(max(abs(xnew-xold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        xold <- xnew
        sgupdate <- matrix(0,p,p)
        sgupdate[row(sgupdate)!=col(sgupdate)] <- subgradh(xold[row(xold)!=col(xold)],1)
        xnew <- xold - t*( - solve(xold) + S - lambda*sgupdate)
        err <- append(err,max(abs(xnew-xold)))
    }
    return(list(omegahat = xnew,error = err,diff = sqrt(sum(omegastar-xnew)^2)))
}

glassoista2 <- function(Y,xold,t,lambda,eps=10^(-3),omegastar){
    n <- dim(Y)[1]
    p <- dim(Y)[2]
    S <- (1/n)*t(Y)%*%Y
    xnew <- xold + t*(-solve(xold)+S)
    xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
    err <- max(abs(xnew-xold))
    itr <- 1
    while(max(abs(xnew-xold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        xold <- xnew
        xnew <- xold + t*(-solve(xold)+S)
        xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
        err <- append(err,max(abs(xnew-xold)))
    }
return(list(omegahat = xnew,error = err,diff = sqrt(sum(omegastar-xnew)^2)))
}

glassofista2 <- function(Y,xold,t,lambda,eps=10^(-3),omegastar){
    n <- dim(Y)[1]
    p <- dim(Y)[2]
    S <- (1/n)*t(Y)%*%Y
    itr <- 1
    zeta0 <- xold
    xnew <- zeta0 + t*(-solve(zeta0)+S)
    xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
    zeta0 <- xnew + ((itr-1)/(itr+2))*(xnew-xold)
    err <- max(abs(xnew-xold))
    while(max(abs(xnew-xold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        xold <- xnew
        zeta0 <- xold
        xnew <- zeta0 + t*(-solve(zeta0)+S)
        xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
        zeta0 <- xnew + ((itr-1)/(itr+2))*(xnew-xold)
        err <- append(err,max(abs(xnew-xold)))
    }
    return(list(omegahat = xnew,error = err,diff = sqrt(sum(omegastar-xnew)^2)))
}

glassosubg3 <- function(Y,xold,t,lambda,eps=10^(-3),omegastar){
n <- dim(Y)[1]
p <- dim(Y)[2]
S <- (1/n)*t(Y)%*%Y
sgupdate <- matrix(0,p,p)
sgupdate[row(sgupdate)!=col(sgupdate)] <- subgradh(xold[row(xold)!=col(xold)],1)
xnew <- xold + t*(solve(xold)-S + lambda*sgupdate)
itr <- 1
err <- max(abs(xnew-xold))
while(max(abs(xnew-xold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        xold <- xnew
        sgupdate <- matrix(0,p,p)
        sgupdate[row(sgupdate)!=col(sgupdate)] <- subgradh(xold[row(xold)!=col(xold)],1)
        xnew <- xold + t*(solve(xold)-S + lambda*sgupdate)
        err <- append(err,max(abs(xnew-xold)))
    }
    return(list(omegahat = xnew,error = err,diff = sqrt(sum(omegastar-xnew)^2)))
}

glassoista3 <- function(Y,xold,t,lambda,eps=10^(-3),omegastar){
    n <- dim(Y)[1]
    p <- dim(Y)[2]
    S <- (1/n)*t(Y)%*%Y
    xnew <- xold + t*(solve(xold)-S)
    xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
    err <- max(abs(xnew-xold))
    itr <- 1
    while(max(abs(xnew-xold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        xold <- xnew
        xnew <- xold + t*(solve(xold)-S)
        xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
        err <- append(err,max(abs(xnew-xold)))
    }
return(list(omegahat = xnew,error = err,diff = sqrt(sum(omegastar-xnew)^2)))
}

glassofista3 <- function(Y,xold,t,lambda,eps=10^(-3),omegastar){
    n <- dim(Y)[1]
    p <- dim(Y)[2]
    S <- (1/n)*t(Y)%*%Y
    itr <- 1
    zeta0 <- xold
    xnew <- zeta0 + t*(solve(zeta0)-S)
    xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
    zeta0 <- xnew + ((itr-1)/(itr+2))*(xnew-xold)
    err <- max(abs(xnew-xold))
    while(max(abs(xnew-xold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        xold <- xnew
        zeta0 <- xold
        xnew <- zeta0 + t*(solve(zeta0)-S)
        xnew[row(xnew)!=col(xnew)] <- softhresh(xnew[row(xnew)!=col(xnew)],lambda*t)
        zeta0 <- xnew + ((itr-1)/(itr+2))*(xnew-xold)
        err <- append(err,max(abs(xnew-xold)))
    }
    return(list(omegahat = xnew,error = err,diff = sqrt(sum(omegastar-xnew)^2)))
}

tprfpr = function(path, theta, verbose = TRUE){
	gcinfo(verbose = FALSE)	
	ROC = list()
	
	theta = as.matrix(theta)
	d = ncol(theta)
	pos.total = sum(theta!=0)
	neg.total = d*(d-1) - pos.total
	
	if(verbose) cat("Computing F1 scores, false positive rates and true positive rates....")
	ROC$tp = rep(0,length(path))
   	ROC$fp = rep(0,length(path))
   	ROC$F1 = rep(0,length(path))
   	for (r in 1:length(path)){
   		tmp = as.matrix(path[[r]]) 
   		tp.all = (theta!=0)*(tmp!=0)
   		diag(tp.all) = 0
		ROC$tp[r] <- sum(tp.all!=0)/pos.total
		fp.all = (theta==0)*(tmp!=0)
		diag(fp.all) = 0
		ROC$fp[r] <- sum(fp.all!=0)/neg.total
		
		fn = 1 - ROC$tp[r]
		precision = ROC$tp[r]/(ROC$tp[r]+ROC$fp[r])
		recall = ROC$tp[r]/(ROC$tp[r]+fn)
		ROC$F1[r] = 2*precision*recall/(precision+recall)
		if(is.na(ROC$F1[r]))	ROC$F1[r] = 0
	}
	if(verbose) cat("done.\n")
		
	rm(precision,recall,tp.all,fp.all,path,theta,fn)
   	gc()	
		
	ord.fp = order(ROC$fp)
	
	tmp1 = ROC$fp[ord.fp]
	tmp2 = ROC$tp[ord.fp]
        return(list(fp = tmp1, tp = tmp2))
}
