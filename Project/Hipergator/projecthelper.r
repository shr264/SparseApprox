g <- function(betanew,y,X){(1/2)*t(y-X%*%betanew)%*%(y-X%*%betanew)}
h <- function(betanew,lambda){lambda*sum(abs(betanew))}
subgradh <- function(beta,lambda){ifelse(abs(beta)>lambda,sign(lambda),0)}
gradg <-function(y,X,beta){t(-X)%*%(y-X%*%beta)}
F <- function(betanew,y,X,lambda){g(betanew,y,X) + h(betanew,lambda)}
Q <- function(betanew,betaold,y,X,lambda,L){g(betaold,y,X) + t(betanew-betaold)%*%gradg(y,X,betaold) + (L/2)*(t(betaold-betanew)%*%(betaold-betanew)) + h(betanew,lambda)}
subgF <- function(beta,y,X,lambda){t(-X)%*%(y-X%*%beta) + lambda*subgradh(betaold,1)}

subgradientconst <- function(y,X,betaold,lambda,eps=10^(-3),t, betastar){
    betanew <- betaold - t*(-t(X)%*%(y-X%*%betaold) + lambda*subgradh(betaold,1))
    itr <- 1
    err <- max(abs(betanew-betaold))
    while(max(abs(betanew-betaold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        betanew <- betaold - t*(-t(X)%*%(y-X%*%betaold) + lambda*subgradh(betaold,1))
        err <- append(err,max(abs(betanew-betaold)))
    }
    return(list(betahat = betanew,error = err,diff = sqrt(sum(betastar-betanew)^2)))
}

subgradientconst2 <- function(y,X,betaold,lambda,eps=10^(-3),t, betastar){
    betanew <- betaold - t*(-t(X)%*%(y-X%*%betaold) + lambda*subgradh(betaold,1))
    itr <- 1
    err <- max(abs(betanew-betaold))
    while((itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        betanew <- betaold - t*(-t(X)%*%(y-X%*%betaold) + lambda*subgradh(betaold,1))
        err <- append(err,max(abs(betanew-betaold)))
    }
    return(list(betahat = betanew,error = err,diff = sqrt(sum(betastar-betanew)^2)))
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
    return(list(betahat = betanew,error = err, diff = sqrt(sum(betastar-betanew)^2)))
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

softhresh <- function(x,lambda){sign(x)*pmax(abs(x)-lambda,0)}

ista <- function(y,X,lambda,betaold, eps = 10^(-3),betastar){
    eigenvalues <- eigen(t(X)%*%X)$values
    t <- 1/max(eigenvalues)
    betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t)
    err <- max(abs(betanew-betaold))
    itr <- 1
    while(max(abs(betanew-betaold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t)
        err <- append(err,max(abs(betanew-betaold)))
    }
    return(list(betahat = betanew,error = err,diff = sqrt(sum(betastar-betanew)^2)))
}

ista2 <- function(y,X,lambda,betaold, eps = 10^(-3),betastar){
    eigenvalues <- eigen(t(X)%*%X)$values
    t <- 1/max(eigenvalues)
    betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t)
    err <- max(abs(betanew-betaold))
    itr <- 1
    while(itr<100){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        betanew <- softhresh(betaold + t*t(X)%*%(y-X%*%betaold),lambda*t)
        err <- append(err,max(abs(betanew-betaold)))
    }
    return(list(betahat = betanew,error = err,diff = sqrt(sum(betastar-betanew)^2)))
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
    while(max(abs(betanew-betaold)>eps)&&(itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
        zeta0 <- betanew + ((itr-1)/(itr+2))*(betanew-betaold)
        err <- append(err,max(abs(betanew-betaold)))
    }
    return(list(betahat = betanew,error = err,diff = sqrt(sum(betastar-betanew)^2)))
}

fista2 <- function(y,X,lambda,betaold, eps = 10^(-3),betastar){
    eigenvalues <- eigen(t(X)%*%X)$values
    t <- 1/max(eigenvalues)
    itr <- 1
    zeta0 <- betaold
    betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
    zeta0 <- betanew + ((itr-1)/(itr+2))*(betanew-betaold)
    err <- max(abs(betanew-betaold))
    while((itr<100)){
        cat('...',itr,'\n')
        itr <-  itr + 1
        betaold <- betanew
        betanew <- softhresh(zeta0 + t*t(X)%*%(y-X%*%zeta0),lambda*t)
        zeta0 <- betanew + ((itr-1)/(itr+2))*(betanew-betaold)
        err <- append(err,max(abs(betanew-betaold)))
    }
    return(list(betahat = betanew,error = err,diff = sqrt(sum(betastar-betanew)^2)))
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
    while(max(abs(betanew-betaold)>eps)&&(itr<100)){
        gammaold <- gammanew
        etaold <- etanew
        betaold <- betanew
        betanew <- solve(t(X)%*%X + rho*diag(p))%*%(rho*gammaold - etaold + t(X)%*%y)
        gammanew <- softhresh(betanew + etaold/rho,lambda/rho)
        etanew <- etaold + rho*(betanew-gammanew)
        err <- append(err,max(abs(betanew-betaold)))
        itr <- itr + 1
    }
    return(list(betahat = betanew,error = err, diff = sqrt(sum(betastar-betanew)^2)))
}

admm2 <- function(y, X, lambda, betaold, eps = 10^(-3), rho, betastar){
    p <- length(betaold) 
    gammaold <- etaold <- betaold
    betanew <- solve(t(X)%*%X + rho*diag(p))%*%(rho*gammaold - etaold + t(X)%*%y)
    gammanew <- softhresh(betanew + etaold/rho,lambda/rho)
    etanew <- etaold + rho*(betanew-gammanew)
    itr <-  1
    err <-max(abs(betanew-betaold))
    while((itr<100)){
        gammaold <- gammanew
        etaold <- etanew
        betaold <- betanew
        betanew <- solve(t(X)%*%X + rho*diag(p))%*%(rho*gammaold - etaold + t(X)%*%y)
        gammanew <- softhresh(betanew + etaold/rho,lambda/rho)
        etanew <- etaold + rho*(betanew-gammanew)
        err <- append(err,max(abs(betanew-betaold)))
        itr <- itr + 1
    }
    return(list(betahat = betanew,error = err, diff = sqrt(sum(betastar-betanew)^2)))
}

psnr <- function(original,noisy){
    MSE <- sum((original-noisy)^2)/(dim(original)[1]*dim(original)[2])
    MAX <- max(original)
    out <- 10*log10((MAX^2)/MSE)
    out
}
