A <- matrix(0,nrow = 3, ncol = 3)
A[1,1] = integrate(function(t){rep(1,length(t))}, 0, 1)$value
A[2,1] = A[1,2] = integrate(function(t){t}, 0, 1)$value
A[3,1] = A[1,3] = integrate(function(t){(t^2)}, 0, 1)$value
A[3,2] = A[2,3] = integrate(function(t){(t^3)}, 0, 1)$value
A[2,2] = integrate(function(t){(t^2)}, 0, 1)$value
A[3,3] = integrate(function(t){(t^4)}, 0, 1)$value

b <- matrix(0,nrow = 3, ncol = 1)
b[1] = integrate(function(t){(exp(t))}, 0, 1)$value
b[2] = integrate(function(t){(exp(t))*t}, 0, 1)$value
b[3] = integrate(function(t){(exp(t))*(t^2)}, 0, 1)$value

a = solve(A,b)
t = seq(0,1,length=200)
xhat = a[1] + a[2]*t + a[3]*t^2
x = exp(t)
xtaylor = 1 + t + t^2/2
pdf('hw1p1.pdf')
plot(t,x, type = 'l', col = 1, ylab = "y", xlab = "t")
lines(t,xhat, col = 2)
lines(t,xtaylor, col = 3)
legend(0,2.5,c(expression(e^x),expression(hat(x)),expression(x[Taylor])),lty = rep(1,3),col = 1:3)
dev.off()
f <- function(t){
    (exp(t)-(a[1] + a[2]*t + a[3]*(t^2)))^2
}
err <- integrate(f,0,1)
ftaylor <- function(t){
    (exp(t)-(1 + 1*t + (1/2)*(t^2)))^2
}
errtaylor <- integrate(ftaylor,0,1)
errtaylor$value - err$value

A1 <- matrix(0,nrow = 3, ncol = 3)
A1[1,1] = integrate(function(t){16*(t-1/2)^2*1}, 0, 1)$value
A1[2,1] = A1[1,2] = integrate(function(t){16*(t-1/2)^2*t}, 0, 1)$value
A1[3,1] = A1[1,3] = integrate(function(t){16*(t-1/2)^2*(t^2)}, 0, 1)$value
A1[3,2] = A1[2,3] = integrate(function(t){16*(t-1/2)^2*(t^3)}, 0, 1)$value
A1[2,2] = integrate(function(t){16*(t-1/2)^2*(t^2)}, 0, 1)$value
A1[3,3] = integrate(function(t){16*(t-1/2)^2*(t^4)}, 0, 1)$value

b1 <- matrix(0,nrow = 3, ncol = 1)
b1[1] = integrate(function(t){16*(t-1/2)^2*(exp(t))}, 0, 1)$value
b1[2] = integrate(function(t){16*(t-1/2)^2*(exp(t))*t}, 0, 1)$value
b1[3] = integrate(function(t){16*(t-1/2)^2*(exp(t))*(t^2)}, 0, 1)$value

a1 = solve(A1,b1)

t = seq(0,1,length=200)
x1hat = a1[1] + a1[2]*t + a1[3]*t^2

pdf('hw1p1-2.pdf')
plot(t,x, type = 'l', col = 1, ylab = "y", xlab = "t")
lines(t,x1hat, col = 2)
lines(t,xhat, col = 3)
legend(0,2.5,c(expression(e^x),expression(hat(x)[1]),expression(hat(x))),lty = rep(1,3),col = 1:3)
dev.off()

f <- function(t){
    (16*(t-1/2)^2)*(exp(t)-(a1[1] + a1[2]*t + a1[3]*(t^2)))^2
}
err1 <- integrate(f,0,1)
err1$value-err$value

