library("EBImage")
library("dtt")


N = 16
vni <- function(n,i,N){
    if(n==0){out <- ifelse(i<=N,sqrt(1/N),NULL)
        } else {
            out <- ifelse(i<=N,sqrt(2/N)*cos(pi*n*(i-1/2)/N),NULL)
        }
    out
}
sum(vni(1,1:16,N)*vni(1,1:16,N))

sum(vni(1,1:16,N)*vni(2,1:16,N))

size = 5

flower <- channel((readImage("flower.jpg")),"gray")
flower2 <- channel((readImage("flower2.jpg")),"gray")
flower <- resize(flower,size,size)
flower2 <- resize(flower2,size,size)
display(flower)
display(flower2)
###angle between flower and flower2
acos(sum(flower*flower2)/(sqrt(sum(flower*flower))*sqrt(sum(flower2*flower2))))


sum(outer(vni(1,1:16,16),vni(1,1:16,16))*outer(vni(1,1:16,16),vni(1,1:16,16)))

kk = 5
Aflower = matrix(0,nrow=5,ncol=5)
for(n in 1:kk){
    for(m in 1:kk){
        Aflower[m,n] = sum(outer(vni(m,1:size,size),vni(n,1:size,size))*flower)
    }
}

### now we do a regression log(amn) = -b log n + log C to get estimates for b and log c. This is mostly to get a good estimate for -b.

Aflowervec = matrix(Aflower,nrow = 25)

an = sort(abs(Aflowervec), decreasing = TRUE)

reg1 <- lm(log(an) ~ log(1:25))

### estimate for -b is -16.01. As b>1/2, set b = 16.

C = max(an*(n^(16)))







