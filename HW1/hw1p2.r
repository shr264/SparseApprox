library("EBImage")

N = 16
vni <- function(n,i,N){
    if(n=0){out <- ifelse(i<=N,sqrt(1/N),NULL)
        } else {
            out <- ifelse(i<=N,sqrt(2/N)*cos(pi*n*(i-1/2)/N),NULL)
        }
    out
}
sum(vni(1,1:16,N)*vni(1,1:16,N))

sum(vni(1,1:16,N)*vni(2,1:16,N))

size = 100

flower <- channel((readImage("flower.jpg")),"gray")
flower2 <- channel((readImage("flower2.jpg")),"gray")
tiger <- channel((readImage("Tiger.jpg")),"gray")
taj <- channel((readImage("taj.jpg")),"gray")
flower <- resize(flower,size,size)
flower2 <- resize(flower2,size,size)
tiger <- resize(tiger,size,size)
tiger <- resize(tiger,size,size)
display(flower)
display(flower2)
display(tiger)
###angle between flower and flower2
acos(sum(flower*flower2)/(sqrt(sum(flower*flower))*sqrt(sum(flower2*flower2))))
###angle between flower and tiger
acos(sum(flower*tiger)/(sqrt(sum(flower*flower))*sqrt(sum(tiger*tiger))))

###angle between flower2 and tiger
acos(sum(flower2*tiger)/(sqrt(sum(flower2*flower2))*sqrt(sum(tiger*tiger))))

###angle between flower2 and taj
acos(sum(flower2*taj)/(sqrt(sum(flower2*flower2))*sqrt(sum(taj*taj))))

###angle between flower and taj
acos(sum(flower*taj)/(sqrt(sum(flower*flower))*sqrt(sum(taj*taj))))

###angle between taj and tiger
acos(sum(taj*tiger)/(sqrt(sum(taj*taj))*sqrt(sum(tiger*tiger))))


sum(outer(vni(1,1:16,16),vni(1,1:16,16))*outer(vni(1,1:16,16),vni(1,1:16,16)))

plot(1:16,outer(vni(0,1:16,16),vni(0,1:16,16))[1,], type = 'l',col = 1,ylim = c(-0.3,0.3),xlab = " ", ylab = " ")
for(m in 0:3){
    for(n in 0:3){
        lines(outer(vni(m,1:16,16),vni(n,1:16,16))[1,], col = m+n)
    }
}

### 25% energy means k = 50
numberk = 10
norm = rep(0,numberk)
k = floor(seq(20,100,length=numberk))
for (kk in 1:numberk){
    Aflower = matrix(0,nrow=k[kk],ncol=k[kk])
    for(n in 1:k[kk]){
        for(m in 1:k[kk]){
            Aflower[m,n] = sum(outer(vni(m,1:size,size),vni(n,1:size,size))*flower)
        }
    }

    xhatflower = matrix(0,nrow=size,ncol=size)
    for(n in 1:k[kk]){
        for(m in 1:k[kk]){
            xhatflower = xhatflower + Aflower[m,n]*outer(vni(m,1:size,size),vni(n,1:size,size))
        }
    }
    display(channel(xhatflower,"gray"))
    norm[kk] = sqrt(sum((flower - xhatflower)^2))
}
pdf('norm.pdf')
plot(k,norm, xlab = "# of basis elements", ylab = "Normed Error", type = 'l', col = 2)
dev.off()

### 25% energy means k = 50
numberk = 10
norm = rep(0,numberk)
k = floor(seq(20,100,length=numberk))
for (kk in 1:numberk){
    Aflower2 = matrix(0,nrow=k[kk],ncol=k[kk])
    for(n in 1:k[kk]){
        for(m in 1:k[kk]){
            Aflower2[m,n] = sum(outer(vni(m,1:size,size),vni(n,1:size,size))*flower2)
        }
    }

    xhatflower2 = matrix(0,nrow=size,ncol=size)
    for(n in 1:k[kk]){
        for(m in 1:k[kk]){
            xhatflower2 = xhatflower2 + Aflower2[m,n]*outer(vni(m,1:size,size),vni(n,1:size,size))
        }
    }
    display(channel(xhatflower2,"gray"))
    norm[kk] = sqrt(sum((flower2 - xhatflower2)^2))
}
pdf('norm2.pdf')
plot(k,norm, xlab = "# of basis elements", ylab = "Normed Error", type = 'l', col = 2)
dev.off()

### 25% energy means k = 50
numberk = 10
norm = rep(0,numberk)
k = floor(seq(20,100,length=numberk))
for (kk in 1:numberk){
    Atiger = matrix(0,nrow=k[kk],ncol=k[kk])
    for(n in 1:k[kk]){
        for(m in 1:k[kk]){
            Atiger[m,n] = sum(outer(vni(m,1:size,size),vni(n,1:size,size))*tiger)
        }
    }

    xhattiger = matrix(0,nrow=size,ncol=size)
    for(n in 1:k[kk]){
        for(m in 1:k[kk]){
            xhattiger = xhattiger + Atiger[m,n]*outer(vni(m,1:size,size),vni(n,1:size,size))
        }
    }
    display(channel(xhattiger,"gray"))
    norm[kk] = sqrt(sum((tiger - xhattiger)^2))
}
pdf('norm3.pdf')
plot(k,norm, xlab = "# of basis elements", ylab = "Normed Error", type = 'l', col = 2)
dev.off()
