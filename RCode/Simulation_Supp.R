#This file contains the simulation results in Figure 1 of our online Supplementary Material.

library(VGAM)
library(ggplot2)


###S1, gamma
delta = c(2:5)/20
p=c(60)
n=c(40,50,60)
j=1
k=1
i=1
m.prop=matrix(nrow=200,ncol=length(delta))
m.sing=matrix(nrow=200,ncol=length(delta))
for(i in 1:length(delta)){
  #  for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    theta = log(1+c(1:p[j])*delta[i])
    b = runif(n[k],1,3)
    
    Theta = matrix(c(rep(theta,n[k]/2),rep(0,(n[k]/2)*p[j])),nrow=n[k],byrow=T)+
      b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,i]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.sing[rounds,i]=max(kendall.tau(1:p[j],order(order(svd(data)$v[,1]))),
                        kendall.tau(p[j]:1,order(order(svd(data)$v[,1]))))
  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.sing)),
                   alpha= rep(rep(as.factor(delta),each=200),2),
                   Method = rep(c("pi.hat","pi.svd"),each=200*length(delta)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=alpha, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S1: n=40, p=60")+scale_fill_grey(start = 0.3, end = 0.8)
p


######p
i=1
j=1
k=1
delta = c(2:5)/20
p=c(45,60,75,90)
n=c(40,50,60)
m.prop=matrix(nrow=200,ncol=4)
m.sing=matrix(nrow=200,ncol=4)
#for(i in 1:length(delta)){
for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    theta = log(1+c(1:p[j])*delta[i])
    b = runif(n[k],1,3)
    c = runif(n[k],1,3)
    Theta = matrix(c(rep(theta,n[k]/2),rep(0,(n[k]/2)*p[j])),nrow=n[k],byrow=T)+ 
      b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,j]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.sing[rounds,j]=max(kendall.tau(1:p[j],order(order(svd(data)$v[,1]))),
                        kendall.tau(p[j]:1,order(order(svd(data)$v[,1]))))
  }
  #    }
  #  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.sing)),
                   p = rep(rep(as.factor(p),each=200),2),
                   Method = rep(c("pi.hat","pi.svd"),each=200*length(p)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=p, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S1: n=40, alpha=0.1")+scale_fill_grey(start = 0.3, end = 0.8)
p

####### n

i=1
j=1
k=1
delta = c(1:5)/10
p=c(60)
n=c(30,40,50,60)
m.prop=matrix(nrow=200,ncol=4)
m.sing=matrix(nrow=200,ncol=4)
for(k in 1:length(n)){
  for(rounds in 1:200){
    theta = log(1+c(1:p[j])*delta[i])
    b = runif(n[k],1,3)
    Theta = matrix(c(rep(theta,n[k]/2),rep(0,(n[k]/2)*p[j])),nrow=n[k],byrow=T)+ 
      b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,k]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.sing[rounds,k]=max(kendall.tau(1:p[j],order(order(svd(data)$v[,1]))),
                        kendall.tau(p[j]:1,order(order(svd(data)$v[,1]))))
    
  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.sing)),
                   n = rep(rep(as.factor(n),each=200),2),
                   Method = rep(c("pi.hat","pi.svd"),each=200*length(n)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=n, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S1: p=60, alpha=0.1")+scale_fill_grey(start = 0.3, end = 0.8)
p


###S2, gamma
delta = c(2:5)/20
p=c(60)
n=c(40,50,60)
j=1
k=1
i=1
m.prop=matrix(nrow=200,ncol=length(delta))
m.sing=matrix(nrow=200,ncol=length(delta))
for(i in 1:length(delta)){
  #  for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    theta = c(1:p[j])*delta[i]
    b = runif(n[k],1,3)
    
    Theta = matrix(c(theta,rep(0,(n[k]-1)*p[j])),nrow=n[k],byrow=T)+
      b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,i]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.sing[rounds,i]=max(kendall.tau(1:p[j],order(order(svd(data)$v[,1]))),
                        kendall.tau(p[j]:1,order(order(svd(data)$v[,1]))))
  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.sing)),
                   alpha = rep(rep(as.factor(delta),each=200),2),
                   Method = rep(c("pi.hat","pi.svd"),each=200*length(delta)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=alpha, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S2: n=40, p=60")+scale_fill_grey(start = 0.3, end = 0.8)
p


######p
i=1
j=1
k=1
delta = c(2:5)/20
p=c(45,60,75,90)
n=c(40,50,60)
m.prop=matrix(nrow=200,ncol=4)
m.sing=matrix(nrow=200,ncol=4)
#for(i in 1:length(delta)){
for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    theta = c(1:p[j])*delta[i]
    b = runif(n[k],1,3)
    
    Theta = matrix(c(theta,rep(0,(n[k]-1)*p[j])),nrow=n[k],byrow=T)+
      b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,j]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.sing[rounds,j]=max(kendall.tau(1:p[j],order(order(svd(data)$v[,1]))),
                        kendall.tau(p[j]:1,order(order(svd(data)$v[,1]))))
  }
  #    }
  #  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.sing)),
                   p = rep(rep(as.factor(p),each=200),2),
                   Method = rep(c("pi.hat","pi.svd"),each=200*length(p)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=p, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S2: n=40, alpha=0.1")+scale_fill_grey(start = 0.3, end = 0.8)
p

####### n

i=1
j=1
k=1
delta = c(1:5)/10
p=c(60)
n=c(30,40,50,60)
m.prop=matrix(nrow=200,ncol=4)
m.sing=matrix(nrow=200,ncol=4)
for(k in 1:length(n)){
  for(rounds in 1:200){
    theta = c(1:p[j])*delta[i]
    b = runif(n[k],1,3)
    Theta = matrix(c(theta,rep(0,(n[k]-1)*p[j])),nrow=n[k],byrow=T)+
      b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,k]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.sing[rounds,k]=max(kendall.tau(1:p[j],order(order(svd(data)$v[,1]))),
                        kendall.tau(p[j]:1,order(order(svd(data)$v[,1]))))
    
  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.sing)),
                   n = rep(rep(as.factor(n),each=200),2),
                   Method = rep(c("pi.hat","pi.svd"), each=200*length(n)) )
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=n, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S2: p=60, alpha=0.1")+scale_fill_grey(start = 0.3, end = 0.8)
p
