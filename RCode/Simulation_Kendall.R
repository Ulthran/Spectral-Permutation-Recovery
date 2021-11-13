#This file contains the simulation codes under the normalized Kendall's tau distance measures.

library(VGAM)
library(ggplot2)
library(corrplot)

###S1, gamma
delta = c(2:5)/20
p=c(75)
n=c(40,50,60)
j=1
k=1
i=1
m.prop=matrix(nrow=200,ncol=length(delta))
m.mean=matrix(nrow=200,ncol=length(delta))
m.max=matrix(nrow=200,ncol=length(delta))
for(i in 1:length(delta)){
  #  for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:n[k]/2){
      Theta[row.i,] =  log(1+c(1:p[j])*runif(1,delta[i]/2,delta[i]))
    }
    for(row.i in (n[k]/2+1):n[k]){
      Theta[row.i,] = log(1+c(1:p[j])*runif(1,0,0.01))
    }
    b = runif(n[k],1,3)
    
    Theta = Theta+ b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,i]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.mean[rounds,i]= max(kendall.tau(1:p[j],order(order(apply(data,2,mean),decreasing=T))),
                          kendall.tau(p[j]:1,order(order(apply(data,2,mean),decreasing=T))))
    m.max[rounds,i]=max(kendall.tau(1:p[j],order(order(apply(data,2,max),decreasing=T))),
                        kendall.tau(p[j]:1,order(order(apply(data,2,max),decreasing=T))))
  
  
    }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.mean),as.vector(m.max)),
                   alpha= rep(rep(as.factor(delta),each=200),3),
                   Method = rep(c("pi.hat","pi.mean","pi.max"),each=200*length(delta)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=alpha, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S1: n=40, p=75")+scale_fill_grey(start = 0.3, end = 0.8)
p

######p

i=1
j=1
k=1
delta = c(2:5)/20
p=c(60,75,90,120)
n=c(40,50,60)
m.prop=matrix(nrow=200,ncol=4)
m.mean=matrix(nrow=200,ncol=4)
m.max=matrix(nrow=200,ncol=4)
#for(i in 1:length(delta)){
for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:n[k]/2){
      Theta[row.i,] =  log(1+c(1:p[j])*runif(1,delta[i]/2,delta[i]))
    }
    for(row.i in (n[k]/2+1):n[k]){
      Theta[row.i,] = log(1+c(1:p[j])*runif(1,0,0.01))
    }
    b = runif(n[k],1,3)
    
    Theta = Theta+ b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,j]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.mean[rounds,j]= max(kendall.tau(1:p[j],order(order(apply(data,2,mean),decreasing=T))),
                          kendall.tau(p[j]:1,order(order(apply(data,2,mean),decreasing=T))))
    m.max[rounds,j]=max(kendall.tau(1:p[j],order(order(apply(data,2,max),decreasing=T))),
                        kendall.tau(p[j]:1,order(order(apply(data,2,max),decreasing=T))))
  }
  #    }
  #  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.mean),as.vector(m.max)),
                   p = rep(rep(as.factor(p),each=200),3),
                   Method = rep(c("pi.hat","pi.mean","pi.max"),each=200*length(p)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=p, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S1: n=40, alpha=0.1")+scale_fill_grey(start = 0.3, end = 0.8)
p

####### n

i=1
j=1
k=1
delta = c(1:5)/10
p=c(75)
n=c(30,40,50,60)
m.prop=matrix(nrow=200,ncol=4)
m.mean=matrix(nrow=200,ncol=4)
m.max=matrix(nrow=200,ncol=4)
for(k in 1:length(n)){
  for(rounds in 1:200){
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:n[k]/2){
      Theta[row.i,] =  log(1+c(1:p[j])*runif(1,delta[i]/2,delta[i]))
    }
    for(row.i in (n[k]/2+1):n[k]){
      Theta[row.i,] = log(1+c(1:p[j])*runif(1,0,0.01))
    }
    b = runif(n[k],1,3)
    
    Theta = Theta+ b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,k]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.mean[rounds,k]= max(kendall.tau(1:p[j],order(order(apply(data,2,mean),decreasing=T))),
                          kendall.tau(p[j]:1,order(order(apply(data,2,mean),decreasing=T))))
    m.max[rounds,k]=max(kendall.tau(1:p[j],order(order(apply(data,2,max),decreasing=T))),
                        kendall.tau(p[j]:1,order(order(apply(data,2,max),decreasing=T))))
    
  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.mean),as.vector(m.max)),
                   n = rep(rep(as.factor(n),each=200),3),
                   Method = rep(c("pi.hat","pi.mean","pi.max"),each=200*length(n)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=n, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S1: p=75, alpha=0.1")+scale_fill_grey(start = 0.3, end = 0.8)
p


###S1 -- sparse, gamma
delta = c(2:5)*200
p=c(75)
n=c(40,50,60)
j=1
k=1
i=1
m.prop=matrix(nrow=200,ncol=length(delta))
m.mean=matrix(nrow=200,ncol=length(delta))
m.max=matrix(nrow=200,ncol=length(delta))
for(i in 1:length(delta)){
  #  for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:3){
      Theta[row.i,] =  log(1+c(1:p[j])*runif(1,delta[i]/2,delta[i]))
    }
    for(row.i in 4:n[k]){
      Theta[row.i,] = log(1+c(1:p[j])*runif(1,0,0.01))
    }
    b = runif(n[k],1,3)
    
    Theta = Theta+ b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,i]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.mean[rounds,i]= max(kendall.tau(1:p[j],order(order(apply(data,2,mean),decreasing=T))),
                          kendall.tau(p[j]:1,order(order(apply(data,2,mean),decreasing=T))))
    m.max[rounds,i]=max(kendall.tau(1:p[j],order(order(apply(data,2,max),decreasing=T))),
                        kendall.tau(p[j]:1,order(order(apply(data,2,max),decreasing=T))))
  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.mean),as.vector(m.max)),
                   alpha= rep(rep(as.factor(delta),each=200),3),
                   Method = rep(c("pi.hat","pi.mean","pi.max"),each=200*length(delta)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=alpha, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S3: n=40, p=75")+scale_fill_grey(start = 0.3, end = 0.8)
p


###p
i=3
j=1
k=1
delta = c(2:5)*200
p=c(60,75,90,120)
n=c(40,50,60)
m.prop=matrix(nrow=200,ncol=4)
m.mean=matrix(nrow=200,ncol=4)
m.max=matrix(nrow=200,ncol=4)
#for(i in 1:length(delta)){
for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:3){
      Theta[row.i,] =  log(1+c(1:p[j])*runif(1,delta[i]/2,delta[i]))
    }
    for(row.i in 4:n[k]){
      Theta[row.i,] = log(1+c(1:p[j])*runif(1,0,0.01))
    }
    b = runif(n[k],1,3)
    
    Theta = Theta+ b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,j]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.mean[rounds,j]= max(kendall.tau(1:p[j],order(order(apply(data,2,mean),decreasing=T))),
                          kendall.tau(p[j]:1,order(order(apply(data,2,mean),decreasing=T))))
    m.max[rounds,j]=max(kendall.tau(1:p[j],order(order(apply(data,2,max),decreasing=T))),
                        kendall.tau(p[j]:1,order(order(apply(data,2,max),decreasing=T))))
  }
  #    }
  #  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.mean),as.vector(m.max)),
                   p = rep(rep(as.factor(p),each=200),3),
                   Method = rep(c("pi.hat","pi.mean","pi.max"),each=200*length(p)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=p, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S3: n=40, alpha=800")+scale_fill_grey(start = 0.3, end = 0.8)
p

####### n

i=3
j=1
k=1
delta = c(2:5)*200
p=c(75)
n=c(30,40,50,60)
m.prop=matrix(nrow=200,ncol=4)
m.mean=matrix(nrow=200,ncol=4)
m.max=matrix(nrow=200,ncol=4)
for(k in 1:length(n)){
  for(rounds in 1:200){
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:3){
      Theta[row.i,] =  log(1+c(1:p[j])*runif(1,delta[i]/2,delta[i]))
    }
    for(row.i in 4:n[k]){
      Theta[row.i,] = log(1+c(1:p[j])*runif(1,0,0.01))
    }
    b = runif(n[k],1,3)
    
    Theta = Theta+ b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,k]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.mean[rounds,k]= max(kendall.tau(1:p[j],order(order(apply(data,2,mean),decreasing=T))),
                          kendall.tau(p[j]:1,order(order(apply(data,2,mean),decreasing=T))))
    m.max[rounds,k]=max(kendall.tau(1:p[j],order(order(apply(data,2,max),decreasing=T))),
                        kendall.tau(p[j]:1,order(order(apply(data,2,max),decreasing=T))))
    
  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.mean),as.vector(m.max)),
                   n = rep(rep(as.factor(n),each=200),3),
                   Method = rep(c("pi.hat","pi.mean","pi.max"),each=200*length(n)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=n, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S3: p=75, alpha=800")+scale_fill_grey(start = 0.3, end = 0.8)
p

###S2 -- dense, gamma
delta = c(2:5)/20
p=c(75)
n=c(40,50,60)
j=1
k=1
i=1
m.prop=matrix(nrow=200,ncol=length(delta))
m.mean=matrix(nrow=200,ncol=length(delta))
m.max=matrix(nrow=200,ncol=length(delta))
for(i in 1:length(delta)){
  #  for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:n[k]/2){
      Theta[row.i,] =  c(1:p[j])*runif(1,delta[i]/2,delta[i])
    }
    for(row.i in (n[k]/2+1):n[k]){
      Theta[row.i,] = c(1:p[j])*runif(1,0,delta[i]/10)
    }
    
    b = runif(n[k],1,3)
    
    Theta = Theta+b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,i]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.mean[rounds,i]= max(kendall.tau(1:p[j],order(order(apply(data,2,mean),decreasing=T))),
                          kendall.tau(p[j]:1,order(order(apply(data,2,mean),decreasing=T))))
    m.max[rounds,i]=max(kendall.tau(1:p[j],order(order(apply(data,2,max),decreasing=T))),
                        kendall.tau(p[j]:1,order(order(apply(data,2,max),decreasing=T))))
  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.mean),as.vector(m.max)),
                   alpha = rep(rep(as.factor(delta),each=200),3),
                   Method = rep(c("pi.hat","pi.mean","pi.max"),each=200*length(delta)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=alpha, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S2: n=40, p=75")+scale_fill_grey(start = 0.3, end = 0.8)
p


######p

i=1
j=1
k=1
delta = c(2:5)/20
p=c(45,60,75,90)
n=c(40,50,60)
m.prop=matrix(nrow=200,ncol=4)
m.mean=matrix(nrow=200,ncol=4)
m.max=matrix(nrow=200,ncol=4)
#for(i in 1:length(delta)){
for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:n[k]/2){
      Theta[row.i,] =  c(1:p[j])*runif(1,delta[i]/2,delta[i])
    }
    for(row.i in (n[k]/2+1):n[k]){
      Theta[row.i,] = c(1:p[j])*runif(1,0,delta[i]/10)
    }
    
    b = runif(n[k],1,3)
    
    Theta = Theta+b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,j]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.mean[rounds,j]= max(kendall.tau(1:p[j],order(order(apply(data,2,mean),decreasing=T))),
                          kendall.tau(p[j]:1,order(order(apply(data,2,mean),decreasing=T))))
    m.max[rounds,j]=max(kendall.tau(1:p[j],order(order(apply(data,2,max),decreasing=T))),
                        kendall.tau(p[j]:1,order(order(apply(data,2,max),decreasing=T))))
  }
  #    }
  #  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.mean),as.vector(m.max)),
                   p = rep(rep(as.factor(p),each=200),3),
                   Method = rep(c("pi.hat","pi.mean","pi.max"),each=200*length(p)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=p, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S2: n=40, alpha=0.1")+scale_fill_grey(start = 0.3, end = 0.8)
p

####### n

i=1
j=1
k=1
delta = c(1:5)/10
p=c(75)
n=c(30,40,50,60)
m.prop=matrix(nrow=200,ncol=4)
m.mean=matrix(nrow=200,ncol=4)
m.max=matrix(nrow=200,ncol=4)
for(k in 1:length(n)){
  for(rounds in 1:200){
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:n[k]/2){
      Theta[row.i,] =  c(1:p[j])*runif(1,delta[i]/2,delta[i])
    }
    for(row.i in (n[k]/2+1):n[k]){
      Theta[row.i,] = c(1:p[j])*runif(1,0,delta[i]/10)
    }
    
    b = runif(n[k],1,3)
    
    Theta = Theta+b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,k]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.mean[rounds,k]= max(kendall.tau(1:p[j],order(order(apply(data,2,mean),decreasing=T))),
                          kendall.tau(p[j]:1,order(order(apply(data,2,mean),decreasing=T))))
    m.max[rounds,k]=max(kendall.tau(1:p[j],order(order(apply(data,2,max),decreasing=T))),
                        kendall.tau(p[j]:1,order(order(apply(data,2,max),decreasing=T))))
    
  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.mean),as.vector(m.max)),
                   n = rep(rep(as.factor(n),each=200),3),
                   Method = rep(c("pi.hat","pi.mean","pi.max"),each=200*length(n)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=n, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S2: p=75, alpha=0.1")+scale_fill_grey(start = 0.3, end = 0.8)
p


###S2 -- sparse, gamma
delta = c(2:5)/20
p=c(75)
n=c(40,50,60)
j=1
k=1
i=1
m.prop=matrix(nrow=200,ncol=length(delta))
m.mean=matrix(nrow=200,ncol=length(delta))
m.max=matrix(nrow=200,ncol=length(delta))
for(i in 1:length(delta)){
  #  for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:3){
      Theta[row.i,] =  c(1:p[j])*runif(1,delta[i]/2,delta[i])
    }
    for(row.i in 4:n[k]){
      Theta[row.i,] = c(1:p[j])*runif(1,0,delta[i]/10)
    }
    
    b = runif(n[k],1,3)
    
    Theta = Theta+b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,i]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.mean[rounds,i]= max(kendall.tau(1:p[j],order(order(apply(data,2,mean),decreasing=T))),
                          kendall.tau(p[j]:1,order(order(apply(data,2,mean),decreasing=T))))
    m.max[rounds,i]=max(kendall.tau(1:p[j],order(order(apply(data,2,max),decreasing=T))),
                        kendall.tau(p[j]:1,order(order(apply(data,2,max),decreasing=T))))
  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.mean),as.vector(m.max)),
                   alpha = rep(rep(as.factor(delta),each=200),3),
                   Method = rep(c("pi.hat","pi.mean","pi.max"),each=200*length(delta)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=alpha, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S4: n=40, p=75")+scale_fill_grey(start = 0.3, end = 0.8)
p


######p

i=1
j=1
k=1
delta = c(2:5)/20
p=c(45,60,75,90)
n=c(40,50,60)
m.prop=matrix(nrow=200,ncol=4)
m.mean=matrix(nrow=200,ncol=4)
m.max=matrix(nrow=200,ncol=4)
#for(i in 1:length(delta)){
for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:3){
      Theta[row.i,] =  c(1:p[j])*runif(1,delta[i]/2,delta[i])
    }
    for(row.i in 4:n[k]){
      Theta[row.i,] = c(1:p[j])*runif(1,0,delta[i]/10)
    }
    
    b = runif(n[k],1,3)
    
    Theta = Theta+b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,j]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.mean[rounds,j]= max(kendall.tau(1:p[j],order(order(apply(data,2,mean),decreasing=T))),
                          kendall.tau(p[j]:1,order(order(apply(data,2,mean),decreasing=T))))
    m.max[rounds,j]=max(kendall.tau(1:p[j],order(order(apply(data,2,max),decreasing=T))),
                        kendall.tau(p[j]:1,order(order(apply(data,2,max),decreasing=T))))
  }
  #    }
  #  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.mean),as.vector(m.max)),
                   p = rep(rep(as.factor(p),each=200),3),
                   Method = rep(c("pi.hat","pi.mean","pi.max"),each=200*length(p)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=p, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S4: n=40, alpha=0.1")+scale_fill_grey(start = 0.3, end = 0.8)
p

####### n

i=1
j=1
k=1
delta = c(1:5)/10
p=c(75)
n=c(30,40,50,60)
m.prop=matrix(nrow=200,ncol=4)
m.mean=matrix(nrow=200,ncol=4)
m.max=matrix(nrow=200,ncol=4)
for(k in 1:length(n)){
  for(rounds in 1:200){
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:3){
      Theta[row.i,] =  c(1:p[j])*runif(1,delta[i]/2,delta[i])
    }
    for(row.i in 4:n[k]){
      Theta[row.i,] = c(1:p[j])*runif(1,0,delta[i]/10)
    }
    
    b = runif(n[k],1,3)
    
    Theta = Theta+b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,k]=max(kendall.tau(1:p[j],order(order(svd(data.c)$v[,1]))),
                         kendall.tau(p[j]:1,order(order(svd(data.c)$v[,1]))))
    m.mean[rounds,k]= max(kendall.tau(1:p[j],order(order(apply(data,2,mean),decreasing=T))),
                          kendall.tau(p[j]:1,order(order(apply(data,2,mean),decreasing=T))))
    m.max[rounds,k]=max(kendall.tau(1:p[j],order(order(apply(data,2,max),decreasing=T))),
                        kendall.tau(p[j]:1,order(order(apply(data,2,max),decreasing=T))))
    
  }
}

results=data.frame(Kendall.tau= c(as.vector(m.prop),as.vector(m.mean),as.vector(m.max)),
                   n = rep(rep(as.factor(n),each=200),3),
                   Method = rep(c("pi.hat","pi.mean","pi.max"),each=200*length(n)))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=n, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("S4: p=75, alpha=0.1")+scale_fill_grey(start = 0.3, end = 0.8)
p

