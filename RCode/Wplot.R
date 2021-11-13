#This file contains the simulation codes for the Figure 2 in our online supplementary material.

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
W=matrix(nrow=200,ncol=n[k])
W.s=matrix(nrow=200,ncol=n[k])
#for(i in 1:length(delta)){
  #  for(j in 1:length(p)){
  #    for(k in 1:length(n)){
  for(rounds in 1:200){
    
    Theta=matrix(nrow=n[k],ncol=p[j])
    for(row.i in 1:n[k]/2){
      Theta[row.i,] =  log(1+c(1:p[j])*runif(1,delta[i]/2,delta[i]))
    }
    for(row.i in (n[k]/2+1):n[k]){
      Theta[row.i,] = log(1+c(1:p[j])*runif(1,0,delta[i]/10))
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
  
    if(svd(data.c)$u[1,1]>0) W[rounds,]= svd(data.c)$u[,1] else
      W[rounds,]= -svd(data.c)$u[,1]
    
    temp=rep(0,n[k])
    for(bb in 1:p[j]){
      temp[ apply(data,2,which.max)[bb]]=temp[ apply(data,2,which.max)[bb]]+1
    }
    W.s[rounds,]=temp/p[j]
    
 # }
}

#n=40, p=75, alpha=0.1, S1
corrplot(t(as.matrix(W)), outline = FALSE, method="color", type = "full",is.corr=FALSE, tl.pos="n")
corrplot(t(as.matrix(W.s)), outline = FALSE, method="color", type = "full",is.corr=FALSE, tl.pos="n")

##S3
delta = c(2:5)*200
p=c(75)
n=c(40,50,60)
j=1
k=1
i=1
m.prop=matrix(nrow=200,ncol=length(delta))
m.mean=matrix(nrow=200,ncol=length(delta))
m.max=matrix(nrow=200,ncol=length(delta))
W=matrix(nrow=200,ncol=n[k])
W.s=matrix(nrow=200,ncol=n[k])
#for(i in 1:length(delta)){
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
    if(svd(data.c)$u[1,1]>0) W[rounds,]= svd(data.c)$u[,1] else
    W[rounds,]= -svd(data.c)$u[,1]
    
    temp=rep(0,n[k])
    for(bb in 1:p[j]){
      temp[ apply(data,2,which.max)[bb]]=temp[ apply(data,2,which.max)[bb]]+1
    }
    W.s[rounds,]=temp/p[j]
  }
#}

#n=40, p=75, alpha=400, S3
corrplot(t(as.matrix(W)), outline = FALSE, method="color", type = "full",is.corr=FALSE, tl.pos="n")
corrplot(t(as.matrix(W.s)), outline = FALSE, method="color", type = "full",is.corr=FALSE, tl.pos="n")

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
W=matrix(nrow=200,ncol=n[k])
W.s=matrix(nrow=200,ncol=n[k])
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
    if(svd(data.c)$u[1,1]>0) W[rounds,]= svd(data.c)$u[,1] else
      W[rounds,]= -svd(data.c)$u[,1]
    
    temp=rep(0,n[k])
    for(bb in 1:p[j]){
      temp[ apply(data,2,which.max)[bb]]=temp[ apply(data,2,which.max)[bb]]+1
    }
    W.s[rounds,]=temp/p[j]
  }
}

#n=40, p=75, alpha=.1, S2
corrplot(t(as.matrix(W)), outline = FALSE, method="color", type = "full",is.corr=FALSE, tl.pos="n")
corrplot(t(as.matrix(W.s)), outline = FALSE, method="color", type = "full",is.corr=FALSE, tl.pos="n")


###S4 
delta = c(2:5)/20
p=c(75)
n=c(40,50,60)
j=1
k=1
i=1
m.prop=matrix(nrow=200,ncol=length(delta))
m.mean=matrix(nrow=200,ncol=length(delta))
m.max=matrix(nrow=200,ncol=length(delta))
W=matrix(nrow=200,ncol=n[k])
W.s=matrix(nrow=200,ncol=n[k])
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
    if(svd(data.c)$u[1,1]>0) W[rounds,]= svd(data.c)$u[,1] else
      W[rounds,]= -svd(data.c)$u[,1]
    
    temp=rep(0,n[k])
    for(bb in 1:p[j]){
      temp[ apply(data,2,which.max)[bb]]=temp[ apply(data,2,which.max)[bb]]+1
    }
    W.s[rounds,]=temp/p[j]
  }
}

#n=40, p=75, alpha=.1, S4
corrplot(t(as.matrix(W)), outline = FALSE, method="color", type = "full",is.corr=FALSE, tl.pos="n")
corrplot(t(as.matrix(W.s)), outline = FALSE, method="color", type = "full",is.corr=FALSE, tl.pos="n")
