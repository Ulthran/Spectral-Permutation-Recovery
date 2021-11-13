#This file contains the simulation codes under the 0-1 loss functions.

###S1--dense, gamma
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
    Theta=matrix(nrow=n[k],ncol=p)
    for(row.i in 1:n[k]/2){
      Theta[row.i,] =  log(1+c(1:p[j])*runif(1,delta[i]/2,delta[i]))
    }
    for(row.i in (n[k]/2+1):n[k]){
      Theta[row.i,] = log(1+c(1:p[j])*runif(1,0,0.01))
    }
    b = runif(n[k],1,3)
    
    Theta = Theta+ b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,0.025)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,i]=min(sum(1:p[j]!=order(order(svd(data.c)$v[,1]))),
                         sum(p[j]:1!=order(order(svd(data.c)$v[,1]))))>0
    m.mean[rounds,i]= min(sum(1:p[j] != order(order(apply(data,2,mean),decreasing=T))),
                          sum(p[j]:1 != order(order(apply(data,2,mean),decreasing=T))))>0
    m.max[rounds,i]=min(sum(1:p[j] != order(order(apply(data,2,max),decreasing=T))),
                        sum(p[j]:1 != order(order(apply(data,2,max),decreasing=T))))>0
  }
}


data.frame(colMeans(m.prop),colMeans(m.mean),colMeans(m.max))


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
      Theta[row.i,] =  log(1+c(1:p[j])*runif(1,delta[i]/2,delta[i]))
    }
    for(row.i in (n[k]/2+1):n[k]){
      Theta[row.i,] = log(1+c(1:p[j])*runif(1,0,0.01))
    }
    b = runif(n[k],1,3)
    
    Theta = Theta+ b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,0.025)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,j]=min(sum(1:p[j]!=order(order(svd(data.c)$v[,1]))),
                         sum(p[j]:1!=order(order(svd(data.c)$v[,1]))))>0
    m.mean[rounds,j]= min(sum(1:p[j] != order(order(apply(data,2,mean),decreasing=T))),
                          sum(p[j]:1 != order(order(apply(data,2,mean),decreasing=T))))>0
    m.max[rounds,j]=min(sum(1:p[j] != order(order(apply(data,2,max),decreasing=T))),
                        sum(p[j]:1 != order(order(apply(data,2,max),decreasing=T))))>0
  }
}


data.frame(colMeans(m.prop),colMeans(m.mean),colMeans(m.max))

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
    Z=rnorm(n[k]*p[j],0,0.025)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,k]=min(sum(1:p[j]!=order(order(svd(data.c)$v[,1]))),
                         sum(p[j]:1!=order(order(svd(data.c)$v[,1]))))>0
    m.mean[rounds,k]= min(sum(1:p[j] != order(order(apply(data,2,mean),decreasing=T))),
                          sum(p[j]:1 != order(order(apply(data,2,mean),decreasing=T))))>0
    m.max[rounds,k]=min(sum(1:p[j] != order(order(apply(data,2,max),decreasing=T))),
                        sum(p[j]:1 != order(order(apply(data,2,max),decreasing=T))))>0
  }
}


data.frame(colMeans(m.prop),colMeans(m.mean),colMeans(m.max))



###S1 -- sparse, gamma
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
      Theta[row.i,] =  log(1+c(1:p[j])*runif(1,delta[i]/2,delta[i]))
    }
    for(row.i in 4:n[k]){
      Theta[row.i,] = log(1+c(1:p[j])*runif(1,0,0.01))
    }
    b = runif(n[k],1,3)
    
    Theta = Theta+ b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,0.0075)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,i]=min(sum(1:p[j]!=order(order(svd(data.c)$v[,1]))),
                         sum(p[j]:1!=order(order(svd(data.c)$v[,1]))))>0
    m.mean[rounds,i]= min(sum(1:p[j] != order(order(apply(data,2,mean),decreasing=T))),
                          sum(p[j]:1 != order(order(apply(data,2,mean),decreasing=T))))>0
    m.max[rounds,i]=min(sum(1:p[j] != order(order(apply(data,2,max),decreasing=T))),
                        sum(p[j]:1 != order(order(apply(data,2,max),decreasing=T))))>0
  }
}

data.frame(colMeans(m.prop),colMeans(m.mean),colMeans(m.max))


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
    for(row.i in 1:3){
      Theta[row.i,] =  log(1+c(1:p[j])*runif(1,delta[i]/2,delta[i]))
    }
    for(row.i in 4:n[k]){
      Theta[row.i,] = log(1+c(1:p[j])*runif(1,0,0.01))
    }
    b = runif(n[k],1,3)
    
    Theta = Theta+ b %o% rep(1,p[j])
    Z=rnorm(n[k]*p[j],0,0.0075)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,j]=min(sum(1:p[j]!=order(order(svd(data.c)$v[,1]))),
                         sum(p[j]:1!=order(order(svd(data.c)$v[,1]))))>0
    m.mean[rounds,j]= min(sum(1:p[j] != order(order(apply(data,2,mean),decreasing=T))),
                          sum(p[j]:1 != order(order(apply(data,2,mean),decreasing=T))))>0
    m.max[rounds,j]=min(sum(1:p[j] != order(order(apply(data,2,max),decreasing=T))),
                        sum(p[j]:1 != order(order(apply(data,2,max),decreasing=T))))>0
  }
  #    }
  #  }
}

data.frame(colMeans(m.prop),colMeans(m.mean),colMeans(m.max))

####### n

i=1
j=1
k=1
delta = c(2:5)/20
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
    Z=rnorm(n[k]*p[j],0,0.0075)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,k]=min(sum(1:p[j]!=order(order(svd(data.c)$v[,1]))),
                         sum(p[j]:1!=order(order(svd(data.c)$v[,1]))))>0
    m.mean[rounds,k]= min(sum(1:p[j] != order(order(apply(data,2,mean),decreasing=T))),
                          sum(p[j]:1 != order(order(apply(data,2,mean),decreasing=T))))>0
    m.max[rounds,k]=min(sum(1:p[j] != order(order(apply(data,2,max),decreasing=T))),
                        sum(p[j]:1 != order(order(apply(data,2,max),decreasing=T))))>0
    
  }
}

data.frame(colMeans(m.prop),colMeans(m.mean),colMeans(m.max))

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
    Z=rnorm(n[k]*p[j],0,0.1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,i]=min(sum(1:p[j]!=order(order(svd(data.c)$v[,1]))),
                         sum(p[j]:1!=order(order(svd(data.c)$v[,1]))))>0
    m.mean[rounds,i]= min(sum(1:p[j] != order(order(apply(data,2,mean),decreasing=T))),
                          sum(p[j]:1 != order(order(apply(data,2,mean),decreasing=T))))>0
    m.max[rounds,i]=min(sum(1:p[j] != order(order(apply(data,2,max),decreasing=T))),
                        sum(p[j]:1 != order(order(apply(data,2,max),decreasing=T))))>0
  }
}

data.frame(colMeans(m.prop),colMeans(m.mean),colMeans(m.max))


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
    Z=rnorm(n[k]*p[j],0,0.1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,j]=min(sum(1:p[j]!=order(order(svd(data.c)$v[,1]))),
                         sum(p[j]:1!=order(order(svd(data.c)$v[,1]))))>0
    m.mean[rounds,j]= min(sum(1:p[j] != order(order(apply(data,2,mean),decreasing=T))),
                          sum(p[j]:1 != order(order(apply(data,2,mean),decreasing=T))))>0
    m.max[rounds,j]=min(sum(1:p[j] != order(order(apply(data,2,max),decreasing=T))),
                        sum(p[j]:1 != order(order(apply(data,2,max),decreasing=T))))>0
  }
  #    }
  #  }
}

data.frame(colMeans(m.prop),colMeans(m.mean),colMeans(m.max))

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
    Z=rnorm(n[k]*p[j],0,0.1)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,k]=min(sum(1:p[j]!=order(order(svd(data.c)$v[,1]))),
                         sum(p[j]:1!=order(order(svd(data.c)$v[,1]))))>0
    m.mean[rounds,k]= min(sum(1:p[j] != order(order(apply(data,2,mean),decreasing=T))),
                          sum(p[j]:1 != order(order(apply(data,2,mean),decreasing=T))))>0
    m.max[rounds,k]=min(sum(1:p[j] != order(order(apply(data,2,max),decreasing=T))),
                        sum(p[j]:1 != order(order(apply(data,2,max),decreasing=T))))>0
    
  }
}

data.frame(colMeans(m.prop),colMeans(m.mean),colMeans(m.max))

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
    Z=rnorm(n[k]*p[j],0,0.025)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,i]=min(sum(1:p[j]!=order(order(svd(data.c)$v[,1]))),
                         sum(p[j]:1!=order(order(svd(data.c)$v[,1]))))>0
    m.mean[rounds,i]= min(sum(1:p[j] != order(order(apply(data,2,mean),decreasing=T))),
                          sum(p[j]:1 != order(order(apply(data,2,mean),decreasing=T))))>0
    m.max[rounds,i]=min(sum(1:p[j] != order(order(apply(data,2,max),decreasing=T))),
                        sum(p[j]:1 != order(order(apply(data,2,max),decreasing=T))))>0
  }
}

data.frame(colMeans(m.prop),colMeans(m.mean),colMeans(m.max))


######p
###p
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
    Z=rnorm(n[k]*p[j],0,0.025)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,j]=min(sum(1:p[j]!=order(order(svd(data.c)$v[,1]))),
                         sum(p[j]:1!=order(order(svd(data.c)$v[,1]))))>0
    m.mean[rounds,j]= min(sum(1:p[j] != order(order(apply(data,2,mean),decreasing=T))),
                          sum(p[j]:1 != order(order(apply(data,2,mean),decreasing=T))))>0
    m.max[rounds,j]=min(sum(1:p[j] != order(order(apply(data,2,max),decreasing=T))),
                        sum(p[j]:1 != order(order(apply(data,2,max),decreasing=T))))>0
  }
  #    }
  #  }
}

data.frame(colMeans(m.prop),colMeans(m.mean),colMeans(m.max))
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
    Z=rnorm(n[k]*p[j],0,0.025)
    data=Theta+Z
    data.c=data-(rowMeans(data) %o% rep(1,p[j]))
    m.prop[rounds,k]=min(sum(1:p[j]!=order(order(svd(data.c)$v[,1]))),
                         sum(p[j]:1!=order(order(svd(data.c)$v[,1]))))>0
    m.mean[rounds,k]= min(sum(1:p[j] != order(order(apply(data,2,mean),decreasing=T))),
                          sum(p[j]:1 != order(order(apply(data,2,mean),decreasing=T))))>0
    m.max[rounds,k]=min(sum(1:p[j] != order(order(apply(data,2,max),decreasing=T))),
                        sum(p[j]:1 != order(order(apply(data,2,max),decreasing=T))))>0
    
  }
}

data.frame(colMeans(m.prop),colMeans(m.mean),colMeans(m.max))

