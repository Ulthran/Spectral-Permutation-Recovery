#This file contains codes for the analysis of the synthetic microbiome datasets.

library(VGAM)

m.prop=c()
m.prop2=c()
m.prop3=c()
m.prop4=c()
m.prop5=c()
m.naive.sum=c()
m.naive.sum2=c()
m.naive.sum3=c()
m.naive.sum4=c()
m.naive.sum5=c()

n.size=35 # set to be one of 15,20,25,30 and 35, then change the variables 'm.prop' and 'm.naive.sum' to one of 'm.prop#' and 'm.naive.sum#' with #=2,3,4,5.
for(i in 1:41){
  filename=paste0("/Data/Subset12info/maxbin.",sprintf("%03d", i),"_ptr.txt.subset1.info") #set as the path of the directory 'Subset12info'
  data=fread(file=filename,header= T)
  data=t(data)
  colnames(data)=data[1,]
  data=data[-1,]
  p=dim(data)[2]
  n=dim(data)[1]
  if(sum(is.na(as.numeric(data[n,])))>0){
    data=data[,-which(is.na(as.numeric(data[n,])))]
    p=dim(data)[2]
    n=dim(data)[1]}
  t.dist=as.numeric(data[n,])
  data=data[1:(n-6),]
  data=as.matrix(data)
  p=dim(data)[2]
  n=dim(data)[1]
  data=matrix(as.numeric(data),ncol=p)
  if(sum(is.na(colSums(data)))>0){
    t.dist=t.dist[-which(is.na(colSums(data)))]
    data=data[,-which(is.na(colSums(data)))]
    p=dim(data)[2]
    n=dim(data)[1]
  }
  #data = rbind(data,rnorm(p))
  data.c=data-(rowMeans(data) %o% rep(1,p))
  sel = sample(n,n.size,replace=T)
  data.c=data.c[sel,]
  data=data[sel,]
  
  m.naive = c()
  for(r in 1:length(sel)){
    m.naive[r] = max(kendall.tau(order(order(t.dist,decreasing=F)),order(order(data.c[r,]))),
                     kendall.tau(order(order(t.dist,decreasing=T)),order(order(data.c[r,]))))
  }
  
  m.prop[i]=max(kendall.tau(order(order(t.dist,decreasing=F)),order(order(svd(data.c)$v[,1]))),
                kendall.tau(order(order(t.dist,decreasing=T)),order(order(svd(data.c)$v[,1]))))
  m.naive.sum[i] = mean(m.naive)
}

results=data.frame(Kendall.tau= c(m.prop,m.naive.sum,m.prop2,m.naive.sum2,
                                  m.prop3,m.naive.sum3,m.prop4,m.naive.sum4,
                                  m.prop5,m.naive.sum5),
                   n = rep(rep(as.factor(c(15,20,25,30,35)),each=2*41)),
                   Method = rep(rep(c("pi.hat","iRep"),each=41),5))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=n, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+ggtitle("Synthetic Data")+scale_fill_grey(start = 0.3, end = 0.8)
p

############# change p


library(VGAM)

m.prop=c()
m.prop2=c()
m.prop3=c()
m.prop4=c()
m.prop5=c()
m.naive.sum=c()
m.naive.sum2=c()
m.naive.sum3=c()
m.naive.sum4=c()
m.naive.sum5=c()


n.size=35
p.size = 70  # set to be one of 70,80,90,100 and 110, then change the variables 'm.prop' and 'm.naive.sum' to one of 'm.prop#' and 'm.naive.sum#' with #=2,3,4,5.

for(i in 1:41){
  filename=paste0("/Data/Subset12info/maxbin.",sprintf("%03d", i),"_ptr.txt.subset2.info") #set as the path of the directory 'Subset12info'
  data=fread(file=filename,header= T)
  data=t(data)
  colnames(data)=data[1,]
  data=data[-1,]
  p=dim(data)[2]
  n=dim(data)[1]
  if(sum(is.na(as.numeric(data[n,])))>0){
    data=data[,-which(is.na(as.numeric(data[n,])))]
    p=dim(data)[2]
    n=dim(data)[1]}
  t.dist=as.numeric(data[n,])
  data=data[1:(n-6),]
  data=as.matrix(data)
  p=dim(data)[2]
  n=dim(data)[1]
  data=matrix(as.numeric(data),ncol=p)
  if(sum(is.na(colSums(data)))>0){
    t.dist=t.dist[-which(is.na(colSums(data)))]
    data=data[,-which(is.na(colSums(data)))]
    p=dim(data)[2]
    n=dim(data)[1]
  }
  # data = rbind(data,rnorm(p,0,0.5))
  data=rbind(data,0.1*order(order(t.dist,decreasing=T)))
  data.c=data-(rowMeans(data) %o% rep(1,p))
  p.sel=sample(p,p.size,replace=T)
  data.c=data.c[,p.sel]
  t.dist=t.dist[p.sel]
  data=data[,p.sel]
  m.naive = c()
  for(r in 1:n){
    m.naive[r] = max(kendall.tau(order(order(t.dist,decreasing=F)),order(order(data[r,]))),
                     kendall.tau(order(order(t.dist,decreasing=T)),order(order(data[r,]))))
  }
  
  m.prop[i]=max(kendall.tau(order(order(t.dist,decreasing=F)),order(order(svd(data.c)$v[,1]))),
                kendall.tau(order(order(t.dist,decreasing=T)),order(order(svd(data.c)$v[,1]))))
  m.naive.sum[i] = mean(m.naive)
}


results=data.frame(Kendall.tau= c(m.prop,m.naive.sum,m.prop2,m.naive.sum2,
                                  m.prop3,m.naive.sum3,m.prop4,m.naive.sum4,
                                  m.prop5,m.naive.sum5),
                   p = rep(rep(as.factor(c(70,80,90,100,110)),each=2*41)),
                   Method = rep(rep(c("pi.hat","iRep"),each=41),5))
results$Kendall.tau=(1-results$Kendall.tau)/2

p<-ggplot(results, aes(x=p, y=Kendall.tau, fill=Method)) +
  geom_boxplot()+scale_fill_grey(start = 0.3, end = 0.8)
p

