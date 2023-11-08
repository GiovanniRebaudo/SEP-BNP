rm(list=ls())

library(DirichletReg)
library(invgamma)
library(tidyr)
library(dplyr)
library(truncnorm)

#read-in and clean data
Dietswap_dataset <- readRDS("Dietswap_dataset.RDS")

data = Dietswap_dataset[Dietswap_dataset$timepoint==1,]
data2=pivot_wider(data[,1:3],names_from = Sample,values_from = Abundance)

na_in_rows=sapply(1:dim(data2)[1],function(i) sum(data2[i,]==0))
#remove all na rows
data3=data.frame(data2[na_in_rows!=38,])

y=data3[,2:39]
rownames(y)=data3[,1]
#y=log(y+1)
N=dim(y)[1]*dim(y)[2]
n=dim(y)[2]
B=dim(y)[1]
# normalize y 
gammas=colMeans(y)
y=sapply(1:n,function(i) as.matrix(y)[,i]/gammas[i])


Update_mu <- function(zl,nl,sig2_l,mu0,sig0){
  var=1/(1/sig0^2 + nl/sig2_l)
  mean=var*(mu0/sig0^2+sum(zl)/sig2_l)
  mu=rtruncnorm(1,0,Inf,mean,sqrt(var))
  return(mu)
}


Update_sig2 <- function(zl,nl,mu_l,a,b){
  sig2=rinvgamma(1,a+nl/2,b+sum((zl-mu_l)^2)/2)
  return(sig2)
}

###################

#cluster number truncated at:
K=3
L=8

# priors and hyperparameters 
alpha=10
beta =10

vk=c(rbeta(K-1,1,alpha),1)
pi=c(vk[1],sapply(2:K,function(i) vk[i]*prod(1-vk[1:(i-1)])))

vkl=matrix(0,K,L)
for(k in 1:K){
  vkl[k,]=c(rbeta(L-1,1,beta),1)
}

w=matrix(0,K,L)
for(k in 1:K){
  w[k,]=c(vkl[k,1],sapply(2:L,function(i) vkl[k,i]*prod(1-vkl[k,1:(i-1)])))
}

#mu_l and sig2_l 
mu0=mean(as.matrix(y))
sig0=10
#mu=rnorm(L,mu0,sig0)

colvar <- function(y){
  return(sapply(1:dim(y)[2],function(j) var(y[,j])))
}

z=as.matrix(y)
rowvar <- function(y){
  return(sapply(1:dim(y)[1],function(j) var(y[j,])))
}

a0=3
b0=2*var(as.vector(as.matrix(y)))
#other test prior settings 
#b0=120
#b0=1000
#b0=2*mean(colvar(as.matrix((y))))
#b0=2*var(as.vector(as.matrix(y)))
#sig2=rinvgamma(L,a0,b0)

#initialize differently
#mu=sapply(1:L, function(l) mean(y[mik_expand==l]))
#sig2=sapply(1:L, function(l) var(y[mik_expand==l]))

###################################################
#use hclust for patients   Sj and pi 
# order by most frequent OTUs 
# hclust for OTU  -- mu_l and sig_l

# in each cluster fix one OTU -- then can use non-informative prior
# sig0 much bigger than now to cover the range of OTUs
# truncate the crazy normal by maximum and minimum of OTUS
### -- when we sample from prior for empty clusters

#initialize Sj with hclust
ty=data.frame(t(y))
ds=dist(ty,method='euclidean')
hc <- hclust(ds, method = "complete" )
sub_grp <- cutree(hc, k = K)
table(sub_grp)
ty %>% mutate(cluster = sub_grp) -> ty
Sj=ty$cluster
Sj0=Sj

### initialize mu with hclust
ydf <- data.frame(y[-c(1),])
d <- dist(ydf, method = "euclidean")
hc1 <- hclust(d, method = "complete" )
sub_grp <- cutree(hc1, k = L-1)
table(sub_grp)
ydf %>% mutate(cluster = sub_grp) -> ydf
mu=rep(0,L)
for(l in 2:L){
  mu[l]= mean(as.matrix(ydf[ydf$cluster == l-1, ]%>% select(-cluster)))
}
mu[1]=mean(as.matrix(y[c(1),]))
sig2=rep(var(as.vector(as.matrix(ydf[ydf$cluster == which.max(table(ydf$cluster)), ]%>% select(-cluster)))),L)

as.data.frame(y) %>% mutate(all= rowSums(y)) -> ydf2
ydf2 %>% arrange(desc(all)) %>% head()

# the first OTU is much bigger than others, 
# which corresponding to mu[1] here in hclust
# fix mik for i=1, mik=1

##############
#initialize mik
nk=sapply(1:K,function(k) sum(Sj==k))
mik=matrix(0,K,B)
for(k in 1:K){
  yk=as.matrix(y[,Sj==k])
  if(nk[k]==0){
    mik[k,]= sample(1:L,B,replace=T,prob=w[k,])
    mik[k,c(1)]=1
  }else{
    mik[k,c(1)]=1
    for(b in 2:B){
      log_prob=rep(0,L)
      for(l in 1:L){
        if(nk[k]==1){
          loglk=dnorm(yk[b],mu[l],sqrt(sig2[l]),log=T)
        }else{
          loglk=dnorm(yk[b,],mu[l],sqrt(sig2[l]),log=T)
        }
        sum_loglk=sum(loglk)
        if(sum_loglk==-Inf){
          sum_loglk=-10^5
        }
        log_prob[l]=log(w[k,l])+sum_loglk
      }
      log_prob=log_prob-max(log_prob)
      mik[k,b]=sample(1:L,1,prob = exp(log_prob))
      #https://en.wikipedia.org/wiki/Categorical_distribution#Sampling_via_the_Gumbel_distribution
      #gumble=-log(-log(runif(L,0,1)))
      #mik[k,b]=which.max(log_prob+gumble)
    }}
}

mik_expand=t(mik[Sj,])


####### create place holders
T=10000
pi_store=matrix(0,T+1,K)
S_store=matrix(0,T+1,n)
w_store=array(0,dim=c(T+1,K,L))
mik_store=array(0,dim=c(T+1,K,B))
mu_store=matrix(0,T+1,L)
sig2_store=matrix(0,T+1,L)
nl_store=matrix(0,T+1,L)
ll_store=rep(0,T+1)

z=matrix(0,B,n)

z=as.matrix(y)

#loglikelihood check functions 
ll_fun1 <- function(i,j){
  return(dnorm(z[i,j],mu[mik[Sj[j],i]],sqrt(sig2[mik[Sj[j],i]]),log = T))
}

ll_fun <- function(){
  return(sum(sapply(1:B, function(x) mapply(ll_fun1,x,1:n))))}

ll_last=ll_fun()

pi_store[1,]=pi
w_store[1,,]=w
S_store[1,]=Sj
mik_store[1,,]=mik
mu_store[1,]=mu
sig2_store[1,]=sig2
nl_store[1,]=sapply(1:L,function(l) sum(mik_expand==l))
ll_store[1]=ll_last

#### MCMC
for(iter in 1:T){
  
  # update vk and pi
  vk=c(sapply(1:(K-1),function(i) rbeta(1,1+nk[i],alpha+sum(nk[(i+1):K]))),1)
  pi=c(vk[1],sapply(2:K,function(i) vk[i]*prod(1-vk[1:(i-1)])))
  pi_store[iter+1,]=pi
  
  
  #update w_k 
  for(k in 1:K){
    nkl=sapply(1:L,function(l) sum(mik[k,]==l))
    vkl[k,]=c(sapply(1:(L-1),function(i) rbeta(1,1+nkl[i],alpha+sum(nkl[(i+1):L]))),1)  
    w[k,]=c(vkl[k,1],sapply(2:L,function(i) vkl[k,i]*prod(1-vkl[k,1:(i-1)])))
  }
  w_store[iter+1,,]=w
  
  ll=ll_fun()
  if(ll<ll_last-1000){
    print("after w")
    browser()
  }
  ll_last=ll
  
  #update Sj
  log_prob_k=rep(0,K)
  for(j in 1:n){
    for(k in 1:K){
      log_lk=rep(0,B)
      for(b in 1:B){
        l=mik[k,b]
        log_lk[b]= dnorm(z[b,j],mu[l],sqrt(sig2[l]),log=T)
      }
      sum_loglk=sum(log_lk)
      if(sum_loglk==-Inf){
        sum_loglk=-10^5
      }
      log_prob_k[k]=log(pi[k])+sum_loglk
    }
    log_prob_k = log_prob_k - max(log_prob_k)
    Sj[j]=sample(1:K,1,prob=exp(log_prob_k))
  }
  S_store[iter+1,]=Sj
  nk=sapply(1:K,function(k) sum(Sj==k))
  
  ll=ll_fun()
  if(ll<ll_last-1000){
    print("after S")
    browser()
  }
  ll_last=ll
  
  
  #update M_ik 
  for(k in 1:K){
    zk=z[,Sj==k]
    if(nk[k]==0){
      mik[k,]= sample(1:L,B,replace=T,prob=w[k,])
      mik[k,c(1)]=1
    }else{
      mik[k,c(1)]=1
      for(b in 2:B){
        log_prob=rep(0,L)
        for(l in 1:L){
          if(nk[k]==1){
            loglk=dnorm(zk[b],mu[l],sqrt(sig2[l]),log=T)
          }else{
            loglk=dnorm(zk[b,],mu[l],sqrt(sig2[l]),log=T)
          }
          sum_loglk=sum(loglk)
          if(sum_loglk==-Inf){
            sum_loglk=-10^5
          }
          log_prob[l]=log(w[k,l])+sum_loglk
        }
        log_prob=log_prob-max(log_prob)
        mik[k,b]=sample(1:L,1,prob = exp(log_prob))
      }}
  }
  
  
  mik_store[iter+1,,]=mik
  
  mik_expand=t(mik[Sj,])
  
  ll=ll_fun()
  if(ll<ll_last-1000){
    print("after mik")
    browser()
  }
  ll_last=ll
  
  #update mul and sigl
  for(l in 1:L){
    zl= z[mik_expand==l]
    nl=sum(mik_expand==l)
    nl_store[iter+1,l]=nl
    if(nl>0){
      mu[l]=Update_mu(zl,nl,sig2[l],mu0,sig0)
      sig2[l]=Update_sig2(zl,nl,mu[l],a0,b0)
    }else{
      mu[l]=rtruncnorm(1,0,quantile(as.matrix(y),0.9),mu0,sig0)
      sig2[l]=rinvgamma(1,a0,b0) #sig2_store[iter,l]#
    }}
  
  mu_store[iter+1,]=mu
  sig2_store[iter+1,]=sig2
  
  ll_store[iter+1]=sum(sapply(1:B, function(x) mapply(ll_fun1,x,1:n)))
  
  ll=ll_fun()
  if(ll<ll_last-1000){
    print("after mu/sig")
    browser()
  }
  ll_last=ll
  
  if(iter>1){
    if(ll_store[iter]<ll_store[iter-1]-1000){
      browser()
    }}
  
  
  print(iter)
}


##### check plotting loglikelihood, cluster y mean and weights
par(mar=c(2,2,1,1))
plot(ll_store,type='l')

par(mar=c(2,2,1,1))
par(mfrow=c(3,1))
hist(mik_store[T,1,])
hist(mik_store[T,2,])
hist(mik_store[T,3,])
par(mfrow=c(1,1))

plot(rowMeans(y[,Sj0==1]),type='l')
lines(rowMeans(y[,Sj0==2]),col='red')
lines(rowMeans(y[,Sj0==4]),col='blue')

y1=y[order(rowSums(y),decreasing = T),]
y_s1=as.matrix(y1[,Sj==1])
y_s2=as.matrix(y1[,Sj==2])
y_s3=as.matrix(y1[,Sj==3])

par(mfrow=c(1,1))
plot(rowMeans(y_s1),type='l',xlab='OTU 1-120',ylab='Abundance',cex.lab=0.6,cex.axis=0.6)
lines(rowMeans(y_s2),type='l',col='red')
lines(rowMeans(y_s3),type='l',col='red')

pdf_s1=rowSums(y_s1)/sum(rowSums(y_s1))
pdf_s2=rowSums(y_s2)/sum(rowSums(y_s2))
pdf_s3=rowSums(y_s3)/sum(rowSums(y_s3))

plot(c(0,cumsum(pdf_s1)),type='l',ylim=c(0,1),xlim=c(0,120),xlab='OTU',ylab="CDF of Relative Abundance")
lines(c(0,cumsum(pdf_s2)),type='l',col='red')
lines(c(0,cumsum(pdf_s3)),type='l',col='blue')

### plot according to country label
data_country=pivot_wider(data[,c(1:2,6)],names_from = Sample,values_from = nationality)
data_country=data_country[1,2:dim(data_country)[2]]

y_af=as.matrix(y1[,data_country=='AFR'])
y_am=as.matrix(y1[,data_country=='AAM'])

as.matrix(data_country[Sj==1])[1,]
as.matrix(data_country[Sj==2])[1,]
as.matrix(data_country[Sj==3])[1,]

pdf_af=rowSums(y_af)/sum(rowSums(y_af))
pdf_am=rowSums(y_am)/sum(rowSums(y_am))

plot(c(0,cumsum(pdf_af)),type='l',ylim=c(0,1),xlim=c(0,120),xlab='OTU',ylab="CDF of Relative Abundance")
lines(c(0,cumsum(pdf_am)),type='l',col='red')


mik_store[T,1,order(rowSums(y),decreasing = T)]
mik_store[T,2,order(rowSums(y),decreasing = T)]
mik_store[T,3,order(rowSums(y),decreasing = T)]
w_store[T,1,]
w_store[T,2,]
w_store[T,3,]

w1_ave=sapply(1:L,function(i) w_store[T,1,i]/sum(mik_store[T,1,]==i))
mik1_plot=w1_ave[mik_store[T,1,order(rowSums(y),decreasing = T)]]

w2_ave=sapply(1:L,function(i) w_store[T,2,i]/sum(mik_store[T,2,]==i))
mik2_plot=w2_ave[mik_store[T,2,order(rowSums(y),decreasing = T)]]

w3_ave=sapply(1:L,function(i) w_store[T,3,i]/sum(mik_store[T,3,]==i))
mik3_plot=w3_ave[mik_store[T,3,order(rowSums(y),decreasing = T)]]

plot(cumsum(mik1_plot),type='l')
lines(cumsum(mik2_plot),type='l',col='red')
lines(cumsum(mik3_plot),type='l',col='blue')


###### Do David Dahl 
idx=1:T
idx2=idx[idx>7500]
TT=length(idx2)

#one_pos means positions of 1s, it is to save space for the binary matrix
one_pos=list()

for(tt in 1:TT){
  ii=idx2[tt]
  s_temp=S_store[[ii]]
  cl_temp=unique(s_temp)
  s_mat_temp=matrix(0,n,n)
  for(ss in cl_temp){
    s_mat_temp[which(s_temp==ss),which(s_temp==ss)]=1
  }
  svec=s_mat_temp[lower.tri(s_mat_temp,diag = F)]
  one_pos[[tt]]=which(svec==1)
  print(tt)
}

hem_dist=rep(0, TT)

hem_dist_mat=matrix(0,TT,TT)
for(ii in 1:TT){
  for(jj in ii:TT){
    print(c(ii,jj))
    s1=rep(0,n*(n-1)/2)
    s2=rep(0,n*(n-1)/2)
    s1[one_pos[[ii]]]=1
    s2[one_pos[[jj]]]=1
    hem_dist_mat[ii,jj]=sum(s1!=s2)*2
  }
}

hem_dist_mat2=t(hem_dist_mat)+hem_dist_mat
hem_dist = rowMeans(hem_dist_mat2)
which.min(hem_dist)

### from the cluster assignmebt and iteration picked from result above, 
### do David Dahl for the nested partitions
idx_star=idx2[which.min(hem_dist)]
s_star=S_store[idx_star,]

idx_star_all=idx[sapply(1:T, function(i) sum(S_store[i,]==s_star)==n)]
idx_star_m=idx_star_all[idx_star_all>8000]
TT=length(idx_star_m)
one_pos_m1=list()
one_pos_m2=list()
one_pos_m3=list()

for(tt in 1:TT){
  ii=idx_star_m[tt]
  s_temp=mik_store[ii,1,]
  cl_temp=unique(s_temp)
  s_mat_temp=matrix(0,B,B)
  for(ss in cl_temp){
    s_mat_temp[which(s_temp==ss),which(s_temp==ss)]=1
  }
  svec=s_mat_temp[lower.tri(s_mat_temp,diag = F)]
  one_pos_m1[[tt]]=which(svec==1)
  
  s_temp=mik_store[ii,2,]
  cl_temp=unique(s_temp)
  s_mat_temp=matrix(0,B,B)
  for(ss in cl_temp){
    s_mat_temp[which(s_temp==ss),which(s_temp==ss)]=1
  }
  svec=s_mat_temp[lower.tri(s_mat_temp,diag = F)]
  one_pos_m2[[tt]]=which(svec==1)
  
  s_temp=mik_store[ii,3,]
  cl_temp=unique(s_temp)
  s_mat_temp=matrix(0,B,B)
  for(ss in cl_temp){
    s_mat_temp[which(s_temp==ss),which(s_temp==ss)]=1
  }
  svec=s_mat_temp[lower.tri(s_mat_temp,diag = F)]
  one_pos_m3[[tt]]=which(svec==1)
  
  print(tt)
  
}


hem_dist_mat_m1=matrix(0,TT,TT)
hem_dist_mat_m2=matrix(0,TT,TT)
hem_dist_mat_m3=matrix(0,TT,TT)

for(ii in 1:TT){
  for(jj in ii:TT){
    print(c(ii,jj))
    s1=rep(0,B*(B-1)/2)
    s2=rep(0,B*(B-1)/2)
    s1[one_pos_m1[[ii]]]=1
    s2[one_pos_m1[[jj]]]=1
    #hem_dist[ii]=hem_dist[ii]+sum(s1!=s2)
    hem_dist_mat_m1[ii,jj]=sum(s1!=s2)*2
    
    s1=rep(0,B*(B-1)/2)
    s2=rep(0,B*(B-1)/2)
    s1[one_pos_m2[[ii]]]=1
    s2[one_pos_m2[[jj]]]=1
    #hem_dist[ii]=hem_dist[ii]+sum(s1!=s2)
    hem_dist_mat_m2[ii,jj]=sum(s1!=s2)*2
    
    s1=rep(0,B*(B-1)/2)
    s2=rep(0,B*(B-1)/2)
    s1[one_pos_m3[[ii]]]=1
    s2[one_pos_m3[[jj]]]=1
    #hem_dist[ii]=hem_dist[ii]+sum(s1!=s2)
    hem_dist_mat_m3[ii,jj]=sum(s1!=s2)*2
  }
  #hem_dist[ii]=hem_dist[ii]/T
}

hem_dist_mat2_m1=t(hem_dist_mat_m1)+hem_dist_mat_m1
hem_dist_m1 = rowMeans(hem_dist_mat2_m1)
which.min(hem_dist_m1)

hem_dist_mat2_m2=t(hem_dist_mat_m2)+hem_dist_mat_m2
hem_dist_m2 = rowMeans(hem_dist_mat2_m2)
which.min(hem_dist_m2)

hem_dist_mat2_m3=t(hem_dist_mat_m3)+hem_dist_mat_m3
hem_dist_m3 = rowMeans(hem_dist_mat2_m3)
which.min(hem_dist_m3)


####### form the cluster assignment, plot pdfs
mik_star=matrix(0,3,B)
mik_star[1,]=mik_store[idx_star_m[which.min(hem_dist_m1)],1,]
mik_star[2,]=mik_store[idx_star_m[which.min(hem_dist_m2)],2,]
mik_star[3,]=mik_store[idx_star_m[which.min(hem_dist_m3)],3,]


y_s1=as.matrix(y1[,s_star==1])
y_s2=as.matrix(y1[,s_star==2])
y_s3=as.matrix(y1[,s_star==3])

pdf_s1=rowSums(y_s1)/sum(rowSums(y_s1))
pdf_s2=rowSums(y_s2)/sum(rowSums(y_s2))
pdf_s3=rowSums(y_s3)/sum(rowSums(y_s3))

pdf('Sjcluster.pdf',width=4,height=3)
plot(c(0,cumsum(pdf_s1)),type='l',ylim=c(0,1),xlim=c(0,120),xlab='OTU',ylab="CDF of Relative Abundance")
lines(c(0,cumsum(pdf_s2)),type='l',col='red')
lines(c(0,cumsum(pdf_s3)),type='l',col='blue')
legend("bottomright",lty=c(1,1,1),col=c("black",'red','blue'),legend=c("cluster 1","cluster 2","cluster 3"),bty = 'n',cex=0.5)
dev.off()

library(ggplot2)
ggplot()+
  geom_line(aes(0:119,c(0,cumsum(pdf_s1)),col='cluster1'))+
  geom_line(aes(0:119,c(0,cumsum(pdf_s2)),col='cluster2'))+
  geom_line(aes(0:119,c(0,cumsum(pdf_s3)),col='cluster3'))+
  labs(x='OTU',y="Cumulative Relative Frequencies")+
  theme(legend.position = 'bottom')

ggsave("Sjclustergg.pdf",width=5,height=3.5,units = 'in')

for(j in 1:dim(y_s1)[2]){
  pdf_s1=y_s1[,j]/sum(y_s1[,j])
  points(c(0,cumsum(pdf_s1)),pch=1,cex=0.1)
}


#### plot heatmaps

data_country=pivot_wider(data[,c(1:2,6)],names_from = Sample,values_from = nationality)
data_country=data_country[1,2:dim(data_country)[2]]


y[order(mik_star[1,]),s_star==1]
y[order(mik_star[2,]),s_star==2]
y[order(mik_star[3,]),s_star==3]

od1=unlist(sapply(1:8,function (i) which(mik_star[1,]==i)))
od2=unlist(sapply(1:8,function (i) which(mik_star[2,]==i)))
od3=unlist(sapply(1:8,function (i) which(mik_star[3,]==i)))

pdf("heatmap.pdf",width=8,height=3)
par(mfrow=c(1,3))
heatmap(y[order(mik_star[1,]),s_star==1],Colv = NA, Rowv = NA)
heatmap(y[order(mik_star[2,]),s_star==2],Colv = NA, Rowv = NA)
heatmap(y[order(mik_star[3,]),s_star==3],Colv = NA, Rowv = NA)
dev.off()

pdf("heatmap1.pdf",width=4,height=3)
heatmap(y[od1,s_star==1],Colv = NA, Rowv = NA)
dev.off()

pdf("heatmap2.pdf",width=4,height=3)
heatmap(y[od2,s_star==2],Colv = NA, Rowv = NA)
dev.off()

pdf("heatmap3.pdf",width=4,height=3)
heatmap(y[od3,s_star==3],Colv = NA, Rowv = NA)
dev.off()


heatmap(y[,s_star==1],Colv = NA, Rowv = NA)
heatmap(y[,s_star==2],Colv = NA, Rowv = NA)
heatmap(y[,s_star==3],Colv = NA, Rowv = NA)


charOTU <- function(s1,s2){
  idx=sort(abs(rowMeans(y[,s_star==s1])-rowMeans(y[,s_star==s2])),decreasing=TRUE,index.return=T)$ix
  names=c(rownames(y)[idx==1],rownames(y)[idx==2],rownames(y)[idx==3])
  return(names)
}
charOTU(1,2)
charOTU(1,3)
charOTU(2,3)

Z=log(y[order(mik_star[3,]),s_star==3])

Zmelt <- reshape2::melt(Z)
glimpse(Zmelt)

ggplot(Zmelt,aes(Var2,Var1,fill=value)) +
  geom_raster() +
  scale_fill_viridis_c()

Zmelts <- rbind(
  reshape2::melt(log(y[od1,s_star==1])) %>%
    mutate(s_star = "1"),
  reshape2::melt(log(y[od2,s_star==2])) %>%
    mutate(s_star = "2"),
  reshape2::melt(log(y[od3,s_star==3])) %>%
    mutate(s_star = "3")
)

ggplot(Zmelts,aes(Var2,Var1,fill=value)) +
  geom_raster() +
  facet_wrap(~s_star)+
  scale_fill_viridis_c()


