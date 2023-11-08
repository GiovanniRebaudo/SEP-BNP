rm(list=ls())
set.seed(135)

load("data_protein.RData")
library(dplyr)
library(tidyr)
library(mvtnorm)
library(invgamma)
library(splines)

# cleaning data and construct design matrix 
rownames(PL)=c()
PL=cbind(c(rep(0,21),rep(1,21)),PL)
colnames(PL)[1]="z"

PL2=as.data.frame(PL) %>% group_by(z,age) %>% summarise_each(funs(mean(., na.rm = TRUE))) %>% as.matrix()
PL2=PL2[,-which(sapply(1:dim(PL2)[2],function(i) all(is.na(PL[,i]))))]

age=PL2[,2]
age_std=(PL2[,2]-mean(PL2[,2]))/sd(PL2[,2])
x=cbind(PL2[,1],age_std,PL2[,1]*age_std)
y= PL2[,3:dim(PL2)[2]]

#total number of proteins
n=dim(y)[2]

# construct splines using bs package
bsx2=bs(x[,2],degree = 3,knots = c(quantile(x[,2],0.33),quantile(x[,2],0.67)),intercept = T)
x=cbind(bsx2,x[,1]*bsx2)


#truncate at N=20 clusters
tot_iter=5000
N=20
#set priors
alpha=1
V=rbeta(N,1,alpha)
pi=c(V[1],sapply(2:N,function(i) V[i]*prod(1-V[1:(i-1)])))
betasig0=diag(dim(x)[2])
beta=rmvnorm(N,mean=rep(0,dim(x)[2]),sigma = betasig0)
beta_store=list()
pi_store=list()
ll_store=NULL
sigma2_h=rinvgamma(N,1,1)
sigma2h_store=list()
alphap=rnorm(n,3,sd = sqrt(5))

#create space holders
alphap_store=list()
s_store=list()
gamma_store=matrix(0,tot_iter,n)

#create scalar age offset 
age_unq=unique(sort(age))
age_match=sapply(1:dim(y)[1], function(i) which(age_unq==PL2[i,2]))
deltasig0=0.01
delta=rnorm(length(age_unq),0,deltasig0)
delta_store=list()
delta_expand=rep(0,dim(y)[1])
for(i in 1:dim(y)[1]){
  delta_expand[i]=delta[age_match[i]]
}

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}

## functions and bins needed for Chi2-test if one wishes to use
# xch=vector()
# R=vector()
# bin=5
# binqt=c(0.2,0.4,0.6,0.8)
# zij <- function(yij,x){
#   if(yij <= x[1]){
#     return(1)
#   }
#   if(yij > x[1] & yij <= x[2]){
#     return(2)
#   }
#   if(yij > x[2] & yij <= x[3]){
#     return(3)
#   }
#   if(yij > x[3] & yij <= x[4]){
#     return(4)
#   }
#   if(yij > x[4]){
#     return(5)
#   }
# }
# this is for normal-QQ test if one wishes to do so
zip=rep(0,tot_iter)

for(iter in 1:tot_iter){
  
  ## update cluster assignment s
  mu=x %*% t(beta)
  
  ll=matrix(0,n,N)
  s=vector(length = n,mode="numeric")
  
  ll_store[iter]=0
  #rseq= sample(n,n,replace = F)
  for(j in 1:n){
    for(nb in 1:N){
      for(mm in 1:dim(y)[1]){
        if(!is.na(y[mm,j])){
          ll[j,nb]= ll[j,nb]+ log(dnorm(y[mm,j]-alphap[j]-delta_expand[mm],mu[mm,nb],sqrt(sigma2_h[nb]) ))
        }
      }
      ll[j,nb]= ll[j,nb]+log(pi[nb])
     }
    #s[j]=sample(1:N,size=1,prob = exp(ll[j,]))
    #https://en.wikipedia.org/wiki/Categorical_distribution#Sampling_via_the_Gumbel_distribution
    gumble=-log(-log(runif(N,0,1)))
    s[j]=which.max(ll[j,]+gumble)
    ll_store[iter]=ll_store[iter]+ll[j,s[j]]
  }
  #note that ll_store[iter] corresponds to beta_store[[iter-1]]
  
  s_store[[iter]]=s
  
  ###### update beta in each cluster
  ns=sapply(1:N,function(i) sum(s==i))
  for(g in 1:N){
    if(sum(s==g)==0){
      next
    }else{
      cc=which(s==g)
      yg=matrix(0,dim(y)[1],length(cc))
      for(rr in 1:length(cc)){
        tt=cc[rr]
        yg[,rr]=y[,tt]-alphap[tt]-delta_expand
      }
      yg=as.vector(yg)
      yg=yg[!is.na(yg)]
      # this part is for formatting designing matrix in each cluster
      xg=matrix(0,length(yg),dim(x)[2])
      rxg=1
      for(tt in cc){
        for(mm in 1:dim(y)[1]){
          if(!is.na(y[mm,tt])){
            xg[rxg,]=x[mm,]
            rxg=rxg+1
          }
        }
      }
      sigma_n_inv=rep(1/sigma2_h[g],length(yg))
      sigtemp= solve(solve(betasig0)+ t(xg) %*% (sigma_n_inv * xg))
      sigtemp= makeSymm(sigtemp)
      #beta prior mean 0, don't write out 
      meantemp= sigtemp %*% (t(xg) %*% (sigma_n_inv * yg))
      beta[g,]=rmvnorm(1,mean=meantemp,sigma = sigtemp)
      sigma2_h[g]=rinvgamma(1,1+ns[g]*dim(y)[1]/2,1+sum((yg-xg %*% beta[g,])^2)/2)
      # gamma calculation
      gamma_store[iter,which(s==g)]=x[32,7:12] %*% meantemp[7:12]-x[17,7:12] %*% meantemp[7:12]
    }
  }
  
  beta_store[[iter]]=beta
  sigma2h_store[[iter]]=sigma2_h
  
  #update alpha_p
  for(j in 1:n){
    yj=y[!is.na(y[,j]),j]
    nj=length(yj)
    xj=x[!is.na(y[,j]),]
    dj=delta_expand[!is.na(y[,j])]
    alphap[j]=rnorm(1,(sum(yj-xj %*% beta[s[j],]-dj)/sigma2_h[s[j]])/(1/1+dim(y)[1]/sigma2_h[s[j]]),sqrt(1/(1/1+dim(y)[1]/sigma2_h[s[j]])))
  }
  alphap_store[[iter]]=alphap
  
  #update delta 
  for(t in 1:length(age_unq)){
    ti=which(age_match==t)
    
    yi_d=list()
    idx=list()
    sigs=list()
    ni= rep(0,length(ti))
    for(tt in 1:length(ti)){
      i=ti[tt]
      yi_vec=y[i,]-alphap- x[i,] %*% t(beta[s,])
      yi_d[[tt]]=yi_vec[!is.na(yi_vec)]
      idx[[tt]]=!is.na(yi_vec)
      ni[tt]=length(yi_d[[tt]])
      s2=s[!is.na(y[i,])]
      sigs[[tt]]=sigma2_h[s2]
    }
    
    yi_d=unlist(yi_d)
    sigs=unlist(sigs)
    delta[t]=rnorm(1,sum(yi_d/sigs)/(1/deltasig0^2+sum(1/sigs)),sqrt(1/(1/deltasig0+sum(1/sigs))))
  }
  delta_store[[iter]]=delta
  
  #expand delta as it was for unique ages -- expand it to all y
  delta_expand=rep(0,dim(y)[1])
  for(i in 1:dim(y)[1]){
    delta_expand[i]=delta[age_match[i]]
  }
  
  # update pi
  V=c(sapply(1:(N-1),function(i) rbeta(1,1+ns[i],alpha+sum(ns[(i+1):N]))),1)
  pi=c(V[1],sapply(2:N,function(i) V[i]*prod(1-V[1:(i-1)])))
  pi_store[[iter]]=pi
  
  
  ## chi2-test if one wishes to do 
  # xch[iter]=rchisq(1,bin-1)
  # nn=length(y)
  # z=matrix(0,nn,bin)
  # tt=1
  # #cl=as.numeric(names(table(s)))
  # #xdelta=sapply(1:dim(y)[1], function(mm) x[mm,] %*% delta_expand[mm,])
  # for(j in 1:n){
  #   yj=y[,j]-alphap[j]-delta_expand
  #   yj=yj[!is.na(yj)]
  #   nj=length(yj)
  #   p=exp(ll[j,])
  #   p=p/sum(p)
  #   mj= rowSums(sapply(1:dim(beta)[1], function(gg) p[gg]*x %*% beta[gg,])) #xx %*% beta[s[j],]
  #   stdj=sqrt(sum(p^2*sigma2_h))
  #   #mj=xx %*% beta[s[j],]
  #   #stdj=sqrt(sigma2_h[s[j]])
  #   for(ii in 1:nj){
  #     yqt= qnorm(binqt,mj[ii],stdj)
  #     z[tt,zij(yj[ii],yqt)]=1
  #     tt=tt+1
  #   }
  # }
  # 
  # mk=colSums(z)
  # R[iter]=sum(((mk-nn*1/bin)/sqrt(nn*1/bin))^2)
  # 
  
  ####### for normal QQ test #######
  zip_i=sample(1:dim(y)[1],1)
  zip_p=sample(1:n,1)
  mip=alphap[zip_p]+x[zip_i,] %*% beta[s[zip_p],]+delta_expand[zip_i]
  zip[iter]=(y[zip_i,zip_p]-mip)/sqrt(sigma2_h[s[zip_p]])
  
  print(iter)
}


zzip=zip[seq(0.5*tot_iter,tot_iter,5)]
zzip=zzip[!is.na(zzip)]


plot(ll_store[1000:5000],type="l")

###########################################
#Do David Dahl
#library(Matrix)
#s_mat=Matrix(data=0, length(idx), ncol=n^2, sparse = T)
K_store=vector()
for(iter in 1:tot_iter){
  K_store[iter]=length(unique(s_store[[iter]]))
}
table(K_store)


idx=which(K_store== as.numeric(names(which.max(table(K_store[2000:3000])))))

idx2=idx[idx>4500]
TT=length(idx2)
one_pos=list()

for(tt in 1:TT){
  ii=idx2[tt]
  s_temp=s_store[[ii]]
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

#######################################
# resample variables based on cluster number result from above
idx_star=idx2[which.min(hem_dist)]
s_star=s_store[[idx_star]]

nc=length(unique(s_star))
tb=table(s_star)
ns=as.vector(tb)
cl_star= as.numeric(names(tb))

TT=1000

alpha_star=alphap_store[[idx_star]]
sigma2h_star=sigma2h_store[[idx_star]][cl_star]

beta_star_post=array(0,dim=c(TT,nc,dim(x)[2]))
beta_star_post[1,,]=beta_store[[idx_star]][cl_star]

sigma2h_star_post=matrix(0,TT,nc)
sigma2h_star_post[1,]=sigma2h_star

alpha_star_post=matrix(0,TT,n)
alpha_star_post[1,]=alpha_star

delta_star_post=matrix(0,TT,length(age_unq))
delta_star_post[1,]=delta_store[[idx_star]]

delta_star_post_expand=matrix(0,TT,dim(y)[1])
for(i in 1:dim(y)[1]){
  delta_star_post_expand[1,i]=delta_star_post[1,age_match[i]]
}

for(qq in 1:(TT-1)){
  
  for(g in 1:nc){
    yg=as.vector(y[,which(s_star==cl_star[g])]-alpha_star_post[qq,which(s_star==cl_star[g])]-delta_star_post_expand[qq,])
    yg=yg[!is.na(yg)]
    xg=matrix(0,length(yg),dim(x)[2])
    cc=which(s_star==cl_star[g])
    rxg=1
    for(tt in cc){
      for(mm in 1:dim(y)[1]){
        if(!is.na(y[mm,tt])){
          xg[rxg,]=x[mm,]
          rxg=rxg+1
        }
      }
    }
    sigma_n_inv=rep(1/sigma2h_star_post[qq,g],length(yg))
    sigtemp= solve(solve(betasig0)+ t(xg) %*% (sigma_n_inv * xg))
    sigtemp= makeSymm(sigtemp)
    meantemp= sigtemp %*% (t(xg) %*% (sigma_n_inv * yg))
    beta_star_post[qq+1,g,]=rmvnorm(1,mean=meantemp,sigma = sigtemp)
    sigma2h_star_post[qq+1,g]=rinvgamma(1,1+length(yg)/2,1+sum((yg-xg %*% beta_star_post[qq+1,g,])^2)/2)
    
  }
  
  for(j in 1:n){
    sj=which(cl_star==s_star[j])
    yj=y[!is.na(y[,j]),j]
    nj=length(yj)
    xj=x[!is.na(y[,j]),]
    dj=delta_star_post_expand[qq,!is.na(y[,j])]
    alpha_star_post[qq+1,j]=rnorm(1,(5/5+sum(yj-xj %*% beta_star_post[qq+1,sj,] -dj)/sigma2h_star_post[qq+1,sj])/(1/5+nj/sigma2h_star_post[qq+1,sj]),sqrt(1/(1/5+nj/sigma2h_star_post[qq+1,sj])))
    
  }
  
  sjs=sapply(1:n, function(j) which(cl_star==s_star[j]))
  for(t in 1:length(age_unq)){
    ti=which(age_match==t)
    yi_d=list()
    idx=list()
    sigs=list()
    ni= rep(0,length(ti))
    for(tt in 1:length(ti)){
      i=ti[tt]
      yi_vec=y[i,]-alpha_star_post[qq+1,]- x[i,] %*% t(beta_star_post[qq+1,sjs,])
      yi_d[[tt]]=yi_vec[!is.na(yi_vec)]
      idx[[tt]]=!is.na(yi_vec)
      ni[tt]=length(yi_d[[tt]])
      s2=sjs[!is.na(y[i,])]
      sigs[[tt]]=sigma2h_star_post[qq+1,s2]
    }
    
    yi_d=unlist(yi_d)
    sigs=unlist(sigs)
    delta_star_post[qq+1,t]=rnorm(1,sum(yi_d/sigs)/(1/deltasig0^2+sum(1/sigs)),sqrt(1/(1/deltasig0^2+sum(1/sigs))))
  }
  
  for(i in 1:dim(y)[1]){
    delta_star_post_expand[qq+1,i]=delta_star_post[qq+1,age_match[i]]
  }
  
  
  print(qq)
}
beta_star_post_mean=apply(beta_star_post,c(2,3),mean)
alpha_star_post_mean=colMeans(alpha_star_post)
sigma2h_star_post_mean=colMeans(sigma2h_star_post)
delta_star_post_mean=colMeans(delta_star_post)
delta_star_post_expand_mean=colMeans(delta_star_post_expand)


##### calculate R^2 in each cluster
R_sq=rep(0,nc)

muc=matrix(0,nc,32)
for(i in 1:nc){
  muc[i,]=x %*% beta_star_post_mean[i,]
}

xdelta=sapply(1:32, function(i) delta_star_post_expand_mean[i])

for(ii in 1:length(cl_star)){
  idx_cl=which(s_star==cl_star[ii])
  
  yhat=matrix(0,32,length(idx_cl))
  for(i in 1:length(idx_cl)){
    yhat[,i]=x %*% beta_star_post_mean[ii,] + xdelta + rep(alpha_star_post_mean[idx_cl[i]],dim(x)[1])
  }
  
  
  ssr_cl=sum((y[,idx_cl]- yhat)^2,na.rm = T)
  sst_cl=sum((y[,idx_cl]- mean(y[,idx_cl],na.rm=T))^2,na.rm=T)
  R_sq[ii]=1-ssr_cl/sst_cl
}
R_sq

#### calculate mu and ymean in each cluster
muc=matrix(0,nc,32)
for(i in 1:nc){
  muc[i,]=x %*% beta_star_post_mean[i,]
}

a_mat=matrix(0,32,n)
for(i in 1:n){
  a_mat[,i]=alpha_star_post_mean[i]
}

ymean0=matrix(0,32,nc)
for(i in 1:nc){
  pts=which(s_star==cl_star[i])
  if(length(pts)>1){
    ymean0[,i]=rowMeans(y[,pts],na.rm = T)
  }else{
    ymean0[,i]=y[,pts]
  }
}

muc0=matrix(0,nc,32)
for(i in 1:nc){
  pts=which(s_star==cl_star[i])
  muc0[i,]=x %*% beta_star_post_mean[i,] + xdelta  +mean(alpha_star_post_mean[pts])
}

###### plotting 
col_list=c(1:8,'orange','purple','dark green','maroon','dark blue','dark gray','brown')
name_list=c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)","(i)","(j)","(k)","(l)")

pdf("fit-splines-nosplit-moreb-deltai-t-scalar-5000-2.pdf",width = 9,height = 12)
par(mfrow=c(4,3))
for(i in 1:nc){
  plot(PL2[1:16,2],muc0[i,1:16],type='l',col=col_list[i],ylim=c(-4,8),xlab='Age',ylab='(log) Expected Protein Abundance',main=substitute(paste(nn,' ',R^2,'=',mm), list(nn=name_list[i],mm=round(R_sq[i],3))),cex.main= 1)
  lines(PL2[17:32,2],muc0[i,17:32],lty=2,col=col_list[i])
  lines(PL2[1:16,2],ymean0[1:16,i],col=col_list[i],lty=3)
  lines(PL2[17:32,2],ymean0[17:32,i],col=col_list[i],lty=4)
  legend("topright", lty=c(1,2,3,4),legend=c("patient(fit)","control(fit)","patient(obs)","control(obs)"),col=col_list[i],bty = 'n',cex=0.8)
  }
dev.off()
par(mfrow=c(1,1))

nc_order=c(11,1,3,4,5,6,7,10,9,8,2)
pdf("fit-splines-nosplit-moreb-deltai-t-scalar-5000-2-adjod-pcflip-2.pdf",width = 7.5,height = 10)
par(mfrow=c(4,3))
for(ii in 1:nc){
  i=nc_order[ii]
  plot(PL2[17:32,2],muc0[i,17:32],type='l',col=col_list[ii],ylim=c(-4,8),xlab='Age',ylab='(log) Expected Protein Abundance',main=substitute(paste(nn,' ',R^2,'=',mm), list(nn=name_list[ii],mm=round(R_sq[i],3))),cex.main= 1)
  lines(PL2[1:16,2],muc0[i,1:16],lty=2,col=col_list[ii])
  lines(PL2[17:32,2],ymean0[17:32,i],col=col_list[ii],lty=3)
  lines(PL2[1:16,2],ymean0[1:16,i],col=col_list[ii],lty=4)
  legend("topright", lty=c(1,2,3,4),legend=c("patient(fit)","control(fit)","patient(obs)","control(obs)"),col=col_list[ii],bty = 'n',cex=0.8)
  }
dev.off()
par(mfrow=c(1,1))

############### rank gamma based on 0-1 Loss ####################
my_gamma_postmean=colMeans(gamma_store[4501:tot_iter,])
top_y = vector()
for(i in 1:100){
  top_y[i]=colnames(y)[which(rank(my_gamma_postmean, ties.method = 'first')==i)]
}

gamma_store_abs= abs(gamma_store)

rank_mat=matrix(0,500,n)
for(i in 4501:tot_iter){
  rank_mat[i-4500,]=rank(-gamma_store_abs[i,], ties.method = 'first')
}

top_y_sel=vector()
for(i in 1:100){
  top_y_sel[i]=colnames(y)[which(rank(colMeans(rank_mat),ties.method = 'first')==i)]
}


rank01ls=1/500*colSums(rank_mat/(n+1)>0.975)
rank(rank01ls,ties.method = 'first')
top_y_01ls=vector()
for(i in 1:100){
  top_y_01ls[i]=colnames(y)[which(rank(rank01ls,ties.method = 'first')==i)]
}

######### Normal- QQtest #####################
yhat=matrix(0,dim(y)[1],n)
zip=matrix(0,dim(y)[1],n)
for(i in 1:dim(y)[1]){
  for(p in 1:n){
    c=which(cl_star==s_star[p])
    yhat[i,p]=alpha_star_post[1,p]+ x[i,]%*%beta_star_post[1,c,]+delta_star_post_expand[1,i]
    zip[i,p]=(y[i,p]-yhat[i,p])/sigma2h_star_post[1,c]
  }
}
zip=as.vector(zip[!is.na(zip)])
NN=length(zip)

pdf("QQ-plot.pdf",width=4,height=3)
qqnorm(zip, pch = 1, frame = FALSE,ylim = c(-4,4),xlim=c(-4,4))
lines(seq(-4,4,0.1),seq(-4,4,0.1))
qqline(zip, col = "steelblue", lwd = 2)
dev.off()

library(car)
pdf("QQ-plot-2.pdf",width=4,height=3)
qqPlot(zip)
dev.off()

pdf("QQ-plot-noref.pdf",width=4,height=3)
qqnorm(zip, pch = 1, frame = FALSE)
dev.off()

