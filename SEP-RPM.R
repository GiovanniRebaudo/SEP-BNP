# Codes accompanying "Separate Exchangeability as 
# Modeling Principle in Bayesian Nonparametrics"

# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.15
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))


library(DirichletReg)
library(invgamma)
library(tidyr)
library(dplyr)
library(truncnorm)
library(scales)
library(viridisLite)
library(ggplot2)
theme_set(theme_bw(base_size = 14))
library(RColorBrewer)
library(salso)


###################################################N
## create global variables for data & field dimensions
# Read and clean the data
###################################################N
## read in data
read.dta  <- function()
{##read-in and clean data
  Dietswap_dataset <- readRDS("Data-and-Results/Dietswap_dataset.RDS")
  
  data = Dietswap_dataset[Dietswap_dataset$timepoint==1,]
  data2=pivot_wider(data[,1:3],names_from = Sample,values_from = Abundance)
  na_in_rows=sapply(1:dim(data2)[1],function(i) sum(data2[i,]==0))
  data3=data.frame(data2[na_in_rows!=38,]) #remove all na rows
  
  
  Y0=data3[,2:39]                ## separate OTU counts
  rownames(Y0)=data3[,1]         ## and OTU names as row names..
  
  n <- dim(Y0)[2]                ## n patients
  drop = apply(Y0==0,1,sum)>n/4  ## drop all OTUS with more than n/4 0's
  Y  <- Y0[!drop,]               ## dropping 20 of the originally 119 OTU's
  ## now total # 0's is only 5.
  gammas=colMeans(Y)
  Y=sapply(1:n,function(i) as.matrix(Y)[,i]/gammas[i])
  y=log(Y+0.05)                  ## log transform (we only have 5 0's)
  return(list(y=y,Y=Y)) ## return both - the log transf "y" and raw "Y"
}
## create global variables for data & fied dimensions
out <- read.dta() ## maybe as.matrix(..)?? See if needed..
y <-   out$y    ## log transf
Y  <-  out$Y    ## abs scale

N <- dim(y)[1]*dim(y)[2]  # global var's
n <- dim(y)[2] # # patients
B <- dim(y)[1] # # OTU's

##### max size of the subject (K) and OTU (L) partition
K=10
L=15
## hyperprior pars
alpha= 1      ## for GEM prior on pi
beta = 1      ## GEM prior on w[k], k = 1, ..., K

## hyperpars (use "empirical Bayes" i.e., match sample moments)
mu0=mean(as.matrix(y))                  ## mu_l ~ N(mu0, sig0)
sig0=3  ## SD!

a0=3                                    ## sig2_l ~ Ga(a0/2, b0/2)
b0=var(c(y))/4

source("SEP_fcts.R")

############ exploratory plots ###########
Exploratory = FALSE
if(Exploratory){
Dietswap_dataset <- readRDS("Data-and-Results/Dietswap_dataset.RDS")
data = Dietswap_dataset[Dietswap_dataset$timepoint==1,]
data2=pivot_wider(data[,1:3],names_from = Sample,values_from = Abundance)
na_in_rows=sapply(1:dim(data2)[1],function(i) sum(data2[i,]==0))
data3=data.frame(data2[na_in_rows!=38,]) #remove all na rows
y=data3[,2:39]
gammas=colMeans(y)
y=sapply(1:n,function(i) as.matrix(y)[,i]/gammas[i])

y1=y[order(rowSums(y),decreasing = T),]

data_country=pivot_wider(data[,c(1:2,6)], names_from = Sample, 
                         values_from = nationality)
data_country=data_country[1,2:dim(data_country)[2]]

data_country=as.matrix(data_country)[1,]
names(data_country)=NULL

data_label=data_country
data_label[data_country=='AFR']=1
data_label[data_country=='AAM']=2

data_label2 = data_country
data_label2[data_country=='AFR'] = 'AF'
data_label2[data_country=='AAM'] = 'AA'

y_af=as.matrix(y1[,data_country=='AFR'])
y_am=as.matrix(y1[,data_country=='AAM'])


pdf_af=rowSums(y_af)/sum(rowSums(y_af))
pdf_am=rowSums(y_am)/sum(rowSums(y_am))

pdf_all=rowSums(y1)/sum(rowSums(y1))

P = ggplot()+
  geom_line(aes(0:119,c(0,cumsum(pdf_af)),col='AF'))+
  geom_line(aes(0:119,c(0,cumsum(pdf_am)),col='AA'))+
  geom_line(aes(0:119,c(0,cumsum(pdf_all)),col='ALL'))+
  labs(x='OTU', y="Cumulative Relative Frequencies")+
  theme(legend.position = 'right')+  labs(color=NULL)
ggsave(plot=P, file ="Image/Emp_CDF_OTU.pdf", 
       width=16, height=12, units = 'cm')

dfplot_af=cbind(1:119,rowMeans(y_af)) %>% data.frame()
colnames(dfplot_af) =  c("OTU","Abundace") 

P1 = ggplot(dfplot_af,aes(OTU,Abundace))+
  geom_bar(stat='identity')+
  labs(x='OTU (AF)',y='Scaled Abundance')
ggsave(plot=P1, file ="Image/normedcount-hist-af.pdf",
       width=5, height=3.5, units = 'in')

dfplot_am=cbind(1:119,rowMeans(y_am)) %>% data.frame()
colnames(dfplot_am) =  c("OTU","Abundace") 

P2 = ggplot(dfplot_am,aes(OTU,Abundace))+
  geom_bar(stat='identity')+
  labs(x='OTU (AA)',y='Scaled Abundance')
ggsave(plot=P2, file ="Image/normedcount-hist-am.pdf", 
       width=5, height=3.5, units = 'in')

#######plot dendograms based on hierarchical clustering 
colvar <- function(y){
  return(sapply(1:dim(y)[2],function(j) var(y[,j])))
}
rowvar <- function(y){
  return(sapply(1:dim(y)[1],function(i) var(y[i,])))
}
rvy=rowvar(y)
yhc=y[order(rvy,decreasing = T)<=10,]
rownames(yhc)=rownames(y)[order(rvy,decreasing = T)<=10]
colnames(yhc)=data_label

d <- dist(yhc, method = "euclidean")
hc2 <- hclust(d, method = "complete" )
plot(hc2,cex=0.5,main='')

d <- dist(t(yhc), method = "euclidean")
hc2 <- hclust(d, method = "complete" )
plot(hc2,cex=0.5,main='')

dend_expr <- as.dendrogram(hc2)
tree_labels <- dendro_data(dend_expr, type = "rectangle")
tree_labels$labels <- cbind(tree_labels$labels, Subject = as.factor(data_label2))

P = ggplot() +
  geom_segment(data =segment(tree_labels), aes(x=x, y=y, xend=xend, yend=yend))+
  geom_segment(data = tree_labels$segments %>%
                 filter(yend == 0) %>%
                 left_join(tree_labels$labels, by = "x"), 
               aes(x=x, y=y.x, xend=xend, yend=yend, color = Subject)) +
  geom_text(data = label(tree_labels), 
            aes(x=x, y=y, label=label, hjust=-1), size=2) +
  coord_flip() +
  scale_y_reverse(expand=c(0.2, 0)) +
  scale_colour_brewer(palette = "Dark2") + 
  theme_dendro()+
  theme(legend.title=element_blank())

ggsave(plot=P, file ="Image/otu-subject-hc.pdf", 
       width=5, height=3.5, units = 'in')
}


# Plot Results
# Run analysis, save results and produce some plots
Run_MCMC = FALSE
if(Run_MCMC){
  ex(niter  = 1e4, #iteration MCMC
     niter0 = 1e3, # estimating Sj after niter0 iterations and stop updating Sj
     niter1 = 2e3 # just for plot 
  )
}

# Load output
pi    <- read.myfile("pi.txt")
Sj    <-  read.myfile("Sj.txt")
w     <- read.myfile("w.txt",K,L) # w[k, iter, l]
mki   <- read.myfile(c(mki), "mki.txt", K*B, append=app)
niter <- dim(w)[2]

mu  <- read.myfile("mu.txt")
sig2  <- read.myfile("sig.txt",L)

logl <-  read.myfile("logl.txt")
logl <-  read.myfile("logl.txt")

#
write(pi, "pi.txt", K, append=app)
write(Sj, "Sj.txt", n, append=app)
write(c(w), "w.txt", L , append=app)
write(c(mki), "mki.txt", K*B, append=app)
write(mu, "mu.txt", L, append=app)
write(sig2, "sig.txt", L, append=app)
nk=sapply(1:K,function(k) sum(Sj==k))
write(nk, "nk.txt", K, append=app)
nl=sapply(1:L,function(l) sum(c(mki)==l))
write(nl, "nl.txt", L, append=app)
write(ll, "logl.txt", 1, append=app)

## write summary of the iteration to "iter.txt"
nk=sapply(1:K,function(k) sum(Sj==k))
nl= sapply(1:L, function(l) sum(c(mki)==l))
nkp=sum(nk>0)
nlp=sum(nl>0)
summ = c(iter, nkp, nlp, nk, nl, ll, pmki)
write(summ, "iter.txt", length(summ), append=app)

##### check plotting log-likelihood, cluster y mean and weights
par(mar=c(2,2,1,1))
plot(logl[-1],type='l')

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



