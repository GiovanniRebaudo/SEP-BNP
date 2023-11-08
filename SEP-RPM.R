# Codes accompanying "Separate Exchangeability as 
# Modeling Principle in Bayesian Nonparametrics"

# Load relevant libraries, functions and data ----------------------------------
rm(list=ls())
# Set the working directory to the current folder 
# Code to set the working directory to the current folder from RStudio
library(rstudioapi) # version 0.15
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

source("SEP_fcts.R")


library(DirichletReg)
library(invgamma)
library(tidyr)
library(dplyr)
library(truncnorm)
library(scales)
library(viridisLite)
library(ggplot2)
library(RColorBrewer)
library(salso)

###################################################N
## create global variables for data & field dimensions
# Read and clean the data
out <- read.dta("Data-and-Results/Dietswap_dataset.RDS")
summary(out)
y <-   out$y    ## log transform
Y <-  out$Y    ## abs scale

N <- dim(y)[1]*dim(y)[2]  # global var's
n <- dim(y)[2] # # patients
B <- dim(y)[1] # # OTU's

##### max size of the subject (K) and OTU (L) partition
K = 10
L = 15
## hyperprior pars
alpha = 1      ## for GEM prior on pi
beta  = 1      ## GEM prior on w[k], k=1...K

## hyperpars (use "empirical bayes" i.e., match sample moments)
mu0  = mean(as.matrix(y))                  ## mu_l ~ N(mu0,sig0)
sig0 = 3  ## SD!
a0   = 3                                    ## sig2_l ~ Ga(a0/2, b0/2)
b0   = var(c(y))/4




## for debugging only -- keepin track of each factor in log posterior
## in the two most recent (old, new) posterior evaluations 
lpot          = c("lik","mu","sig","Sj","mki","pi","w")
lpo           = matrix(0,nrow=2,ncol=length(lpot)) ## to save most recent 2 eval of log post
colnames(lpo) = lpot
rownames(lpo) = c("old","new")


niter = 500
iter  = 0



#######################################################
## make figures,  also called at the end of ex(),
## but could be called separately

niter=1000; niter0=250; niter1=750
fig4  <- function(niter=NULL, niter0=NULL, niter1=NULL)
{# create plots 4,5 & 6
    ## Fig 4 ###############
    SjMC  <- read.myfile("Sj.txt",n)
    mu = read.myfile("mu.txt", L)
    Sj  <- SjMC[nrow(SjMC),]
    nk=sapply(1:K,function(k) sum(Sj==k))
    K1 = sum(nk>1)
    klist = which(nk>1)
    s  <- apply(Y,1,sum)
    idx <- order(s,decreasing=T)
    kdx  <- which( nk>0 )
    sk = matrix(0,nrow=K,ncol=B)
    for (k in klist){ # non-empty clusters; clusters k=1,2,3
        sk[k,] = apply(Y[,Sj==k],1,sum)
    }
    Sk = apply(sk,1,sum)
    sk = sk/Sk
    csk = apply(sk[,idx],1,cumsum) # cum sum for each k, OTU's ordered by total frequ
    lwd = nk/max(nk)*3
    matplot(1:B, csk[,klist],type="l", xlab="OTU", ylab="CUMSUM",lwd=lwd[klist],bty="l")
    legend(60,0.6, col=klist,lty=klist,legend=klist,bty="n")
    return(klist)
}

fig5 <- function(klist, niter=NULL, niter0=NULL, niter1=NULL)
{## klist = list of subject clusters for which to plot nested clusters
  ## get klist as return value from fig4()
  ## prepares Fig 5 ##
      mki2 <- read.myfile("mki.txt",K)
    summ = read.myfile("iter.txt")
    it = summ[,1]
    col = brewer.pal(9, name="YlOrRd")
    ## col = viridis(10)
    
    M  <- nrow(mki2)
    M1 = which(it >= niter1)[1]
    mki  <- array(0,dim=c(K,B,M))
    for(m in 1:M)
      for(b in 1:B)
        mki[,b,m] = mki2[m, (b-1)*K+(1:K)]
       
    pk  <- array(0,dim=c(B,B,K))
    for(k in klist){
        for(m in M1:M)
            for(b in 2:B)
                for(b2 in 1:(b-1)){
                    pk[b,b2,k] = pk[b,b2,k]+(mki[k,b,m]==mki[k,b2,m])
                }#k
        pk[,,k] = pk[,,k]/(M-M1+1)
        pk[,,k] = pk[,,k] + t(pk[,,k])
        C = t(mki[k,,M1:M])
        ck = salso(C)
        idx = order(ck)
        ## for plotting
        pkp = (pk[,,k])
        image(pkp[idx,idx])
    }# k  
}

fmap = function(x)
{
  1/(1+exp(-20*(x-0.25)))
  }
fn <- "w.txt"; p <- K; q  <-  L
read.myfile  <- function(fn, p=NULL, q=NULL)
{ # read matrix (q=NULL) or array (q>0)
    X  <- as.matrix(read.table(fn))
    if (!is.null(q)){ # read array
        niter = nrow(X)/p
        X2  <- array(c(X), dim=c(p,niter,q))
        X <- X2
    }
    return(X)
}
                     
        
    
plt.Gk  <- function(w,mu,sig2,add=F,Sj,M=100)
{ ## plots current Gk
    nk=sapply(1:K,function(k) sum(Sj==k))
    xx <- seq(from=-3,to=3,length=M)
    yy  <- matrix(0,nrow=K,ncol=M)
    for(h in 1:K){
        for(l in 1:L)
            yy[h,]  <- yy[h,] + w[h,l]*dnorm(xx,m=mu[l],sd=sqrt(sig2[l]))
    }
    lwd = nk/max(nk)*3
    if (!add){
        matplot(xx,t(yy),bty="l",type="l",lwd=lwd)
        legend(2,max(yy),legend=1:K,col=1:K,lty=1:K,bty="n",adj=1)
    } else {
    ## text(rep(xx[M],K), yy[,M], text=(1:K),adj=0)
        matlines(xx, t(yy),type="l", lwd=lwd,col=1)
    }
}


plt.Gkbar  <- function()
        {
    w  <- read.myfile("w.txt",K,L) # w[k, iter, l]
    niter  <- dim(w)[2]
    pi <- read.myfile("pi.txt",K)
    mu  <- read.myfile("mu.txt",L)
    sig2  <- read.myfile("sig.txt",L)
    Sj  <-  read.myfile("Sj.txt",K)

    ## clustering of patients
    pij  <- matrix(0,n,n)
    for(iter in 1:niter){
        Sj1  <-  Sj[1,]
    }
}

